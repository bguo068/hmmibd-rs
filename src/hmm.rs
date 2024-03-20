use crate::{
    data::{FracRecord, InputData, OutputBuffer, SegRecord},
    matrix::AsOption,
    model::*,
};

pub struct HmmRunner<'a> {
    data: &'a InputData,
}

impl<'a> HmmRunner<'a> {
    pub fn new(data: &'a InputData) -> Self {
        Self { data }
    }

    pub fn run_hmm_on_pair(
        &self,
        pair: (usize, usize),
        a: &mut Vec<[[f64; 2]; 2]>,
        b: &mut Vec<[f64; 2]>,
        out: &mut OutputBuffer<'a>,
        suppress_frac: bool,
    ) {
        let mut ms = ModelParamState::new(self.data);
        let mut cv = PerChrModelVariables::new();
        let mut rs = RunningStats::new();
        let total_nsites = self.data.sites.get_pos_cm_slice().len();
        a.clear();
        b.clear();
        a.resize(total_nsites, [[0.0; 2]; 2]);
        b.resize(total_nsites, [0.0; 2]);

        for iiter in 0..self.data.args.max_iter {
            ms.iiter = iiter as usize;
            ms.model.max_phi = 0.0;
            rs = RunningStats::new();

            self.run_over_fit_iterations(pair, &mut ms, &mut rs, &mut cv, a, b, out);
            if ms.is_fit_finished() {
                break;
            }
            // - update ms and finish_fit
            self.update_model(&rs, &mut ms);
        }
        if suppress_frac {
            return;
        }
        // print fraction of allele IBD for a given pair
        self.print_hmm_frac(&mut ms, pair, &mut rs, out);
    }

    pub fn run_over_fit_iterations(
        &self,
        pair: (usize, usize),
        ms: &mut ModelParamState,
        rs: &mut RunningStats,
        cv: &mut PerChrModelVariables,
        a: &mut Vec<[[f64; 2]; 2]>,
        b: &mut Vec<[f64; 2]>,
        out: &mut OutputBuffer<'a>,
    ) {
        let nchrom = self.data.genome.get_nchrom();
        for chrid in 0..nchrom as usize {
            self.run_over_single_chromsome(pair, chrid, ms, rs, cv, a, b, out);
        }
    }

    pub fn run_over_single_chromsome(
        &self,
        pair: (usize, usize),
        chrid: usize,
        ms: &mut ModelParamState,
        rs: &mut RunningStats,
        cv: &mut PerChrModelVariables,
        a: &mut Vec<[[f64; 2]; 2]>,
        b: &mut Vec<[f64; 2]>,
        out: &mut OutputBuffer<'a>,
    ) {
        // skip empty chromosomes
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);
        if start_chr == end_chr {
            return;
        }
        cv.resize_and_clear(self.data, chrid);
        // iterator over snp forward to calculate
        // - a/b
        // - alpha
        // - delta
        // - psi
        // - max_phiL
        self.run_over_snp_1_fwd(pair, chrid, ms, cv, a, b);

        // iterator over snp backward to calcuate
        // - beta
        self.run_over_snp_2_back(chrid, cv, a, b);

        // - traj
        self.run_over_snp_3_back(chrid, cv);

        // iterator over snp forward again to calcualte
        // - xi
        // - gamma
        // - update rs
        self.run_over_snp_4_fwd(chrid, rs, cv, a, b);

        if (ms.iiter == (self.data.args.max_iter - 1) as usize) || ms.finish_fit {
            self.print_final_viterbi_trajectory(pair, chrid, rs, cv, out, ms);
            self.tabulate_sites_by_state_for_viterbi_traj(chrid, rs, cv);
        }
    }
    pub fn run_over_snp_1_fwd(
        &self,
        pair: (usize, usize),
        chrid: usize,
        ms: &mut ModelParamState,
        cv: &mut PerChrModelVariables,
        a_gw: &mut Vec<[[f64; 2]; 2]>,
        b_gw: &mut Vec<[f64; 2]>,
    ) {
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);
        let geno_x = &self.data.geno[pair.0][start_chr..end_chr];
        let geno_y = &self.data.geno[pair.1][start_chr..end_chr];
        let args = &self.data.args;
        let freq1_all = &self.data.freq1;
        let freq2_all = match self.data.freq2.as_ref() {
            Some(freq2) => freq2,
            None => freq1_all,
        };
        let freq1 = freq1_all.get_row_chunk_view(start_chr, end_chr);
        let freq2 = freq2_all.get_row_chunk_view(start_chr, end_chr);
        let pi = &ms.model.pi;
        let nall = &self.data.nall[start_chr..end_chr];
        let majall = &self.data.majall[start_chr..end_chr];
        let cm = &self.data.sites.get_pos_cm_slice()[start_chr..end_chr];
        let nsites = end_chr - start_chr;
        let a_chr = &mut a_gw[start_chr..end_chr];
        let b_chr = &mut b_gw[start_chr..end_chr];

        for t in 0..nsites {
            let g_xt = geno_x[t];
            let g_yt = geno_y[t];
            let eps = args.eps;
            // let mut sv = PerSnpModelVariables::new();

            // get num of informative sites
            if (ms.iiter == 0) && (g_xt.is_some()) && (g_yt.is_some()) {
                if (g_xt == g_yt) && (g_yt == majall[t]) {
                } else {
                    ms.num_info_site += 1;
                }
                if g_xt != g_yt {
                    ms.num_discord += 1;
                }
            }

            // update b for a given sites (full loop will update whole chromosome)
            let b = &mut b_chr[t];
            if ms.iiter == 0 {
                let pright = 1.0 - eps * (nall[t] - 1) as f64;
                let (b0t, b1t) = if (g_xt == u8::MAX) || (g_yt == u8::MAX) {
                    // missing
                    (1.0, 1.0)
                } else {
                    let g_xt = g_xt as usize;
                    let g_yt = g_yt as usize;
                    let f1_x = freq1[t][g_xt].as_option().unwrap();
                    let f1_y = freq1[t][g_yt].as_option().unwrap();
                    let f2_x = freq2[t][g_xt].as_option().unwrap();
                    let f2_y = freq2[t][g_yt].as_option().unwrap();
                    if g_xt == g_yt {
                        // concordant genotype
                        // Schaffner's Additional File 1: Eq 2 and 3
                        cv.nsites += 1;
                        let fmean = (f1_x + f2_y) / 2.0;
                        let b0t = pright * pright * fmean + eps * eps * (1.0 - fmean);
                        let b1t = pright * pright * f1_x * f2_y
                            + pright * eps * f1_x * (1.0 - f2_y)
                            + pright * eps * (1.0 - f1_x) * f2_y
                            + eps * eps * (1.0 - f1_x) * (1.0 - f2_y);
                        (b0t, b1t)
                    } else {
                        // discordant genotype
                        // Schaffner's Additional File 1: Eq 4 and 5
                        cv.nsites += 1;
                        let fmeani = (f1_x + f2_x) / 2.0;
                        let fmeanj = (f1_y + f2_y) / 2.0;
                        let b0t =
                            pright * eps * (fmeani + fmeanj) + eps * eps * (1.0 - fmeani - fmeanj);
                        let b1t = pright * pright * f1_x * f2_y
                            + pright * eps * (f1_x * (1.0 - f2_y) + f2_y * (1.0 - f1_x))
                            + eps * eps * (1.0 - f1_x) * (1.0 - f2_y);
                        (b0t, b1t)
                    }
                };
                b[0] = b0t;
                b[1] = b1t;
            }

            // calcuate alpha, delta(phi) and psi
            if t == 0 {
                // initiation
                // Rabiner's Eq 32b
                cv.psi[0][t] = 0;
                cv.psi[1][t] = 0;
                // Rabiner's Eq 32a (in log scale)
                cv.phi[0][t] = pi[0].ln() + b[0].ln();
                cv.phi[1][t] = pi[1].ln() + b[1].ln();
                // Rabiner's Eq 19
                cv.alpha[0][t] = pi[0] * b[0];
                cv.alpha[1][t] = pi[1] * b[1];
                // Scale init
                cv.scale[0] = 1.0;
            } else {
                // only recalculate for each  pair and iteration in the 1st fwd loop

                let a = &mut a_chr[t - 1];
                // induction
                // Schaffer's Additional File 1 Eq 1
                let ptrans = ms.model.k_rec * (cm[t] - cm[t - 1]) as f64 / 100.0;
                a[0][1] = pi[1] - pi[1] * (-ptrans).exp();
                a[1][0] = pi[0] - pi[0] * (-ptrans).exp();
                a[0][0] = 1.0 - a[0][1];
                a[1][1] = 1.0 - a[1][0];

                // Rabiner's Eq 20 for (alpha)
                // Rabiner's Eq 33a for (phi)
                // Rabiner's Eq 33b for (psi)
                for i_o in 0..2 {
                    let mut max_val = f64::MIN;
                    let mut alpha_it = 0.0;

                    cv.scale[t] = 0.0;
                    for i_i in 0..2 {
                        // Rabinar's Eq 33b (target)
                        let score = cv.phi[i_i][t - 1] + a[i_i][i_o].ln();
                        if score > max_val {
                            max_val = score;
                            // Rabiner's Eq 33b for (psi)
                            cv.psi[i_o][t] = i_i as u8;
                        }
                        // Rabiner's Eq 33a for (phi)
                        cv.phi[i_o][t] = max_val + b[i_o].ln();
                        // Rabiner's Eq 20 for (alpha, part)
                        alpha_it += cv.alpha[i_i][t - 1] * a[i_i][i_o] * b[i_o];
                    }
                    // Rabiner's Eq 20 for (alpha, sum of part)
                    cv.alpha[i_o][t] = alpha_it;
                    // Scale addup
                    cv.scale[t] += cv.alpha[i_o][t];
                }
                // scaling
                for i in 0..2 {
                    cv.alpha[i][t] = cv.alpha[i][t] / cv.scale[t];
                }
            }
        }

        // Rabiner's Eq 34a and 34b
        let mut max_phi_local = cv.phi[1][nsites - 1];
        let mut max = 1;
        if cv.phi[1][nsites - 1] < cv.phi[0][nsites - 1] {
            // Rabiner's Eq 34a
            max_phi_local = cv.phi[0][nsites - 1];
            // Rabiner's Eq 34b
            max = 0;
        }
        // Rabiner's Eq 34b
        cv.traj[nsites - 1] = max;

        // Schaffner Additional File 1, page 4: "... positions on different
        // chromosomes are considered infinitely separated. In the code, the
        // latter is achieved by fitting data for different chromosomes
        // separately under a given iteration of the model (in other words, data
        // from all chromosomes are fit using common values of π and k) ..."
        ms.model.max_phi += max_phi_local;
    }

    /// Cacluating Beta
    pub fn run_over_snp_2_back(
        &self,
        chrid: usize,
        cv: &mut PerChrModelVariables,
        a_gw: &Vec<[[f64; 2]; 2]>,
        b_gw: &mut Vec<[f64; 2]>,
    ) {
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);
        let nsites = end_chr - start_chr;
        let a_chr = &a_gw[start_chr..end_chr];
        let b_chr = &mut b_gw[start_chr..end_chr];

        // init beta
        // Rabiner's Eq 24
        cv.beta[0][nsites - 1] = 1.0;
        cv.beta[1][nsites - 1] = 1.0;
        // induction
        for t in (0..nsites - 1).rev() {
            // let mut sv = PerSnpModelVariables::new();
            let a = a_chr[t];
            let b = b_chr[t + 1];

            // Rabiner's Eq 25
            let mut sum = [0.0, 0.0];
            for (i_i, i_o) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
                sum[i_i] += cv.beta[i_o][t + 1] * a[i_i][i_o] * b[i_o];
            }
            cv.beta[0][t] = sum[0] / cv.scale[t];
            cv.beta[1][t] = sum[1] / cv.scale[t];
        }
    }
    /// backtracking (finding traj/q sequence)
    pub fn run_over_snp_3_back(&self, chrid: usize, cv: &mut PerChrModelVariables) {
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);
        let nsites = end_chr - start_chr;

        // Rabiner's Eq 35
        let mut max = cv.traj[nsites - 1] as usize;
        for t in (0..nsites - 1).rev() {
            cv.traj[t] = cv.psi[max][t + 1];
            max = cv.traj[t] as usize;
        }
    }

    pub fn run_over_snp_4_fwd(
        &self,
        chrid: usize,
        // ms: &ModelParamState,
        rs: &mut RunningStats,
        cv: &PerChrModelVariables,
        a_gw: &Vec<[[f64; 2]; 2]>,
        b_gw: &mut Vec<[f64; 2]>,
    ) {
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);
        let pos = &self.data.sites.get_pos_slice()[start_chr..end_chr];
        let nsites = end_chr - start_chr;
        let mut sv = PerSnpModelVariables::new();
        let a_chr = &a_gw[start_chr..end_chr];
        let b_chr = &b_gw[start_chr..end_chr];

        for t in 0..nsites {
            // gamma => Rabiner's Eq 28
            let mut p_ibd = cv.alpha[0][t] * cv.beta[0][t];
            let p_dbd = cv.alpha[1][t] * cv.beta[1][t];
            p_ibd = p_ibd / (p_ibd + p_dbd);

            // Part of Schaffner's Additional File 1: model update via Baum-Welch
            // Similar to Rabiner's Eq 39a and 40a, also related to Eq 43a
            rs.count_ibd_fb += p_ibd;
            rs.count_dbd_fb += 1.0 - p_ibd;

            // not the last snp for following code
            if t == nsites - 1 {
                break;
            }

            let a = a_chr[t];
            // Rabiner's Eq 37
            let b = b_chr[t + 1];
            let mut xisum = 0.0;
            for (i_i, i_o) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
                sv.xi[i_i][i_o] = cv.alpha[i_i][t] * a[i_i][i_o] * b[i_o] * cv.beta[i_o][t + 1];
                xisum += sv.xi[i_i][i_o];
            }
            for (is, js) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
                sv.xi[is][js] /= xisum;
            }

            // Rabiner's Eq 38
            sv.gamma[0] = sv.xi[0][1] + sv.xi[0][0];
            sv.gamma[1] = sv.xi[1][1] + sv.xi[1][0];

            let delpos = (pos[t + 1] - pos[t]) as f64;
            rs.seq_ibd_fb += delpos * sv.xi[0][0];
            rs.seq_dbd_fb += delpos * sv.xi[1][1];

            // Part of Schaffner's Additional File 1: model update via Baum-Welch
            // Also similar to Rabiner's 40b (trans_obs) and 40c (trans_pred)
            rs.trans_obs += sv.xi[0][1] + sv.xi[1][0];
            rs.trans_pred += sv.gamma[0] * a[0][1] + sv.gamma[1] * a[1][0];
        }
    }

    pub fn print_final_viterbi_trajectory(
        &self,
        pair: (usize, usize),
        chrid: usize,
        rs: &mut RunningStats,
        cv: &mut PerChrModelVariables,
        out: &mut OutputBuffer<'a>,
        ms: &ModelParamState,
    ) {
        // if no useful sites, no segment generated
        if cv.nsites == 0 {
            return;
        }
        // skip if too old
        if let Some(max_tmrca) = self.data.args.filt_max_tmrca {
            if ms.model.k_rec > max_tmrca {
                return;
            }
        }
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);
        let pos = &self.data.sites.get_pos_slice()[start_chr..end_chr];
        let cm = &self.data.sites.get_pos_cm_slice()[start_chr..end_chr];
        let nsites = end_chr - start_chr;

        let chrname = self.data.genome.get_chrname(chrid);
        let sample1 = &self.data.samples.v()[pair.0];
        let sample2 = &self.data.samples.v()[pair.1];
        let gw_chr_start_bp = self.data.genome.get_gwchrstarts()[chrid];

        let mut t_seg_start = 0;
        let mut add_seq;
        //
        //      |-----|==========|-------|===========|
        //          start       -1 \isnp
        for t in 1..nsites {
            if cv.traj[t] == cv.traj[t - 1] {
                continue;
            }
            rs.ntrans += 1;

            add_seq = (pos[t - 1] - pos[t_seg_start] + 1) as f64;
            match cv.traj[t - 1] == 0 {
                true => rs.seq_ibd += add_seq,
                false => rs.seq_dbd += add_seq,
            };
            // print
            let ibd = cv.traj[t - 1];
            let seg = SegRecord {
                sample1,
                sample2,
                chrname,
                start_pos: pos[t_seg_start] - gw_chr_start_bp,
                end_pos: pos[t - 1] - gw_chr_start_bp,
                ibd,
                n_snp: t - t_seg_start,
            };

            if self.data.args.filt_ibd_only && (ibd == 1) {
                t_seg_start = t;
                continue;
            }

            // skip if too short
            if let Some(min_seg_cm) = self.data.args.filt_min_seg_cm {
                let cm = cm[t - 1] - cm[t_seg_start];
                if cm < min_seg_cm as f64 {
                    t_seg_start = t;
                    continue;
                }
            }
            out.add_seg(seg);
            t_seg_start = t;
        }
        // Process the "hanging" segments
        add_seq = (nsites + 1) as f64;
        let ibd = cv.traj[nsites - 1];
        match ibd == 0 {
            true => rs.seq_ibd += add_seq,
            false => rs.seq_dbd += add_seq,
        };

        // if output IBD only and segment is not IBD
        if self.data.args.filt_ibd_only && (ibd == 1) {
            return;
        }
        // skip if too short
        if let Some(min_seg_cm) = self.data.args.filt_min_seg_cm {
            let cm = cm[nsites - 1] - cm[t_seg_start];
            if cm < min_seg_cm as f64 {
                return;
            }
        }
        let seg = SegRecord {
            sample1,
            sample2,
            chrname,
            start_pos: pos[t_seg_start] - gw_chr_start_bp,
            end_pos: pos[nsites - 1] - gw_chr_start_bp,
            ibd,
            n_snp: nsites - 1 - t_seg_start + 1,
        };
        out.add_seg(seg);
    }

    pub fn tabulate_sites_by_state_for_viterbi_traj(
        &self,
        chrid: usize,
        rs: &mut RunningStats,
        cv: &mut PerChrModelVariables,
    ) {
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);

        for (t, _) in (start_chr..end_chr).enumerate() {
            match cv.traj[t] == 0 {
                true => rs.count_ibd_vit += 1,
                false => rs.count_dbd_vit += 1,
            }
        }
    }
    pub fn update_model(&self, rs: &RunningStats, ms: &mut ModelParamState) {
        let args = &self.data.args;
        let pi = &mut ms.model.pi;
        if (rs.count_ibd_fb + rs.count_dbd_fb) == 0.0 {
            eprintln!("Insufficient info to estimate parameters");
            eprintln!("Do you have only one variant per chromosome?");
        }

        // reestimate parameters
        //
        // Part of Schaffner's Additional File 1: model update via Baum-Welch;
        // trans_obs/pred: similar to Rabiner's 40b and 40c ;
        // count_ibd_fb/count_ibd_fb: Similar to Rabiner's Eq 39a and 40a, also
        // related to Eq 43a );
        pi[0] = rs.count_ibd_fb / (rs.count_ibd_fb + rs.count_dbd_fb);
        ms.model.k_rec *= rs.trans_obs / rs.trans_pred;

        // cap k_rec
        if ms.model.k_rec > self.data.args.k_rec_max {
            ms.model.k_rec = self.data.args.k_rec_max;
        }

        // prevent model being trapped
        if (ms.iiter < args.max_iter as usize) && (!ms.finish_fit) {
            if pi[0] < 1e-5 {
                pi[0] = 1e-5;
            } else if pi[0] > 1.0 - 1e-5 {
                pi[0] = 1.0 - 1e-5;
            } else {
            }
            if ms.model.k_rec < 1e-5 {
                ms.model.k_rec = 1e-5;
            }
        }
        pi[1] = 1.0 - pi[0];

        // calculate delta
        let delpi = pi[0] - ms.last_model.pi[0];
        let delk = ms.model.k_rec - ms.last_model.k_rec;
        let delrelk = delk / ms.model.k_rec;

        // update last model
        ms.last_model = ms.model;

        // check fit
        if (delpi.abs() < args.fit_thresh_dpi)
            && ((delk.abs() < args.fit_thresh_dk) || (delrelk < args.fit_thresh_drelk))
        {
            ms.finish_fit = true;
        }
    }
    pub fn print_hmm_frac(
        &self,
        ms: &ModelParamState,
        pair: (usize, usize),
        rs: &RunningStats,
        out: &mut OutputBuffer<'a>,
    ) {
        let sample1 = &self.data.samples.v()[pair.0];
        let sample2 = &self.data.samples.v()[pair.1];

        let sum = ms.num_info_site;
        let discord = ms.num_discord as f64 / ms.num_info_site as f64;
        let max_phi = ms.model.max_phi;
        let iter = ms.iiter;
        let k_rec = ms.model.k_rec;
        let ntrans = rs.ntrans;
        let seq_ibd_ratio = rs.seq_ibd / (rs.seq_ibd + rs.seq_dbd);
        let count_ibd_fb_ratio =
            rs.count_ibd_fb as f64 / (rs.count_dbd_fb + rs.count_ibd_fb) as f64;
        let count_ibd_vit_ratio =
            rs.count_ibd_vit as f64 / (rs.count_dbd_vit + rs.count_ibd_vit) as f64;
        let frac = FracRecord {
            sample1,
            sample2,
            sum,
            discord,
            max_phi,
            iter,
            k_rec,
            ntrans,
            seq_ibd_ratio,
            count_ibd_fb_ratio,
            count_ibd_vit_ratio,
        };
        out.add_frac(frac);
    }
}

pub struct RunningStats {
    /// Running sum of probility of observing transition (over the genome)
    /// See Shanffner 2018, Additional File 1, page 5.
    pub trans_obs: f64,
    /// Running sum of predicted probility of transition (over the genome)
    pub trans_pred: f64,
    /// Counter of transisions
    pub ntrans: usize,
    /// Running sum of length of segments in BP that are IBD
    pub seq_ibd: f64,
    /// Running sum of length of segments in BP that are not-IBD
    pub seq_dbd: f64,
    /// Running sum of probability of site IBD given the model and oberservations
    /// (\sum_t gamma_t(i))
    pub count_ibd_fb: f64,
    /// Running sum of probability of site not-IBD given the model and oberservations
    /// (\sum_t gamma_t(i))
    pub count_dbd_fb: f64,
    /// number of sites IBD from the best path (Viterbi)
    pub seq_ibd_fb: f64,
    pub seq_dbd_fb: f64,
    pub count_ibd_vit: usize,
    /// number of sites not-IBD from the best path (Viterbi)
    pub count_dbd_vit: usize,
}

impl RunningStats {
    pub fn new() -> Self {
        Self {
            trans_obs: 0.0,
            trans_pred: 0.0,
            ntrans: 0,
            seq_ibd: 0.0,
            seq_dbd: 0.0,
            count_ibd_fb: 0.0,
            count_dbd_fb: 0.0,
            seq_ibd_fb: 0.0,
            seq_dbd_fb: 0.0,
            count_ibd_vit: 0,
            count_dbd_vit: 0,
        }
    }
}

#[test]
fn test_hmm() {
    use crate::args::Arguments;
    use crate::data::OutputFiles;
    let args = Arguments::new_for_test();

    let input = InputData::from_args(&args);

    let runner = HmmRunner::new(&input);

    {
        let out = OutputFiles::new_from_args(&args, None, None);
        let mut a = vec![];
        let mut b = vec![];
        for pair in input.pairs.iter() {
            let pair = (pair.0 as usize, pair.1 as usize);
            let mut out = OutputBuffer::new(&out, 1, 1);
            runner.run_hmm_on_pair(pair, &mut a, &mut b, &mut out, false);
            out.flush_frac();
            out.flush_segs();
        }
    }

    // run hmmibd
    use std::process::{Command, Stdio};
    if !std::path::Path::new("hmmIBD").exists() {
        if !std::path::Path::new("hmmIBD.c").exists() {
            eprintln!("./hmmIBD or ./hmmIBD.c can be found in current folder");
            std::process::exit(-1);
        }
        Command::new("gcc")
            .args(["hmmIBD.c", "-lm", "-O2", "-o", "hmmIBD"])
            .status()
            .unwrap();
    }
    Command::new("./hmmIBD")
        .args([
            "-i",
            "samp_data/pf3k_Cambodia_13.txt",
            "-I",
            "samp_data/pf3k_Ghana_13.txt",
            "-f",
            "samp_data/freqs_pf3k_Cambodia_13.txt",
            "-F",
            "samp_data/freqs_pf3k_Ghana_13.txt",
            "-o",
            "tmp_hmmibd",
        ])
        .stdout(Stdio::null())
        .status()
        .unwrap();

    // load data from hmmibd-rs results
    {
        use std::collections::HashMap;
        #[derive(Default, Eq, PartialEq, PartialOrd, Ord, Debug)]
        struct Seg {
            sample1: u32,
            sample2: u32,
            chr: u32,
            start: u32,
            end: u32,
            different: u32,
            nsnp: u32,
        }
        fn load_data_from_hmm_seg_res(
            path: &str,
            sample_name_map: &HashMap<String, u32>,
        ) -> Vec<Seg> {
            let mut res = vec![];
            for line in std::fs::read_to_string(path)
                .unwrap()
                .trim()
                .split("\n")
                .skip(1)
            {
                let mut seg = Seg::default();
                for (ifield, field) in line.split("\t").enumerate() {
                    match ifield {
                        0 => seg.sample1 = sample_name_map[field],
                        1 => seg.sample2 = sample_name_map[field],
                        2 => seg.chr = field.parse::<u32>().unwrap(),
                        3 => seg.start = field.parse::<u32>().unwrap(),
                        4 => seg.end = field.parse::<u32>().unwrap(),
                        5 => seg.different = field.parse::<u32>().unwrap(),
                        6 => seg.nsnp = field.parse::<u32>().unwrap(),
                        _ => {}
                    }
                }
                res.push(seg);
                res.sort();
            }
            res
        }

        #[derive(Default, PartialEq, PartialOrd, Debug)]
        struct Frac {
            pub sample1: u32,
            pub sample2: u32,
            pub num_info_sites: u32,
            pub discord: f64,
            pub max_phi: f64,
            pub iter: u32,
            pub k_rec: f64,
            pub ntrans: u32,
            pub seq_ibd_ratio: f64,
            pub count_ibd_fb_ratio: f64,
            pub count_ibd_vit_ratio: f64,
        }

        fn load_data_from_hmm_frac_res(
            path: &str,
            sample_name_map: &HashMap<String, u32>,
        ) -> Vec<Frac> {
            let mut res = vec![];
            for line in std::fs::read_to_string(path)
                .unwrap()
                .trim()
                .split("\n")
                .skip(1)
            {
                let mut frac = Frac::default();
                for (ifield, field) in line.split("\t").enumerate() {
                    match ifield {
                        0 => frac.sample1 = sample_name_map[field],
                        1 => frac.sample2 = sample_name_map[field],
                        2 => frac.num_info_sites = field.parse::<u32>().unwrap(),
                        3 => frac.discord = field.parse::<f64>().unwrap(),
                        4 => frac.max_phi = field.parse::<f64>().unwrap(),
                        5 => frac.iter = field.parse::<u32>().unwrap(),
                        6 => frac.k_rec = field.parse::<f64>().unwrap(),
                        7 => frac.ntrans = field.parse::<u32>().unwrap(),
                        8 => frac.seq_ibd_ratio = field.parse::<f64>().unwrap(),
                        9 => frac.count_ibd_fb_ratio = field.parse::<f64>().unwrap(),
                        10 => frac.count_ibd_vit_ratio = field.parse::<f64>().unwrap(),
                        _ => {}
                    }
                }
                res.push(frac);
                res.sort_by(|a, b| a.partial_cmp(&b).unwrap());
            }
            res
        }

        let sample_name_map = input.samples.m();

        assert_eq!(
            load_data_from_hmm_seg_res("tmp_hmmibdrs.hmm.txt", sample_name_map),
            load_data_from_hmm_seg_res("tmp_hmmibd.hmm.txt", sample_name_map)
        );
        assert_eq!(
            load_data_from_hmm_frac_res("tmp_hmmibdrs.hmm.txt", sample_name_map),
            load_data_from_hmm_frac_res("tmp_hmmibd.hmm.txt", sample_name_map)
        );
    }
}

#[test]
fn test_hmm_with_bcf() {
    use crate::args::Arguments;
    use crate::data::OutputFiles;
    let args = Arguments::new_for_test_bcf();

    let out = OutputFiles::new_from_args(&args, None, None);
    let input = InputData::from_args(&args);

    let runner = HmmRunner::new(&input);

    for pair in input.pairs.iter().take(2) {
        let pair = (pair.0 as usize, pair.1 as usize);
        let mut out = OutputBuffer::new(&out, 1, 1);
        runner.run_hmm_on_pair(pair, &mut out, false);
    }
}
