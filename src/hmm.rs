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
        out: &mut OutputBuffer<'a>,
        suppress_frac: bool,
    ) {
        let mut ms = ModelParamState::new(self.data);
        let mut cv = PerChrModelVariables::new();
        let mut rs = RunningStats::new();

        for iiter in 0..self.data.args.max_iter {
            ms.iiter = iiter as usize;
            ms.model.max_phi = 0.0;
            rs = RunningStats::new();

            self.run_over_fit_iterations(pair, &mut ms, &mut rs, &mut cv, out);
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
        out: &mut OutputBuffer<'a>,
    ) {
        let nchrom = self.data.genome.get_nchrom();
        for chrid in 0..nchrom as usize {
            self.run_over_single_chromsome(pair, chrid, ms, rs, cv, out);
            // if (chrid > 0) || (ms.iiter > 0) {
            //     continue;
            // }
            //     println!(
            //         "iter: {}, chrid {}: \ntraj: {:?}, \nalpha: {:?}, \nbeta: {:?}, \nb={:?} \nphi={:?} \npsi={:?}",
            //         ms.iiter,
            //         chrid,
            //         &cv.traj[0..4],
            //         &cv.alpha.get_row_raw_slice(0)[0..6],
            //         &cv.beta.get_row_raw_slice(0)[0..6],
            //         &cv.b.get_row_raw_slice(0)[0..6],
            //         &cv.phi.get_row_raw_slice(0)[0..6],
            //         &cv.psi.get_row_raw_slice(0)[0..6],
            //     );
        }
    }

    pub fn run_over_single_chromsome(
        &self,
        pair: (usize, usize),
        chrid: usize,
        ms: &mut ModelParamState,
        rs: &mut RunningStats,
        cv: &mut PerChrModelVariables,
        out: &mut OutputBuffer<'a>,
    ) {
        // skip empty chromosomes
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);
        if start_chr == end_chr {
            return;
        }
        cv.resize_and_clear(self.data, chrid);
        // iterator over snp forward to calculate
        // - b
        // - alpha
        // - delta
        // - psi
        // - max_phiL
        self.run_over_snp_1_fwd(pair, chrid, ms, cv);

        // iterator over snp backward to calcuate
        // - beta
        self.run_over_snp_2_back(chrid, ms, cv);

        // - traj
        self.run_over_snp_3_back(chrid, cv);

        // iterator over snp forward again to calcualte
        // - xi
        // - gamma
        // - update rs
        self.run_over_snp_4_fwd(chrid, ms, rs, cv);

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
        // rs: &mut RunningStats,
        cv: &mut PerChrModelVariables,
    ) {
        let cm = self.data.sites.get_pos_cm_slice();
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);
        let geno = &self.data.geno;
        let args = &self.data.args;
        let freq1 = &self.data.freq1;
        let freq2 = match self.data.freq2.as_ref() {
            Some(freq2) => freq2,
            None => freq1,
        };
        let pi = &ms.model.pi;
        let nall = &self.data.nall;

        for isnp in start_chr..end_chr {
            let snp_ind = (isnp - start_chr) as usize;
            let isnp = isnp as usize;
            let gi = geno[pair.0][isnp];
            let gj = geno[pair.1][isnp];
            let eps = args.eps;
            let mut sv = PerSnpModelVariables::new();

            // get num of informative sites
            if (ms.iiter == 0) && (gi.is_some()) && (gj.is_some()) {
                // println!("self.data.majall={:?}", &self.data.majall);
                if (gi == gj) && (gj == self.data.majall[isnp]) {
                } else {
                    ms.num_info_site += 1;
                }
                if gi != gj {
                    ms.num_discord += 1;
                }
            }

            // update b for a given sites (full loop will update whole chromosome)
            let pright = 1.0 - eps * (nall[isnp] - 1) as f64;
            let (b0t, b1t) = if (gi == u8::MAX) || (gj == u8::MAX) {
                // missing
                (1.0, 1.0)
            } else {
                let f1i = freq1[isnp][gi as usize].as_option().unwrap() as f64;
                let f1j = freq1[isnp][gj as usize].as_option().unwrap() as f64;
                let f2i = freq2[isnp][gi as usize].as_option().unwrap() as f64;
                let f2j = freq2[isnp][gj as usize].as_option().unwrap() as f64;
                if gi == gj {
                    // concordant genotype
                    cv.nsites += 1;
                    let fmean = ((f1i + f2j) / 2.0) as f64;
                    let b0t = pright * pright * fmean + eps * eps * (1.0 - fmean);
                    let b1t = pright * pright * f1i * f2j
                        + pright * eps * f1i * (1.0 - f2j)
                        + pright * eps * (1.0 - f1i) * f2j
                        + eps * eps * (1.0 - f1i) * (1.0 - f2j);
                    (b0t, b1t)
                } else {
                    // discordant genotype
                    cv.nsites += 1;
                    let fmeani = (f1i + f2i) / 2.0;
                    let fmeanj = (f1j + f2j) / 2.0;
                    let b0t =
                        pright * eps * (fmeani + fmeanj) + eps * eps * (1.0 - fmeani - fmeanj);
                    let b1t = pright * pright * f1i * f2j
                        + pright * eps * (f1i * (1.0 - f2j) + f2j * (1.0 - f1i))
                        + eps * eps * (1.0 - f1i) * (1.0 - f2j);
                    (b0t, b1t)
                }
            };
            cv.b[0][snp_ind] = b0t;
            cv.b[1][snp_ind] = b1t;

            // println!("snp_ind={snp_ind}\tpos={}\tb[0]={b0t:.5}\tb[1]={b1t:.5}\tgi={gi}\tgj={gj}\tpright={pright}\t{}", pos[isnp],nall[snp_ind]);

            // println!("++++++++++++++++++ snp_ind={snp_ind}");

            // calcuate alpha, delta(phi) and psi
            if snp_ind == 0 {
                // initiation
                cv.psi[0][snp_ind] = 0;
                cv.psi[1][snp_ind] = 0;
                cv.phi[0][snp_ind] = pi[0].ln() + b0t.ln();
                cv.phi[1][snp_ind] = pi[1].ln() + b1t.ln();
                cv.alpha[0][snp_ind] = pi[0] * b0t;
                cv.alpha[1][snp_ind] = pi[1] * b1t;
                // assert!(cv.alpha[0][snp_ind].gt(&0.0));
                // assert!(cv.alpha[1][snp_ind].gt(&0.0));
                // TODO: scalre value is very confusing
                cv.scale[0] = 1.0;
            } else {
                // induction
                let ptrans = ms.model.k_rec * (cm[isnp] - cm[isnp - 1]) as f64 / 100.0;
                sv.a[0][1] = 1.0 - pi[0] - (1.0 - pi[0]) * (-ptrans).exp();
                sv.a[1][0] = 1.0 - pi[1] - (1.0 - pi[1]) * (-ptrans).exp();
                sv.a[0][0] = 1.0 - sv.a[0][1];
                sv.a[1][1] = 1.0 - sv.a[1][0];

                for js in 0..2 {
                    let mut max_val = f64::MIN;
                    let mut alpha_jt = 0.0;

                    // TODO: scalre value is very confusing, should move it out of for js in 0..2 loop
                    cv.scale[snp_ind] = 0.0;
                    for is in 0..2 {
                        let score = cv.phi[is][snp_ind - 1] + sv.a[is][js].ln();
                        if score > max_val {
                            max_val = score;
                            cv.psi[js][snp_ind] = is as u8;
                        }
                        cv.phi[js][snp_ind] = max_val + cv.b[js][snp_ind].ln();
                        alpha_jt += cv.alpha[is][snp_ind - 1] * sv.a[is][js] * cv.b[js][snp_ind];
                        // println!(
                        //     "chrid={chrid}, isnp={isnp} snp_ind={snp_ind} js={js} is={is} alpha={} a={} b={} ",
                        //     cv.alpha[is][snp_ind - 1],
                        //     sv.a[is][js],
                        //     cv.b[js][snp_ind],
                        // );
                    }
                    cv.alpha[js][snp_ind] = alpha_jt;

                    // TODO: scalre value is very confusing
                    cv.scale[snp_ind] += cv.alpha[js][snp_ind];
                }
                // scaling
                // println!("scale: {} ", cv.scale[snp_ind]);
                for js in 0..2 {
                    // assert!(cv.scale[snp_ind] > 0.0);
                    cv.alpha[js][snp_ind] = cv.alpha[js][snp_ind] / cv.scale[snp_ind];
                }
            }

            // println!("snp_ind={snp_ind}\tpos={}\tpsi0={:.4}\tpsi1={:.4}\tphi0={:.4}\tphi1={:.4}\talpha0={:.4}\talpha1={:.4}",
            //     pos[isnp],
            //     cv.psi[0][snp_ind],
            //     cv.psi[1][snp_ind],
            //     cv.phi[0][snp_ind],
            //     cv.phi[1][snp_ind],
            //     cv.alpha[0][snp_ind],
            //     cv.alpha[1][snp_ind]);
        }
        // let last_snp = (end_chr - 1) as usize;
        let last_snp_indx = (end_chr - start_chr - 1) as usize;
        let mut max_phi_local = cv.phi[1][last_snp_indx];
        let mut max = 1;
        if cv.phi[1][last_snp_indx] < cv.phi[0][last_snp_indx] {
            max_phi_local = cv.phi[0][last_snp_indx];
            max = 0;
        }
        cv.traj[last_snp_indx] = max;
        ms.model.max_phi += max_phi_local;

        // println!(
        //     "last_snp={last_snp_indx}\tmax={max}\tmax_phi{}",
        //     ms.model.max_phi
        // );
    }
    pub fn run_over_snp_2_back(
        &self,
        chrid: usize,
        ms: &mut ModelParamState,
        cv: &mut PerChrModelVariables,
    ) {
        let cm = self.data.sites.get_pos_cm_slice();
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);
        let pi = &ms.model.pi;

        // Cacluating Beta
        ///////////////////////////////////////

        // init beta
        let last = end_chr - start_chr - 1;
        cv.beta[0][last] = 1.0;
        cv.beta[1][last] = 1.0;
        // induction
        for isnp in (start_chr..end_chr).rev().skip(1) {
            let snp_ind = isnp - start_chr;
            let mut sv = PerSnpModelVariables::new();

            // let ptrans = ms.model.k_rec * args.rec_rate * (pos[isnp + 1] - pos[isnp]) as f64;
            let ptrans = ms.model.k_rec * (cm[isnp + 1] - cm[isnp]) as f64 / 100.0;
            sv.a[0][1] = 1.0 - pi[0] - (1.0 - pi[0]) * (-ptrans).exp();
            sv.a[1][0] = 1.0 - pi[1] - (1.0 - pi[1]) * (-ptrans).exp();
            sv.a[0][0] = 1.0 - sv.a[0][1];
            sv.a[1][1] = 1.0 - sv.a[1][0];

            let mut sum = [0.0, 0.0];
            for (is, js) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
                sum[is] += cv.beta[js][snp_ind + 1] * sv.a[is][js] * cv.b[js][snp_ind + 1];
            }
            cv.beta[0][snp_ind] = sum[0] / cv.scale[snp_ind];
            cv.beta[1][snp_ind] = sum[1] / cv.scale[snp_ind];

            // println!(
            //     "snp_ind={snp_ind}\tpos={}\tbeta0={:.4}\tbeta1={:.4}",
            //     pos[isnp], cv.beta[0][snp_ind], cv.beta[1][snp_ind]
            // )
        }
    }
    pub fn run_over_snp_3_back(&self, chrid: usize, cv: &mut PerChrModelVariables) {
        // backtracking
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);
        // let pos = self.data.sites.get_pos_slice();
        // let gw_chr_start = self.data.genome.get_gwchrstarts()[chrid];

        let last_snp = end_chr - start_chr - 1;
        let mut max = cv.traj[last_snp] as usize;
        // println!("last_snp max: {}", max);
        for isnp in (start_chr..end_chr).rev().skip(1) {
            let snp_ind = isnp - start_chr;
            cv.traj[snp_ind] = cv.psi[max][snp_ind + 1];
            max = cv.traj[snp_ind] as usize;
            // println!(
            //     "snp_ind={snp_ind}\tpos={}\ttraj={}",
            //     pos[isnp] - gw_chr_start,
            //     cv.traj[snp_ind]
            // )
        }
    }

    pub fn run_over_snp_4_fwd(
        &self,
        chrid: usize,
        ms: &mut ModelParamState,
        rs: &mut RunningStats,
        cv: &mut PerChrModelVariables,
    ) {
        let pos = self.data.sites.get_pos_slice();
        let cm = self.data.sites.get_pos_cm_slice();
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);
        let pi = &ms.model.pi;
        let mut sv = PerSnpModelVariables::new();

        for (snp_ind, isnp) in (start_chr..end_chr).enumerate() {
            let mut p_ibd = cv.alpha[0][snp_ind] * cv.beta[0][snp_ind];
            let p_dbd = cv.alpha[1][snp_ind] * cv.beta[1][snp_ind];
            p_ibd = p_ibd / (p_ibd + p_dbd);

            rs.count_ibd_fb += p_ibd;
            rs.count_dbd_fb += 1.0 - p_ibd;

            // not the last snp for following code
            if snp_ind == end_chr - start_chr - 1 {
                break;
            }

            let delpos = (pos[isnp + 1] - pos[isnp]) as f64;
            // let ptrans = ms.model.k_rec * args.rec_rate * delpos;
            let ptrans = ms.model.k_rec * (cm[isnp + 1] - cm[isnp]) as f64 / 100.0;

            sv.a[0][1] = 1.0 - pi[0] - (1.0 - pi[0]) * (-ptrans).exp();
            sv.a[1][0] = 1.0 - pi[1] - (1.0 - pi[1]) * (-ptrans).exp();
            sv.a[0][0] = 1.0 - sv.a[0][1];
            sv.a[1][1] = 1.0 - sv.a[1][0];

            let mut xisum = 0.0;
            for (is, js) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
                sv.xi[is][js] = cv.alpha[is][snp_ind]
                    * sv.a[is][js]
                    * cv.b[js][snp_ind + 1]
                    * cv.beta[js][snp_ind + 1];
                xisum += sv.xi[is][js];
            }
            for (is, js) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
                sv.xi[is][js] /= xisum;
            }
            sv.gamma[0] = sv.xi[0][1] + sv.xi[0][0];
            sv.gamma[1] = sv.xi[1][1] + sv.xi[1][0];

            // println!("snp_ind={snp_ind}\tpos={}\tdelpos={delpos:.0}\txi={:.3},{:.3},{:.3},{:.3}\tgamma={:.3},{:.3}\ta01/10={:.6},{:.6}",
            //     pos[isnp] - self.data.genome.get_gwchrstarts()[chrid],
            //     sv.xi[0][0],
            //     sv.xi[0][1],
            //     sv.xi[1][0],
            //     sv.xi[1][1],
            //     sv.gamma[0],
            //     sv.gamma[1],
            //     sv.a[0][1],
            //     sv.a[1][0],
            // );

            rs.seq_ibd_fb += delpos * sv.xi[0][0];
            rs.seq_dbd_fb += delpos * sv.xi[1][1];
            rs.trans_obs += sv.xi[0][1] + sv.xi[1][0];
            rs.trans_pred += sv.gamma[0] * sv.a[0][1] + sv.gamma[1] * sv.a[1][0];
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
        if cv.nsites == 0 {
            return;
        }
        // skip if too old
        if let Some(max_tmrca) = self.data.args.filt_max_tmrca {
            let tmrca = ms.model.k_rec;
            if tmrca > max_tmrca {
                return;
            }
        }
        let pos = self.data.sites.get_pos_slice();
        let cm = self.data.sites.get_pos_cm_slice();
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);

        let chrname = self.data.genome.get_chrname(chrid);
        let sample1 = &self.data.samples.v()[pair.0];
        let sample2 = &self.data.samples.v()[pair.1];
        let gw_chr_starts = self.data.genome.get_gwchrstarts()[chrid];

        let mut start_snp = 0;
        let mut add_seq;
        //
        //      |-----|==========|-------|===========|
        //          start       -1 \isnp
        for (isnp, ipos) in (start_chr..end_chr).enumerate().skip(1) {
            if cv.traj[isnp] == cv.traj[isnp - 1] {
                continue;
            }

            rs.ntrans += 1;

            add_seq = (pos[ipos - 1] - pos[start_chr + start_snp] + 1) as f64;
            match cv.traj[isnp - 1] == 0 {
                true => rs.seq_ibd += add_seq,
                false => rs.seq_dbd += add_seq,
            };
            // print
            let start_pos = pos[start_chr + start_snp] - gw_chr_starts;
            let end_pos = pos[start_chr + isnp - 1] - gw_chr_starts;
            let ibd = cv.traj[isnp - 1];
            let n_snp = isnp - start_snp;
            let seg = SegRecord {
                sample1,
                sample2,
                chrname,
                start_pos,
                end_pos,
                ibd,
                n_snp,
            };

            if self.data.args.filt_ibd_only && (ibd == 1) {
                start_snp = isnp;
                continue;
            }

            // skip if too short
            if let Some(min_seg_cm) = self.data.args.filt_min_seg_cm {
                let cm = cm[start_chr + isnp - 1] - cm[start_chr + start_snp];
                if cm < min_seg_cm as f64 {
                    start_snp = isnp;
                    continue;
                }
            }
            out.add_seg(seg);
            // write!(
            //     &mut out.seg_file,
            //     "{sample1}\t{sample2}\t{chrname}\t{start_pos}\t{end_pos}\t{ibd}\t{n_snp}\n"
            // )
            // .unwrap();
            // // write!(&mut out.seg_file, "{:?}", cv.traj).unwrap();

            start_snp = isnp;
        }
        //
        let isnp = end_chr - start_chr - 1;
        let start_pos = pos[start_chr + start_snp] - gw_chr_starts;
        let end_pos = pos[end_chr - 1] - gw_chr_starts;
        add_seq = (end_pos - start_pos + 1) as f64;
        let ibd = cv.traj[isnp];
        let n_snp = isnp - start_snp + 1;
        match ibd == 0 {
            true => rs.seq_ibd += add_seq,
            false => rs.seq_dbd += add_seq,
        };
        if self.data.args.filt_ibd_only && (ibd == 1) {
            return;
        }

        // skip if too short
        if let Some(min_seg_cm) = self.data.args.filt_min_seg_cm {
            let cm = cm[end_chr - 1] - cm[start_chr + start_snp];
            if cm < min_seg_cm as f64 {
                return;
            }
        }
        let seg = SegRecord {
            sample1,
            sample2,
            chrname,
            start_pos,
            end_pos,
            ibd,
            n_snp,
        };
        out.add_seg(seg);
        // write!(
        //     &mut out.seg_file,
        //     "{sample1}\t{sample2}\t{chrname}\t{start_pos}\t{end_pos}\t{ibd}\t{n_snp}\n"
        // )
        // .unwrap();
    }

    pub fn tabulate_sites_by_state_for_viterbi_traj(
        &self,
        chrid: usize,
        rs: &mut RunningStats,
        cv: &mut PerChrModelVariables,
    ) {
        // let pos = self.data.sites.get_pos_slice();
        let (start_chr, end_chr) = self.data.sites.get_chrom_pos_idx_ranges(chrid);

        for (isnp, _) in (start_chr..end_chr).enumerate() {
            match cv.traj[isnp] == 0 {
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
        pi[0] = rs.count_ibd_fb / (rs.count_ibd_fb + rs.count_dbd_fb);
        ms.model.k_rec *= rs.trans_obs / rs.trans_pred;
        // println!(
        //     "trans_obs={:0.5}, trans_pred={:.5}",
        //     rs.trans_obs, rs.trans_pred,
        // );
        // assert!(!ms.model.k_rec.is_nan());
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
        // let delprob = ms.model.max_phi - ms.last_model.max_phi;
        let delrelk = delk / ms.model.k_rec;

        // update last model
        ms.last_model = ms.model;

        // check fit
        if (delpi.abs() < args.fit_thresh_dpi)
            && ((delk.abs() < args.fit_thresh_dk) || (delrelk < args.fit_thresh_drelk))
        {
            ms.finish_fit = true;
        }
        // println!(
        //     "iter={}\tk_rec={:.4}\tpi=[0]={:.4}\tmax_phi={:.4}\tfinish_fit={}\tcount_ibd_fb={:.3}\tcount_dbd_fb={:.3}\ttrans_obs={:.3}\ttrans_pred={:.3}",
        //     ms.iiter, ms.model.k_rec, ms.model.pi[0], ms.model.max_phi, ms.finish_fit, rs.count_ibd_fb, rs.count_dbd_fb, rs.trans_obs, rs.trans_pred
        // );
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

        // use std::io::Write;
        // write!(
        //     &mut out.frac_file,
        //     "{sample1}\t{sample2}\t{sum}\t{discord:.4}\t{max_phi:0.5e}\t{iter}\t{k_rec:.3}\t{ntrans}\t{seq_ibd_ratio:.5}\t{count_ibd_fb_ratio:.5}\t{count_ibd_vit_ratio:.5}\n"
        // )
        // .unwrap();
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

    let out = OutputFiles::new_from_args(&args, None, None);
    let input = InputData::from_args(&args);

    let runner = HmmRunner::new(&input);

    for pair in input.pairs.iter().take(2) {
        let pair = (pair.0 as usize, pair.1 as usize);
        let mut out = OutputBuffer::new(&out, 1, 1);
        runner.run_hmm_on_pair(pair, &mut out, false);
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

        assert_ne!(
            load_data_from_hmm_seg_res("tmp_hmmibdrs.hmm.txt", sample_name_map),
            load_data_from_hmm_seg_res("tmp_hmmibd.hmm.txt", sample_name_map)
        );
        assert_ne!(
            load_data_from_hmm_frac_res("tmp_hmmibdrs.hmm.txt", sample_name_map),
            load_data_from_hmm_frac_res("tmp_hmmibd.hmm.txt", sample_name_map)
        );
    }
}
