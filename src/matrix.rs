use std::ops::Index;
use std::ops::IndexMut;

pub trait AsOption: Sized {
    fn as_option(&self) -> Option<Self>;
    fn is_none(&self) -> bool;
    fn is_some(&self) -> bool;
    fn none() -> Self;
}
impl AsOption for u8 {
    fn as_option(&self) -> Option<Self> {
        match *self {
            u8::MAX => None,
            _ => Some(*self),
        }
    }

    fn is_none(&self) -> bool {
        u8::MAX == *self
    }
    fn is_some(&self) -> bool {
        u8::MAX != *self
    }

    fn none() -> Self {
        u8::MAX
    }
}
impl AsOption for u32 {
    fn as_option(&self) -> Option<Self> {
        match *self {
            u32::MAX => None,
            _ => Some(*self),
        }
    }
    fn is_none(&self) -> bool {
        u32::MAX == *self
    }
    fn is_some(&self) -> bool {
        u32::MAX != *self
    }
    fn none() -> Self {
        u32::MAX
    }
}
impl AsOption for f32 {
    fn as_option(&self) -> Option<Self> {
        match self.is_nan() {
            true => None,
            _ => Some(*self),
        }
    }
    fn is_none(&self) -> bool {
        self.is_nan()
    }
    fn is_some(&self) -> bool {
        !self.is_nan()
    }
    fn none() -> Self {
        f32::NAN
    }
}
impl AsOption for f64 {
    fn as_option(&self) -> Option<Self> {
        match self.is_nan() {
            true => None,
            _ => Some(*self),
        }
    }
    fn is_none(&self) -> bool {
        self.is_nan()
    }
    fn is_some(&self) -> bool {
        !self.is_nan()
    }
    fn none() -> Self {
        f64::NAN
    }
}

pub struct Matrix<T>
where
    T: Copy + Default + AsOption + PartialOrd,
{
    data: Vec<T>,
    ncols: usize,
    nrows: usize,
}

impl<T> Matrix<T>
where
    T: Copy + Default + AsOption + PartialOrd,
{
    pub fn from_shape(nrows: usize, ncols: usize, init_val: T) -> Self {
        let data = vec![init_val; nrows * ncols];
        Self { data, ncols, nrows }
    }

    pub fn resize_and_clear(&mut self, nrows: usize, ncols: usize, init_val: T) {
        self.data.clear();
        self.data.resize(nrows * ncols, init_val);
        self.ncols = ncols;
        self.nrows = nrows;
    }

    pub fn transpose(&mut self) {
        let mut data2 = Vec::<T>::with_capacity(self.data.len());
        for col in 0..self.ncols {
            for row in 0..self.nrows {
                data2.push(self.get_at(row, col));
            }
        }
        *self = Self {
            data: data2,
            ncols: self.nrows,
            nrows: self.ncols,
        };
    }

    pub fn get_nrows(&self) -> usize {
        self.nrows
    }
    pub fn get_ncols(&self) -> usize {
        self.ncols
    }

    pub fn get_row_raw_slice(&self, row: usize) -> &[T] {
        &self.data[(row * self.ncols)..((row + 1) * self.ncols)]
    }

    pub fn get_row_iter(&self, row: usize) -> impl Iterator<Item = Option<T>> + '_ {
        self.data[(row * self.ncols)..((row + 1) * self.ncols)]
            .into_iter()
            .map(|x| x.as_option())
    }

    pub fn get_at(&self, row: usize, col: usize) -> T {
        self.data[self.ncols * row + col]
    }

    pub fn set_at(&mut self, row: usize, col: usize, val: T) {
        self.data[self.ncols * row + col] = val;
    }
    pub fn merge(&mut self, other: &Self) {
        assert_eq!(self.ncols, other.ncols);
        self.data.extend(other.data.iter());
        self.nrows += other.nrows;
    }
}

impl<T> Index<usize> for Matrix<T>
where
    T: Copy + Default + AsOption + PartialOrd,
{
    type Output = [T];
    fn index(&self, index: usize) -> &Self::Output {
        let s = index * self.ncols;
        let e = s + self.ncols;
        &self.data[s..e]
    }
}

impl<T> IndexMut<usize> for Matrix<T>
where
    T: Copy + Default + AsOption + PartialOrd,
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        let s = index * self.ncols;
        let e = s + self.ncols;
        &mut self.data[s..e]
    }
}
pub struct MatrixBuilder<T: Copy + Default + AsOption + PartialOrd> {
    data: Vec<T>,
    ncols: usize,
}
impl<T> MatrixBuilder<T>
where
    T: Copy + Default + AsOption + PartialOrd,
{
    pub fn new(ncols: usize) -> Self {
        Self {
            data: vec![],
            ncols,
        }
    }

    pub fn push(&mut self, val: Option<T>) {
        match val {
            Some(val) => self.data.push(val),
            None => self.data.push(T::none()),
        }
    }

    pub fn finish(&mut self) -> Matrix<T> {
        assert_eq!(self.data.len() % self.ncols, 0);
        let nrows = self.data.len() / self.ncols;
        Matrix {
            data: std::mem::take(&mut self.data),
            ncols: self.ncols,
            nrows,
        }
    }
}
