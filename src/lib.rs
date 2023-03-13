extern crate nalgebra as na;
use na::{DMatrix, DVector};
use rand::{Rng, thread_rng};

pub fn createA_2dgeom(sig: fn((f64,f64)) -> f64 , d: fn((f64,f64)) -> f64 , ps: &DMatrix<(f64,f64)>) -> DMatrix<f64>{
    eprintln!("Creating 2d A matrix");
    let cols = ps.ncols();
    let rows = ps.nrows();
    let n = rows*cols;
    let mut a = DMatrix::from_element(n, n, 0.0 as f64);
    let del: f64=1.0;

    for (i,mut r) in a.row_iter_mut().enumerate(){
        if i>0 {
            r[i-1]+=-d(ps[i])/del.powf(2.0);
        }
        if i%cols<cols-1 {
            r[i+1]+=-d(ps[i])/del.powf(2.0);
        }
        r[i]+=(2.0*d(ps[i])/del.powf(2.0)) + sig(ps[i]);

        if i>cols {
            r[i-cols]+=-d(ps[i])/del.powf(2.0);
        }
        if i<n-cols-1 {
            r[i+cols]+=-d(ps[i])/del.powf(2.0);
        }
        r[i]+=(2.0*d(ps[i])/del.powf(2.0)) + sig(ps[i]);
    }
    a
}


pub fn createA_geom(sig: fn(f64) -> f64 , d: fn(f64) -> f64 , xs: &Vec<f64>) -> DMatrix<f64>{
    let n = xs.len();
    let mut a = DMatrix::from_element(n, n, 0.0 as f64);
    let diff: Vec<f64>= xs.windows(2).map(|e| e[1] - e[0]).collect();
    //let del = a/(n as f64);
    for (i,mut r) in a.row_iter_mut().enumerate(){
        if i>0 {
            r[i-1]=-d(xs[i])/diff[i%(n-1)].powi(2);
        }
        if i<n-1 {
            r[i+1]=-d(xs[i])/diff[i%(n-1)].powi(2);
        }
        r[i]=(2.0*d(xs[i])/diff[i%(n-1)].powi(2)) + sig(xs[i]);
    }
    a
}

pub fn sol_geom(sig: fn(f64) -> f64 , d: fn(f64) -> f64, vsig: fn(f64) -> f64 , xs: &Vec<f64>) -> (f64, DVector<f64>) {
    let a = createA_geom(sig, d, xs);
    let b = DMatrix::from_diagonal(&DVector::from_iterator(xs.len(), xs.clone().into_iter().map(vsig)));
    //let ainv = a.try_inverse().unwrap();
    //let qr = QR::new(ainv * b);
    eigen(a, b)
}

pub fn sol_2dgeom(sig: fn((f64,f64)) -> f64 , d: fn((f64,f64)) -> f64, vsig: fn(&(f64,f64)) -> f64, ps: &DMatrix<(f64,f64)>) -> (f64, DMatrix<f64>) {
    eprintln!("Solving 2d geometry");
    let a = createA_2dgeom(sig, d, ps);
    let b = DMatrix::from_diagonal(&DVector::from_iterator(ps.len(), ps.clone().into_iter().map(vsig)));
    let (k,p) = inv_pow(a, b);
    let mut phi: DMatrix<f64> = DMatrix::zeros(ps.nrows(), ps.ncols());
    for i in 0..phi.nrows()*phi.ncols(){
        phi[i] = p[i];
    }
    (k,phi)
}

fn eigen(a: DMatrix<f64>, b: DMatrix<f64>) -> (f64, DVector<f64>) {
    let M = a.try_inverse().unwrap() * b;
    let eigen = M.clone().schur().eigenvalues().unwrap();
    let k = eigen[0];
    let A = M.clone() - k * DMatrix::identity(M.nrows(), M.ncols());
    let decomp = A.lu();
    let zero: DVector<f64> = DVector::zeros(M.nrows());
    let phi = decomp.solve(&zero).unwrap();
    (k,phi)
    
}

fn inv_pow(a: DMatrix<f64>, b: DMatrix<f64>) -> (f64, DVector<f64>) {
    eprintln!("Running inv pow algo");
    let mut rng = thread_rng();
    let mut phi: DVector<f64> = DVector::from_vec((0..b.ncols()).map(|_| rng.gen()).collect());
    eprintln!("Calculating inverse!");
    let c = a.try_inverse().unwrap() * b;
    eprintln!("Done");
    let mut converged = false;
    let mut kold: f64 =0.0;
    let mut k: f64 = 0.0;
    while !converged {
        phi = c.clone()*phi;
        k = *phi.iter().max_by(|a,b| a.abs().total_cmp(&b.abs())).unwrap();
        let phim = phi.iter().max_by(|a,b| a.total_cmp(&b)).unwrap();
        phi = phi.map(|x| x/phim);
        if (kold-k).abs() < 1e-7 {
            converged = true;
        }
        kold=k;
    }
    (k, phi)
}
