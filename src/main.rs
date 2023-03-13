#![feature(iter_intersperse)]
#![allow(non_snake_case)]
#![allow(non_snake_case)]

extern crate nalgebra as na;
use na::DMatrix;
use neutron::*;

fn main(){
    //Input data
    let a=100.0;
    let n=100;

    const W: f64=1.0;

    let sigma=|p: (f64,f64)| -> f64 {if p.0.powi(2)+p.1.powi(2)<W.powi(2){0.9}else{0.1532}};
    let vsigma=|p: &(f64,f64)| -> f64 {0.157};
    let d=|p: (f64,f64)| -> f64 {3.85};

    let xs: Vec<f64> = (0..n).map(|i| { (a/(n as f64))*(i as f64) - a/2.0}).collect();
    let ys: Vec<f64> = (0..n).map(|i| { (a/(n as f64))*(i as f64) - a/2.0}).collect();

    let ps: DMatrix<(f64,f64)> =DMatrix::from_iterator(n,n,(0..n*n).map(|i| {(xs[i%n],ys[i/n])}));
    let (k,phi) = sol_2dgeom(sigma, d, vsigma, &ps);

    for i in 0..n*n{
        let p = ps[i];
        println!("{}\t{}\t{}",p.0,p.1,phi[i])
    }
    //dbg!(k);
    ////let A = createA_geom(sigma, d, &x);

    ////let s=|x: &f64| -> f64 {x*0.1};
    ////let mut S = DVector::from(x.clone().iter().map(s).collect::<Vec<f64>>());
    ////let decomp = A.clone().lu();
    ////let x = decomp.solve(&S).expect("Linear resolution failed.");
    //for (xval,pval) in zip(x,phi.iter()) {
    //    println!("{}\t{}",xval,pval);
    //}
}
