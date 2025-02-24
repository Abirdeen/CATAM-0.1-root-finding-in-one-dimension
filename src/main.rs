#![allow(dead_code)]

extern crate num;
use num::signum;
use num::traits::pow;

extern crate composing;
use composing::compose_fn;

mod root_search;
mod functional;
mod test_function;
mod test_func_derivative;

type ContinuousFunction<'a> = dyn Fn(f64) -> f64 + 'a;

fn main() -> Result<(), Box<&'static str>> {

    let initial_func: &ContinuousFunction = &(test_function::trig as fn(f64) -> f64);
    let deriv: &ContinuousFunction = &(test_func_derivative::trig as fn(f64) -> f64);
    let initial_val: f64 = -4.8;
    let trunc_err: f64 = 1.0/pow(10.0, 5);
    let max_iter: usize = 100;

    let transformed_func: Box<ContinuousFunction> = functional::x_minus(functional::newton_raphson(&initial_func, deriv));
    let output: Result<(f64, Vec<f64>), Vec<f64>> = root_search::fixed_point(&transformed_func, initial_val, trunc_err, max_iter);
    match output {
        Ok(v) => {
            let (mut res, seq) = v;
            res = (res/trunc_err).round()*trunc_err;
            println!("Root is at {} ± {}", res, trunc_err);
            println!("{}", seq.len());
        }
        Err(_e) => {println!("Sequence did not converge in {} iterations", max_iter)}
    } 

    // let round: fn(f64) -> f64 = |x| (x * pow(10.0, 6)).round()/pow(10.0, 6);
    // let mut prev_diff: f64 = 1.0;
    // for (index, value) in seq.iter().enumerate() {
    //     let diff: f64 = *value-4.0;
    //     if index == 0 {
    //         println!("| {} | {} | {} |||", index, round(*value), round(diff));
    //         prev_diff = diff;
    //         continue
    //     };
    //     let ratio: f64 = diff/prev_diff;
    //     let approx_one: f64 = diff*(index as f64)*7.0/40.0;
    //     if index < 6 || index > 730 {
    //         println!("| {} | {} | {} | {} | {} |", index, round(*value), round(diff), round(ratio), round(approx_one));
    //     }
    // }
    Ok(())
}