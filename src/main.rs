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

    Ok(())
}