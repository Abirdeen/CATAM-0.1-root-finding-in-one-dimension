extern crate num;
use num::signum;
use num::traits::pow;

extern crate composing;
use composing::compose_fn;

type ContinuousFunction<'a> = dyn Fn(f64) -> f64 + 'a;

fn main() {
    let trunc_err: f64 = 0.00000001;
    let res:Result<f64, &str>  = root_search::binary(&test_function::trig, (-10.0, 15.0), trunc_err);

    // let initial_func: &ContinuousFunction = &(test_function::trig as fn(f64) -> f64);
    // let k: f64 = 0.0;
    // let initial_val: f64 = -2.0;
    // let trunc_err: f64 = 1.0/pow(10.0, 5);
    // let max_iter: usize = 10;

    // let transformed_func: Box<ContinuousFunction> = functional::x_minus(functional::frac(&initial_func, &k));
    // let (mut res, _seq) = root_search::fixed_point(&transformed_func, initial_val, trunc_err, max_iter);
    // res = (res/trunc_err).round()*trunc_err;
    // println!("Root is at {} ± {}", res, trunc_err);
}

#[allow(dead_code)]
mod root_search {
    use ContinuousFunction;
    use {compose_fn, signum};

    /// Find a root of a continuous function using binary search.
    ///
    /// Parameters
    /// ----------
    /// * `func` : A continuous function with a sign change over the given domain.
    /// * `domain` : The start and end points of the search interval. Note that `func` must be computable over the entire domain, including end points.
    /// * `trunc_err` : A float representing the acceptable truncation error for the search; e.g. `trunc_err=1` will result in finding the root +-1.
    ///
    /// Returns
    /// -------
    /// * `f64` : A floating point representing the root of the function.
    ///
    /// Examples
    /// --------
    /// This example should output 0.0:
    /// ```rust
    /// fn main() {
    ///     let mut res:f64 = root_search::binary(&test_function::identity, (-1,2), 0.1)
    ///     res = res.round();
    ///     println!({}, res);
    /// }
    /// ```
    /// ---
    /// This example should output 0.5:
    /// ```rust
    ///fn main() {
    ///     let mut res:f64 = root_search::binary(&test_function::identity, (-1,2), 0.1)
    ///     res = (res*10.0).round()/10.0;
    ///     println!({}, res);
    /// }
    /// ```
    pub fn binary(
        func: &ContinuousFunction,
        domain: (f64, f64),
        trunc_err: f64,
    ) -> Result<f64, &'static str> {
        // Since we only care about signs, we get a tiny speedup by working with the signs of the function
        // values, rather than the values themselves.
        let func_sgn = compose_fn!(func => signum);

        let (mut start, mut end): (f64, f64) = domain;
        let (mut start_val, end_val): (f64, f64) = (func_sgn(start), func_sgn(end));

        if start_val == 0.0 {
            return Ok(start);
        }
        if end_val == 0.0 {
            return Ok(end);
        }
        if start_val * end_val > 0.0 {
            return Err("Error: no sign change at endpoints!");
        }

        let val: f64 = loop {
            let midpoint: f64 = (end + start) / 2.0;
            let test_val: f64 = func_sgn(midpoint);

            if test_val == 0.0 || (end - start) < trunc_err {
                break midpoint;
            }
            if test_val * start_val > 0.0 {
                start = midpoint;
                start_val = test_val;
            } else {
                end = midpoint;
            }
        };

        return Ok(val);
    }

    /// Return the fixed point of a function where one exists.
    ///
    /// Parameters
    /// ----------
    /// * `func` : A continuous function with a fixed point, e.g. a contraction mapping.
    /// * `initial_val` : An initial guess for the location of the fixed point.
    /// * `trunc_err` : A float representing the acceptable truncation error for the search; e.g. `trunc_err=1` will result in finding the fixed point +-1.
    /// * `max_iter` : The maximum number of iterations the algorithm will use before it declares there is no fixed point. This is highly dependent on both the rate of convergence and the truncation error.
    ///
    /// Returns
    /// -------
    /// * `f64` : A float representing the fixed point.
    /// * `func_vals` : A vector of all computed iterations.
    pub fn fixed_point(
        func: &ContinuousFunction,
        initial_val: f64,
        trunc_err: f64,
        max_iter: usize,
    ) -> (f64, Vec<f64>) {
        let mut func_vals: Vec<f64> = Vec::with_capacity(max_iter);
        let mut current_val: f64 = initial_val;
        for _ in 1..max_iter {
            func_vals.push(current_val);
            let next_val: f64 = func(current_val);
            if (next_val - current_val).abs() < trunc_err {
                current_val = next_val;
                break;
            }
            current_val = next_val;
        }
        return (current_val, func_vals);
    }
}

/// A collection of 'functionals' for use in fixed-point iteration.
///
/// Functionals
/// -----------
/// * `x_minus` : id-F.
/// * `identity` : F
/// * `frac` : F/(2+k) unless k=-2, in which case returns the constant function 1; since x -> x-1 has no fixed point, this behaviour is fine.
#[allow(dead_code)]
mod functional {
    use ContinuousFunction;

    pub fn x_minus<'a>(func: Box<ContinuousFunction<'a>>) -> Box<ContinuousFunction<'a>> {
        Box::new(move |x: f64| -> f64 { x - func(x) })
    }

    pub fn identity<'a>(func: &'a ContinuousFunction) -> Box<ContinuousFunction<'a>> {
        Box::new(move |x: f64| -> f64 { func(x) })
    }

    pub fn frac<'a, 'b: 'a>(
        func: &'a ContinuousFunction,
        k: &'b f64,
    ) -> Box<ContinuousFunction<'a>> {
        if *k == -2.0 {
            Box::new(|_x: f64| -> f64 { 1.0 })
        } else {
            Box::new(move |x: f64| -> f64 { &func(x) / (2.0 + *k) })
        }
    }
}

/// A collection of 'test functions' for analysing the root-finders.
///
/// Functions
/// ---------
/// * `identity` : Returns the input. Has a root at `x=0`.
/// * `polynom` : A cubic polynomial. Has roots at `x=0.5` and `x=4`.
/// * `trig` : The sum of a linear polynomial and a sinosoidal function. Has a root at `x=-2.88...`
#[allow(dead_code)]
mod test_function {
    use pow;
    /// Test
    pub fn identity(x: f64) -> f64 {
        return x;
    }

    pub fn polynom(x: f64) -> f64 {
        let res: f64 = pow(x, 3) - 8.5 * pow(x, 2) + 20.0 * x - 8.0;
        return res;
    }

    pub fn trig(x: f64) -> f64 {
        let res: f64 = 2.0 * x - 3.0 * x.sin() + 5.0;
        return res;
    }
}

/// The derivatives of 'test functions', for use in e.g. the Raphson-Newton algorithm.
///
/// Functions
/// ---------
/// * `identity` : Returns the input. Has a root at `x=0`.
/// * `polynom` : A cubic polynomial. Has roots at `x=0.5` and `x=4`.
/// * `trig` : The sum of a linear polynomial and a sinosoidal function. Has a root at `x=-2.88...`
#[allow(dead_code)]
mod test_func_derivatives {
    use pow;

    pub fn identity(_x: f64) -> f64 {
        return 1.0;
    }

    pub fn polynom(x: f64) -> f64 {
        let res: f64 = 3.0 * pow(x, 2) - 17.0 * x + 20.0;
        return res;
    }

    pub fn trig(x: f64) -> f64 {
        let res: f64 = 2.0 - 3.0 * x.cos();
        return res;
    }
}
