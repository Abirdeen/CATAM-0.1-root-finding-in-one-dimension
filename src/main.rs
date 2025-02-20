extern crate num;
use num::signum;
use num::traits::pow;

extern crate composing;
use composing::compose_fn;

fn main() {
    let trunc_err: f64 = 0.00000001;
    let res: f64 = root_search::binary(&test_functions::trig, (-10.0, 15.0), trunc_err);

    let rounded_res = (res / trunc_err).round() * trunc_err;

    println!("A zero for your function was found at {rounded_res}+-{trunc_err}");

}

#[allow(dead_code)]
mod root_search {
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
    ///     let res:f64 = root_search::binary(&test_functions::identity, (-1,2), 0.1).round();
    ///     println!({}, res);
    /// }
    /// ```
    /// ---
    /// This example should output 0.5:
    /// ```rust
    ///fn main() {
    ///     let res:f64 = root_search::binary(&test_functions::identity, (-1,2), 0.1).round();
    ///     println!({}, res);
    /// }
    /// ```
    pub fn binary(func: &dyn Fn(f64) -> f64, domain: (f64, f64), trunc_err: f64) -> f64 {
        // Since we only care about signs, we get a tiny speedup by working with the signs of the function
        // values, rather than the values themselves.
        let func_sgn = compose_fn!(func => signum);

        let (mut start, mut end): (f64, f64) = domain;
        let (mut start_val, end_val): (f64, f64) = (func_sgn(start), func_sgn(end));

        if start_val == 0.0 {
            return start;
        }
        if end_val == 0.0 {
            return end;
        }
        if start_val * end_val > 0.0 {
            // This should actually return an error, need to figure out error types later.
            return 1.0;
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

        return val;
    }

    /// Return the fixed point of a function where one exists.
    ///
    /// Parameters
    /// ----------
    /// * `func` : A continuous function with a fixed point, e.g. a contraction mapping.
    /// * `trunc_err` : A float representing the acceptable truncation error for the search; e.g. `trunc_err=1` will result in finding the fixed point +-1.
    /// * `max_iter` : The maximum number of iterations the algorithm will use before it declares there is no fixed point. This is highly dependent on both the rate of convergence and the truncation error.
    ///
    /// Returns
    /// -------
    /// * `f64` : A float representing the fixed point.
    /// * `func_vals` : A vector of all computed iterations.
    pub fn fixed_point(func: &dyn Fn(f64) -> f64, trunc_err: f64, max_iter: usize) -> (f64, Vec<f64>) {
        let mut func_vals = Vec::with_capacity(max_iter);
        let mut current_val: f64 = 1.0;
        for _ in 1..max_iter {
            func_vals.push(current_val);
            let next_val: f64 = func(current_val);
            if (next_val-current_val).abs() < trunc_err {
                current_val = next_val;
                break;
            }
            current_val = next_val;
        }
        return (current_val,func_vals);
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
mod test_functions {
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

// trait Function {
//     fn call(&self, x: f64) -> f64;
//   }
// trait Derivable: Function {
//     fn d_dx(&self) -> Box<dyn Function>;
// }

// macro_rules! def_derivable {
//     ($func_name:ident, $var:ident, $func:expr, $deriv_name:ident, $dvar:ident, $deriv:expr) => {
        
//         struct $func_name;
//         impl Function for $func_name {
//         fn call(&self, $var: f64) -> f64  {
//             $func
//         }
//         }
//         impl Derivable for $func_name {
//         fn d_dx(&self) -> Box<dyn Function> {
//             Box::new($deriv_name)
//         }
//         }

//         struct $deriv_name;
//         impl Function for $deriv_name {
//         fn call(&self, $dvar: f64) -> f64  {
//             $deriv
//         }
//         }

//     };
// }

// impl Function for Box<dyn Function> {
//     fn call(&self, x: f64) -> f64 {
//         self.downcast()
//     }
// }

// def_derivable!(Identity, x, x, DIdentity, _x, 1.0);

// fn quotient<F: Derivable>(func: F) -> impl Function {
//     struct Quotient<F1, F2>(F1, F2);
//     impl<F1: Function, F2: Function> Function for Quotient<F1, F2> {
//         fn call(&self, x: f64) -> f64 {
//           self.0.call(x) / self.1.call(x)
//         }
//     }
//     let derivative = func.d_dx();
//     Quotient(func, derivative)
// }