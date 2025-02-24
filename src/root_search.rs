/*!
A collection of algorithms for use in root-finding.

Functions
---------
* `binary` : Binary search, a.k.a interval bisection.
* `fixed_point` : Fixed point iteration.
*/

use ContinuousFunction;
use {compose_fn, signum};

/**
Find a root of a continuous function using binary search.

Parameters
----------
* `func` : A continuous function with a sign change over the given domain.
* `domain` : The start and end points of the search interval. Note that `func` must be computable over the entire domain, including end points.
* `trunc_err` : A float representing the acceptable truncation error for the search; e.g. `trunc_err=1` will result in finding the root +-1.

Returns
-------
* `f64` : A floating point representing the root of the function.

Errors
------
* If the function doesn't change sign at the endpoints

Examples
--------
In this example, `root_search::binary` should return 0.0 ± 0.1.
```rust
let res:f64 = root_search::binary(&test_function::identity, (-1.0,2.0), 0.1)?;
Ok(assert_eq!(0.0, (res*10.0).round()/10.0))
```
---
In this example, `root_search::binary` should return -2.9 ± 0.1.
```rust
let res:f64 = root_search::binary(&test_function::trig, (-3.0,-2.0), 0.1)?;
Ok(assert_eq!(-2.9, (res*10.0).round()/10.0))
```
*/
pub fn binary(
    func: &ContinuousFunction,
    domain: (f64, f64),
    trunc_err: f64,
) -> Result<f64, Box<&'static str>> {
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
        return Err(Box::new("Error: no sign change at endpoints!"));
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

/**
Return the fixed point of a function where one exists.

Parameters
----------
* `func` : A continuous function with a fixed point, e.g. a contraction mapping.
* `initial_val` : An initial guess for the location of the fixed point.
* `trunc_err` : A float representing the acceptable truncation error for the search; e.g. `trunc_err=1` will result in finding the fixed point +-1.
* `max_iter` : The maximum number of iterations the algorithm will use before it declares there is no fixed point. This is highly dependent on both the rate of convergence and the truncation error.

Returns
-------
* `f64` : A float representing the fixed point.
* `Vec<f64>` : A vector of all computed iterations.

Errors
------
* `Vec<f64>` : If function fails to converge in `max_iter`, returns the current sequence of computed iterations as an error.

Examples
--------
In this example, `root_search::fixed_point` should return 0.8 ± 0.1, the fixed point of the cos function.
```rust
let (res, _seq) = root_search::fixed_point(&(|x:f64| -> f64 {x.cos()}), 2.5, 0.1, 10)?;
Ok(assert_eq!(0.8, (res*10.0).round()/10.0))
```
*/
pub fn fixed_point(
    func: &ContinuousFunction,
    initial_val: f64,
    trunc_err: f64,
    max_iter: usize,
) -> Result<(f64, Vec<f64>), Vec<f64>> {
    let mut func_vals: Vec<f64> = Vec::with_capacity(max_iter);
    let mut current_val: f64 = initial_val;
    for _ in 1..max_iter {
        func_vals.push(current_val);
        let next_val: f64 = func(current_val);
        if (next_val - current_val).abs() < trunc_err {
            current_val = next_val;
            return Ok((current_val, func_vals))
        }
        current_val = next_val;
    }
    func_vals.push(current_val);
    Err(func_vals)
}
