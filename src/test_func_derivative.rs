/*!
The derivatives of 'test functions', for use in e.g. the Raphson-Newton algorithm.

Functions
---------
* `identity` : Returns the input. Has a root at `x=0`.
* `polynom` : A cubic polynomial. Has roots at `x=0.5` and `x=4`.
* `trig` : The sum of a linear polynomial and a sinosoidal function. Has a root at `x=-2.88...`
*/

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
