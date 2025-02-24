/*!
A collection of 'test functions' for analysing the root-finders.

Functions
---------
* `identity` : Returns the input. Has a root at `x=0`.
* `polynom` : A cubic polynomial. Has roots at `x=0.5` and `x=4`.
* `trig` : The sum of a linear polynomial and a sinosoidal function. Has a root at `x=-2.88...`
*/

use pow;

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
