/*!
A collection of 'functionals' for use in fixed-point iteration.

Functionals
-----------
* `x_minus` : id-F.
* `identity` : F
* `frac` : F/(2+k) unless k=-2, in which case returns the constant function 1; since x -> x-1 has no fixed point, this behaviour is fine and avoids complicated return types.
* `newton_raphson` : F/F'
*/

use ContinuousFunction;

/**
Inverts a function, i.e. applies the transform `f(x) -> 1-f(x)`.

Parameters
----------
* `func` : A continuous function.

Returns
-------
* `ContinuousFunction` : The inversion of `func`. 

Examples
--------
In this example, `functional::x_minus` should return f(x) = x-cos(x), which we test at x=pi.
```rust
let res:f64 = functional::x_minus(Box::new(|x:f64| -> f64 {x.cos()}))(3.1415);
assert_eq!(4.1, (res*10.0).round()/10.0);
```
*/
pub fn x_minus<'a>(func: Box<ContinuousFunction<'a>>) -> Box<ContinuousFunction<'a>> {
    Box::new(move |x: f64| -> f64 { x - func(x) })
}

/**
Applies the identity transform to a function.

Parameters
----------
* `func` : A continuous function.

Returns
-------
* `ContinuousFunction` : The inversion of `func`. 

Examples
--------
In this example, `functional::identity` should return f(x) = cos(x), which we test at x=pi.
```rust
let res:f64 = functional::identity(&(|x:f64| -> f64 {x.cos()}))(3.1415);
assert_eq!(-1.0, (res*10.0).round()/10.0);
```
*/
pub fn identity<'a>(func: &'a ContinuousFunction) -> Box<ContinuousFunction<'a>> {
    Box::new(func)
}

/**
Applies the identity transform to a function.

Parameters
----------
* `func` : A continuous function.

Returns
-------
* `func`

Examples
--------
In this example, `functional::identity` should return f(x) = cos(x), which we test at x=pi.
```rust
let res:f64 = functional::identity(&(|x:f64| -> f64 {x.cos()}))(3.1415);
assert_eq!(-1.0, (res*10.0).round()/10.0);
```
*/
pub fn frac<'a, 'b: 'a>(func: &'a ContinuousFunction, k: &'b f64) -> Box<ContinuousFunction<'a>> {
    if *k == -2.0 {
        Box::new(|_x: f64| -> f64 { 1.0 })
    } else {
        Box::new(move |x: f64| -> f64 { &func(x) / (2.0 + *k) })
    }
}

/**
Applies the Newton-Raphson transform to a function.

Parameters
----------
* `func` : A continuous function.
* `deriv` : The derivative of `func`.

Returns
-------
* `ContinuousFunction` : Input function divided by its derivative, except where the derivative is 0.

Examples
--------
In this example, `functional::newton_raphson` should return f(x) = tan(x), which we test at x=pi/4.
```rust
let res:f64 = functional::newton_raphson(&(|x:f64| -> f64 {x.sin()}), &(|x:f64| -> f64 {x.cos()}))(3.1415/4.0);
assert_eq!(1.0, (res*10.0).round()/10.0);
```
*/
pub fn newton_raphson<'a, 'b: 'a>(
    func: &'a ContinuousFunction,
    deriv: &'b ContinuousFunction,
) -> Box<ContinuousFunction<'a>> {
    Box::new(move |x: f64| -> f64 {
        if deriv(x) == 0.0 {
            1.0
        } else {
            &func(x) / deriv(x)
        }
    })
}
