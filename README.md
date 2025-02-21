# CATAM-0.1-root-finding-in-one-dimension

A rewrite of my Cambridge CATAM undergraduate project, "Root finding in one dimension".

## Rust

This project is written in [Rust](https://www.rust-lang.org/), which is a highly performant compiled language with memory management. For this project, where we may manipulate large volumes of data, performance is prioritised. Rust is a multi-paradigm language, supporting object-oriented, functional and data-oriented designs. This project mostly follows a functional design - after all, we're attempting to study properties of functions - with some object-oriented aspects.

## The project

The aim of this project is to study iteration methods for the numerical solution of an algebraic
or transcendental equation $F(x) = 0$.

### Test functions

We focus on three main example functions:
- $F(x)=x$, with a simple root at $x=0$;
- $F(x)=2x-3\sin(x)+5$, with a simple root at $x\approx -2.88$; and
- $F(x)=x^3-8.5x^2+20x-8=(x-\frac{1}{2})(x-4)^2$, with a simple root at $x=\frac{1}{2}$ and a repeated root at $x=4$.

These functions are very well-behaved (they are [smooth](https://en.wikipedia.org/wiki/Smoothness), have bounded 3rd derivative, etc.), but only some of these nice properties are used in each iteration method.

These functions are implemented in the ```test_functions``` module as ```identity```, ```trig``` and ```polynom``` respectively.

### Binary search

The binary search algorithm is based on Bolzano's Theorem, a special case of the [Intermediate Value Theorem](https://en.wikipedia.org/wiki/Intermediate_value_theorem). In short, if a continuous function $F$ changes sign in an interval $[a,b]$, then it must have a root in the interval. Taking the midpoint $m=\frac{a+b}{2}$, and considering the sign of $F(m)$, we can determine if the root is in $[a,m]$ or $[m,b]$. Repeating this process, we can locate the root to arbitrary precision.

A binary search will locate any simple root on any continuous function with a linear [rate of convergence](https://en.wikipedia.org/wiki/Rate_of_convergence).

We implement this algorithm as ```root_search::binary```.

### Fixed point iteration

[Fixed point iteration](https://en.wikipedia.org/wiki/Fixed-point_iteration) is a method for finding fixed points for some function $f$, that is, values for which $f(x)=x$. 

The sequence $x_{N+1}=f(x_N)$ may be convergent or divergent for a given starting value $x_0$. If it is convergent, then by [sequential continuity of f](https://en.wikipedia.org/wiki/Continuous_function#Definition_in_terms_of_limits_of_sequences), the limit $x_*$ will be a fixed point of $f$.

When this method is convergent, it is at least linearly convergent for simple roots.

We implement a fixed point iteration algorithm as ```root_search::fixed_point```.

### A nice class of functionals for fixed-point iteration

If $\Gamma$ is any functional such that $\Gamma(F)(x)=0$ exactly when $F(x)=0$, then roots of $F$ correspond to fixed points of $f(x) = x - \Gamma(F)(x)$. We can thus leverage our fixed-point iteration algorithm to find roots of $F$. Picking $\Gamma$ carefully is key to both avoiding divergence, and to having a fast rate of convergence.

Two choices of $\Gamma$ are implemented as functionals in our program.

- $\Gamma(F) = \frac{F}{2+k}$ is implemented as `functional::frac`.

- $\Gamma(F) = \frac{F}{F'}$, where $F'$ is the derivative of $F$, is implemented as `functional::newton_raphson`. $F'$ must be manually defined and passed as an input.

Some specific methods, like the Newton-Raphson method, have better convergence. These often depend on properties like differentiability, so are not always an appropriate choice.

