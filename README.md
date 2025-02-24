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

The binary search algorithm is based on Bolzano's Theorem, a special case of the [Intermediate Value Theorem](https://en.wikipedia.org/wiki/Intermediate_value_theorem). In short, if a continuous function $F$ changes sign in an interval $[a,b]$, then it must have a root in the interval. By taking the midpoint $m=\frac{a+b}{2}$, and considering the sign of $F(m)$, we can determine if the root is in $[a,m]$ or $[m,b]$. Repeating this process, we can locate the root to arbitrary precision.

A binary search will locate any simple root on any continuous function with a linear [rate of convergence](https://en.wikipedia.org/wiki/Rate_of_convergence).

We implement this algorithm as ```root_search::binary```.

### Fixed point iteration

[Fixed point iteration](https://en.wikipedia.org/wiki/Fixed-point_iteration) is a method for finding fixed points for some function $f$, that is, values for which $f(x)=x$. 

The sequence $x_{N+1}=f(x_N)$ may be convergent or divergent for a given starting value $x_0$. If it is convergent, then by [sequential continuity](https://en.wikipedia.org/wiki/Continuous_function#Definition_in_terms_of_limits_of_sequences) of $f$, the limit $x_*$ will be a fixed point of $f$.

When this method is convergent, it is at least linearly convergent for simple roots.

We implement a fixed point iteration algorithm as ```root_search::fixed_point```.

### A nice class of functionals for fixed-point iteration

If $\Gamma$ is any functional such that $\Gamma(F)(x)=0$ exactly when $F(x) = 0$, then roots of $F$ correspond to fixed points of $f(x) = x - \Gamma(F)(x)$. We can thus leverage our fixed-point iteration algorithm to find roots of $F$. Picking $\Gamma$ carefully is key to both avoiding divergence, and to having a fast rate of convergence.

Two choices of $\Gamma$ are implemented as functionals in our program.

- $\Gamma(F) = \frac{F}{2+k}$ is implemented as `functional::frac`.

- $\Gamma(F) = \frac{F}{F'}$, where $F'$ is the derivative of $F$, is implemented as `functional::newton_raphson`. $F'$ must be manually defined and passed as an input.

Some specific methods, like the Newton-Raphson method, have better convergence. These often depend on properties like differentiability, so are not always an appropriate choice.

## Problems

The original CATAM project involved certain explicit questions and problems, which are reproduced (and solved) here.

### Problem one: 

Show, with the help of a graph, that $F(x) = 2x - 3sin(x) + 5$ has exactly one root.

#### Solution: 

We can see that $2x + 8 = 2x + 3 + 5 \ge F(x) \ge 2x - 3 + 5 = 2x + 2$, since $1 \ge sin(x) \ge -1$. So $F(x)<0$ when $x<-4$, and $F(x)>0$ when $x>-1$. Graphing $F$ in the range $-4 \le x \le -1$, we can see there is exactly one root:

<img src="./images/trig-graph-desmos.png" alt="Plot of F" width="400"/>

---

### Problem two:

Write an implementation of binary search

#### Solution:

Implemented in `root_search::binary`.

---

### Problem three:

Write an implementation of fixed-point iteration.

#### Solution:

Implemented in `root_search::fixed_point`.

---

### Problem four:

Use fixed-point iteration to find the root of $F$ by taking the transform $f = x - \Gamma(F) = x - \frac{F}{2+k}$.

1. First, run the program with $k = 0$, a truncation error of $10^{-5}$, $x_0 = -2$, and $N_{max}=10$. Plot $y = f(x)$ and $y = x$ on the same graph, and use these plots to show why convergence should not occur. Explain the divergence by identifying a theoretical criterion that has been violated.

2. Determine the values of $k$ for which convergence is guaranteed if $x_N$ remains in the range $(−\pi, −\pi/2)$.

3. Choose, giving reasons, a value of $k$ for which monotonic convergence should occur near the root, and also a value for which oscillatory convergence should occur near the root. Verify that these two values of $k$ give the expected behaviour, by running the program with $N_{max} = 20$. 

4. Also run the case $k = 16$. This should converge only slowly, so set $N_{max} = 50$.

5. Discuss whether your results are consistent with first-order convergence.

#### Solution:

For the whole of this question, we use the following code in `main.rs`, varying our parameters and making minor modifications as needed:
```rust
fn main() -> Result<(), Vec<64>> {
    let initial_func: &ContinuousFunction = &(test_function::trig as fn(f64) -> f64);
    let k: f64 = 0.0;
    let initial_val: f64 = -2.0;
    let trunc_err: f64 = 1.0/pow(10.0, 5);
    let max_iter: usize = 10;

    let transformed_func: Box<ContinuousFunction> = functional::x_minus(functional::frac(&initial_func, &k));
    let (mut res, _seq) = root_search::fixed_point(&transformed_func, initial_val, trunc_err, max_iter)?;
    res = (res/trunc_err).round()*trunc_err;
    println!("Root is at {} ± {}", res, trunc_err);
}
```

1. Running our code with the given parameters, a root is reported at -3.9749 ± 0.00001. From our previous work, we know this isn't correct. The true root should be between -3 and -2. In fact, a binary search reveals the root is at $x \approx -2.8832$.

    For $x_N$ close to $x_*$, we can consider the Taylor expansion of $f$ to see that $x_N = f(x_{N-1}) = f(x_*) + f'(x_*)(x_{N-1} - x_*) + ... = x_* + f'(x_*)(x_{N-1} - x_*) + ...$

    Writing $\epsilon_N = x_N-x_*$, we find that 

    $\epsilon_N \approx f'(x_*)\epsilon_{N-1}$. 

    So if $|f'(x_*)| \cong |\frac{3cos(x_*) + k}{2 + k}| > 1$, then the iteration will diverge. For $k = 0$, $|f'(-2.8832)| \approx 1.4502 > 1$, so divergence is expected.

2. The [mean-value theorem](https://en.wikipedia.org/wiki/Mean_value_theorem) says that, for a function continuous on $[a,b]$ and differentiable on $(a,b)$, there's some $\xi \in (a,b)$ such that $f(b)-f(a) = f'(\xi)(b-a)$. 

    In particular, $\epsilon_N = f'(\xi)\epsilon_{N-1}$ for some $\xi \in (x_{N-1}, x_*)$. From this, we can deduce that if $|f'(x)| < 1$ for all $x$ in some interval $[a,b]$, then $f$ is a contraction mapping on $[a,b]$, so the iterations will converge here.

    With $[a,b] = [-\pi, -\frac{\pi}{2}]$, and $k > \frac{1}{2}$, we can explicitly compute $|f'(x)| = |\frac{3cos(x) + k}{2 + k}| < \min(|\frac{k - 3}{k + 2}|, |\frac{k}{k+2}|) < 1$.

3. We expect oscillatory convergence when $f'(x) < 0$ near the root, and monotonic convergence when $f'(x) > 0$ near the root. Plotting $f'(x_*) = \frac{3cos(x_*)+k}{2+k}$ as $k$ varies, we find a sign change at the critical value $k_c = 2.9 ± 0.01$. Picking values of $k$ slightly larger or smaller than $k_c$, we should expect rapid monotonic and oscillatory convergence, respectively.

    For $k = 2.85$, a truncation error of $10^{-10}$, $x_0 = -2$, and $N_{max}=20$, we find $x_* = -2.8832368726 ± 0.0000000001$ and the following table of iterates:

    | $N$ | $x_N$        |
    | --- | ------------ |
    | 0   | -2           |
    | 1   | -2.768637583 |
    | 2   | -2.883242049 |
    | 3   | -2.883236818 |
    | 4   | -2.883236873 |
    | 5   | -2.883236872 |

    For $k = 2.95$ and other parameters the same, we find $x_* = -2.8832368726 ± 0.0000000001$ and the following table of iterates:

    | $N$ | $x_N$        |
    | --- | ------------ |
    | 0   | -2           |
    | 1   | -2.753109551 |
    | 2   | -2.880409724 |
    | 3   | -2.883207942 |
    | 4   | -2.883236582 |
    | 5   | -2.883236869 |
    | 6   | -2.883236872 |

4. Running the case $k = 16$ with $N_{max} = 50$, a truncation error of $10^{-5}$ and $x_0 = -2$, we get the following table of iterates:

    | $N$ | $x_N$     | $N$ | $x_N$     | $N$ | $x_N$     |
    | --- | --------- | --- | --------- | --- | --------- |
    | 1   | -2.207105 | 12  | -2.860369 | 23  | -2.882541 |
    | 2   | -2.373698 | 13  | -2.866583 | 24  | -2.882730 |
    | 3   | -2.503502 | 14  | -2.871111 | 25  | -2.882868 |
    | 4   | -2.602390 | 15  | -2.874409 | 26  | -2.882968 |
    | 5   | -2.676588 | 16  | -2.876810 | 27  | -2.883041 |
    | 6   | -2.731705 | 17  | -2.878559 | 28  | -2.883094 |
    | 7   | -2.772378 | 18  | -2.879832 | 29  | -2.883133 |
    | 8   | -2.802260 | 19  | -2.880758 | 30  | -2.883161 |
    | 9   | -2.824152 | 20  | -2.881433 | 31  | -2.883182 |
    | 10  | -2.840158 | 21  | -2.881924 | 32  | -2.883197 |
    | 11  | -2.851844 | 22  | -2.882281 | 33  | -2.883207 |

5. Since $\epsilon_N \approx f'(x_*)\epsilon_{N-1}$ by our work in part 1, we know that $\lim\limits_{N \to \infty}(|\frac{\epsilon_N}{\epsilon_{N-1}}|) = |f'(x_*)|$, so when $0 < |f'(x_*)| < 1$, we get linear convergence. This is consistent with our results above.

    Moreover, when $|f'(x_*)| < \frac{1}{2}$, this method is faster than interval bisection.

---

### Problem five:

Use your program to find the double root of equation $G(x) = x^3 - 8.5x^2 + 20x - 8$ by taking $g = x-\Gamma(G) = x-\frac{G}{20}$. 

By considering $g'(x_*)$ explain why convergence will be slow at a multiple root for any
choice of functional $\Gamma$ given by post-composition with a differentiable function.

In your calculations some care may be needed over the choice of $x_0$. Since convergence will be slow, take $N_{max} = 1000$. Is this an example of first-order convergence?

#### Solution:

For the whole of this question, we use the following code in `main.rs`, varying our parameters and making minor modifications as needed:
```rust
fn main() -> Result<(), Vec<f64>> {
    let initial_func: &ContinuousFunction = &(test_function::polynom as fn(f64) -> f64);
    let initial_val: f64 = 4.5;
    let trunc_err: f64 = 1.0/pow(10.0, 5);
    let max_iter: usize = 1000;

    let transformed_func: Box<ContinuousFunction> = functional::x_minus(functional::frac(&initial_func, &18.0));
    let (mut res, _seq) = root_search::fixed_point(&transformed_func, initial_val, trunc_err, max_iter)?;
    res = (res/trunc_err).round()*trunc_err;
    println!("Root is at {} ± {}", res, trunc_err);
}
```

We use this code to construct the following table:

| $N$ | $x_N$    | $\epsilon_N$ | $\epsilon_N/\epsilon_{N-1}$     | $7N\epsilon_N/40$ |
| --- | -------- | ------------ | ------------------------------- | ----------------- |
| 0   | 4.5      | 0.5          |                                 |                   |
| 1   | 4.45     | 0.45         | 0.9                             | 0.07875           |
| 2   | 4.410006 | 0.410006     | 0.820013                        | 0.143502          |
| 3   | 4.377142 | 0.377142     | 0.754283                        | 0.197999          |
| 4   | 4.349568 | 0.349568     | 0.699136                        | 0.244698          |
| 5   | 4.326048 | 0.326048     | 0.652096                        | 0.285292          |
| ⋮   | ⋮         | ⋮            | ⋮                                | ⋮                 |
| 731 | 4.007585 | 0.007585     | 0.01517                         | 0.970318          |
| 732 | 4.007575 | 0.007575     | 0.01515                         | 0.970353          |
| 733 | 4.007565 | 0.007565     | 0.01513                         | 0.970388          |
| 734 | 4.007555 | 0.007555     | 0.01511                         | 0.970423          |
| 735 | 4.007545 | 0.007545     | 0.01509                         | 0.970457          |

We see from our table that $x_N$ converges very slowly to 4, that the ratio $\epsilon_N/\epsilon_{N-1}$ slowly shrinks over (suggesting better than linear convergence), and that the ratio $7N\epsilon_N/40$ slowly tends to $1$.

Now, if $\Gamma(G) = h\circ G$ for some differentiable function $h$, then we can compute $g' = 1 - \frac{dh}{dG} * \frac{dG}{dx}$. When $x_*$ is a repeated root, $G'(x_*) = 0$, so $g'(x_*) = 1$. Since $\epsilon_N \approx g'(x_*)\epsilon_{N-1}$, it follows that the convergence or divergence of the iterations will be small, and depend on higher-order powers of $\epsilon_{N-1}$.

We can examine the Taylor expansion further:

$x_N = x_{N-1} + \frac{1}{2}g''(x_*)\epsilon_{N-1}^2 + O(\epsilon_{N-1}^3)$,

so $\epsilon_N = \epsilon_{N-1} + \frac{1}{2} g''(x_*) \epsilon_{N-1}^2 + O(\epsilon_{N-1}^3)$.

We can compute $g''(x) = -(6x + 17)/20$, so $g''(x_*) = -7/20$. So if $\epsilon_N \sim \frac{\kappa}{N}$ for some constant $\kappa$, then $\kappa = 40/7$, explaining the behaviour of the last column in our previous table.

---
