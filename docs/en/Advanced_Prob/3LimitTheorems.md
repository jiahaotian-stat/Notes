# Chapter 3: Limit Theorems, Stochastic Convergence, and Infinitely Divisible Distributions

This chapter begins with the classical limit theorems for Bernoulli sums, then systematically introduces almost sure convergence, convergence in probability, and convergence in distribution. We discuss the continuous mapping theorem, Slutsky's theorem, Lévy's continuity theorem, the Borel-Cantelli lemmas, and subsequence characterizations. Finally, we introduce infinitely divisible distributions, their basic properties and standard examples, the Lévy-Khintchine representation, and their connection with limits of triangular arrays.

---

## 1. Classical Limit Theorems for Bernoulli Sums

Let

\[
S_n=\sum_{i=1}^n\xi_i,
\]

where $\xi_1,\ldots,\xi_n$ are independent and identically distributed with

\[
P(\xi_i=1)=p,
\qquad
P(\xi_i=0)=1-p.
\]

Then $S_n\sim\operatorname{Binomial}(n,p)$.

### 1.1 Bernoulli Weak Law of Large Numbers

!!! success "Theorem 1.1 (Bernoulli Weak Law of Large Numbers)"

    As $n\to\infty$,

    \[
    \frac{S_n}{n}\xrightarrow{P}p.
    \]

??? proof "Proof of Theorem 1.1 (click to expand)"

    Since

    \[
    E\left(\frac{S_n}{n}\right)=p,
    \qquad
    \operatorname{Var}\left(\frac{S_n}{n}\right)
    =\frac{p(1-p)}{n},
    \]

    for every $\varepsilon>0$, Chebyshev's inequality gives

    \[
    P\left(\left|\frac{S_n}{n}-p\right|>\varepsilon\right)
    \leq
    \frac{p(1-p)}{n\varepsilon^2}
    \longrightarrow0.
    \]

    Therefore, $S_n/n\xrightarrow{P}p$. $\square$

### 1.2 Bernoulli Central Limit Theorem

!!! success "Theorem 1.2 (De Moivre-Laplace Central Limit Theorem)"

    If $0<p<1$, then

    \[
    Z_n
    =\frac{S_n-np}{\sqrt{np(1-p)}}
    \xrightarrow{d}N(0,1).
    \]

Let $F_n$ be the distribution function of $Z_n$ and let $\Phi$ be the standard normal distribution function. Since $\Phi$ is continuous, Pólya's theorem further yields

\[
\Delta_n
=\sup_{x\in\mathbb{R}}|F_n(x)-\Phi(x)|
\longrightarrow0.
\]

Thus, for large $n$, standardized binomial tail probabilities can be approximated by normal tail probabilities:

\[
P(Z_n>x)
=1-F_n(x)
\approx1-\Phi(x).
\]

### 1.3 Poisson Limit Theorem

For each $n$, let $\xi_{n,1},\ldots,\xi_{n,n}$ be independent and identically distributed with

\[
P(\xi_{n,i}=1)=p_n,
\qquad
P(\xi_{n,i}=0)=1-p_n.
\]

Define

\[
S_n=\sum_{i=1}^n\xi_{n,i}.
\]

!!! success "Theorem 1.3 (Poisson Limit Theorem)"

    If

    \[
    np_n\longrightarrow\lambda>0,
    \]

    then

    \[
    S_n\xrightarrow{d}S,
    \qquad
    S\sim\operatorname{Poisson}(\lambda).
    \]

??? proof "Proof of Theorem 1.3 (click to expand)"

    For every fixed nonnegative integer $k$,

    \[
    P(S_n=k)
    =\binom{n}{k}p_n^k(1-p_n)^{n-k}.
    \]

    Since $np_n\to\lambda$,

    \[
    \binom{n}{k}p_n^k
    =\frac{n(n-1)\cdots(n-k+1)}{k!}p_n^k
    \longrightarrow\frac{\lambda^k}{k!},
    \]

    and

    \[
    (1-p_n)^{n-k}\longrightarrow e^{-\lambda}.
    \]

    Therefore,

    \[
    P(S_n=k)
    \longrightarrow
    e^{-\lambda}\frac{\lambda^k}{k!}
    =P(S=k).
    \]

    The probabilities on the right sum to $1$ over $k\geq0$, so $S_n\xrightarrow{d}S$. $\square$

!!! note "Scope and Extensions"

    The independent and identically distributed assumption can be relaxed. For example, laws of large numbers and central limit theorems can be studied for dependent sequences, Markov chains, and processes satisfying mixing conditions. Extreme-value statistics often converge to extreme-value laws such as the Gumbel distribution. Random variables may also take values in a general Banach space, where limit theorems are closely related to the geometry of that space.

---

## 2. Three Modes of Convergence for Random Variables

Let $X_n$ and $X$ be defined on the same probability space $(\Omega,\mathcal{A},P)$.

### 2.1 Almost Sure Convergence

!!! info "Definition 2.1 (Almost Sure Convergence)"

    If there exists $\Omega_0\in\mathcal{A}$ with $P(\Omega_0)=1$ such that, for every $\omega\in\Omega_0$,

    \[
    X_n(\omega)\longrightarrow X(\omega),
    \]

    then $X_n$ is said to converge almost surely to $X$, written

    \[
    X_n\xrightarrow{a.s.}X.
    \]

Equivalently, for every $\varepsilon>0$,

\[
P\left(\limsup_{n\to\infty}
\{|X_n-X|>\varepsilon\}\right)=0.
\]

The limit superior of events is

\[
\limsup_{n\to\infty}A_n
=\bigcap_{N=1}^{\infty}\bigcup_{n\geq N}A_n,
\]

which is the event that $A_n$ occurs infinitely often.

### 2.2 Convergence in Probability

!!! info "Definition 2.2 (Convergence in Probability)"

    If, for every $\varepsilon>0$,

    \[
    P(|X_n-X|>\varepsilon)\longrightarrow0,
    \]

    then $X_n$ is said to converge in probability to $X$, written

    \[
    X_n\xrightarrow{P}X.
    \]

Markov's inequality provides a common sufficient condition. If there exists $r>0$ such that

\[
E|X_n-X|^r\longrightarrow0,
\]

then, for every $\varepsilon>0$,

\[
P(|X_n-X|>\varepsilon)
\leq\frac{E|X_n-X|^r}{\varepsilon^r}
\longrightarrow0,
\]

and hence $X_n\xrightarrow{P}X$.

### 2.3 Convergence in Distribution

Let $F_n$ and $F$ be the distribution functions of $X_n$ and $X$, respectively, and define

\[
D_F=\{x\in\mathbb{R}:F\text{ is continuous at }x\}.
\]

!!! info "Definition 2.3 (Convergence in Distribution)"

    If, for every $x\in D_F$,

    \[
    F_n(x)\longrightarrow F(x),
    \]

    then $X_n$ is said to converge in distribution to $X$, written

    \[
    X_n\xrightarrow{d}X.
    \]

Convergence in distribution depends only on the marginal distributions and does not require $X_n$ and $X$ to be defined on the same probability space. If

\[
X_n\xrightarrow{d}X
\quad\text{and}\quad
X_n\xrightarrow{d}Y,
\]

then $X$ and $Y$ have the same distribution.

!!! example "Example 2.4 (Why Only Continuity Points Are Required)"

    Let $X_n\equiv1/n$ and $X\equiv0$. Then $X_n\to X$, but at $x=0$,

    \[
    F_n(0)=0,
    \qquad
    F(0)=1.
    \]

    Since $0$ is a discontinuity point of $F$, this does not prevent $X_n\xrightarrow{d}X$.

---

## 3. Operational Rules for Convergence in Probability and Distribution

### 3.1 Algebra of Convergence in Probability

!!! success "Theorem 3.1 (Algebraic Properties of Convergence in Probability)"

    If

    \[
    X_n\xrightarrow{P}X,
    \qquad
    Y_n\xrightarrow{P}Y,
    \]

    then

    \[
    X_n\pm Y_n\xrightarrow{P}X\pm Y,
    \]

    \[
    X_nY_n\xrightarrow{P}XY,
    \]

    and, if $P(Y\neq0)=1$,

    \[
    \frac{X_n}{Y_n}\xrightarrow{P}\frac{X}{Y}.
    \]

### 3.2 Continuous Mapping Theorem

!!! success "Theorem 3.2 (Continuous Mapping Theorem)"

    If $X_n\xrightarrow{P}X$ and $f:\mathbb{R}\to\mathbb{R}$ is continuous, then

    \[
    f(X_n)\xrightarrow{P}f(X).
    \]

    If $X_n\xrightarrow{d}X$ and the discontinuity set $D_f$ of $f$ satisfies

    \[
    P(X\in D_f)=0,
    \]

    then

    \[
    f(X_n)\xrightarrow{d}f(X).
    \]

### 3.3 Convergence of Types and Slutsky's Theorem

!!! success "Theorem 3.3 (Convergence of Types)"

    Suppose that $X_n\xrightarrow{d}X$, where $X$ is nondegenerate. If

    \[
    a_nX_n+b_n\xrightarrow{d}Y,
    \]

    and $Y$ is nondegenerate, then there exist constants $a>0$ and $b\in\mathbb{R}$ such that

    \[
    a_n\longrightarrow a,
    \qquad
    b_n\longrightarrow b,
    \]

    and

    \[
    Y\overset{d}{=}aX+b.
    \]

!!! success "Theorem 3.4 (Slutsky's Theorem)"

    If

    \[
    X_n\xrightarrow{d}X,
    \qquad
    Y_n\xrightarrow{P}c,
    \]

    where $c$ is a constant, then

    \[
    X_n+Y_n\xrightarrow{d}X+c,
    \]

    \[
    X_nY_n\xrightarrow{d}cX,
    \]

    and, if $c\neq0$,

    \[
    \frac{X_n}{Y_n}\xrightarrow{d}\frac{X}{c}.
    \]

!!! warning "Marginal Convergence in Distribution Alone Is Insufficient"

    In general, $X_n\xrightarrow{d}X$ and $Y_n\xrightarrow{d}Y$ do not imply

    \[
    X_n+Y_n\xrightarrow{d}X+Y.
    \]

    If one additionally has joint convergence

    \[
    (X_n,Y_n)\xrightarrow{d}(X,Y),
    \]

    then the continuous mapping theorem gives

    \[
    X_n+Y_n\xrightarrow{d}X+Y.
    \]

---

## 4. Characteristic Functions and Weak Convergence

Let the characteristic functions of $X_n$ and $X$ be

\[
\varphi_n(t)=E(e^{itX_n}),
\qquad
\varphi(t)=E(e^{itX}).
\]

### 4.1 Lévy's Continuity Theorem

!!! success "Theorem 4.1 (Lévy's Continuity Theorem)"

    The following statements are equivalent:

    1. $X_n\xrightarrow{d}X$;
    2. for every $t\in\mathbb{R}$,

    \[
    \varphi_n(t)\longrightarrow\varphi(t).
    \]

    More generally, if $\varphi_n(t)$ converges pointwise to a function $\varphi(t)$ that is continuous at $0$, then $\varphi$ is the characteristic function of a probability distribution, and the corresponding distributions converge weakly to that distribution.

This theorem converts weak convergence of distribution functions into pointwise convergence of characteristic functions. It is a fundamental tool for proving central limit theorems and studying sums of independent random variables.

### 4.2 Relations Among Modes of Convergence

The three modes of convergence satisfy

\[
X_n\xrightarrow{a.s.}X
\quad\Longrightarrow\quad
X_n\xrightarrow{P}X
\quad\Longrightarrow\quad
X_n\xrightarrow{d}X.
\]

In general, neither converse holds. However, if the limit is a constant $c$, then

\[
X_n\xrightarrow{d}c
\quad\Longleftrightarrow\quad
X_n\xrightarrow{P}c.
\]

---

## 5. Almost Sure Convergence, Borel-Cantelli, and Subsequences

### 5.1 Borel-Cantelli Criterion

!!! success "Theorem 5.1 (Complete Convergence Implies Almost Sure Convergence)"

    If, for every $\varepsilon>0$,

    \[
    \sum_{n=1}^{\infty}
    P(|X_n-X|>\varepsilon)<\infty,
    \]

    then

    \[
    X_n\xrightarrow{a.s.}X.
    \]

??? proof "Proof of Theorem 5.1 (click to expand)"

    Fix $\varepsilon>0$ and let

    \[
    A_n^{(\varepsilon)}
    =\{|X_n-X|>\varepsilon\}.
    \]

    By the first Borel-Cantelli lemma,

    \[
    \sum_{n=1}^{\infty}P(A_n^{(\varepsilon)})<\infty
    \quad\Longrightarrow\quad
    P(A_n^{(\varepsilon)}\ \text{i.o.})=0.
    \]

    Taking a countable intersection over all positive rational $\varepsilon$ gives $X_n\to X$ almost surely. $\square$

!!! example "Example 5.2 (A Fourth-Moment Proof of the Bernoulli Strong Law)"

    For the Bernoulli sum $S_n$,

    \[
    E\left|\frac{S_n}{n}-p\right|^4
    =O(n^{-2}).
    \]

    Hence, by Markov's inequality,

    \[
    P\left(\left|\frac{S_n}{n}-p\right|>\varepsilon\right)
    \leq
    \frac{E|S_n/n-p|^4}{\varepsilon^4}
    =O(n^{-2}).
    \]

    Since $\sum_n n^{-2}<\infty$, Theorem 5.1 yields

    \[
    \frac{S_n}{n}\xrightarrow{a.s.}p.
    \]

### 5.2 Subsequence Characterization of Convergence in Probability

!!! success "Theorem 5.3 (Riesz Subsequence Principle)"

    $X_n\xrightarrow{P}X$ if and only if every subsequence $\{X_{n_k}\}$ has a further subsequence $\{X_{n_{k_j}}\}$ such that

    \[
    X_{n_{k_j}}\xrightarrow{a.s.}X.
    \]

??? proof "Construction for the Sufficient Direction (click to expand)"

    Assume $X_n\xrightarrow{P}X$. From any subsequence, recursively choose $n_{k_j}$ such that

    \[
    P\left(|X_{n_{k_j}}-X|>2^{-j}\right)<2^{-j}.
    \]

    Consequently,

    \[
    \sum_{j=1}^{\infty}
    P\left(|X_{n_{k_j}}-X|>2^{-j}\right)<\infty.
    \]

    The Borel-Cantelli lemma gives $X_{n_{k_j}}\to X$ almost surely. The reverse direction follows directly by contradiction. $\square$

### 5.3 Skorohod Representation Theorem

!!! success "Theorem 5.4 (Skorohod Representation Theorem)"

    If

    \[
    X_n\xrightarrow{d}X,
    \]

    then there exists another probability space carrying random variables $X_n'$ and $X'$ such that

    \[
    X_n'\overset{d}{=}X_n,
    \qquad
    X'\overset{d}{=}X,
    \]

    and

    \[
    X_n'\xrightarrow{a.s.}X'.
    \]

!!! warning "The Probability Space Has Changed"

    The Skorohod representation theorem does not say that the original $X_n$ converge almost surely on their original probability space. It only constructs versions with the same marginal distributions on a new common probability space where almost sure convergence holds.

---

## 6. Infinitely Divisible Distributions

### 6.1 Definition and Equivalent Characterizations

!!! info "Definition 6.1 (Infinitely Divisible Distribution)"

    If, for every positive integer $n$, there exist independent and identically distributed random variables

    \[
    X_{n,1},\ldots,X_{n,n}
    \]

    such that

    \[
    X\overset{d}{=}
    X_{n,1}+\cdots+X_{n,n},
    \]

    then the distribution of $X$ is called infinitely divisible.

If $F$ is the distribution function of $X$ and $F_n$ is the distribution function of $X_{n,1}$, then the definition is equivalent to

\[
F=F_n^{*n},
\]

where $*$ denotes convolution. If $\varphi$ and $\varphi_n$ are the characteristic functions of $X$ and $X_{n,1}$, respectively, then equivalently,

\[
\varphi(t)=\varphi_n(t)^n,
\qquad t\in\mathbb{R}.
\]

### 6.2 Standard Examples

!!! example "Example 6.2 (Degenerate Distributions)"

    If $X\equiv c$, then

    \[
    X=\frac{c}{n}+\cdots+\frac{c}{n},
    \]

    so every degenerate distribution is infinitely divisible. Its characteristic function satisfies

    \[
    e^{itc}
    =\left(e^{itc/n}\right)^n.
    \]

!!! example "Example 6.3 (Normal Distribution)"

    If $X\sim N(\mu,\sigma^2)$, then

    \[
    \varphi(t)
    =\exp\left(i\mu t-\frac{\sigma^2t^2}{2}\right)
    =\left[
    \exp\left(
    i\frac{\mu}{n}t-
    \frac{\sigma^2}{2n}t^2
    \right)
    \right]^n.
    \]

    Thus one may take

    \[
    X_{n,k}\sim N\left(\frac{\mu}{n},\frac{\sigma^2}{n}\right).
    \]

!!! example "Example 6.4 (Poisson Distribution)"

    If $X\sim\operatorname{Poisson}(\lambda)$, then

    \[
    \varphi(t)
    =\exp\{\lambda(e^{it}-1)\}
    =\left[
    \exp\left\{\frac{\lambda}{n}(e^{it}-1)\right\}
    \right]^n.
    \]

    Thus one may take

    \[
    X_{n,k}\sim\operatorname{Poisson}\left(\frac{\lambda}{n}\right).
    \]

!!! example "Example 6.5 (Cauchy Distribution)"

    The characteristic function of the standard Cauchy distribution is

    \[
    \varphi(t)=e^{-|t|}
    =\left(e^{-|t|/n}\right)^n.
    \]

    Hence the standard Cauchy distribution is infinitely divisible.

### 6.3 Examples That Are Not Infinitely Divisible

The characteristic function of the Rademacher distribution is

\[
\varphi(t)=\cos t,
\]

whereas the characteristic function of $\operatorname{Uniform}(-1,1)$ is

\[
\varphi(t)=\frac{\sin t}{t}.
\]

Both characteristic functions have zeros. The next section proves that the characteristic function of an infinitely divisible distribution cannot vanish. Therefore, neither distribution is infinitely divisible.

---

## 7. Basic Properties of Infinitely Divisible Distributions

### 7.1 Convolution and Symmetrization

!!! success "Proposition 7.1 (Closure Properties)"

    If $f$ and $g$ are characteristic functions of infinitely divisible distributions, then the following are also characteristic functions of infinitely divisible distributions:

    \[
    f(t)g(t),
    \]

    \[
    f(-t)=\overline{f(t)},
    \]

    and

    \[
    |f(t)|^2=f(t)f(-t).
    \]

These correspond, respectively, to the sum of independent infinitely divisible random variables, reflection, and the difference of two independent identically distributed random variables.

### 7.2 The Characteristic Function Has No Zeros

!!! success "Theorem 7.2 (Nonvanishing Property)"

    If $f$ is the characteristic function of an infinitely divisible distribution, then

    \[
    f(t)\neq0,
    \qquad \forall t\in\mathbb{R}.
    \]

??? proof "Proof of Theorem 7.2 (click to expand)"

    For every $n$, there exists a characteristic function $f_n$ such that

    \[
    f(t)=f_n(t)^n.
    \]

    Suppose that $f(t_0)=0$ for some $t_0$. Since $f(0)=1$ and $f$ is continuous, choose the positive zero $t_0$ closest to the origin. Then $f$ has no zeros on $[-\delta,\delta]$ for some $\delta>0$.

    From

    \[
    |f_n(t)|=|f(t)|^{1/n},
    \]

    it follows that, for fixed $t$ with $f(t)\neq0$,

    \[
    |f_n(t)|\longrightarrow1.
    \]

    Every characteristic function $h$ satisfies

    \[
    0\leq1-|h(2t)|^2
    \leq4\{1-|h(t)|^2\}.
    \]

    Starting from $|f_n(t)|\to1$ near the origin and repeatedly applying this inequality propagates nonvanishing from a neighborhood of the origin to the entire real line. This implies $f(t_0)\neq0$, a contradiction. $\square$

### 7.3 A Bounded Infinitely Divisible Distribution Is Degenerate

!!! success "Theorem 7.3 (Bounded Infinitely Divisible Distributions Are Degenerate)"

    Suppose that $X$ is infinitely divisible and that there exists $M<\infty$ such that

    \[
    |X|\leq M
    \quad\text{almost surely},
    \]

    Then $X$ must be degenerate.

??? proof "Proof of Theorem 7.3 (click to expand)"

    For every $n$, write

    \[
    X\overset{d}{=}
    X_{n,1}+\cdots+X_{n,n},
    \]

    where $X_{n,1},\ldots,X_{n,n}$ are independent and identically distributed. Since the sum lies almost surely in $[-M,M]$ and the summands are independent, the essential oscillation of each summand is at most $2M/n$. After translating the summands while keeping the total translation unchanged, we may arrange that

    \[
    |X_{n,k}|\leq\frac{M}{n}
    \quad\text{almost surely}.
    \]

    Therefore,

    \[
    \operatorname{Var}(X_{n,1})
    \leq E(X_{n,1}^2)
    \leq\frac{M^2}{n^2}.
    \]

    By independence,

    \[
    \operatorname{Var}(X)
    =n\operatorname{Var}(X_{n,1})
    \leq\frac{M^2}{n}.
    \]

    Letting $n\to\infty$ gives $\operatorname{Var}(X)=0$, so $X$ is degenerate. $\square$

### 7.4 Closure Under Weak Convergence

!!! success "Theorem 7.4 (Weak Closure of Infinite Divisibility)"

    If every $F_m$ is infinitely divisible and

    \[
    F_m\xrightarrow{w}F,
    \]

    then $F$ is also infinitely divisible.

??? proof "Proof Idea for Theorem 7.4 (click to expand)"

    Let $f_m$ and $f$ be the characteristic functions of $F_m$ and $F$, respectively. By Lévy's continuity theorem,

    \[
    f_m(t)\longrightarrow f(t).
    \]

    Fix $n$. Since $F_m$ is infinitely divisible, there exists a characteristic function $f_{n,m}$ such that

    \[
    f_{n,m}(t)^n=f_m(t).
    \]

    Tightness permits the selection of a subsequence along which $f_{n,m}$ converges to a characteristic function $f_n$. Passing to the limit gives

    \[
    f_n(t)^n=f(t).
    \]

    Since this holds for every $n$, $F$ is infinitely divisible. The nonvanishing property guarantees that a continuous choice of the $n$th root can be made consistently. $\square$

---

## 8. Lévy-Khintchine Representation

!!! success "Theorem 8.1 (Lévy-Khintchine Formula)"

    A characteristic function $\varphi$ corresponds to an infinitely divisible distribution if and only if there exist $\gamma\in\mathbb{R}$, a Lévy measure $\nu$ satisfying

    \[
    \nu(\{0\})=0,
    \qquad
    \int_{\mathbb{R}}(1\wedge x^2)\,\nu(dx)<\infty,
    \]

    and some $\sigma^2\geq0$ such that

    \[
    \varphi(t)
    =\exp\left\{
    i\gamma t-\frac{\sigma^2t^2}{2}
    +\int_{\mathbb{R}}
    \left(e^{itx}-1-
    \frac{itx}{1+x^2}\right)
    \nu(dx)
    \right\}.
    \]

Equivalently, one may absorb the Gaussian component into a bounded nondecreasing function $G$ and write the Kolmogorov canonical form

\[
\varphi(t)
=\exp\left\{
i\gamma t+
\int_{-\infty}^{\infty}
\left(e^{itx}-1-
\frac{itx}{1+x^2}\right)
\frac{1+x^2}{x^2}\,dG(x)
\right\}.
\]

Here $G$ may be chosen left-continuous, bounded, and nondecreasing, with

\[
G(-\infty)=0,
\qquad
G(+\infty)<\infty.
\]

The three terms describe drift, Gaussian fluctuation, and jump behavior. Once the truncation function is fixed, $\gamma$, $\sigma^2$, and $\nu$ are unique.

---

## 9. Limits of Infinitesimal Triangular Arrays

For each $n$, let

\[
\xi_{n,1},\ldots,\xi_{n,k_n}
\]

be mutually independent, and define

\[
S_n=\sum_{k=1}^{k_n}\xi_{n,k}.
\]

!!! info "Definition 9.1 (Infinitesimal Triangular Array)"

    If, for every $\varepsilon>0$,

    \[
    \max_{1\leq k\leq k_n}
    P(|\xi_{n,k}|>\varepsilon)
    \longrightarrow0.
    \]

    then the triangular array is called infinitesimal.

!!! success "Theorem 9.2 (Limits of Infinitesimal Arrays Are Infinitely Divisible)"

    If the triangular array is infinitesimal and

    \[
    S_n\xrightarrow{d}S,
    \]

    then the distribution of $S$ is infinitely divisible.

??? proof "Proof Idea for Theorem 9.2 (click to expand)"

    Fix a positive integer $m$. Partition each row of independent variables into $m$ consecutive blocks so that the blocks contribute as equally as possible to the total characteristic exponent, and sum the variables within each block. The infinitesimal condition ensures that the contribution of any single variable tends to zero, so the block sums have the same limiting distribution along a suitable subsequence.

    Denote the block sums by $Y_{n,1},\ldots,Y_{n,m}$. They are independent and satisfy

    \[
    S_n=Y_{n,1}+\cdots+Y_{n,m}+o_P(1).
    \]

    Passing to a subsequential limit produces independent and identically distributed random variables $Y_1,\ldots,Y_m$ such that

    \[
    S\overset{d}{=}Y_1+\cdots+Y_m.
    \]

    Since $m$ is arbitrary, the distribution of $S$ is infinitely divisible. $\square$

!!! note "Role in Limit Theory"

    Infinitely divisible distributions are precisely the possible limits of infinitesimal triangular arrays of independent random variables. Normal, Poisson, and more general stable distributions all belong to this class.

---

## 10. Chapter Summary

- Under different normalizations, Bernoulli sums yield the law of large numbers, the central limit theorem, and the Poisson limit.
- Almost sure convergence implies convergence in probability, which in turn implies convergence in distribution; the converses generally fail.
- The continuous mapping theorem, Slutsky's theorem, and Lévy's continuity theorem are central tools for weak convergence.
- The Borel-Cantelli lemma converts summability of probabilities into almost sure convergence, and convergence in probability admits a characterization through almost surely convergent subsequences.
- An infinitely divisible distribution can be decomposed into a sum of any number of independent identically distributed random variables, and its characteristic function has no zeros.
- Infinitely divisible distributions are closed under convolution and weak convergence and are completely characterized by the Lévy-Khintchine formula.
- Every weak limit of an infinitesimal independent triangular array is infinitely divisible.
