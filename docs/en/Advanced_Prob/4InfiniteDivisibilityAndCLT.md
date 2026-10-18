# Lecture 4: Infinite Divisibility, Infinitesimal Arrays, and the Central Limit Theorem

This lecture first completes the Lévy–Khintchine representation of infinitely divisible distributions and illustrates its structure with the normal and Poisson distributions. It then introduces infinitesimal triangular arrays and explains why infinitely divisible distributions are precisely the possible limits of their row sums. Finally, it turns to the central limit theorem and presents the characteristic-function method, the Lindeberg replacement technique, and Stein's method.

---

## 1. Infinitely Divisible Distributions

### 1.1 Definition

!!! info "Definition 1.1 (Infinitely divisible distribution)"

    Let the distribution function and characteristic function of a random variable $X$ be $F$ and $f$, respectively. If, for every positive integer $n$, there exist independent and identically distributed random variables

    \[
    X_{n,1},\ldots,X_{n,n}
    \]

    such that

    \[
    X\overset{d}{=}X_{n,1}+\cdots+X_{n,n},
    \]

    then $F$ is called an **infinitely divisible distribution function**, and $f$ is called an **infinitely divisible characteristic function**, abbreviated as i.d.d.f. and i.d.c.f.

In terms of characteristic functions, this is equivalent to the following statement: for every $n\geq1$, there exists a characteristic function $f_n$ such that

\[
f(t)=\bigl(f_n(t)\bigr)^n,
\qquad t\in\mathbb{R}.
\]

!!! note "Basic examples"

    Normal, Poisson, and compound Poisson distributions are all infinitely divisible.

### 1.2 The Lévy–Khintchine Representation

!!! success "Theorem 1.2 (Canonical Lévy–Khintchine representation)"

    A characteristic function $f$ is infinitely divisible if and only if there exist a constant $\gamma\in\mathbb{R}$ and a bounded, right-continuous, nondecreasing function $G$ satisfying $G(-\infty)=0$ such that

    \[
    f(t)=\exp\{\psi(t)\},
    \]

    where

    \[
    \psi(t)
    =i\gamma t
    +\int_{-\infty}^{\infty}
    \left(
    e^{itx}-1-\frac{itx}{1+x^2}
    \right)
    \frac{1+x^2}{x^2}\,dG(x).
    \]

    At $x=0$, the integrand is interpreted by continuous extension as

    \[
    \lim_{x\to0}
    \left(
    e^{itx}-1-\frac{itx}{1+x^2}
    \right)
    \frac{1+x^2}{x^2}
    =-\frac{t^2}{2}.
    \]

The ordered pair $(\gamma,G)$ is called a set of **canonical parameters** of the infinitely divisible distribution. The jump of $G$ at the origin produces the Gaussian component, whereas its increments away from the origin produce the jump component.

??? proof "Constructive proof of sufficiency (click to expand)"

    Suppose that $\psi$ has the displayed form. The linear term $i\gamma t$ is the characteristic exponent of a constant shift. If $G$ has a jump at the origin of size

    \[
    \sigma^2=G(0+)-G(0-),
    \]

    its contribution to the exponent is

    \[
    -\frac{\sigma^2t^2}{2},
    \]

    which is exactly the characteristic exponent of the centered normal distribution $N(0,\sigma^2)$.

    Next, consider a finite interval $[a,b]$ not containing the origin. Take a partition

    \[
    a=x_0<x_1<\cdots<x_m=b
    \]

    and points $\xi_j\in(x_{j-1},x_j]$. The $j$th term in the corresponding Riemann–Stieltjes sum can be written as

    \[
    \lambda_j\left(e^{it\xi_j}-1\right)-itc_j,
    \]

    where

    \[
    \lambda_j
    =\frac{1+\xi_j^2}{\xi_j^2}
    \bigl(G(x_j)-G(x_{j-1})\bigr)\geq0,
    \]

    and

    \[
    c_j
    =\frac{\xi_j}{1+\xi_j^2}\lambda_j.
    \]

    Hence the exponential of this term factors as

    \[
    \exp\left\{\lambda_j(e^{it\xi_j}-1)\right\}
    \exp\{-itc_j\}.
    \]

    The first factor is the characteristic function of a compound Poisson distribution with jump size $\xi_j$ and intensity $\lambda_j$, while the second is the characteristic function of a constant shift. Their product is therefore infinitely divisible.

    Let the mesh of the partition tend to zero and then let $a\downarrow-\infty$ and $b\uparrow\infty$. By the continuity theorem for characteristic functions, the limit $e^{\psi(t)}$ is again a characteristic function. Moreover, for every $n$,

    \[
    e^{\psi(t)}
    =\left(e^{\psi(t)/n}\right)^n,
    \]

    and $e^{\psi(t)/n}$ is given by the parameters $(\gamma/n,G/n)$. Thus $f$ is infinitely divisible. $\square$

??? proof "Idea of the necessity proof (textbook supplement; click to expand)"

    If $f$ is infinitely divisible, then $f=f_n^n$ for every $n$. Regard the corresponding $n$ i.i.d. variables as one row of an infinitesimal triangular array. Analyze separately the truncated means, truncated variances, and distributions of large jumps in each row, and extract convergent subsequences. This yields limiting parameters $\gamma$ and $G$. A uniform small-term expansion of

    \[
    n\log f_n(t)
    \]

    then gives precisely the exponent $\psi(t)$ above. Uniqueness of the canonical parameters ensures that all subsequences yield the same pair $(\gamma,G)$, which establishes the Lévy–Khintchine representation for the full sequence. $\square$

### 1.3 The Normal Distribution

If $X\sim N(\mu,\sigma^2)$, then

\[
f(t)=\exp\left\{i\mu t-\frac{\sigma^2t^2}{2}\right\}.
\]

In the canonical representation, one may take

\[
\gamma=\mu
\]

and let $G$ have a single jump of size $\sigma^2$ at the origin; that is,

\[
G(x)=
\begin{cases}
0, & x<0,\\
\sigma^2, & x\geq0.
\end{cases}
\]

Thus, the infinite divisibility of a normal distribution comes from its Gaussian component.

### 1.4 The Poisson Distribution

If $X\sim\operatorname{Pois}(\lambda)$, then

\[
f(t)=\exp\left\{\lambda(e^{it}-1)\right\}.
\]

One may take

\[
G(x)=
\begin{cases}
0, & x<1,\\
\lambda/2, & x\geq1,
\end{cases}
\qquad
\gamma=\frac{\lambda}{2}.
\]

Indeed, the jump of $G$ at $1$ gives

\[
2\cdot\frac{\lambda}{2}
\left(e^{it}-1-\frac{it}{2}\right),
\]

which, when combined with $i\gamma t=i\lambda t/2$, yields $\lambda(e^{it}-1)$.

---

## 2. Infinitesimal Triangular Arrays and Limit Distributions

### 2.1 Infinitesimal Triangular Arrays

Consider a triangular array

\[
\{X_{n,k}:1\leq k\leq k_n,\ n\geq1\}
\]

and suppose that the random variables in each row are mutually independent. Write the sum of the $n$th row as

\[
S_n=\sum_{k=1}^{k_n}X_{n,k}.
\]

!!! info "Definition 2.1 (Infinitesimal array)"

    The triangular array is called **infinitesimal** if, for every $\varepsilon>0$,

    \[
    \max_{1\leq k\leq k_n}
    P\bigl(|X_{n,k}|>\varepsilon\bigr)
    \longrightarrow0.
    \]

This condition says that, as $n\to\infty$, no individual term in a row has a macroscopic effect on the sum.

### 2.2 Every Limit Is Infinitely Divisible

!!! success "Theorem 2.2 (Limit of an infinitesimal array)"

    If the variables in each row are independent, the triangular array is infinitesimal, and

    \[
    S_n\overset{d}{\longrightarrow}S,
    \]

    then the distribution of $S$ is infinitely divisible.

??? proof "Proof idea (click to expand)"

    Fix a positive integer $m$. Infinitesimality allows the $n$th row to be divided consecutively into $m$ independent blocks, each contributing approximately $1/m$ of the total characteristic exponent. Denote the block sums by

    \[
    Y_{n,1}^{(m)},\ldots,Y_{n,m}^{(m)}.
    \]

    A tightness and diagonal-subsequence argument makes all block sums converge along a common subsequence. After adjusting negligible centering constants, the $m$ limiting blocks can be chosen to have the same distribution. Hence

    \[
    S\overset{d}{=}Y_1^{(m)}+\cdots+Y_m^{(m)},
    \]

    where $Y_1^{(m)},\ldots,Y_m^{(m)}$ are independent and identically distributed. Since $m$ is arbitrary, the distribution of $S$ is infinitely divisible. $\square$

### 2.3 Converse

!!! success "Theorem 2.3 (Array representation of an infinitely divisible distribution)"

    If $S$ has an infinitely divisible distribution, then there exists a row-wise independent infinitesimal triangular array whose row sums converge in distribution to $S$. In fact, every row sum can be chosen to have the same distribution as $S$.

??? proof "Proof (click to expand)"

    By infinite divisibility, for every $n$ there exist independent and identically distributed variables

    \[
    X_{n,1},\ldots,X_{n,n}
    \]

    such that

    \[
    S\overset{d}{=}\sum_{k=1}^nX_{n,k}.
    \]

    If $f$ is the characteristic function of $S$, choose its $n$th convolution root $f_n$ so that

    \[
    \bigl(f_n(t)\bigr)^n=f(t).
    \]

    Since $f_n(t)\to1$, Lévy's continuity theorem gives

    \[
    X_{n,1}\overset{d}{\longrightarrow}0.
    \]

    When the limit is constant, convergence in distribution is equivalent to convergence in probability. Therefore,

    \[
    P\bigl(|X_{n,1}|>\varepsilon\bigr)\longrightarrow0.
    \]

    The variables in the same row have the same distribution, so

    \[
    \max_{1\leq k\leq n}
    P\bigl(|X_{n,k}|>\varepsilon\bigr)
    =P\bigl(|X_{n,1}|>\varepsilon\bigr)
    \longrightarrow0.
    \]

    Thus the array is infinitesimal, while every row sum has the same distribution as $S$. $\square$

### 2.4 Three Limits for Bernoulli Sums

Let $\xi_1,\xi_2,\ldots$ be i.i.d. random variables satisfying

\[
P(\xi_i=1)=p,
\qquad
P(\xi_i=0)=q=1-p.
\]

Write

\[
T_n=\sum_{i=1}^n\xi_i.
\]

#### Law-of-Large-Numbers Scale

By the weak law of large numbers,

\[
\frac{T_n}{n}\overset{P}{\longrightarrow}p.
\]

Regard $X_{n,i}=\xi_i/n$ as a triangular-array element. For every fixed $\varepsilon>0$, when $n>1/\varepsilon$,

\[
P\bigl(|X_{n,i}|>\varepsilon\bigr)=0,
\]

so the array is infinitesimal.

#### Central-Limit-Theorem Scale

We have

\[
\frac{T_n-np}{\sqrt{npq}}
=\sum_{i=1}^n\frac{\xi_i-p}{\sqrt{npq}}
\overset{d}{\longrightarrow}N(0,1).
\]

The corresponding array elements

\[
X_{n,i}=\frac{\xi_i-p}{\sqrt{npq}}
\]

satisfy $|X_{n,i}|\leq1/\sqrt{npq}\to0$ and hence also form an infinitesimal array.

#### Poisson Rare-Event Scale

Suppose that $\xi_{n,1},\ldots,\xi_{n,n}$ are i.i.d. within the $n$th row and satisfy

\[
P(\xi_{n,k}=1)=p_n,
\qquad
P(\xi_{n,k}=0)=1-p_n.
\]

If

\[
np_n\longrightarrow\lambda\in(0,\infty),
\]

then $p_n\to0$, the array is infinitesimal, and

\[
\sum_{k=1}^n\xi_{n,k}
\overset{d}{\longrightarrow}
\operatorname{Pois}(\lambda).
\]

The characteristic-function proof is

\[
\mathbb{E}\left[
\exp\left\{
it\sum_{k=1}^{n}\xi_{n,k}
\right\}
\right]
=
\left(
1+p_n\left(e^{it}-1\right)
\right)^n
\longrightarrow
\exp\left\{
\lambda\left(e^{it}-1\right)
\right\}.
\]

---

## 3. Lévy's Continuity Theorem

Let the characteristic functions of $Z_n$ and $Z$ be $\varphi_n$ and $\varphi$, respectively.

!!! success "Theorem 3.1 (Lévy's continuity theorem)"

    If

    \[
    Z_n\overset{d}{\longrightarrow}Z,
    \]

    then, for every $t\in\mathbb{R}$,

    \[
    \varphi_n(t)\longrightarrow\varphi(t).
    \]

    Conversely, if $\varphi_n(t)$ converges for every $t$ to a function $\varphi(t)$ that is continuous at $0$, then $\varphi$ is the characteristic function of some random variable $Z$, and

    \[
    Z_n\overset{d}{\longrightarrow}Z.
    \]

??? proof "Forward implication (click to expand)"

    For fixed $t$, the function $x\mapsto e^{itx}$ is bounded and continuous. By the Portmanteau theorem for convergence in distribution,

    \[
    Ee^{itZ_n}\longrightarrow Ee^{itZ}.
    \]

    Thus $\varphi_n(t)\to\varphi(t)$. $\square$

??? proof "Idea of the converse proof (click to expand)"

    Continuity of $\varphi$ at the origin, together with $\varphi_n(0)=1$, implies tightness of the distributions of $\{Z_n\}$. By Helly's selection theorem, every subsequence has a further subsequence converging in distribution to some $Z$. The forward implication shows that the limiting characteristic function of this further subsequence must be $\varphi$. Uniqueness of characteristic functions implies that all convergent subsequences have the same limiting distribution. Hence the whole sequence converges in distribution, and $\varphi$ is the characteristic function of the limit. $\square$

---

## 4. The i.i.d. Central Limit Theorem

Let $X_1,X_2,\ldots$ be i.i.d. and satisfy

\[
EX_1=\mu,
\qquad
\operatorname{Var}(X_1)=\sigma^2\in(0,\infty).
\]

Define

\[
Z_n=\frac{\sum_{i=1}^n(X_i-\mu)}{\sigma\sqrt{n}}.
\]

!!! success "Theorem 4.1 (Lindeberg–Lévy central limit theorem)"

    Under the preceding assumptions,

    \[
    Z_n\overset{d}{\longrightarrow}N(0,1).
    \]

### 4.1 Characteristic-Function Proof

??? proof "Proof (click to expand)"

    Let

    \[
    Y_i=\frac{X_i-\mu}{\sigma}.
    \]

    Then $EY_i=0$ and $EY_i^2=1$. Let $g$ be the characteristic function of $Y_1$. Since the second moment is finite, $g$ has the expansion

    \[
    g(u)=1-\frac{u^2}{2}+o(u^2),
    \qquad u\to0.
    \]

    By independence, the characteristic function of $Z_n=n^{-1/2}\sum_{i=1}^nY_i$ is

    \[
    \varphi_n(t)
    =\left[g\left(\frac{t}{\sqrt n}\right)\right]^n.
    \]

    Hence

    \[
    g\left(\frac{t}{\sqrt n}\right)
    =1-\frac{t^2}{2n}+o\left(\frac1n\right),
    \]

    and therefore

    \[
    \varphi_n(t)
    \longrightarrow e^{-t^2/2}.
    \]

    The right-hand side is the characteristic function of the standard normal distribution. Lévy's continuity theorem now yields

    \[
    Z_n\overset{d}{\longrightarrow}N(0,1).
    \]

    $\square$

### 4.2 The Lindeberg Replacement Technique

After standardization, suppose that $EX_i=0$ and $EX_i^2=1$. Let $Y_i\sim N(0,1)$ be i.i.d., and suppose that the two sequences are independent of one another. For a sufficiently smooth bounded function $f$, define

\[
W_{n,r}
=\frac{Y_1+\cdots+Y_{r-1}+X_{r+1}+\cdots+X_n}{\sqrt n}.
\]

The telescoping decomposition gives

\[
\begin{aligned}
&E f\left(\frac{X_1+\cdots+X_n}{\sqrt n}\right)
-E f\left(\frac{Y_1+\cdots+Y_n}{\sqrt n}\right)\\
&=\sum_{r=1}^n
E\left[
f\left(W_{n,r}+\frac{X_r}{\sqrt n}\right)
-f\left(W_{n,r}+\frac{Y_r}{\sqrt n}\right)
\right].
\end{aligned}
\]

??? proof "Lindeberg replacement argument (click to expand)"

    Take a second-order Taylor expansion at $W_{n,r}$:

    \[
    f(W_{n,r}+u)
    =f(W_{n,r})+uf'(W_{n,r})
    +\frac{u^2}{2}f''(W_{n,r})+R_{n,r}(u).
    \]

    Since $X_r$ and $Y_r$ are independent of $W_{n,r}$ and

    \[
    EX_r=EY_r=0,
    \qquad
    EX_r^2=EY_r^2=1,
    \]

    the zeroth-, first-, and second-order terms cancel in every difference, leaving only the Taylor remainders. Truncate $X_r$ first, and then use $EX_r^2<\infty$ and uniform continuity of $f''$ to obtain

    \[
    \sum_{r=1}^n
    E\left|R_{n,r}\left(\frac{X_r}{\sqrt n}\right)
    -R_{n,r}\left(\frac{Y_r}{\sqrt n}\right)\right|
    \longrightarrow0.
    \]

    Consequently,

    \[
    E f\left(\frac{X_1+\cdots+X_n}{\sqrt n}\right)
    -E f\left(\frac{Y_1+\cdots+Y_n}{\sqrt n}\right)
    \longrightarrow0.
    \]

    Since

    \[
    \frac{Y_1+\cdots+Y_n}{\sqrt n}\sim N(0,1),
    \]

    expectations converge for a convergence-determining class of smooth test functions. Smooth approximation extends the conclusion to every $f\in C_b(\mathbb R)$, proving the central limit theorem. $\square$

### 4.3 The Lindeberg Condition in the Textbook

The preceding replacement idea extends to independent triangular arrays whose variables need not be identically distributed. Suppose that $EX_{n,k}=0$ and define

\[
B_n=\sum_{k=1}^{k_n}EX_{n,k}^2.
\]

!!! info "Definition 4.2 (Lindeberg condition)"

    The array satisfies the Lindeberg condition if, for every $\varepsilon>0$,

    \[
    \frac1{B_n}
    \sum_{k=1}^{k_n}
    E\left[
    X_{n,k}^2
    \mathbf 1_{\{|X_{n,k}|>\varepsilon\sqrt{B_n}\}}
    \right]
    \longrightarrow0.
    \]

!!! success "Theorem 4.3 (Sufficiency in the Lindeberg–Feller theorem)"

    If the variables in each row are independent and centered, $B_n>0$, and the Lindeberg condition holds, then

    \[
    \frac{\sum_{k=1}^{k_n}X_{n,k}}{\sqrt{B_n}}
    \overset{d}{\longrightarrow}N(0,1).
    \]

The Lindeberg condition excludes the possibility that one large jump dominates the sum and implies the Feller infinitesimal condition

\[
\max_{1\leq k\leq k_n}
\frac{EX_{n,k}^2}{B_n}
\longrightarrow0.
\]

---

## 5. Stein's Method

### 5.1 The Normal Integration-by-Parts Identity

Let $Z\sim N(0,1)$. For every suitably smooth and integrable function $f$,

\[
E\bigl[Zf(Z)\bigr]=E\bigl[f'(Z)\bigr].
\]

??? proof "Proof (click to expand)"

    Write the standard normal density as

    \[
    \phi(z)=\frac1{\sqrt{2\pi}}e^{-z^2/2}.
    \]

    Since $\phi'(z)=-z\phi(z)$, integration by parts gives

    \[
    \begin{aligned}
    E[Zf(Z)]
    &=\int_{-\infty}^{\infty}zf(z)\phi(z)\,dz\\
    &=-\int_{-\infty}^{\infty}f(z)\,d\phi(z)\\
    &=\int_{-\infty}^{\infty}f'(z)\phi(z)\,dz\\
    &=E[f'(Z)].
    \end{aligned}
    \]

    The boundary term vanishes. $\square$

### 5.2 Stein's Characterization

!!! success "Theorem 5.1 (Stein's characterization of the standard normal distribution)"

    A random variable $W$ has distribution $N(0,1)$ if and only if, for a sufficiently rich class of absolutely continuous functions $f$,

    \[
    E\bigl[f'(W)-Wf(W)\bigr]=0.
    \]

??? proof "Proof of the converse (click to expand)"

    Fix $z\in\mathbb R$ and define

    \[
    h_z(x)=\mathbf 1_{(-\infty,z]}(x)-\Phi(z).
    \]

    Consider Stein's equation

    \[
    f_z'(x)-xf_z(x)=h_z(x).
    \]

    It has a bounded absolutely continuous solution. If $W$ satisfies Stein's identity, then

    \[
    0=E\bigl[f_z'(W)-Wf_z(W)\bigr]
    =E h_z(W).
    \]

    Hence

    \[
    P(W\leq z)=\Phi(z)
    \]

    for every $z$, so $W\sim N(0,1)$. $\square$

### 5.3 Proving the Central Limit Theorem with Stein's Identity

Again suppose that $EX_i=0$ and $EX_i^2=1$, and define

\[
Z_n=\frac1{\sqrt n}\sum_{i=1}^nX_i,
\qquad
Z_{n,i}=Z_n-\frac{X_i}{\sqrt n}.
\]

Observe that $Z_{n,i}$ and $X_i$ are independent. For a smooth function $f$,

\[
\begin{aligned}
E[Z_nf(Z_n)]
&=\frac1{\sqrt n}\sum_{i=1}^nE[X_if(Z_n)]\\
&=\frac1{\sqrt n}\sum_{i=1}^n
E\left[X_i\bigl(f(Z_n)-f(Z_{n,i})\bigr)\right].
\end{aligned}
\]

The second equality uses

\[
E[X_if(Z_{n,i})]=EX_i\,Ef(Z_{n,i})=0.
\]

Expand the difference to first order:

\[
f(Z_n)-f(Z_{n,i})
=\frac{X_i}{\sqrt n}f'(Z_{n,i})+R_{n,i}.
\]

It follows that

\[
E[Z_nf(Z_n)]
=\frac1n\sum_{i=1}^nE[f'(Z_{n,i})]
+\frac1{\sqrt n}\sum_{i=1}^nE[X_iR_{n,i}].
\]

Using $EX_i^2=1$, a truncation argument, and uniform continuity of $f'$, one obtains

\[
\frac1{\sqrt n}\sum_{i=1}^nE[X_iR_{n,i}]
\longrightarrow0,
\]

and

\[
\frac1n\sum_{i=1}^nE[f'(Z_{n,i})]
-E[f'(Z_n)]
\longrightarrow0.
\]

Therefore,

\[
E\bigl[f'(Z_n)-Z_nf(Z_n)\bigr]
\longrightarrow0.
\]

Thus $Z_n$ asymptotically satisfies Stein's identity for the standard normal distribution. Stein's equation and tightness then yield

\[
Z_n\overset{d}{\longrightarrow}N(0,1).
\]

---

## 6. Logical Structure of the Lecture

!!! note "Summary of the main line"

    **Infinite divisibility:** a distribution can be decomposed, for every $n$, into the sum of $n$ i.i.d. random variables.

    **Lévy–Khintchine representation:** translation, a Gaussian component, and a jump component completely characterize infinitely divisible characteristic functions.

    **Infinitesimal triangular arrays:** every weak limit of their row sums is infinitely divisible, and every infinitely divisible distribution can be obtained from such row sums.

    **Lévy's continuity theorem:** pointwise convergence of characteristic functions is converted into convergence in distribution.

    **Central limit theorem:** it can be proved by characteristic functions, Lindeberg replacement, or Stein's identity; the three approaches emphasize product structure, moment matching, and distributional characterization, respectively.
