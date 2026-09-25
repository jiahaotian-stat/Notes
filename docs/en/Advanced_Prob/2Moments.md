# Chapter 2: Independence, Numerical Characteristics, Characteristic Functions, and Convergence in Distribution

This chapter first discusses the independence of random variables and its preservation under measurable transformations, then introduces expectation, variance, moments, and Chebyshev's inequality. On this basis, we introduce moment generating functions and characteristic functions, and use characteristic functions to study sums of independent random variables and the uniqueness of distributions. Finally, we introduce convergence in probability, convergence in distribution, limit theorems for sums of Bernoulli random variables, and the central limit theorem.

---

## 1. Independence of Random Variables

### 1.1 Independence of Two Random Variables

Let $X,Y$ be random variables on the probability space $(\Omega,\mathcal{A},P)$, with joint distribution function

\[
F_{X,Y}(x,y)=P(X\leq x,Y\leq y).
\]

!!! info "Definition 1.1 (Independence of Random Variables)"

    If for any Borel sets $A,B\in\mathcal{B}(\mathbb{R})$,

    \[
    P(X\in A,Y\in B)
    =P(X\in A)P(Y\in B),
    \]

    then the random variables $X$ and $Y$ are called mutually independent.

It suffices to verify this equality on a family of sets that generates the Borel $\sigma$-algebra. Therefore, $X$ and $Y$ are independent if and only if

\[
F_{X,Y}(x,y)=F_X(x)F_Y(y),
\qquad \forall x,y\in\mathbb{R}.
\]

For discrete random variables, if the possible values of $X$ are $x_i$ and the possible values of $Y$ are $y_j$, then independence is equivalent to

\[
P(X=x_i,Y=y_j)
=P(X=x_i)P(Y=y_j),
\qquad \forall i,j.
\]

If the joint distribution has a density, then independence is equivalent to

\[
f_{X,Y}(x,y)=f_X(x)f_Y(y)
\]

for Lebesgue almost every $(x,y)$.

### 1.2 Independence Is Preserved under Measurable Transformations

!!! success "Proposition 1.2 (Functions of Independent Random Variables Are Still Independent)"

    If $X$ and $Y$ are independent, and $f,g:\mathbb{R}\to\mathbb{R}$ are Borel measurable functions, then

    \[
    U=f(X),
    \qquad
    V=g(Y)
    \]

    are also mutually independent.

??? proof "Proof of Proposition 1.2 (click to expand)"

    For any Borel sets $A,B$,

    \[
    \begin{aligned}
    P(U\in A,V\in B)
    &=P\bigl(X\in f^{-1}(A),Y\in g^{-1}(B)\bigr)\\
    &=P\bigl(X\in f^{-1}(A)\bigr)
      P\bigl(Y\in g^{-1}(B)\bigr)\\
    &=P(U\in A)P(V\in B).
    \end{aligned}
    \]

    Here $f^{-1}(A)$ and $g^{-1}(B)$ are both Borel sets. $\square$

### 1.3 Independence of Multiple Random Variables

!!! info "Definition 1.3 (Multiple Random Variables Are Mutually Independent)"

    Random variables $X_1,\ldots,X_m$ are called mutually independent if for any Borel sets $A_1,\ldots,A_m$,

    \[
    P(X_1\in A_1,\ldots,X_m\in A_m)
    =\prod_{j=1}^mP(X_j\in A_j).
    \]

If the joint density exists, then mutual independence is equivalent to

\[
f_{X_1,\ldots,X_m}(x_1,\ldots,x_m)
=\prod_{j=1}^mf_{X_j}(x_j)
\]

almost everywhere.

!!! success "Proposition 1.4 (Functions of Independent Subgroups)"

    Let $X_1,\ldots,X_m$ be mutually independent. If

    \[
    I,J\subseteq\{1,\ldots,m\},
    \qquad I\cap J=\varnothing,
    \]

    and $f$ and $g$ are Borel measurable functions of the corresponding dimensions, then

    \[
    f\bigl((X_i)_{i\in I}\bigr)
    \quad\text{and}\quad
    g\bigl((X_j)_{j\in J}\bigr)
    \]

    are mutually independent.

---

## 2. Numerical Characteristics of Random Variables

### 2.1 Measure-Theoretic Definition of Expectation

Let the distribution of the random variable $X$ be $P_X=P\circ X^{-1}$. The expectation can be written as an integral with respect to the distribution:

\[
EX=\int_{\mathbb{R}}x\,P_X(dx)
=\int_{\mathbb{R}}x\,dF_X(x).
\]

The last integral can be understood as a Lebesgue-Stieltjes integral. It unifies the discrete and continuous cases:

\[
EX=\sum_i x_iP(X=x_i)
\]

and

\[
EX=\int_{-\infty}^{\infty}xf_X(x)\,dx.
\]

Equivalently, one may also integrate directly on the original probability space:

\[
EX=\int_\Omega X(\omega)\,P(d\omega).
\]

### 2.2 From Simple Random Variables to General Random Variables

If $A\in\mathcal{A}$, then the expectation of the indicator function $\mathbf{1}_A$ is

\[
E\mathbf{1}_A=P(A).
\]

If $X$ is a simple random variable

\[
X=\sum_{i=1}^na_i\mathbf{1}_{A_i},
\]

where $A_1,\ldots,A_n$ are pairwise disjoint, then define

\[
EX=\sum_{i=1}^na_iP(A_i).
\]

For a nonnegative random variable $X\geq0$, one can take a sequence of nonnegative simple random variables $X_n$ such that

\[
0\leq X_n\uparrow X,
\]

and define

\[
EX=\lim_{n\to\infty}EX_n.
\]

For a general real-valued random variable, define the positive and negative parts

\[
X^+=\max\{X,0\},
\qquad
X^-=\max\{-X,0\},
\]

then

\[
X=X^+-X^-,
\qquad
|X|=X^++X^-.
\]

When $EX^+$ and $EX^-$ are not both infinite, one can define

\[
EX=EX^+-EX^-.
\]

In particular, $X$ is integrable if and only if

\[
E|X|=EX^++EX^-<\infty.
\]

In this case, $EX$ is a finite real number.

### 2.3 Basic Properties of Expectation

!!! success "Proposition 2.1 (Basic Properties of Expectation)"

    Assume the relevant expectations exist. Then:

    (i) If $a\leq X\leq b$, then

    \[
    a\leq EX\leq b.
    \]

    (ii) For constants $a,b$,

    \[
    E(aX+b)=aEX+b.
    \]

    (iii) If $X,Y$ are integrable, then

    \[
    E(X+Y)=EX+EY.
    \]

    (iv) If $X\leq Y$ almost everywhere, then

    \[
    EX\leq EY.
    \]

!!! success "Theorem 2.2 (Monotone Convergence Theorem)"

    If $X_n\geq0$ and

    \[
    X_n\uparrow X
    \quad\text{almost everywhere},
    \]

    then

    \[
    EX_n\uparrow EX.
    \]

### 2.4 Law of the Unconscious Statistician

!!! success "Theorem 2.3 (LOTUS)"

    Let $Y=g(X)$. If the expectation exists, then

    \[
    E\{g(X)\}
    =\int_{\mathbb{R}}g(x)\,P_X(dx)
    =\int_{\mathbb{R}}g(x)\,dF_X(x).
    \]

    In the discrete case,

    \[
    E\{g(X)\}
    =\sum_i g(x_i)P(X=x_i),
    \]

    and in the continuous case with a density,

    \[
    E\{g(X)\}
    =\int_{-\infty}^{\infty}g(x)f_X(x)\,dx.
    \]

This formula does not require first finding the distribution of $Y=g(X)$.

### 2.5 Variance, Covariance, and Moments

!!! info "Definition 2.4 (Variance)"

    If $EX^2<\infty$, then the variance of $X$ is defined as

    \[
    \operatorname{Var}(X)
    =E(X-EX)^2
    =EX^2-(EX)^2.
    \]

!!! info "Definition 2.5 (Moments and Central Moments)"

    If the expectation exists, then $EX^k$ is called the $k$-th raw moment of $X$, and $E(X-EX)^k$ is called the $k$-th central moment of $X$.

!!! info "Definition 2.6 (Covariance)"

    If $X,Y$ are square integrable, then

    \[
    \operatorname{Cov}(X,Y)
    =E\{(X-EX)(Y-EY)\}.
    \]

By direct expansion, we obtain

\[
\operatorname{Var}(X+Y)
=\operatorname{Var}(X)
+\operatorname{Var}(Y)
+2\operatorname{Cov}(X,Y).
\]

If $X_1,\ldots,X_n$ are pairwise uncorrelated, then

\[
\operatorname{Var}\left(\sum_{i=1}^nX_i\right)
=\sum_{i=1}^n\operatorname{Var}(X_i).
\]

In particular, mutually independent and square integrable random variables are pairwise uncorrelated, so the above identity holds.

### 2.6 Chebyshev's Inequality

!!! success "Theorem 2.7 (Chebyshev's Inequality)"

    If $EX^2<\infty$, then for any $\varepsilon>0$,

    \[
    P(|X-EX|\geq\varepsilon)
    \leq\frac{\operatorname{Var}(X)}{\varepsilon^2}.
    \]

??? proof "Proof of Theorem 2.7 (click to expand)"

    On the event $\{|X-EX|\geq\varepsilon\}$,

    \[
    \mathbf{1}_{\{|X-EX|\geq\varepsilon\}}
    \leq\frac{|X-EX|^2}{\varepsilon^2}.
    \]

    Taking expectations of both sides gives

    \[
    \begin{aligned}
    P(|X-EX|\geq\varepsilon)
    &=E\mathbf{1}_{\{|X-EX|\geq\varepsilon\}}\\
    &\leq\frac{E|X-EX|^2}{\varepsilon^2}\\
    &=\frac{\operatorname{Var}(X)}{\varepsilon^2}.
    \end{aligned}
    \]

    This completes the proof. $\square$

### 2.7 Moments of Common Distributions

!!! example "Example 2.8 (Cauchy Distribution)"

    If $X$ follows the standard Cauchy distribution, then

    \[
    f_X(x)=\frac{1}{\pi(1+x^2)}.
    \]

    Since $E|X|=\infty$, $EX$ does not exist.

!!! example "Example 2.9 (Normal Distribution)"

    If $X\sim N(0,\sigma^2)$, then for $k\geq0$,

    \[
    EX^{2k+1}=0,
    \]

    \[
    EX^{2k}=(2k-1)!!\,\sigma^{2k}.
    \]

!!! example "Example 2.10 (Poisson Distribution)"

    If $X\sim\operatorname{Poisson}(\lambda)$, then the falling factorial moments satisfy

    \[
    E\{X(X-1)\cdots(X-m+1)\}
    =\lambda^m.
    \]

---

## 3. Moment Generating Functions and Characteristic Functions

### 3.1 Moment Generating Function

!!! info "Definition 3.1 (Moment Generating Function)"

    The moment generating function of a random variable $X$ is defined as

    \[
    M_X(t)=Ee^{tX},
    \]

    where $t$ is taken over real numbers for which this expectation is finite.

If $M_X(t)$ is finite in some neighborhood of $t=0$, then it uniquely determines the distribution of $X$ in that neighborhood, and under the corresponding differentiability conditions

\[
M_X^{(k)}(0)=EX^k.
\]

!!! warning "Moments Do Not Always Uniquely Determine a Distribution"

    Merely knowing that all moments are equal generally does not guarantee that two random variables have the same distribution. If the moment generating function exists in a neighborhood of the origin, then the uniqueness of the moment generating function can be used to obtain uniqueness of the distribution.

### 3.2 Definition of the Characteristic Function

!!! info "Definition 3.2 (Characteristic Function)"

    The characteristic function of a random variable $X$ is defined as

    \[
    \varphi_X(t)=Ee^{itX},
    \qquad t\in\mathbb{R}.
    \]

    Equivalently,

    \[
    \varphi_X(t)
    =E\cos(tX)+iE\sin(tX)
    =\int_{\mathbb{R}}e^{itx}\,P_X(dx).
    \]

Unlike the moment generating function, the characteristic function exists for every real-valued random variable and every $t\in\mathbb{R}$, because

\[
|e^{itX}|=1.
\]

### 3.3 Characteristic Functions of Common Distributions

!!! example "Example 3.3 (Poisson Distribution)"

    If $X\sim\operatorname{Poisson}(\lambda)$, then

    \[
    \varphi_X(t)
    =\exp\{\lambda(e^{it}-1)\}.
    \]

!!! example "Example 3.4 (Normal Distribution)"

    If $X\sim N(\mu,\sigma^2)$, then

    \[
    \varphi_X(t)
    =\exp\left(i\mu t-\frac{\sigma^2t^2}{2}\right).
    \]

!!! example "Example 3.5 (Standard Cauchy Distribution)"

    If $X$ follows the standard Cauchy distribution, then

    \[
    \varphi_X(t)=e^{-|t|}.
    \]

### 3.4 Basic Properties of Characteristic Functions

!!! success "Proposition 3.6 (Basic Properties of Characteristic Functions)"

    Any characteristic function satisfies:

    (i)

    \[
    \varphi_X(0)=1.
    \]

    (ii)

    \[
    |\varphi_X(t)|\leq1,
    \qquad \forall t\in\mathbb{R}.
    \]

    (iii)

    \[
    \varphi_X(-t)=\overline{\varphi_X(t)}.
    \]

    (iv) $\varphi_X$ is uniformly continuous on $\mathbb{R}$.

??? proof "Proof of Uniform Continuity (click to expand)"

    For any $t,h\in\mathbb{R}$,

    \[
    \begin{aligned}
    |\varphi_X(t+h)-\varphi_X(t)|
    &=\left|E\left[e^{itX}(e^{ihX}-1)\right]\right|\\
    &\leq E|e^{ihX}-1|.
    \end{aligned}
    \]

    As $h\to0$, $|e^{ihX}-1|\to0$ almost everywhere, and

    \[
    |e^{ihX}-1|\leq2.
    \]

    By the dominated convergence theorem,

    \[
    E|e^{ihX}-1|\to0.
    \]

    The right-hand side is independent of $t$, so $\varphi_X$ is uniformly continuous. $\square$

### 3.5 Nonnegative Definiteness

!!! success "Proposition 3.7 (Bochner Nonnegative Definiteness)"

    For any $n\geq1$, complex numbers $a_1,\ldots,a_n$, and real numbers $t_1,\ldots,t_n$,

    \[
    \sum_{j=1}^n\sum_{k=1}^n
    a_j\overline{a_k}\,
    \varphi_X(t_j-t_k)
    \geq0.
    \]

??? proof "Proof of Proposition 3.7 (click to expand)"

    By the definition of the characteristic function,

    \[
    \begin{aligned}
    \sum_{j,k}a_j\overline{a_k}\varphi_X(t_j-t_k)
    &=E\left(\sum_ja_je^{it_jX}\right)
      \overline{\left(\sum_ka_ke^{it_kX}\right)}\\
    &=E\left|\sum_{j=1}^na_je^{it_jX}\right|^2
    \geq0.
    \end{aligned}
    \]

    This completes the proof. $\square$

### 3.6 Derivatives of Characteristic Functions and Moments

!!! success "Theorem 3.8 (Derivatives of Characteristic Functions)"

    If $E|X|^k<\infty$, then $\varphi_X$ is at least $k$ times differentiable, and

    \[
    \varphi_X^{(j)}(t)
    =E\left[(iX)^je^{itX}\right],
    \qquad j=1,\ldots,k.
    \]

    In particular,

    \[
    \varphi_X^{(j)}(0)=i^jEX^j.
    \]

Therefore, when the corresponding moments exist, the characteristic function has the expansion near the origin

\[
\varphi_X(t)
=\sum_{j=0}^k\frac{i^jEX^j}{j!}t^j+o(t^k).
\]

### 3.7 Sums of Independent Random Variables

!!! success "Proposition 3.9 (Characteristic Function of a Sum of Independent Random Variables)"

    If $X$ and $Y$ are independent, then

    \[
    \varphi_{X+Y}(t)
    =\varphi_X(t)\varphi_Y(t).
    \]

    More generally, if $X_1,\ldots,X_n$ are mutually independent, then

    \[
    \varphi_{\sum_{j=1}^nX_j}(t)
    =\prod_{j=1}^n\varphi_{X_j}(t).
    \]

??? proof "Proof of Proposition 3.9 (click to expand)"

    By independence,

    \[
    \begin{aligned}
    \varphi_{X+Y}(t)
    &=E e^{it(X+Y)}\\
    &=E(e^{itX}e^{itY})\\
    &=Ee^{itX}\,Ee^{itY}.
    \end{aligned}
    \]

    The case of multiple random variables follows by induction. $\square$

!!! example "Example 3.10 (Sum of Independent Normal Variables)"

    If

    \[
    X_i\sim N(\mu_i,\sigma_i^2),
    \qquad i=1,\ldots,n,
    \]

    and $X_1,\ldots,X_n$ are mutually independent, then

    \[
    \sum_{i=1}^nX_i
    \sim N\left(\sum_{i=1}^n\mu_i,
    \sum_{i=1}^n\sigma_i^2\right).
    \]

    This is because the product of the corresponding characteristic functions is

    \[
    \exp\left\{
    it\sum_{i=1}^n\mu_i
    -\frac{t^2}{2}\sum_{i=1}^n\sigma_i^2
    \right\},
    \]

    which is exactly the characteristic function of the above normal distribution.

---

## 4. Characteristic Functions and Distributions

### 4.1 Uniqueness Theorem

!!! success "Theorem 4.1 (Uniqueness of Characteristic Functions)"

    The characteristic function uniquely determines the probability distribution. That is,

    \[
    X\overset d=Y
    \quad\Longleftrightarrow\quad
    \varphi_X(t)=\varphi_Y(t),
    \qquad \forall t\in\mathbb{R}.
    \]

### 4.2 Inversion Formula

!!! success "Theorem 4.2 (Lévy Inversion Formula)"

    Let $F$ be the distribution function of $X$, let $x_1<x_2$, and suppose $x_1,x_2$ are both continuity points of $F$. Then

    \[
    F(x_2)-F(x_1)
    =\lim_{T\to\infty}\frac{1}{2\pi}
    \int_{-T}^T
    \frac{e^{-itx_1}-e^{-itx_2}}{it}
    \varphi_X(t)\,dt.
    \]

This formula shows that one can recover the probability mass of the distribution function between continuity points from the characteristic function, thereby obtaining the uniqueness theorem.

!!! note "Discontinuity Points of Distribution Functions"

    The jump size of the distribution function at point $x$ is

    \[
    F(x)-F(x-)=P(X=x).
    \]

    Therefore, the discontinuity points of a distribution function are at most countable.

### 4.3 Fourier Inversion of Densities

Under suitable integrability conditions, if $X$ has density $f_X$, then

\[
\varphi_X(t)
=\int_{-\infty}^{\infty}e^{itx}f_X(x)\,dx,
\]

and by Fourier inversion one obtains

\[
f_X(x)
=\frac{1}{2\pi}
\int_{-\infty}^{\infty}e^{-itx}\varphi_X(t)\,dt.
\]

!!! example "Example 4.3 (Symmetric Two-Point Distribution)"

    If

    \[
    P(X=1)=P(X=-1)=\frac12,
    \]

    then

    \[
    \varphi_X(t)
    =\frac12e^{it}+\frac12e^{-it}
    =\cos t.
    \]

---

## 5. Convergence of Sequences of Random Variables

### 5.1 Convergence in Probability

!!! info "Definition 5.1 (Convergence in Probability)"

    If for every $\varepsilon>0$,

    \[
    P(|X_n-X|>\varepsilon)\to0,
    \]

    then $X_n$ is said to converge to $X$ in probability, denoted by

    \[
    X_n\xrightarrow{P}X.
    \]

### 5.2 Convergence in Distribution

!!! info "Definition 5.2 (Convergence in Distribution)"

    Let $F_n$ and $F$ be the distribution functions of $X_n$ and $X$, respectively. If for every continuity point $x$ of $F$,

    \[
    F_n(x)\to F(x),
    \]

    then $X_n$ is said to converge to $X$ in distribution, denoted by

    \[
    X_n\xrightarrow{d}X.
    \]

!!! warning "Why Only Convergence at Continuity Points of the Limiting Distribution Is Required"

    Let $x_n=1/n$, and regard $X_n\equiv x_n$ and $X\equiv0$ as degenerate random variables. Clearly

    \[
    X_n\to X
    \]

    holds pointwise, and hence $X_n\xrightarrow{d}X$. However, at $x=0$,

    \[
    F_n(0)=0,
    \qquad
    F(0)=1.
    \]

    Therefore, the distribution functions need not converge at discontinuity points of the limiting distribution.

Convergence in probability implies convergence in distribution, but the converse generally does not hold.

### 5.3 Lévy Continuity Theorem

!!! success "Theorem 5.3 (Lévy Continuity Theorem)"

    If the characteristic functions of $X_n$ and $X$ are $\varphi_n$ and $\varphi$, respectively, then

    \[
    X_n\xrightarrow{d}X
    \quad\Longleftrightarrow\quad
    \varphi_n(t)\to\varphi(t),
    \qquad \forall t\in\mathbb{R}.
    \]

More generally, if $\varphi_n(t)$ converges pointwise to a function $\varphi(t)$, and $\varphi$ is continuous at $0$, then $\varphi$ is the characteristic function of some probability distribution, and the corresponding distributions converge weakly to that distribution.

---

## 6. Limit Theorems for Sums of Bernoulli Random Variables

Let event $A$ have probability $p$, and let

\[
\xi_i(\omega)=
\begin{cases}
1,&\omega_i\in A,\\
0,&\omega_i\notin A,
\end{cases}
\]

where $\xi_1,\ldots,\xi_n$ are mutually independent, and set

\[
S_n=\sum_{i=1}^n\xi_i.
\]

Then

\[
S_n\sim\operatorname{Binomial}(n,p),
\]

that is,

\[
P(S_n=k)=\binom nkp^k(1-p)^{n-k}.
\]

### 6.1 Bernoulli Weak Law of Large Numbers

!!! success "Theorem 6.1 (Bernoulli Weak Law of Large Numbers)"

    For fixed $0<p<1$,

    \[
    \frac{S_n}{n}\xrightarrow{P}p.
    \]

??? proof "Proof of Theorem 6.1 (click to expand)"

    Because

    \[
    E\left(\frac{S_n}{n}\right)=p,
    \qquad
    \operatorname{Var}\left(\frac{S_n}{n}\right)
    =\frac{p(1-p)}{n},
    \]

    by Chebyshev's inequality, for every $\varepsilon>0$,

    \[
    P\left(\left|\frac{S_n}{n}-p\right|>\varepsilon\right)
    \leq\frac{p(1-p)}{n\varepsilon^2}
    \to0.
    \]

    This completes the proof. $\square$

### 6.2 Poisson Limit Theorem

!!! success "Theorem 6.2 (Poisson Limit Theorem)"

    Let

    \[
    S_n\sim\operatorname{Binomial}(n,p_n),
    \qquad np_n\to\lambda\in(0,\infty).
    \]

    Then for every fixed nonnegative integer $k$,

    \[
    P(S_n=k)
    \to e^{-\lambda}\frac{\lambda^k}{k!}.
    \]

    Therefore

    \[
    S_n\xrightarrow{d}\operatorname{Poisson}(\lambda).
    \]

### 6.3 De Moivre-Laplace Theorem

Let $q=1-p$.

!!! success "Theorem 6.3 (De Moivre-Laplace Central Limit Theorem)"

    If $S_n\sim\operatorname{Binomial}(n,p)$, where $0<p<1$ is fixed, then

    \[
    \frac{S_n-np}{\sqrt{npq}}
    \xrightarrow{d}N(0,1).
    \]

    Equivalently, for every $x\in\mathbb{R}$,

    \[
    P\left(
    \frac{S_n-np}{\sqrt{npq}}\leq x
    \right)
    \to
    \Phi(x)
    =\frac{1}{\sqrt{2\pi}}
    \int_{-\infty}^xe^{-t^2/2}\,dt.
    \]

This theorem can be proved starting from the binomial probability mass function by using Stirling's formula

\[
n!\sim\sqrt{2\pi n}\left(\frac ne\right)^n.
\]

---

## 7. Characteristic Function Proof of the Central Limit Theorem

!!! success "Theorem 7.1 (Lindeberg-Lévy Central Limit Theorem)"

    Let $X_1,X_2,\ldots$ be independent and identically distributed, and

    \[
    EX_1=\mu,
    \qquad
    \operatorname{Var}(X_1)=\sigma^2\in(0,\infty).
    \]

    Let

    \[
    S_n=\sum_{i=1}^nX_i.
    \]

    Then

    \[
    \frac{S_n-n\mu}{\sigma\sqrt n}
    \xrightarrow{d}N(0,1).
    \]

??? proof "Characteristic Function Proof (click to expand)"

    Let

    \[
    Y_i=\frac{X_i-\mu}{\sigma},
    \]

    then $EY_i=0$ and $EY_i^2=1$. Denote the characteristic function of $Y_i$ by $\varphi_Y$. From the second-order expansion near the origin,

    \[
    \varphi_Y(t)
    =1-\frac{t^2}{2}+o(t^2),
    \qquad t\to0.
    \]

    By independence, the characteristic function of the standardized sum is

    \[
    \varphi_{n}(t)
    =\left[\varphi_Y\left(\frac{t}{\sqrt n}\right)\right]^n.
    \]

    Therefore

    \[
    \begin{aligned}
    \varphi_n(t)
    &=\left[
    1-\frac{t^2}{2n}+o\left(\frac1n\right)
    \right]^n\\
    &\to e^{-t^2/2}.
    \end{aligned}
    \]

    The function $e^{-t^2/2}$ is the characteristic function of the standard normal distribution. By the Lévy continuity theorem,

    \[
    \frac{S_n-n\mu}{\sigma\sqrt n}
    \xrightarrow{d}N(0,1).
    \]

    This completes the proof. $\square$

---

## 8. Summary of This Chapter

The main conclusions of this chapter are as follows:

1. The independence of random variables can be characterized by the product form of the joint distribution function, probability mass function, or joint density, and is preserved after applying measurable transformations separately.
2. Expectation can be uniformly represented as a Lebesgue-Stieltjes integral with respect to the distribution; integrability is equivalent to $E|X|<\infty$.
3. Variance, covariance, and moments characterize the numerical characteristics of random variables, and Chebyshev's inequality converts variance into an upper bound for tail probabilities.
4. Characteristic functions always exist and uniquely determine the distribution; the characteristic function of a sum of independent random variables equals the product of their respective characteristic functions.
5. The Lévy continuity theorem connects pointwise convergence of characteristic functions with convergence in distribution.
6. The Bernoulli weak law of large numbers, the Poisson limit theorem, the De Moivre-Laplace theorem, and the Lindeberg-Lévy central limit theorem describe, respectively, the limiting behavior of frequency stabilization, rare events, and standardized sums.