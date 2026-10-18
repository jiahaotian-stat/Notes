# Chapter 8: Consistency, Moment Estimators, and the Delta Method

This chapter introduces the basic workflow of large-sample statistical inference. We first review modes of convergence, the continuous mapping theorem, and Slutsky's theorem, and then discuss identifiability, separation, and consistency of minimum-distance estimators. We next study empirical raw and central moments, consistency and asymptotic normality of moment estimators, and finally the multivariate delta method together with two supplementary proofs.

---

## 1. The Large-Sample Workflow

Let

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}
P_\theta,
\]

and let the statistic $T_n$ estimate a target $g(\theta)$. Large-sample analysis asks three successive questions:

| Question | Statistical Concept |
| --- | --- |
| Does $T_n$ approach $g(\theta)$? | Consistency |
| At what scale does $T_n-g(\theta)$ fluctuate? | Rate of convergence |
| What is the limiting distribution? | Inference |

Consistency is therefore not the final goal, but it is usually the first condition needed to justify a local expansion.

---

## 2. Modes of Convergence and Two Transfer Tools

Let $Z_n,Z$ be random vectors in a metric space with metric $d$.

!!! info "Definition 2.1 (Four Modes of Convergence)"

    **Convergence in probability:**

    \[
    Z_n\overset{p}{\longrightarrow}Z
    \]

    means that, for every $\epsilon>0$,

    \[
    P\{d(Z_n,Z)>\epsilon\}\longrightarrow0.
    \]

    **Almost sure convergence:**

    \[
    Z_n\overset{\mathrm{a.s.}}{\longrightarrow}Z
    \]

    means

    \[
    P\left(\lim_{n\to\infty}Z_n=Z\right)=1.
    \]

    **Convergence in $L^r$:**

    \[
    Z_n\overset{L^r}{\longrightarrow}Z
    \]

    means

    \[
    E\,d(Z_n,Z)^r\longrightarrow0.
    \]

    **Convergence in distribution:**

    \[
    Z_n\overset{d}{\longrightarrow}Z
    \]

    means that, for every bounded continuous function $f$,

    \[
    Ef(Z_n)\longrightarrow Ef(Z).
    \]

The standard implications are

\[
Z_n\overset{\mathrm{a.s.}}{\longrightarrow}Z
\quad\Longrightarrow\quad
Z_n\overset{p}{\longrightarrow}Z
\quad\Longrightarrow\quad
Z_n\overset{d}{\longrightarrow}Z,
\]

and

\[
Z_n\overset{L^r}{\longrightarrow}Z
\quad\Longrightarrow\quad
Z_n\overset{p}{\longrightarrow}Z.
\]

!!! success "Theorem 2.2 (Continuous Mapping and Slutsky)"

    If

    \[
    Z_n\overset{p}{\longrightarrow}Z
    \]

    and $g$ is continuous at $Z$ almost surely, then

    \[
    g(Z_n)\overset{p}{\longrightarrow}g(Z).
    \]

    Corresponding continuous mapping results hold for almost sure convergence and convergence in distribution.

    In particular, if

    \[
    Z_n\overset{d}{\longrightarrow}Z,
    \qquad
    W_n\overset{p}{\longrightarrow}c,
    \]

    then

    \[
    Z_n+W_n\overset{d}{\longrightarrow}Z+c,
    \]

    \[
    Z_nW_n\overset{d}{\longrightarrow}cZ,
    \]

    and, when $c\neq0$,

    \[
    \frac{Z_n}{W_n}
    \overset{d}{\longrightarrow}
    \frac{Z}{c}.
    \]

??? proof "Proof idea (click to expand)"

    Continuity converts the event

    \[
    d\{g(Z_n),g(Z)\}>\epsilon
    \]

    into the event that $Z_n$ lies outside a sufficiently small neighborhood of $Z$.

    For Slutsky's theorem, first obtain the joint convergence

    \[
    (Z_n,W_n)
    \overset{d}{\longrightarrow}
    (Z,c),
    \]

    because the second component converges to a constant. Then apply the continuous mapping theorem to addition, multiplication, and division. $\square$

---

## 3. Consistency and Identifiability

!!! info "Definition 3.1 (Weak and Strong Consistency)"

    If, for every $\theta\in\Theta$,

    \[
    T_n\overset{p}{\longrightarrow}g(\theta),
    \]

    then $T_n$ is weakly consistent for $g(\theta)$.

    If, for every $\theta\in\Theta$,

    \[
    T_n\overset{\mathrm{a.s.}}{\longrightarrow}g(\theta),
    \]

    then $T_n$ is strongly consistent for $g(\theta)$.

!!! info "Definition 3.2 (Identifiability)"

    The model is identifiable if

    \[
    P_{\theta_1}=P_{\theta_2}
    \quad\Longrightarrow\quad
    \theta_1=\theta_2.
    \]

    If only a functional $g(\theta)$ is to be estimated, the necessary condition is

    \[
    P_{\theta_1}=P_{\theta_2}
    \quad\Longrightarrow\quad
    g(\theta_1)=g(\theta_2).
    \]

!!! success "Proposition 3.3 (Identifiability Is Necessary)"

    If

    \[
    P_{\theta_1}=P_{\theta_2},
    \qquad
    g(\theta_1)\neq g(\theta_2),
    \]

    then no estimator can consistently estimate $g(\theta)$ at both parameter values.

??? proof "Proof of Proposition 3.3 (click to expand)"

    Every statistic $T_n$ has the same distribution under $\theta_1$ and $\theta_2$. Thus, if it converges in probability to $g(\theta_1)$ under $P_{\theta_1}$, the same convergence holds under $P_{\theta_2}$.

    On the other hand, consistency at $\theta_2$ would require convergence in probability to $g(\theta_2)$. A sequence of random variables cannot converge in probability to two distinct constants. $\square$

### 3.1 Why Identifiability Is Not Sufficient

Identifiability is only a pointwise uniqueness statement. Consistency additionally requires the data distribution to separate parameter values in a stable way.

Consider the following intentionally pathological example:

\[
X\sim\operatorname{Bernoulli}\{h(\theta)\},
\qquad
h(\theta)=
\begin{cases}
\theta, & \theta\in\mathbb{Q},\\
1-\theta, & \theta\notin\mathbb{Q},
\end{cases}
\qquad
0<\theta<1.
\tag{1}
\]

The map $h$ is one-to-one, so the model is identifiable, yet no estimator is consistent at every $\theta$. The complete proof is given in Section 9. It uses the fact that a pointwise limit of continuous functions must have many continuity points, whereas $h$ is discontinuous everywhere except at $1/2$.

!!! note "Practical Interpretation"

    This example is not a model that students will encounter in applications. Its role is to separate two logical statements:

    \[
    \text{identifiability}
    =\text{a unique population label},
    \]

    \[
    \text{consistency}
    =\text{stable recovery from empirical data}.
    \]

    In regular parametric models, compactness, continuity, and separation usually provide the missing stability.

---

## 4. Separation and Minimum-Distance Consistency

Let $F_\theta$ be the distribution function under $P_\theta$, let $F_n$ be the empirical distribution function, and define

\[
d(F,G)
=\sup_x|F(x)-G(x)|.
\]

At the true parameter $\theta_0$, suppose that, for every $\epsilon>0$,

\[
\eta(\theta_0,\epsilon)
:=
\inf_{\|\theta-\theta_0\|\geq\epsilon}
d(F_\theta,F_{\theta_0})
>0.
\tag{2}
\]

Condition (2) is often called strong separation or well-separated identification. It is stronger than pointwise identifiability.

Let $\widehat{\theta}_n$ be an approximate minimum-distance estimator satisfying

\[
d(F_{\widehat{\theta}_n},F_n)
\leq
\inf_{\theta\in\Theta}d(F_\theta,F_n)
+r_n,
\qquad
r_n\overset{\mathrm{a.s.}}{\longrightarrow}0.
\tag{3}
\]

!!! success "Theorem 4.1 (Strong Consistency of Minimum Distance)"

    Under conditions (2)-(3),

    \[
    \widehat{\theta}_n
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \theta_0.
    \]

    If empirical convergence and optimization error hold only in probability, the same argument gives weak consistency.

??? proof "Proof of Theorem 4.1 (click to expand)"

    The Glivenko-Cantelli theorem gives

    \[
    \Delta_n
    :=d(F_n,F_{\theta_0})
    \overset{\mathrm{a.s.}}{\longrightarrow}0.
    \]

    Since $\theta_0$ is an admissible competitor in (3),

    \[
    d(F_{\widehat{\theta}_n},F_n)
    \leq\Delta_n+r_n.
    \]

    The triangle inequality then yields

    \[
    \begin{aligned}
    d(F_{\widehat{\theta}_n},F_{\theta_0})
    &\leq
    d(F_{\widehat{\theta}_n},F_n)
    +d(F_n,F_{\theta_0})\\
    &\leq
    2\Delta_n+r_n
    \overset{\mathrm{a.s.}}{\longrightarrow}0.
    \end{aligned}
    \]

    Fix $\epsilon>0$. If

    \[
    \|\widehat{\theta}_n-\theta_0\|
    \geq\epsilon,
    \]

    then condition (2) forces

    \[
    d(F_{\widehat{\theta}_n},F_{\theta_0})
    \geq\eta(\theta_0,\epsilon)>0.
    \]

    This is eventually impossible on the probability-one event where $2\Delta_n+r_n\to0$. Therefore,

    \[
    \widehat{\theta}_n
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \theta_0.
    \]

    $\square$

### 4.1 Where Strong Convexity Fits

For a general criterion estimator

\[
\widehat{\theta}_n
\in\operatorname*{arg\,min}_\theta Q_n(\theta),
\]

the population analogue of condition (2) is

\[
\inf_{\|\theta-\theta_0\|\geq\epsilon}
\{Q(\theta)-Q(\theta_0)\}
>0.
\tag{4}
\]

If $Q$ is $m$-strongly convex and minimized at $\theta_0$, then

\[
Q(\theta)-Q(\theta_0)
\geq
\frac{m}{2}\|\theta-\theta_0\|^2,
\tag{5}
\]

so (4) follows. Strong convexity is sufficient but not necessary; distribution distances may be nonsmooth and still satisfy separation.

Strong convexity also gives a quantitative result. Let

\[
\Delta_n
=\sup_\theta|Q_n(\theta)-Q(\theta)|,
\]

and suppose

\[
Q_n(\widehat{\theta}_n)
\leq Q_n(\theta_0)+r_n.
\]

Then

\[
\|\widehat{\theta}_n-\theta_0\|
\leq
\left(
\frac{4\Delta_n+2r_n}{m}
\right)^{1/2}.
\tag{6}
\]

Indeed,

\[
Q(\widehat{\theta}_n)-Q(\theta_0)
\leq2\Delta_n+r_n,
\]

and combining this with (5) gives (6). Thus almost sure uniform convergence yields strong consistency, whereas uniform convergence in probability yields weak consistency.

---

## 5. Raw and Central Moments

Let

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}P,
\]

and suppose the required moments exist. Define the population raw and central moments

\[
\alpha_r=EX^r,
\qquad
\mu_r=E(X-\mu)^r,
\qquad
\mu=\alpha_1,
\]

and their empirical versions

\[
A_{r,n}
=\frac{1}{n}\sum_{i=1}^nX_i^r,
\qquad
M_{r,n}
=\frac{1}{n}\sum_{i=1}^n(X_i-\overline{X}_n)^r.
\]

A central moment is a moment about the mean, whereas a raw moment is a moment about zero.

!!! success "Theorem 5.1 (Strong Consistency of Empirical Moments)"

    If

    \[
    E|X|^r<\infty,
    \]

    then

    \[
    A_{r,n}
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \alpha_r,
    \qquad
    M_{r,n}
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \mu_r.
    \tag{7}
    \]

    The conclusion also holds jointly for any finite collection of raw and central moments whose required absolute moments exist.

??? proof "Proof of Theorem 5.1 (click to expand)"

    For every $j\leq r$, $E|X|^j<\infty$, so the strong law gives

    \[
    A_{j,n}
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \alpha_j.
    \]

    Let $A_{0,n}=\alpha_0=1$. The binomial identity gives

    \[
    M_{r,n}
    =\sum_{j=0}^r
    \binom{r}{j}
    (-\overline{X}_n)^{r-j}A_{j,n},
    \]

    and

    \[
    \mu_r
    =\sum_{j=0}^r
    \binom{r}{j}
    (-\mu)^{r-j}\alpha_j.
    \tag{8}
    \]

    Since

    \[
    \overline{X}_n
    =A_{1,n}
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \mu,
    \]

    the continuous mapping theorem applied to the polynomial in (8) gives

    \[
    M_{r,n}
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \mu_r.
    \]

    Applying the same argument to a finite-dimensional vector proves joint convergence. $\square$

Bias is different from consistency. $A_{r,n}$ is exactly unbiased for $\alpha_r$, whereas $M_{r,n}$ is generally biased for $\mu_r$. For example,

\[
EM_{2,n}
=\frac{n-1}{n}\mu_2.
\]

Under

\[
E|X|^{2r}<\infty,
\]

the calculation in Section 10 gives

\[
EM_{r,n}-\mu_r
=O(n^{-1}).
\tag{9}
\]

Equation (9) is useful for finite-sample bias, but it is not needed for the simpler strong-consistency proof above.

---

## 6. Both Versions of the Method of Moments

Suppose the parameter is characterized by a vector containing raw moments, central moments, or both:

\[
m(\theta)
=\bigl(
\alpha_{r_1}(\theta),\ldots,\alpha_{r_k}(\theta),
\mu_{t_1}(\theta),\ldots,\mu_{t_\ell}(\theta)
\bigr).
\]

The empirical vector $\widehat{m}_n$ replaces each population moment by the corresponding $A_{r,n}$ or $M_{t,n}$.

A method-of-moments estimator may be obtained by solving

\[
m(\widehat{\theta}_n)
=\widehat{m}_n,
\]

or, when a local inverse $h$ is available, by defining

\[
\widehat{\theta}_n
=h(\widehat{m}_n),
\qquad
h\{m(\theta)\}=\theta.
\]

!!! success "Corollary 6.1 (Consistency of Moment Estimators)"

    If the required absolute moments exist and $h$ is continuous at $m(\theta)$, then

    \[
    \widehat{\theta}_n
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \theta.
    \]

    If only weak convergence of $\widehat{m}_n$ is available, the conclusion is weak consistency.

??? proof "Proof of Corollary 6.1 (click to expand)"

    Theorem 5.1 gives

    \[
    \widehat{m}_n
    \overset{\mathrm{a.s.}}{\longrightarrow}
    m(\theta).
    \]

    Apply the continuous mapping theorem to the continuous map $h$. $\square$

!!! warning "The Choice of Moments Matters"

    The classical method of moments is often introduced using the first $k$ raw moments. Central moments are equally legitimate when they identify the desired scale or shape.

    However, central moments alone cannot identify a location parameter: $\mu_1=0$, and every central moment is invariant under location shifts. In a location-scale problem, a natural choice is the first raw moment for location and the second central moment for scale.

---

## 7. Asymptotic Normality of Empirical Moments

For fixed positive integers $r_1,\ldots,r_k$, suppose

\[
E|X|^{2\max_jr_j}<\infty.
\]

The multivariate central limit theorem gives

\[
\sqrt{n}
\begin{bmatrix}
A_{r_1,n}-\alpha_{r_1}\\
\vdots\\
A_{r_k,n}-\alpha_{r_k}
\end{bmatrix}
\overset{d}{\longrightarrow}
N_k(0,B),
\]

where

\[
B_{ij}
=\alpha_{r_i+r_j}
-\alpha_{r_i}\alpha_{r_j}.
\tag{10}
\]

For a central moment, the polynomial representation (8) and the delta method yield a particularly useful asymptotic linear expansion.

!!! success "Theorem 7.1 (Central-Moment Expansion)"

    If

    \[
    E|X|^{2r}<\infty,
    \]

    then

    \[
    \sqrt{n}(M_{r,n}-\mu_r)
    =\frac{1}{\sqrt{n}}
    \sum_{i=1}^n\psi_r(X_i)
    +o_p(1),
    \tag{11}
    \]

    where

    \[
    \psi_r(x)
    =(x-\mu)^r-\mu_r
    -r\mu_{r-1}(x-\mu).
    \]

    Consequently,

    \[
    \sqrt{n}(M_{r,n}-\mu_r)
    \overset{d}{\longrightarrow}
    N\{0,\operatorname{Var}(\psi_r(X))\}.
    \]

??? proof "Proof of Theorem 7.1 (click to expand)"

    Let

    \[
    Y_i=X_i-\mu,
    \qquad
    \overline{Y}_n=\overline{X}_n-\mu.
    \]

    A first-order expansion around $\overline{Y}_n=0$ gives

    \[
    \begin{aligned}
    M_{r,n}
    &=\frac{1}{n}\sum_{i=1}^n
    (Y_i-\overline{Y}_n)^r\\
    &=\frac{1}{n}\sum_{i=1}^nY_i^r
    -r\overline{Y}_n
    \frac{1}{n}\sum_{i=1}^nY_i^{r-1}
    +O_p(\overline{Y}_n^2).
    \end{aligned}
    \]

    Here,

    \[
    \overline{Y}_n=O_p(n^{-1/2}),
    \qquad
    \frac{1}{n}\sum_{i=1}^nY_i^{r-1}
    \overset{p}{\longrightarrow}\mu_{r-1}.
    \]

    Thus, after subtracting $\mu_r$ and multiplying by $\sqrt{n}$,

    \[
    \begin{aligned}
    \sqrt{n}(M_{r,n}-\mu_r)
    &=
    \frac{1}{\sqrt{n}}
    \sum_{i=1}^n(Y_i^r-\mu_r)\\
    &\quad
    -r\mu_{r-1}
    \frac{1}{\sqrt{n}}
    \sum_{i=1}^nY_i
    +o_p(1).
    \end{aligned}
    \]

    This is (11), and the scalar central limit theorem completes the proof. $\square$

For $r=2$, $\mu_1=0$, so

\[
\psi_2(x)
=(x-\mu)^2-\sigma^2.
\]

Therefore,

\[
\sqrt{n}(M_{2,n}-\sigma^2)
\overset{d}{\longrightarrow}
N(0,\mu_4-\sigma^4).
\]

---

## 8. The Delta Method

!!! success "Theorem 8.1 (Multivariate Delta Method)"

    Let $\xi_n\in\mathbb{R}^k$ satisfy

    \[
    \sqrt{n}(\xi_n-a)
    \overset{d}{\longrightarrow}
    N_k(0,B).
    \]

    If $g:\mathbb{R}^k\to\mathbb{R}^m$ is differentiable at $a$, with Jacobian

    \[
    C=Dg(a),
    \]

    then

    \[
    \sqrt{n}\{g(\xi_n)-g(a)\}
    \overset{d}{\longrightarrow}
    N_m(0,CBC^{\mathsf T}).
    \tag{12}
    \]

    More generally, if $s_n\to\infty$ and

    \[
    s_n(\xi_n-a)
    \overset{d}{\longrightarrow}Z,
    \]

    then

    \[
    s_n\{g(\xi_n)-g(a)\}
    \overset{d}{\longrightarrow}
    Dg(a)Z.
    \]

??? proof "Proof of Theorem 8.1 (click to expand)"

    Differentiability means

    \[
    g(x)-g(a)
    =C(x-a)+R(x),
    \]

    with

    \[
    \frac{\|R(x)\|}{\|x-a\|}
    \longrightarrow0
    \qquad
    (x\to a).
    \]

    Since

    \[
    \sqrt{n}(\xi_n-a)
    =O_p(1),
    \]

    we have

    \[
    \xi_n\overset{p}{\longrightarrow}a,
    \]

    and

    \[
    \sqrt{n}R(\xi_n)
    =
    \frac{\|R(\xi_n)\|}{\|\xi_n-a\|}
    \sqrt{n}\|\xi_n-a\|
    =o_p(1).
    \]

    Hence,

    \[
    \sqrt{n}\{g(\xi_n)-g(a)\}
    =
    C\sqrt{n}(\xi_n-a)
    +o_p(1).
    \]

    Slutsky's theorem proves (12). Replacing $\sqrt{n}$ by $s_n$ gives the general-rate statement in the same way. $\square$

### 8.1 Example: The Logarithm of a Poisson Mean

If

\[
X_i
\overset{\mathrm{iid}}{\sim}
\operatorname{Poisson}(\lambda),
\]

then

\[
\sqrt{n}(\overline{X}_n-\lambda)
\overset{d}{\longrightarrow}
N(0,\lambda).
\]

Take

\[
g(x)=\log x,
\qquad
g'(\lambda)=\frac{1}{\lambda}.
\]

The delta method gives

\[
\sqrt{n}
\{\log\overline{X}_n-\log\lambda\}
\overset{d}{\longrightarrow}
N(0,1/\lambda).
\]

---

## 9. Why the Identifiable Model Has No Globally Consistent Estimator

We now prove the claim about example (1) in Section 3.1.

First, $h$ is injective. It is clearly injective within the rationals and within the irrationals. If a rational number $q$ and an irrational number $u$ satisfied

\[
q=1-u,
\]

then $u=1-q$ would be rational, a contradiction.

Suppose that an estimator $T_n$ were consistent for every $\theta\in(0,1)$. Replacing it by

\[
0\vee T_n\wedge1
\]

preserves consistency, so assume

\[
0\leq T_n\leq1.
\]

For $p\in(0,1)$, define

\[
H_n(p)=E_pT_n,
\]

where the expectation is computed from $n$ independent and identically distributed $\operatorname{Bernoulli}(p)$ observations. For fixed $n$, $H_n(p)$ is a polynomial in $p$ and is therefore continuous.

Boundedness and convergence in probability imply convergence in $L^1$. Therefore:

If $p$ is rational, consistency at $\theta=p$ gives

\[
H_n(p)\longrightarrow p.
\]

If $p$ is irrational, set $\theta=1-p$. Then $\theta$ is also irrational, and the data distribution under this parameter is $\operatorname{Bernoulli}(p)$. Consistency therefore gives

\[
H_n(p)\longrightarrow1-p.
\]

Thus the continuous functions $H_n$ converge pointwise to $h$. A pointwise limit of continuous functions on an interval is a Baire class-one function, and every Baire class-one function has a dense set of continuity points.

However, for every $p\neq1/2$, rational and irrational sequences approaching $p$ give respective limits $p$ and $1-p$ for $h$. Thus $h$ is discontinuous at every $p\neq1/2$, contradicting the Baire theorem.

Therefore, no globally consistent estimator exists.

---

## 10. Bias Order for an Empirical Central Moment

Let

\[
Y_i=X_i-\mu,
\qquad
\overline{Y}_n
=\frac{1}{n}\sum_{i=1}^nY_i,
\qquad
A_{s,n}^Y
=\frac{1}{n}\sum_{i=1}^nY_i^s.
\]

Then

\[
M_{r,n}
=\sum_{j=0}^r
\binom{r}{j}
(-\overline{Y}_n)^j
A_{r-j,n}^Y.
\tag{13}
\]

The $j=0$ term has expectation $\mu_r$. For $j=1$, independence and $EY=0$ give

\[
E(\overline{Y}_nA_{r-1,n}^Y)
=\frac{\mu_r}{n}.
\]

For $j\geq2$, the Cauchy-Schwarz inequality gives

\[
\begin{aligned}
E\left|
\overline{Y}_n^jA_{r-j,n}^Y
\right|
&\leq
\{E|\overline{Y}_n|^{2j}\}^{1/2}
\{E|A_{r-j,n}^Y|^2\}^{1/2}\\
&=O(n^{-j/2})\\
&=O(n^{-1}).
\end{aligned}
\]

The moment inequalities for sums of independent centered variables and the condition

\[
E|X|^{2r}<\infty
\]

justify these bounds. Taking expectations in (13) and summing finitely many terms proves

\[
EM_{r,n}-\mu_r
=O(n^{-1}).
\]

---

## 11. Chapter Summary

!!! success "Key Conclusions"

    1. Large-sample analysis studies consistency, rate of convergence, and limiting distribution in sequence.

    2. Identifiability is necessary for consistent estimation, but stable recovery usually also requires continuity, compactness, or strong separation.

    3. Empirical raw and central moments are strongly consistent when the corresponding absolute moments exist, so continuously invertible moment estimators are also strongly consistent.

    4. Empirical central moments are generally biased, but their bias is $O(n^{-1})$ when the absolute moment of order $2r$ exists.

    5. Empirical raw moments satisfy a multivariate central limit theorem, while empirical central moments admit asymptotically linear influence-function expansions.

    6. The delta method transfers the limiting distribution of an estimator through a smooth transformation using the first derivative.
