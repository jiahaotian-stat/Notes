# Chapter 3: Completeness, Sufficiency, and Optimal Unbiased Estimation

This chapter studies how to find optimal estimators within the restricted class of unbiased estimators. The core tools are **completeness**, the **Rao-Blackwell theorem**, and the **Lehmann-Scheffé theorem**. Finally, we introduce Basu's theorem and illustrate the concrete applications of these results through Bernoulli, Poisson, Normal, and Uniform models.

---

## 1. Why Study Unbiased Estimators

Let the target parameter be $g(\theta)$, the estimator be $\delta(X)$, and the loss function be $L\{g(\theta),a\}$. The risk of the estimator is defined as

\[
R(\theta,\delta)
=E_\theta L\{g(\theta),\delta(X)\}.
\]

Under squared error loss, the risk can be decomposed as

\[
R(\theta,\delta)
=\operatorname{Var}_\theta\{\delta(X)\}
+\operatorname{Bias}_\theta(\delta)^2,
\]

where

\[
\operatorname{Bias}_\theta(\delta)
=E_\theta\delta(X)-g(\theta).
\]

### 1.1 Why an Unconstrained Uniformly Optimal Estimator Usually Does Not Exist

!!! success "Proposition 1.1 (An Unconstrained Uniformly Optimal Estimator Usually Does Not Exist)"

    Suppose the family of probability measures $\{P_\theta:\theta\in\Theta\}$ is pairwise mutually absolutely continuous, and $g(\theta)$ is not constant. Under squared error loss, there is no estimator $\delta_0$ whose risk is no greater than that of every other estimator at every parameter point.

??? proof "Proof of Proposition 1.1 (click to expand)"

    Fix $\theta_0$, and compare $\delta_0$ with the constant estimator

    \[
    \delta(X)\equiv g(\theta_0).
    \]

    The latter has zero risk at $\theta_0$. If $\delta_0$ is uniformly optimal among all estimators, then

    \[
    0\leq R(\theta_0,\delta_0)\leq0.
    \]

    Hence

    \[
    \delta_0=g(\theta_0),
    \qquad P_{\theta_0}\text{-a.s.}
    \]

    Since the $P_\theta$ are mutually absolutely continuous, the above holds almost everywhere under every $P_\theta$. Choosing another parameter value $\theta_1$, the same argument gives

    \[
    \delta_0=g(\theta_1),
    \qquad P_{\theta_1}\text{-a.s.}
    \]

    This contradicts the fact that $g$ is not constant. $\square$

To obtain meaningful comparisons, one may restrict the class of estimators, for example by considering invariant estimators, Bayes estimators, minimax estimators, or unbiased estimators. This chapter studies the last option.

### 1.2 Unbiased Estimators, MRU, and UMVU

!!! info "Definition 1.2 (Unbiased Estimator)"

    If an estimator $\delta(X)$ satisfies

    \[
    E_\theta\delta(X)=g(\theta),
    \qquad \forall\theta\in\Theta,
    \]

    then $\delta$ is called an unbiased estimator of $g(\theta)$.

!!! info "Definition 1.3 (MRU and UMVU)"

    If an unbiased estimator $\delta_0$ has risk no greater than that of any other unbiased estimator at every parameter point, then $\delta_0$ is called a **minimum-risk unbiased estimator (MRU)**.

    Under squared error loss, an MRU is called a **uniformly minimum-variance unbiased estimator (UMVU)**.

For unbiased estimators, the bias term is zero, so

\[
R(\theta,\delta)=\operatorname{Var}_\theta(\delta).
\]

The optimality here is restricted optimality: a UMVU estimator is optimal only within the class of unbiased estimators; after enlarging the comparison class, it may still be dominated by some biased estimator.

### 1.3 Does an Unbiased Estimator Exist?

Unbiased estimation is essentially an inverse problem: find a function $\delta$ such that

\[
g(\theta)
=\int\delta(x)p_\theta(x)\,d\mu(x),
\qquad \forall\theta\in\Theta.
\]

Not every target function $g$ belongs to the range of this expectation operator.

!!! success "Proposition 1.4 (Estimable Targets in Bernoulli Samples)"

    Let

    \[
    X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(p),
    \qquad 0<p<1.
    \]

    A real-valued function $g(p)$ has an unbiased estimator based on this sample if and only if $g$ is a polynomial in $p$ of degree at most $n$.

??? proof "Proof of Proposition 1.4 (click to expand)"

    Let

    \[
    T=\sum_{i=1}^nX_i.
    \]

    For any statistic $D(X)$, let

    \[
    d_t=E\{D(X)\mid T=t\}.
    \]

    Grouping sample points according to the value of $T$, we obtain

    \[
    E_pD(X)
    =\sum_{t=0}^nd_t\binom ntp^t(1-p)^{n-t}.
    \]

    The right-hand side is a polynomial of degree at most $n$, so necessity holds.

    Conversely, the Bernstein polynomials

    \[
    \binom ntp^t(1-p)^{n-t},
    \qquad t=0,\ldots,n,
    \]

    form a basis of the space of polynomials of degree at most $n$. Expanding $g$ in this basis and defining the corresponding coefficients as $d_t$, then taking $D=d_T$, yields an unbiased estimator of $g(p)$. $\square$

Therefore, when $r\leq n$, $p^r$ can be estimated unbiasedly; whereas $1/p$, $\log p$, and most nonlinear standardized shape parameters cannot be estimated unbiasedly from a Bernoulli sample of fixed size $n$.

---

## 2. Completeness: Uniqueness from Expectations

### 2.1 Complete and Boundedly Complete Statistics

!!! info "Definition 2.1 (Complete Statistic)"

    A statistic $T$ is called **complete** with respect to the model family $\mathcal{P}$ if for any function $a(T)$ satisfying integrability, whenever

    \[
    E_\theta a(T)=0,
    \qquad \forall\theta\in\Theta,
    \]

    it follows that

    \[
    a(T)=0,
    \qquad P_\theta\text{-a.s.},\quad \forall\theta\in\Theta.
    \]

!!! info "Definition 2.2 (Boundedly Complete Statistic)"

    If the above conclusion is required only for bounded functions $a$, then $T$ is called **boundedly complete**.

Completeness implies bounded completeness. It and sufficiency describe different properties: sufficiency means that outside $T$ there is no longer any parameter information; completeness means that expectations under the model family can uniquely identify functions of $T$.

!!! success "Lemma 2.3 (Uniqueness of Unbiased Functions)"

    If $T$ is complete, and $a_1(T)$ and $a_2(T)$ are both unbiased estimators of the same target function $g(\theta)$, then

    \[
    a_1(T)=a_2(T),
    \qquad P_\theta\text{-a.s.},\quad \forall\theta\in\Theta.
    \]

??? proof "Proof of Lemma 2.3 (click to expand)"

    The difference $a_1(T)-a_2(T)$ is integrable, and for all $\theta$,

    \[
    E_\theta\{a_1(T)-a_2(T)\}=0.
    \]

    By the definition of completeness, the difference is almost everywhere zero. $\square$

### 2.2 Completeness of the Binomial Family

!!! example "Example 2.4 (Binomial Completeness)"

    Let

    \[
    T\sim\operatorname{Binomial}(n,p),
    \qquad 0<p<1.
    \]

    If $E_pa(T)=0$ for all $p$, then

    \[
    0=\sum_{t=0}^na(t)\binom ntp^t(1-p)^{n-t}.
    \]

    Dividing both sides by $(1-p)^n$ and letting

    \[
    z=\frac{p}{1-p}>0,
    \]

    we obtain

    \[
    0=\sum_{t=0}^na(t)\binom ntz^t,
    \qquad \forall z>0.
    \]

    If a polynomial is identically zero on an interval, then all its coefficients are zero. Hence $a(t)=0$ for $t=0,\ldots,n$, so the Binomial family is complete.

### 2.3 Completeness of the Poisson Family

!!! example "Example 2.5 (Poisson Completeness)"

    Let

    \[
    T\sim\operatorname{Poisson}(n\lambda),
    \qquad \lambda>0.
    \]

    If $E_\lambda a(T)=0$, then

    \[
    0=e^{n\lambda}E_\lambda a(T)
    =\sum_{t=0}^{\infty}a(t)\frac{(n\lambda)^t}{t!}.
    \]

    Integrability guarantees that this power series converges absolutely for every $\lambda>0$. Since it is identically zero on an interval, all coefficients are zero, that is, $a(t)=0$ for all $t$. Hence the Poisson family is complete.

---

## 3. Complete Sufficiency of Full Exponential Families

### 3.1 Complete Sufficiency Theorem

!!! success "Theorem 3.1 (Complete Sufficient Statistic in a Full Exponential Family)"

    Suppose the model has the exponential family form

    \[
    p_\eta(x)
    =\exp\{\eta^TT(x)-A(\eta)\}h(x),
    \qquad \eta\in H\subseteq\mathbb{R}^k,
    \]

    and the natural parameter space $H$ contains a nonempty open subset of $\mathbb{R}^k$. Then $T$ is sufficient and complete for the model family indexed by this open set. Therefore, for any larger model family that contains this subfamily, $T$ remains complete.

??? proof "Proof of Theorem 3.1 (click to expand)"

    By the factorization theorem, $T$ is a sufficient statistic.

    Let $\nu$ be the image measure of the base measure $h(x)\mu(dx)$ under the mapping $x\mapsto T(x)$. Then the induced distribution of $T$ is

    \[
    P_\eta^T(dt)
    =e^{\eta^Tt-A(\eta)}\nu(dt).
    \]

    Suppose $a(T)$ is integrable and its expectation is identically zero on an open set. Choose an interior point $\eta_0$ of this open set. For any $u$ near the origin,

    \[
    \int a(t)e^{(\eta_0+u)^Tt}\nu(dt)=0.
    \]

    Define the finite signed measure

    \[
    M(dt)=a(t)e^{\eta_0^Tt}\nu(dt).
    \]

    Then the two-sided Laplace transform of $M$ satisfies

    \[
    L_M(u)=\int e^{u^Tt}M(dt)=0
    \]

    for all $u$ in some neighborhood of the origin. By the local uniqueness lemma for the Laplace transform below, $M=0$. Since $e^{\eta_0^Tt}>0$, we obtain

    \[
    a(t)=0,
    \qquad \nu\text{-a.e.}
    \]

    Every $P_\eta^T$ has a strictly positive density with respect to $\nu$, so $a(T)=0$ almost surely under every $P_\eta$. Hence $T$ is complete. $\square$

!!! success "Lemma 3.2 (Local Uniqueness of the Laplace Transform)"

    Let $M$ be a finite signed measure on $\mathbb{R}^k$ whose exponential moments are finite in some neighborhood of the origin. If

    \[
    \int e^{u^Tt}M(dt)=0
    \]

    in this neighborhood, then $M=0$.

??? proof "Proof of Lemma 3.2 (click to expand)"

    For a complex vector $z$ whose real part lies in the neighborhood where the exponential moments exist, define

    \[
    L_M(z)=\int e^{z^Tt}M(dt).
    \]

    By dominating the integrand with slightly larger real exponential moments, one can show that $L_M$ is analytic and that differentiation under the integral sign is permitted.

    Starting from $L_M$ being identically zero on a real open set, applying the identity theorem for one complex variable coordinate by coordinate extends the zero identity to a connected complex tube domain containing the imaginary axis. Therefore, for every $s\in\mathbb{R}^k$,

    \[
    \widehat M(s)
    =L_M(is)
    =\int e^{is^Tt}M(dt)
    =0.
    \]

    By uniqueness of the Fourier transform of a finite signed measure, $M=0$. $\square$

### 3.2 Why the “Open Set” Condition Is Needed

A natural parameter space with full affine span is enough to guarantee minimal sufficiency, because finitely many likelihood ratios can recover the natural statistic. However, completeness requires that the Laplace transforms agree on a neighborhood. A curved parameter space may have full affine span but contain no open set. Therefore, minimal sufficiency alone does not imply completeness.

### 3.3 Full Normal Family

!!! example "Example 3.3 (Full Normal Family)"

    Let

    \[
    X_1,\ldots,X_n
    \overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2).
    \]

    The natural parameters are

    \[
    \eta_1=\frac{\mu}{\sigma^2},
    \qquad
    \eta_2=-\frac{1}{2\sigma^2},
    \]

    and the natural parameter space is

    \[
    H=\mathbb{R}\times(-\infty,0),
    \]

    which is an open set. Therefore

    \[
    \left(\sum_{i=1}^nX_i,\sum_{i=1}^nX_i^2\right)
    \]

    is a complete sufficient statistic for $(\mu,\sigma^2)$.

    When $n\geq2$, it is in one-to-one correspondence with

    \[
    (\overline X,S^2),
    \qquad
    S^2=\frac{1}{n-1}\sum_{i=1}^n(X_i-\overline X)^2.
    \]

    Hence $(\overline X,S^2)$ is also a complete sufficient statistic.

!!! note "A Common Misunderstanding"

    When both $\mu$ and $\sigma^2$ are unknown, $(\overline X,S^2)$ remains complete for the full normal family. What may lose completeness is a curved normal submodel whose natural parameter space contains no open set. “Multiple parameters” does not mean “not complete”.

---

## 4. Rao-Blackwell Improvement

### 4.1 Rao-Blackwell Theorem under Convex Loss

!!! success "Theorem 4.1 (Rao-Blackwell Theorem under Convex Loss)"

    Let $T$ be a sufficient statistic, let $\delta(X)$ be an integrable estimator, and define

    \[
    \delta_T(T)=E_\theta\{\delta(X)\mid T\}.
    \]

    Sufficiency guarantees that the right-hand side can be chosen as a function of $T$ that does not depend on the unknown parameter $\theta$. If $a\mapsto L\{g(\theta),a\}$ is convex, then

    \[
    R(\theta,\delta_T)
    \leq R(\theta,\delta),
    \qquad \forall\theta\in\Theta.
    \]

    If $\delta$ is an unbiased estimator of $g(\theta)$, then $\delta_T$ remains unbiased.

??? proof "Proof of Theorem 4.1 (click to expand)"

    Sufficiency provides a conditional distribution kernel independent of $\theta$, so the conditional expectation can be regarded as a measurable function of the same $T$. By the conditional Jensen inequality,

    \[
    L\{g(\theta),\delta_T(T)\}
    \leq E_\theta\left[L\{g(\theta),\delta(X)\}\mid T\right].
    \]

    Taking expectations of both sides gives the risk inequality.

    In addition, by the tower property,

    \[
    E_\theta\delta_T
    =E_\theta\{E_\theta(\delta\mid T)\}
    =E_\theta\delta.
    \]

    Hence unbiasedness is preserved. $\square$

### 4.2 Exact Improvement under Squared Error Loss

!!! success "Corollary 4.2 (Exact Risk Improvement under Squared Error Loss)"

    Under squared error loss,

    \[
    R(\theta,\delta)-R(\theta,\delta_T)
    =E_\theta\{\operatorname{Var}_\theta(\delta\mid T)\}
    \geq0.
    \]

    Unless $\delta$ is already a function of $T$ (ignoring $P_\theta$-null sets), the improvement is strict.

??? proof "Proof of Corollary 4.2 (click to expand)"

    Write

    \[
    \delta-g(\theta)
    =(\delta-\delta_T)+(\delta_T-g(\theta)).
    \]

    The first term has conditional expectation zero given $T$, so after expanding the square, the expectation of the cross term is zero. Hence

    \[
    E_\theta(\delta-g(\theta))^2
    =E_\theta(\delta-\delta_T)^2
    +E_\theta(\delta_T-g(\theta))^2.
    \]

    Here

    \[
    E_\theta(\delta-\delta_T)^2
    =E_\theta\{\operatorname{Var}_\theta(\delta\mid T)\}.
    \]

    Rearranging gives the result. $\square$

This equality is exactly the conditional variance decomposition, and it also explains the connection between Rao-Blackwellization, conditional Monte Carlo, and modern variance reduction methods: it averages out extra randomness that does not affect the inferential target.

### 4.3 Rao-Blackwellization of Bernoulli Power Functions

!!! example "Example 4.3 (Bernoulli Power Functions)"

    Let

    \[
    X_1,\ldots,X_n
    \overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(p),
    \qquad
    T=\sum_{i=1}^nX_i.
    \]

    When $r\leq n$, $X_1\cdots X_r$ is an unbiased estimator of $p^r$. Given $T=t$, all permutations containing $t$ ones are equally likely, so

    \[
    E(X_1\cdots X_r\mid T=t)
    =\frac{(t)_r}{(n)_r},
    \]

    where the falling factorial is defined as

    \[
    (x)_r=x(x-1)\cdots(x-r+1).
    \]

    Therefore $(T)_r/(n)_r$ remains unbiased, and its variance is no greater than that of the original product estimator.

---

## 5. Lehmann-Scheffé Theorem

!!! success "Theorem 5.1 (Lehmann-Scheffé Theorem under Convex Loss)"

    Let $T$ be a complete sufficient statistic. If $\delta_0(T)$ is an unbiased estimator of $g(\theta)$ with finite risk, then $\delta_0$ is an MRU estimator under any loss that is convex in the action $a$.

    Under squared error loss, $\delta_0$ is the unique UMVU estimator. More generally, under strictly convex loss, $\delta_0$ is the unique MRU estimator, where uniqueness ignores sets that are null for all models.

??? proof "Proof of Theorem 5.1 (click to expand)"

    Take any unbiased estimator $D(X)$ of $g(\theta)$. Rao-Blackwellize it:

    \[
    D_T(T)=E(D\mid T).
    \]

    Then $D_T$ remains unbiased, and its risk is no greater than that of $D$. Since $D_T(T)$ and $\delta_0(T)$ are both unbiased functions of the complete statistic $T$, Lemma 2.3 gives

    \[
    D_T(T)=\delta_0(T),
    \qquad P_\theta\text{-a.s.},\quad \forall\theta.
    \]

    Therefore

    \[
    R(\theta,\delta_0)
    =R(\theta,D_T)
    \leq R(\theta,D),
    \]

    which proves MRU optimality.

    Under squared error loss, Corollary 4.2 shows that equality can occur only when $D=D_T=\delta_0$ almost everywhere. Under strictly convex loss, the equality condition in Jensen's inequality gives the same conclusion. $\square$

!!! note "Standard Construction Steps for UMVU"

    (i) Find a complete sufficient statistic $T$;

    (ii) Find any unbiased estimator $D(X)$ of the target function;

    (iii) Compute $E(D\mid T)$, or directly find an unbiased function of $T$;

    (iv) Use the Lehmann-Scheffé theorem to obtain uniqueness and UMVU optimality.

---

## 6. Four Representative UMVU Constructions

### 6.1 Bernoulli Power Functions

The statistic

\[
T=\sum_{i=1}^nX_i\sim\operatorname{Binomial}(n,p)
\]

is sufficient by the factorization theorem and complete by Example 2.4. Therefore

\[
\widehat{p^r}_{\mathrm{UMVU}}
=\frac{(T)_r}{(n)_r},
\qquad r=1,\ldots,n.
\]

### 6.2 Poisson Point Probability

Let

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}\operatorname{Poisson}(\lambda),
\qquad
T=\sum_{i=1}^nX_i.
\]

For fixed $k\geq0$, the indicator function $\mathbf{1}\{X_1=k\}$ unbiasedly estimates

\[
g_k(\lambda)
=P_\lambda(X_1=k)
=e^{-\lambda}\frac{\lambda^k}{k!}.
\]

Given $T=t$, $(X_1,\ldots,X_n)$ follows a multinomial distribution with total $t$ and cell probabilities $1/n$. Therefore

\[
\widehat g_k(T)
=\binom Tk
\left(\frac1n\right)^k
\left(1-\frac1n\right)^{T-k}
\mathbf{1}\{T\geq k\}.
\]

By Poisson completeness and the Lehmann-Scheffé theorem, $\widehat g_k(T)$ is the UMVU estimator of $g_k(\lambda)$.

### 6.3 Mean and Variance of a Normal Population

Let

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2).
\]

By Example 3.3, $(\overline X,S^2)$ is a complete sufficient statistic for the full normal family. Also,

\[
E_{\mu,\sigma^2}\overline X=\mu,
\qquad
E_{\mu,\sigma^2}S^2=\sigma^2,
\]

so $\overline X$ and $S^2$ are the UMVU estimators of $\mu$ and $\sigma^2$, respectively.

Moreover,

\[
E\overline X^2
=\mu^2+\frac{\sigma^2}{n},
\]

so

\[
\overline X^2-\frac{S^2}{n}
\]

is an unbiased estimator of $\mu^2$. It is a function of the complete sufficient statistic, so it is also the UMVU estimator of $\mu^2$.

### 6.4 Uniform One-Endpoint Model

Let

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}\operatorname{Uniform}(0,\theta),
\qquad
M=X_{(n)}.
\]

By the factorization theorem, $M$ is a sufficient statistic. If $E_\theta a(M)=0$ for every $\theta>0$, then

\[
0=\frac{n}{\theta^n}
\int_0^\theta a(m)m^{n-1}\,dm.
\]

Therefore, for each $\theta$, the indefinite integral in the above expression is zero. By differentiation in the Lebesgue sense,

\[
a(\theta)\theta^{n-1}=0
\]

almost everywhere, so $M$ is complete.

On the other hand, when $r>-n$,

\[
E_\theta M^r
=\frac{n}{n+r}\theta^r.
\]

Therefore, as long as the target function is defined,

\[
\widehat{\theta^r}_{\mathrm{UMVU}}
=\frac{n+r}{n}M^r.
\]

This nonregular example shows that complete sufficiency is not limited to exponential families with fixed support.

---

## 7. Basu's Theorem: Completeness and Ancillary Statistics

### 7.1 Ancillary Statistics

!!! info "Definition 7.1 (Ancillary Statistic)"

    If the distribution of a statistic $A(X)$ does not depend on the parameter, then $A$ is called an **ancillary statistic**.

Ancillary statistics can still reflect estimation precision, design imbalance, or sampling geometry. “Distribution independent of the parameter” does not mean “useless”.

### 7.2 Basu's Theorem

!!! success "Theorem 7.2 (Basu's Theorem)"

    If $T$ is a sufficient and boundedly complete statistic and $A$ is an ancillary statistic, then under every $P_\theta$, $T$ and $A$ are independent.

??? proof "Proof of Theorem 7.2 (click to expand)"

    Fix a measurable set $B$ in the range of $A$. By sufficiency, one can choose a version

    \[
    q_B(T)=P_\theta(A\in B\mid T)
    \]

    that is independent of $\theta$.

    By ancillarity,

    \[
    P_\theta(A\in B)=c_B
    \]

    is also independent of $\theta$. Therefore

    \[
    E_\theta\{q_B(T)-c_B\}=0,
    \qquad \forall\theta.
    \]

    The function $q_B-c_B$ is bounded. By bounded completeness,

    \[
    q_B(T)=c_B,
    \qquad P_\theta\text{-a.s.}
    \]

    For any measurable set $C$ in the range of $T$,

    \[
    \begin{aligned}
    P_\theta(A\in B,T\in C)
    &=E_\theta\left[\mathbf{1}\{T\in C\}q_B(T)\right]\\
    &=c_BP_\theta(T\in C)\\
    &=P_\theta(A\in B)P_\theta(T\in C).
    \end{aligned}
    \]

    Hence $T$ and $A$ are independent. $\square$

### 7.3 Why $\overline X$ and $S^2$ Are Independent

!!! example "Example 7.3 (Independence of the Normal Sample Mean and Sample Variance)"

    Fix $\sigma^2$, and consider the submodel indexed only by $\mu$:

    \[
    X_1,\ldots,X_n
    \overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2).
    \]

    In this case, $\overline X$ is a complete sufficient statistic for $\mu$, while

    \[
    \frac{(n-1)S^2}{\sigma^2}
    \sim\chi_{n-1}^2
    \]

    has a distribution independent of $\mu$, so $S^2$ is an ancillary statistic in this submodel. By Basu's theorem,

    \[
    \overline X\perp S^2.
    \]

    Since this conclusion holds for each fixed $\sigma^2$, it holds in the full normal family.

---

## 8. Functional Target: Estimating the Normal Distribution Function

Let

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2),
\]

where $\sigma^2$ is known, $n\geq2$, and $y\in\mathbb{R}$ is fixed. The indicator function

\[
D=\mathbf{1}\{X_1\leq y\}
\]

unbiasedly estimates

\[
F_\mu(y)
=\Phi\left(\frac{y-\mu}{\sigma}\right).
\]

Here the complete sufficient statistic is $\overline X$. By the conditional distribution formula for the joint normal distribution,

\[
X_1\mid\overline X=x
\sim N\left(x,\sigma^2\left(1-\frac1n\right)\right).
\]

Therefore, Rao-Blackwellization gives

\[
\widehat F_{\mathrm{UMVU}}(y)
=\Phi\left(
\frac{y-\overline X}
{\sigma\sqrt{1-1/n}}
\right).
\]

A similar Gaussian convolution identity gives the UMVU estimator of the population density at $y$:

\[
\widehat f_{\mathrm{UMVU}}(y)
=\frac{1}{\sigma\sqrt{1-1/n}}
\phi\left(
\frac{y-\overline X}
{\sigma\sqrt{1-1/n}}
\right).
\]

Taking the expectation of the above expression is equivalent to convolving the density of $N(\mu,\sigma^2/n)$ with the density of $N(0,\sigma^2(1-1/n))$, and the result is exactly the density of $N(\mu,\sigma^2)$ at $y$. Completeness further guarantees that this estimator is the unique UMVU estimator.

---

## 9. Uniform Two-Endpoint Model

Let

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}\operatorname{Uniform}(a,b),
\qquad a<b,\quad n\geq2,
\]

and denote

\[
U=X_{(1)},
\qquad
V=X_{(n)}.
\]

### 9.1 Complete Sufficiency of Extreme Order Statistics

!!! success "Theorem 9.1 (Complete Sufficiency of Extreme Order Statistics)"

    The statistic $(U,V)$ is a complete sufficient statistic for the parameters $(a,b)$.

??? proof "Proof of Theorem 9.1 (click to expand)"

    The joint density of the sample is

    \[
    (b-a)^{-n}\mathbf{1}\{a<U<V<b\},
    \]

    so by the factorization theorem, $(U,V)$ is sufficient.

    The joint density of $(U,V)$ is

    \[
    f_{a,b}(u,v)
    =\frac{n(n-1)(v-u)^{n-2}}{(b-a)^n}
    \mathbf{1}\{a<u<v<b\}.
    \]

    Suppose that for all $a<b$,

    \[
    E_{a,b}c(U,V)=0.
    \]

    Multiplying both sides by $(b-a)^n$, we obtain

    \[
    G(a,b)
    =\int_a^b\int_u^b
    c(u,v)n(n-1)(v-u)^{n-2}\,dv\,du
    =0.
    \]

    Let

    \[
    q(u,v)=c(u,v)n(n-1)(v-u)^{n-2}.
    \]

    Absolute integrability permits differentiation in the Lebesgue sense. For almost every $a<b$,

    \[
    \frac{\partial G}{\partial a}(a,b)
    =-\int_a^bq(a,v)\,dv,
    \]

    and

    \[
    \frac{\partial^2G}{\partial b\,\partial a}(a,b)
    =-q(a,b).
    \]

    Since $G\equiv0$, we have $q(a,b)=0$ almost everywhere. The weight $n(n-1)(b-a)^{n-2}$ is positive when $a<b$, so $c(a,b)=0$ almost everywhere. From the joint density of $(U,V)$, it follows that

    \[
    c(U,V)=0,
    \qquad P_{a,b}\text{-a.s.},
    \]

    for every parameter pair $(a,b)$. Hence $(U,V)$ is complete. $\square$

### 9.2 UMVU Estimates of the Two Endpoints

For a standard Uniform$(0,1)$ sample,

\[
EU=\frac{1}{n+1},
\qquad
EV=\frac{n}{n+1}.
\]

By a linear transformation, we obtain

\[
E_{a,b}U=\frac{na+b}{n+1},
\qquad
E_{a,b}V=\frac{a+nb}{n+1}.
\]

Solving these two linear equations gives

\[
\widehat a_{\mathrm{UMVU}}
=\frac{nU-V}{n-1},
\qquad
\widehat b_{\mathrm{UMVU}}
=\frac{nV-U}{n-1}.
\]

By Theorem 9.1 and the Lehmann-Scheffé theorem, they are the UMVU estimators of $a$ and $b$, respectively.

---

## 10. No Randomization Is Needed under Convex Loss

Suppose a randomized estimator selects an action according to the distribution $Q_x(da)$ after observing $X=x$. Define the barycenter of this randomized action as

\[
\overline\delta(x)=\int a\,Q_x(da).
\]

If the randomized estimator is unbiased, then

\[
E_\theta\overline\delta(X)
=E_\theta\int a\,Q_X(da)
=g(\theta),
\]

so $\overline\delta$ is an ordinary nonrandomized unbiased estimator.

By convexity and Jensen's inequality, for each $x$,

\[
L\{g(\theta),\overline\delta(x)\}
\leq\int L\{g(\theta),a\}\,Q_x(da).
\]

Integrating over $X$ then shows that the risk of the nonrandomized barycenter estimator is no greater than that of the original randomized estimator. Therefore, under convex loss, randomization provides no additional benefit for unbiased estimation.

This differs from discrete hypothesis testing: in discrete tests, the action space and constraint structure are different, and randomization may still be very important.

---

## 11. Proof Checklist for This Chapter

Before asserting that an estimator is UMVU, or using an independence theorem, one should check the following questions:

1. Does an unbiased estimator of the target function exist over the entire parameter space?
2. Is the statistic used for conditioning sufficient, so that $E_\theta(\delta\mid T)$ can be chosen as the same statistic independent of the parameter?
3. Is the statistic complete, or merely minimal sufficient?
4. When using completeness of exponential families, does the natural parameter space contain an open set, rather than merely having full affine span?
5. Do the integrability conditions required for expectations, power series, or Laplace transforms hold?
6. Is the loss function used in the Rao-Blackwell theorem convex? Only after excluding the equality case can one claim that the improvement is strict.
7. Before invoking the Lehmann-Scheffé theorem, is the candidate estimator both unbiased and a function of a complete sufficient statistic?
8. When using Basu's theorem, is the second statistic truly ancillary for the same parameter family, and is $T$ boundedly complete for that parameter family?
9. Is the final conclusion only optimality within the class of unbiased estimators, or does it incorrectly claim admissibility among all estimators?

---

## 12. Summary of This Chapter

The logic of this chapter can be summarized as:

\[
\text{Sufficiency}
+\text{Rao-Blackwell improvement}
+\text{uniqueness guaranteed by completeness}
\Longrightarrow
\text{Lehmann-Scheffé optimality}.
\]

A complete sufficient statistic provides both a dimension-reduction structure and uniqueness of unbiased functions. The Rao-Blackwell theorem improves any unbiased estimator into a function of a sufficient statistic, and the Lehmann-Scheffé theorem further identifies this function as the unique UMVU estimator. Basu's theorem reveals the independence relationship between a complete sufficient statistic and an ancillary statistic.