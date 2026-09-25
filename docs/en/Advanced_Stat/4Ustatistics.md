# Chapter 4: MVUE, U-Statistics, and Information Lower Bounds

This chapter first discusses minimum-risk unbiased estimation and minimum-variance unbiased estimation, and gives the covariance orthogonality characterization of the MVUE. It then extends the conclusions to vector targets and introduces U-statistics, a systematic method for constructing unbiased estimators. On this basis, we derive the scalar and multiparameter Cramér-Rao lower bounds, study the equality conditions and attainability of the lower bounds, and finally introduce the Hammersley-Chapman-Robbins bound, which applies to nonregular models.

---

## 1. MRU, MVUE, and the Comparison Class

Let the parameter be $\theta\in\Theta$, where possibly $\Theta\subseteq\mathbb{R}^k$, and let the scalar target to be estimated be $g(\theta)$. Under squared error loss,

\[
L\{g(\theta),a\}=\{a-g(\theta)\}^2.
\]

The risk of an estimator $\delta$ can be decomposed as

\[
R(\theta,\delta)
=\operatorname{Var}_\theta(\delta)
+\operatorname{Bias}_\theta(\delta)^2.
\]

Therefore, within the class of unbiased estimators,

\[
R(\theta,\delta)=\operatorname{Var}_\theta(\delta).
\]

!!! info "Definition 1.1 (MRU and MVUE)"

    If an unbiased estimator $\widehat g$ has risk, under a given loss function, no greater than that of any other unbiased estimator for every $\theta\in\Theta$, then $\widehat g$ is called a **minimum-risk unbiased estimator**, abbreviated MRU.

    Under squared error loss, an MRU is called a **minimum-variance unbiased estimator**, abbreviated MVUE, and is also often called a UMVU estimator.

Under squared loss, MRU and MVUE are completely equivalent. Note that the dimension of the parameter $\theta$ does not affect this conclusion: even if the parameter belongs to a multiparameter model, as long as the target $g(\theta)$ is scalar, one can still compare unbiased estimators using scalar variance.

---

## 2. Orthogonality Characterization of the MVUE

Let

\[
\mathcal{Z}
=\left\{Z:E_\theta Z=0,\ \forall\theta\in\Theta\right\}
\]

denote the linear space consisting of estimators that have mean zero over the entire parameter space.

!!! success "Theorem 2.1 (Covariance Characterization)"

    Let $\widehat g$ be an unbiased estimator of $g(\theta)$ with finite variance for every $\theta$. Then $\widehat g$ is an MVUE if and only if, for every $Z\in\mathcal{Z}$ and every $\theta_0$ at which $Z$ has finite variance,

    \[
    \operatorname{Cov}_{\theta_0}(\widehat g,Z)=0.
    \]

??? proof "Proof of Theorem 2.1 (click to expand)"

    First suppose $\widehat g$ is an MVUE. For any $Z\in\mathcal{Z}$ and any $a\in\mathbb{R}$, the perturbed estimator

    \[
    \widehat g+aZ
    \]

    is still an unbiased estimator of $g(\theta)$. At a fixed $\theta_0$,

    \[
    \begin{aligned}
    &\operatorname{Var}_{\theta_0}(\widehat g+aZ)
    -\operatorname{Var}_{\theta_0}(\widehat g)\\
    &\qquad
    =a^2\operatorname{Var}_{\theta_0}(Z)
    +2a\operatorname{Cov}_{\theta_0}(\widehat g,Z).
    \end{aligned}
    \]

    By the optimality of the MVUE, the above expression is nonnegative for every $a\in\mathbb{R}$. If the covariance were nonzero, choosing $a$ with the opposite sign and sufficiently small absolute value would make the right-hand side negative, a contradiction. Therefore

    \[
    \operatorname{Cov}_{\theta_0}(\widehat g,Z)=0.
    \]

    Conversely, suppose the above orthogonality condition holds. Let $D$ be any other finite-variance unbiased estimator of $g(\theta)$, and let

    \[
    Z=D-\widehat g.
    \]

    Then $Z\in\mathcal{Z}$, so

    \[
    \begin{aligned}
    \operatorname{Var}_\theta(D)
    &=\operatorname{Var}_\theta(\widehat g+Z)\\
    &=\operatorname{Var}_\theta(\widehat g)
      +\operatorname{Var}_\theta(Z)
      +2\operatorname{Cov}_\theta(\widehat g,Z)\\
    &=\operatorname{Var}_\theta(\widehat g)
      +\operatorname{Var}_\theta(Z)\\
    &\geq\operatorname{Var}_\theta(\widehat g).
    \end{aligned}
    \]

    Hence $\widehat g$ is an MVUE. $\square$

!!! note "Geometric Interpretation"

    At each fixed parameter value, the square-integrable statistics form a Hilbert space. The MVUE is orthogonal to all admissible perturbations whose mean is zero over the entire model.

    The phrase “over the entire model” is very important: if a statistic has mean zero only at a single parameter value, adding it to the original estimator generally does not preserve global unbiasedness.

### 2.1 Complete Sufficient Statistics and Orthogonality

!!! success "Corollary 2.2 (Complete Sufficiency Implies Orthogonality)"

    Let $T$ be a complete sufficient statistic. If $\widehat g=a(T)$ is an unbiased estimator of $g(\theta)$, then $\widehat g$ satisfies the orthogonality condition of Theorem 2.1 and is therefore an MVUE.

??? proof "Proof of Corollary 2.2 (click to expand)"

    For any $Z\in\mathcal{Z}$, by sufficiency, one can choose a common version

    \[
    q(T)=E_\theta(Z\mid T)
    \]

    that does not depend on $\theta$.

    By the law of iterated expectations,

    \[
    E_\theta q(T)=E_\theta Z=0,
    \qquad \forall\theta\in\Theta.
    \]

    By completeness, $q(T)=0$ almost everywhere. Therefore

    \[
    \begin{aligned}
    \operatorname{Cov}_\theta\{a(T),Z\}
    &=E_\theta\{a(T)Z\}\\
    &=E_\theta\left[a(T)E_\theta(Z\mid T)\right]\\
    &=E_\theta\{a(T)q(T)\}=0.
    \end{aligned}
    \]

    Then by Theorem 2.1, $a(T)$ is an MVUE. This is exactly the orthogonality form of the Lehmann-Scheffé theorem. $\square$

---

## 3. Vector Targets and Multiparameter Models

One must distinguish the dimension of the parameter from the dimension of the target:

\[
\theta\in\mathbb{R}^k,
\qquad
g(\theta)\in\mathbb{R}^r.
\]

The covariance characterization in the previous section already applies to scalar targets under any $k$. Now consider

\[
\widehat g
=(\widehat g_1,\ldots,\widehat g_r)^\mathsf{T}
\]

to be an unbiased estimator of

\[
g=(g_1,\ldots,g_r)^\mathsf{T}.
\]

!!! info "Definition 3.1 (Vector MVUE)"

    If each component $\widehat g_j$ is a scalar MVUE of $g_j$, then $\widehat g$ is called a vector MVUE.

For symmetric matrices, write

\[
A\preceq B
\]

to mean that $B-A$ is positive semidefinite.

!!! success "Theorem 3.2 (Equivalent Characterizations of a Vector MVUE)"

    Assume that all relevant second moments are finite. Then the following statements are equivalent:

    1. Each $\widehat g_j$ is an MVUE of $g_j$;
    2. For every $c\in\mathbb{R}^r$, $c^\mathsf{T}\widehat g$ is an MVUE of $c^\mathsf{T}g(\theta)$;
    3. For every unbiased vector estimator $D$ of $g$,

    \[
    \operatorname{Cov}_\theta(\widehat g)
    \preceq
    \operatorname{Cov}_\theta(D),
    \qquad \forall\theta.
    \]

??? proof "Proof of Theorem 3.2 (click to expand)"

    Suppose condition 1 holds. For any $Z\in\mathcal{Z}$, by Theorem 2.1,

    \[
    \operatorname{Cov}_\theta(\widehat g_j,Z)=0,
    \qquad j=1,\ldots,r.
    \]

    Hence for any $c\in\mathbb{R}^r$,

    \[
    \operatorname{Cov}_\theta(c^\mathsf{T}\widehat g,Z)
    =\sum_{j=1}^r c_j
    \operatorname{Cov}_\theta(\widehat g_j,Z)
    =0.
    \]

    Applying Theorem 2.1 again yields condition 2.

    If condition 2 holds and $D$ is any unbiased vector estimator, then for every $c$,

    \[
    \begin{aligned}
    c^\mathsf{T}\operatorname{Cov}_\theta(\widehat g)c
    &=\operatorname{Var}_\theta(c^\mathsf{T}\widehat g)\\
    &\leq\operatorname{Var}_\theta(c^\mathsf{T}D)\\
    &=c^\mathsf{T}\operatorname{Cov}_\theta(D)c.
    \end{aligned}
    \]

    This is exactly the covariance matrix inequality in the positive semidefinite order, so condition 3 holds.

    Finally suppose condition 3 holds. Fix $j$, and let $d$ be any scalar unbiased estimator of $g_j$. Replace only the $j$-th component of $\widehat g$ with $d$ to obtain an unbiased vector estimator $D$. Taking the quadratic form of the covariance matrix inequality on both sides with respect to the $j$-th coordinate vector $e_j$ gives

    \[
    \operatorname{Var}_\theta(\widehat g_j)
    \leq\operatorname{Var}_\theta(d).
    \]

    Hence each component is an MVUE, that is, condition 1 holds. $\square$

!!! success "Corollary 3.3 (Vector Orthogonality)"

    A vector estimator $\widehat g$ is a vector MVUE if and only if, for every scalar $Z\in\mathcal{Z}$ with finite variance,

    \[
    \operatorname{Cov}_\theta(\widehat g,Z)
    =0\in\mathbb{R}^r.
    \]

For Euclidean squared loss, within the class of unbiased estimators,

\[
E_\theta\|\widehat g-g(\theta)\|_2^2
=\operatorname{tr}\left\{
\operatorname{Cov}_\theta(\widehat g)
\right\}.
\]

However, comparison of covariance matrices in the positive semidefinite order is stronger than comparison of traces, because it controls the variance of every linear combination $c^\mathsf{T}\widehat g$.

!!! example "Example 3.4 (Mean and Variance of a Normal Population)"

    Let

    \[
    X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}
    N(\mu,\sigma^2).
    \]

    From the previous chapter, $(\overline X,S^2)$ is a complete sufficient statistic for $(\mu,\sigma^2)$. Hence it is a vector MVUE, and

    \[
    \operatorname{Cov}_{\mu,\sigma^2}
    \begin{pmatrix}
    \overline X\\
    S^2
    \end{pmatrix}
    =
    \begin{pmatrix}
    \sigma^2/n & 0\\
    0 & 2\sigma^4/(n-1)
    \end{pmatrix}.
    \]

    The off-diagonal entries are zero because

    \[
    \overline X\perp S^2.
    \]

    However, being a vector MVUE does not mean that the multiparameter Cramér-Rao lower bound is necessarily attainable; this example will be discussed again later.

---

## 4. U-Statistics: Systematic Construction of Unbiased Estimators

Let

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}F.
\]

If the target functional can be written as

\[
q(F)=E_Fh(X_1,\ldots,X_m),
\qquad m\leq n,
\]

where the kernel $h$ is integrable, then one can use $h$ to construct an unbiased estimator.

### 4.1 Symmetrization of the Kernel

!!! success "Proposition 4.1 (Symmetrization Does Not Increase Variance)"

    Define the symmetrized kernel

    \[
    h_s(x_1,\ldots,x_m)
    =\frac{1}{m!}\sum_\pi
    h(x_{\pi(1)},\ldots,x_{\pi(m)}),
    \]

    where the sum runs over all permutations $\pi$. Then

    \[
    E_Fh_s=E_Fh.
    \]

    When second moments exist, we also have

    \[
    \operatorname{Var}_F\{h_s(X_1,\ldots,X_m)\}
    \leq
    \operatorname{Var}_F\{h(X_1,\ldots,X_m)\}.
    \]

??? proof "Proof of Proposition 4.1 (click to expand)"

    Let $\Pi$ be a uniform random permutation independent of the sample. By exchangeability,

    \[
    h(X_{\Pi(1)},\ldots,X_{\Pi(m)})
    \]

    has the same distribution as $h(X_1,\ldots,X_m)$. Conditional on the sample,

    \[
    \begin{aligned}
    h_s(X_1,\ldots,X_m)
    =E\bigl[&h(X_{\Pi(1)},\ldots,X_{\Pi(m)})\\
    &\mid X_1,\ldots,X_m\bigr].
    \end{aligned}
    \]

    The law of iterated expectations gives equality of expectations, and the conditional variance decomposition gives the variance inequality. $\square$

Therefore, in what follows we may assume without loss of generality that the kernel $h$ is symmetric.

### 4.2 Definition and Unbiasedness of U-Statistics

!!! info "Definition 4.2 (U-Statistic)"

    The U-statistic generated by an $m$-th order symmetric kernel $h$ is defined as

    \[
    U_n
    =\binom{n}{m}^{-1}
    \sum_{1\leq i_1<\cdots<i_m\leq n}
    h(X_{i_1},\ldots,X_{i_m}).
    \]

!!! success "Proposition 4.3 (Unbiasedness of U-Statistics)"

    If $q(F)=E_Fh(X_1,\ldots,X_m)$ is finite, then

    \[
    E_FU_n=q(F).
    \]

Each summand has the same expectation as $h(X_1,\ldots,X_m)$, and there are $\binom{n}{m}$ summands, so the conclusion follows immediately.

### 4.3 Rao-Blackwell Representation of U-Statistics

!!! success "Theorem 4.4 (U-Statistics as Rao-Blackwell Estimators)"

    Suppose $F$ ranges over all distributions that have a density with respect to some fixed nonatomic measure, and suppose $h$ has finite second moment. Then

    \[
    U_n
    =E_F\left\{
    h(X_1,\ldots,X_m)
    \mid X_{(1)},\ldots,X_{(n)}
    \right\}.
    \]

    Therefore the variance of $U_n$ is no greater than that of the crude estimator $h(X_1,\ldots,X_m)$. In this completely nonparametric model, $U_n$ is the unique MVUE of $q(F)$.

??? proof "Proof of Theorem 4.4 (click to expand)"

    Conditional on the order statistics, the $n!$ ways of assigning the observed values to the original labels are equally likely. Taking the conditional average of $h$ over these assignments is equivalent to averaging equally over all $m$-element subsets. Since $h$ is symmetric, the result is exactly $U_n$.

    Hence $U_n$ is the conditional expectation of $h(X_1,\ldots,X_m)$ given the order statistics. The Rao-Blackwell theorem gives the variance inequality.

    Appendix A proves that the vector of order statistics is a complete sufficient statistic in this completely nonparametric model. Since $U_n$ is a symmetric function of it and is unbiased, by the Lehmann-Scheffé theorem, $U_n$ is the unique MVUE. $\square$

### 4.4 Typical Examples

!!! example "Example 4.5 (Estimating $\mu^2$)"

    Let

    \[
    \mu=E_FX,
    \qquad E_FX^2<\infty.
    \]

    Because

    \[
    E_F(X_1X_2)=\mu^2,
    \]

    take the second-order kernel $h(x_1,x_2)=x_1x_2$. The corresponding U-statistic is

    \[
    \begin{aligned}
    U_n
    &=\binom{n}{2}^{-1}\sum_{i<j}X_iX_j\\
    &=\frac{\left(\sum_{i=1}^nX_i\right)^2
    -\sum_{i=1}^nX_i^2}{n(n-1)}.
    \end{aligned}
    \]

    It is the nonparametric MVUE of $\mu^2$.

!!! example "Example 4.6 (Estimating the Variance)"

    Take the symmetric kernel

    \[
    h(x_1,x_2)=\frac{(x_1-x_2)^2}{2}.
    \]

    Then

    \[
    E_Fh(X_1,X_2)=\operatorname{Var}_F(X).
    \]

    Using the identity

    \[
    \sum_{i<j}(X_i-X_j)^2
    =n\sum_{i=1}^n(X_i-\overline X)^2,
    \]

    we obtain

    \[
    \begin{aligned}
    U_n
    &=\binom{n}{2}^{-1}
    \sum_{i<j}\frac{(X_i-X_j)^2}{2}\\
    &=\frac{1}{n-1}
    \sum_{i=1}^n(X_i-\overline X)^2.
    \end{aligned}
    \]

    Therefore the usual unbiased sample variance is a second-order U-statistic. The degree of the variance functional is $2$; Appendix B will show that in a completely nonparametric model there is no unbiased estimator of the variance that uses only one observation.

!!! example "Example 4.7 (Estimating the Distribution Function)"

    For fixed $t$, take the first-order kernel

    \[
    h_t(x)=\mathbf{1}_{\{x\leq t\}}.
    \]

    The corresponding first-order U-statistic is

    \[
    F_n(t)=\frac{1}{n}\sum_{i=1}^n
    \mathbf{1}_{\{X_i\leq t\}},
    \]

    that is, the value of the empirical distribution function at $t$.

U-statistics are also widely used in rank tests, Kendall's tau, kernel two-sample statistics, and various modern pairwise estimators.

---

## 5. Scalar Cramér-Rao Inequality

Let $p_\theta$ be a density with respect to a common dominating measure, with the parameter $\theta$ lying in some open interval. Consider the following regularity conditions:

!!! note "Regularity Conditions"

    **(R1)** The support of the distribution does not depend on $\theta$, and $p_\theta>0$ on that support.

    **(R2)** $p_\theta(x)$ is differentiable with respect to $\theta$, and the integrals of $p_\theta$ and $\delta p_\theta$ permit interchange of differentiation and integration.

    **(R3)** The score function

    \[
    S_\theta(X)
    =\frac{\partial}{\partial\theta}
    \log p_\theta(X)
    \]

    is square integrable, and the Fisher information satisfies

    \[
    0<I(\theta)
    =E_\theta S_\theta^2<\infty.
    \]

### 5.1 Score Identity

!!! success "Lemma 5.1 (Score Identity)"

    Under conditions (R1)-(R3),

    \[
    E_\theta S_\theta=0.
    \]

    If $\delta$ is an unbiased estimator of a differentiable target $g(\theta)$, and condition (R2) also applies to $\delta$, then

    \[
    \operatorname{Cov}_\theta(\delta,S_\theta)
    =g'(\theta).
    \]

??? proof "Proof of Lemma 5.1 (click to expand)"

    Because

    \[
    S_\theta
    =\frac{\partial p_\theta/\partial\theta}{p_\theta},
    \]

    we have

    \[
    \begin{aligned}
    E_\theta S_\theta
    &=\int\frac{\partial p_\theta(x)}{\partial\theta}
    \,d\mu(x)\\
    &=\frac{d}{d\theta}\int p_\theta(x)\,d\mu(x)=0.
    \end{aligned}
    \]

    Similarly, from $E_\theta\delta=g(\theta)$,

    \[
    \begin{aligned}
    g'(\theta)
    &=\frac{d}{d\theta}
    \int\delta(x)p_\theta(x)\,d\mu(x)\\
    &=\int\delta(x)S_\theta(x)p_\theta(x)\,d\mu(x)\\
    &=E_\theta(\delta S_\theta).
    \end{aligned}
    \]

    And since $E_\theta S_\theta=0$, the last expression is exactly $\operatorname{Cov}_\theta(\delta,S_\theta)$. $\square$

### 5.2 Cramér-Rao Lower Bound

!!! success "Theorem 5.2 (Scalar Cramér-Rao Lower Bound)"

    Under the conditions of Lemma 5.1, every unbiased estimator $\delta$ of $g(\theta)$ satisfies

    \[
    \operatorname{Var}_\theta(\delta)
    \geq
    \frac{\{g'(\theta)\}^2}{I(\theta)}.
    \]

    For an independent and identically distributed sample of size $n$,

    \[
    I_n(\theta)=nI_1(\theta).
    \]

??? proof "Proof of Theorem 5.2 (click to expand)"

    By the Cauchy-Schwarz inequality and the score identity,

    \[
    \begin{aligned}
    \{g'(\theta)\}^2
    &=\operatorname{Cov}_\theta(\delta,S_\theta)^2\\
    &\leq\operatorname{Var}_\theta(\delta)
    \operatorname{Var}_\theta(S_\theta)\\
    &=\operatorname{Var}_\theta(\delta)I(\theta).
    \end{aligned}
    \]

    Rearranging gives the conclusion. For independent samples, the sample score is the sum of the individual observation scores. The terms are independent and have mean zero, so the variance of the total score is $nI_1(\theta)$. $\square$

!!! success "Lemma 5.3 (Expected Hessian Identity)"

    If $p_\theta$ is twice differentiable and the second derivative also permits interchange with integration, then

    \[
    I(\theta)
    =-E_\theta\left\{
    \frac{\partial^2}{\partial\theta^2}
    \log p_\theta(X)
    \right\}.
    \]

??? proof "Proof of Lemma 5.3 (click to expand)"

    Direct differentiation gives

    \[
    \frac{\partial^2}{\partial\theta^2}\log p_\theta
    =\frac{p_\theta''}{p_\theta}
    -\left(\frac{p_\theta'}{p_\theta}\right)^2.
    \]

    Taking expectations, the first term is

    \[
    \int p_\theta''\,d\mu
    =\frac{d^2}{d\theta^2}1=0,
    \]

    while the second term is $-I(\theta)$, so the conclusion holds. $\square$

!!! info "Definition 5.4 (Finite-Sample Efficiency)"

    An unbiased estimator is called **efficient** at a parameter value $\theta$ if it makes the Cramér-Rao inequality an equality at $\theta$. If it attains equality for every $\theta$ in the model, it is called **globally efficient**.

Every globally efficient unbiased estimator is an MVUE, because the Cramér-Rao lower bound applies to all unbiased competitors. But the converse does not hold: an MVUE need not attain the Cramér-Rao lower bound.

---

## 6. Equality Conditions and Attainability of the Scalar Lower Bound

!!! success "Theorem 6.1 (Equality Condition for the Cramér-Rao Lower Bound)"

    Under the scalar Cramér-Rao regularity conditions, the inequality becomes an equality at $\theta$ if and only if

    \[
    \delta(X)-g(\theta)
    =\frac{g'(\theta)}{I(\theta)}S_\theta(X),
    \qquad P_\theta\text{-a.s.}
    \]

??? proof "Proof of Theorem 6.1 (click to expand)"

    The Cauchy-Schwarz inequality becomes an equality if and only if the two centered random variables are linearly related. Therefore

    \[
    \delta-g(\theta)=a(\theta)S_\theta
    \]

    almost surely. Taking the covariance of both sides with the score function and using the score identity gives

    \[
    g'(\theta)=a(\theta)I(\theta).
    \]

    Therefore

    \[
    a(\theta)=\frac{g'(\theta)}{I(\theta)}.
    \]

    Conversely, if the equality condition holds, then

    \[
    \operatorname{Var}_\theta(\delta)
    =\frac{\{g'(\theta)\}^2}{I(\theta)^2}
    \operatorname{Var}_\theta(S_\theta)
    =\frac{\{g'(\theta)\}^2}{I(\theta)}.
    \]

    Hence the Cramér-Rao lower bound is attained. $\square$

!!! warning "The Attainability Condition Is Very Strict"

    The left-hand side of the equality condition must be a statistic that does not depend on the unknown parameter, while the right-hand side contains the score function, which depends on the parameter. Hence unbiasedness alone is not sufficient to guarantee attainability of the lower bound.

### 6.1 Global Attainability and Exponential Families

!!! success "Theorem 6.2 (Global Equality Implies Exponential Family Structure)"

    Suppose the regularity conditions hold on some interval, $I(\theta)>0$, and $g'(\theta)\neq0$. If there exists a single estimator $\delta$ that attains the Cramér-Rao lower bound at every $\theta$, then the model must have the form

    \[
    p_\theta(x)
    =\exp\{q(\theta)\delta(x)+r(\theta)+c(x)\}.
    \]

    Hence global attainability requires the model to be a one-parameter exponential family with the efficient estimator as the natural statistic.

??? proof "Proof of Theorem 6.2 (click to expand)"

    Rewrite the equality condition as

    \[
    \frac{\partial}{\partial\theta}
    \log p_\theta(x)
    =\frac{I(\theta)}{g'(\theta)}\delta(x)
    -\frac{I(\theta)g(\theta)}{g'(\theta)}.
    \]

    Let the two coefficients on the right-hand side that depend only on the parameter be $q'(\theta)$ and $r'(\theta)$, respectively, and integrate with respect to $\theta$ to obtain

    \[
    \log p_\theta(x)
    =q(\theta)\delta(x)+r(\theta)+c(x).
    \]

    Exponentiating yields the conclusion. $\square$

The above theorem is only a necessary result. For regular exponential families, one can further precisely characterize which targets can attain the lower bound.

!!! success "Theorem 6.3 (Efficient Affine Statistics in One-Parameter Exponential Families)"

    Let the regular one-parameter exponential family be

    \[
    p_\theta(x)
    =C(\theta)\exp\{Q(\theta)T(x)\}h(x),
    \qquad Q'(\theta)\neq0,
    \]

    with

    \[
    0<\operatorname{Var}_\theta(T)<\infty.
    \]

    For any constants $a,b$, the estimator

    \[
    \delta=aT+b
    \]

    is globally efficient for the target

    \[
    g(\theta)=E_\theta(aT+b).
    \]

    Conversely, under the condition that $T$ is nondegenerate, any single estimator that is efficient throughout the entire family must have the affine form $aT+b$.

??? proof "Proof of Theorem 6.3 (click to expand)"

    Write

    \[
    m(\theta)=E_\theta T.
    \]

    Differentiating the normalization identity gives

    \[
    \frac{C'(\theta)}{C(\theta)}
    =-Q'(\theta)m(\theta).
    \]

    Therefore

    \[
    S_\theta
    =Q'(\theta)\{T-m(\theta)\},
    \]

    and

    \[
    I(\theta)
    =Q'(\theta)^2\operatorname{Var}_\theta(T).
    \]

    Applying the score identity to $T$ gives

    \[
    m'(\theta)
    =Q'(\theta)\operatorname{Var}_\theta(T).
    \]

    For $g(\theta)=am(\theta)+b$, we have

    \[
    g'(\theta)
    =aQ'(\theta)\operatorname{Var}_\theta(T).
    \]

    Hence

    \[
    \frac{g'(\theta)}{I(\theta)}S_\theta
    =a\{T-m(\theta)\}
    =\delta-g(\theta).
    \]

    By Theorem 6.1, $\delta$ is globally efficient.

    Conversely, by the equality condition,

    \[
    \delta(x)-g(\theta)
    =\frac{g'(\theta)Q'(\theta)}{I(\theta)}
    \{T(x)-m(\theta)\}.
    \]

    Therefore for each $\theta$, $\delta$ is an affine function of $T$. Comparing the expressions at two different parameter values and using the nondegeneracy of $T$, the slope and intercept must be constants that do not depend on the parameter. $\square$

!!! warning "Exponential Family Structure Is Still Not Enough to Guarantee Efficiency for Arbitrary Targets"

    Even within an exponential family, a nonlinear target need not have an efficient unbiased estimator. The exponential family structure only makes equality possible; whether the lower bound is actually attainable is determined by the affine compatibility between the target and the natural statistic.

!!! example "Example 6.4 (UMVU but Not Efficient)"

    Let

    \[
    X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}
    N(\mu,\sigma^2),
    \]

    where $\sigma^2$ is known, and write

    \[
    v=\frac{\sigma^2}{n}.
    \]

    The statistic

    \[
    \delta=\overline X^2-v
    \]

    is the UMVU estimator of $\mu^2$. Since $\overline X\sim N(\mu,v)$,

    \[
    \operatorname{Var}_\mu(\delta)
    =4\mu^2v+2v^2.
    \]

    The corresponding Cramér-Rao lower bound is

    \[
    \frac{\{(\mu^2)'\}^2}{n/\sigma^2}
    =4\mu^2v.
    \]

    There is a strictly positive difference $2v^2$, so the lower bound for this target is not attainable. This example shows that MVUE does not imply efficiency.

---

## 7. Multiparameter Cramér-Rao Inequality

Let

\[
\theta=(\theta_1,\ldots,\theta_k)^\mathsf{T}
\]

belong to an open subset of $\mathbb{R}^k$. Define the score vector and Fisher information matrix as

\[
S_\theta=\nabla_\theta\log p_\theta(X),
\qquad
I(\theta)=E_\theta(S_\theta S_\theta^\mathsf{T}),
\]

and assume that $I(\theta)$ is positive definite. For a differentiable target $g:\mathbb{R}^k\to\mathbb{R}^r$, denote its Jacobian matrix by

\[
G(\theta)
=\frac{\partial g(\theta)}{\partial\theta^\mathsf{T}}
\in\mathbb{R}^{r\times k}.
\]

!!! success "Theorem 7.1 (Multiparameter Cramér-Rao Lower Bound)"

    When the vector versions of (R1)-(R3) hold, every unbiased vector estimator $\delta\in\mathbb{R}^r$ of $g(\theta)$ satisfies

    \[
    \operatorname{Cov}_\theta(\delta)
    \succeq
    G(\theta)I(\theta)^{-1}G(\theta)^\mathsf{T}.
    \]

    Equality holds if and only if

    \[
    \delta-g(\theta)
    =G(\theta)I(\theta)^{-1}S_\theta,
    \qquad P_\theta\text{-a.s.}
    \]

??? proof "Proof of Theorem 7.1 (click to expand)"

    Differentiating each unbiasedness identity gives the score identity in matrix form:

    \[
    E_\theta S_\theta=0,
    \qquad
    \operatorname{Cov}_\theta(\delta,S_\theta)=G(\theta).
    \]

    Define the residual

    \[
    R
    =\delta-g(\theta)-G(\theta)I(\theta)^{-1}S_\theta.
    \]

    Expanding its covariance matrix directly gives

    \[
    \begin{aligned}
    \operatorname{Cov}_\theta(R)
    &=\operatorname{Cov}_\theta(\delta)
    -GI^{-1}\operatorname{Cov}_\theta(S_\theta,\delta)\\
    &\quad
    -\operatorname{Cov}_\theta(\delta,S_\theta)I^{-1}G^\mathsf{T}
    +GI^{-1}II^{-1}G^\mathsf{T}\\
    &=\operatorname{Cov}_\theta(\delta)
    -GI^{-1}G^\mathsf{T}.
    \end{aligned}
    \]

    Any covariance matrix is positive semidefinite, so the lower bound follows. Equality holds if and only if $\operatorname{Cov}_\theta(R)=0$, and since $ER=0$, this is equivalent to $R=0$ almost surely. $\square$

### 7.1 Scalar Targets with Nuisance Parameters

If the target is scalar, that is, $r=1$, then

\[
G(\theta)=\nabla g(\theta)^\mathsf{T},
\]

so

\[
\operatorname{Var}_\theta(\delta)
\geq
\nabla g(\theta)^\mathsf{T}
I(\theta)^{-1}
\nabla g(\theta).
\]

Here one must use the full Fisher information matrix including all parameters, especially the nuisance parameters.

### 7.2 Natural Exponential Families

!!! success "Corollary 7.2 (Efficient Estimators in Natural Exponential Families)"

    Let the regular natural exponential family be

    \[
    p_\eta(x)
    =\exp\{\eta^\mathsf{T}T(x)-A(\eta)\}h(x),
    \qquad \eta\in\mathcal{H}\subseteq\mathbb{R}^k,
    \]

    with

    \[
    I(\eta)=\nabla^2A(\eta)
    \]

    positive definite. For a constant matrix $C\in\mathbb{R}^{r\times k}$ and a constant vector $b\in\mathbb{R}^r$, the estimator

    \[
    \delta=CT+b
    \]

    is efficient for the target

    \[
    g(\eta)
    =CE_\eta T+b
    =C\nabla A(\eta)+b.
    \]

    If the affine support of $T$ is $\mathbb{R}^k$, then these are exactly the targets that have globally efficient unbiased estimators.

??? proof "Proof of Corollary 7.2 (click to expand)"

    In a natural exponential family,

    \[
    S_\eta=T-\nabla A(\eta),
    \]

    and

    \[
    I(\eta)
    =\operatorname{Cov}_\eta(T)
    =\nabla^2A(\eta).
    \]

    For the target $g(\eta)=C\nabla A(\eta)+b$, we have

    \[
    G= C\nabla^2A=CI.
    \]

    Therefore

    \[
    GI^{-1}S_\eta
    =C\{T-E_\eta T\}
    =\delta-g(\eta),
    \]

    which satisfies the equality condition of the multiparameter lower bound.

    Conversely, the equality condition implies that for each $\eta$, the same estimator is an affine function of $T$, with coefficients $G(\eta)I(\eta)^{-1}$. Comparing the affine representations at two parameter values and using the full affine support of $T$, the coefficients and intercept must be constants. Hence $\delta=CT+b$, and taking expectations yields the corresponding target. $\square$

### 7.3 Normal Mean and Variance: MVUE but Not Fully Efficient

Let

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}N(\mu,\tau),
\qquad \tau=\sigma^2,
\]

where both $\mu$ and $\tau$ are unknown. The two components of the score function are

\[
S_\mu
=\frac{1}{\tau}\sum_{i=1}^n(X_i-\mu)
\]

and

\[
S_\tau
=-\frac{n}{2\tau}
+\frac{1}{2\tau^2}
\sum_{i=1}^n(X_i-\mu)^2.
\]

The third central moment of the normal distribution is zero, so the two score components are uncorrelated. Also, since

\[
\operatorname{Var}\{(X_i-\mu)^2\}=2\tau^2,
\]

the Fisher information matrix is

\[
I(\mu,\tau)
=
\begin{pmatrix}
n/\tau & 0\\
0 & n/(2\tau^2)
\end{pmatrix}.
\]

For the identity target

\[
g(\mu,\tau)=(\mu,\tau)^\mathsf{T},
\]

the matrix lower bound is

\[
I(\mu,\tau)^{-1}
=
\begin{pmatrix}
\tau/n & 0\\
0 & 2\tau^2/n
\end{pmatrix}.
\]

The covariance matrix of the vector MVUE $(\overline X,S^2)^\mathsf{T}$ is

\[
\operatorname{Cov}_{\mu,\tau}
\begin{pmatrix}
\overline X\\
S^2
\end{pmatrix}
=
\begin{pmatrix}
\tau/n & 0\\
0 & 2\tau^2/(n-1)
\end{pmatrix}.
\]

Therefore $\overline X$ attains the lower bound for the mean component, but $S^2$ does not attain the lower bound for the variance component. Indeed, the equality condition for estimating $\tau$ requires

\[
\begin{aligned}
\delta_\tau
&=\tau+\frac{2\tau^2}{n}S_\tau\\
&=\frac{1}{n}\sum_{i=1}^n(X_i-\mu)^2.
\end{aligned}
\]

The right-hand side depends on the unknown nuisance parameter $\mu$, so it is not a statistic. There is no single unbiased estimator that can satisfy this equality condition for all $(\mu,\tau)$.

If $\mu$ is known, then the above expression is indeed a legitimate estimator and attains the lower bound. This example clearly distinguishes between the MVUE and the attainability of the Cramér-Rao lower bound.

---

## 8. Nonregular Models: Hammersley-Chapman-Robbins Bound

The proof of the Cramér-Rao lower bound requires differentiation of the model. When the model is not differentiable, or when the support varies with the parameter, a finite-difference form of the lower bound can be used.

!!! success "Theorem 8.1 (Hammersley-Chapman-Robbins Bound)"

    Suppose

    \[
    P_{\theta+\Delta}\ll P_\theta,
    \]

    and define the likelihood ratio

    \[
    L_\Delta(X)
    =\frac{p_{\theta+\Delta}(X)}{p_\theta(X)}.
    \]

    If $L_\Delta$ has finite second moment under $P_\theta$, and $\delta$ is an unbiased estimator of $g(\theta)$, then

    \[
    \operatorname{Var}_\theta(\delta)
    \geq
    \frac{\{g(\theta+\Delta)-g(\theta)\}^2}
    {E_\theta\{L_\Delta(X)-1\}^2}.
    \]

    Taking the supremum over all admissible $\Delta$ yields the strongest bound of this form.

??? proof "Proof of Theorem 8.1 (click to expand)"

    By the change-of-measure formula,

    \[
    E_\theta L_\Delta=1,
    \]

    and

    \[
    E_\theta(\delta L_\Delta)
    =E_{\theta+\Delta}\delta
    =g(\theta+\Delta).
    \]

    Let

    \[
    Z=L_\Delta-1.
    \]

    Then

    \[
    E_\theta Z=0,
    \]

    and

    \[
    \operatorname{Cov}_\theta(\delta,Z)
    =g(\theta+\Delta)-g(\theta).
    \]

    Applying the Cauchy-Schwarz inequality,

    \[
    \operatorname{Cov}_\theta(\delta,Z)^2
    \leq
    \operatorname{Var}_\theta(\delta)
    \operatorname{Var}_\theta(Z),
    \]

    which gives the desired lower bound. $\square$

Under the condition of differentiability in quadratic mean, dividing both the numerator and denominator by $\Delta^2$ and letting $\Delta\to0$ recovers the scalar Cramér-Rao lower bound.

!!! warning "Models Whose Support Varies with the Parameter"

    For $\operatorname{Unif}(0,\theta)$, the support varies with $\theta$, so conditions (R1)-(R2) fail, and the usual score computation should not be used directly. In this case, the maximum order statistic and finite-difference-type information lower bounds are more appropriate.

---

## Appendix A: Complete Order Statistics in the Completely Nonparametric Model

Let $\mu$ be a nonatomic $\sigma$-finite measure, and let $\mathcal{F}_\mu$ be the family of all probability densities with respect to $\mu$. If

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}f,
\]

denote the vector of order statistics by

\[
O=(X_{(1)},\ldots,X_{(n)}).
\]

!!! success "Theorem A.1 (Complete Sufficiency of Order Statistics)"

    In the completely nonparametric family

    \[
    \left\{f^n:f\in\mathcal{F}_\mu\right\}
    \]

    the vector of order statistics $O$ is a complete sufficient statistic.

??? proof "Proof of Theorem A.1 (click to expand)"

    The joint density

    \[
    \prod_{i=1}^n f(x_i)
    \]

    is a symmetric function of the sample and hence can be written as a function of the order statistics. By the factorization theorem, $O$ is sufficient.

    To prove completeness, suppose $a(O)$ is integrable and that for every $f\in\mathcal{F}_\mu$,

    \[
    E_fa(O)=0.
    \]

    Define the symmetric function

    \[
    H(x_1,\ldots,x_n)
    =a\{x_{(1)},\ldots,x_{(n)}\}.
    \]

    Then for every probability density $f$,

    \[
    \int H(x_1,\ldots,x_n)
    \prod_{i=1}^nf(x_i)
    \,d\mu(x_1)\cdots d\mu(x_n)=0.
    \]

    Take an arbitrary base density $f_0$ and densities $f_1,\ldots,f_n$. For $\alpha_j>0$, define the normalized positive mixture

    \[
    f_\alpha
    =\frac{f_0+\sum_{j=1}^n\alpha_jf_j}
    {1+\sum_{j=1}^n\alpha_j}.
    \]

    Substituting $f_\alpha$ into the previous expression and multiplying by

    \[
    \left(1+\sum_{j=1}^n\alpha_j\right)^n,
    \]

    yields a polynomial in the free positive variables $\alpha_1,\ldots,\alpha_n$. This polynomial is identically zero on an open subset of $\mathbb{R}^n$, so every coefficient is zero.

    In particular, the coefficient of $\alpha_1\cdots\alpha_n$ is the sum of the integrals corresponding to all permutations. Since $H$ is symmetric, these integrals are all equal, so

    \[
    \int H(x_1,\ldots,x_n)
    \prod_{i=1}^nf_i(x_i)
    \,d\mu^n(x)=0
    \]

    for arbitrary densities $f_1,\ldots,f_n$.

    Let each $f_i$ be the normalized indicator function of an arbitrary finite-measure set $A_i$. Then $H$ has zero integral over every finite-measure rectangle

    \[
    A_1\times\cdots\times A_n.
    \]

    By the monotone class theorem, the corresponding signed measure is zero on every finite-measure measurable set; by $\sigma$-finiteness,

    \[
    H=0,
    \qquad \mu^n\text{-a.e.}
    \]

    Hence for every $f$, $a(O)=0$ almost surely, and completeness is proved. $\square$

!!! note "Why Is Nonatomicity Needed?"

    Nonatomicity guarantees that ties among sample values occur with probability zero, so that conditional on the order statistics, the conditional distribution of the label permutations is uniform. This is precisely the fact used in the conditional-expectation representation of U-statistics.

    If the distribution contains atoms, the counting of permutation orbits must be corrected according to repeated observed values, but the symmetrization idea still holds.

---

## Appendix B: Degree and Nonexistence in Nonparametric Unbiased Estimation

If $m$ is the smallest positive integer satisfying

\[
q(F)=E_Fh(X_1,\ldots,X_m),
\]

then the **degree** of the unbiasedly estimable functional $q(F)$ is $m$.

!!! success "Proposition B.1 (Necessary Mixture Polynomial Condition)"

    If the degree of $q$ is at most $m$, then for any two distributions $F,G$ in the convex model, the map

    \[
    \alpha
    \longmapsto
    q\{\alpha F+(1-\alpha)G\}
    \]

    must be a polynomial in $\alpha$ of degree at most $m$.

??? proof "Proof of Proposition B.1 (click to expand)"

    Let

    \[
    H_\alpha=\alpha F+(1-\alpha)G.
    \]

    Then

    \[
    q(H_\alpha)
    =\int h(x_1,\ldots,x_m)
    \prod_{j=1}^m
    \{\alpha\,dF(x_j)+(1-\alpha)\,dG(x_j)\}.
    \]

    Expanding the product yields a polynomial of degree at most $m$. $\square$

For the variance functional,

\[
\begin{aligned}
\operatorname{Var}_{\alpha F+(1-\alpha)G}(X)
&=\alpha E_FX^2+(1-\alpha)E_GX^2\\
&\quad
-\{\alpha E_FX+(1-\alpha)E_GX\}^2.
\end{aligned}
\]

When $F$ and $G$ have different means, this is a genuine quadratic polynomial in $\alpha$. Hence the variance cannot have degree $1$; the second-order kernel in Example 4.6 shows that the degree of the variance is exactly $2$.

The standard deviation generally does not have finite degree, because the square root of the above quadratic polynomial is generally not a polynomial. This shows that U-statistics can construct unbiased estimators only for functionals that have a finite-order expectation representation.

---

## 9. Checklist for Verifying Information Lower Bounds

Before claiming that an estimator attains the Cramér-Rao lower bound, one should check each of the following:

1. Is the estimator unbiased over the entire parameter space?
2. Is the support independent of the parameter? Are the required interchanges of differentiation and integration valid?
3. Has the correct score function for the entire sample been used?
4. In multiparameter problems, has the full Fisher information matrix including nuisance parameters been used?
5. Is the current result merely a lower bound? Does the estimator actually satisfy the scalar or vector equality condition?
6. Does the equality formula define a single statistic that does not depend on unknown parameters?
7. Are the three concepts “efficient”, “MVUE”, and “admissible” clearly distinguished?

---

## 10. Summary of This Chapter

- Under squared loss, the risk of an unbiased estimator is just its variance, so MRU and MVUE are equivalent.
- The core characterization of the MVUE is: it is orthogonal in the covariance sense to every globally zero-mean unbiased perturbation.
- A vector MVUE is equivalent to every linear combination being a scalar MVUE, and also to the covariance matrix being minimal in the Loewner order.
- U-statistics systematically construct unbiased estimators by averaging a kernel over all subsamples; in completely nonparametric models, they are also the unique MVUE of the corresponding functional.
- The Cramér-Rao inequality is derived from the score identity and the Cauchy-Schwarz inequality, but attaining the lower bound also requires a strict equality condition.
- Global efficiency is closely related to exponential family structure; even an MVUE need not attain the information lower bound.
- In multiparameter models, the full information matrix must be used, and nuisance parameters cannot be ignored.
- When the regularity conditions fail, the Hammersley-Chapman-Robbins bound provides a finite-difference alternative lower bound.