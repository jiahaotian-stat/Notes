# Chapter 6: Bayes Rules, Generalized Bayes, and Empirical Bayes

This chapter moves from frequentist risk to Bayesian decision theory. In the previous chapter, admissibility compared the entire risk function
\[
\theta\longmapsto R(\theta,\delta),
\]
but risk curves of two estimators often cross, so they cannot be uniformly ordered. Bayesian decision theory introduces a prior distribution on the parameter space, averages the risk function into a scalar Bayes risk, and turns the global optimization problem into pointwise minimization of posterior expected loss.

We will study Bayes risk and posterior risk, Bayes estimators under common loss functions, conjugate priors, generalized Bayes procedures, Jeffreys priors, MAP and regularization, and the basic ideas of empirical Bayes.

---

## 1. From frequentist risk to Bayes risk

### 1.1 Rule, estimator, and estimate

Let the sample space be $\mathcal X$ and the action space be $\mathcal A$. A **decision rule** is a measurable map

\[
\delta:\mathcal X\to\mathcal A.
\]

When the action $a$ is intended to estimate some target $g(\theta)$, the same object is also called an **estimator**.

After observing a particular data value $X=x$,

\[
\delta(x)
\]

is the realized number or vector and is called an **estimate**.

!!! note "Terminology"

    - **Bayes rule**: a Bayes-optimal decision rule in a general decision problem;
    - **Bayes estimator**: a Bayes rule when the action is used to estimate a parameter or a function $g(\theta)$;
    - **Bayes estimate**: the realized value $\delta_\pi(x)$ after observing $x$.

### 1.2 Frequentist risk

Given a loss function $L(\theta,a)$, the frequentist risk is

\[
R(\theta,\delta)
=
E_\theta L\{\theta,\delta(X)\},
\]

where $\theta$ is treated as a fixed unknown quantity and the expectation is taken only over repeated sampling

\[
X\sim P_\theta.
\]

Thus frequentist theory compares the entire function

\[
\theta\mapsto R(\theta,\delta).
\]

### 1.3 Introducing a prior distribution

Bayesian decision theory introduces a prior distribution $\pi$ on the parameter space $\Theta$:

\[
\theta\sim\pi,
\qquad
X\mid\theta\sim P_\theta.
\]

The joint distribution of $(\theta,X)$ can therefore be written as

\[
H(d\theta,dx)
=
\pi(d\theta)P_\theta(dx).
\]

The marginal distribution of $X$ is

\[
m_\pi(dx)
=
\int_\Theta P_\theta(dx)\,\pi(d\theta).
\]

Whenever a regular conditional distribution exists, the posterior distribution is denoted by

\[
\pi(d\theta\mid x).
\]

The same joint law admits the factorization

\[
\pi(d\theta)P_\theta(dx)
=
m_\pi(dx)\pi(d\theta\mid x).
\]

If $P_\theta$ has density $p(x\mid\theta)$ with respect to $\mu$ and the prior has density $\pi(\theta)$ with respect to $\nu$, then

\[
m_\pi(x)
=
\int_\Theta
p(x\mid\theta)\pi(\theta)\,\nu(d\theta),
\]

and

\[
\pi(\theta\mid x)
=
\frac{p(x\mid\theta)\pi(\theta)}
{m_\pi(x)}.
\]

Since the normalizing constant $m_\pi(x)$ does not depend on the action $a$, a Bayes rule may equivalently be found by minimizing

\[
\int_\Theta
L(\theta,a)p(x\mid\theta)\pi(\theta)\,\nu(d\theta).
\]

---

## 2. Bayes risk and posterior risk

!!! info "Definition 2.1 (Bayes risk)"

    For a prior $\pi$, the **Bayes risk**, or **integrated risk**, of a decision rule $\delta$ is

    \[
    r(\pi,\delta)
    =
    \int_\Theta
    R(\theta,\delta)\,\pi(d\theta).
    \]

    Equivalently,

    \[
    r(\pi,\delta)
    =
    \int_\Theta
    \int_{\mathcal X}
    L\{\theta,\delta(x)\}
    P_\theta(dx)\pi(d\theta).
    \]

    A rule $\delta_\pi$ is called a **Bayes rule** with respect to $\pi$ if it minimizes $r(\pi,\delta)$ over all decision rules.

Define the **posterior expected loss**

\[
C(a,x)
=
\int_\Theta
L(\theta,a)\pi(d\theta\mid x).
\]

Using the two factorizations of the joint distribution,

\[
\begin{aligned}
r(\pi,\delta)
&=
\int_{\mathcal X}
\int_\Theta
L\{\theta,\delta(x)\}
\pi(d\theta\mid x)
m_\pi(dx)\\
&=
\int_{\mathcal X}
C\{\delta(x),x\}\,m_\pi(dx).
\end{aligned}
\]

Thus Bayes risk minimization can be performed separately for each value of $x$.

!!! success "Theorem 2.2 (Posterior-risk characterization)"

    Assume that at least one decision rule has finite Bayes risk. If, for $m_\pi$-almost every $x$,

    \[
    C\{\delta_\pi(x),x\}
    =
    \inf_{a\in\mathcal A}C(a,x),
    \]

    then $\delta_\pi$ is a Bayes rule.

??? proof "Proof idea for Theorem 2.2 (click to expand)"

    Since

    \[
    r(\pi,\delta)
    =
    \int_{\mathcal X}
    C\{\delta(x),x\}\,m_\pi(dx),
    \]

    choosing, for each observed value $x$, an action that minimizes $C(a,x)$ also minimizes the entire integral.

    Therefore,

    \[
    \boxed{
    \text{minimize Bayes risk}
    \iff
    \text{minimize posterior expected loss pointwise}
    }.
    \]

---

## 3. Bayes estimators under common loss functions

Let the target be $g(\theta)$ and let the action satisfy $a\in\mathbb R$.

### 3.1 Squared-error loss

Suppose

\[
L(\theta,a)
=
\{g(\theta)-a\}^2.
\]

The posterior risk is

\[
C(a,x)
=
E\left[
\{g(\theta)-a\}^2
\mid X=x
\right].
\]

Expanding,

\[
C(a,x)
=
E\{g(\theta)^2\mid x\}
-2aE\{g(\theta)\mid x\}
+a^2.
\]

Differentiating with respect to $a$ and setting the derivative equal to zero gives

\[
a
=
E\{g(\theta)\mid X=x\}.
\]

!!! success "Conclusion: squared-error loss"

    Under squared-error loss, the Bayes estimator is the **posterior mean**:

    \[
    \boxed{
    \delta_\pi(x)
    =
    E\{g(\theta)\mid X=x\}
    }.
    \]

The minimum posterior risk is

\[
\inf_a C(a,x)
=
\operatorname{Var}\{g(\theta)\mid X=x\}.
\]

### 3.2 Weighted squared-error loss

Consider

\[
L(\theta,a)
=
w(\theta)\{g(\theta)-a\}^2,
\qquad
w(\theta)\ge 0.
\]

Define

\[
A_j(x)
=
E\{w(\theta)g(\theta)^j\mid X=x\},
\qquad
j=0,1,2.
\]

Then

\[
C(a,x)
=
A_2(x)-2aA_1(x)+a^2A_0(x).
\]

When

\[
0<A_0(x)<\infty
\]

and the required moments are finite,

\[
\boxed{
\delta_\pi(x)
=
\frac{
E\{w(\theta)g(\theta)\mid X=x\}
}{
E\{w(\theta)\mid X=x\}
}
}.
\]

Moreover,

\[
C(a,x)
=
C\{\delta_\pi(x),x\}
+
A_0(x)\{a-\delta_\pi(x)\}^2.
\]

!!! example "Example 3.1 (Relative squared error for a positive parameter)"

    Suppose the posterior, in rate parameterization, is

    \[
    \theta\mid x
    \sim
    \operatorname{Gamma}(k,r),
    \]

    and the loss is

    \[
    L(\theta,a)
    =
    \left(\frac a\theta-1\right)^2
    =
    \theta^{-2}(a-\theta)^2.
    \]

    Hence

    \[
    g(\theta)=\theta,
    \qquad
    w(\theta)=\theta^{-2}.
    \]

    For $k>2$,

    \[
    \delta_\pi(x)
    =
    \frac{E(\theta^{-1}\mid x)}
    {E(\theta^{-2}\mid x)}
    =
    \frac{k-2}{r}.
    \]

    Under ordinary squared-error loss, however, the Bayes estimator is

    \[
    E(\theta\mid x)
    =
    \frac{k}{r}.
    \]

!!! note "Important idea"

    **The posterior distribution itself has not changed.**

    What changes is the loss function, and therefore the optimal action changes as well. A Bayesian decision is determined not only by the posterior distribution, but also by how the cost of error is defined.

### 3.3 Absolute-error loss

If

\[
L(\theta,a)
=
|g(\theta)-a|,
\]

then any posterior median is a Bayes estimator; that is, $a$ satisfies

\[
P\{g(\theta)\le a\mid X=x\}
\ge \frac12,
\]

and

\[
P\{g(\theta)\ge a\mid X=x\}
\ge \frac12.
\]

Hence

\[
\boxed{
\text{absolute loss}
\longrightarrow
\text{posterior median}
}.
\]

### 3.4 0–1 loss and MAP

For a discrete parameter, consider exact 0–1 loss

\[
L(\theta,a)
=
\mathbf 1\{a\ne\theta\}.
\]

The posterior risk is

\[
C(a,x)
=
P(\theta\ne a\mid X=x)
=
1-P(\theta=a\mid X=x).
\]

Thus the Bayes rule is the posterior mode:

\[
\boxed{
\delta_\pi(x)
\in
\arg\max_a
P(\theta=a\mid X=x)
}.
\]

For a continuous parameter,

\[
P(\theta=a\mid X=x)=0
\]

typically holds, so exact 0–1 loss is not useful. A common alternative is interval 0–1 loss,

\[
L(\theta,a)
=
\mathbf 1\{|g(\theta)-a|>c\}.
\]

The optimal action is the center of a length-$2c$ interval having the largest posterior probability:

\[
\delta_\pi(x)
\in
\arg\max_a
P\{a-c\le g(\theta)\le a+c\mid X=x\}.
\]

As $c\downarrow0$, this informally leads to the posterior mode, or MAP estimator, when a posterior density exists.

!!! note "Common losses and Bayes point estimators"

    \[
    \begin{array}{c|c}
    \text{Loss} & \text{Bayes point estimator}\\
    \hline
    \text{Squared error} & \text{Posterior mean}\\
    \text{Absolute error} & \text{Posterior median}\\
    \text{0--1 / local 0--1} & \text{Posterior mode / MAP}
    \end{array}
    \]

---

## 4. Conjugate priors

!!! info "Definition 4.1 (Conjugate family)"

    For a sampling model $p(x\mid\theta)$, a family of priors $\mathcal F$ is called **conjugate** if

    \[
    \pi\in\mathcal F
    \Longrightarrow
    \pi(\cdot\mid x)\in\mathcal F
    \]

    for every possible observation $x$.

The key advantage of conjugacy is **algebraic closure**: the prior and posterior have the same functional form, and only the hyperparameters are updated.

---

## 5. Beta–Binomial model

Let

\[
X_1,\ldots,X_n
\overset{\text{iid}}{\sim}
\operatorname{Bernoulli}(p),
\]

and define

\[
T=\sum_{i=1}^nX_i.
\]

The likelihood satisfies

\[
L(p;x)
\propto
p^T(1-p)^{n-T}.
\]

Take a Beta prior

\[
p\sim\operatorname{Beta}(a,b),
\]

with density

\[
\pi(p)
=
\frac{\Gamma(a+b)}
{\Gamma(a)\Gamma(b)}
p^{a-1}(1-p)^{b-1},
\qquad
0<p<1.
\]

Then

\[
\pi(p\mid x)
\propto
p^{a+T-1}
(1-p)^{b+n-T-1},
\]

so

\[
\boxed{
p\mid X
\sim
\operatorname{Beta}(a+T,b+n-T)
}.
\]

Under squared-error loss,

\[
\delta_\pi(X)
=
E(p\mid X)
=
\frac{a+T}{a+b+n}.
\]

This can be rewritten in shrinkage form:

\[
\frac{a+T}{a+b+n}
=
\frac{n}{n+a+b}\bar X
+
\frac{a+b}{n+a+b}
\frac{a}{a+b}.
\]

Thus the posterior mean shrinks the sample proportion $\bar X$ toward the prior mean

\[
\frac{a}{a+b},
\]

while

\[
a+b
\]

acts as a prior strength.

!!! example "Example 5.1 (Uniform prior and Laplace smoothing)"

    If

    \[
    a=b=1,
    \]

    then

    \[
    \delta_\pi(X)
    =
    \frac{T+1}{n+2}.
    \]

    This is Laplace smoothing. Even when $T=0$ or $T=n$, the estimate is not exactly $0$ or $1$.

---

## 6. Gamma–Poisson model (optional enrichment)

Let

\[
X_1,\ldots,X_n
\overset{\text{iid}}{\sim}
\operatorname{Poisson}(\lambda),
\qquad
T=\sum_{i=1}^nX_i.
\]

The likelihood satisfies

\[
L(\lambda;x)
\propto
e^{-n\lambda}\lambda^T.
\]

Use the Gamma prior in rate parameterization:

\[
\lambda
\sim
\operatorname{Gamma}(\alpha,\beta),
\]

\[
\pi(\lambda)
=
\frac{\beta^\alpha}{\Gamma(\alpha)}
\lambda^{\alpha-1}e^{-\beta\lambda},
\qquad
\lambda>0.
\]

Then

\[
\pi(\lambda\mid x)
\propto
\lambda^{\alpha+T-1}
e^{-(\beta+n)\lambda},
\]

and therefore

\[
\boxed{
\lambda\mid X
\sim
\operatorname{Gamma}(\alpha+T,\beta+n)
}.
\]

Under squared-error loss,

\[
\delta_\pi(X)
=
E(\lambda\mid X)
=
\frac{\alpha+T}{\beta+n}.
\]

It may also be written as

\[
\frac{\alpha+T}{\beta+n}
=
\frac{n}{\beta+n}
\frac{T}{n}
+
\frac{\beta}{\beta+n}
\frac{\alpha}{\beta}.
\]

Hence it is a weighted average of the sample mean and the prior mean, with $\beta$ behaving like a prior sample size.

---

## 7. Gamma–Exponential model: comparing three estimation principles

Let $n\ge2$ and

\[
X_1,\ldots,X_n
\overset{\text{iid}}{\sim}
f(x\mid\theta)
=
\theta e^{-\theta x}\mathbf1\{x>0\},
\qquad
\theta>0.
\]

Define

\[
T=\sum_{i=1}^nX_i.
\]

The likelihood is

\[
L(\theta;x)
\propto
\theta^n e^{-\theta T}.
\]

Take the Gamma prior

\[
\theta
\sim
\operatorname{Gamma}(\beta,\alpha),
\]

whose density is proportional to

\[
\theta^{\beta-1}e^{-\alpha\theta}.
\]

Then

\[
\boxed{
\theta\mid X
\sim
\operatorname{Gamma}(n+\beta,\alpha+T)
}.
\]

We wish to estimate the reliability probability that a new component survives beyond $t_0$:

\[
g(\theta)
=
P_\theta(X_{\mathrm{new}}>t_0)
=
e^{-\theta t_0}.
\]

### 7.1 Plug-in estimator

The maximum likelihood estimator is

\[
\widehat\theta_{\mathrm{MLE}}
=
\frac nT.
\]

By invariance of the MLE,

\[
\boxed{
\widehat g_{\mathrm{plug}}
=
e^{-nt_0/T}
}.
\]

### 7.2 Minimum-risk unbiased estimator

Because

\[
\mathbf1\{X_1>t_0\}
\]

satisfies

\[
E_\theta\mathbf1\{X_1>t_0\}
=
e^{-\theta t_0}
=
g(\theta),
\]

it is an unbiased estimator of $g(\theta)$.

Conditional on the complete sufficient statistic $T$,

\[
\left(
\frac{X_1}{T},
\ldots,
\frac{X_n}{T}
\right)
\]

is uniform on the simplex. Rao–Blackwellization therefore gives

\[
\boxed{
\widehat g_{\mathrm{MRU}}
=
E\{\mathbf1(X_1>t_0)\mid T\}
=
\left(1-\frac{t_0}{T}\right)_+^{n-1}
}.
\]

Under convex loss, this is the minimum-risk unbiased rule; in particular, under squared-error loss it is the UMVU estimator.

### 7.3 Bayes estimator

Under squared-error loss, the Bayes estimator is

\[
E(e^{-\theta t_0}\mid X).
\]

Using the Laplace transform of the Gamma posterior,

\[
\boxed{
\delta_\pi(X)
=
\left(
\frac{\alpha+T}
{\alpha+T+t_0}
\right)^{n+\beta}
=
\left(
1+\frac{t_0}{\alpha+T}
\right)^{-(n+\beta)}
}.
\]

!!! note "The three estimators solve different optimization problems"

    - Plug-in estimator: substitute $\widehat\theta_{\mathrm{MLE}}$ into the target function;
    - MRU / UMVU estimator: optimize within the class of unbiased estimators;
    - Bayes estimator: minimize Bayes risk under a specified prior and loss.

    There is no reason for them to coincide.

---

## 8. Normal–Normal model

Let

\[
X_1,\ldots,X_n
\overset{\text{iid}}{\sim}
N(\theta,\sigma^2),
\]

where $\sigma^2$ is known, and take the prior

\[
\theta
\sim
N(\mu,\tau^2).
\]

The likelihood depends on the data only through $\bar X$. Completing the square gives

\[
\theta\mid X
\sim
N(m_n,v_n),
\]

where

\[
v_n
=
\frac{1}
{n/\sigma^2+1/\tau^2},
\]

and

\[
m_n
=
v_n
\left(
\frac{n\bar X}{\sigma^2}
+
\frac{\mu}{\tau^2}
\right).
\]

Equivalently,

\[
m_n
=
\frac{n\sigma^{-2}}
{n\sigma^{-2}+\tau^{-2}}
\bar X
+
\frac{\tau^{-2}}
{n\sigma^{-2}+\tau^{-2}}
\mu.
\]

Under squared-error loss,

\[
\boxed{
\delta_\pi(X)=m_n
}.
\]

The posterior risk of this Bayes rule is

\[
\operatorname{Var}(\theta\mid X)
=
v_n.
\]

In this conjugate normal model, $v_n$ does not depend on the realized data $X$.

!!! note "Shrinkage interpretation"

    The posterior mean is a **precision-weighted average** of the sample mean and the prior mean.

    The data precision is

    \[
    \frac{n}{\sigma^2},
    \]

    while the prior precision is

    \[
    \frac1{\tau^2}.
    \]

    As $n$ increases, the data receive more weight; as $\tau^2$ becomes small, the prior becomes more concentrated and receives more weight.

---

## 9. Conjugate priors in exponential families (optional enrichment)

Consider the canonical exponential family

\[
p_\eta(t)
=
\exp\{\eta^Tt-A(\eta)\}K(t),
\qquad
\eta\in\mathcal H.
\]

For i.i.d. data $T_1,\ldots,T_n$, the likelihood is proportional to

\[
\exp\left\{
\eta^T\sum_{i=1}^nT_i
-
nA(\eta)
\right\}.
\]

A conjugate prior may be written as

\[
\pi_{n_0,t_0}(\eta)
\propto
\exp\{
n_0t_0^T\eta-n_0A(\eta)
\}.
\]

The posterior is

\[
\pi(\eta\mid T_1,\ldots,T_n)
\propto
\exp\left\{
\left(
n_0t_0+\sum_{i=1}^nT_i
\right)^T\eta
-
(n_0+n)A(\eta)
\right\}.
\]

Hence the hyperparameters update according to

\[
n_0
\longmapsto
n_0+n,
\]

and

\[
t_0
\longmapsto
\frac{
n_0t_0+\sum_{i=1}^nT_i
}{
n_0+n
}.
\]

!!! note "Pseudo-observation interpretation"

    The quantity $n_0$ behaves like an effective prior sample size, while $t_0$ behaves like a prior mean sufficient statistic.

---

## 10. Improper priors and generalized Bayes

### 10.1 Improper prior

If

\[
\int_\Theta\pi(d\theta)=\infty,
\]

then $\pi$ is called an **improper prior**.

It is not a probability distribution, but it may still yield a proper posterior. If

\[
m(x)
=
\int_\Theta
p(x\mid\theta)\pi(d\theta)
<
\infty,
\]

then the formal posterior

\[
\pi(d\theta\mid x)
=
\frac{
p(x\mid\theta)\pi(d\theta)
}{
m(x)
}
\]

is a genuine probability distribution.

!!! info "Definition 10.1 (Generalized Bayes rule)"

    A decision rule is called a **generalized Bayes rule** with respect to an improper prior if it minimizes posterior expected loss under the formal posterior generated by that prior.

### 10.2 Flat prior for the normal mean

Let

\[
X_1,\ldots,X_n
\overset{\text{iid}}{\sim}
N(\theta,\sigma^2),
\]

with $\sigma^2$ known.

Use the improper flat prior

\[
\pi(\theta)\propto1.
\]

Then

\[
\theta\mid X
\sim
N\left(
\bar X,
\frac{\sigma^2}{n}
\right).
\]

Therefore, under squared-error loss,

\[
\boxed{
\delta(X)=\bar X
}
\]

is a generalized Bayes estimator.

### 10.3 Why is $\bar X$ not proper Bayes?

!!! success "Proposition 10.2 ($\bar X$ is not proper Bayes)"

    In the unrestricted normal location model, $\bar X$ cannot be the posterior mean under any proper prior on $\mathbb R$, although it is generalized Bayes under Lebesgue measure.

??? proof "Proof of Proposition 10.2 (click to expand)"

    Let

    \[
    s^2=\frac{\sigma^2}{n},
    \]

    and suppose $\pi$ is a proper prior. The marginal density of $\bar X$ is

    \[
    m_\pi(x)
    =
    \int
    \phi_s(x-\theta)\,\pi(d\theta),
    \]

    where $\phi_s$ is the normal density with variance $s^2$.

    Differentiating with respect to $x$ gives Tweedie's identity:

    \[
    E_\pi(\theta\mid\bar X=x)
    =
    x
    +
    s^2
    \frac{m_\pi'(x)}{m_\pi(x)}.
    \]

    If

    \[
    E_\pi(\theta\mid\bar X=x)=x
    \]

    for every $x$, then necessarily

    \[
    m_\pi'(x)=0
    \]

    for every $x$.

    Thus $m_\pi(x)$ must be a positive constant on $\mathbb R$. But a positive constant cannot integrate to $1$ over $\mathbb R$.

    Therefore no proper prior can have posterior mean identically equal to $x$. $\square$

!!! warning "Proper Bayes and generalized Bayes are different"

    The estimator

    \[
    \bar X
    \]

    is generalized Bayes because a flat measure produces a proper posterior;

    but it is not proper Bayes because no genuine probability prior produces it as the posterior mean.

    Moreover, an improper prior generally does not define a finite prior-averaged Bayes risk, so being generalized Bayes does not mean minimizing a finite Bayes risk under some proper prior.

---

## 11. Jeffreys prior and reference priors

A so-called “non-informative prior” should not be interpreted as literally containing no information. A density that is flat in one parameterization will generally not remain flat after reparameterization.

A common invariant default in one-dimensional problems is Jeffreys prior.

!!! info "Definition 11.1 (Jeffreys prior)"

    If $I(\theta)$ is the Fisher information, define

    \[
    \boxed{
    \pi_J(\theta)
    \propto
    \sqrt{I(\theta)}
    }.
    \]

Its main motivation is reparameterization invariance. If

\[
\phi=h(\theta)
\]

is a smooth one-to-one transformation, then Jeffreys prior satisfies

\[
\pi_J(\theta)d\theta
=
\pi_J(\phi)d\phi.
\]

### 11.1 Bernoulli parameter

For $n$ Bernoulli observations,

\[
I_n(p)
=
\frac{n}{p(1-p)}.
\]

Therefore

\[
\pi_J(p)
\propto
\{p(1-p)\}^{-1/2}.
\]

This is exactly the

\[
\operatorname{Beta}\left(\frac12,\frac12\right)
\]

prior.

If

\[
T=\sum_{i=1}^nX_i,
\]

then the posterior is

\[
p\mid X
\sim
\operatorname{Beta}
\left(
T+\frac12,
n-T+\frac12
\right).
\]

### 11.2 Location family

For

\[
p(x\mid\theta)
=
f(x-\theta),
\]

the Fisher information is constant in $\theta$, so

\[
\boxed{
\pi_J(\theta)\propto1
}.
\]

### 11.3 Scale family

For

\[
p(x\mid\sigma)
=
\frac1\sigma
f\left(\frac x\sigma\right),
\]

the Fisher information is proportional to $1/\sigma^2$, so

\[
\boxed{
\pi_J(\sigma)\propto\frac1\sigma
}.
\]

This prior is usually improper but is invariant under scale transformations.

### 11.4 Multivariate parameter

If

\[
\theta\in\mathbb R^d,
\]

then Jeffreys prior takes the form

\[
\pi_J(\theta)
\propto
\sqrt{\det I(\theta)},
\]

where $I(\theta)$ is the Fisher information matrix.

### 11.5 Reference priors

The motivation of reference priors is related to Jeffreys priors: one seeks a prior that asymptotically maximizes the amount of information gained from prior to posterior, commonly measured through Kullback–Leibler divergence.

In regular one-parameter problems, the reference prior often agrees with Jeffreys prior.

In multiparameter problems, especially with nuisance parameters, the resulting reference prior may depend on

- which component is the parameter of interest;
- the ordering used in the construction.

!!! warning "Required caution"

    Invariance alone does not guarantee

    - posterior propriety;
    - good finite-sample risk;
    - consistency in models with many nuisance parameters.

    In applications one should still check posterior propriety, prior sensitivity, and robustness.

---

## 12. MAP and regularization (optional enrichment)

Consider the Gaussian linear model

\[
y=A\beta+\sigma Z,
\qquad
Z\sim N(0,I_n).
\]

Suppose the prior density has the form

\[
\pi(\beta)
\propto
\exp\{-\lambda P(\beta)\}.
\]

The log posterior satisfies

\[
\log\pi(\beta\mid y)
=
-\frac1{2\sigma^2}
\|y-A\beta\|^2
-\lambda P(\beta)
+\text{constant}.
\]

Hence the MAP estimator is equivalent to minimizing

\[
\boxed{
\frac1{2\sigma^2}
\|y-A\beta\|^2
+
\lambda P(\beta)
}.
\]

Therefore:

- if
  \[
  P(\beta)=\frac12\|\beta\|_2^2,
  \]
  we obtain ridge regression;

- if
  \[
  P(\beta)=\|\beta\|_1,
  \]
  we obtain lasso-type soft thresholding;

- an $\ell_0$ penalty yields hard thresholding but is computationally harder.

### 12.1 Spike-and-slab prior

In a Gaussian sequence model,

\[
X_i\mid\mu_i
\sim
N(\mu_i,1),
\]

one may use

\[
\mu_i
\sim
(1-w)\delta_0+w\gamma,
\]

where $\delta_0$ is a point mass at zero and $\gamma$ is a continuous distribution.

Such a spike-and-slab prior explicitly encodes sparsity and can produce shrinkage or thresholding through the posterior median or MAP.

---

## 13. Empirical Bayes: learning a prior from an ensemble

In ordinary Bayesian analysis, the prior $\nu$ is specified before using the data for the current decision.

In empirical Bayes (EB), an ensemble of related observations is used to estimate

- the prior distribution $\nu$;
- its hyperparameters;
- or the marginal distribution required by the Bayes rule.

The estimated object is then plugged into the oracle Bayes rule.

Suppose

\[
(\theta_i,X_i),
\qquad
i=1,\ldots,n,
\]

are i.i.d. from the hierarchical model

\[
\theta_i\sim\nu,
\qquad
X_i\mid\theta_i\sim P_{\theta_i}.
\]

### 13.1 Historical-data EB

Use historical observations

\[
X_1,\ldots,X_m
\]

to estimate $\nu$, and then apply the fitted Bayes rule to a new observation $X_{m+1}$.

Conditional on the historical sample, this resembles an ordinary Bayes analysis using a fitted prior.

### 13.2 Compound / parallel EB

Observe many current units

\[
X_1,\ldots,X_m,
\]

estimate the common mixing distribution from the entire ensemble, and then estimate all of

\[
\theta_1,\ldots,\theta_m
\]

using that fitted distribution.

This is the large-scale setting behind many modern shrinkage and multiple-testing methods.

### 13.3 Parametric and nonparametric EB

Empirical Bayes may be

- **parametric EB**: for example, assume a normal prior and estimate only its mean and variance;
- **nonparametric EB**: estimate an unrestricted mixing distribution.

An important idea is that estimating the entire prior $\nu$ may be unnecessary. If the oracle Bayes rule depends only on the marginal distribution $f_\nu$, one may estimate $f_\nu$ directly.

### 13.4 Poisson empirical Bayes

Suppose

\[
X\mid\theta
\sim
\operatorname{Poisson}(\theta),
\qquad
\theta\sim\nu.
\]

The marginal mass function is

\[
f_\nu(x)
=
\int_0^\infty
e^{-\theta}
\frac{\theta^x}{x!}
\nu(d\theta).
\]

The posterior mean is

\[
\delta_\nu(x)
=
E(\theta\mid X=x).
\]

Since

\[
\begin{aligned}
(x+1)f_\nu(x+1)
&=
(x+1)
\int_0^\infty
e^{-\theta}
\frac{\theta^{x+1}}{(x+1)!}
\nu(d\theta)\\
&=
\int_0^\infty
\theta
e^{-\theta}
\frac{\theta^x}{x!}
\nu(d\theta),
\end{aligned}
\]

we obtain Robbins' identity:

\[
\boxed{
\delta_\nu(x)
=
\frac{(x+1)f_\nu(x+1)}
{f_\nu(x)}
}.
\]

Thus empirical Bayes can estimate

\[
f_\nu(x)
\quad\text{and}\quad
f_\nu(x+1)
\]

directly from a large ensemble and plug them into this ratio.

!!! note "What does empirical Bayes change?"

    Because

    \[
    \widehat\nu
    \]

    is itself data-dependent,

    \[
    \delta_{\widehat\nu}
    \]

    is generally not a Bayes rule with respect to any fixed prior.

    Reusing the same data both to learn the prior and to estimate the current parameters can introduce

    - hyperparameter uncertainty;
    - overconfidence;
    - additional risk-analysis issues.

    Possible responses include sample splitting, cross-fitting, explicit risk analysis, or full hierarchical Bayes.

The benefit is equally important: EB can learn the degree of pooling or shrinkage adaptively from the observed population.

---

## 14. Why does Bayes enter naturally here?

The logic of the preceding lectures can be summarized as

\[
\boxed{
\text{UMVU compares only unbiased estimators}
}
\]

\[
\Downarrow
\]

\[
\boxed{
\text{Admissibility compares full risk functions}
}
\]

\[
\Downarrow
\]

\[
\boxed{
\text{risk curves of different procedures often cross}
}
\]

\[
\Downarrow
\]

\[
\boxed{
\text{weight the parameter space by a prior }\pi
}
\]

\[
\Downarrow
\]

\[
\boxed{
\text{minimize Bayes risk / posterior risk}
}.
\]

Admissibility removes only rules that are uniformly dominated, but it usually does not identify a unique procedure among the many remaining admissible rules.

A prior $\pi$ assigns weights to the parameter space and converts the function-valued comparison

\[
\theta\mapsto R(\theta,\delta)
\]

into the scalar criterion

\[
r(\pi,\delta)
=
\int_\Theta
R(\theta,\delta)\pi(d\theta).
\]

### 14.1 Two viewpoints on Bayes

#### Decision-theoretic viewpoint

The prior need not initially be interpreted as a subjective belief. It may instead be regarded as a weighting distribution over the parameter space.

From this viewpoint, Bayesian procedures

- turn a partial ordering of risk functions into scalar optimization;
- generate candidates for admissible and complete classes;
- provide lower bounds for minimax risk and connect to least-favorable priors;
- explain why shrinkage and regularization arise naturally.

#### Statistical-inference viewpoint

The prior represents information available before the current experiment:

\[
\theta\sim\pi,
\qquad
X\mid\theta\sim P_\theta.
\]

The likelihood supplies information from the current data, and Bayes' formula combines the two into the posterior

\[
\pi(d\theta\mid x).
\]

The posterior may then be used for

- point estimation;
- uncertainty quantification;
- prediction;
- subsequent scientific decision making.

### 14.2 Why does weighted risk become posterior risk?

By the tower property,

\[
\begin{aligned}
r(\pi,\delta)
&=
\int_\Theta
R(\theta,\delta)\pi(d\theta)\\
&=
E_{\theta,X}
L\{\theta,\delta(X)\}\\
&=
E_X
\left[
E\{L(\theta,\delta(X))\mid X\}
\right].
\end{aligned}
\]

Therefore,

\[
\boxed{
\text{weighted frequentist risk}
\Longleftrightarrow
\text{posterior expected loss}
}.
\]

This is the bridge between the decision-theoretic and Bayesian-inference viewpoints.

---

## 15. What Bayesian methods contribute to statistics

Bayesian ideas play several distinct roles in statistics:

1. **Prior information**: incorporate previous studies, expert knowledge, physical restrictions, or external data.
2. **Regularization and pooling**: priors naturally produce shrinkage and allow related parameters to share information.
3. **Prediction and uncertainty**: posterior and posterior-predictive distributions propagate parameter uncertainty.
4. **Hierarchical and latent-variable models**: random effects, mixtures, missing data, and latent states can be handled in a unified probability model.
5. **Scientific decision making**: the posterior may be combined with asymmetric medical, economic, or scientific losses.
6. **Frequentist theory**: Bayes and generalized-Bayes constructions are useful for proving admissibility, minimaxity, and complete-class results.
7. **Empirical Bayes**: an ensemble of related problems can be used to learn prior structure and adapt the amount of pooling.

!!! warning "Bayesian analysis is not automatically correct"

    Conclusions may be sensitive to

    - the prior;
    - the likelihood;
    - hierarchical modeling assumptions.

    Posterior propriety, prior sensitivity, calibration, and robustness should therefore be checked.

---

## 16. Chapter summary

The main conclusions of this chapter are:

1. Bayes risk is the frequentist risk averaged over a prior distribution.
2. A Bayes rule can be obtained by minimizing posterior expected loss pointwise.
3. Squared-error, absolute-error, and 0–1 losses lead respectively to the posterior mean, median, and mode/MAP.
4. Conjugate priors preserve functional form from prior to posterior and can be interpreted in terms of pseudo-observations.
5. The Beta–Binomial, Gamma–Poisson, and Normal–Normal models all illustrate the shrinkage structure of posterior means.
6. The Gamma–Exponential example clearly distinguishes plug-in, MRU/UMVU, and Bayes optimization principles.
7. Improper priors can generate generalized Bayes rules but must not be confused with proper Bayes procedures.
8. In the normal location model, $\bar X$ is generalized Bayes under a flat prior but is not Bayes under any proper prior.
9. Jeffreys prior is constructed from Fisher information and is motivated primarily by invariance under reparameterization.
10. MAP is closely connected to penalized optimization, and ridge, lasso, and sparse priors all admit Bayesian interpretations.
11. Empirical Bayes learns a prior or marginal distribution from an ensemble of related observations and thereby adapts the amount of shrinkage.
12. Bayesian decision theory is not a replacement for frequentist risk analysis; it is a natural extension when risk curves cannot be uniformly ordered.
