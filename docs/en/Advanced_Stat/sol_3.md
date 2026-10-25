---
tags:
  - Advanced Mathematical Statistics
  - Admissibility
  - Shrinkage Estimation
  - Bayes Estimation
  - James-Stein Estimation
  - Homework
---

# Advanced Mathematical Statistics: Assignment 3


## Problem 1: Risk of a Linear Shrinkage Estimator

**Problem.** Let

\[
X\sim N(\theta,1).
\]

Under squared-error loss, compute the risk of $aX$ for estimating $\theta$, and show that no fixed $a<1$ dominates $X$ over all $\theta\in\mathbb R$.

??? success "Solution (click to expand)"

    Since

    \[
    E_\theta(X)=\theta,
    \qquad
    \operatorname{Var}_\theta(X)=1,
    \]

    we have

    \[
    E_\theta(X^2)=1+\theta^2.
    \]

    Let $\delta_a(X)=aX$. Its squared-error risk is

    \[
    R(\theta,\delta_a)
    =E_\theta[(aX-\theta)^2].
    \]

    Expanding gives

    \[
    R(\theta,\delta_a)
    =a^2E_\theta(X^2)-2a\theta E_\theta(X)+\theta^2.
    \]

    Therefore,

    \[
    \boxed{
    R(\theta,\delta_a)
    =a^2+(a-1)^2\theta^2
    }.
    \]

    For the unbiased estimator $\delta_0(X)=X$,

    \[
    R(\theta,\delta_0)
    =E_\theta[(X-\theta)^2]
    =1.
    \]

    The difference between the two risks is

    \[
    R(\theta,\delta_a)-R(\theta,\delta_0)
    =a^2-1+(a-1)^2\theta^2.
    \]

    For every fixed $a<1$, $(a-1)^2>0$, and hence

    \[
    \lim_{|\theta|\to\infty}
    \{R(\theta,\delta_a)-R(\theta,\delta_0)\}
    =+\infty.
    \]

    Consequently, there is always some $\theta\in\mathbb R$ for which

    \[
    R(\theta,\delta_a)>R(\theta,\delta_0).
    \]

    Thus no fixed $a<1$ dominates $X$ over the entire parameter space.

---

## Problem 2: Inadmissibility by Projection onto the Parameter Space

**Problem.** Suppose that $\vartheta\in[c,d]$ and the loss is $L(|\vartheta-a|)$, where $L$ is increasing. Show that every rule taking values outside $[c,d]$ with positive probability is inadmissible.

??? success "Solution (click to expand)"

    Let $\delta(X)$ be a decision rule. Define its projection onto the parameter interval by

    \[
    g(a)=
    \begin{cases}
    c, & a<c,\\
    a, & c\leq a\leq d,\\
    d, & a>d,
    \end{cases}
    \]

    and let

    \[
    \delta^*(X)=g(\delta(X)).
    \]

    For each fixed $\vartheta\in[c,d]$, compare the losses of the actions $a$ and $g(a)$.

    If $a<c$, then $\vartheta\geq c>a$, and therefore

    \[
    |\vartheta-g(a)|
    =|\vartheta-c|
    =\vartheta-c
    <\vartheta-a
    =|\vartheta-a|.
    \]

    If $a>d$, then $\vartheta\leq d<a$, and therefore

    \[
    |\vartheta-g(a)|
    =|\vartheta-d|
    =d-\vartheta
    <a-\vartheta
    =|\vartheta-a|.
    \]

    If $c\leq a\leq d$, then $g(a)=a$, so the two distances are equal. Hence

    \[
    |\vartheta-\delta^*(X)|
    \leq|\vartheta-\delta(X)|
    \qquad\text{for every }\vartheta\in[c,d].
    \]

    Since $L$ is increasing,

    \[
    L(|\vartheta-\delta^*(X)|)
    \leq L(|\vartheta-\delta(X)|).
    \]

    Taking expectations yields

    \[
    R(\vartheta,\delta^*)
    \leq R(\vartheta,\delta)
    \qquad\text{for all }\vartheta\in[c,d].
    \]

    If $L$ is strictly increasing and, for some $\vartheta_0$,

    \[
    P_{\vartheta_0}\{\delta(X)\notin[c,d]\}>0,
    \]

    then the loss inequality is strict on an event of positive probability. Therefore,

    \[
    R(\vartheta_0,\delta^*)
    <R(\vartheta_0,\delta).
    \]

    Thus $\delta^*$ dominates $\delta$, and $\delta$ is inadmissible.

---

## Problem 3: Pitman Estimation in Location Families

**Problem.** Let $X_1,\ldots,X_n$ be a random sample from a location family with density $f(x-\theta)$. The Pitman estimator is

\[
d_P(X)
=\frac{
\displaystyle\int_{-\infty}^{\infty}
t\prod_{i=1}^nf(X_i-t)\,dt
}{
\displaystyle\int_{-\infty}^{\infty}
\prod_{i=1}^nf(X_i-t)\,dt
}.
\]

(a) Show that $d_P$ is location equivariant.

(b) Find $d_P$ when $f$ is the standard normal density.

(c) Find $d_P$ when $f$ is uniform on $(-1/2,1/2)$.

??? success "Solution (click to expand)"

    **(a) Location equivariance**

    For every $b\in\mathbb R$, write

    \[
    X+b\mathbf 1=(X_1+b,\ldots,X_n+b).
    \]

    Then

    \[
    d_P(X+b\mathbf 1)
    =\frac{
    \displaystyle\int_{-\infty}^{\infty}
    t\prod_{i=1}^nf(X_i+b-t)\,dt
    }{
    \displaystyle\int_{-\infty}^{\infty}
    \prod_{i=1}^nf(X_i+b-t)\,dt
    }.
    \]

    Make the change of variables $u=t-b$, so that $t=u+b$ and $dt=du$. The numerator becomes

    \[
    \int_{-\infty}^{\infty}
    (u+b)\prod_{i=1}^nf(X_i-u)\,du.
    \]

    Splitting this expression gives

    \[
    \int_{-\infty}^{\infty}
    u\prod_{i=1}^nf(X_i-u)\,du
    +b\int_{-\infty}^{\infty}
    \prod_{i=1}^nf(X_i-u)\,du.
    \]

    The denominator becomes

    \[
    \int_{-\infty}^{\infty}
    \prod_{i=1}^nf(X_i-u)\,du.
    \]

    Therefore, whenever the defining integrals exist,

    \[
    \boxed{d_P(X+b\mathbf 1)=d_P(X)+b},
    \]

    proving location equivariance.

    **(b) Standard normal case**

    When

    \[
    f(z)=\frac1{\sqrt{2\pi}}e^{-z^2/2}
    \]

    we have

    \[
    \prod_{i=1}^nf(X_i-t)
    =(2\pi)^{-n/2}
    \exp\left\{-\frac12\sum_{i=1}^n(X_i-t)^2\right\}.
    \]

    Completing the square gives

    \[
    \sum_{i=1}^n(X_i-t)^2
    =n(t-\overline X)^2
    +\sum_{i=1}^n(X_i-\overline X)^2.
    \]

    The factor independent of $t$ cancels from the ratio, so

    \[
    d_P(X)
    =\frac{
    \displaystyle\int_{-\infty}^{\infty}
    t\exp\{-n(t-\overline X)^2/2\}\,dt
    }{
    \displaystyle\int_{-\infty}^{\infty}
    \exp\{-n(t-\overline X)^2/2\}\,dt
    }.
    \]

    This ratio is the mean of a $N(\overline X,1/n)$ density in $t$. Hence

    \[
    \boxed{d_P(X)=\overline X}.
    \]

    **(c) Uniform case**

    When

    \[
    f(z)=\mathbf 1_{\{-1/2<z<1/2\}}
    \]

    we have

    \[
    \prod_{i=1}^nf(X_i-t)
    =\prod_{i=1}^n
    \mathbf 1_{\{X_i-1/2<t<X_i+1/2\}}.
    \]

    Therefore,

    \[
    \prod_{i=1}^nf(X_i-t)
    =\mathbf 1_{\{X_{(n)}-1/2<t<X_{(1)}+1/2\}}.
    \]

    Define

    \[
    a=X_{(n)}-\frac12,
    \qquad
    b=X_{(1)}+\frac12.
    \]

    Under the model, $a<b$ almost surely. Hence

    \[
    d_P(X)
    =\frac{\int_a^bt\,dt}{\int_a^b1\,dt}
    =\frac{a+b}{2}.
    \]

    Thus

    \[
    \boxed{
    d_P(X)=\frac{X_{(1)}+X_{(n)}}2
    },
    \]

    the midrange of the sample.

---

## Problem 4: An Improper Prior and a Generalized Bayes Action

**Problem.** Let

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2),
\]

and use the improper prior

\[
\pi(\mu,\sigma^2)
=\sigma^{-2}\mathbf 1_{\{\sigma^2>0\}}.
\]

Show that the marginal posterior distribution of $\mu$ is a Student $t$ distribution centered at $\overline X$, and obtain the generalized Bayes action for estimating $\mu/\sigma$ under squared-error loss.

??? success "Solution (click to expand)"

    Let

    \[
    v=\sigma^2,
    \qquad
    S=\sum_{i=1}^n(X_i-\overline X)^2=(n-1)s^2,
    \]

    where $s^2$ is the usual sample variance. Assume $n\geq2$ and $S>0$; these conditions hold almost surely under a nondegenerate normal model.

    Using the likelihood and the prior $\pi(\mu,v)\propto v^{-1}$, the joint posterior kernel is

    \[
    \pi(\mu,v\mid X)
    \propto
    v^{-(n/2+1)}
    \exp\left\{-\frac1{2v}
    \sum_{i=1}^n(X_i-\mu)^2\right\}.
    \]

    The identity

    \[
    \sum_{i=1}^n(X_i-\mu)^2
    =S+n(\mu-\overline X)^2
    \]

    yields

    \[
    \pi(\mu,v\mid X)
    \propto
    v^{-(n/2+1)}
    \exp\left\{-\frac{S+n(\mu-\overline X)^2}{2v}\right\}.
    \]

    Integrating out $v$ gives the marginal posterior kernel of $\mu$:

    \[
    \pi(\mu\mid X)
    \propto
    \int_0^\infty
    v^{-(n/2+1)}
    \exp\left\{-\frac{S+n(\mu-\overline X)^2}{2v}\right\}dv.
    \]

    Set

    \[
    A=S+n(\mu-\overline X)^2
    \]

    and make the substitution $w=A/(2v)$. This gives

    \[
    \pi(\mu\mid X)
    \propto A^{-n/2}.
    \]

    Therefore,

    \[
    \pi(\mu\mid X)
    \propto
    \left\{1+
    \frac{n(\mu-\overline X)^2}{(n-1)s^2}
    \right\}^{-n/2}.
    \]

    The location-scale density kernel of a Student $t_\nu(m,a)$ distribution is

    \[
    \left\{1+\frac1\nu
    \left(\frac{\mu-m}{a}\right)^2
    \right\}^{-(\nu+1)/2}.
    \]

    Comparing the two expressions gives

    \[
    \nu=n-1,
    \qquad
    m=\overline X,
    \qquad
    a=\frac{s}{\sqrt n}.
    \]

    Hence

    \[
    \boxed{
    \left.
    \frac{\mu-\overline X}{s/\sqrt n}
    \right|X
    \sim t_{n-1}
    }.
    \]

    Under squared-error loss, the generalized Bayes action for estimating

    \[
    \eta=\frac\mu\sigma=\frac\mu{\sqrt v}
    \]

    is the posterior mean

    \[
    \delta_{\mathrm{GB}}(X)
    =E\left(\left.\frac\mu{\sqrt v}\right|X\right).
    \]

    From the joint posterior distribution,

    \[
    \mu\mid v,X
    \sim N\left(\overline X,\frac vn\right),
    \]

    so

    \[
    E(\mu\mid v,X)=\overline X.
    \]

    Integrating the conditional normal density over $\mu$ also gives

    \[
    \pi(v\mid X)
    \propto
    v^{-\{(n-1)/2+1\}}
    \exp\left(-\frac S{2v}\right).
    \]

    Thus

    \[
    v\mid X\sim\operatorname{Inv-Gamma}(\alpha,\beta),
    \]

    where

    \[
    \alpha=\frac{n-1}{2},
    \qquad
    \beta=\frac S2.
    \]

    We use the parameterization

    \[
    p(v)=\frac{\beta^\alpha}{\Gamma(\alpha)}
    v^{-\alpha-1}e^{-\beta/v},
    \qquad v>0.
    \]

    By iterated expectation,

    \[
    \delta_{\mathrm{GB}}(X)
    =\overline X E(v^{-1/2}\mid X).
    \]

    The negative half-moment of the inverse-gamma distribution is

    \[
    E(v^{-1/2}\mid X)
    =\beta^{-1/2}
    \frac{\Gamma(\alpha+1/2)}{\Gamma(\alpha)}.
    \]

    Substituting $\alpha=(n-1)/2$ and $\beta=S/2$ gives

    \[
    E(v^{-1/2}\mid X)
    =\sqrt{\frac2S}
    \frac{\Gamma(n/2)}{\Gamma((n-1)/2)}.
    \]

    Therefore,

    \[
    \boxed{
    \delta_{\mathrm{GB}}(X)
    =\overline X
    \sqrt{\frac{2}{(n-1)s^2}}
    \frac{\Gamma(n/2)}{\Gamma((n-1)/2)}
    }.
    \]

    When $n\geq2$ and $S>0$, the posterior is proper and the displayed posterior moment is finite, so the generalized Bayes action is well defined.

---

## Problem 5: A Modified James-Stein Estimator

**Problem.** Let

\[
X\sim N_d(\theta,I_d),
\qquad
\theta\in\mathbb R^d,
\]

under squared-error loss. Write

\[
\mathbf 1=(1,\ldots,1)^\top,
\qquad
\overline x=\frac{\langle x,\mathbf 1\rangle}{d}.
\]

Consider the modified James-Stein estimator

\[
\widehat\theta(x)
=\overline x\mathbf 1
+\left\{1-h(\|x-\overline x\mathbf 1\|^2)\right\}
(x-\overline x\mathbf 1).
\]

Equivalently,

\[
\widehat\theta(x)
=x-h(\|x-\overline x\mathbf 1\|^2)
(x-\overline x\mathbf 1).
\]

(a) Use SURE to give an expression for the risk of this estimator, assuming all regularity conditions required for SURE.

(b) Specialize the expression to $h(u)=C/u$, where $C$ is a constant.

(c) For $d\geq4$, choose $C$ so that $\widehat\theta$ dominates the unbiased estimator $\widehat\theta_u(x)=x$. Explain what happens when $d=3$.

(d) Give an empirical Bayes interpretation of $\widehat\theta$.

??? success "Solution (click to expand)"

    Define

    \[
    z=x-\overline x\mathbf 1,
    \qquad
    u=\|z\|^2=\sum_{i=1}^dz_i^2.
    \]

    **(a) The SURE risk expression**

    Write the estimator as

    \[
    \widehat\theta(x)=x+g(x),
    \qquad
    g(x)=-h(u)z.
    \]

    For $X\sim N_d(\theta,I_d)$, Stein's unbiased risk identity gives

    \[
    R(\theta,\widehat\theta)
    =d+E_\theta\left[
    \|g(X)\|^2+2\operatorname{div}g(X)
    \right].
    \]

    First,

    \[
    \|g(x)\|^2
    =h(u)^2\|z\|^2
    =uh(u)^2.
    \]

    Since

    \[
    z_i=x_i-\overline x,
    \qquad
    \overline x=\frac1d\sum_{j=1}^dx_j,
    \]

    we have

    \[
    \frac{\partial z_i}{\partial x_i}=1-\frac1d,
    \qquad
    \sum_{i=1}^dz_i=0.
    \]

    Moreover,

    \[
    \frac{\partial u}{\partial x_i}
    =2\sum_{j=1}^dz_j
    \left(\mathbf 1_{\{i=j\}}-\frac1d\right)
    =2z_i.
    \]

    Therefore, by the product and chain rules,

    \[
    \operatorname{div}g(x)
    =\sum_{i=1}^d
    \frac\partial{\partial x_i}[-h(u)z_i].
    \]

    More explicitly,

    \[
    \operatorname{div}g(x)
    =-2h'(u)\sum_{i=1}^dz_i^2
    -h(u)\sum_{i=1}^d\left(1-\frac1d\right).
    \]

    Thus

    \[
    \operatorname{div}g(x)
    =-2uh'(u)-(d-1)h(u).
    \]

    Substituting into SURE gives

    \[
    \boxed{
    R(\theta,\widehat\theta)
    =d+E_\theta\left[
    uh(u)^2-4uh'(u)-2(d-1)h(u)
    \right]
    },
    \]

    where

    \[
    u=\|X-\overline X\mathbf 1\|^2.
    \]

    **(b) Taking $h(u)=C/u$**

    In this case,

    \[
    h'(u)=-\frac C{u^2}.
    \]

    Hence

    \[
    uh(u)^2=\frac{C^2}{u},
    \qquad
    -4uh'(u)=\frac{4C}{u},
    \qquad
    -2(d-1)h(u)=-\frac{2(d-1)C}{u}.
    \]

    Therefore,

    \[
    \boxed{
    R(\theta,\widehat\theta)
    =d+C\{C-2(d-3)\}
    E_\theta\left[
    \frac1{\|X-\overline X\mathbf 1\|^2}
    \right]
    }.
    \]

    This expression requires the relevant expectation to be finite and the integrability conditions for SURE to hold; in particular, they hold for $d\geq4$.

    **(c) Dominance and the case $d=3$**

    The unbiased estimator $\widehat\theta_u(X)=X$ has constant risk

    \[
    R(\theta,\widehat\theta_u)
    =E_\theta\|X-\theta\|^2
    =d.
    \]

    For $d\geq4$, the residual vector

    \[
    Z=X-\overline X\mathbf 1
    \]

    lies in a subspace of dimension $d-1\geq3$. In an orthonormal basis for that subspace,

    \[
    u=\|Z\|^2
    \sim\chi_{d-1}^{\prime 2}(\lambda),
    \qquad
    \lambda=\|\theta-\overline\theta\mathbf 1\|^2.
    \]

    A noncentral chi-square variable with more than two degrees of freedom satisfies

    \[
    0<E_\theta(1/u)<\infty.
    \]

    Hence the risk difference is

    \[
    R(\theta,\widehat\theta)
    -R(\theta,\widehat\theta_u)
    =C\{C-2(d-3)\}E_\theta(1/u).
    \]

    Choosing

    \[
    0<C<2(d-3)
    \]

    makes the risk difference strictly negative for every $\theta$. The quadratic factor is minimized at

    \[
    C=d-3
    \]

    giving

    \[
    \boxed{
    \widehat\theta(X)
    =\overline X\mathbf 1
    +\left(
    1-\frac{d-3}{\|X-\overline X\mathbf 1\|^2}
    \right)
    (X-\overline X\mathbf 1)
    }.
    \]

    Thus the modified James-Stein rule strictly dominates $X$ under total squared-error risk.

    When $d=3$, the residual subspace has dimension only $2$. In this case,

    \[
    u\sim\chi_2^{\prime 2}(\lambda)
    \]

    and $E_\theta(1/u)=\infty$. For every $C\neq0$,

    \[
    E_\theta\|g(X)\|^2
    =C^2E_\theta(1/u)
    =\infty.
    \]

    Therefore, the formal SURE expression in part (b) cannot be applied when $d=3$ because its regularity and integrability assumptions fail. The only finite-risk member of this particular $C/u$ family is $C=0$, which gives $\widehat\theta=X$ and no strict improvement.

    **(d) Empirical Bayes interpretation**

    Consider the normal-normal hierarchical model

    \[
    X_i\mid\theta_i\sim N(\theta_i,1),
    \qquad
    \theta_i\mid\mu,\tau^2\sim N(\mu,\tau^2),
    \qquad i=1,\ldots,d.
    \]

    Here the $\theta_i$ are assigned a working prior even though they are treated as fixed unknown parameters when the frequentist risk is evaluated. Bayes' rule gives the posterior mean

    \[
    E(\theta_i\mid X_i,\mu,\tau^2)
    =\frac{\tau^2X_i+\mu}{1+\tau^2}.
    \]

    Equivalently,

    \[
    E(\theta_i\mid X_i,\mu,\tau^2)
    =\mu+left(1-\frac1{1+\tau^2}\right)(X_i-\mu).
    \]

    Marginally,

    \[
    X_i\sim N(\mu,1+\tau^2).
    \]

    Estimate the prior center by $\widehat\mu=\overline X$ and define

    \[
    u=\sum_{i=1}^d(X_i-\overline X)^2.
    \]

    Use the plug-in estimate of the marginal variance

    \[
    1+\widehat\tau^2=\frac{u}{d-3},
    \qquad d\geq4.
    \]

    Then

    \[
    \widehat\theta_i^{\mathrm{EB}}
    =\overline X
    +\left(1-\frac{d-3}{u}\right)(X_i-\overline X).
    \]

    In vector form,

    \[
    \boxed{
    \widehat\theta^{\mathrm{EB}}(X)
    =\overline X\mathbf 1
    +\left(
    1-\frac{d-3}{\|X-\overline X\mathbf 1\|^2}
    \right)
    (X-\overline X\mathbf 1)
    }.
    \]

    This is exactly the modified James-Stein estimator with $C=d-3$.
