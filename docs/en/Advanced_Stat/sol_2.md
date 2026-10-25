---
tags:
  - Advanced Mathematical Statistics
  - Completeness
  - UMVUE
  - U-Statistics
  - Information Inequality
  - Homework
---

# Advanced Mathematical Statistics: Assignment 2


## Problem 1: Binomial UMVUEs

**Problem.** Let

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(p),
\qquad p\in(0,1),
\]

and define

\[
T=\sum_{i=1}^nX_i.
\]

(i) Find the UMVUE of $p^m$, where $m\leq n$.

(ii) Find the UMVUE of $P(X_1+\cdots+X_m=k)$, where $k\leq m\leq n$.

??? success "Solution (click to expand)"

    The joint probability mass function is

    \[
    f_p(x_1,\ldots,x_n)
    =p^{\sum_{i=1}^n x_i}(1-p)^{n-\sum_{i=1}^n x_i}.
    \]

    By the factorization theorem, $T=\sum_{i=1}^nX_i$ is sufficient for $p$. Moreover,

    \[
    T\sim\operatorname{Binomial}(n,p)
    \]

    The binomial family is a full one-parameter exponential family, so $T$ is complete. Therefore, by the Lehmann-Scheffé theorem, every unbiased function of $T$ is the UMVUE of its expectation.

    **(i) UMVUE of $p^m$**

    Let

    \[
    W=X_1X_2\cdots X_m.
    \]

    Then

    \[
    E_p(W)=p^m.
    \]

    Rao-Blackwellize $W$ with respect to $T$. For $t=0,\ldots,n$,

    \[
    E_p(W\mid T=t)
    =P(X_1=\cdots=X_m=1\mid T=t).
    \]

    Conditional on a total of $t$ successes, all placements of the successes among the $n$ positions are equally likely. Hence

    \[
    E_p(W\mid T=t)
    =\frac{\binom{n-m}{t-m}}{\binom nt}
    =\frac{\binom tm}{\binom nm}.
    \]

    A binomial coefficient is understood to be zero when its lower index is outside its admissible range. Thus

    \[
    \boxed{
    \widehat{p^m}_{\mathrm{UMVUE}}
    =\frac{\binom Tm}{\binom nm}
    }.
    \]

    **(ii) UMVUE of $P(X_1+\cdots+X_m=k)$**

    Let

    \[
    W=\mathbf 1_{\{X_1+\cdots+X_m=k\}}.
    \]

    Then

    \[
    E_p(W)
    =P(X_1+\cdots+X_m=k)
    =\binom mkp^k(1-p)^{m-k}.
    \]

    Conditional on $T=t$, the number of successes among the first $m$ observations has a hypergeometric distribution. Therefore,

    \[
    E_p(W\mid T=t)
    =\frac{\binom mk\binom{n-m}{t-k}}{\binom nt}.
    \]

    Consequently,

    \[
    \boxed{
    \widehat{P(X_1+\cdots+X_m=k)}_{\mathrm{UMVUE}}
    =\frac{\binom mk\binom{n-m}{T-k}}{\binom nT}
    }.
    \]

---

## Problem 2: A Poisson UMVUE for a Transformed Parameter

**Problem.** Let

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Pois}(\theta),
\]

define

\[
S=\sum_{i=1}^nX_i,
\]

and fix $a>0$. Find the UMVUE of

\[
    g(\theta)=e^{-a\theta}
\]

??? success "Solution (click to expand)"

    The joint probability mass function can be written as

    \[
    f_\theta(x_1,\ldots,x_n)
    =\left(\prod_{i=1}^n\frac1{x_i!}\right)
    e^{-n\theta}\theta^{\sum_{i=1}^nx_i}.
    \]

    Thus $S$ is sufficient. Moreover,

    \[
    S\sim\operatorname{Pois}(n\theta),
    \]

    and the Poisson family is complete. Hence $S$ is complete and sufficient.

    For every real number $q$, the probability generating function of a Poisson random variable gives

    \[
    E_\theta(q^S)=\exp\{n\theta(q-1)\}.
    \]

    Choose

    \[
    q=1-\frac an
    \]

    Then

    \[
    E_\theta\left[\left(1-\frac an\right)^S\right]
    =e^{-a\theta}.
    \]

    Therefore, by the Lehmann-Scheffé theorem,

    \[
    \boxed{
    \widehat{e^{-a\theta}}_{\mathrm{UMVUE}}
    =\left(1-\frac an\right)^S
    }.
    \]

    When $a=n$ and $S=0$, the value $0^0$ is interpreted as $1$, consistently with the power-series representation of the probability generating function.

---

## Problem 3: Best Unbiased Estimators in Endpoint Families

**Problem.** Let $X_1,\ldots,X_n$ be i.i.d. with density

\[
f(x\mid\theta)=c(\theta)m(x)\mathbf 1_{\{a<x<\theta\}},
\]

where

\[
M(\theta)=\int_a^\theta m(x)\,dx=\frac1{c(\theta)}.
\]

Let $Y=X_{(n)}$. Its density is

\[
f_Y(y\mid\theta)
=nm(y)c(\theta)^nM(y)^{n-1}\mathbf 1_{\{a<y<\theta\}}.
\]

(i) Use the factorization theorem to show that $Y$ is sufficient. Then prove directly from the definition of completeness that $Y$ is complete.

(ii) Let $h$ be differentiable. Write the unbiasedness equation for an estimator $\delta(Y)$ of $h(\theta)$, differentiate that equation with respect to $\theta$, and thereby obtain the UMVUE of $h(\theta)$.

(iii) Specialize the result to $X_i\sim\operatorname{Unif}(0,\theta)$ and $h(\theta)=\theta^r$, where $r>0$. Explain explicitly why $(n+r)Y^r/n$ is a special case of the general formula.

??? success "Solution (click to expand)"

    **(i) Sufficiency and completeness**

    The joint density is

    \[
    f_\theta(x_1,\ldots,x_n)
    =c(\theta)^n\prod_{i=1}^nm(x_i)
    \prod_{i=1}^n\mathbf 1_{\{a<x_i<\theta\}}.
    \]

    Rewrite it as

    \[
    f_\theta(x_1,\ldots,x_n)
    =\left(\prod_{i=1}^nm(x_i)\mathbf 1_{\{x_i>a\}}\right)
    c(\theta)^n\mathbf 1_{\{x_{(n)}<\theta\}}.
    \]

    By the factorization theorem, $Y=X_{(n)}$ is sufficient for $\theta$.

    To prove completeness, suppose that a measurable function $q$ satisfies

    \[
    E_\theta[q(Y)]=0
    \qquad\text{for every }\theta>a.
    \]

    Using the density of $Y$, we obtain

    \[
    nc(\theta)^n
    \int_a^\theta q(y)m(y)M(y)^{n-1}\,dy=0.
    \]

    Since $c(\theta)>0$, this is equivalent to

    \[
    \int_a^\theta q(y)m(y)M(y)^{n-1}\,dy=0
    \qquad\text{for every }\theta>a.
    \]

    Differentiating with respect to $\theta$ gives

    \[
    q(\theta)m(\theta)M(\theta)^{n-1}=0
    \]

    almost everywhere. On the support, $m(\theta)>0$ and $M(\theta)>0$, so $q(\theta)=0$ almost everywhere. Hence $Y$ is complete.

    **(ii) General UMVUE formula**

    If $\delta(Y)$ is unbiased for $h(\theta)$, then

    \[
    h(\theta)
    =nc(\theta)^n
    \int_a^\theta\delta(y)m(y)M(y)^{n-1}\,dy.
    \]

    Since $c(\theta)=1/M(\theta)$,

    \[
    h(\theta)M(\theta)^n
    =n\int_a^\theta\delta(y)m(y)M(y)^{n-1}\,dy.
    \]

    Differentiating both sides with respect to $\theta$ yields

    \[
    h'(\theta)M(\theta)^n
    +nh(\theta)M(\theta)^{n-1}m(\theta)
    =n\delta(\theta)m(\theta)M(\theta)^{n-1}.
    \]

    Therefore,

    \[
    \delta(\theta)
    =h(\theta)+\frac{M(\theta)}{nm(\theta)}h'(\theta).
    \]

    Replacing $\theta$ by $Y$ gives the unbiased estimator

    \[
    \delta(Y)
    =h(Y)+\frac{M(Y)}{nm(Y)}h'(Y).
    \]

    Since $Y$ is complete and sufficient, the Lehmann-Scheffé theorem gives

    \[
    \boxed{
    \widehat{h(\theta)}_{\mathrm{UMVUE}}
    =h(Y)+\frac{M(Y)}{nm(Y)}h'(Y)
    }.
    \]

    **(iii) The uniform special case**

    For $X_i\sim\operatorname{Unif}(0,\theta)$,

    \[
    a=0,
    \qquad
    m(x)=1,
    \qquad
    M(y)=\int_0^y1\,dx=y.
    \]

    When $h(\theta)=\theta^r$,

    \[
    h'(y)=ry^{r-1}.
    \]

    Substituting into the general formula gives

    \[
    \delta(Y)
    =Y^r+\frac Yn rY^{r-1}
    =\left(1+\frac rn\right)Y^r.
    \]

    Hence

    \[
    \boxed{
    \widehat{\theta^r}_{\mathrm{UMVUE}}
    =\frac{n+r}{n}Y^r
    }.
    \]

---

## Problem 4: Bernoulli Sample Mean and the Cramér-Rao Lower Bound

**Problem.** Let

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(p).
\]

Show that $\overline X$ is unbiased for $p$ and attains the Cramér-Rao lower bound.

??? success "Solution (click to expand)"

    Since $E_p(X_i)=p$,

    \[
    E_p(\overline X)
    =\frac1n\sum_{i=1}^nE_p(X_i)
    =p.
    \]

    Thus $\overline X$ is unbiased. By independence,

    \[
    \operatorname{Var}_p(\overline X)
    =\frac1{n^2}\sum_{i=1}^n\operatorname{Var}_p(X_i)
    =\frac{p(1-p)}n.
    \]

    The log-likelihood is

    \[
    \ell(p)
    =\left(\sum_{i=1}^nX_i\right)\log p
    +\left(n-\sum_{i=1}^nX_i\right)\log(1-p).
    \]

    Its second derivative is

    \[
    \ell''(p)
    =-\frac{\sum_{i=1}^nX_i}{p^2}
    -\frac{n-\sum_{i=1}^nX_i}{(1-p)^2}.
    \]

    The Fisher information in the sample is

    \[
    I_n(p)
    =-E_p[\ell''(p)]
    =\frac{np}{p^2}+\frac{n(1-p)}{(1-p)^2}
    =\frac{n}{p(1-p)}.
    \]

    For every unbiased estimator of $p$, the Cramér-Rao lower bound is

    \[
    \frac1{I_n(p)}=\frac{p(1-p)}n.
    \]

    Since

    \[
    \operatorname{Var}_p(\overline X)
    =\frac{p(1-p)}n
    =\frac1{I_n(p)},
    \]

    we conclude that

    \[
    \boxed{\overline X\text{ attains the Cramér-Rao lower bound}}.
    \]

---

## Problem 5: A Concrete U-Statistic for a Product of Moments

**Problem.** Let $X_1,\ldots,X_n$ be i.i.d. with

\[
E(X_1^2)<\infty,
\qquad
E(X_1^{-2})<\infty.
\]

Define

\[
\mu=E(X_1),
\qquad
\overline\mu=E(X_1^{-1}).
\]

Find a U-statistic that is unbiased for $\mu\overline\mu$.

??? success "Solution (click to expand)"

    For two independent observations, define the symmetric kernel

    \[
    h(x_1,x_2)
    =\frac12\left(\frac{x_1}{x_2}+\frac{x_2}{x_1}\right).
    \]

    The moment assumptions ensure that the required expectations are finite. By independence,

    \[
    E[h(X_1,X_2)]
    =\frac12\left\{E(X_1)E(X_2^{-1})
    +E(X_2)E(X_1^{-1})\right\}.
    \]

    Hence

    \[
    E[h(X_1,X_2)]=\mu\overline\mu.
    \]

    Therefore, the order-two U-statistic associated with $h$,

    \[
    \boxed{
    U_n
    =\binom n2^{-1}
    \sum_{1\leq i<j\leq n}
    \frac12\left(\frac{X_i}{X_j}+\frac{X_j}{X_i}\right)
    }
    \]

    is unbiased for $\mu\overline\mu$.

---

## Problem 6: Why the Cramér-Rao Bound Fails in a Nonregular Model

**Problem.** Let

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Unif}(0,\theta),
\qquad \theta>0,
\]

and let $Y=X_{(n)}$. Recall that

\[
f_Y(y\mid\theta)
=\frac{ny^{n-1}}{\theta^n}\mathbf 1_{\{0<y<\theta\}}.
\]

(i) Show why the usual interchange of differentiation and integration fails when applied to $E_\theta(Y)$.

(ii) Verify that

\[
\widehat\theta=\frac{n+1}{n}Y
\]

is unbiased. Using the completeness and sufficiency result for the endpoint family from Problem 3, explain why it is the UMVUE. Compute its variance and show that the formal regular Fisher-information inequality does not hold.

(iii) Identify precisely which regularity condition in the usual Cramér-Rao proof fails, and explain the role of the moving support.

??? success "Solution (click to expand)"

    **(i) Failure of the interchange of differentiation and integration**

    Direct calculation gives

    \[
    E_\theta(Y)
    =\int_0^\theta y\frac{ny^{n-1}}{\theta^n}\,dy
    =\frac{n}{n+1}\theta.
    \]

    Therefore,

    \[
    \frac d{d\theta}E_\theta(Y)=\frac{n}{n+1}.
    \]

    If one ignores the dependence of the upper integration limit on $\theta$ and differentiates only the density, then

    \[
    \int_0^\theta y\frac\partial{\partial\theta}
    \left(\frac{ny^{n-1}}{\theta^n}\right)dy
    =-\frac{n^2}{n+1},
    \]

    which is not equal to $n/(n+1)$. The missing term is the boundary contribution. By Leibniz's rule,

    \[
    \frac d{d\theta}
    \int_0^\theta yf_Y(y\mid\theta)\,dy
    =\int_0^\theta y\frac\partial{\partial\theta}
    f_Y(y\mid\theta)\,dy
    +\theta f_Y(\theta\mid\theta).
    \]

    The right-hand side equals

    \[
    -\frac{n^2}{n+1}+n
    =\frac{n}{n+1}.
    \]

    **(ii) UMVUE, variance, and failure of the formal CRLB**

    Since

    \[
    E_\theta(Y)=\frac{n}{n+1}\theta,
    \]

    we have

    \[
    E_\theta\left(\frac{n+1}{n}Y\right)=\theta.
    \]

    Thus $\widehat\theta=(n+1)Y/n$ is unbiased. By Problem 3, $Y$ is complete and sufficient in this endpoint family, so the Lehmann-Scheffé theorem shows that $\widehat\theta$ is the UMVUE.

    Next,

    \[
    E_\theta(Y^2)
    =\int_0^\theta y^2\frac{ny^{n-1}}{\theta^n}\,dy
    =\frac{n}{n+2}\theta^2.
    \]

    Therefore,

    \[
    \operatorname{Var}_\theta(Y)
    =\frac{n}{n+2}\theta^2
    -\left(\frac{n}{n+1}\theta\right)^2
    =\frac{n}{(n+2)(n+1)^2}\theta^2.
    \]

    Hence

    \[
    \boxed{
    \operatorname{Var}_\theta(\widehat\theta)
    =\frac{\theta^2}{n(n+2)}
    }.
    \]

    If the regular-model Fisher-information calculation is applied formally to one observation, then inside the support,

    \[
    \frac\partial{\partial\theta}\log f(X\mid\theta)
    =-\frac1\theta.
    \]

    Thus the formal Fisher information in one observation is

    \[
    I_1^{\mathrm{formal}}(\theta)
    =E_\theta\left[
    \left\{\frac\partial{\partial\theta}
    \log f(X\mid\theta)\right\}^2
    \right]
    =\frac1{\theta^2}.
    \]

    Formal additivity of information for a regular model would give

    \[
    I_n^{\mathrm{formal}}(\theta)=\frac n{\theta^2},
    \]

    and hence the formal Cramér-Rao lower bound

    \[
    \frac1{I_n^{\mathrm{formal}}(\theta)}
    =\frac{\theta^2}{n}.
    \]

    However,

    \[
    \frac{\theta^2}{n(n+2)}<\frac{\theta^2}{n}.
    \]

    Thus the formal regular Fisher-information inequality is “violated.” This is not a contradiction because the model is nonregular.

    **(iii) The failed regularity condition**

    The key failure is that the support of the density depends on the parameter, invalidating the usual interchange of differentiation and integration. Here the support is $(0,\theta)$ and moves with $\theta$.

    For one observation,

    \[
    \int_0^\theta f(x\mid\theta)\,dx=1.
    \]

    Differentiating only the integrand gives

    \[
    \int_0^\theta
    \frac\partial{\partial\theta}f(x\mid\theta)\,dx
    =\int_0^\theta\left(-\frac1{\theta^2}\right)dx
    =-\frac1\theta\neq0.
    \]

    The missing boundary term is

    \[
    f(\theta\mid\theta)=\frac1\theta,
    \]

    which restores the derivative of the total probability to zero. Equivalently, the usual identity

    \[
    E_\theta\left[
    \frac\partial{\partial\theta}\log f(X\mid\theta)
    \right]=0
    \]

    fails here. Therefore, the standard Cramér-Rao proof does not apply.
