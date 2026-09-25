---
tags:
  - Advanced Mathematical Statistics
  - Statistical Inference
  - Homework Solutions
---

# Advanced Mathematical Statistics Assignment 1

---

## Problem 1: Unbiasedness and Squared-Error Risk

**Problem.** Let

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(\theta).
\]

Show that $\overline X$ is unbiased for $\theta$ and compute its risk under squared-error loss.

??? success "Solution (click to expand)"

    Since

    \[
    E_\theta(X_i)=\theta,
    \]

    we have

    \[
    E_\theta(\overline X)
    =E_\theta\left(\frac{1}{n}\sum_{i=1}^nX_i\right)
    =\frac{1}{n}\sum_{i=1}^nE_\theta(X_i)
    =\theta.
    \]

    Hence $\overline X$ is unbiased for $\theta$.

    Under squared-error loss,

    \[
    L(\theta,a)=(a-\theta)^2,
    \]

    the risk of $\overline X$ is

    \[
    \begin{aligned}
    R(\theta,\overline X)
    &=E_\theta(\overline X-\theta)^2\\
    &=\operatorname{Var}_\theta(\overline X)
      +\{E_\theta(\overline X)-\theta\}^2\\
    &=\operatorname{Var}_\theta(\overline X).
    \end{aligned}
    \]

    By independence,

    \[
    \begin{aligned}
    \operatorname{Var}_\theta(\overline X)
    &=\frac{1}{n^2}\sum_{i=1}^n
    \operatorname{Var}_\theta(X_i)\\
    &=\frac{1}{n^2}n\theta(1-\theta)\\
    &=\frac{\theta(1-\theta)}{n}.
    \end{aligned}
    \]

    Therefore,

    \[
    \boxed{R(\theta,\overline X)
    =\frac{\theta(1-\theta)}{n}}.
    \]

---

## Problem 2: Student's (t) Distribution

**Problem.** Let $X$ have Student's $t$ distribution with $p$ degrees of freedom.

(a) Derive the mean and variance of $X$.

(b) Show that $X^2$ has an $F$ distribution with $1$ and $p$ degrees of freedom.

(c) Let $f_X(x\mid p)$ denote the density of $X$. Show that, for each fixed $x$,

\[
\lim_{p\to\infty}f_X(x\mid p)
=\frac{1}{\sqrt{2\pi}}e^{-x^2/2}.
\]

??? success "Solution (click to expand)"

    Use the standard representation

    \[
    X=\frac{Z}{\sqrt{W/p}},
    \qquad
    Z\sim N(0,1),
    \qquad
    W\sim\chi_p^2,
    \qquad
    Z\perp W.
    \]

    **(a) Mean and variance**

    For $r<p/2$ and $W\sim\chi_p^2$,

    \[
    E(W^{-r})
    =2^{-r}\frac{\Gamma(p/2-r)}{\Gamma(p/2)}.
    \]

    In particular,

    \[
    E(W^{-1/2})<\infty
    \iff p>1,
    \qquad
    E(W^{-1})<\infty
    \iff p>2.
    \]

    Thus, when $p>1$, $X$ is integrable and

    \[
    E(X)
    =\sqrt{p}\,E(Z)E(W^{-1/2})
    =0.
    \]

    Therefore,

    \[
    \boxed{E(X)=0\quad(p>1)}.
    \]

    For $p\leq1$, the mean does not exist. For $p>2$,

    \[
    \begin{aligned}
    E(X^2)
    &=pE(Z^2)E(W^{-1})\\
    &=p\cdot1\cdot
    \frac{1}{2}\frac{\Gamma(p/2-1)}{\Gamma(p/2)}\\
    &=p\cdot\frac{1}{p-2}\\
    &=\frac{p}{p-2}.
    \end{aligned}
    \]

    Since $E(X)=0$,

    \[
    \boxed{\operatorname{Var}(X)
    =\frac{p}{p-2}\quad(p>2)}.
    \]

    For $1<p\leq2$, the mean exists but the variance is infinite. For $p\leq1$, the usual variance is not defined because the mean does not exist.

    **(b) Distribution of (X^2)**

    Since $Z^2\sim\chi_1^2$ and $Z^2\perp W$,

    \[
    X^2
    =\frac{Z^2}{W/p}
    =\frac{Z^2/1}{W/p}.
    \]

    By the definition of the $F$ distribution,

    \[
    \boxed{X^2\sim F_{1,p}}.
    \]

    **(c) Limiting density**

    The $t_p$ density is

    \[
    f_X(x\mid p)
    =\frac{\Gamma\{(p+1)/2\}}
    {\sqrt{p\pi}\,\Gamma(p/2)}
    \left(1+\frac{x^2}{p}\right)^{-(p+1)/2}.
    \]

    Write

    \[
    A_p
    =\frac{\Gamma\{(p+1)/2\}}
    {\sqrt{p\pi}\,\Gamma(p/2)},
    \qquad
    B_p
    =\left(1+\frac{x^2}{p}\right)^{-(p+1)/2}.
    \]

    Let $z=p/2$. By the standard Gamma-ratio asymptotic formula

    \[
    \frac{\Gamma(z+a)}{\Gamma(z+b)}
    \sim z^{a-b},
    \qquad z\to\infty,
    \]

    we obtain

    \[
    A_p
    \sim\frac{(p/2)^{1/2}}{\sqrt{p\pi}}
    =\frac{1}{\sqrt{2\pi}}.
    \]

    Moreover,

    \[
    \begin{aligned}
    \log B_p
    &=-\frac{p+1}{2}
    \log\left(1+\frac{x^2}{p}\right)\\
    &=-\frac{1}{2}\left(1+\frac{1}{p}\right)
    p\log\left(1+\frac{x^2}{p}\right)\\
    &\longrightarrow-\frac{x^2}{2},
    \end{aligned}
    \]

    because

    \[
    p\log\left(1+\frac{x^2}{p}\right)
    \longrightarrow x^2.
    \]

    Hence

    \[
    B_p\longrightarrow e^{-x^2/2}.
    \]

    Therefore, for every fixed $x$,

    \[
    \boxed{\lim_{p\to\infty}f_X(x\mid p)
    =\frac{1}{\sqrt{2\pi}}e^{-x^2/2}}.
    \]

---

## Problem 3: Sample-Variance Identity and Its Covariance with the Sample Mean

**Problem.** Let $X_1,\ldots,X_n$ be a random sample, and define $\overline X$ and $S^2$ in the usual way.

(a) Show that

\[
S^2
=\frac{1}{2n(n-1)}
\sum_{i=1}^n\sum_{j=1}^n(X_i-X_j)^2.
\]

(b) Assume that the $X_i$ have finite fourth moment. Let

\[
\theta_1=E(X_i),
\qquad
\theta_j=E\{(X_i-\theta_1)^j\},
\qquad j=2,3,4.
\]

Find $\operatorname{Cov}(\overline X,S^2)$ in terms of $\theta_1,\ldots,\theta_4$. Under what conditions is $\operatorname{Cov}(\overline X,S^2)=0$?

??? success "Solution (click to expand)"

    **(a) Pairwise-difference representation**

    Expanding the double sum gives

    \[
    \begin{aligned}
    \sum_{i=1}^n\sum_{j=1}^n(X_i-X_j)^2
    &=\sum_{i,j}(X_i^2-2X_iX_j+X_j^2)\\
    &=2n\sum_{i=1}^nX_i^2
      -2\left(\sum_{i=1}^nX_i\right)^2\\
    &=2n\left(\sum_{i=1}^nX_i^2-n\overline X^2\right)\\
    &=2n\sum_{i=1}^n(X_i-\overline X)^2.
    \end{aligned}
    \]

    Since

    \[
    S^2=\frac{1}{n-1}
    \sum_{i=1}^n(X_i-\overline X)^2,
    \]

    it follows that

    \[
    \boxed{S^2
    =\frac{1}{2n(n-1)}
    \sum_{i=1}^n\sum_{j=1}^n(X_i-X_j)^2}.
    \]

    **(b) Covariance with the sample mean**

    Set

    \[
    Y_i=X_i-\theta_1,
    \qquad
    \overline Y=\overline X-\theta_1.
    \]

    Then

    \[
    E(Y_i)=0,
    \qquad
    E(Y_i^j)=\theta_j,
    \qquad j=2,3,4,
    \]

    and

    \[
    S^2
    =\frac{1}{n-1}
    \left(\sum_{i=1}^nY_i^2-n\overline Y^2\right).
    \]

    Since $E(\overline Y)=0$,

    \[
    \operatorname{Cov}(\overline X,S^2)
    =\operatorname{Cov}(\overline Y,S^2)
    =E(\overline YS^2).
    \]

    First,

    \[
    \begin{aligned}
    E\left(\overline Y\sum_{i=1}^nY_i^2\right)
    &=\frac{1}{n}
    \sum_{j=1}^n\sum_{i=1}^nE(Y_jY_i^2)\\
    &=\frac{1}{n}\sum_{i=1}^nE(Y_i^3)\\
    &=\theta_3,
    \end{aligned}
    \]

    where all terms with $i\neq j$ vanish by independence and centering. Also,

    \[
    \begin{aligned}
    E(\overline Y^3)
    &=\frac{1}{n^3}E\left(\sum_{i=1}^nY_i\right)^3\\
    &=\frac{1}{n^3}\sum_{i=1}^nE(Y_i^3)\\
    &=\frac{\theta_3}{n^2},
    \end{aligned}
    \]

    because every mixed term contains a factor with mean zero. Therefore,

    \[
    \begin{aligned}
    \operatorname{Cov}(\overline X,S^2)
    &=\frac{1}{n-1}
    \left\{\theta_3-n\left(\frac{\theta_3}{n^2}\right)\right\}\\
    &=\frac{\theta_3}{n}.
    \end{aligned}
    \]

    Thus,

    \[
    \boxed{\operatorname{Cov}(\overline X,S^2)
    =\frac{\theta_3}{n}}.
    \]

    Consequently,

    \[
    \operatorname{Cov}(\overline X,S^2)=0
    \iff\theta_3=0.
    \]

    In particular, every distribution symmetric about its mean has $\theta_3=0$, so $\overline X$ and $S^2$ are uncorrelated. Symmetry is sufficient but not necessary.

---

## Problem 4: Sufficiency of |X|

**Problem.** Let $X$ be one observation from a $N(0,\sigma^2)$ population. Is $|X|$ a sufficient statistic for $\sigma^2$?

??? success "Solution (click to expand)"

    The density is

    \[
    f_\sigma(x)
    =\frac{1}{\sqrt{2\pi}\sigma}
    \exp\left(-\frac{x^2}{2\sigma^2}\right).
    \]

    Since $x^2=|x|^2$,

    \[
    f_\sigma(x)
    =\underbrace{
    \frac{1}{\sqrt{2\pi}\sigma}
    \exp\left(-\frac{|x|^2}{2\sigma^2}\right)
    }_{g_\sigma(|x|)}
    \underbrace{1}_{h(x)}.
    \]

    The support is $\mathbb{R}$ and does not depend on $\sigma^2$. Therefore, by the Fisher-Neyman factorization theorem,

    \[
    \boxed{|X|\text{ is sufficient for }\sigma^2}.
    \]

---

## Problem 5: Exponential-Family Representation and Minimal Sufficiency

**Problem.** Let $X_1,\ldots,X_n$ have density

\[
f_{\mu,\sigma}(x)
=\exp\left\{-\left(\frac{x-\mu}{\sigma}\right)^4
-\xi(\mu,\sigma)\right\},
\qquad
\mu\in\mathbb{R},
\quad \sigma>0.
\]

Show that the joint family is an exponential family and identify a minimal sufficient statistic.

??? success "Solution (click to expand)"

    First expand

    \[
    -\left(\frac{x-\mu}{\sigma}\right)^4
    =-\frac{x^4}{\sigma^4}
    +\frac{4\mu x^3}{\sigma^4}
    -\frac{6\mu^2x^2}{\sigma^4}
    +\frac{4\mu^3x}{\sigma^4}
    -\frac{\mu^4}{\sigma^4}.
    \]

    Hence the joint density is

    \[
    \begin{aligned}
    f_{\mu,\sigma}(\mathbf{x})
    &=\prod_{i=1}^nf_{\mu,\sigma}(x_i)\\
    &=\exp\Biggl\{
    -\frac{1}{\sigma^4}\sum_{i=1}^nx_i^4
    +\frac{4\mu}{\sigma^4}\sum_{i=1}^nx_i^3
    -\frac{6\mu^2}{\sigma^4}\sum_{i=1}^nx_i^2\\
    &\qquad
    +\frac{4\mu^3}{\sigma^4}\sum_{i=1}^nx_i
    -n\left(\frac{\mu^4}{\sigma^4}
    +\xi(\mu,\sigma)\right)
    \Biggr\}.
    \end{aligned}
    \]

    Thus the family is a curved exponential family with natural-statistic vector

    \[
    T(\mathbf{X})
    =\left(
    \sum_{i=1}^nX_i,
    \sum_{i=1}^nX_i^2,
    \sum_{i=1}^nX_i^3,
    \sum_{i=1}^nX_i^4
    \right).
    \]

    One corresponding natural-parameter vector is

    \[
    \eta(\mu,\sigma)
    =\left(
    \frac{4\mu^3}{\sigma^4},
    -\frac{6\mu^2}{\sigma^4},
    \frac{4\mu}{\sigma^4},
    -\frac{1}{\sigma^4}
    \right).
    \]

    To prove minimal sufficiency, use the likelihood-ratio criterion. For two sample points $\mathbf{x}$ and $\mathbf{y}$,

    \[
    \begin{aligned}
    \log\frac{f_{\mu,\sigma}(\mathbf{x})}
    {f_{\mu,\sigma}(\mathbf{y})}
    &=-\frac{1}{\sigma^4}
    \sum_{i=1}^n\{(x_i-\mu)^4-(y_i-\mu)^4\}\\
    &=-\frac{1}{\sigma^4}
    \left(D_4-4\mu D_3+6\mu^2D_2-4\mu^3D_1\right),
    \end{aligned}
    \]

    where

    \[
    D_k=\sum_{i=1}^nx_i^k-
    \sum_{i=1}^ny_i^k,
    \qquad k=1,2,3,4.
    \]

    The likelihood ratio is independent of $(\mu,\sigma)$ if and only if

    \[
    D_4-4\mu D_3+6\mu^2D_2-4\mu^3D_1=0
    \]

    for every $\mu\in\mathbb{R}$. A polynomial is identically zero if and only if each coefficient is zero, so

    \[
    D_1=D_2=D_3=D_4=0.
    \]

    This is equivalent to

    \[
    T(\mathbf{x})=T(\mathbf{y}).
    \]

    Therefore, by the likelihood-ratio characterization of minimal sufficiency,

    \[
    \boxed{
    T(\mathbf{X})
    =\left(
    \sum_{i=1}^nX_i,
    \sum_{i=1}^nX_i^2,
    \sum_{i=1}^nX_i^3,
    \sum_{i=1}^nX_i^4
    \right)
    }
    \]

    is minimal sufficient for $(\mu,\sigma)$.

---

## Problem 6: Minimal Sufficient but Not Complete

**Problem.** For a sample from $\operatorname{Uniform}(\theta,\theta+1)$, show that the minimal sufficient statistic $(X_{(1)},X_{(n)})$ is not complete.

??? success "Solution (click to expand)"

    The joint density is

    \[
    \begin{aligned}
    f_\theta(\mathbf{x})
    &=\prod_{i=1}^n
    \mathbf{1}_{\{\theta<x_i<\theta+1\}}\\
    &=\mathbf{1}_{\{\theta<X_{(1)},\ X_{(n)}<\theta+1\}}\\
    &=\mathbf{1}_{\{X_{(n)}-1<\theta<X_{(1)}\}}.
    \end{aligned}
    \]

    Hence, by the Fisher-Neyman factorization theorem,

    \[
    T=(X_{(1)},X_{(n)})
    \]

    is sufficient.

    To prove minimal sufficiency, for a sample point $\mathbf{x}$ define

    \[
    I_\mathbf{x}=(x_{(n)}-1,x_{(1)}).
    \]

    On the union of the supports, namely for sample points satisfying

    \[
    x_{(n)}-x_{(1)}\leq1,
    \]

    we have

    \[
    f_\theta(\mathbf{x})
    =\mathbf{1}_{\{\theta\in I_\mathbf{x}\}}.
    \]

    For two sample points $\mathbf{x}$ and $\mathbf{y}$, the likelihood ratio is independent of $\theta$ precisely when

    \[
    I_\mathbf{x}=I_\mathbf{y},
    \]

    which is equivalent to

    \[
    x_{(1)}=y_{(1)},
    \qquad
    x_{(n)}=y_{(n)}.
    \]

    Thus the likelihood-ratio criterion shows that $(X_{(1)},X_{(n)})$ is minimal sufficient.

    It remains to show that it is not complete. First suppose $n\geq2$ and consider the range

    \[
    R=X_{(n)}-X_{(1)},
    \]

    which is a measurable function of $T$. Write

    \[
    X_i=\theta+U_i,
    \qquad
    U_i\overset{\mathrm{iid}}{\sim}\operatorname{Uniform}(0,1).
    \]

    Then

    \[
    R=U_{(n)}-U_{(1)},
    \]

    so its distribution does not depend on $\theta$. Moreover,

    \[
    E\{U_{(n)}\}=\frac{n}{n+1},
    \qquad
    E\{U_{(1)}\}=\frac{1}{n+1},
    \]

    and hence

    \[
    E_\theta(R)=\frac{n-1}{n+1}.
    \]

    Define

    \[
    g(T)
    =X_{(n)}-X_{(1)}-\frac{n-1}{n+1}.
    \]

    Then, for every $\theta$,

    \[
    E_\theta\{g(T)\}=0,
    \]

    while $g(T)$ is not almost surely zero because the range is nondegenerate. Therefore $T$ is not complete. Consequently,

    \[
    \boxed{(X_{(1)},X_{(n)})
    \text{ is minimal sufficient but not complete}}.
    \]

    If $n=1$, noncompleteness can be shown separately by taking

    \[
    g(X_1)=\sin(2\pi X_1).
    \]

    Indeed,

    \[
    E_\theta\{g(X_1)\}
    =\int_\theta^{\theta+1}\sin(2\pi x)\,dx
    =0
    \]

    for every $\theta$, whereas $g(X_1)$ is not almost surely zero.
