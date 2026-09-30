# Chapter 5: Admissibility, Shrinkage, and the James–Stein Phenomenon

This chapter revisits optimality in estimation from a decision-theoretic perspective. Earlier chapters focused mainly on optimality within the class of unbiased estimators, such as MVUE and UMVU theory. Here we remove the unbiasedness restriction and compare estimators through their full risk functions. This leads to dominance, admissibility, Rao–Blackwell risk improvement, shrinkage estimation, invariance, and the celebrated James–Stein phenomenon.

---

## 1. From Unbiasedness to Risk Optimality

Under squared-error loss, the risk of an estimator $\delta(X)$ decomposes as

\[
R(\theta,\delta)
=\operatorname{Var}_\theta\{\delta(X)\}
+\operatorname{Bias}_\theta\{\delta(X)\}^2.
\]

If $\delta$ is unbiased, its risk equals its variance. Thus MVUE/UMVU theory compares variances within the class of unbiased estimators.

Unbiasedness, however, is only a constraint and is not the same as global risk optimality. An estimator with a small bias may have a substantially smaller variance and therefore a smaller mean squared error.

!!! note "Central question of this chapter"

    Instead of asking

    > Which unbiased estimator has the smallest variance?

    we ask

    > Is there another estimator whose risk is never larger and is strictly smaller for at least some parameter values?

This is the starting point of admissibility and shrinkage theory.

---

## 2. Decision-Theoretic Framework

Let

\[
X\sim P_\theta,\qquad \theta\in\Theta,
\]

and let the action space be $\mathcal D$. A decision rule is a measurable function

\[
\delta:\mathcal X\to\mathcal D.
\]

Given a loss function $L(\theta,a)$, the risk function of $\delta$ is

\[
R(\theta,\delta)
=E_\theta L\{\theta,\delta(X)\}.
\]

For point estimation of $g(\theta)$ under squared-error loss,

\[
L(\theta,a)=\|a-g(\theta)\|^2,
\]

so

\[
R(\theta,\delta)
=E_\theta\|\delta(X)-g(\theta)\|^2.
\]

The object being compared is therefore not a single number but the entire risk function

\[
\theta\longmapsto R(\theta,\delta).
\]

### 2.1 Dominance and admissibility

!!! info "Definition 2.1 (Dominance and admissibility)"

    An estimator $\delta_1$ **dominates** another estimator $\delta_0$ if

    \[
    R(\theta,\delta_1)
    \leq R(\theta,\delta_0),
    \qquad \forall\theta\in\Theta,
    \]

    with strict inequality for at least one value of $\theta$.

    An estimator is **inadmissible** if it is dominated by another estimator.

    If no estimator dominates it, it is **admissible**.

!!! note "Operational meaning"

    To prove that $\delta_0$ is inadmissible, it is enough to construct one estimator $\delta_1$ whose risk is no larger everywhere and strictly smaller somewhere.

    Proving admissibility is generally harder because every possible dominating estimator must be ruled out.

---

## 3. Sufficiency as a Decision-Theoretic Improvement

Lecture 4 used Rao–Blackwellization to reduce the variance of an unbiased estimator. The same idea is more general: under a convex loss, conditioning can improve the risk of a decision rule whether or not the original rule is unbiased.

!!! success "Proposition 3.1 (Rao–Blackwell improvement under convex loss)"

    Let $T=T(X)$ be sufficient for $\theta$. Assume that the action space is convex and that, for each fixed $\theta$,

    \[
    a\longmapsto L(\theta,a)
    \]

    is convex. For any decision rule $\delta(X)$ with finite risk, define

    \[
    \delta^*(T)
    =E_\theta\{\delta(X)\mid T\}.
    \]

    Because $T$ is sufficient, the conditional law of $X$ given $T$ does not depend on $\theta$, so $\delta^*(T)$ is a genuine statistic rather than a parameter-dependent expression.

    Then

    \[
    R(\theta,\delta^*)
    \leq R(\theta,\delta),
    \qquad \forall\theta\in\Theta.
    \]

    If the loss is strictly convex and, for some $\theta_0$, $\delta(X)$ is not almost surely a function of $T$, then the inequality is strict at $\theta_0$ and $\delta$ is inadmissible.

??? proof "Proof of Proposition 3.1 (click to expand)"

    Conditional Jensen's inequality gives

    \[
    L\!\left(\theta,E_\theta\{\delta(X)\mid T\}\right)
    \leq
    E_\theta\{L(\theta,\delta(X))\mid T\}.
    \]

    Taking expectations yields

    \[
    R(\theta,\delta^*)
    \leq R(\theta,\delta).
    \]

    Under strict convexity, the inequality is strict whenever the conditional distribution of $\delta(X)$ given $T$ is nondegenerate on a set of positive probability. $\square$

!!! example "Example 3.2 (An inadmissible rule that ignores sufficiency)"

    Let

    \[
    X_1,X_2\overset{\mathrm{iid}}\sim N(\theta,1),
    \]

    and estimate $\theta$ under squared-error loss.

    The statistic

    \[
    T=\bar X
    \]

    is sufficient, while $\delta(X)=X_1$ is not a function of $T$. Moreover,

    \[
    E_\theta(X_1\mid\bar X)=\bar X.
    \]

    Therefore $\bar X$ dominates $X_1$. Directly,

    \[
    R(\theta,X_1)=1,
    \qquad
    R(\theta,\bar X)=\frac12.
    \]

    The main point is not the arithmetic: under a strictly convex loss, a rule that retains irrelevant conditional noise cannot be admissible.

!!! note "Why this matters beyond classical theory"

    Conditioning away ancillary or algorithmic randomness reappears in conditional Monte Carlo, variance reduction, and the construction of stabilized estimators.

---

## 4. Why UMVU Does Not Imply Admissibility

An estimator may be optimal among all unbiased estimators and still be dominated when biased estimators are allowed.

### 4.1 The classical normal-variance example

!!! example "Example 4.1 (Estimating a normal variance)"

    Let

    \[
    X_1,\ldots,X_n
    \overset{\mathrm{iid}}\sim N(\mu,\sigma^2),
    \]

    where both $\mu$ and $\sigma^2$ are unknown. The usual unbiased variance estimator is

    \[
    S^2
    =\frac{1}{n-1}
    \sum_{i=1}^n(X_i-\bar X)^2.
    \]

    It is the UMVU estimator of $\sigma^2$.

Under squared-error loss for estimating $\sigma^2$,

\[
L(\sigma^2,a)=(a-\sigma^2)^2,
\]

we have

\[
\frac{(n-1)S^2}{\sigma^2}
\sim\chi^2_{n-1},
\]

and

\[
\operatorname{Var}(S^2)
=\frac{2\sigma^4}{n-1}.
\]

Since $S^2$ is unbiased,

\[
R(\sigma^2,S^2)
=\frac{2\sigma^4}{n-1}.
\]

Now consider the biased estimator

\[
cS^2,\qquad c>0.
\]

Its risk is

\[
\begin{aligned}
R(\sigma^2,cS^2)
&=E(cS^2-\sigma^2)^2\\
&=c^2\operatorname{Var}(S^2) + \{E(cS^2)-\sigma^2\}^2\\
&=\sigma^4
\left\{
\frac{2c^2}{n-1}+(c-1)^2
\right\}.
\end{aligned}
\]

Minimizing this quadratic in $c$ gives

\[
c^*=\frac{n-1}{n+1}.
\]

Hence

\[
\frac{n-1}{n+1}S^2
=\frac{1}{n+1}
\sum_{i=1}^n(X_i-\bar X)^2
\]

has smaller squared-error risk than $S^2$ for every $\sigma^2>0$.

!!! warning "Conclusion"

    Although $S^2$ is UMVU for $\sigma^2$, it is inadmissible under squared-error loss.

    Thus

    \[
    \boxed{\text{UMVU optimality is only optimality within the unbiased class.}}
    \]

The improvement comes from shrinkage toward zero: a small bias is introduced in exchange for a larger variance reduction.

### 4.2 Stein's second improvement

!!! example "Example 4.2 (Stein's 1964 refinement)"

    The estimator

    \[
    \delta_0
    =\frac{1}{n+1}
    \sum_{i=1}^n(X_i-\bar X)^2
    \]

    is itself inadmissible. Stein (1964) showed that

    \[
    \delta_{\mathrm{St}}
    =\min\left\{
    \frac{1}{n+1}\sum_{i=1}^n(X_i-\bar X)^2,
    \frac{1}{n+2}\sum_{i=1}^nX_i^2
    \right\}
    \]

    has smaller squared-error risk for estimating $\sigma^2$.

The second term uses evidence that the total second moment is small and prevents the variance estimate from being implausibly large relative to it.

!!! note "Teaching priority"

    The first domination calculation for $S^2$ is core material.

    Stein's 1964 refinement mainly illustrates that one improvement need not be the endpoint of the admissibility story.

!!! note "Meaning of shrinkage toward zero"

    Shrinkage toward zero does not claim that the true $\sigma^2$ is close to zero. It means that the estimator is pulled toward a lower-variance target to reduce mean squared error. The cost is bias; the benefit is variance reduction.

---

## 5. Shrinkage as a Bias–Variance Tradeoff

Consider an estimator $\delta(X)$ of a scalar target $g(\theta)$ and the shrinkage rule

\[
\delta_c(X)=c\delta(X),
\qquad 0<c<1.
\]

If $\delta$ is unbiased, then

\[
\operatorname{Bias}_\theta(\delta_c)
=(c-1)g(\theta),
\]

and

\[
\operatorname{Var}_\theta(\delta_c)
=c^2\operatorname{Var}_\theta(\delta).
\]

Therefore

\[
R(\theta,\delta_c)-R(\theta,\delta)
=(c^2-1)\operatorname{Var}_\theta(\delta)
+(c-1)^2g(\theta)^2.
\]

Shrinkage helps when the variance reduction dominates the squared-bias cost.

### 5.1 Fixed shrinkage in the one-dimensional normal mean problem

!!! example "Example 5.1 (One-dimensional normal mean)"

    Let

    \[
    X\sim N(\theta,1),
    \]

    and estimate $\theta$ under squared-error loss.

    The usual estimator $X$ has constant risk

    \[
    R(\theta,X)=1.
    \]

    For fixed $0<a<1$, the shrinkage estimator $aX$ has risk

    \[
    R(\theta,aX)
    =a^2+(1-a)^2\theta^2.
    \]

The estimator $aX$ improves on $X$ near $\theta=0$, but for large $|\theta|$ the squared-bias term dominates.

Hence no fixed $a<1$ uniformly dominates $X$ in one dimension.

!!! note "Lesson"

    Uniform improvement generally requires the amount of shrinkage to adapt to the observed data rather than remain fixed.

### 5.2 Affine rules in the scalar normal model

Let

\[
X\sim N(\theta,\sigma^2),
\]

where $\sigma^2$ is known, and estimate $\theta\in\mathbb R$ under squared-error loss. Consider

\[
\delta_{a,\beta}(X)=aX+\beta.
\]

!!! success "Proposition 5.2 (Classification of affine rules in the scalar normal model)"

    For $\delta_{a,\beta}(X)=aX+\beta$:

    | Parameter values | Admissibility |
    |---|---|
    | $a>1$, $a<0$, or $a=1,\beta\neq0$ | Inadmissible |
    | $a=0$, $0<a<1$, or $a=1,\beta=0$ | Admissible |

Its risk is

\[
R(\theta,aX+\beta)
=a^2\sigma^2+igl\{(a-1)\theta+\beta\bigr\}^2.
\]

For example, when $a=1$ and $\beta\neq0$, $X+\beta$ is dominated by $X$.

For $0<a<1$, the affine rule is the posterior mean under a proper normal prior. Specifically, if

\[
m=\frac{\beta}{1-a},
\qquad
\tau^2=\frac{a}{1-a}\sigma^2,
\]

then $aX+\beta$ can be interpreted as shrinkage toward the prior center $m$.

!!! note "Contraction versus expansion"

    Shrinkage with $0<a<1$ often admits a Bayes interpretation, whereas expansion with $a>1$ is difficult to justify uniformly.

---

## 6. Equivariance and Invariant Risk

Many statistical problems possess symmetries. In such settings, it is natural to consider procedures that respect the same symmetries.

### 6.1 Location equivariance

Let

\[
X\sim N(\theta,\sigma^2).
\]

If the data are shifted by $a$, the parameter is shifted by the same amount. A location-equivariant estimator satisfies

\[
\delta(x+a)=\delta(x)+a.
\]

The usual estimator

\[
\delta(x)=x
\]

is location equivariant.

!!! info "Definition 6.1 (Equivariant estimator)"

    Suppose a group $G$ acts on both the sample space and the parameter space. An estimator $\delta$ is equivariant if

    \[
    \delta(gx)=g\delta(x),
    \qquad g\in G.
    \]

### 6.2 Invariant loss and risk

A loss function is invariant if

\[
L(g\theta,ga)=L(\theta,a),
\qquad g\in G.
\]

When both the model and the loss are invariant, equivariant estimators often have risk functions that are constant over group orbits, greatly simplifying the decision problem.

!!! example "Example 6.2 (Normal location model)"

    If

    \[
    X\sim N(\theta,\sigma^2),
    \qquad
    L(\theta,a)=(a-\theta)^2,
    \]

    then the equivariant estimator $\delta(X)=X$ has risk

    \[
    R(\theta,X)=\sigma^2,
    \]

    which is constant in $\theta$.

!!! warning "Equivariance is not admissibility"

    A rule may be optimal within the equivariant class and still be inadmissible among all decision rules.

    The James–Stein phenomenon is the canonical example.

---

## 7. The Multivariate Normal Mean Problem

Consider

\[
X\sim N_p(\theta,\sigma^2I_p),
\qquad
\theta\in\mathbb R^p,
\]

where $\sigma^2$ is known, with total squared-error loss

\[
L(\theta,a)=\|a-\theta\|^2.
\]

The usual estimator is

\[
\delta_0(X)=X.
\]

It has many strong optimality properties:

1. Each coordinate $X_j$ is UMVU for $\theta_j$;
2. it has minimum total risk among unbiased vector estimators;
3. it is the MLE;
4. it is minimax;
5. it is translation equivariant and best equivariant.

Its risk is

\[
R(\theta,X)
=E_\theta\|X-\theta\|^2
=p\sigma^2.
\]

!!! note "The surprise"

    Even though $X$ is simultaneously coordinatewise UMVU, MLE, best equivariant, and minimax, it is still inadmissible when $p\geq3$.

---

## 8. The James–Stein Estimator

For $p\geq3$, define

\[
\delta_c(X)
=
\left(
1-
\frac{c\sigma^2(p-2)}{\|X\|^2}
\right)X,
\qquad 0<c<2.
\]

This estimator shrinks $X$ toward the origin. Its shrinkage factor is

\[
1-\frac{c\sigma^2(p-2)}{\|X\|^2}.
\]

When $\|X\|^2$ is large, shrinkage is mild; when $\|X\|^2$ is small, shrinkage is strong. Thus the procedure is data-adaptive rather than a fixed contraction.

!!! success "Theorem 8.1 (James–Stein phenomenon)"

    Let

    \[
    X\sim N_p(\theta,\sigma^2I_p),
    \qquad p\geq3.
    \]

    Under squared-error loss, every

    \[
    0<c<2
    \]

    satisfies

    \[
    R(\theta,\delta_c)
    <R(\theta,X)
    =p\sigma^2,
    \qquad \forall\theta\in\mathbb R^p.
    \]

    The canonical choice $c=1$ gives

    \[
    \boxed{
    \delta_{\mathrm{JS}}(X)
    =
    \left(
    1-\frac{\sigma^2(p-2)}{\|X\|^2}
    \right)X
    }.
    \]

    Hence the usual estimator $X$ is inadmissible for $p\geq3$.

!!! warning "Dimension matters"

    For $p=1$ and $p=2$, the usual estimator $X$ is admissible.

    The James–Stein inadmissibility phenomenon begins only at

    \[
    p\geq3.
    \]

---

## 9. Risk Calculation Through Stein's Identity and SURE

The James–Stein result can be understood elegantly through Stein's identity.

Write an estimator as

\[
\delta(X)=X+g(X),
\]

where

\[
g:\mathbb R^p\to\mathbb R^p
\]

is weakly differentiable and satisfies suitable integrability conditions.

!!! success "Theorem 9.1 (Stein's identity)"

    For

    \[
    X\sim N_p(\theta,\sigma^2I_p),
    \]

    \[
    E_\theta\{(X-\theta)^Tg(X)\}
    =\sigma^2E_\theta\{\nabla\cdot g(X)\},
    \]

    where

    \[
    \nabla\cdot g(x)
    =\sum_{j=1}^p
    \frac{\partial g_j(x)}{\partial x_j}
    \]

    is the divergence.

Therefore

\[
\begin{aligned}
R(\theta,X+g)
&=E_\theta\|X+g(X)-\theta\|^2\\
&=E_\theta\|X-\theta\|^2
  +E_\theta\|g(X)\|^2
  +2E_\theta\{(X-\theta)^Tg(X)\}\\
&=p\sigma^2
  +E_\theta\left[
  \|g(X)\|^2
  +2\sigma^2\nabla\cdot g(X)
  \right].
\end{aligned}
\]

Thus

\[
\boxed{
\operatorname{SURE}(X+g)
=p\sigma^2
+\|g(X)\|^2
+2\sigma^2\nabla\cdot g(X)
}
\]

is an unbiased estimator of the risk. This is Stein's Unbiased Risk Estimate (SURE).

### 9.1 Applying SURE to the James–Stein form

Take

\[
g(x)=-\frac{a}{\|x\|^2}x,
\qquad a>0.
\]

Then

\[
\|g(x)\|^2
=\frac{a^2}{\|x\|^2}.
\]

Also,

\[
\nabla\cdot g(x)
=-a\nabla\cdot\left(\frac{x}{\|x\|^2}\right).
\]

A direct calculation gives

\[
\nabla\cdot\left(\frac{x}{\|x\|^2}\right)
=\frac{p-2}{\|x\|^2}.
\]

Hence

\[
R\left(\theta,X-\frac{aX}{\|X\|^2}\right)
-p\sigma^2
=
E_\theta\left[
\frac{a^2-2a\sigma^2(p-2)}{\|X\|^2}
\right].
\]

The quadratic term is minimized at

\[
a^*=\sigma^2(p-2).
\]

Writing

\[
a=c\sigma^2(p-2)
\]

gives

\[
\boxed{
R(\theta,\delta_c)-R(\theta,X)
=-c(2-c)\sigma^4(p-2)^2
E_\theta\left(\frac{1}{\|X\|^2}\right)
}.
\]

For $0<c<2$, the right-hand side is strictly negative, proving the James–Stein improvement.

!!! note "Why is the critical dimension equal to 3?"

    The key term is

    \[
    p-2
    \]

    in the divergence calculation. It becomes positive precisely when $p\geq3$.

---

## 10. Positive-Part James–Stein

The original James–Stein shrinkage factor

\[
1-\frac{\sigma^2(p-2)}{\|X\|^2}
\]

can become negative. When

\[
\|X\|^2<\sigma^2(p-2),
\]

the estimator reverses the direction of $X$, which is intuitively undesirable.

The positive-part James–Stein estimator is therefore defined as

\[
\boxed{
\delta_{\mathrm{JS}+}(X)
=
\left(
1-\frac{\sigma^2(p-2)}{\|X\|^2}
\right)_+X
}
\]

where

\[
u_+=\max(u,0).
\]

!!! success "Important fact"

    The positive-part James–Stein estimator dominates the original James–Stein estimator.

    Hence the original James–Stein estimator is itself inadmissible.

The positive-part rule is nevertheless not fully admissible either. Thus finding one improvement does not mean that the admissibility problem has been solved completely.

---

## 11. Interpretation and Shrinkage Targets

The James–Stein phenomenon is not magic. When many noisy coordinates are estimated jointly, the collection of coordinates can be used to infer how much signal is present and therefore how strongly one should shrink.

### 11.1 Hierarchical and empirical-Bayes interpretation

Consider the hierarchical model

\[
Y_i\mid\theta_i
\sim N(\theta_i,\sigma_0^2),
\]

with

\[
\theta_i\sim N(\mu,\tau^2).
\]

The posterior mean of $\theta_i$ is a weighted average of $Y_i$ and $\mu$, naturally producing shrinkage toward $\mu$.

Estimating the unknown amount of shrinkage from the ensemble leads to empirical-Bayes rules of James–Stein type.

### 11.2 The shrinkage center need not be zero

Zero is only one convenient shrinkage center. For a fixed scientifically meaningful

\[
m\in\mathbb R^p,
\]

define

\[
\delta_m(X)
=m+
\left(
1-\frac{\sigma^2(p-2)}{\|X-m\|^2}
\right)(X-m).
\]

This shrinks $X$ toward $m$.

!!! warning "A random shrinkage target requires a new calculation"

    If $m$ is estimated from the data, for example by the grand mean, then the effective dimension and the risk calculation change.

    One cannot simply substitute a random target into the fixed-target formula without rechecking the degrees of freedom and the risk derivation.

---

## 12. Affine Exponential-Family Rules

Consider the natural exponential family

\[
dP_\eta(x)
=C(\eta)\exp\{\eta T(x)\}\,d\mu(x),
\qquad
\eta\in(\underline\eta,\overline\eta),
\]

and estimate

\[
g(\eta)=E_\eta T
\]

under squared-error loss.

Consider the affine shrinkage rule

\[
\delta_{k,\lambda}(T)
=\frac{T+k\lambda}{1+\lambda},
\qquad \lambda\geq0.
\]

It shrinks $T$ toward $k$.

!!! success "Theorem 12.1 (Karlin's admissibility criterion)"

    Suppose there exists an interior point $\eta_0$ such that

    \[
    \int_{\underline\eta}^{\eta_0}
    C(\eta)^{-\lambda}e^{-k\lambda\eta}\,d\eta
    =\infty
    \]

    and

    \[
    \int_{\eta_0}^{\overline\eta}
    C(\eta)^{-\lambda}e^{-k\lambda\eta}\,d\eta
    =\infty.
    \]

    Then $\delta_{k,\lambda}$ is admissible for estimating $g(\eta)$.

For $\lambda=0$, this includes the unshrunk rule $T$ in standard full exponential families.

!!! note "Course-level message"

    The long proof is not required here. The main lesson is that admissible affine shrinkage rules in exponential families are closely tied to conjugate-Bayes or generalized-Bayes structure.

---

## 13. Why the James–Stein Result Is Surprising

There are four especially important aspects of the phenomenon.

### 13.1 The coordinates need not be scientifically related

Even if the $p$ means represent unrelated scientific quantities, joint estimation under total squared error

\[
\sum_{j=1}^p(\hat\theta_j-\theta_j)^2
\]

allows information sharing across coordinates.

### 13.2 Each coordinate may become biased while total risk improves

The James–Stein estimator tolerates a small bias in individual coordinates in exchange for a sufficiently large aggregate variance reduction.

### 13.3 An excellent estimator can still be inadmissible

The usual estimator $X$ is simultaneously

- unbiased;
- the MLE;
- equivariant;
- minimax;
- coordinatewise UMVU.

Yet it is dominated when $p\geq3$.

### 13.4 The dimension threshold is genuine

The improvement comes directly from the factor

\[
p-2
\]

in the divergence calculation, so the critical dimension is exactly $p=3$.

---

## 14. Comparison of Finite-Sample Optimality Criteria

| Criterion | Meaning | Caveat |
|---|---|---|
| Unbiasedness | Correct on average for every parameter value | May have unnecessarily large variance |
| UMVU | Minimum variance among unbiased estimators | May be inadmissible among all estimators |
| Admissibility | No estimator uniformly dominates it | Weak criterion; many admissible rules may be incomparable |
| Minimaxity | Minimizes worst-case risk | Can be conservative and does not imply admissibility |
| Equivariance | Respects the symmetries of the problem | Best equivariant may be inadmissible outside the equivariant class |
| Bayes optimality | Minimizes average risk under a prior | Depends on the prior; improper priors require care |

---

## 15. Chapter Summary

The main conclusions of this chapter are:

1. Decision theory compares procedures through the risk function $R(\theta,\delta)$ rather than only through the variance of unbiased estimators.
2. An estimator is inadmissible if another estimator has no larger risk for every parameter value and strictly smaller risk somewhere.
3. Under a convex loss, Rao–Blackwell conditioning on a sufficient statistic reduces risk; under strict convexity, it often gives a strict improvement when the original rule is not already a function of the sufficient statistic.
4. UMVU does not imply global risk optimality. In the normal variance problem, $S^2$ is UMVU but is dominated by a suitable biased shrinkage estimator.
5. Shrinkage trades a controlled amount of bias for variance reduction. Fixed shrinkage cannot uniformly dominate the usual estimator in the one-dimensional normal mean problem.
6. Equivariance exploits symmetry and simplifies risk calculations, but the best equivariant rule may still be inadmissible in the full decision class.
7. For $X\sim N_p(\theta,\sigma^2I_p)$ with $p\geq3$, the James–Stein estimator uniformly dominates the usual estimator $X$.
8. Stein's identity and SURE give a concise derivation of the James–Stein risk improvement; the key quantity is the factor $p-2$ in the divergence.
9. The positive-part James–Stein estimator further dominates the original James–Stein estimator, showing that one risk improvement does not automatically establish admissibility.
10. Shrinkage is deeply connected with hierarchical Bayes, empirical Bayes, and generalized-Bayes structures in exponential families.
