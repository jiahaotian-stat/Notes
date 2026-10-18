# Chapter 7: Admissibility, MRU, and Minimaxity

This chapter studies four optimality criteria in statistical decision theory: Bayes optimality, admissibility, minimum risk unbiasedness (MRU), and minimaxity. We then present three routes from Bayes optimality to admissibility, Blyth's method, and least favorable priors, and use Bernoulli and normal mean problems to distinguish these criteria.

---

## 1. Four Optimality Criteria

Let the frequentist risk be

\[
R(\theta,\delta)
=E_\theta L\{\theta,\delta(X)\}.
\]

For a prior probability distribution $\pi$, define the Bayes risk

\[
r(\pi,\delta)
=\int_\Theta R(\theta,\delta)\,\pi(d\theta),
\qquad
B(\pi)=\inf_\delta r(\pi,\delta).
\]

For a nonzero, possibly improper, measure $\nu$, the same integral is denoted by $r(\nu,\delta)$. A rule that minimizes this integral among rules with finite integrated risk is called a generalized Bayes rule for $\nu$.

| Criterion | Definition | What It Protects Against |
| --- | --- | --- |
| Bayes | Minimizes prior-average risk | Average performance under a chosen prior |
| Admissible | Is not uniformly dominated | Avoidable improvement at every parameter |
| MRU | Minimizes risk among unbiased rules | Variance within the unbiased class |
| Minimax | Minimizes $\sup_\theta R(\theta,\delta)$ | The worst-case parameter value |

These criteria are related, but none should be used as a synonym for another. In particular, MRU is defined relative to the restricted class of unbiased rules, whereas admissibility and minimaxity are normally defined by comparison with all decision rules.

---

## 2. When a Bayes Rule Is Admissible

We say that $\delta_1$ dominates $\delta_0$ if

\[
R(\theta,\delta_1)
\leq R(\theta,\delta_0)
\qquad
\text{for every }\theta
\]

and the inequality is strict at least at one parameter value. A rule is admissible if no other rule dominates it.

!!! success "Theorem 2.1 (Three Routes from Bayes Optimality to Admissibility)"

    Let $\delta_\nu$ minimize $r(\nu,\delta)$ for a nonzero measure $\nu$. Any one of the following sets of conditions is sufficient for the admissibility of $\delta_\nu$.

    **(i) Continuity and full support:** Every risk function is continuous; every nonempty open subset of $\Theta$ has positive $\nu$-measure; and at least one rule has finite integrated risk.

    **(ii) Strong uniqueness:** $r(\nu,\delta_\nu)<\infty$, and $\delta_\nu$ is the unique Bayes rule in the following decision-theoretic sense: every other minimizing rule $\delta$ satisfies

    \[
    P_\theta\{\delta(X)=\delta_\nu(X)\}=1
    \qquad
    \text{for every }\theta.
    \]

    **(iii) Countable parameter space:** $\Theta$ is countable,

    \[
    0<\nu\{\theta\}<\infty
    \qquad
    \text{for every }\theta,
    \]

    and $r(\nu,\delta_\nu)<\infty$.

    In condition (i), $\nu$ may be an improper $\sigma$-finite measure, provided that the integrated risks are well defined.

??? proof "Proof of Theorem 2.1 (click to expand)"

    **Case (i).** Optimality and the existence of a finite-risk competitor imply

    \[
    r(\nu,\delta_\nu)<\infty.
    \]

    Suppose $\delta$ dominates $\delta_\nu$, with strict improvement at $\theta_0$. The continuous function

    \[
    g(\theta)
    =R(\theta,\delta_\nu)-R(\theta,\delta)
    \]

    is nonnegative and positive at $\theta_0$. Hence, for some $\epsilon>0$ and some nonempty open neighborhood $U$ of $\theta_0$,

    \[
    g(\theta)\geq\epsilon.
    \]

    Full support gives $\nu(U)>0$, so

    \[
    r(\nu,\delta_\nu)-r(\nu,\delta)
    =\int_\Theta g(\theta)\,\nu(d\theta)>0,
    \]

    contradicting the optimality of $\delta_\nu$.

    **Case (ii).** Dominance gives

    \[
    r(\nu,\delta)\leq r(\nu,\delta_\nu),
    \]

    while optimality gives the reverse inequality. Thus $\delta$ is also a Bayes rule. Strong uniqueness says that the two rules agree $P_\theta$-almost surely for every $\theta$, making strict risk improvement impossible.

    **Case (iii).** If strict improvement occurs at $\theta_0$, then the term corresponding to $\theta_0$ in the countable sum

    \[
    r(\nu,\delta_\nu)-r(\nu,\delta)
    =\sum_{\theta\in\Theta}
    \{R(\theta,\delta_\nu)-R(\theta,\delta)\}\nu\{\theta\}
    \]

    is positive and all other terms are nonnegative. This again contradicts Bayes optimality. $\square$

!!! note "A Point Prior Is Not Automatically Dangerous"

    Full support is only a sufficient condition, not a necessary one. A Bayes rule under a point-mass prior may still be admissible if the strong uniqueness condition in Theorem 2.1(ii) can be verified.

---

## 3. Bernoulli Sampling and an Improper Prior

Let

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}
\operatorname{Bernoulli}(\theta),
\qquad
0\leq\theta\leq1,
\qquad
T=\sum_{i=1}^nX_i,
\]

and use squared-error loss. Consider the improper measure

\[
\nu(d\theta)
=\frac{d\theta}{\theta(1-\theta)}.
\tag{1}
\]

For $1\leq t\leq n-1$, the formal posterior is

\[
\theta\mid T=t
\sim\operatorname{Beta}(t,n-t),
\qquad
E(\theta\mid T=t)=\frac{t}{n}.
\tag{2}
\]

At $t=0$ and $t=n$, however, the formal posterior is improper. Therefore, the statement "the posterior mean is $T/n$" is incomplete unless the two endpoints are handled separately.

!!! success "Proposition 3.1 (A Rigorous Generalized Bayes Calculation)"

    Under the prior measure (1), the unique rule with minimum finite integrated risk is

    \[
    \delta_0(T)=\frac{T}{n}.
    \]

    This rule is admissible and is also the MRU estimator under squared-error loss.

??? proof "Proof of Proposition 3.1 (click to expand)"

    For a rule $a_t=\delta(t)$, Tonelli's theorem gives

    \[
    r(\nu,\delta)
    =\sum_{t=0}^n
    \binom{n}{t}
    \int_0^1
    (a_t-\theta)^2
    \theta^{t-1}(1-\theta)^{n-t-1}\,d\theta.
    \tag{3}
    \]

    For $1\leq t\leq n-1$, the $t$th term is uniquely minimized by $t/n$, the mean of the finite Beta measure.

    When $t=0$, the integral in (3) diverges near zero unless $a_0=0$; it is finite when $a_0=0$. Similarly, the $t=n$ term is finite if and only if $a_n=1$. Hence $a_t=t/n$ is the unique finite-risk minimizer, including the two endpoints.

    Its frequentist risk and integrated risk are

    \[
    R(\theta,T/n)
    =\frac{\theta(1-\theta)}{n},
    \qquad
    r(\nu,T/n)
    =\int_0^1\frac{1}{n}\,d\theta
    =\frac{1}{n}.
    \]

    The risk is continuous, and the measure in (1) assigns positive mass to every open subinterval. Theorem 2.1(i) therefore proves that $T/n$ is admissible.

    Finally, $T/n$ is unbiased, and $T$ is complete and sufficient for the binomial family. By the Lehmann-Scheffe theorem, $T/n$ is UMVU; equivalently, it is MRU under squared-error loss. $\square$

!!! warning "Lesson"

    An improper prior is not a probability distribution, and a formal posterior need not exist at every sample point. When exceptional samples occur, generalized Bayes claims should be checked through integrated risk.

---

## 4. Normal Mean: The Affine Rules $c\overline{X}+d$

Suppose

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}
N(\theta,\sigma^2),
\qquad
v=\frac{\sigma^2}{n},
\qquad
Y=\overline{X}\sim N(\theta,v),
\]

where $\sigma^2$ is known. For

\[
\delta_{c,d}(Y)=cY+d,
\]

the risk is

\[
R(\theta,\delta_{c,d})
=c^2v+\{(c-1)\theta+d\}^2.
\tag{4}
\]

!!! success "Proposition 4.1 (Complete Classification of Affine Rules)"

    Among the affine rules, and in fact against arbitrary competitors, the following rules are admissible:

    \[
    0\leq c<1,
    \qquad
    d\in\mathbb{R},
    \qquad
    \text{and }(c,d)=(1,0).
    \]

    Every other affine rule is inadmissible.

??? proof "Proof of Proposition 4.1 (click to expand)"

    **When $0<c<1$.** Choose the prior

    \[
    \theta\sim N(\mu,\tau^2),
    \qquad
    \tau^2=\frac{cv}{1-c},
    \qquad
    \mu=\frac{d}{1-c}.
    \]

    The posterior mean is

    \[
    \frac{\tau^2}{\tau^2+v}Y
    +\frac{v}{\tau^2+v}\mu
    =cY+d.
    \]

    This is the unique Bayes rule for a proper prior with a positive density on the whole real line, and its risk is continuous and finite. It is therefore admissible by Theorem 2.1.

    **When $c=0$.** The constant rule $d$ is Bayes under the point-mass prior at $\theta=d$. Any other Bayes rule has zero risk at $\theta=d$, so it equals $d$ $P_d$-almost surely. All scalar normal distributions have positive densities and are mutually absolutely continuous, so the rule equals $d$ $P_\theta$-almost surely for every $\theta$. Strong uniqueness proves admissibility.

    **When $c>1$.** By (4),

    \[
    R(\theta,cY+d)
    \geq c^2v
    >v
    =R(\theta,Y),
    \]

    so $Y$ dominates $cY+d$.

    **When $c<0$.** Let

    \[
    m=\frac{d}{1-c}.
    \]

    Then

    \[
    R(\theta,cY+d)
    =c^2v+(1-c)^2(\theta-m)^2
    >(\theta-m)^2
    =R(\theta,m),
    \]

    so the constant rule $m$ dominates $cY+d$.

    **When $c=1$ and $d\neq0$.**

    \[
    R(\theta,Y+d)
    =v+d^2
    >v,
    \]

    so $Y$ dominates $Y+d$.

    It remains to justify the case $(c,d)=(1,0)$. The rule $Y$ is generalized Bayes for Lebesgue measure, but its integrated risk is infinite, so Theorem 2.1(i) cannot be applied directly. Its admissibility will be proved by Blyth's method in Section 5. $\square$

The statistic $Y=\overline{X}$ is also MRU: it is unbiased and complete sufficient for $\theta$, so the Lehmann-Scheffe theorem applies. This fact alone does not yet prove either admissibility or minimaxity.

---

## 5. Blyth's Method: Admissibility from a Prior Sequence

Blyth's method is useful when the target rule is generalized Bayes for an improper measure but ordinary Bayes-risk arguments fail.

!!! success "Theorem 5.1 (Blyth's Method)"

    Let $\mathcal{D}_c$ be a complete class of rules with continuous risk functions; that is, every rule outside $\mathcal{D}_c$ is weakly dominated by a rule in $\mathcal{D}_c$. Let $\delta_0\in\mathcal{D}_c$.

    Suppose that, for every $j$, there exists a nonzero measure $\nu_j$ such that $r(\nu_j,\delta_0)<\infty$, and for every $\theta_0$ and $\eta>0$,

    \[
    0<\nu_j\{B(\theta_0,\eta)\}<\infty,
    \]

    and

    \[
    \frac{
    r(\nu_j,\delta_0)
    -\inf_{\delta\in\mathcal{D}_c}r(\nu_j,\delta)
    }{
    \nu_j\{B(\theta_0,\eta)\}
    }
    \longrightarrow0.
    \tag{5}
    \]

    Then $\delta_0$ is admissible.

??? proof "Proof of Theorem 5.1 (click to expand)"

    Suppose $\delta_1\in\mathcal{D}_c$ dominates $\delta_0$, with strict improvement at $\theta_0$. By continuity, there exist $\eta>0$ and $\epsilon>0$ such that

    \[
    R(\theta,\delta_0)-R(\theta,\delta_1)
    \geq\epsilon
    \qquad
\text{on the set }B(\theta_0,\eta)\text{ itself}.
    \]

    Therefore,

    \[
    \begin{aligned}
    r(\nu_j,\delta_0)
    -\inf_{\delta\in\mathcal{D}_c}r(\nu_j,\delta)
    &\geq
    r(\nu_j,\delta_0)-r(\nu_j,\delta_1)\\
    &\geq
    \epsilon\,\nu_j\{B(\theta_0,\eta)\}.
    \end{aligned}
    \]

    Thus the ratio in (5) is at least $\epsilon$ for every $j$, a contradiction. Completeness of $\mathcal{D}_c$ rules out a dominator outside the class. $\square$

### 5.1 Application: Admissibility of the Scalar Normal Mean

Return to

\[
Y\sim N(\theta,v),
\qquad
\delta_0(Y)=Y.
\]

Let

\[
\nu_j=N(0,j).
\]

The corresponding Bayes rule and Bayes risk are

\[
\delta_j(Y)=\frac{j}{j+v}Y,
\qquad
B(\nu_j)=\frac{vj}{v+j}.
\]

Since $r(\nu_j,Y)=v$, the Bayes-risk gap is

\[
r(\nu_j,Y)-B(\nu_j)
=\frac{v^2}{v+j}
=O(j^{-1}).
\tag{6}
\]

For the fixed interval

\[
B(\theta_0,\eta)
=(\theta_0-\eta,\theta_0+\eta),
\]

we have

\[
\begin{aligned}
\nu_j\{B(\theta_0,\eta)\}
&=
\Phi\left(\frac{\theta_0+\eta}{\sqrt{j}}\right)
-\Phi\left(\frac{\theta_0-\eta}{\sqrt{j}}\right)\\
&\sim\frac{2\eta\phi(0)}{\sqrt{j}}.
\end{aligned}
\tag{7}
\]

Dividing (6) by (7) gives

\[
O(j^{-1/2})\longrightarrow0.
\]

Under the standard complete-class regularity conditions for this scalar normal problem, Theorem 5.1 proves that

\[
Y=\overline{X}
\]

is admissible.

!!! note "The Core of Blyth's Argument"

    The essential step is the rate comparison: the Bayes-risk gap vanishes as $j^{-1}$, faster than the prior mass $j^{-1/2}$ of any fixed neighborhood. Merely saying that the proper priors "approach the flat prior" is not enough.

---

## 6. Minimaxity and Least Favorable Priors

Define

\[
M(\delta)
=\sup_{\theta\in\Theta}R(\theta,\delta),
\qquad
V=\inf_\delta M(\delta).
\]

For every proper prior $\pi$ and every rule $\delta$,

\[
B(\pi)
\leq r(\pi,\delta)
\leq M(\delta).
\tag{8}
\]

Thus every Bayes risk is a lower bound for the minimax risk.

!!! success "Theorem 6.1 (Two Bayes Proofs of Minimaxity)"

    **(i)** If $\delta_\pi$ is Bayes for a proper prior $\pi$ and

    \[
    B(\pi)
    =r(\pi,\delta_\pi)
    =M(\delta_\pi)
    <\infty,
    \tag{9}
    \]

    then $\delta_\pi$ is minimax and $\pi$ is least favorable. In particular, a proper Bayes rule with constant risk is minimax.

    **(ii)** If $\{\pi_j\}$ is a sequence of proper priors and

    \[
    B(\pi_j)
    \longrightarrow M(\delta_0)
    <\infty,
    \tag{10}
    \]

    then $\delta_0$ is minimax and $\{\pi_j\}$ is a least favorable sequence.

??? proof "Proof of Theorem 6.1 (click to expand)"

    Under (9), for every $\delta$, (8) gives

    \[
    M(\delta)
    \geq B(\pi)
    =M(\delta_\pi),
    \]

    which proves minimaxity.

    Moreover, for every prior $\pi'$,

    \[
    B(\pi')
    \leq V
    =M(\delta_\pi)
    =B(\pi),
    \]

    so $\pi$ is least favorable.

    Under (10), (8) gives, for every $j$ and every $\delta$,

    \[
    B(\pi_j)\leq M(\delta).
    \]

    Letting $j\to\infty$ yields

    \[
    M(\delta_0)\leq M(\delta),
    \]

    proving that $\delta_0$ is minimax and that $\{\pi_j\}$ is a least favorable sequence. $\square$

---

## 7. Binomial Example: MRU Need Not Be Minimax

Let

\[
T\sim\operatorname{Bin}(n,\theta),
\]

and estimate $\theta$ under squared-error loss. The MRU estimator $T/n$ has risk

\[
R(\theta,T/n)
=\frac{\theta(1-\theta)}{n},
\qquad
M(T/n)=\frac{1}{4n}.
\]

It is admissible by Proposition 3.1, but it is not minimax.

Take the proper prior

\[
\theta\sim\operatorname{Beta}(a,a).
\]

Its Bayes rule is

\[
\delta_a(T)
=\frac{T+a}{n+2a}.
\]

For every $a>0$,

\[
R(\theta,\delta_a)
=\frac{
n\theta(1-\theta)+a^2(1-2\theta)^2
}{
(n+2a)^2
}.
\tag{11}
\]

Choose

\[
a=\frac{\sqrt{n}}{2}.
\]

Then $a^2=n/4$, and the numerator of (11) becomes

\[
n\left\{
\theta(1-\theta)
+\frac{1}{4}(1-2\theta)^2
\right\}
=\frac{n}{4}.
\]

Therefore,

\[
R(\theta,\delta_a)
=\frac{1}{4(\sqrt{n}+1)^2}
\qquad
\text{for every }\theta.
\]

This rule is proper Bayes and has constant risk. Theorem 6.1(i) therefore proves that it is minimax. Moreover,

\[
\frac{1}{4(\sqrt{n}+1)^2}
<\frac{1}{4n}
=M(T/n).
\]

!!! info "Interpretation"

    The MRU rule refuses bias and has zero risk at $\theta=0,1$, but pays a larger risk near $\theta=1/2$. The minimax rule accepts bias near the boundary to flatten and reduce the worst-case risk. Both rules are admissible.

---

## 8. Normal Example: The Same Priors, a Different Conclusion

For

\[
Y\sim N(\theta,v),
\]

the rule $Y$ has constant risk $v$. Take

\[
\pi_j=N(0,j).
\]

Then

\[
B(\pi_j)
=\frac{vj}{v+j}
\longrightarrow v
=M(Y).
\]

Theorem 6.1(ii) therefore proves that $Y$ is minimax.

The same prior sequence serves two different purposes:

- For minimaxity, compare the absolute Bayes risk $B(\pi_j)$ with the maximum risk $M(Y)$.

- For admissibility, compare the Bayes-risk gap $r(\pi_j,Y)-B(\pi_j)$ with the local prior mass of every neighborhood.

Minimaxity does not imply admissibility. For a multivariate normal mean of dimension at least three, the usual estimator is minimax but is dominated by James-Stein shrinkage. That result belongs to Chapter 5 and is not repeated here.

---

## 9. Relations Worth Remembering

!!! success "Proposition 9.1 (Constant Risk and Uniqueness)"

    **(i)** An admissible rule with finite constant risk is minimax.

    **(ii)** A unique minimax rule is admissible.

    **(iii)** Under a strictly convex loss and a convex action space, an admissible constant-risk rule is the unique minimax rule, up to almost-sure equality.

??? proof "Proof of Proposition 9.1 (click to expand)"

    Let

    \[
    R(\theta,\delta_0)=C.
    \]

    If $\delta_0$ were not minimax, there would exist $\delta_1$ such that

    \[
    M(\delta_1)<C.
    \]

    Then $\delta_1$ would dominate $\delta_0$, contradicting admissibility. This proves (i).

    If a unique minimax rule were inadmissible, a rule dominating it would have no larger maximum risk and would therefore also be minimax, contradicting uniqueness. This proves (ii).

    For (iii), suppose $\delta_1$ is another minimax rule. By convexity of the action space, the midpoint rule

    \[
    \delta_m
    =\frac{\delta_0+\delta_1}{2}
    \]

    is allowed. Strict convexity of the loss gives

    \[
    R(\theta,\delta_m)
    \leq
    \frac{1}{2}R(\theta,\delta_0)
    +\frac{1}{2}R(\theta,\delta_1)
    \leq C,
    \]

    with strict inequality wherever the two rules differ with positive $P_\theta$-probability. If they differ for some $\theta$, then $\delta_m$ dominates $\delta_0$, contradicting admissibility. Hence the minimax rule is unique up to almost-sure equality. $\square$

---

## 10. Chapter Summary

!!! success "Key Conclusions"

    1. Bayes optimality, admissibility, MRU, and minimaxity are distinct notions of optimality with different comparison classes and risk objectives.

    2. Continuity and full support, strong uniqueness, or a countable parameter space can each promote Bayes optimality to admissibility.

    3. Improper priors require rigorous treatment through integrated risk, especially at sample points where the formal posterior fails.

    4. Blyth's method proves admissibility of generalized Bayes rules by comparing the convergence rate of the Bayes-risk gap with local prior mass.

    5. Bayes risks provide lower bounds for minimax risk; constant-risk Bayes rules and least favorable prior sequences can prove minimaxity.

    6. An MRU rule need not be minimax, and a minimax rule need not be admissible; constant risk, uniqueness, and convexity can establish stronger relations.
