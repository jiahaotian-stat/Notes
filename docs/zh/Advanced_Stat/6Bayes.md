# 第六章：Bayes 规则、广义 Bayes 与经验 Bayes

本章从频率学派的风险函数进一步走向 Bayes 决策论。前一章讨论可容许性时，我们比较整个风险函数
\[
\theta\longmapsto R(\theta,\delta),
\]
但两个估计量的风险曲线常常会交叉，因此无法给出统一排序。Bayes 方法通过在参数空间上引入先验分布，将风险函数加权平均为一个标量的 Bayes 风险，并把全局优化问题转化为逐点最小化后验期望损失。

本章将依次介绍 Bayes 风险与后验风险、常见损失下的 Bayes 估计量、共轭先验、广义 Bayes、Jeffreys 先验、MAP 与正则化，以及经验 Bayes 的基本思想。

---

## 1. 从频率学派风险到 Bayes 风险

### 1.1 Rule、estimator 与 estimate

设样本空间为 $\mathcal X$，行动空间为 $\mathcal A$。一个**决策规则**（decision rule）是可测映射

\[
\delta:\mathcal X\to\mathcal A.
\]

当行动 $a$ 的目的是估计某个目标 $g(\theta)$ 时，同一个对象也称为**估计量**（estimator）。

在观测到具体数据 $X=x$ 后，

\[
\delta(x)
\]

是一个实际的数值或向量，称为**估计值**（estimate）。

!!! note "术语区别"

    - **Bayes rule**：一般决策问题中的 Bayes 最优决策规则；
    - **Bayes estimator**：当行动是估计某个参数或函数 $g(\theta)$ 时的 Bayes rule；
    - **Bayes estimate**：给定具体观测 $x$ 后得到的数值 $\delta_\pi(x)$。

### 1.2 频率学派风险

给定损失函数 $L(\theta,a)$，频率学派风险为

\[
R(\theta,\delta)
=
E_\theta L\{\theta,\delta(X)\},
\]

其中 $\theta$ 被视为固定未知量，期望仅关于重复抽样

\[
X\sim P_\theta
\]

计算。

因此，频率学派比较的是整个函数

\[
\theta\mapsto R(\theta,\delta).
\]

### 1.3 引入先验分布

Bayes 决策论在参数空间 $\Theta$ 上引入先验分布 $\pi$：

\[
\theta\sim\pi,
\qquad
X\mid\theta\sim P_\theta.
\]

于是 $(\theta,X)$ 的联合分布可以写为

\[
H(d\theta,dx)
=
\pi(d\theta)P_\theta(dx).
\]

$X$ 的边缘分布为

\[
m_\pi(dx)
=
\int_\Theta P_\theta(dx)\,\pi(d\theta).
\]

若存在正则条件分布，则后验分布记为

\[
\pi(d\theta\mid x).
\]

同一个联合分布可以分解为

\[
\pi(d\theta)P_\theta(dx)
=
m_\pi(dx)\pi(d\theta\mid x).
\]

若 $P_\theta$ 相对于测度 $\mu$ 有密度 $p(x\mid\theta)$，先验相对于测度 $\nu$ 有密度 $\pi(\theta)$，则

\[
m_\pi(x)
=
\int_\Theta
p(x\mid\theta)\pi(\theta)\,\nu(d\theta),
\]

并且

\[
\pi(\theta\mid x)
=
\frac{p(x\mid\theta)\pi(\theta)}
{m_\pi(x)}.
\]

由于归一化常数 $m_\pi(x)$ 与行动 $a$ 无关，因此在寻找 Bayes rule 时，可以直接最小化

\[
\int_\Theta
L(\theta,a)p(x\mid\theta)\pi(\theta)\,\nu(d\theta).
\]

---

## 2. Bayes 风险与后验风险

!!! info "定义 2.1（Bayes 风险）"

    对给定先验 $\pi$，决策规则 $\delta$ 的 **Bayes 风险**（Bayes risk）或 **积分风险**（integrated risk）定义为

    \[
    r(\pi,\delta)
    =
    \int_\Theta
    R(\theta,\delta)\,\pi(d\theta).
    \]

    等价地，

    \[
    r(\pi,\delta)
    =
    \int_\Theta
    \int_{\mathcal X}
    L\{\theta,\delta(x)\}
    P_\theta(dx)\pi(d\theta).
    \]

    若 $\delta_\pi$ 在所有决策规则中最小化 $r(\pi,\delta)$，则称 $\delta_\pi$ 为关于先验 $\pi$ 的 **Bayes rule**。

定义给定 $X=x$ 时的**后验期望损失**

\[
C(a,x)
=
\int_\Theta
L(\theta,a)\pi(d\theta\mid x).
\]

利用联合分布的两种分解，

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

这说明 Bayes 风险的最小化可以逐个 $x$ 进行。

!!! success "定理 2.2（后验风险刻画）"

    假设至少存在一个决策规则具有有限 Bayes 风险。若对 $m_\pi$-几乎处处的 $x$，

    \[
    C\{\delta_\pi(x),x\}
    =
    \inf_{a\in\mathcal A}C(a,x),
    \]

    则 $\delta_\pi$ 是 Bayes rule。

??? proof "定理 2.2 的证明思路（点击展开）"

    因为

    \[
    r(\pi,\delta)
    =
    \int_{\mathcal X}
    C\{\delta(x),x\}\,m_\pi(dx),
    \]

    对每个观测值 $x$ 分别选择使 $C(a,x)$ 最小的行动，就会使整个积分达到最小。

    因此，

    \[
    \boxed{
    \text{最小化 Bayes 风险}
    \iff
    \text{逐点最小化后验期望损失}
    }.
    \]

---

## 3. 常见损失函数下的 Bayes 估计量

设目标为 $g(\theta)$，行动 $a\in\mathbb R$。

### 3.1 平方误差损失

若

\[
L(\theta,a)
=
\{g(\theta)-a\}^2,
\]

则后验风险为

\[
C(a,x)
=
E\left[
\{g(\theta)-a\}^2
\mid X=x
\right].
\]

展开得

\[
C(a,x)
=
E\{g(\theta)^2\mid x\}
-2aE\{g(\theta)\mid x\}
+a^2.
\]

对 $a$ 求导并令其为零，

\[
a
=
E\{g(\theta)\mid X=x\}.
\]

!!! success "结论：平方误差损失"

    在平方误差损失下，Bayes 估计量是**后验均值**：

    \[
    \boxed{
    \delta_\pi(x)
    =
    E\{g(\theta)\mid X=x\}
    }.
    \]

最小后验风险为

\[
\inf_a C(a,x)
=
\operatorname{Var}\{g(\theta)\mid X=x\}.
\]

### 3.2 加权平方误差损失

考虑

\[
L(\theta,a)
=
w(\theta)\{g(\theta)-a\}^2,
\qquad
w(\theta)\ge 0.
\]

令

\[
A_j(x)
=
E\{w(\theta)g(\theta)^j\mid X=x\},
\qquad
j=0,1,2.
\]

则

\[
C(a,x)
=
A_2(x)-2aA_1(x)+a^2A_0(x).
\]

当

\[
0<A_0(x)<\infty
\]

且相应矩有限时，

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

同时，

\[
C(a,x)
=
C\{\delta_\pi(x),x\}
+
A_0(x)\{a-\delta_\pi(x)\}^2.
\]

!!! example "例 3.1（正参数的相对平方误差）"

    假设后验分布采用 rate 参数化：

    \[
    \theta\mid x
    \sim
    \operatorname{Gamma}(k,r),
    \]

    并使用损失

    \[
    L(\theta,a)
    =
    \left(\frac a\theta-1\right)^2
    =
    \theta^{-2}(a-\theta)^2.
    \]

    因此

    \[
    g(\theta)=\theta,
    \qquad
    w(\theta)=\theta^{-2}.
    \]

    当 $k>2$ 时，

    \[
    \delta_\pi(x)
    =
    \frac{E(\theta^{-1}\mid x)}
    {E(\theta^{-2}\mid x)}
    =
    \frac{k-2}{r}.
    \]

    而普通平方误差损失下的 Bayes 估计量为

    \[
    E(\theta\mid x)
    =
    \frac{k}{r}.
    \]

!!! note "重要思想"

    **后验分布本身并没有改变。**

    改变的是损失函数，因此最优行动也会改变。Bayes 决策不仅由后验分布决定，还由我们如何定义“错误的代价”决定。

### 3.3 绝对误差损失

若

\[
L(\theta,a)
=
|g(\theta)-a|,
\]

则任意后验中位数都是 Bayes 估计量，即 $a$ 满足

\[
P\{g(\theta)\le a\mid X=x\}
\ge \frac12,
\]

以及

\[
P\{g(\theta)\ge a\mid X=x\}
\ge \frac12.
\]

因此

\[
\boxed{
\text{absolute loss}
\longrightarrow
\text{posterior median}
}.
\]

### 3.4 0–1 损失与 MAP

若参数是离散的，考虑精确 0–1 损失

\[
L(\theta,a)
=
\mathbf 1\{a\ne\theta\}.
\]

后验风险为

\[
C(a,x)
=
P(\theta\ne a\mid X=x)
=
1-P(\theta=a\mid X=x).
\]

因此 Bayes rule 为后验众数：

\[
\boxed{
\delta_\pi(x)
\in
\arg\max_a
P(\theta=a\mid X=x)
}.
\]

对连续参数，通常

\[
P(\theta=a\mid X=x)=0,
\]

所以精确 0–1 损失并不实用。可以改用区间 0–1 损失

\[
L(\theta,a)
=
\mathbf 1\{|g(\theta)-a|>c\}.
\]

此时应选择长度为 $2c$ 的区间中后验概率最大的区间中心：

\[
\delta_\pi(x)
\in
\arg\max_a
P\{a-c\le g(\theta)\le a+c\mid X=x\}.
\]

当 $c\downarrow0$ 时，在存在后验密度的情况下，这在直觉上导向后验众数，即 MAP。

!!! note "常见损失与 Bayes 点估计"

    \[
    \begin{array}{c|c}
    \text{损失函数} & \text{Bayes 点估计}\\
    \hline
    \text{平方误差} & \text{后验均值}\\
    \text{绝对误差} & \text{后验中位数}\\
    \text{0--1 / 局部 0--1} & \text{后验众数 / MAP}
    \end{array}
    \]

---

## 4. 共轭先验

!!! info "定义 4.1（共轭先验族）"

    对抽样模型 $p(x\mid\theta)$，如果一类先验分布 $\mathcal F$ 满足

    \[
    \pi\in\mathcal F
    \Longrightarrow
    \pi(\cdot\mid x)\in\mathcal F
    \]

    对所有可能的观测 $x$ 都成立，则称 $\mathcal F$ 是该模型的**共轭先验族**。

共轭性的关键优点是**代数闭合性**：先验与后验具有相同的函数形式，只更新超参数。

---

## 5. Beta–Binomial 模型

设

\[
X_1,\ldots,X_n
\overset{\text{iid}}{\sim}
\operatorname{Bernoulli}(p),
\]

并记

\[
T=\sum_{i=1}^nX_i.
\]

似然函数满足

\[
L(p;x)
\propto
p^T(1-p)^{n-T}.
\]

取 Beta 先验

\[
p\sim\operatorname{Beta}(a,b),
\]

其密度为

\[
\pi(p)
=
\frac{\Gamma(a+b)}
{\Gamma(a)\Gamma(b)}
p^{a-1}(1-p)^{b-1},
\qquad
0<p<1.
\]

于是

\[
\pi(p\mid x)
\propto
p^{a+T-1}
(1-p)^{b+n-T-1},
\]

因此

\[
\boxed{
p\mid X
\sim
\operatorname{Beta}(a+T,b+n-T)
}.
\]

平方误差损失下，

\[
\delta_\pi(X)
=
E(p\mid X)
=
\frac{a+T}{a+b+n}.
\]

将其写成收缩形式：

\[
\frac{a+T}{a+b+n}
=
\frac{n}{n+a+b}\bar X
+
\frac{a+b}{n+a+b}
\frac{a}{a+b}.
\]

因此后验均值把样本比例 $\bar X$ 收缩到先验均值

\[
\frac{a}{a+b},
\]

而

\[
a+b
\]

可以理解为先验强度。

!!! example "例 5.1（均匀先验与 Laplace smoothing）"

    若

    \[
    a=b=1,
    \]

    则

    \[
    \delta_\pi(X)
    =
    \frac{T+1}{n+2}.
    \]

    这就是 Laplace smoothing。即使 $T=0$ 或 $T=n$，估计值也不会恰好等于 $0$ 或 $1$。

---

## 6. Gamma–Poisson 模型（拓展）

设

\[
X_1,\ldots,X_n
\overset{\text{iid}}{\sim}
\operatorname{Poisson}(\lambda),
\qquad
T=\sum_{i=1}^nX_i.
\]

似然函数满足

\[
L(\lambda;x)
\propto
e^{-n\lambda}\lambda^T.
\]

采用 rate 参数化的 Gamma 先验：

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

则

\[
\pi(\lambda\mid x)
\propto
\lambda^{\alpha+T-1}
e^{-(\beta+n)\lambda},
\]

因此

\[
\boxed{
\lambda\mid X
\sim
\operatorname{Gamma}(\alpha+T,\beta+n)
}.
\]

在平方误差损失下，

\[
\delta_\pi(X)
=
E(\lambda\mid X)
=
\frac{\alpha+T}{\beta+n}.
\]

又可写为

\[
\frac{\alpha+T}{\beta+n}
=
\frac{n}{\beta+n}
\frac{T}{n}
+
\frac{\beta}{\beta+n}
\frac{\alpha}{\beta}.
\]

因此它是样本均值与先验均值的加权平均，而 $\beta$ 起到类似“先验样本量”的作用。

---

## 7. Gamma–Exponential 模型：三种估计思想的比较

设 $n\ge2$，且

\[
X_1,\ldots,X_n
\overset{\text{iid}}{\sim}
f(x\mid\theta)
=
\theta e^{-\theta x}\mathbf1\{x>0\},
\qquad
\theta>0.
\]

记

\[
T=\sum_{i=1}^nX_i.
\]

似然为

\[
L(\theta;x)
\propto
\theta^n e^{-\theta T}.
\]

取 Gamma 先验

\[
\theta
\sim
\operatorname{Gamma}(\beta,\alpha),
\]

其中密度正比于

\[
\theta^{\beta-1}e^{-\alpha\theta}.
\]

于是

\[
\boxed{
\theta\mid X
\sim
\operatorname{Gamma}(n+\beta,\alpha+T)
}.
\]

我们要估计一个新元件存活超过 $t_0$ 的可靠性概率：

\[
g(\theta)
=
P_\theta(X_{\mathrm{new}}>t_0)
=
e^{-\theta t_0}.
\]

### 7.1 Plug-in 估计量

最大似然估计为

\[
\widehat\theta_{\mathrm{MLE}}
=
\frac nT.
\]

利用 MLE 的不变性，

\[
\boxed{
\widehat g_{\mathrm{plug}}
=
e^{-nt_0/T}
}.
\]

### 7.2 最小风险无偏估计量

由于

\[
\mathbf1\{X_1>t_0\}
\]

满足

\[
E_\theta\mathbf1\{X_1>t_0\}
=
e^{-\theta t_0}
=
g(\theta),
\]

它是 $g(\theta)$ 的无偏估计量。

条件于完全充分统计量 $T$ 后，

\[
\left(
\frac{X_1}{T},
\ldots,
\frac{X_n}{T}
\right)
\]

在单纯形上均匀分布，因此 Rao–Blackwell 化得到

\[
\boxed{
\widehat g_{\mathrm{MRU}}
=
E\{\mathbf1(X_1>t_0)\mid T\}
=
\left(1-\frac{t_0}{T}\right)_+^{n-1}
}.
\]

在凸损失下，这是 minimum-risk unbiased rule；特别是在平方误差损失下，它是 UMVU 估计量。

### 7.3 Bayes 估计量

平方误差损失下，Bayes 估计量是

\[
E(e^{-\theta t_0}\mid X).
\]

利用 Gamma 后验的 Laplace transform，

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

!!! note "三种估计量解决的是不同优化问题"

    - Plug-in estimator：把 $\widehat\theta_{\mathrm{MLE}}$ 代入目标函数；
    - MRU / UMVU estimator：在无偏估计量类中追求最优；
    - Bayes estimator：在给定先验与损失下最小化 Bayes 风险。

    它们没有理由必须相同。

---

## 8. Normal–Normal 模型

设

\[
X_1,\ldots,X_n
\overset{\text{iid}}{\sim}
N(\theta,\sigma^2),
\]

其中 $\sigma^2$ 已知，并取先验

\[
\theta
\sim
N(\mu,\tau^2).
\]

似然只通过 $\bar X$ 依赖数据。配方后得到

\[
\theta\mid X
\sim
N(m_n,v_n),
\]

其中

\[
v_n
=
\frac{1}
{n/\sigma^2+1/\tau^2},
\]

以及

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

等价地，

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

平方误差损失下，

\[
\boxed{
\delta_\pi(X)=m_n
}.
\]

Bayes rule 的后验风险为

\[
\operatorname{Var}(\theta\mid X)
=
v_n.
\]

在这个共轭模型中，$v_n$ 与具体观测值 $X$ 无关。

!!! note "收缩解释"

    后验均值是样本均值与先验均值的**精度加权平均**。

    数据精度为

    \[
    \frac{n}{\sigma^2},
    \]

    先验精度为

    \[
    \frac1{\tau^2}.
    \]

    当 $n$ 增大时，数据权重增加；当 $\tau^2$ 很小时，先验更集中，先验权重增加。

---

## 9. 指数族中的共轭先验（拓展）

考虑规范指数族

\[
p_\eta(t)
=
\exp\{\eta^Tt-A(\eta)\}K(t),
\qquad
\eta\in\mathcal H.
\]

对于 i.i.d. 数据 $T_1,\ldots,T_n$，似然正比于

\[
\exp\left\{
\eta^T\sum_{i=1}^nT_i
-
nA(\eta)
\right\}.
\]

一类共轭先验可以写为

\[
\pi_{n_0,t_0}(\eta)
\propto
\exp\{
n_0t_0^T\eta-n_0A(\eta)
\}.
\]

后验为

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

因此超参数更新为

\[
n_0
\longmapsto
n_0+n,
\]

以及

\[
t_0
\longmapsto
\frac{
n_0t_0+\sum_{i=1}^nT_i
}{
n_0+n
}.
\]

!!! note "Pseudo-observation 解释"

    $n_0$ 可以看作先验的“有效样本量”，而 $t_0$ 可以看作先验中的“平均充分统计量”。

---

## 10. Improper prior 与 generalized Bayes

### 10.1 Improper prior

如果

\[
\int_\Theta\pi(d\theta)=\infty,
\]

则称 $\pi$ 为**非正常先验**（improper prior）。

它不是一个概率分布，但仍可能产生正常后验。若

\[
m(x)
=
\int_\Theta
p(x\mid\theta)\pi(d\theta)
<
\infty,
\]

则形式后验

\[
\pi(d\theta\mid x)
=
\frac{
p(x\mid\theta)\pi(d\theta)
}{
m(x)
}
\]

是一个真正的概率分布。

!!! info "定义 10.1（广义 Bayes rule）"

    如果一个决策规则最小化由 improper prior 所得到的形式后验下的后验期望损失，则称它为关于该 improper prior 的 **generalized Bayes rule**。

### 10.2 正态均值的平坦先验

设

\[
X_1,\ldots,X_n
\overset{\text{iid}}{\sim}
N(\theta,\sigma^2),
\]

其中 $\sigma^2$ 已知。

取 improper flat prior

\[
\pi(\theta)\propto1.
\]

则

\[
\theta\mid X
\sim
N\left(
\bar X,
\frac{\sigma^2}{n}
\right).
\]

因此在平方误差损失下，

\[
\boxed{
\delta(X)=\bar X
}
\]

是 generalized Bayes estimator。

### 10.3 为什么 $\bar X$ 不是 proper Bayes？

!!! success "命题 10.2（$\bar X$ 不是 proper Bayes）"

    在无约束的正态位置模型中，$\bar X$ 不可能是任何 $\mathbb R$ 上 proper prior 所产生的后验均值；但它是 Lebesgue 平坦测度下的 generalized Bayes estimator。

??? proof "命题 10.2 的证明（点击展开）"

    令

    \[
    s^2=\frac{\sigma^2}{n},
    \]

    并设 $\pi$ 是某个 proper prior。$\bar X$ 的边缘密度为

    \[
    m_\pi(x)
    =
    \int
    \phi_s(x-\theta)\,\pi(d\theta),
    \]

    其中 $\phi_s$ 是方差为 $s^2$ 的正态密度。

    对 $x$ 求导可得到 Tweedie identity：

    \[
    E_\pi(\theta\mid\bar X=x)
    =
    x
    +
    s^2
    \frac{m_\pi'(x)}{m_\pi(x)}.
    \]

    如果对所有 $x$ 都有

    \[
    E_\pi(\theta\mid\bar X=x)=x,
    \]

    那么必须有

    \[
    m_\pi'(x)=0
    \]

    对所有 $x$ 成立。

    因此 $m_\pi(x)$ 必须是 $\mathbb R$ 上的正常数。但一个正常数不可能在 $\mathbb R$ 上积分为 $1$。

    所以不存在 proper prior 使后验均值恒等于 $x$。$\square$

!!! warning "Proper Bayes 与 generalized Bayes 不能混淆"

    \[
    \bar X
    \]

    是 generalized Bayes estimator，因为 flat measure 可以产生正常后验；

    但它不是 proper Bayes estimator，因为不存在真正的概率先验使其成为后验均值。

    此外，improper prior 下通常不存在有限的 prior-averaged Bayes risk，因此 generalized Bayes 并不意味着它在某个 proper prior 下最小化有限的 Bayes 风险。

---

## 11. Jeffreys 先验与 reference prior

所谓“non-informative prior”不应理解为真正没有信息。一个在参数 $\theta$ 下看起来平坦的密度，在重新参数化后通常不再平坦。

在一维问题中，一个常用的具有不变性的默认先验是 Jeffreys prior。

!!! info "定义 11.1（Jeffreys prior）"

    若 $I(\theta)$ 是 Fisher information，则

    \[
    \boxed{
    \pi_J(\theta)
    \propto
    \sqrt{I(\theta)}
    }.
    \]

其核心动机是重新参数化不变性。若

\[
\phi=h(\theta)
\]

是一一光滑变换，则 Jeffreys prior 满足

\[
\pi_J(\theta)d\theta
=
\pi_J(\phi)d\phi.
\]

### 11.1 Bernoulli 参数

对 $n$ 个 Bernoulli 观测，

\[
I_n(p)
=
\frac{n}{p(1-p)}.
\]

因此

\[
\pi_J(p)
\propto
\{p(1-p)\}^{-1/2}.
\]

这正是

\[
\operatorname{Beta}\left(\frac12,\frac12\right)
\]

先验。

若

\[
T=\sum_{i=1}^nX_i,
\]

则后验为

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

若

\[
p(x\mid\theta)
=
f(x-\theta),
\]

则 Fisher information 与 $\theta$ 无关，因此

\[
\boxed{
\pi_J(\theta)\propto1
}.
\]

### 11.3 Scale family

若

\[
p(x\mid\sigma)
=
\frac1\sigma
f\left(\frac x\sigma\right),
\]

则 Fisher information 与 $1/\sigma^2$ 成正比，因此

\[
\boxed{
\pi_J(\sigma)\propto\frac1\sigma
}.
\]

它通常是 improper 的，但具有尺度不变性。

### 11.4 多维参数

若

\[
\theta\in\mathbb R^d,
\]

则 Jeffreys prior 写为

\[
\pi_J(\theta)
\propto
\sqrt{\det I(\theta)},
\]

其中 $I(\theta)$ 是 Fisher information matrix。

### 11.5 Reference prior

Reference prior 的思想与 Jeffreys prior 相关：它试图选择一个先验，使得从先验到后验所获得的信息在渐近意义下最大，这种信息通常通过 Kullback–Leibler divergence 描述。

在规则的一参数问题中，reference prior 往往与 Jeffreys prior 一致。

但在多参数问题中，尤其存在 nuisance parameter 时，reference prior 可能依赖于：

- 关注的是哪个参数分量；
- 参数进入构造的顺序。

!!! warning "必要的谨慎"

    不变性本身并不能保证：

    - 后验一定 proper；
    - 有良好的有限样本风险；
    - 在高维 nuisance parameter 模型中一定一致。

    实际应用时仍需检查 posterior propriety、先验敏感性与模型稳健性。

---

## 12. MAP 与正则化（拓展）

考虑 Gaussian 线性模型

\[
y=A\beta+\sigma Z,
\qquad
Z\sim N(0,I_n).
\]

假设先验密度具有形式

\[
\pi(\beta)
\propto
\exp\{-\lambda P(\beta)\}.
\]

后验对数密度满足

\[
\log\pi(\beta\mid y)
=
-\frac1{2\sigma^2}
\|y-A\beta\|^2
-\lambda P(\beta)
+\text{constant}.
\]

因此 MAP 估计量等价于最小化

\[
\boxed{
\frac1{2\sigma^2}
\|y-A\beta\|^2
+
\lambda P(\beta)
}.
\]

所以：

- 若
  \[
  P(\beta)=\frac12\|\beta\|_2^2,
  \]
  则得到 ridge regression；

- 若
  \[
  P(\beta)=\|\beta\|_1,
  \]
  则得到 lasso 型软阈值；

- 若使用 $\ell_0$ penalty，则对应 hard thresholding，但优化更困难。

### 12.1 Spike-and-slab prior

在 Gaussian sequence model 中，

\[
X_i\mid\mu_i
\sim
N(\mu_i,1),
\]

可以使用

\[
\mu_i
\sim
(1-w)\delta_0+w\gamma,
\]

其中 $\delta_0$ 是 $0$ 点的点质量，$\gamma$ 是连续分布。

这种 spike-and-slab prior 显式编码稀疏性，并可通过 posterior median 或 MAP 产生收缩或阈值估计。

---

## 13. Empirical Bayes：从一组问题中学习先验

普通 Bayes 分析中，先验 $\nu$ 在用于当前决策的数据之前给定。

经验 Bayes（empirical Bayes, EB）则利用一组相关观测来估计：

- 先验分布 $\nu$；
- 先验的超参数；
- 或 Bayes rule 所需的边缘分布。

然后将估计结果代入 oracle Bayes rule。

设

\[
(\theta_i,X_i),
\qquad
i=1,\ldots,n,
\]

独立同分布于层次模型

\[
\theta_i\sim\nu,
\qquad
X_i\mid\theta_i\sim P_{\theta_i}.
\]

### 13.1 Historical-data EB

用历史样本

\[
X_1,\ldots,X_m
\]

估计 $\nu$，然后对一个新的观测 $X_{m+1}$ 使用拟合后的 Bayes rule。

条件于历史数据，这看起来就像使用了一个已经估计好的先验做普通 Bayes 推断。

### 13.2 Compound / parallel EB

同时观察很多当前单位

\[
X_1,\ldots,X_m,
\]

利用整个 ensemble 估计共同的 mixing distribution，再同时估计所有

\[
\theta_1,\ldots,\theta_m.
\]

这是现代大规模 shrinkage 与 multiple testing 方法中的典型场景。

### 13.3 Parametric 与 nonparametric EB

经验 Bayes 可以是：

- **parametric EB**：例如假设先验为正态分布，只估计其均值与方差；
- **nonparametric EB**：直接估计一个不受参数化限制的 mixing distribution。

一个重要思想是：有时并不需要估计完整的先验 $\nu$。如果 oracle Bayes rule 只依赖边缘分布 $f_\nu$，那么直接估计 $f_\nu$ 即可。

### 13.4 Poisson empirical Bayes

设

\[
X\mid\theta
\sim
\operatorname{Poisson}(\theta),
\qquad
\theta\sim\nu.
\]

边缘概率质量函数为

\[
f_\nu(x)
=
\int_0^\infty
e^{-\theta}
\frac{\theta^x}{x!}
\nu(d\theta).
\]

后验均值为

\[
\delta_\nu(x)
=
E(\theta\mid X=x).
\]

利用

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

得到 Robbins identity：

\[
\boxed{
\delta_\nu(x)
=
\frac{(x+1)f_\nu(x+1)}
{f_\nu(x)}
}.
\]

因此，经验 Bayes 可以直接利用大量观测估计

\[
f_\nu(x)
\quad\text{和}\quad
f_\nu(x+1),
\]

再代入上述比例。

!!! note "Empirical Bayes 改变了什么？"

    因为

    \[
    \widehat\nu
    \]

    本身是数据依赖的，所以

    \[
    \delta_{\widehat\nu}
    \]

    一般不再是关于某个固定先验的严格 Bayes rule。

    同一批数据既用于学习先验，又用于当前估计，可能带来：

    - hyperparameter uncertainty；
    - 过度自信；
    - 数据重复使用带来的风险分析问题。

    可采用 sample splitting、cross-fitting、显式风险分析，或者完整 hierarchical Bayes 来处理这些问题。

另一方面，EB 的优势也很明显：它可以根据观测到的整体人群结构，自适应地学习 pooling / shrinkage 的强度。

---

## 14. 为什么 Bayes 会自然出现在这里？

前几章的逻辑可以写成：

\[
\boxed{
\text{UMVU 只比较无偏估计量}
}
\]

\[
\Downarrow
\]

\[
\boxed{
\text{Admissibility 比较完整风险函数}
}
\]

\[
\Downarrow
\]

\[
\boxed{
\text{不同估计量的风险曲线往往交叉}
}
\]

\[
\Downarrow
\]

\[
\boxed{
\text{用先验 }\pi\text{ 对参数空间加权}
}
\]

\[
\Downarrow
\]

\[
\boxed{
\text{最小化 Bayes risk / posterior risk}
}
\]

可容许性只能排除那些被统一支配的规则，却通常无法在剩余的大量可容许规则中选出唯一方案。

先验 $\pi$ 为参数空间赋予权重，从而把函数值比较

\[
\theta\mapsto R(\theta,\delta)
\]

转化为标量准则

\[
r(\pi,\delta)
=
\int_\Theta
R(\theta,\delta)\pi(d\theta).
\]

### 14.1 两种理解 Bayes 的视角

#### 决策论视角

先验不一定必须首先解释为“主观信念”，也可以理解为参数空间上的权重分布。

从这个角度，Bayes 方法可以：

- 将风险函数之间的偏序变为一个标量优化问题；
- 产生 admissible / complete class 理论中的候选规则；
- 为 minimax risk 提供下界并连接 least-favorable prior；
- 自然解释 shrinkage 与 regularization。

#### 统计推断视角

先验则表示当前实验之前已有的信息：

\[
\theta\sim\pi,
\qquad
X\mid\theta\sim P_\theta.
\]

似然描述当前数据的信息，Bayes formula 将两者结合成后验

\[
\pi(d\theta\mid x).
\]

随后后验可以用于：

- 点估计；
- 不确定性量化；
- 预测；
- 后续科学决策。

### 14.2 加权风险为什么变成后验风险？

利用 tower property，

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

因此

\[
\boxed{
\text{weighted frequentist risk}
\Longleftrightarrow
\text{posterior expected loss}
}.
\]

这就是决策论观点与 Bayes 推断观点之间的桥梁。

---

## 15. Bayes 方法对统计学的主要贡献

Bayesian 方法在统计学中有多个不同用途：

1. **先验信息**：整合已有研究、专家知识、物理约束或外部数据。
2. **正则化与 pooling**：先验自然产生 shrinkage，并帮助相关参数共享信息。
3. **预测与不确定性**：后验及 posterior predictive distribution 能传播参数不确定性。
4. **层次模型与潜变量模型**：统一处理 random effects、mixtures、missing data 与 latent states。
5. **科学决策**：后验可以与非对称的医学、经济或科学损失结合选择行动。
6. **频率学派理论工具**：Bayes 与 generalized Bayes 构造常用于证明 admissibility、minimaxity 与 complete-class 结果。
7. **Empirical Bayes**：从一组相关问题中学习先验结构，并自适应 pooling 强度。

!!! warning "Bayes 分析并非自动正确"

    结论可能对以下因素敏感：

    - 先验；
    - 似然模型；
    - 层次结构假设。

    因此需要关注 posterior propriety、prior sensitivity、calibration 与 robustness。

---

## 16. 本章小结

本章的主要结论如下：

1. Bayes risk 是频率学派风险关于先验分布的加权平均。
2. Bayes rule 可以通过逐点最小化 posterior expected loss 得到。
3. 平方误差、绝对误差与 0–1 损失分别对应 posterior mean、median 与 mode/MAP。
4. 共轭先验使先验与后验保持同一函数形式，并可解释为 pseudo-observations。
5. Beta–Binomial、Gamma–Poisson 与 Normal–Normal 模型都展示了 posterior mean 的 shrinkage 结构。
6. Gamma–Exponential 例子清楚区分了 plug-in、MRU/UMVU 与 Bayes 三类优化目标。
7. Improper prior 可以产生 generalized Bayes rule，但不能与 proper Bayes 混为一谈。
8. 正态位置模型中的 $\bar X$ 是 flat prior 下的 generalized Bayes estimator，却不是任何 proper prior 下的 Bayes estimator。
9. Jeffreys prior 通过 Fisher information 构造，核心动机是参数变换下的不变性。
10. MAP 与 penalized optimization 紧密相连，ridge、lasso 与 sparse priors 都可以从 Bayes 角度理解。
11. Empirical Bayes 通过一组相关观测学习先验或边缘分布，并据此自适应 shrinkage。
12. Bayes 方法不是对频率学派风险分析的替代，而是当风险曲线无法统一排序时的一种自然决策论延伸。
