# 第八章：相合性、矩估计与 Delta 方法

本章介绍大样本统计推断的基本工作流程。我们先回顾随机变量的收敛方式、连续映射定理与 Slutsky 定理，再讨论可识别性、分离条件和最小距离估计的相合性。随后研究经验原点矩与中心矩、矩估计的相合性和渐近正态性，最后给出多元 Delta 方法以及两个补充证明。

---

## 1. 大样本分析的工作流程

设

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}
P_\theta,
\]

并用统计量 $T_n$ 估计目标 $g(\theta)$。大样本分析依次回答三个问题：

| 问题 | 统计概念 |
| --- | --- |
| $T_n$ 是否趋近于 $g(\theta)$？ | 相合性 |
| $T_n-g(\theta)$ 以什么尺度波动？ | 收敛速度 |
| 极限分布是什么？ | 推断 |

因此，相合性不是最终目标，但它通常是进行局部展开的第一个必要条件。

---

## 2. 收敛方式与两个传递工具

设 $Z_n,Z$ 是某个度量空间中的随机向量，距离为 $d$。

!!! info "定义 2.1（四种收敛）"

    **依概率收敛：**

    \[
    Z_n\overset{p}{\longrightarrow}Z
    \]

    是指对每个 $\epsilon>0$，

    \[
    P\{d(Z_n,Z)>\epsilon\}\longrightarrow0.
    \]

    **几乎处处收敛：**

    \[
    Z_n\overset{\mathrm{a.s.}}{\longrightarrow}Z
    \]

    是指

    \[
    P\left(\lim_{n\to\infty}Z_n=Z\right)=1.
    \]

    **$L^r$ 收敛：**

    \[
    Z_n\overset{L^r}{\longrightarrow}Z
    \]

    是指

    \[
    E\,d(Z_n,Z)^r\longrightarrow0.
    \]

    **依分布收敛：**

    \[
    Z_n\overset{d}{\longrightarrow}Z
    \]

    是指对每个有界连续函数 $f$，

    \[
    Ef(Z_n)\longrightarrow Ef(Z).
    \]

标准蕴含关系为

\[
Z_n\overset{\mathrm{a.s.}}{\longrightarrow}Z
\quad\Longrightarrow\quad
Z_n\overset{p}{\longrightarrow}Z
\quad\Longrightarrow\quad
Z_n\overset{d}{\longrightarrow}Z,
\]

以及

\[
Z_n\overset{L^r}{\longrightarrow}Z
\quad\Longrightarrow\quad
Z_n\overset{p}{\longrightarrow}Z.
\]

!!! success "定理 2.2（连续映射定理与 Slutsky 定理）"

    若

    \[
    Z_n\overset{p}{\longrightarrow}Z
    \]

    且 $g$ 在 $Z$ 处几乎处处连续，则

    \[
    g(Z_n)\overset{p}{\longrightarrow}g(Z).
    \]

    对几乎处处收敛与依分布收敛，也有相应的连续映射结论。

    特别地，若

    \[
    Z_n\overset{d}{\longrightarrow}Z,
    \qquad
    W_n\overset{p}{\longrightarrow}c,
    \]

    则

    \[
    Z_n+W_n\overset{d}{\longrightarrow}Z+c,
    \]

    \[
    Z_nW_n\overset{d}{\longrightarrow}cZ,
    \]

    并且当 $c\neq0$ 时，

    \[
    \frac{Z_n}{W_n}
    \overset{d}{\longrightarrow}
    \frac{Z}{c}.
    \]

??? proof "证明思路（点击展开）"

    连续性把事件

    \[
    d\{g(Z_n),g(Z)\}>\epsilon
    \]

    转化为 $Z_n$ 落在 $Z$ 的某个充分小邻域之外的事件。

    对于 Slutsky 定理，先得到联合收敛

    \[
    (Z_n,W_n)
    \overset{d}{\longrightarrow}
    (Z,c),
    \]

    因为第二个分量收敛到常数。再分别对加法、乘法与除法应用连续映射定理。 $\square$

---

## 3. 相合性与可识别性

!!! info "定义 3.1（弱相合与强相合）"

    如果对每个 $\theta\in\Theta$，

    \[
    T_n\overset{p}{\longrightarrow}g(\theta),
    \]

    则称 $T_n$ 是 $g(\theta)$ 的弱相合估计量。

    如果对每个 $\theta\in\Theta$，

    \[
    T_n\overset{\mathrm{a.s.}}{\longrightarrow}g(\theta),
    \]

    则称 $T_n$ 是 $g(\theta)$ 的强相合估计量。

!!! info "定义 3.2（可识别性）"

    如果

    \[
    P_{\theta_1}=P_{\theta_2}
    \quad\Longrightarrow\quad
\theta_1=\theta_2.
    \]

    则称模型可识别。

    若只估计泛函 $g(\theta)$，所需的条件为

    \[
    P_{\theta_1}=P_{\theta_2}
    \quad\Longrightarrow\quad
    g(\theta_1)=g(\theta_2).
    \]

!!! success "命题 3.3（可识别性是必要条件）"

    若

    \[
    P_{\theta_1}=P_{\theta_2},
    \qquad
    g(\theta_1)\neq g(\theta_2),
    \]

    则不存在能在两个参数值处同时相合地估计 $g(\theta)$ 的估计量。

??? proof "命题 3.3 的证明（点击展开）"

    任意统计量 $T_n$ 在 $\theta_1$ 与 $\theta_2$ 下的分布相同。因此，若它在 $P_{\theta_1}$ 下依概率收敛到 $g(\theta_1)$，则在 $P_{\theta_2}$ 下也依概率收敛到同一个常数。

    另一方面，在 $\theta_2$ 处的相合性要求它依概率收敛到 $g(\theta_2)$。一个随机变量序列不可能依概率收敛到两个不同的常数。 $\square$

### 3.1 为什么可识别性不充分

可识别性只是逐点唯一性陈述。相合性还要求数据分布能够稳定地区分不同参数。

考虑一个刻意构造的病态例子：

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

映射 $h$ 是一一映射，因此模型可识别，但不存在对每个 $\theta$ 都相合的估计量。完整证明见第 9 节：证明使用了连续函数列的逐点极限必须具有许多连续点，而 $h$ 除了 $1/2$ 以外处处不连续。

!!! note "实际含义"

    这个例子并非实际应用中常见的模型。它的作用是区分两个逻辑命题：

    \[
    \text{可识别性}
    =\text{总体层面的唯一标签},
    \]

    \[
    \text{相合性}
    =\text{由经验数据进行稳定恢复}.
    \]

    在正则参数模型中，紧致性、连续性与分离条件通常提供所缺少的稳定性。

---

## 4. 分离条件与最小距离估计的相合性

设 $F_\theta$ 是 $P_\theta$ 下的分布函数，$F_n$ 是经验分布函数，并定义

\[
d(F,G)
=\sup_x|F(x)-G(x)|.
\]

在真实参数 $\theta_0$ 处，假设对每个 $\epsilon>0$，

\[
\eta(\theta_0,\epsilon)
:=
\inf_{\|\theta-\theta_0\|\geq\epsilon}
d(F_\theta,F_{\theta_0})
>0.
\tag{2}
\]

条件 (2) 常称为强分离或良好分离的识别条件。它强于逐点可识别性。

设 $\widehat{\theta}_n$ 是近似最小距离估计量，满足

\[
d(F_{\widehat{\theta}_n},F_n)
\leq
\inf_{\theta\in\Theta}d(F_\theta,F_n)
+r_n,
\qquad
r_n\overset{\mathrm{a.s.}}{\longrightarrow}0.
\tag{3}
\]

!!! success "定理 4.1（最小距离估计的强相合性）"

    在条件 (2)-(3) 下，

    \[
    \widehat{\theta}_n
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \theta_0.
    \]

    如果经验收敛与优化误差只依概率成立，则同一论证给出弱相合性。

??? proof "定理 4.1 的证明（点击展开）"

    由 Glivenko-Cantelli 定理，

    \[
    \Delta_n
    :=d(F_n,F_{\theta_0})
    \overset{\mathrm{a.s.}}{\longrightarrow}0.
    \]

    因为 $\theta_0$ 是 (3) 中允许的竞争者，

    \[
    d(F_{\widehat{\theta}_n},F_n)
    \leq\Delta_n+r_n.
    \]

    三角不等式给出

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

    固定 $\epsilon>0$。如果

    \[
    \|\widehat{\theta}_n-\theta_0\|
    \geq\epsilon,
    \]

    则条件 (2) 强迫

    \[
    d(F_{\widehat{\theta}_n},F_{\theta_0})
    \geq\eta(\theta_0,\epsilon)>0.
    \]

    在 $2\Delta_n+r_n\to0$ 的概率为一事件上，这最终不可能发生，因此

    \[
    \widehat{\theta}_n
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \theta_0.
    \]

    $\square$

### 4.1 强凸性在何处起作用

对一般准则估计量

\[
\widehat{\theta}_n
\in\operatorname*{arg\,min}_\theta Q_n(\theta),
\]

条件 (2) 的总体版本为

\[
\inf_{\|\theta-\theta_0\|\geq\epsilon}
\{Q(\theta)-Q(\theta_0)\}
>0.
\tag{4}
\]

若 $Q$ 是 $m$-强凸函数并在 $\theta_0$ 处达到最小值，则

\[
Q(\theta)-Q(\theta_0)
\geq
\frac{m}{2}\|\theta-\theta_0\|^2,
\tag{5}
\]

所以 (4) 成立。强凸性是充分条件，但并非必要条件；分布距离即使不光滑，也可能满足分离条件。

强凸性还给出定量结论。令

\[
\Delta_n
=\sup_\theta|Q_n(\theta)-Q(\theta)|,
\]

并假设

\[
Q_n(\widehat{\theta}_n)
\leq Q_n(\theta_0)+r_n.
\]

则

\[
\|\widehat{\theta}_n-\theta_0\|
\leq
\left(
\frac{4\Delta_n+2r_n}{m}
\right)^{1/2}.
\tag{6}
\]

事实上，

\[
Q(\widehat{\theta}_n)-Q(\theta_0)
\leq2\Delta_n+r_n,
\]

再结合 (5) 即得 (6)。因此，几乎处处一致收敛产生强相合性，而依概率一致收敛产生弱相合性。

---

## 5. 原点矩与中心矩

设

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}P,
\]

并假设所需矩存在。定义总体原点矩与中心矩

\[
\alpha_r=EX^r,
\qquad
\mu_r=E(X-\mu)^r,
\qquad
\mu=\alpha_1,
\]

以及相应的经验矩

\[
A_{r,n}
=\frac{1}{n}\sum_{i=1}^nX_i^r,
\qquad
M_{r,n}
=\frac{1}{n}\sum_{i=1}^n(X_i-\overline{X}_n)^r.
\]

“中心矩”是相对于均值的矩，原点矩是相对于零点的矩。

!!! success "定理 5.1（经验矩的强相合性）"

    若

    \[
    E|X|^r<\infty,
    \]

    则

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

    对任意有限个原点矩与中心矩，只要相应绝对矩存在，结论也联合成立。

??? proof "定理 5.1 的证明（点击展开）"

    对每个 $j\leq r$，有 $E|X|^j<\infty$，所以强大数定律给出

    \[
    A_{j,n}
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \alpha_j.
    \]

    令 $A_{0,n}=\alpha_0=1$。二项式恒等式给出

    \[
    M_{r,n}
    =\sum_{j=0}^r
    \binom{r}{j}
    (-\overline{X}_n)^{r-j}A_{j,n},
    \]

    以及

    \[
    \mu_r
    =\sum_{j=0}^r
    \binom{r}{j}
    (-\mu)^{r-j}\alpha_j.
    \tag{8}
    \]

    由于

    \[
    \overline{X}_n
    =A_{1,n}
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \mu,
    \]

    对 (8) 中的多项式应用连续映射定理，得到

    \[
    M_{r,n}
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \mu_r.
    \]

    对有限维向量应用同一论证，即得联合收敛。 $\square$

偏差与相合性是不同概念。$A_{r,n}$ 对 $\alpha_r$ 严格无偏，而 $M_{r,n}$ 通常对 $\mu_r$ 有偏。例如，

\[
EM_{2,n}
=\frac{n-1}{n}\mu_2.
\]

在

\[
E|X|^{2r}<\infty,
\]

的条件下，第 10 节的计算给出

\[
EM_{r,n}-\mu_r
=O(n^{-1}).
\tag{9}
\]

公式 (9) 对研究有限样本偏差有用，但对上面更简单的强相合性证明并不需要。

---

## 6. 两种矩估计

假设参数由一个同时包含原点矩、中心矩或二者的向量刻画：

\[
m(\theta)
=\bigl(
\alpha_{r_1}(\theta),\ldots,\alpha_{r_k}(\theta),
\mu_{t_1}(\theta),\ldots,\mu_{t_\ell}(\theta)
\bigr).
\]

经验向量 $\widehat{m}_n$ 用 $A_{r,n}$ 或 $M_{t,n}$ 替换相应总体矩。

矩估计量可以通过解方程

\[
m(\widehat{\theta}_n)
=\widehat{m}_n,
\]

得到；当存在局部逆映射 $h$ 时，也可以定义

\[
\widehat{\theta}_n
=h(\widehat{m}_n),
\qquad
h\{m(\theta)\}=\theta.
\]

!!! success "推论 6.1（矩估计量的相合性）"

    若所需绝对矩存在，且 $h$ 在 $m(\theta)$ 处连续，则

    \[
    \widehat{\theta}_n
    \overset{\mathrm{a.s.}}{\longrightarrow}
    \theta.
    \]

    如果只有 $\widehat{m}_n$ 的弱收敛，则结论相应为弱相合性。

??? proof "推论 6.1 的证明（点击展开）"

    由定理 5.1，

    \[
    \widehat{m}_n
    \overset{\mathrm{a.s.}}{\longrightarrow}
    m(\theta).
    \]

    对连续映射 $h$ 应用连续映射定理即可。 $\square$

!!! warning "矩的选择很重要"

    经典矩估计通常使用前 $k$ 个原点矩。当中心矩能够识别尺度或形状时，使用中心矩同样合理。

    但是，仅使用中心矩无法识别位置参数：$\mu_1=0$，而且所有中心矩都对位置平移不变。在位置-尺度问题中，一个自然选择是用一阶原点矩识别位置，用二阶中心矩识别尺度。

---

## 7. 经验矩的渐近正态性

对固定正整数 $r_1,\ldots,r_k$，假设

\[
E|X|^{2\max_jr_j}<\infty.
\]

多元中心极限定理给出

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

其中

\[
B_{ij}
=\alpha_{r_i+r_j}
-\alpha_{r_i}\alpha_{r_j}.
\tag{10}
\]

对于中心矩，利用多项式表示 (8) 与 Delta 方法可以得到特别有用的渐近线性展开。

!!! success "定理 7.1（中心矩的渐近展开）"

    若

    \[
    E|X|^{2r}<\infty,
    \]

    则

    \[
    \sqrt{n}(M_{r,n}-\mu_r)
    =\frac{1}{\sqrt{n}}
    \sum_{i=1}^n\psi_r(X_i)
    +o_p(1),
    \tag{11}
    \]

    其中

    \[
    \psi_r(x)
    =(x-\mu)^r-\mu_r
    -r\mu_{r-1}(x-\mu).
    \]

    因此，

    \[
    \sqrt{n}(M_{r,n}-\mu_r)
    \overset{d}{\longrightarrow}
    N\{0,\operatorname{Var}(\psi_r(X))\}.
    \]

??? proof "定理 7.1 的证明（点击展开）"

    令

    \[
    Y_i=X_i-\mu,
    \qquad
    \overline{Y}_n=\overline{X}_n-\mu.
    \]

    在 $\overline{Y}_n=0$ 附近作一阶展开：

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

    这里

    \[
    \overline{Y}_n=O_p(n^{-1/2}),
    \qquad
    \frac{1}{n}\sum_{i=1}^nY_i^{r-1}
    \overset{p}{\longrightarrow}\mu_{r-1}.
    \]

    因而减去 $\mu_r$ 并乘以 $\sqrt{n}$ 后，

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

    这就是 (11)，再应用一维中心极限定理即得结论。 $\square$

当 $r=2$ 时，$\mu_1=0$，所以

\[
\psi_2(x)
=(x-\mu)^2-\sigma^2.
\]

因此

\[
\sqrt{n}(M_{2,n}-\sigma^2)
\overset{d}{\longrightarrow}
N(0,\mu_4-\sigma^4).
\]

---

## 8. Delta 方法

!!! success "定理 8.1（多元 Delta 方法）"

    设 $\xi_n\in\mathbb{R}^k$ 满足

    \[
    \sqrt{n}(\xi_n-a)
    \overset{d}{\longrightarrow}
    N_k(0,B).
    \]

    若 $g:\mathbb{R}^k\to\mathbb{R}^m$ 在 $a$ 处可微，Jacobian 矩阵为

    \[
    C=Dg(a),
    \]

    则

    \[
    \sqrt{n}\{g(\xi_n)-g(a)\}
    \overset{d}{\longrightarrow}
    N_m(0,CBC^{\mathsf T}).
    \tag{12}
    \]

    更一般地，若 $s_n\to\infty$ 且

    \[
    s_n(\xi_n-a)
    \overset{d}{\longrightarrow}Z,
    \]

    则

    \[
    s_n\{g(\xi_n)-g(a)\}
    \overset{d}{\longrightarrow}
    Dg(a)Z.
    \]

??? proof "定理 8.1 的证明（点击展开）"

    可微性意味着

    \[
    g(x)-g(a)
    =C(x-a)+R(x),
    \]

    且

    \[
    \frac{\|R(x)\|}{\|x-a\|}
    \longrightarrow0
    \qquad
    (x\to a).
    \]

    因为

    \[
    \sqrt{n}(\xi_n-a)
    =O_p(1),
    \]

    所以

    \[
    \xi_n\overset{p}{\longrightarrow}a,
    \]

    并且

    \[
    \sqrt{n}R(\xi_n)
    =
    \frac{\|R(\xi_n)\|}{\|\xi_n-a\|}
    \sqrt{n}\|\xi_n-a\|
    =o_p(1).
    \]

    因而

    \[
    \sqrt{n}\{g(\xi_n)-g(a)\}
    =
    C\sqrt{n}(\xi_n-a)
    +o_p(1).
    \]

    由 Slutsky 定理得到 (12)。把 $\sqrt{n}$ 换成 $s_n$，同样得到一般速度的结论。 $\square$

### 8.1 例子：Poisson 均值的对数

若

\[
X_i
\overset{\mathrm{iid}}{\sim}
\operatorname{Poisson}(\lambda),
\]

则

\[
\sqrt{n}(\overline{X}_n-\lambda)
\overset{d}{\longrightarrow}
N(0,\lambda).
\]

取

\[
g(x)=\log x,
\qquad
g'(\lambda)=\frac{1}{\lambda}.
\]

Delta 方法给出

\[
\sqrt{n}
\{\log\overline{X}_n-\log\lambda\}
\overset{d}{\longrightarrow}
N(0,1/\lambda).
\]

---

## 9. 可识别模型为何不存在全局相合估计量

这里证明第 3.1 节中例子 (1) 的结论。

首先，$h$ 是单射。它在有理数集合内部和无理数集合内部显然分别是单射。若某个有理数 $q$ 与某个无理数 $u$ 满足

\[
q=1-u,
\]

则 $u=1-q$ 将是有理数，矛盾。

假设存在对每个 $\theta\in(0,1)$ 都相合的估计量 $T_n$。把它替换为

\[
0\vee T_n\wedge1
\]

不会破坏相合性，因此可以假设

\[
0\leq T_n\leq1.
\]

对 $p\in(0,1)$，定义

\[
H_n(p)=E_pT_n,
\]

其中期望基于 $n$ 个独立同分布的 $\operatorname{Bernoulli}(p)$ 观测计算。对固定 $n$，$H_n(p)$ 是 $p$ 的多项式，因此连续。

有界性与依概率收敛共同推出 $L^1$ 收敛。因此：

若 $p$ 为有理数，则在 $\theta=p$ 处的相合性给出

\[
H_n(p)\longrightarrow p.
\]

若 $p$ 为无理数，令 $\theta=1-p$，则 $\theta$ 仍为无理数，并且该参数下的数据分布为 $\operatorname{Bernoulli}(p)$，所以相合性给出

\[
H_n(p)\longrightarrow1-p.
\]

因此，连续函数列 $H_n$ 逐点收敛到 $h$。连续函数列在区间上的逐点极限是 Baire 第一类函数，而每个 Baire 第一类函数都有稠密的连续点集合。

但是，对任意 $p\neq1/2$，分别取趋近于 $p$ 的有理数列与无理数列，$h$ 的相应极限分别为 $p$ 与 $1-p$。所以 $h$ 在每个 $p\neq1/2$ 处不连续，与 Baire 定理矛盾。

因此，不存在全局相合估计量。

---

## 10. 经验中心矩的偏差阶

令

\[
Y_i=X_i-\mu,
\qquad
\overline{Y}_n
=\frac{1}{n}\sum_{i=1}^nY_i,
\qquad
A_{s,n}^Y
=\frac{1}{n}\sum_{i=1}^nY_i^s.
\]

则

\[
M_{r,n}
=\sum_{j=0}^r
\binom{r}{j}
(-\overline{Y}_n)^j
A_{r-j,n}^Y.
\tag{13}
\]

其中 $j=0$ 项的期望为 $\mu_r$。对于 $j=1$，由独立性与 $EY=0$，

\[
E(\overline{Y}_nA_{r-1,n}^Y)
=\frac{\mu_r}{n}.
\]

对于 $j\geq2$，由 Cauchy-Schwarz 不等式，

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

其中，独立中心随机变量之和的矩不等式与条件

\[
E|X|^{2r}<\infty
\]

保证上述两个界成立。对 (13) 取期望并对有限个项求和，得到

\[
EM_{r,n}-\mu_r
=O(n^{-1}).
\]

---

## 11. 本章总结

!!! success "核心结论"

    1. 大样本分析依次研究相合性、收敛速度与极限分布。

    2. 可识别性是相合估计存在的必要条件，但稳定恢复通常还需要连续性、紧致性或强分离条件。

    3. 经验原点矩与中心矩在相应绝对矩存在时强相合，因而连续可逆的矩估计量也强相合。

    4. 经验中心矩通常有偏，但在 $2r$ 阶绝对矩存在时，其偏差为 $O(n^{-1})$。

    5. 经验原点矩满足多元中心极限定理，经验中心矩具有影响函数形式的渐近线性展开。

    6. Delta 方法通过函数的一阶导数把估计量的极限分布传递给光滑变换。
