# 第四章：MVUE、U-统计量与信息下界

本章首先讨论最小风险无偏估计与最小方差无偏估计，并给出 MVUE 的协方差正交刻画。随后将结论推广至向量目标，介绍 U-统计量这一系统构造无偏估计量的方法。在此基础上，我们推导标量与多参数 Cramér-Rao 下界，研究下界的等号条件与可达性，最后介绍适用于非正则模型的 Hammersley-Chapman-Robbins 界。

---

## 1. MRU、MVUE 与比较类

设参数为 $\theta\in\Theta$，其中可能有 $\Theta\subseteq\mathbb{R}^k$，需要估计的标量目标为 $g(\theta)$。在平方误差损失下，

\[
L\{g(\theta),a\}=\{a-g(\theta)\}^2.
\]

估计量 $\delta$ 的风险可以分解为

\[
R(\theta,\delta)
=\operatorname{Var}_\theta(\delta)
+\operatorname{Bias}_\theta(\delta)^2.
\]

因此，在无偏估计量类中，

\[
R(\theta,\delta)=\operatorname{Var}_\theta(\delta).
\]

!!! info "定义 1.1（MRU 与 MVUE）"

    若无偏估计量 $\widehat g$ 在给定损失函数下，对任意其他无偏估计量和任意 $\theta\in\Theta$，其风险均不更大，则称 $\widehat g$ 为**最小风险无偏估计量**，简称 MRU。

    在平方误差损失下，MRU 称为**最小方差无偏估计量**，简称 MVUE，也常称 UMVU 估计量。

平方损失下，MRU 与 MVUE 完全等价。需要注意的是，参数 $\theta$ 的维数并不影响这一结论：即使参数属于多参数模型，只要目标 $g(\theta)$ 是标量，仍然可以用标量方差比较无偏估计量。

---

## 2. MVUE 的正交刻画

令

\[
\mathcal{Z}
=\left\{Z:E_\theta Z=0,\ \forall\theta\in\Theta\right\}
\]

表示在整个参数空间上均为零均值的估计量所构成的线性空间。

!!! success "定理 2.1（协方差刻画）"

    设 $\widehat g$ 是 $g(\theta)$ 的无偏估计量，并且对每个 $\theta$ 都具有有限方差。则 $\widehat g$ 是 MVUE，当且仅当对任意 $Z\in\mathcal{Z}$ 以及 $Z$ 具有有限方差的任意 $\theta_0$，均有

    \[
    \operatorname{Cov}_{\theta_0}(\widehat g,Z)=0.
    \]

??? proof "定理 2.1 的证明（点击展开）"

    首先假设 $\widehat g$ 是 MVUE。对任意 $Z\in\mathcal{Z}$ 和任意 $a\in\mathbb{R}$，扰动估计量

    \[
    \widehat g+aZ
    \]

    仍然是 $g(\theta)$ 的无偏估计量。在固定的 $\theta_0$ 处，

    \[
    \begin{aligned}
    &\operatorname{Var}_{\theta_0}(\widehat g+aZ)
    -\operatorname{Var}_{\theta_0}(\widehat g)\\
    &\qquad
    =a^2\operatorname{Var}_{\theta_0}(Z)
    +2a\operatorname{Cov}_{\theta_0}(\widehat g,Z).
    \end{aligned}
    \]

    由 MVUE 的最优性，上式对所有 $a\in\mathbb{R}$ 都非负。如果协方差不为零，取符号相反且绝对值充分小的 $a$，右端将为负，产生矛盾。因此

    \[
    \operatorname{Cov}_{\theta_0}(\widehat g,Z)=0.
    \]

    反之，假设上述正交条件成立。令 $D$ 是 $g(\theta)$ 的任意其他有限方差无偏估计量，并令

    \[
    Z=D-\widehat g.
    \]

    则 $Z\in\mathcal{Z}$，故

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

    因此 $\widehat g$ 是 MVUE。$\square$

!!! note "几何解释"

    在每个固定的参数值处，平方可积统计量构成一个 Hilbert 空间。MVUE 与所有在整个模型上期望为零的允许扰动正交。

    这里“在整个模型上”非常重要：如果某个统计量只在单个参数值处均值为零，将它加到原估计量上通常不能保持全局无偏性。

### 2.1 完备充分统计量与正交性

!!! success "推论 2.2（完备充分性蕴含正交性）"

    设 $T$ 是完备充分统计量。若 $\widehat g=a(T)$ 是 $g(\theta)$ 的无偏估计量，则 $\widehat g$ 满足定理 2.1 的正交条件，因而是 MVUE。

??? proof "推论 2.2 的证明（点击展开）"

    对任意 $Z\in\mathcal{Z}$，由充分性，可以选择一个不依赖于 $\theta$ 的共同版本

    \[
    q(T)=E_\theta(Z\mid T).
    \]

    由重期望公式，

    \[
    E_\theta q(T)=E_\theta Z=0,
    \qquad \forall\theta\in\Theta.
    \]

    由完备性，$q(T)=0$ 几乎处处。因此

    \[
    \begin{aligned}
    \operatorname{Cov}_\theta\{a(T),Z\}
    &=E_\theta\{a(T)Z\}\\
    &=E_\theta\left[a(T)E_\theta(Z\mid T)\right]\\
    &=E_\theta\{a(T)q(T)\}=0.
    \end{aligned}
    \]

    再由定理 2.1，$a(T)$ 是 MVUE。这正是 Lehmann-Scheffé 定理的正交形式。$\square$

---

## 3. 向量目标与多参数模型

需要区分参数的维数与目标的维数：

\[
\theta\in\mathbb{R}^k,
\qquad
g(\theta)\in\mathbb{R}^r.
\]

上一节的协方差刻画已经适用于任意 $k$ 下的标量目标。现在考虑

\[
\widehat g
=(\widehat g_1,\ldots,\widehat g_r)^\mathsf{T}
\]

是

\[
g=(g_1,\ldots,g_r)^\mathsf{T}
\]

的无偏估计量。

!!! info "定义 3.1（向量 MVUE）"

    若每个分量 $\widehat g_j$ 都是 $g_j$ 的标量 MVUE，则称 $\widehat g$ 为向量 MVUE。

对于对称矩阵，记

\[
A\preceq B
\]

表示 $B-A$ 为半正定矩阵。

!!! success "定理 3.2（向量 MVUE 的等价刻画）"

    假设所有相关的二阶矩均有限，则下列命题等价：

    1. 每个 $\widehat g_j$ 都是 $g_j$ 的 MVUE；
    2. 对每个 $c\in\mathbb{R}^r$，$c^\mathsf{T}\widehat g$ 都是 $c^\mathsf{T}g(\theta)$ 的 MVUE；
    3. 对 $g$ 的任意无偏向量估计量 $D$，均有

    \[
    \operatorname{Cov}_\theta(\widehat g)
    \preceq
    \operatorname{Cov}_\theta(D),
    \qquad \forall\theta.
    \]

??? proof "定理 3.2 的证明（点击展开）"

    假设条件 1 成立。对任意 $Z\in\mathcal{Z}$，由定理 2.1，

    \[
    \operatorname{Cov}_\theta(\widehat g_j,Z)=0,
    \qquad j=1,\ldots,r.
    \]

    因而对任意 $c\in\mathbb{R}^r$，

    \[
    \operatorname{Cov}_\theta(c^\mathsf{T}\widehat g,Z)
    =\sum_{j=1}^r c_j
    \operatorname{Cov}_\theta(\widehat g_j,Z)
    =0.
    \]

    再用定理 2.1，得到条件 2。

    若条件 2 成立，$D$ 是任意无偏向量估计量，则对每个 $c$，

    \[
    \begin{aligned}
    c^\mathsf{T}\operatorname{Cov}_\theta(\widehat g)c
    &=\operatorname{Var}_\theta(c^\mathsf{T}\widehat g)\\
    &\leq\operatorname{Var}_\theta(c^\mathsf{T}D)\\
    &=c^\mathsf{T}\operatorname{Cov}_\theta(D)c.
    \end{aligned}
    \]

    这正是半正定序意义下的协方差矩阵不等式，因此条件 3 成立。

    最后假设条件 3 成立。固定 $j$，令 $d$ 是 $g_j$ 的任意标量无偏估计量。只将 $\widehat g$ 的第 $j$ 个分量替换为 $d$，得到无偏向量估计量 $D$。在协方差矩阵不等式两侧取关于第 $j$ 个坐标向量 $e_j$ 的二次型，得到

    \[
    \operatorname{Var}_\theta(\widehat g_j)
    \leq\operatorname{Var}_\theta(d).
    \]

    故每个分量都是 MVUE，即条件 1 成立。$\square$

!!! success "推论 3.3（向量正交性）"

    向量估计量 $\widehat g$ 是向量 MVUE，当且仅当对每个具有有限方差的标量 $Z\in\mathcal{Z}$，均有

    \[
    \operatorname{Cov}_\theta(\widehat g,Z)
    =0\in\mathbb{R}^r.
    \]

对于 Euclidean 平方损失，在无偏估计量类中，

\[
E_\theta\|\widehat g-g(\theta)\|_2^2
=\operatorname{tr}\left\{
\operatorname{Cov}_\theta(\widehat g)
\right\}.
\]

但半正定序下的协方差矩阵比较比迹的比较更强，因为它控制了每个线性组合 $c^\mathsf{T}\widehat g$ 的方差。

!!! example "例 3.4（正态总体的均值与方差）"

    设

    \[
    X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}
    N(\mu,\sigma^2).
    \]

    由上一章可知，$(\overline X,S^2)$ 是 $(\mu,\sigma^2)$ 的完备充分统计量。因此它是向量 MVUE，并且

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

    非对角元为零，是因为

    \[
    \overline X\perp S^2.
    \]

    但向量 MVUE 并不意味着多参数 Cramér-Rao 下界一定可达，后文将再次讨论这个例子。

---

## 4. U-统计量：无偏估计量的系统构造

设

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}F.
\]

若目标泛函可以写成

\[
q(F)=E_Fh(X_1,\ldots,X_m),
\qquad m\leq n,
\]

其中核函数 $h$ 可积，则可以利用 $h$ 构造无偏估计量。

### 4.1 核函数的对称化

!!! success "命题 4.1（对称化不会增加方差）"

    定义对称化核函数

    \[
    h_s(x_1,\ldots,x_m)
    =\frac{1}{m!}\sum_\pi
    h(x_{\pi(1)},\ldots,x_{\pi(m)}),
    \]

    其中求和遍历所有排列 $\pi$。则

    \[
    E_Fh_s=E_Fh.
    \]

    当二阶矩存在时，还有

    \[
    \operatorname{Var}_F\{h_s(X_1,\ldots,X_m)\}
    \leq
    \operatorname{Var}_F\{h(X_1,\ldots,X_m)\}.
    \]

??? proof "命题 4.1 的证明（点击展开）"

    令 $\Pi$ 是独立于样本的均匀随机排列。由可交换性，

    \[
    h(X_{\Pi(1)},\ldots,X_{\Pi(m)})
    \]

    与 $h(X_1,\ldots,X_m)$ 同分布。给定样本后，

    \[
    \begin{aligned}
    h_s(X_1,\ldots,X_m)
    =E\bigl[&h(X_{\Pi(1)},\ldots,X_{\Pi(m)})\\
    &\mid X_1,\ldots,X_m\bigr].
    \end{aligned}
    \]

    重期望公式给出期望相等，条件方差分解给出方差不增。$\square$

因此，下文不妨假设核函数 $h$ 是对称的。

### 4.2 U-统计量的定义与无偏性

!!! info "定义 4.2（U-统计量）"

    由 $m$ 阶对称核 $h$ 生成的 U-统计量定义为

    \[
    U_n
    =\binom{n}{m}^{-1}
    \sum_{1\leq i_1<\cdots<i_m\leq n}
    h(X_{i_1},\ldots,X_{i_m}).
    \]

!!! success "命题 4.3（U-统计量的无偏性）"

    若 $q(F)=E_Fh(X_1,\ldots,X_m)$ 有限，则

    \[
    E_FU_n=q(F).
    \]

每个求和项与 $h(X_1,\ldots,X_m)$ 具有相同的期望，而求和共有 $\binom{n}{m}$ 项，因此结论立即成立。

### 4.3 U-统计量的 Rao-Blackwell 表示

!!! success "定理 4.4（U-统计量作为 Rao-Blackwell 估计量）"

    假设 $F$ 遍历所有关于某个固定非原子测度具有密度的分布，并且 $h$ 具有有限二阶矩。则

    \[
    U_n
    =E_F\left\{
    h(X_1,\ldots,X_m)
    \mid X_{(1)},\ldots,X_{(n)}
    \right\}.
    \]

    因此 $U_n$ 的方差不大于粗糙估计量 $h(X_1,\ldots,X_m)$ 的方差。在这个完全非参数模型中，$U_n$ 是 $q(F)$ 的唯一 MVUE。

??? proof "定理 4.4 的证明（点击展开）"

    给定次序统计量后，将观测值分配给原始标签的 $n!$ 种方式等可能。对这些分配求 $h$ 的条件平均，等价于对所有 $m$ 元子集等权平均。由于 $h$ 对称，所得结果恰好是 $U_n$。

    因而 $U_n$ 是 $h(X_1,\ldots,X_m)$ 关于次序统计量的条件期望。Rao-Blackwell 定理给出方差不增。

    附录 A 证明了次序统计量向量在这一完全非参数模型中是完备充分统计量。$U_n$ 是它的对称函数且无偏，因此由 Lehmann-Scheffé 定理，$U_n$ 是唯一 MVUE。$\square$

### 4.4 典型例子

!!! example "例 4.5（估计 $\mu^2$）"

    设

    \[
    \mu=E_FX,
    \qquad E_FX^2<\infty.
    \]

    因为

    \[
    E_F(X_1X_2)=\mu^2,
    \]

    所以取二阶核 $h(x_1,x_2)=x_1x_2$。相应的 U-统计量为

    \[
    \begin{aligned}
    U_n
    &=\binom{n}{2}^{-1}\sum_{i<j}X_iX_j\\
    &=\frac{\left(\sum_{i=1}^nX_i\right)^2
    -\sum_{i=1}^nX_i^2}{n(n-1)}.
    \end{aligned}
    \]

    它是 $\mu^2$ 的非参数 MVUE。

!!! example "例 4.6（估计方差）"

    取对称核

    \[
    h(x_1,x_2)=\frac{(x_1-x_2)^2}{2}.
    \]

    则

    \[
    E_Fh(X_1,X_2)=\operatorname{Var}_F(X).
    \]

    利用恒等式

    \[
    \sum_{i<j}(X_i-X_j)^2
    =n\sum_{i=1}^n(X_i-\overline X)^2,
    \]

    得到

    \[
    \begin{aligned}
    U_n
    &=\binom{n}{2}^{-1}
    \sum_{i<j}\frac{(X_i-X_j)^2}{2}\\
    &=\frac{1}{n-1}
    \sum_{i=1}^n(X_i-\overline X)^2.
    \end{aligned}
    \]

    因此通常的无偏样本方差是二阶 U-统计量。方差泛函的次数为 $2$；附录 B 将说明，在完全非参数模型中不存在只使用一个观测的方差无偏估计量。

!!! example "例 4.7（估计分布函数）"

    对固定的 $t$，取一阶核

    \[
    h_t(x)=\mathbf{1}_{\{x\leq t\}}.
    \]

    相应的一阶 U-统计量为

    \[
    F_n(t)=\frac{1}{n}\sum_{i=1}^n
    \mathbf{1}_{\{X_i\leq t\}},
    \]

    即经验分布函数在 $t$ 处的取值。

U-统计量还广泛用于秩检验、Kendall's tau、核两样本统计量以及各种现代成对估计量。

---

## 5. 标量 Cramér-Rao 不等式

设 $p_\theta$ 是关于共同控制测度的密度，参数 $\theta$ 位于某个开区间内。考虑以下正则条件：

!!! note "正则条件"

    **(R1)** 分布的支撑集不依赖于 $\theta$，并且在该支撑集上 $p_\theta>0$。

    **(R2)** $p_\theta(x)$ 关于 $\theta$ 可微，并且对 $p_\theta$ 和 $\delta p_\theta$ 的积分均允许交换微分与积分次序。

    **(R3)** 得分函数

    \[
    S_\theta(X)
    =\frac{\partial}{\partial\theta}
    \log p_\theta(X)
    \]

    平方可积，且 Fisher 信息满足

    \[
    0<I(\theta)
    =E_\theta S_\theta^2<\infty.
    \]

### 5.1 得分恒等式

!!! success "引理 5.1（得分恒等式）"

    在条件 (R1)-(R3) 下，

    \[
    E_\theta S_\theta=0.
    \]

    若 $\delta$ 是可微目标 $g(\theta)$ 的无偏估计量，并且条件 (R2) 也适用于 $\delta$，则

    \[
    \operatorname{Cov}_\theta(\delta,S_\theta)
    =g'(\theta).
    \]

??? proof "引理 5.1 的证明（点击展开）"

    因为

    \[
    S_\theta
    =\frac{\partial p_\theta/\partial\theta}{p_\theta},
    \]

    所以

    \[
    \begin{aligned}
    E_\theta S_\theta
    &=\int\frac{\partial p_\theta(x)}{\partial\theta}
    \,d\mu(x)\\
    &=\frac{d}{d\theta}\int p_\theta(x)\,d\mu(x)=0.
    \end{aligned}
    \]

    类似地，由 $E_\theta\delta=g(\theta)$，

    \[
    \begin{aligned}
    g'(\theta)
    &=\frac{d}{d\theta}
    \int\delta(x)p_\theta(x)\,d\mu(x)\\
    &=\int\delta(x)S_\theta(x)p_\theta(x)\,d\mu(x)\\
    &=E_\theta(\delta S_\theta).
    \end{aligned}
    \]

    又因为 $E_\theta S_\theta=0$，最后一个式子正是 $\operatorname{Cov}_\theta(\delta,S_\theta)$。$\square$

### 5.2 Cramér-Rao 下界

!!! success "定理 5.2（标量 Cramér-Rao 下界）"

    在引理 5.1 的条件下，$g(\theta)$ 的任意无偏估计量 $\delta$ 都满足

    \[
    \operatorname{Var}_\theta(\delta)
    \geq
    \frac{\{g'(\theta)\}^2}{I(\theta)}.
    \]

    对大小为 $n$ 的独立同分布样本，

    \[
    I_n(\theta)=nI_1(\theta).
    \]

??? proof "定理 5.2 的证明（点击展开）"

    由 Cauchy-Schwarz 不等式和得分恒等式，

    \[
    \begin{aligned}
    \{g'(\theta)\}^2
    &=\operatorname{Cov}_\theta(\delta,S_\theta)^2\\
    &\leq\operatorname{Var}_\theta(\delta)
    \operatorname{Var}_\theta(S_\theta)\\
    &=\operatorname{Var}_\theta(\delta)I(\theta).
    \end{aligned}
    \]

    移项即得结论。对于独立样本，样本得分是各观测得分之和。各项独立且均值为零，所以总得分的方差为 $nI_1(\theta)$。$\square$

!!! success "引理 5.3（期望 Hessian 恒等式）"

    若 $p_\theta$ 二阶可微，并且二阶导数也可以与积分交换，则

    \[
    I(\theta)
    =-E_\theta\left\{
    \frac{\partial^2}{\partial\theta^2}
    \log p_\theta(X)
    \right\}.
    \]

??? proof "引理 5.3 的证明（点击展开）"

    直接求导可得

    \[
    \frac{\partial^2}{\partial\theta^2}\log p_\theta
    =\frac{p_\theta''}{p_\theta}
    -\left(\frac{p_\theta'}{p_\theta}\right)^2.
    \]

    取期望后，第一项为

    \[
    \int p_\theta''\,d\mu
    =\frac{d^2}{d\theta^2}1=0,
    \]

    第二项为 $-I(\theta)$，故结论成立。$\square$

!!! info "定义 5.4（有限样本有效性）"

    若无偏估计量在参数值 $\theta$ 处使 Cramér-Rao 不等式取等号，则称它在 $\theta$ 处是**有效的**。若它在模型中每个 $\theta$ 处都取到等号，则称其为**全局有效的**。

每个全局有效的无偏估计量都是 MVUE，因为 Cramér-Rao 下界适用于所有无偏竞争者。但反过来不成立：MVUE 不一定达到 Cramér-Rao 下界。

---

## 6. 标量下界的等号条件与可达性

!!! success "定理 6.1（Cramér-Rao 下界的等号条件）"

    在标量 Cramér-Rao 正则条件下，不等式在 $\theta$ 处取等号，当且仅当

    \[
    \delta(X)-g(\theta)
    =\frac{g'(\theta)}{I(\theta)}S_\theta(X),
    \qquad P_\theta\text{-a.s.}
    \]

??? proof "定理 6.1 的证明（点击展开）"

    Cauchy-Schwarz 不等式取等号，当且仅当两个中心化随机变量线性相关。因此

    \[
    \delta-g(\theta)=a(\theta)S_\theta
    \]

    几乎处处成立。两侧与得分函数取协方差，并使用得分恒等式，得到

    \[
    g'(\theta)=a(\theta)I(\theta).
    \]

    因此

    \[
    a(\theta)=\frac{g'(\theta)}{I(\theta)}.
    \]

    反之，若等号条件成立，则

    \[
    \operatorname{Var}_\theta(\delta)
    =\frac{\{g'(\theta)\}^2}{I(\theta)^2}
    \operatorname{Var}_\theta(S_\theta)
    =\frac{\{g'(\theta)\}^2}{I(\theta)}.
    \]

    故 Cramér-Rao 下界取到等号。$\square$

!!! warning "可达性条件非常严格"

    等号条件左侧必须是一个不依赖未知参数的统计量，而右侧包含依赖参数的得分函数。因此，满足无偏性并不足以保证下界可达。

### 6.1 全局可达性与指数族

!!! success "定理 6.2（全局取等号蕴含指数族结构）"

    假设正则条件在某个区间上成立，$I(\theta)>0$ 且 $g'(\theta)\neq0$。若存在同一个估计量 $\delta$ 在每个 $\theta$ 处都达到 Cramér-Rao 下界，则模型必具有形式

    \[
    p_\theta(x)
    =\exp\{q(\theta)\delta(x)+r(\theta)+c(x)\}.
    \]

    因此，全局可达性要求模型是以有效估计量为自然统计量的一参数指数族。

??? proof "定理 6.2 的证明（点击展开）"

    将等号条件改写为

    \[
    \frac{\partial}{\partial\theta}
    \log p_\theta(x)
    =\frac{I(\theta)}{g'(\theta)}\delta(x)
    -\frac{I(\theta)g(\theta)}{g'(\theta)}.
    \]

    令右侧两个只依赖参数的系数分别为 $q'(\theta)$ 与 $r'(\theta)$，再关于 $\theta$ 积分，即得

    \[
    \log p_\theta(x)
    =q(\theta)\delta(x)+r(\theta)+c(x).
    \]

    指数化后得到结论。$\square$

上述定理只是必要性结论。对于正则指数族，可以进一步精确刻画哪些目标可达到下界。

!!! success "定理 6.3（一参数指数族中的有效仿射统计量）"

    设正则一参数指数族为

    \[
    p_\theta(x)
    =C(\theta)\exp\{Q(\theta)T(x)\}h(x),
    \qquad Q'(\theta)\neq0,
    \]

    且

    \[
    0<\operatorname{Var}_\theta(T)<\infty.
    \]

    对任意常数 $a,b$，估计量

    \[
    \delta=aT+b
    \]

    对目标

    \[
    g(\theta)=E_\theta(aT+b)
    \]

    是全局有效的。反之，在 $T$ 非退化的条件下，任何在整个族中均有效的同一个估计量都必须具有仿射形式 $aT+b$。

??? proof "定理 6.3 的证明（点击展开）"

    记

    \[
    m(\theta)=E_\theta T.
    \]

    对归一化恒等式求导，得到

    \[
    \frac{C'(\theta)}{C(\theta)}
    =-Q'(\theta)m(\theta).
    \]

    因此

    \[
    S_\theta
    =Q'(\theta)\{T-m(\theta)\},
    \]

    且

    \[
    I(\theta)
    =Q'(\theta)^2\operatorname{Var}_\theta(T).
    \]

    将得分恒等式用于 $T$，得到

    \[
    m'(\theta)
    =Q'(\theta)\operatorname{Var}_\theta(T).
    \]

    对 $g(\theta)=am(\theta)+b$，有

    \[
    g'(\theta)
    =aQ'(\theta)\operatorname{Var}_\theta(T).
    \]

    因而

    \[
    \frac{g'(\theta)}{I(\theta)}S_\theta
    =a\{T-m(\theta)\}
    =\delta-g(\theta).
    \]

    由定理 6.1，$\delta$ 全局有效。

    反之，由等号条件，

    \[
    \delta(x)-g(\theta)
    =\frac{g'(\theta)Q'(\theta)}{I(\theta)}
    \{T(x)-m(\theta)\}.
    \]

    因此对于每个 $\theta$，$\delta$ 都是 $T$ 的仿射函数。比较两个不同参数值处的表达式，并利用 $T$ 非退化，可知斜率和截距必须是不依赖参数的常数。$\square$

!!! warning "指数族结构仍不足以保证任意目标有效"

    即使在指数族中，非线性目标也未必存在有效无偏估计量。指数族结构只使等号成为可能；目标与自然统计量之间的仿射相容性才决定下界是否真正可达。

!!! example "例 6.4（UMVU 但不有效）"

    设

    \[
    X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}
    N(\mu,\sigma^2),
    \]

    其中 $\sigma^2$ 已知，并记

    \[
    v=\frac{\sigma^2}{n}.
    \]

    统计量

    \[
    \delta=\overline X^2-v
    \]

    是 $\mu^2$ 的 UMVU 估计量。因为 $\overline X\sim N(\mu,v)$，

    \[
    \operatorname{Var}_\mu(\delta)
    =4\mu^2v+2v^2.
    \]

    相应的 Cramér-Rao 下界为

    \[
    \frac{\{(\mu^2)'\}^2}{n/\sigma^2}
    =4\mu^2v.
    \]

    两者存在严格为正的差 $2v^2$，所以该目标的下界不可达。这个例子说明：MVUE 并不蕴含有效性。

---

## 7. 多参数 Cramér-Rao 不等式

设

\[
\theta=(\theta_1,\ldots,\theta_k)^\mathsf{T}
\]

属于 $\mathbb{R}^k$ 的一个开集。定义得分向量和 Fisher 信息矩阵为

\[
S_\theta=\nabla_\theta\log p_\theta(X),
\qquad
I(\theta)=E_\theta(S_\theta S_\theta^\mathsf{T}),
\]

并假设 $I(\theta)$ 正定。对于可微目标 $g:\mathbb{R}^k\to\mathbb{R}^r$，将其 Jacobian 矩阵记为

\[
G(\theta)
=\frac{\partial g(\theta)}{\partial\theta^\mathsf{T}}
\in\mathbb{R}^{r\times k}.
\]

!!! success "定理 7.1（多参数 Cramér-Rao 下界）"

    在 (R1)-(R3) 的向量版本成立时，$g(\theta)$ 的任意无偏向量估计量 $\delta\in\mathbb{R}^r$ 都满足

    \[
    \operatorname{Cov}_\theta(\delta)
    \succeq
    G(\theta)I(\theta)^{-1}G(\theta)^\mathsf{T}.
    \]

    等号成立，当且仅当

    \[
    \delta-g(\theta)
    =G(\theta)I(\theta)^{-1}S_\theta,
    \qquad P_\theta\text{-a.s.}
    \]

??? proof "定理 7.1 的证明（点击展开）"

    对每个无偏性恒等式求导，得到矩阵形式的得分恒等式

    \[
    E_\theta S_\theta=0,
    \qquad
    \operatorname{Cov}_\theta(\delta,S_\theta)=G(\theta).
    \]

    定义残差

    \[
    R
    =\delta-g(\theta)-G(\theta)I(\theta)^{-1}S_\theta.
    \]

    直接展开其协方差矩阵，得到

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

    任意协方差矩阵均为半正定矩阵，故得到下界。取等号当且仅当 $\operatorname{Cov}_\theta(R)=0$，又因为 $ER=0$，这等价于 $R=0$ 几乎处处。$\square$

### 7.1 含干扰参数时的标量目标

若目标为标量，即 $r=1$，则

\[
G(\theta)=\nabla g(\theta)^\mathsf{T},
\]

从而

\[
\operatorname{Var}_\theta(\delta)
\geq
\nabla g(\theta)^\mathsf{T}
I(\theta)^{-1}
\nabla g(\theta).
\]

这里必须使用包含所有参数，尤其是干扰参数的完整 Fisher 信息矩阵。

### 7.2 自然指数族

!!! success "推论 7.2（自然指数族中的有效估计量）"

    设正则自然指数族为

    \[
    p_\eta(x)
    =\exp\{\eta^\mathsf{T}T(x)-A(\eta)\}h(x),
    \qquad \eta\in\mathcal{H}\subseteq\mathbb{R}^k,
    \]

    且

    \[
    I(\eta)=\nabla^2A(\eta)
    \]

    正定。对常矩阵 $C\in\mathbb{R}^{r\times k}$ 和常向量 $b\in\mathbb{R}^r$，估计量

    \[
    \delta=CT+b
    \]

    对目标

    \[
    g(\eta)
    =CE_\eta T+b
    =C\nabla A(\eta)+b
    \]

    是有效的。若 $T$ 的仿射支撑为 $\mathbb{R}^k$，则这些恰好是具有全局有效无偏估计量的目标。

??? proof "推论 7.2 的证明（点击展开）"

    在自然指数族中，

    \[
    S_\eta=T-\nabla A(\eta),
    \]

    并且

    \[
    I(\eta)
    =\operatorname{Cov}_\eta(T)
    =\nabla^2A(\eta).
    \]

    对目标 $g(\eta)=C\nabla A(\eta)+b$，有

    \[
    G= C\nabla^2A=CI.
    \]

    因此

    \[
    GI^{-1}S_\eta
    =C\{T-E_\eta T\}
    =\delta-g(\eta),
    \]

    满足多参数下界的等号条件。

    反之，等号条件说明对每个 $\eta$，同一个估计量都是 $T$ 的仿射函数，其系数为 $G(\eta)I(\eta)^{-1}$。比较两个参数值下的仿射表示，并使用 $T$ 的全仿射支撑，可知系数和截距必须为常数。因此 $\delta=CT+b$，取期望便得到相应目标。$\square$

### 7.3 正态均值与方差：MVUE 但未完全有效

设

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}N(\mu,\tau),
\qquad \tau=\sigma^2,
\]

其中 $\mu$ 与 $\tau$ 均未知。得分函数的两个分量为

\[
S_\mu
=\frac{1}{\tau}\sum_{i=1}^n(X_i-\mu)
\]

以及

\[
S_\tau
=-\frac{n}{2\tau}
+\frac{1}{2\tau^2}
\sum_{i=1}^n(X_i-\mu)^2.
\]

中心正态分布的三阶矩为零，所以两个得分分量不相关。又因为

\[
\operatorname{Var}\{(X_i-\mu)^2\}=2\tau^2,
\]

故 Fisher 信息矩阵为

\[
I(\mu,\tau)
=
\begin{pmatrix}
n/\tau & 0\\
0 & n/(2\tau^2)
\end{pmatrix}.
\]

对恒等目标

\[
g(\mu,\tau)=(\mu,\tau)^\mathsf{T},
\]

矩阵下界为

\[
I(\mu,\tau)^{-1}
=
\begin{pmatrix}
\tau/n & 0\\
0 & 2\tau^2/n
\end{pmatrix}.
\]

而向量 MVUE $(\overline X,S^2)^\mathsf{T}$ 的协方差矩阵为

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

因此 $\overline X$ 达到了均值分量的下界，但 $S^2$ 没有达到方差分量的下界。事实上，估计 $\tau$ 时的等号条件要求

\[
\begin{aligned}
\delta_\tau
&=\tau+\frac{2\tau^2}{n}S_\tau\\
&=\frac{1}{n}\sum_{i=1}^n(X_i-\mu)^2.
\end{aligned}
\]

右端依赖未知的干扰参数 $\mu$，因此它不是一个统计量。不存在同一个无偏估计量能对所有 $(\mu,\tau)$ 满足该等号条件。

若 $\mu$ 已知，则上式确实是合法的估计量，并能达到下界。这个例子清楚地区分了 MVUE 与 Cramér-Rao 下界的可达性。

---

## 8. 非正则模型：Hammersley-Chapman-Robbins 界

Cramér-Rao 下界的证明需要对模型求导。当模型不可微，或者支撑集随参数变化时，可以使用有限差分形式的下界。

!!! success "定理 8.1（Hammersley-Chapman-Robbins 界）"

    假设

    \[
    P_{\theta+\Delta}\ll P_\theta,
    \]

    并定义似然比

    \[
    L_\Delta(X)
    =\frac{p_{\theta+\Delta}(X)}{p_\theta(X)}.
    \]

    若 $L_\Delta$ 在 $P_\theta$ 下具有有限二阶矩，且 $\delta$ 是 $g(\theta)$ 的无偏估计量，则

    \[
    \operatorname{Var}_\theta(\delta)
    \geq
    \frac{\{g(\theta+\Delta)-g(\theta)\}^2}
    {E_\theta\{L_\Delta(X)-1\}^2}.
    \]

    对所有允许的 $\Delta$ 取上确界，可得到这一形式下最强的下界。

??? proof "定理 8.1 的证明（点击展开）"

    由换测度公式，

    \[
    E_\theta L_\Delta=1,
    \]

    且

    \[
    E_\theta(\delta L_\Delta)
    =E_{\theta+\Delta}\delta
    =g(\theta+\Delta).
    \]

    令

    \[
    Z=L_\Delta-1.
    \]

    则

    \[
    E_\theta Z=0,
    \]

    并且

    \[
    \operatorname{Cov}_\theta(\delta,Z)
    =g(\theta+\Delta)-g(\theta).
    \]

    应用 Cauchy-Schwarz 不等式，

    \[
    \operatorname{Cov}_\theta(\delta,Z)^2
    \leq
    \operatorname{Var}_\theta(\delta)
    \operatorname{Var}_\theta(Z),
    \]

    即得所需下界。$\square$

在二次均值可微条件下，将分子与分母同时除以 $\Delta^2$，再令 $\Delta\to0$，即可恢复标量 Cramér-Rao 下界。

!!! warning "支撑集随参数变化的模型"

    对 $\operatorname{Unif}(0,\theta)$，支撑集随 $\theta$ 改变，条件 (R1)-(R2) 失效，因此不应直接使用通常的得分计算。此时最大次序统计量和有限差分型信息下界更合适。

---

## 附录 A：完全非参数模型中的完备次序统计量

设 $\mu$ 是非原子 $\sigma$-有限测度，$\mathcal{F}_\mu$ 是所有关于 $\mu$ 的概率密度构成的族。若

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}f,
\]

记次序统计量向量为

\[
O=(X_{(1)},\ldots,X_{(n)}).
\]

!!! success "定理 A.1（次序统计量的完备充分性）"

    在完全非参数族

    \[
    \left\{f^n:f\in\mathcal{F}_\mu\right\}
    \]

    中，次序统计量向量 $O$ 是完备充分统计量。

??? proof "定理 A.1 的证明（点击展开）"

    联合密度

    \[
    \prod_{i=1}^n f(x_i)
    \]

    是样本的对称函数，因此可以写成次序统计量的函数。由因子分解定理，$O$ 充分。

    为证明完备性，假设 $a(O)$ 可积，并且对所有 $f\in\mathcal{F}_\mu$ 都有

    \[
    E_fa(O)=0.
    \]

    定义对称函数

    \[
    H(x_1,\ldots,x_n)
    =a\{x_{(1)},\ldots,x_{(n)}\}.
    \]

    则对每个概率密度 $f$，

    \[
    \int H(x_1,\ldots,x_n)
    \prod_{i=1}^nf(x_i)
    \,d\mu(x_1)\cdots d\mu(x_n)=0.
    \]

    任取基准密度 $f_0$ 和密度 $f_1,\ldots,f_n$。对 $\alpha_j>0$，定义归一化正混合

    \[
    f_\alpha
    =\frac{f_0+\sum_{j=1}^n\alpha_jf_j}
    {1+\sum_{j=1}^n\alpha_j}.
    \]

    将 $f_\alpha$ 代入前式，并乘以

    \[
    \left(1+\sum_{j=1}^n\alpha_j\right)^n,
    \]

    得到关于自由正变量 $\alpha_1,\ldots,\alpha_n$ 的多项式。该多项式在 $\mathbb{R}^n$ 的一个开集上恒为零，所以每个系数都为零。

    其中 $\alpha_1\cdots\alpha_n$ 的系数为所有排列对应积分之和。由于 $H$ 对称，这些积分全部相等，故

    \[
    \int H(x_1,\ldots,x_n)
    \prod_{i=1}^nf_i(x_i)
    \,d\mu^n(x)=0
    \]

    对任意密度 $f_1,\ldots,f_n$ 成立。

    令每个 $f_i$ 为任意有限测度集合 $A_i$ 上的归一化示性函数，可知 $H$ 在每个有限测度矩形

    \[
    A_1\times\cdots\times A_n
    \]

    上的积分都为零。由单调类定理，相应的有符号测度在每个有限测度可测集上均为零；再由 $\sigma$-有限性，

    \[
    H=0,
    \qquad \mu^n\text{-a.e.}
    \]

    因而对每个 $f$，$a(O)=0$ 几乎处处，完备性得证。$\square$

!!! note "为什么需要非原子性？"

    非原子性保证样本中出现并列值的概率为零，从而给定次序统计量后，标签排列的条件分布是均匀的，这正是 U-统计量条件期望表示所使用的事实。

    若分布含有原子，则需要根据重复观测值修正排列轨道的计数，但对称化思想仍然成立。

---

## 附录 B：非参数无偏估计中的次数与不存在性

若 $m$ 是满足下式的最小正整数，

\[
q(F)=E_Fh(X_1,\ldots,X_m),
\]

则称无偏可估泛函 $q(F)$ 的**次数**为 $m$。

!!! success "命题 B.1（必要的混合多项式条件）"

    若 $q$ 的次数至多为 $m$，则对凸模型中的任意两个分布 $F,G$，映射

    \[
    \alpha
    \longmapsto
    q\{\alpha F+(1-\alpha)G\}
    \]

    必须是关于 $\alpha$ 的次数不超过 $m$ 的多项式。

??? proof "命题 B.1 的证明（点击展开）"

    令

    \[
    H_\alpha=\alpha F+(1-\alpha)G.
    \]

    则

    \[
    q(H_\alpha)
    =\int h(x_1,\ldots,x_m)
    \prod_{j=1}^m
    \{\alpha\,dF(x_j)+(1-\alpha)\,dG(x_j)\}.
    \]

    展开乘积后即得到一个次数不超过 $m$ 的多项式。$\square$

对方差泛函，有

\[
\begin{aligned}
\operatorname{Var}_{\alpha F+(1-\alpha)G}(X)
&=\alpha E_FX^2+(1-\alpha)E_GX^2\\
&\quad
-\{\alpha E_FX+(1-\alpha)E_GX\}^2.
\end{aligned}
\]

当 $F$ 与 $G$ 的均值不同时，它是 $\alpha$ 的真正二次多项式。因此，方差不可能具有次数 $1$；而例 4.6 的二阶核说明方差的次数恰为 $2$。

标准差通常不具有有限次数，因为上述二次多项式的平方根一般不是多项式。这说明 U-统计量只能为具有有限阶期望表示的泛函构造无偏估计量。

---

## 9. 信息下界的证明核查清单

在声称某个估计量达到 Cramér-Rao 下界前，应逐项检查：

1. 估计量是否在整个参数空间上无偏？
2. 支撑集是否与参数无关？所需的导数与积分换序是否合法？
3. 是否使用了整个样本的正确得分函数？
4. 在多参数问题中，是否使用了包含干扰参数的完整 Fisher 信息矩阵？
5. 当前结果是否只是一个下界？估计量是否真正满足标量或向量等号条件？
6. 等号公式是否定义了一个不依赖未知参数的统一统计量？
7. 是否清楚地区分了“有效”“MVUE”与“可容许”这三个概念？

---

## 10. 本章小结

- 在平方损失下，无偏估计量的风险就是方差，因此 MRU 与 MVUE 等价。
- MVUE 的核心刻画是：它与每个全局零均值无偏扰动在协方差意义下正交。
- 向量 MVUE 等价于每个线性组合均为标量 MVUE，也等价于协方差矩阵在 Loewner 序下最小。
- U-统计量通过对核函数在所有子样本上平均，系统地构造无偏估计量；在完全非参数模型中，它还是相应泛函的唯一 MVUE。
- Cramér-Rao 不等式由得分恒等式与 Cauchy-Schwarz 不等式导出，但达到下界还必须满足严格的等号条件。
- 全局有效性与指数族结构密切相关；即便是 MVUE，也未必能达到信息下界。
- 多参数模型中必须使用完整信息矩阵，不能忽略干扰参数。
- 当正则条件失效时，Hammersley-Chapman-Robbins 界提供了有限差分型的替代下界。
