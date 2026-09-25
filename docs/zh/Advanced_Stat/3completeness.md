# 第三章：完备充分性与最优无偏估计

本章研究如何在无偏估计量这一受限类别中寻找最优估计量。核心工具是**完备性（completeness）**、**Rao-Blackwell 定理**与 **Lehmann-Scheffé 定理**。最后介绍 Basu 定理，并通过 Bernoulli、Poisson、Normal 与 Uniform 模型说明这些结果的具体应用。

---

## 1. 为什么研究无偏估计量

设目标参数为 $g(\theta)$，估计量为 $\delta(X)$，损失函数为 $L\{g(\theta),a\}$。估计量的风险定义为

\[
R(\theta,\delta)
=E_\theta L\{g(\theta),\delta(X)\}.
\]

在平方损失下，风险可以分解为

\[
R(\theta,\delta)
=\operatorname{Var}_\theta\{\delta(X)\}
+\operatorname{Bias}_\theta(\delta)^2,
\]

其中

\[
\operatorname{Bias}_\theta(\delta)
=E_\theta\delta(X)-g(\theta).
\]

### 1.1 为什么通常不存在无约束的一致最优估计量

!!! success "命题 1.1（无约束的一致最优估计量通常不存在）"

    假设概率测度族 $\{P_\theta:\theta\in\Theta\}$ 两两相互绝对连续，并且 $g(\theta)$ 不是常数。在平方损失下，不存在估计量 $\delta_0$，使其风险在每个参数点都不大于所有其他估计量的风险。

??? proof "命题 1.1 的证明（点击展开）"

    固定 $\theta_0$，并将 $\delta_0$ 与常数估计量

    \[
    \delta(X)\equiv g(\theta_0)
    \]

    比较。后者在 $\theta_0$ 处的风险为零。若 $\delta_0$ 对所有估计量都一致最优，则

    \[
    0\leq R(\theta_0,\delta_0)\leq0.
    \]

    因此

    \[
    \delta_0=g(\theta_0),
    \qquad P_{\theta_0}\text{-a.s.}
    \]

    由于各 $P_\theta$ 相互绝对连续，上式在每个 $P_\theta$ 下都几乎处处成立。再选择另一个参数值 $\theta_1$，同理可得

    \[
    \delta_0=g(\theta_1),
    \qquad P_{\theta_1}\text{-a.s.}
    \]

    这与 $g$ 不是常数矛盾。$\square$

为了得到有意义的比较，可以限制估计量的类别，例如考虑不变估计、Bayes 估计、minimax 估计或无偏估计。本章研究最后一种选择。

### 1.2 无偏估计、MRU 与 UMVU

!!! info "定义 1.2（无偏估计）"

    如果估计量 $\delta(X)$ 满足

    \[
    E_\theta\delta(X)=g(\theta),
    \qquad \forall\theta\in\Theta,
    \]

    则称 $\delta$ 是 $g(\theta)$ 的无偏估计量。

!!! info "定义 1.3（MRU 与 UMVU）"

    若无偏估计量 $\delta_0$ 在每个参数点的风险都不大于任何其他无偏估计量的风险，则称 $\delta_0$ 是**最小风险无偏估计量（minimum-risk unbiased estimator, MRU）**。

    在平方损失下，MRU 称为**一致最小方差无偏估计量（uniformly minimum-variance unbiased estimator, UMVU）**。

对无偏估计量，偏差项为零，因此

\[
R(\theta,\delta)=\operatorname{Var}_\theta(\delta).
\]

这里的最优性是受限最优性：UMVU 估计量只在无偏估计量类中最优；扩大比较范围后，它仍可能被某个有偏估计量支配。

### 1.3 无偏估计量是否存在

无偏估计本质上是一个逆问题：寻找函数 $\delta$，使得

\[
g(\theta)
=\int\delta(x)p_\theta(x)\,d\mu(x),
\qquad \forall\theta\in\Theta.
\]

并非每个目标函数 $g$ 都属于这个期望算子的值域。

!!! success "命题 1.4（Bernoulli 样本中的可估目标）"

    设

    \[
    X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(p),
    \qquad 0<p<1.
    \]

    实值函数 $g(p)$ 存在基于该样本的无偏估计量，当且仅当 $g$ 是次数不超过 $n$ 的关于 $p$ 的多项式。

??? proof "命题 1.4 的证明（点击展开）"

    记

    \[
    T=\sum_{i=1}^nX_i.
    \]

    对任意统计量 $D(X)$，令

    \[
    d_t=E\{D(X)\mid T=t\}.
    \]

    按照 $T$ 的取值对样本点分组，可得

    \[
    E_pD(X)
    =\sum_{t=0}^nd_t\binom ntp^t(1-p)^{n-t}.
    \]

    右端是次数不超过 $n$ 的多项式，因此必要性成立。

    反之，Bernstein 多项式

    \[
    \binom ntp^t(1-p)^{n-t},
    \qquad t=0,\ldots,n,
    \]

    构成次数不超过 $n$ 的多项式空间的一组基。将 $g$ 在这组基下展开，并把相应系数定义为 $d_t$，取 $D=d_T$，即可得到 $g(p)$ 的无偏估计量。$\square$

因此，当 $r\leq n$ 时，$p^r$ 可以被无偏估计；而 $1/p$、$\log p$ 以及多数非线性标准化形状参数都不能由固定大小为 $n$ 的 Bernoulli 样本无偏估计。

---

## 2. 完备性：由期望得到唯一性

### 2.1 完备与有界完备统计量

!!! info "定义 2.1（完备统计量）"

    统计量 $T$ 称为关于模型族 $\mathcal{P}$ **完备**，如果对任意满足可积性的函数 $a(T)$，只要

    \[
    E_\theta a(T)=0,
    \qquad \forall\theta\in\Theta,
    \]

    就有

    \[
    a(T)=0,
    \qquad P_\theta\text{-a.s.},\quad \forall\theta\in\Theta.
    \]

!!! info "定义 2.2（有界完备统计量）"

    如果上述结论只要求对有界函数 $a$ 成立，则称 $T$ **有界完备（boundedly complete）**。

完备性蕴含有界完备性。它与充分性描述的是不同性质：充分性表示 $T$ 之外不再含有参数信息；完备性表示模型族下的期望能够唯一识别 $T$ 的函数。

!!! success "引理 2.3（无偏函数的唯一性）"

    若 $T$ 完备，并且 $a_1(T)$ 与 $a_2(T)$ 都是同一目标函数 $g(\theta)$ 的无偏估计量，则

    \[
    a_1(T)=a_2(T),
    \qquad P_\theta\text{-a.s.},\quad \forall\theta\in\Theta.
    \]

??? proof "引理 2.3 的证明（点击展开）"

    差值 $a_1(T)-a_2(T)$ 可积，并且对所有 $\theta$ 都满足

    \[
    E_\theta\{a_1(T)-a_2(T)\}=0.
    \]

    由完备性的定义，差值几乎处处为零。$\square$

### 2.2 Binomial 分布族的完备性

!!! example "例 2.4（Binomial 完备性）"

    设

    \[
    T\sim\operatorname{Binomial}(n,p),
    \qquad 0<p<1.
    \]

    若对所有 $p$ 都有 $E_pa(T)=0$，则

    \[
    0=\sum_{t=0}^na(t)\binom ntp^t(1-p)^{n-t}.
    \]

    两边除以 $(1-p)^n$，并令

    \[
    z=\frac{p}{1-p}>0,
    \]

    得到

    \[
    0=\sum_{t=0}^na(t)\binom ntz^t,
    \qquad \forall z>0.
    \]

    一个多项式若在区间上恒为零，则所有系数均为零。因此 $a(t)=0$ 对 $t=0,\ldots,n$ 成立，故 Binomial 分布族完备。

### 2.3 Poisson 分布族的完备性

!!! example "例 2.5（Poisson 完备性）"

    设

    \[
    T\sim\operatorname{Poisson}(n\lambda),
    \qquad \lambda>0.
    \]

    若 $E_\lambda a(T)=0$，则

    \[
    0=e^{n\lambda}E_\lambda a(T)
    =\sum_{t=0}^{\infty}a(t)\frac{(n\lambda)^t}{t!}.
    \]

    可积性保证该幂级数对每个 $\lambda>0$ 绝对收敛。因为它在一个区间上恒为零，所以所有系数都为零，即 $a(t)=0$ 对所有 $t$ 成立。因此 Poisson 分布族完备。

---

## 3. 满指数族的完备充分性

### 3.1 完备充分性定理

!!! success "定理 3.1（满指数族中的完备充分统计量）"

    设模型具有指数族形式

    \[
    p_\eta(x)
    =\exp\{\eta^TT(x)-A(\eta)\}h(x),
    \qquad \eta\in H\subseteq\mathbb{R}^k,
    \]

    并且自然参数空间 $H$ 包含 $\mathbb{R}^k$ 中的一个非空开集。则 $T$ 对由该开集索引的模型族是充分且完备的。因此，对于任何包含该子族的更大模型族，$T$ 仍然完备。

??? proof "定理 3.1 的证明（点击展开）"

    由因子分解定理，$T$ 是充分统计量。

    令 $\nu$ 为基准测度 $h(x)\mu(dx)$ 在映射 $x\mapsto T(x)$ 下的像测度，则 $T$ 的诱导分布为

    \[
    P_\eta^T(dt)
    =e^{\eta^Tt-A(\eta)}\nu(dt).
    \]

    假设 $a(T)$ 可积，并且其期望在一个开集上恒为零。选择该开集的内点 $\eta_0$。对原点附近的任意 $u$，

    \[
    \int a(t)e^{(\eta_0+u)^Tt}\nu(dt)=0.
    \]

    定义有限符号测度

    \[
    M(dt)=a(t)e^{\eta_0^Tt}\nu(dt).
    \]

    则 $M$ 的双边 Laplace 变换满足

    \[
    L_M(u)=\int e^{u^Tt}M(dt)=0
    \]

    对原点某个邻域中的所有 $u$ 成立。由下面的局部 Laplace 变换唯一性引理，$M=0$。因为 $e^{\eta_0^Tt}>0$，可得

    \[
    a(t)=0,
    \qquad \nu\text{-a.e.}
    \]

    每个 $P_\eta^T$ 关于 $\nu$ 都有严格正的密度，因此 $a(T)=0$ 在每个 $P_\eta$ 下几乎处处成立。故 $T$ 完备。$\square$

!!! success "引理 3.2（局部 Laplace 变换的唯一性）"

    设 $M$ 是 $\mathbb{R}^k$ 上的有限符号测度，其指数矩在原点的某个邻域内有限。若

    \[
    \int e^{u^Tt}M(dt)=0
    \]

    在该邻域内成立，则 $M=0$。

??? proof "引理 3.2 的证明（点击展开）"

    对实部落在指数矩存在邻域内的复向量 $z$，定义

    \[
    L_M(z)=\int e^{z^Tt}M(dt).
    \]

    用稍大的实指数矩控制被积函数，可以证明 $L_M$ 是解析函数，并允许在积分号下求导。

    从 $L_M$ 在实开集上恒为零出发，逐个坐标应用单复变量的恒等定理，可以把零恒等式延拓到包含虚轴的连通复管状区域。因此，对每个 $s\in\mathbb{R}^k$，

    \[
    \widehat M(s)
    =L_M(is)
    =\int e^{is^Tt}M(dt)
    =0.
    \]

    由有限符号测度的 Fourier 变换唯一性，得到 $M=0$。$\square$

### 3.2 为什么需要“开集”条件

自然参数空间具有满仿射张成足以保证最小充分性，因为有限个似然比就可以恢复自然统计量。但是，完备性要求 Laplace 变换在一个邻域上相等。弯曲参数空间可能具有满仿射张成，却不包含任何开集。因此，最小充分性本身并不能推出完备性。

### 3.3 完整正态分布族

!!! example "例 3.3（完整正态分布族）"

    设

    \[
    X_1,\ldots,X_n
    \overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2).
    \]

    自然参数为

    \[
    \eta_1=\frac{\mu}{\sigma^2},
    \qquad
    \eta_2=-\frac{1}{2\sigma^2},
    \]

    自然参数空间为

    \[
    H=\mathbb{R}\times(-\infty,0),
    \]

    它是开集。因此

    \[
    \left(\sum_{i=1}^nX_i,\sum_{i=1}^nX_i^2\right)
    \]

    是 $(\mu,\sigma^2)$ 的完备充分统计量。

    当 $n\geq2$ 时，它与

    \[
    (\overline X,S^2),
    \qquad
    S^2=\frac{1}{n-1}\sum_{i=1}^n(X_i-\overline X)^2,
    \]

    一一对应。因此 $(\overline X,S^2)$ 也是完备充分统计量。

!!! note "一个常见误解"

    当 $\mu$ 与 $\sigma^2$ 都未知时，$(\overline X,S^2)$ 对完整正态分布族仍然完备。可能失去完备性的是自然参数空间不含开集的弯曲正态子模型。“多参数”并不意味着“不完备”。

---

## 4. Rao-Blackwell 改进

### 4.1 凸损失下的 Rao-Blackwell 定理

!!! success "定理 4.1（凸损失下的 Rao-Blackwell 定理）"

    设 $T$ 是充分统计量，$\delta(X)$ 是可积估计量，并定义

    \[
    \delta_T(T)=E_\theta\{\delta(X)\mid T\}.
    \]

    充分性保证右端可以选为一个不依赖于未知参数 $\theta$ 的 $T$ 的函数。若 $a\mapsto L\{g(\theta),a\}$ 是凸函数，则

    \[
    R(\theta,\delta_T)
    \leq R(\theta,\delta),
    \qquad \forall\theta\in\Theta.
    \]

    如果 $\delta$ 是 $g(\theta)$ 的无偏估计量，则 $\delta_T$ 仍然无偏。

??? proof "定理 4.1 的证明（点击展开）"

    充分性提供了一个与 $\theta$ 无关的条件分布核，因此条件期望可以视为同一个 $T$ 的可测函数。由条件 Jensen 不等式，

    \[
    L\{g(\theta),\delta_T(T)\}
    \leq E_\theta\left[L\{g(\theta),\delta(X)\}\mid T\right].
    \]

    两边取期望便得到风险不等式。

    另外，由塔式法则，

    \[
    E_\theta\delta_T
    =E_\theta\{E_\theta(\delta\mid T)\}
    =E_\theta\delta.
    \]

    因此无偏性得以保持。$\square$

### 4.2 平方损失下的精确改进量

!!! success "推论 4.2（平方损失下的精确风险改进）"

    在平方损失下，

    \[
    R(\theta,\delta)-R(\theta,\delta_T)
    =E_\theta\{\operatorname{Var}_\theta(\delta\mid T)\}
    \geq0.
    \]

    除非 $\delta$ 已经是 $T$ 的函数（忽略 $P_\theta$-零测集），否则改进是严格的。

??? proof "推论 4.2 的证明（点击展开）"

    写成

    \[
    \delta-g(\theta)
    =(\delta-\delta_T)+(\delta_T-g(\theta)).
    \]

    第一项关于 $T$ 的条件期望为零，因此展开平方后交叉项的期望为零。于是

    \[
    E_\theta(\delta-g(\theta))^2
    =E_\theta(\delta-\delta_T)^2
    +E_\theta(\delta_T-g(\theta))^2.
    \]

    其中

    \[
    E_\theta(\delta-\delta_T)^2
    =E_\theta\{\operatorname{Var}_\theta(\delta\mid T)\}.
    \]

    整理即可。$\square$

这个等式正是条件方差分解，也解释了 Rao-Blackwell 化与条件 Monte Carlo、现代方差缩减方法之间的联系：它平均掉了不影响推断目标的额外随机性。

### 4.3 Bernoulli 幂函数的 Rao-Blackwell 化

!!! example "例 4.3（Bernoulli 幂函数）"

    设

    \[
    X_1,\ldots,X_n
    \overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(p),
    \qquad
    T=\sum_{i=1}^nX_i.
    \]

    当 $r\leq n$ 时，$X_1\cdots X_r$ 是 $p^r$ 的无偏估计量。给定 $T=t$ 后，含有 $t$ 个 $1$ 的所有排列等可能，因此

    \[
    E(X_1\cdots X_r\mid T=t)
    =\frac{(t)_r}{(n)_r},
    \]

    其中下降阶乘定义为

    \[
    (x)_r=x(x-1)\cdots(x-r+1).
    \]

    因而 $(T)_r/(n)_r$ 仍然无偏，并且方差不大于原始乘积估计量。

---

## 5. Lehmann-Scheffé 定理

!!! success "定理 5.1（凸损失下的 Lehmann-Scheffé 定理）"

    设 $T$ 是完备充分统计量。若 $\delta_0(T)$ 是 $g(\theta)$ 的无偏估计量，并且风险有限，则 $\delta_0$ 在任何关于行动 $a$ 为凸函数的损失下都是 MRU 估计量。

    在平方损失下，$\delta_0$ 是唯一的 UMVU 估计量。更一般地，在严格凸损失下，$\delta_0$ 是唯一的 MRU 估计量，其中唯一性忽略各模型共同的零测集。

??? proof "定理 5.1 的证明（点击展开）"

    任取 $g(\theta)$ 的无偏估计量 $D(X)$。对其进行 Rao-Blackwell 化：

    \[
    D_T(T)=E(D\mid T).
    \]

    则 $D_T$ 仍然无偏，并且风险不大于 $D$。由于 $D_T(T)$ 与 $\delta_0(T)$ 都是完备统计量 $T$ 的无偏函数，引理 2.3 给出

    \[
    D_T(T)=\delta_0(T),
    \qquad P_\theta\text{-a.s.},\quad \forall\theta.
    \]

    因此

    \[
    R(\theta,\delta_0)
    =R(\theta,D_T)
    \leq R(\theta,D),
    \]

    从而证明 MRU 最优性。

    在平方损失下，推论 4.2 表明等号成立只能发生在 $D=D_T=\delta_0$ 几乎处处时。严格凸损失下，由 Jensen 不等式的等号条件可得同样结论。$\square$

!!! note "UMVU 的标准构造步骤"

    (i) 找到一个完备充分统计量 $T$；

    (ii) 找到目标函数的任意无偏估计量 $D(X)$；

    (iii) 计算 $E(D\mid T)$，或者直接找到一个关于 $T$ 的无偏函数；

    (iv) 使用 Lehmann-Scheffé 定理得到唯一性与 UMVU 最优性。

---

## 6. 四类代表性 UMVU 构造

### 6.1 Bernoulli 幂函数

统计量

\[
T=\sum_{i=1}^nX_i\sim\operatorname{Binomial}(n,p)
\]

由因子分解定理知是充分统计量，由例 2.4 知是完备统计量。因此

\[
\widehat{p^r}_{\mathrm{UMVU}}
=\frac{(T)_r}{(n)_r},
\qquad r=1,\ldots,n.
\]

### 6.2 Poisson 点概率

设

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}\operatorname{Poisson}(\lambda),
\qquad
T=\sum_{i=1}^nX_i.
\]

对固定的 $k\geq0$，示性函数 $\mathbf{1}\{X_1=k\}$ 无偏估计

\[
g_k(\lambda)
=P_\lambda(X_1=k)
=e^{-\lambda}\frac{\lambda^k}{k!}.
\]

给定 $T=t$ 后，$(X_1,\ldots,X_n)$ 服从总数为 $t$、各格概率为 $1/n$ 的 multinomial 分布。因此

\[
\widehat g_k(T)
=\binom Tk
\left(\frac1n\right)^k
\left(1-\frac1n\right)^{T-k}
\mathbf{1}\{T\geq k\}.
\]

由 Poisson 完备性和 Lehmann-Scheffé 定理，$\widehat g_k(T)$ 是 $g_k(\lambda)$ 的 UMVU 估计量。

### 6.3 正态总体的均值与方差

设

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2).
\]

由例 3.3，$(\overline X,S^2)$ 是完整正态分布族的完备充分统计量。又因为

\[
E_{\mu,\sigma^2}\overline X=\mu,
\qquad
E_{\mu,\sigma^2}S^2=\sigma^2,
\]

所以 $\overline X$ 与 $S^2$ 分别是 $\mu$ 与 $\sigma^2$ 的 UMVU 估计量。

此外，

\[
E\overline X^2
=\mu^2+\frac{\sigma^2}{n},
\]

因此

\[
\overline X^2-\frac{S^2}{n}
\]

是 $\mu^2$ 的无偏估计量。它是完备充分统计量的函数，所以也是 $\mu^2$ 的 UMVU 估计量。

### 6.4 Uniform 单端点模型

设

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}\operatorname{Uniform}(0,\theta),
\qquad
M=X_{(n)}.
\]

由因子分解定理，$M$ 是充分统计量。若对每个 $\theta>0$ 都有 $E_\theta a(M)=0$，则

\[
0=\frac{n}{\theta^n}
\int_0^\theta a(m)m^{n-1}\,dm.
\]

因此，对每个 $\theta$，上式中的不定积分均为零。由 Lebesgue 意义下的微分，

\[
a(\theta)\theta^{n-1}=0
\]

几乎处处成立，所以 $M$ 完备。

另一方面，当 $r>-n$ 时，

\[
E_\theta M^r
=\frac{n}{n+r}\theta^r.
\]

因此，只要目标函数有定义，

\[
\widehat{\theta^r}_{\mathrm{UMVU}}
=\frac{n+r}{n}M^r.
\]

这个非正则例子表明，完备充分性并不局限于具有固定支撑集的指数族。

---

## 7. Basu 定理：完备性与辅助统计量

### 7.1 辅助统计量

!!! info "定义 7.1（辅助统计量）"

    如果统计量 $A(X)$ 的分布不依赖于参数，则称 $A$ 是**辅助统计量（ancillary statistic）**。

辅助统计量仍然可以反映估计精度、设计不平衡或抽样几何结构。“分布与参数无关”并不等于“没有用”。

### 7.2 Basu 定理

!!! success "定理 7.2（Basu 定理）"

    若 $T$ 是充分且有界完备的统计量，$A$ 是辅助统计量，则在每个 $P_\theta$ 下，$T$ 与 $A$ 相互独立。

??? proof "定理 7.2 的证明（点击展开）"

    固定 $A$ 的值域中的可测集 $B$。由充分性，可以选择一个与 $\theta$ 无关的版本

    \[
    q_B(T)=P_\theta(A\in B\mid T).
    \]

    由辅助性，

    \[
    P_\theta(A\in B)=c_B
    \]

    也与 $\theta$ 无关。因此

    \[
    E_\theta\{q_B(T)-c_B\}=0,
    \qquad \forall\theta.
    \]

    函数 $q_B-c_B$ 有界。由有界完备性，

    \[
    q_B(T)=c_B,
    \qquad P_\theta\text{-a.s.}
    \]

    对 $T$ 的值域中的任意可测集 $C$，

    \[
    \begin{aligned}
    P_\theta(A\in B,T\in C)
    &=E_\theta\left[\mathbf{1}\{T\in C\}q_B(T)\right]\\
    &=c_BP_\theta(T\in C)\\
    &=P_\theta(A\in B)P_\theta(T\in C).
    \end{aligned}
    \]

    因此 $T$ 与 $A$ 独立。$\square$

### 7.3 为什么 $\overline X$ 与 $S^2$ 独立

!!! example "例 7.3（正态样本均值与样本方差的独立性）"

    固定 $\sigma^2$，考虑只由 $\mu$ 索引的子模型

    \[
    X_1,\ldots,X_n
    \overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2).
    \]

    此时 $\overline X$ 是关于 $\mu$ 的完备充分统计量，而

    \[
    \frac{(n-1)S^2}{\sigma^2}
    \sim\chi_{n-1}^2
    \]

    的分布与 $\mu$ 无关，所以 $S^2$ 是该子模型中的辅助统计量。由 Basu 定理，

    \[
    \overline X\perp S^2.
    \]

    由于该结论对每个固定的 $\sigma^2$ 都成立，所以它在完整正态分布族中成立。

---

## 8. 函数型目标：估计正态分布函数

设

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2),
\]

其中 $\sigma^2$ 已知、$n\geq2$，并固定 $y\in\mathbb{R}$。示性函数

\[
D=\mathbf{1}\{X_1\leq y\}
\]

无偏估计

\[
F_\mu(y)
=\Phi\left(\frac{y-\mu}{\sigma}\right).
\]

此时完备充分统计量为 $\overline X$。由联合正态分布的条件分布公式，

\[
X_1\mid\overline X=x
\sim N\left(x,\sigma^2\left(1-\frac1n\right)\right).
\]

因此，Rao-Blackwell 化给出

\[
\widehat F_{\mathrm{UMVU}}(y)
=\Phi\left(
\frac{y-\overline X}
{\sigma\sqrt{1-1/n}}
\right).
\]

类似的 Gaussian 卷积恒等式给出总体密度在 $y$ 处的 UMVU 估计量：

\[
\widehat f_{\mathrm{UMVU}}(y)
=\frac{1}{\sigma\sqrt{1-1/n}}
\phi\left(
\frac{y-\overline X}
{\sigma\sqrt{1-1/n}}
\right).
\]

对上式取期望，相当于把 $N(\mu,\sigma^2/n)$ 的密度与 $N(0,\sigma^2(1-1/n))$ 的密度卷积，结果正是 $N(\mu,\sigma^2)$ 在 $y$ 处的密度。完备性进一步保证该估计量是唯一的 UMVU 估计量。

---

## 9. Uniform 双端点模型

设

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}\operatorname{Uniform}(a,b),
\qquad a<b,\quad n\geq2,
\]

并记

\[
U=X_{(1)},
\qquad
V=X_{(n)}.
\]

### 9.1 极端次序统计量的完备充分性

!!! success "定理 9.1（极端次序统计量的完备充分性）"

    统计量 $(U,V)$ 是参数 $(a,b)$ 的完备充分统计量。

??? proof "定理 9.1 的证明（点击展开）"

    样本联合密度为

    \[
    (b-a)^{-n}\mathbf{1}\{a<U<V<b\},
    \]

    因此由因子分解定理，$(U,V)$ 充分。

    $(U,V)$ 的联合密度为

    \[
    f_{a,b}(u,v)
    =\frac{n(n-1)(v-u)^{n-2}}{(b-a)^n}
    \mathbf{1}\{a<u<v<b\}.
    \]

    假设对所有 $a<b$ 都有

    \[
    E_{a,b}c(U,V)=0.
    \]

    两边乘以 $(b-a)^n$，得到

    \[
    G(a,b)
    =\int_a^b\int_u^b
    c(u,v)n(n-1)(v-u)^{n-2}\,dv\,du
    =0.
    \]

    令

    \[
    q(u,v)=c(u,v)n(n-1)(v-u)^{n-2}.
    \]

    绝对可积性允许在 Lebesgue 意义下求导。对几乎处处的 $a<b$，

    \[
    \frac{\partial G}{\partial a}(a,b)
    =-\int_a^bq(a,v)\,dv,
    \]

    并且

    \[
    \frac{\partial^2G}{\partial b\,\partial a}(a,b)
    =-q(a,b).
    \]

    因为 $G\equiv0$，所以 $q(a,b)=0$ 几乎处处。权重 $n(n-1)(b-a)^{n-2}$ 在 $a<b$ 时为正，因此 $c(a,b)=0$ 几乎处处。再由 $(U,V)$ 的联合密度可知

    \[
    c(U,V)=0,
    \qquad P_{a,b}\text{-a.s.},
    \]

    对每个参数对 $(a,b)$ 都成立。因此 $(U,V)$ 完备。$\square$

### 9.2 两个端点的 UMVU 估计

对标准 Uniform$(0,1)$ 样本，

\[
EU=\frac{1}{n+1},
\qquad
EV=\frac{n}{n+1}.
\]

经线性变换可得

\[
E_{a,b}U=\frac{na+b}{n+1},
\qquad
E_{a,b}V=\frac{a+nb}{n+1}.
\]

解这两个线性方程，得到

\[
\widehat a_{\mathrm{UMVU}}
=\frac{nU-V}{n-1},
\qquad
\widehat b_{\mathrm{UMVU}}
=\frac{nV-U}{n-1}.
\]

由定理 9.1 与 Lehmann-Scheffé 定理，它们分别是 $a$ 与 $b$ 的 UMVU 估计量。

---

## 10. 凸损失下不需要随机化

设随机化估计量在观察到 $X=x$ 后，按照分布 $Q_x(da)$ 选择行动。定义该随机行动的重心

\[
\overline\delta(x)=\int a\,Q_x(da).
\]

若随机化估计量无偏，则

\[
E_\theta\overline\delta(X)
=E_\theta\int a\,Q_X(da)
=g(\theta),
\]

所以 $\overline\delta$ 是普通的非随机化无偏估计量。

由凸性和 Jensen 不等式，对每个 $x$，

\[
L\{g(\theta),\overline\delta(x)\}
\leq\int L\{g(\theta),a\}\,Q_x(da).
\]

再对 $X$ 积分，得到非随机化重心估计量的风险不大于原随机化估计量。因此，在凸损失下，随机化不会给无偏估计带来额外收益。

这与离散假设检验不同：在离散检验中，行动空间与约束结构不同，随机化仍可能十分重要。

---

## 11. 本章证明核查清单

在断言某个估计量为 UMVU，或者使用独立性定理之前，应检查以下问题：

1. 目标函数在整个参数空间上是否存在无偏估计量？
2. 用来条件化的统计量是否充分，从而使 $E_\theta(\delta\mid T)$ 能够选成同一个与参数无关的统计量？
3. 该统计量是完备的，还是仅仅最小充分？
4. 使用指数族完备性时，自然参数空间是否包含开集，而不仅仅具有满仿射张成？
5. 期望、幂级数或 Laplace 变换所需的可积性条件是否成立？
6. Rao-Blackwell 定理所用的损失函数是否凸？只有排除等号情形后，才能声称改进是严格的。
7. 在引用 Lehmann-Scheffé 定理之前，候选估计量是否既无偏，又是完备充分统计量的函数？
8. 使用 Basu 定理时，第二个统计量是否确实关于同一个参数族辅助，而 $T$ 是否关于该参数族有界完备？
9. 最终结论只是无偏估计量类中的最优性，还是错误地声称了在所有估计量中的可容许性？

---

## 12. 本章小结

本章的逻辑可以概括为：

\[
\text{充分性}
+\text{Rao-Blackwell 改进}
+\text{完备性所保证的唯一性}
\Longrightarrow
\text{Lehmann-Scheffé 最优性}.
\]

完备充分统计量既提供了降维结构，又保证无偏函数的唯一性。Rao-Blackwell 定理把任意无偏估计量改进为充分统计量的函数，Lehmann-Scheffé 定理进一步把该函数识别为唯一的 UMVU 估计量。Basu 定理则揭示了完备充分统计量与辅助统计量之间的独立关系。
