# 第五章：可容许性、收缩与 James–Stein 现象

本章从决策论视角重新审视估计量的“最优性”。前几章主要研究无偏估计量中的最优性，例如 MVUE 与 UMVU；本章则去掉无偏性限制，直接比较估计量的完整风险函数。由此引出支配、可容许性、Rao–Blackwell 风险改进、收缩估计、不变性以及著名的 James–Stein 现象。

---

## 1. 从无偏性走向风险最优性

在平方误差损失下，估计量 $\delta(X)$ 的风险可写成

\[
R(\theta,\delta)
=\operatorname{Var}_\theta\{\delta(X)\}
+\operatorname{Bias}_\theta\{\delta(X)\}^2.
\]

如果 $\delta$ 无偏，则风险就是方差，因此前面的 MVUE/UMVU 理论本质上是在无偏估计量类中比较方差。

但是，无偏性只是一个约束，并不等同于全局风险最优。一个带有少量偏差的估计量，可能通过显著降低方差而获得更小的均方误差。

!!! note "本章的核心问题"

    不再问：

    > 哪一个无偏估计量具有最小方差？

    而是问：

    > 是否存在另一个估计量，使得它对所有参数值的风险都不更大，并且在某些参数值上严格更小？

这正是可容许性与收缩估计理论的出发点。

---

## 2. 决策论框架

设

\[
X\sim P_\theta,\qquad \theta\in\Theta,
\]

行动空间为 $\mathcal D$。一个决策规则是可测函数

\[
\delta:\mathcal X\to\mathcal D.
\]

给定损失函数 $L(\theta,a)$，决策规则 $\delta$ 的风险函数定义为

\[
R(\theta,\delta)
=E_\theta L\{\theta,\delta(X)\}.
\]

对于目标 $g(\theta)$ 的点估计，在平方误差损失下，

\[
L(\theta,a)=\|a-g(\theta)\|^2,
\]

因此

\[
R(\theta,\delta)
=E_\theta\|\delta(X)-g(\theta)\|^2.
\]

这里真正比较的对象不是单个数，而是整个风险函数

\[
\theta\longmapsto R(\theta,\delta).
\]

### 2.1 支配与可容许性

!!! info "定义 2.1（支配与可容许性）"

    若估计量 $\delta_1$ 满足

    \[
    R(\theta,\delta_1)
    \leq R(\theta,\delta_0),
    \qquad \forall\theta\in\Theta,
    \]

    且至少对某一个 $\theta$ 严格成立，则称 $\delta_1$ **支配**（dominates）$\delta_0$。

    如果某估计量被另一个估计量支配，则称它是 **不可容许的**（inadmissible）。

    如果不存在任何估计量能够支配它，则称它是 **可容许的**（admissible）。

!!! note "操作意义"

    要证明 $\delta_0$ 不可容许，只需构造一个估计量 $\delta_1$，使其风险处处不大于 $\delta_0$，并且至少某处严格更小。

    相反，要证明可容许性通常更难，因为需要排除所有可能的支配估计量。

---

## 3. 充分统计量带来的决策论改进

Lecture 4 中，Rao–Blackwell 定理通过条件期望降低无偏估计量的方差。这个思想其实更加一般：只要损失函数关于行动是凸的，条件化就可以改善风险，而不要求原估计量无偏。

!!! success "命题 3.1（凸损失下的 Rao–Blackwell 改进）"

    设 $T=T(X)$ 是关于 $\theta$ 的充分统计量，行动空间为凸集，并且对每个固定的 $\theta$，函数

    \[
    a\longmapsto L(\theta,a)
    \]

    是凸函数。对任意具有有限风险的决策规则 $\delta(X)$，定义

    \[
    \delta^*(T)
    =E_\theta\{\delta(X)\mid T\}.
    \]

    由于 $T$ 充分，给定 $T$ 后 $X$ 的条件分布与 $\theta$ 无关，因此 $\delta^*(T)$ 确实是一个统计量，而不是依赖未知参数的表达式。

    此时

    \[
    R(\theta,\delta^*)
    \leq R(\theta,\delta),
    \qquad \forall\theta\in\Theta.
    \]

    如果损失严格凸，并且对某个 $\theta_0$，$\delta(X)$ 不是 $T$ 的几乎处处函数，则在 $\theta_0$ 处严格改进，因此 $\delta$ 不可容许。

??? proof "命题 3.1 的证明（点击展开）"

    由条件 Jensen 不等式，

    \[
    L\!\left(\theta,E_\theta\{\delta(X)\mid T\}\right)
    \leq
    E_\theta\{L(\theta,\delta(X))\mid T\}.
    \]

    两边再取期望，得到

    \[
    R(\theta,\delta^*)
    \leq R(\theta,\delta).
    \]

    如果损失严格凸，并且给定 $T$ 后 $\delta(X)$ 的条件分布在一个正概率集合上非退化，那么 Jensen 不等式严格成立。$\square$

!!! example "例 3.2（忽略充分性的不可容许估计量）"

    设

    \[
    X_1,X_2\overset{\mathrm{iid}}\sim N(\theta,1),
    \]

    在平方误差损失下估计 $\theta$。

    统计量

    \[
    T=\bar X
    \]

    是充分统计量，而估计量 $\delta(X)=X_1$ 不是 $T$ 的函数。并且

    \[
    E_\theta(X_1\mid \bar X)=\bar X.
    \]

    因此 $\bar X$ 支配 $X_1$。事实上，

    \[
    R(\theta,X_1)=1,
    \qquad
    R(\theta,\bar X)=\frac12.
    \]

    这里真正重要的并不是数值 $1$ 和 $1/2$，而是：在严格凸损失下，一个保留了与参数无关的条件噪声的估计量不可能是可容许的。

!!! note "更广泛的意义"

    通过条件化去除辅助随机性或算法随机性这一思想，也会出现在 conditional Monte Carlo、方差缩减和稳定化估计量的构造中。

---

## 4. 为什么 UMVU 不蕴含可容许性

一个估计量即使在所有无偏估计量中最优，也仍可能在所有估计量组成的更大类中被支配。

### 4.1 正态方差估计的经典例子

!!! example "例 4.1（正态总体方差估计）"

    设

    \[
    X_1,\ldots,X_n
    \overset{\mathrm{iid}}\sim N(\mu,\sigma^2),
    \]

    其中 $\mu$ 与 $\sigma^2$ 都未知。通常的无偏方差估计量为

    \[
    S^2
    =\frac{1}{n-1}
    \sum_{i=1}^n(X_i-\bar X)^2.
    \]

    它是 $\sigma^2$ 的 UMVU 估计量。

在估计 $\sigma^2$ 的平方误差损失

\[
L(\sigma^2,a)=(a-\sigma^2)^2
\]

下，有

\[
\frac{(n-1)S^2}{\sigma^2}
\sim\chi^2_{n-1},
\]

以及

\[
\operatorname{Var}(S^2)
=\frac{2\sigma^4}{n-1}.
\]

由于 $S^2$ 无偏，

\[
R(\sigma^2,S^2)
=\frac{2\sigma^4}{n-1}.
\]

现在考虑带偏估计量

\[
cS^2,\qquad c>0.
\]

其风险为

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

对 $c$ 最小化这个二次函数，得到

\[
c^*=\frac{n-1}{n+1}.
\]

因此

\[
\frac{n-1}{n+1}S^2
=\frac{1}{n+1}
\sum_{i=1}^n(X_i-\bar X)^2
\]

对每个 $\sigma^2>0$ 都具有比 $S^2$ 更小的平方误差风险。

!!! warning "结论"

    $S^2$ 虽然是 $\sigma^2$ 的 UMVU 估计量，但在平方误差损失下却是不可容许的。

    这说明：

    \[
    \boxed{\text{UMVU 最优性只是在无偏估计量类中的最优性。}}
    \]

这种改进本质上来自向零方向的收缩：引入一定偏差，以换取更显著的方差下降。

### 4.2 Stein 的第二次改进

!!! example "例 4.2（Stein 1964 的进一步改进）"

    上面得到的估计量

    \[
    \delta_0
    =\frac{1}{n+1}
    \sum_{i=1}^n(X_i-\bar X)^2
    \]

    本身仍然不可容许。Stein (1964) 证明

    \[
    \delta_{\mathrm{St}}
    =\min\left\{
    \frac{1}{n+1}\sum_{i=1}^n(X_i-\bar X)^2,
    \frac{1}{n+2}\sum_{i=1}^nX_i^2
    \right\}
    \]

    对 $\sigma^2$ 的平方误差风险更小。

第二项利用了总体二阶矩的信息：如果 $\sum_iX_i^2$ 本身很小，那么方差估计值也不应相对它大得不合理。

!!! note "教学重点"

    第一层改进，即证明 $S^2$ 被 $\frac{n-1}{n+1}S^2$ 支配，是本章应掌握的核心计算。

    Stein (1964) 的第二次改进主要说明：找到一个更好的估计量，并不意味着已经到达可容许估计量。

!!! note "“向零收缩”的正确理解"

    “shrinkage toward zero” 并不是说真实的 $\sigma^2$ 接近零，而只是说估计量被拉向一个方差更低的目标，从而降低均方误差。代价是偏差，收益是方差下降。

---

## 5. 收缩：偏差—方差权衡

考虑标量目标 $g(\theta)$ 的估计量 $\delta(X)$，以及收缩规则

\[
\delta_c(X)=c\delta(X),
\qquad 0<c<1.
\]

如果 $\delta$ 无偏，则

\[
\operatorname{Bias}_\theta(\delta_c)
=(c-1)g(\theta),
\]

并且

\[
\operatorname{Var}_\theta(\delta_c)
=c^2\operatorname{Var}_\theta(\delta).
\]

因此风险差为

\[
R(\theta,\delta_c)-R(\theta,\delta)
=(c^2-1)\operatorname{Var}_\theta(\delta)
+(c-1)^2g(\theta)^2.
\]

收缩有利的条件是：方差下降带来的收益超过平方偏差带来的损失。

### 5.1 一维正态均值中的固定收缩

!!! example "例 5.1（一维正态均值）"

    设

    \[
    X\sim N(\theta,1),
    \]

    在平方误差损失下估计 $\theta$。

    通常估计量 $X$ 的风险为

    \[
    R(\theta,X)=1.
    \]

    对固定 $0<a<1$，收缩估计量 $aX$ 的风险为

    \[
    R(\theta,aX)
    =a^2+(1-a)^2\theta^2.
    \]

当 $\theta$ 接近零时，$aX$ 的确可能优于 $X$；但当 $|\theta|$ 很大时，偏差项 $(1-a)^2\theta^2$ 会占主导。

因此，没有任何固定的 $a<1$ 可以在一维中处处支配 $X$。

!!! note "启示"

    要实现统一风险改进，收缩程度通常必须根据观测数据自适应，而不能使用固定的收缩比例。

### 5.2 一维正态模型中的仿射规则

考虑

\[
X\sim N(\theta,\sigma^2),
\]

其中 $\sigma^2$ 已知，在平方误差下估计 $\theta\in\mathbb R$。令

\[
\delta_{a,\beta}(X)=aX+\beta.
\]

!!! success "命题 5.2（一维正态模型中仿射规则的可容许性分类）"

    对 $\delta_{a,\beta}(X)=aX+\beta$：

    | 参数范围 | 可容许性 |
    |---|---|
    | $a>1$、$a<0$，或 $a=1,\beta\neq0$ | 不可容许 |
    | $a=0$、$0<a<1$，或 $a=1,\beta=0$ | 可容许 |

其风险为

\[
R(\theta,aX+\beta)
=a^2\sigma^2+\bigl\{(a-1)\theta+\beta\bigr\}^2.
\]

例如当 $a=1$ 且 $\beta\neq0$ 时，$X+\beta$ 被 $X$ 支配。

对于 $0<a<1$，这个仿射规则可以写成某个 proper normal prior 下的后验均值。具体地，若

\[
m=\frac{\beta}{1-a},
\qquad
\tau^2=\frac{a}{1-a}\sigma^2,
\]

则 $aX+\beta$ 对应向先验中心 $m$ 的收缩。

!!! note "收缩与扩张"

    $0<a<1$ 的收缩往往可以被 Bayes 结构解释；而 $a>1$ 的扩张则很难在所有参数值上统一合理。

---

## 6. 等变性与不变风险

许多统计问题存在对称性。此时，一个自然的思想是要求估计程序尊重这种对称性。

### 6.1 位置等变性

设

\[
X\sim N(\theta,\sigma^2).
\]

如果数据整体平移 $a$，参数也相应平移 $a$。一个位置等变估计量满足

\[
\delta(x+a)=\delta(x)+a.
\]

通常估计量

\[
\delta(x)=x
\]

就是位置等变的。

!!! info "定义 6.1（等变估计量）"

    设群 $G$ 同时作用于样本空间和参数空间。若

    \[
    \delta(gx)=g\delta(x),
    \qquad g\in G,
    \]

    则称 $\delta$ 是等变的（equivariant）。

### 6.2 不变损失与风险

如果损失满足

\[
L(g\theta,ga)=L(\theta,a),
\qquad g\in G,
\]

则称损失函数是不变的。

当模型和损失都具有相同对称性时，等变估计量的风险通常在群轨道上保持常数，从而显著简化决策问题。

!!! example "例 6.2（正态位置模型）"

    若

    \[
    X\sim N(\theta,\sigma^2),
    \qquad
    L(\theta,a)=(a-\theta)^2,
    \]

    则对位置等变估计量 $\delta(X)=X$，

    \[
    R(\theta,X)=\sigma^2,
    \]

    与 $\theta$ 无关。

!!! warning "等变性不等于可容许性"

    一个估计量可能在等变规则类中最优，却在所有估计量组成的更大类中不可容许。

    James–Stein 现象正是这一点的经典例子。

---

## 7. 多元正态均值问题

考虑

\[
X\sim N_p(\theta,\sigma^2I_p),
\qquad
\theta\in\mathbb R^p,
\]

其中 $\sigma^2$ 已知，损失为总平方误差

\[
L(\theta,a)=\|a-\theta\|^2.
\]

通常估计量为

\[
\delta_0(X)=X.
\]

它具有很多非常强的性质：

1. 每个坐标 $X_j$ 都是 $\theta_j$ 的 UMVU 估计量；
2. 在无偏向量估计量中具有最小总风险；
3. 它是 MLE；
4. 它是 minimax；
5. 它是平移等变的，并且是最佳等变规则。

其风险为

\[
R(\theta,X)
=E_\theta\|X-\theta\|^2
=p\sigma^2.
\]

!!! note "令人惊讶之处"

    即使一个估计量同时是 coordinatewise UMVU、MLE、best equivariant 和 minimax，它仍可能不是可容许的。

    当 $p\geq3$ 时，James–Stein 估计量就会统一改进 $X$。

---

## 8. James–Stein 估计量

当 $p\geq3$ 时，定义

\[
\delta_c(X)
=
\left(
1-
\frac{c\sigma^2(p-2)}{\|X\|^2}
\right)X,
\qquad 0<c<2.
\]

这个估计量把 $X$ 向原点收缩，其收缩因子为

\[
1-\frac{c\sigma^2(p-2)}{\|X\|^2}.
\]

当 $\|X\|^2$ 很大时，收缩很弱；当 $\|X\|^2$ 较小时，收缩更强。因此它与前面的一维固定收缩不同，是 **数据自适应收缩**。

!!! success "定理 8.1（James–Stein 现象）"

    设

    \[
    X\sim N_p(\theta,\sigma^2I_p),
    \qquad p\geq3.
    \]

    在平方误差损失下，对任意

    \[
    0<c<2,
    \]

    都有

    \[
    R(\theta,\delta_c)
    <R(\theta,X)
    =p\sigma^2,
    \qquad \forall\theta\in\mathbb R^p.
    \]

    特别地，当 $c=1$ 时得到标准 James–Stein 估计量

    \[
    \boxed{
    \delta_{\mathrm{JS}}(X)
    =
    \left(
    1-\frac{\sigma^2(p-2)}{\|X\|^2}
    \right)X
    }.
    \]

    因此，当 $p\geq3$ 时，通常估计量 $X$ 是不可容许的。

!!! warning "维数阈值"

    当 $p=1$ 或 $p=2$ 时，通常估计量 $X$ 是可容许的。

    James–Stein 不可容许性只从

    \[
    p\geq3
    \]

    开始出现。

---

## 9. Stein 恒等式与 SURE 风险计算

James–Stein 结果可以通过 Stein 恒等式非常简洁地理解。

设估计量写成

\[
\delta(X)=X+g(X),
\]

其中

\[
g:\mathbb R^p\to\mathbb R^p
\]

弱可微，并满足适当可积条件。

!!! success "定理 9.1（Stein 恒等式）"

    对

    \[
    X\sim N_p(\theta,\sigma^2I_p),
    \]

    有

    \[
    E_\theta\{(X-\theta)^Tg(X)\}
    =\sigma^2E_\theta\{\nabla\cdot g(X)\},
    \]

    其中散度定义为

    \[
    \nabla\cdot g(x)
    =\sum_{j=1}^p
    \frac{\partial g_j(x)}{\partial x_j}.
    \]

于是

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

因此

\[
\boxed{
\operatorname{SURE}(X+g)
=p\sigma^2
+\|g(X)\|^2
+2\sigma^2\nabla\cdot g(X)
}
\]

是风险的无偏估计，这就是 Stein's Unbiased Risk Estimate (SURE)。

### 9.1 对 James–Stein 形式应用 SURE

取

\[
g(x)=-\frac{a}{\|x\|^2}x,
\qquad a>0.
\]

首先，

\[
\|g(x)\|^2
=\frac{a^2}{\|x\|^2}.
\]

其次，

\[
\nabla\cdot g(x)
=-a\nabla\cdot\left(\frac{x}{\|x\|^2}\right).
\]

直接计算得到

\[
\nabla\cdot\left(\frac{x}{\|x\|^2}\right)
=\frac{p-2}{\|x\|^2}.
\]

因此

\[
R\left(\theta,X-\frac{aX}{\|X\|^2}\right)
-p\sigma^2
=
E_\theta\left[
\frac{a^2-2a\sigma^2(p-2)}{\|X\|^2}
\right].
\]

关于 $a$ 最小化括号中的二次式，得到

\[
a^*=\sigma^2(p-2).
\]

令

\[
a=c\sigma^2(p-2),
\]

则

\[
\boxed{
R(\theta,\delta_c)-R(\theta,X)
=-c(2-c)\sigma^4(p-2)^2
E_\theta\left(\frac{1}{\|X\|^2}\right)
}.
\]

当 $0<c<2$ 时，右侧严格小于零，于是得到 James–Stein 定理。

!!! note "为什么临界维数是 3？"

    关键就在散度公式中的

    \[
    p-2.
    \]

    当 $p\geq3$ 时，这个量为正，从而可以通过适当收缩获得严格风险改进。

---

## 10. Positive-part James–Stein 估计量

原始 James–Stein 的收缩因子

\[
1-\frac{\sigma^2(p-2)}{\|X\|^2}
\]

可能为负。当

\[
\|X\|^2<\sigma^2(p-2)
\]

时，估计量甚至会把 $X$ 的方向反转，这在直觉上并不理想。

因此定义 positive-part James–Stein 估计量

\[
\boxed{
\delta_{\mathrm{JS}+}(X)
=
\left(
1-\frac{\sigma^2(p-2)}{\|X\|^2}
\right)_+X
}
\]

其中

\[
u_+=\max(u,0).
\]

!!! success "重要事实"

    positive-part James–Stein 估计量支配原始 James–Stein 估计量。

    因此，原始 James–Stein 估计量本身也是不可容许的。

但是，positive-part 版本也并非完全可容许。这说明可容许性理论比“找到一次风险改进”更加微妙。

---

## 11. James–Stein 的解释与收缩中心

James–Stein 现象并不是“魔法”。其核心思想是：同时估计许多带噪坐标时，可以利用整个坐标集合来判断信号规模，并据此决定收缩强度。

### 11.1 层次模型与经验 Bayes 解释

考虑层次模型

\[
Y_i\mid\theta_i
\sim N(\theta_i,\sigma_0^2),
\]

以及

\[
\theta_i\sim N(\mu,\tau^2).
\]

在这个模型下，$\theta_i$ 的后验均值是 $Y_i$ 与 $\mu$ 的加权平均，即自然地产生向 $\mu$ 收缩的估计规则。

如果再利用整个样本估计未知的收缩强度，就会得到 James–Stein 类型的经验 Bayes 规则。

### 11.2 收缩中心不一定是零

零只是一个方便的收缩中心。对于固定的、具有科学意义的

\[
m\in\mathbb R^p,
\]

可以定义

\[
\delta_m(X)
=m+
\left(
1-\frac{\sigma^2(p-2)}{\|X-m\|^2}
\right)(X-m).
\]

它把 $X$ 向 $m$ 收缩。

!!! warning "随机收缩中心需要重新计算"

    如果 $m$ 不是固定常数，而是由数据估计得到，例如取总体坐标平均，那么有效维数与风险推导都会发生变化。

    不能简单地把固定中心公式中的 $m$ 换成随机估计量而不重新检查自由度与风险计算。

---

## 12. 指数族中的仿射收缩规则

考虑自然指数族

\[
dP_\eta(x)
=C(\eta)\exp\{\eta T(x)\}\,d\mu(x),
\qquad
\eta\in(\underline\eta,\overline\eta),
\]

并在平方误差损失下估计

\[
g(\eta)=E_\eta T.
\]

考虑仿射收缩规则

\[
\delta_{k,\lambda}(T)
=\frac{T+k\lambda}{1+\lambda},
\qquad \lambda\geq0.
\]

该估计量把 $T$ 向 $k$ 收缩。

!!! success "定理 12.1（Karlin 可容许性判据）"

    若存在内部点 $\eta_0$，使得

    \[
    \int_{\underline\eta}^{\eta_0}
    C(\eta)^{-\lambda}e^{-k\lambda\eta}\,d\eta
    =\infty
    \]

    且

    \[
    \int_{\eta_0}^{\overline\eta}
    C(\eta)^{-\lambda}e^{-k\lambda\eta}\,d\eta
    =\infty,
    \]

    则 $\delta_{k,\lambda}$ 对 $g(\eta)$ 的估计是可容许的。

当 $\lambda=0$ 时，在标准的 full exponential family 中，该结果包含未收缩规则 $T$。

!!! note "课程层面的理解"

    这个定理的长证明不是本课程重点。需要掌握的是：指数族中的可容许仿射收缩规则，与 conjugate Bayes 或 generalized Bayes 结构存在紧密联系。

---

## 13. 为什么 James–Stein 结果令人惊讶

James–Stein 现象有四个特别值得记住的地方。

### 13.1 坐标不需要具有科学关联

即使 $p$ 个均值代表完全不同的科学量，只要损失函数是总平方误差

\[
\sum_{j=1}^p(\hat\theta_j-\theta_j)^2,
\]

联合估计就允许在坐标之间进行信息共享。

### 13.2 每个坐标可以有偏，但总风险更低

James–Stein 估计量允许每个坐标产生一定偏差，但这种偏差换来了整体更明显的方差下降。

### 13.3 一个“非常优秀”的估计量仍可能不可容许

通常估计量 $X$ 同时具有

- 无偏性；
- MLE 性质；
- 等变性；
- minimax 性质；
- coordinatewise UMVU 性质。

但是，当 $p\geq3$ 时，它仍然被 James–Stein 估计量支配。

### 13.4 维数阈值是真实存在的

风险改进来自散度计算中的

\[
p-2,
\]

因此临界维数恰好出现在 $p=3$。

---

## 14. 有限样本最优性标准的比较

| 标准 | 含义 | 局限 |
|---|---|---|
| 无偏性 | 对每个参数值平均意义下正确 | 可能具有不必要的大方差 |
| UMVU | 在无偏估计量中方差最小 | 在所有估计量中仍可能不可容许 |
| 可容许性 | 不存在统一支配它的估计量 | 条件较弱；许多可容许规则彼此无法比较 |
| Minimax | 最小化最坏情形风险 | 可能保守，也不蕴含可容许性 |
| 等变性 | 尊重统计问题的对称结构 | 最佳等变规则在全体规则中可能不可容许 |
| Bayes 最优性 | 最小化先验下的平均风险 | 依赖先验；improper prior 需要谨慎处理 |

---

## 15. 本章小结

本章的主要结论如下：

1. 决策论用风险函数 $R(\theta,\delta)$ 比较估计规则，而不是只比较无偏估计量的方差。
2. 一个估计量若被另一个估计量在所有参数值上风险不增、并在某处严格降低，则它是不可容许的。
3. 对凸损失函数，基于充分统计量做 Rao–Blackwell 条件化可以降低风险；严格凸时通常能够严格改进非充分统计量函数的规则。
4. UMVU 不等于全局风险最优。正态方差估计中的 $S^2$ 虽为 UMVU，却被适当的带偏收缩估计量支配。
5. 收缩估计通过“增加少量偏差、显著降低方差”实现均方误差改进；固定收缩在一维中无法统一支配通常估计量。
6. 等变性利用模型对称性简化风险，但最佳等变估计量仍可能在更大的决策类中不可容许。
7. 对 $X\sim N_p(\theta,\sigma^2I_p)$，当 $p\geq3$ 时，James–Stein 估计量统一支配通常估计量 $X$。
8. Stein 恒等式与 SURE 给出了 James–Stein 风险改进的简洁推导，关键量是散度中的 $p-2$。
9. positive-part James–Stein 进一步支配原始 James–Stein，说明风险改进并不自动意味着已经达到可容许性。
10. 收缩思想与层次 Bayes、经验 Bayes 以及指数族中的广义 Bayes 结构有深刻联系。
