---
tags:
  - 高等数理统计
  - 可容许性
  - 收缩估计
  - Bayes 估计
  - James-Stein 估计
  - 课后习题
---

# 高等数理统计第三次作业


## 习题 1：线性收缩估计量的风险

**题目。** 设

\[
X\sim N(\theta,1).
\]

在平方误差损失下，计算 $aX$ 估计 $\theta$ 的风险，并证明不存在固定的 $a<1$ 能够在所有 $\theta\in\mathbb R$ 上优于 $X$。

??? success "解答（点击展开）"

    因为

    \[
    E_\theta(X)=\theta,
    \qquad
    \operatorname{Var}_\theta(X)=1,
    \]

    所以

    \[
    E_\theta(X^2)=1+\theta^2.
    \]

    记 $\delta_a(X)=aX$。其平方误差风险为

    \[
    R(\theta,\delta_a)
    =E_\theta[(aX-\theta)^2].
    \]

    展开可得

    \[
    R(\theta,\delta_a)
    =a^2E_\theta(X^2)-2a\theta E_\theta(X)+\theta^2.
    \]

    因此

    \[
    \boxed{
    R(\theta,\delta_a)
    =a^2+(a-1)^2\theta^2
    }.
    \]

    对无偏估计量 $\delta_0(X)=X$，

    \[
    R(\theta,\delta_0)
    =E_\theta[(X-\theta)^2]
    =1.
    \]

    两个风险之差为

    \[
    R(\theta,\delta_a)-R(\theta,\delta_0)
    =a^2-1+(a-1)^2\theta^2.
    \]

    对任意固定的 $a<1$，有 $(a-1)^2>0$，所以

    \[
    \lim_{|\theta|\to\infty}
    \{R(\theta,\delta_a)-R(\theta,\delta_0)\}
    =+\infty.
    \]

    因而总存在某个 $\theta\in\mathbb R$，使得

    \[
    R(\theta,\delta_a)>R(\theta,\delta_0).
    \]

    所以，不存在固定的 $a<1$ 能够在整个参数空间上优于 $X$。

---

## 习题 2：投影到参数空间与不可容许性

**题目。** 设参数 $\vartheta\in[c,d]$，损失函数为 $L(|\vartheta-a|)$，其中 $L$ 单调递增。证明：任何以正概率取到区间 $[c,d]$ 之外数值的决策规则都是不可容许的。

??? success "解答（点击展开）"

    设 $\delta(X)$ 是一个决策规则。定义到参数区间的投影

    \[
    g(a)=
    \begin{cases}
    c, & a<c,\\
    a, & c\leq a\leq d,\\
    d, & a>d,
    \end{cases}
    \]

    并令

    \[
    \delta^*(X)=g(\delta(X)).
    \]

    对每个固定的 $\vartheta\in[c,d]$，比较行动 $a$ 与 $g(a)$ 的损失。

    当 $a<c$ 时，$\vartheta\geq c>a$，故

    \[
    |\vartheta-g(a)|
    =|\vartheta-c|
    =\vartheta-c
    <\vartheta-a
    =|\vartheta-a|.
    \]

    当 $a>d$ 时，$\vartheta\leq d<a$，故

    \[
    |\vartheta-g(a)|
    =|\vartheta-d|
    =d-\vartheta
    <a-\vartheta
    =|\vartheta-a|.
    \]

    当 $c\leq a\leq d$ 时，$g(a)=a$，两个距离相等。因此

    \[
    |\vartheta-\delta^*(X)|
    \leq|\vartheta-\delta(X)|
    \qquad\text{对每个 }\vartheta\in[c,d].
    \]

    由于 $L$ 单调递增，

    \[
    L(|\vartheta-\delta^*(X)|)
    \leq L(|\vartheta-\delta(X)|).
    \]

    取期望可得

    \[
    R(\vartheta,\delta^*)
    \leq R(\vartheta,\delta)
    \qquad\text{对所有 }\vartheta\in[c,d].
    \]

    如果 $L$ 严格递增，并且对某个 $\vartheta_0$，

    \[
    P_{\vartheta_0}\{\delta(X)\notin[c,d]\}>0,
    \]

    那么在一个正概率事件上损失不等式严格成立，从而

    \[
    R(\vartheta_0,\delta^*)
    <R(\vartheta_0,\delta).
    \]

    因此 $\delta^*$ 支配 $\delta$，所以 $\delta$ 不可容许。

---

## 习题 3：位置分布族中的 Pitman 估计量

**题目。** 设 $X_1,\ldots,X_n$ 是来自位置分布族的随机样本，其密度为 $f(x-\theta)$。Pitman 估计量定义为

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

(a) 证明 $d_P$ 具有位置等变性。

(b) 当 $f$ 为标准正态密度时，求 $d_P$。

(c) 当 $f$ 为 $(-1/2,1/2)$ 上的均匀密度时，求 $d_P$。

??? success "解答（点击展开）"

    **(a) 位置等变性**

    对任意 $b\in\mathbb R$，记

    \[
    X+b\mathbf 1=(X_1+b,\ldots,X_n+b).
    \]

    则

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

    作变量替换 $u=t-b$，则 $t=u+b$ 且 $dt=du$。分子变为

    \[
    \int_{-\infty}^{\infty}
    (u+b)\prod_{i=1}^nf(X_i-u)\,du.
    \]

    将其拆开，

    \[
    \int_{-\infty}^{\infty}
    u\prod_{i=1}^nf(X_i-u)\,du
    +b\int_{-\infty}^{\infty}
    \prod_{i=1}^nf(X_i-u)\,du.
    \]

    分母变为

    \[
    \int_{-\infty}^{\infty}
    \prod_{i=1}^nf(X_i-u)\,du.
    \]

    所以，只要定义中的积分存在，

    \[
    \boxed{d_P(X+b\mathbf 1)=d_P(X)+b},
    \]

    即 $d_P$ 具有位置等变性。

    **(b) 标准正态情形**

    当

    \[
    f(z)=\frac1{\sqrt{2\pi}}e^{-z^2/2}
    \]

    时，

    \[
    \prod_{i=1}^nf(X_i-t)
    =(2\pi)^{-n/2}
    \exp\left\{-\frac12\sum_{i=1}^n(X_i-t)^2\right\}.
    \]

    配方得到

    \[
    \sum_{i=1}^n(X_i-t)^2
    =n(t-\overline X)^2
    +\sum_{i=1}^n(X_i-\overline X)^2.
    \]

    与 $t$ 无关的因子在比值中消去，所以

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

    该比值是关于 $t$ 的 $N(\overline X,1/n)$ 密度的均值，因此

    \[
    \boxed{d_P(X)=\overline X}.
    \]

    **(c) 均匀分布情形**

    当

    \[
    f(z)=\mathbf 1_{\{-1/2<z<1/2\}}
    \]

    时，

    \[
    \prod_{i=1}^nf(X_i-t)
    =\prod_{i=1}^n
    \mathbf 1_{\{X_i-1/2<t<X_i+1/2\}}.
    \]

    所以

    \[
    \prod_{i=1}^nf(X_i-t)
    =\mathbf 1_{\{X_{(n)}-1/2<t<X_{(1)}+1/2\}}.
    \]

    记

    \[
    a=X_{(n)}-\frac12,
    \qquad
    b=X_{(1)}+\frac12.
    \]

    在该模型下，几乎必然有 $a<b$。因此

    \[
    d_P(X)
    =\frac{\int_a^bt\,dt}{\int_a^b1\,dt}
    =\frac{a+b}{2}.
    \]

    故

    \[
    \boxed{
    d_P(X)=\frac{X_{(1)}+X_{(n)}}2
    },
    \]

    即样本的中程数。

---

## 习题 4：非正常先验与广义 Bayes 行动

**题目。** 设

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2),
\]

并采用非正常先验

\[
\pi(\mu,\sigma^2)
=\sigma^{-2}\mathbf 1_{\{\sigma^2>0\}}.
\]

证明 $\mu$ 的边际后验分布是以 $\overline X$ 为中心的 Student $t$ 分布，并求在平方误差损失下估计 $\mu/\sigma$ 的广义 Bayes 行动。

??? success "解答（点击展开）"

    令

    \[
    v=\sigma^2,
    \qquad
    S=\sum_{i=1}^n(X_i-\overline X)^2=(n-1)s^2,
    \]

    其中 $s^2$ 是通常的样本方差。假设 $n\geq2$ 且 $S>0$；在非退化正态模型下，这些条件几乎必然成立。

    利用似然函数和先验 $\pi(\mu,v)\propto v^{-1}$，联合后验核为

    \[
    \pi(\mu,v\mid X)
    \propto
    v^{-(n/2+1)}
    \exp\left\{-\frac1{2v}
    \sum_{i=1}^n(X_i-\mu)^2\right\}.
    \]

    由恒等式

    \[
    \sum_{i=1}^n(X_i-\mu)^2
    =S+n(\mu-\overline X)^2
    \]

    可得

    \[
    \pi(\mu,v\mid X)
    \propto
    v^{-(n/2+1)}
    \exp\left\{-\frac{S+n(\mu-\overline X)^2}{2v}\right\}.
    \]

    对 $v$ 积分，得到 $\mu$ 的边际后验核

    \[
    \pi(\mu\mid X)
    \propto
    \int_0^\infty
    v^{-(n/2+1)}
    \exp\left\{-\frac{S+n(\mu-\overline X)^2}{2v}\right\}dv.
    \]

    令

    \[
    A=S+n(\mu-\overline X)^2
    \]

    并作变量替换 $w=A/(2v)$，可得

    \[
    \pi(\mu\mid X)
    \propto A^{-n/2}.
    \]

    因而

    \[
    \pi(\mu\mid X)
    \propto
    \left\{1+
    \frac{n(\mu-\overline X)^2}{(n-1)s^2}
    \right\}^{-n/2}.
    \]

    Student $t_\nu(m,a)$ 的位置-尺度密度核为

    \[
    \left\{1+\frac1\nu
    \left(\frac{\mu-m}{a}\right)^2
    \right\}^{-(\nu+1)/2}.
    \]

    比较可知

    \[
    \nu=n-1,
    \qquad
    m=\overline X,
    \qquad
    a=\frac{s}{\sqrt n}.
    \]

    因此

    \[
    \boxed{
    \left.
    \frac{\mu-\overline X}{s/\sqrt n}
    \right|X
    \sim t_{n-1}
    }.
    \]

    在平方误差损失下，估计

    \[
    \eta=\frac\mu\sigma=\frac\mu{\sqrt v}
    \]

    的广义 Bayes 行动是后验均值

    \[
    \delta_{\mathrm{GB}}(X)
    =E\left(\left.\frac\mu{\sqrt v}\right|X\right).
    \]

    由联合后验分布，

    \[
    \mu\mid v,X
    \sim N\left(\overline X,\frac vn\right),
    \]

    所以

    \[
    E(\mu\mid v,X)=\overline X.
    \]

    将条件正态密度对 $\mu$ 积分，还可得到

    \[
    \pi(v\mid X)
    \propto
    v^{-\{(n-1)/2+1\}}
    \exp\left(-\frac S{2v}\right).
    \]

    因此

    \[
    v\mid X\sim\operatorname{Inv-Gamma}(\alpha,\beta),
    \]

    其中

    \[
    \alpha=\frac{n-1}{2},
    \qquad
    \beta=\frac S2.
    \]

    这里采用的参数化为

    \[
    p(v)=\frac{\beta^\alpha}{\Gamma(\alpha)}
    v^{-\alpha-1}e^{-\beta/v},
    \qquad v>0.
    \]

    由迭代期望，

    \[
    \delta_{\mathrm{GB}}(X)
    =\overline X E(v^{-1/2}\mid X).
    \]

    对逆 Gamma 分布计算负半阶矩，得到

    \[
    E(v^{-1/2}\mid X)
    =\beta^{-1/2}
    \frac{\Gamma(\alpha+1/2)}{\Gamma(\alpha)}.
    \]

    代入 $\alpha=(n-1)/2$、$\beta=S/2$，

    \[
    E(v^{-1/2}\mid X)
    =\sqrt{\frac2S}
    \frac{\Gamma(n/2)}{\Gamma((n-1)/2)}.
    \]

    最终

    \[
    \boxed{
    \delta_{\mathrm{GB}}(X)
    =\overline X
    \sqrt{\frac{2}{(n-1)s^2}}
    \frac{\Gamma(n/2)}{\Gamma((n-1)/2)}
    }.
    \]

    当 $n\geq2$ 且 $S>0$ 时，后验分布正常化且上述后验矩有限，所以该广义 Bayes 行动定义良好。

---

## 习题 5：修正的 James-Stein 估计量

**题目。** 设

\[
X\sim N_d(\theta,I_d),
\qquad
\theta\in\mathbb R^d,
\]

损失函数为平方误差损失。记

\[
\mathbf 1=(1,\ldots,1)^\top,
\qquad
\overline x=\frac{\langle x,\mathbf 1\rangle}{d}.
\]

考虑修正的 James-Stein 估计量

\[
\widehat\theta(x)
=\overline x\mathbf 1
+\left\{1-h(\|x-\overline x\mathbf 1\|^2)\right\}
(x-\overline x\mathbf 1).
\]

等价地，

\[
\widehat\theta(x)
=x-h(\|x-\overline x\mathbf 1\|^2)
(x-\overline x\mathbf 1).
\]

(a) 利用 SURE 给出该估计量的风险表达式，并假设应用 SURE 所需的正则性条件均成立。

(b) 对 $h(u)=C/u$ 化简该表达式，其中 $C$ 为常数。

(c) 当 $d\geq4$ 时，选择 $C$ 使 $\widehat\theta$ 支配无偏估计量 $\widehat\theta_u(x)=x$，并解释 $d=3$ 时的情况。

(d) 给出 $\widehat\theta$ 的经验 Bayes 解释。

??? success "解答（点击展开）"

    令

    \[
    z=x-\overline x\mathbf 1,
    \qquad
    u=\|z\|^2=\sum_{i=1}^dz_i^2.
    \]

    **(a) SURE 风险表达式**

    将估计量写成

    \[
    \widehat\theta(x)=x+g(x),
    \qquad
    g(x)=-h(u)z.
    \]

    对 $X\sim N_d(\theta,I_d)$，Stein 无偏风险恒等式给出

    \[
    R(\theta,\widehat\theta)
    =d+E_\theta\left[
    \|g(X)\|^2+2\operatorname{div}g(X)
    \right].
    \]

    首先，

    \[
    \|g(x)\|^2
    =h(u)^2\|z\|^2
    =uh(u)^2.
    \]

    因为

    \[
    z_i=x_i-\overline x,
    \qquad
    \overline x=\frac1d\sum_{j=1}^dx_j,
    \]

    所以

    \[
    \frac{\partial z_i}{\partial x_i}=1-\frac1d,
    \qquad
    \sum_{i=1}^dz_i=0.
    \]

    此外，

    \[
    \frac{\partial u}{\partial x_i}
    =2\sum_{j=1}^dz_j
    \left(\mathbf 1_{\{i=j\}}-\frac1d\right)
    =2z_i.
    \]

    因此，由乘积法则和链式法则，

    \[
    \operatorname{div}g(x)
    =\sum_{i=1}^d
    \frac\partial{\partial x_i}[-h(u)z_i].
    \]

    进一步，

    \[
    \operatorname{div}g(x)
    =-2h'(u)\sum_{i=1}^dz_i^2
    -h(u)\sum_{i=1}^d\left(1-\frac1d\right).
    \]

    所以

    \[
    \operatorname{div}g(x)
    =-2uh'(u)-(d-1)h(u).
    \]

    代入 SURE，得到

    \[
    \boxed{
    R(\theta,\widehat\theta)
    =d+E_\theta\left[
    uh(u)^2-4uh'(u)-2(d-1)h(u)
    \right]
    },
    \]

    其中

    \[
    u=\|X-\overline X\mathbf 1\|^2.
    \]

    **(b) 取 $h(u)=C/u$**

    此时

    \[
    h'(u)=-\frac C{u^2}.
    \]

    所以

    \[
    uh(u)^2=\frac{C^2}{u},
    \qquad
    -4uh'(u)=\frac{4C}{u},
    \qquad
    -2(d-1)h(u)=-\frac{2(d-1)C}{u}.
    \]

    因此

    \[
    \boxed{
    R(\theta,\widehat\theta)
    =d+C\{C-2(d-3)\}
    E_\theta\left[
    \frac1{\|X-\overline X\mathbf 1\|^2}
    \right]
    }.
    \]

    该表达式要求相应期望有限且满足 SURE 的可积性条件；特别地，在 $d\geq4$ 时这些条件成立。

    **(c) 支配条件与 $d=3$**

    无偏估计量 $\widehat\theta_u(X)=X$ 的风险恒为

    \[
    R(\theta,\widehat\theta_u)
    =E_\theta\|X-\theta\|^2
    =d.
    \]

    当 $d\geq4$ 时，残差向量

    \[
    Z=X-\overline X\mathbf 1
    \]

    位于维数为 $d-1\geq3$ 的子空间中。在该子空间的一组标准正交基下，

    \[
    u=\|Z\|^2
    \sim\chi_{d-1}^{\prime 2}(\lambda),
    \qquad
    \lambda=\|\theta-\overline\theta\mathbf 1\|^2.
    \]

    自由度大于 $2$ 的非中心卡方变量满足

    \[
    0<E_\theta(1/u)<\infty.
    \]

    所以风险差为

    \[
    R(\theta,\widehat\theta)
    -R(\theta,\widehat\theta_u)
    =C\{C-2(d-3)\}E_\theta(1/u).
    \]

    选择

    \[
    0<C<2(d-3)
    \]

    即可使风险差对每个 $\theta$ 都严格为负。二次因子在

    \[
    C=d-3
    \]

    处达到最小值，从而得到

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

    因此，该修正的 James-Stein 估计量在总平方误差风险下严格支配 $X$。

    当 $d=3$ 时，残差空间只有 $2$ 维。这时

    \[
    u\sim\chi_2^{\prime 2}(\lambda)
    \]

    且 $E_\theta(1/u)=\infty$。对任意 $C\neq0$，

    \[
    E_\theta\|g(X)\|^2
    =C^2E_\theta(1/u)
    =\infty.
    \]

    因此，(b) 中的形式 SURE 表达式不能用于 $d=3$；它所要求的正则性和可积性条件不成立。在这个特定的 $C/u$ 估计量族中，唯一具有有限风险的选择是 $C=0$，这时 $\widehat\theta=X$，没有严格改进。

    **(d) 经验 Bayes 解释**

    考虑正态-正态层次模型

    \[
    X_i\mid\theta_i\sim N(\theta_i,1),
    \qquad
    \theta_i\mid\mu,\tau^2\sim N(\mu,\tau^2),
    \qquad i=1,\ldots,d.
    \]

    这里，即使在计算频率风险时 $\theta_i$ 被视为固定未知参数，我们仍给它们指定一个工作先验。由 Bayes 公式，后验均值为

    \[
    E(\theta_i\mid X_i,\mu,\tau^2)
    =\frac{\tau^2X_i+\mu}{1+\tau^2}.
    \]

    等价地，

    \[
    E(\theta_i\mid X_i,\mu,\tau^2)
    =\mu+left(1-\frac1{1+\tau^2}\right)(X_i-\mu).
    \]

    边际上，

    \[
    X_i\sim N(\mu,1+\tau^2).
    \]

    用数据估计先验中心 $\widehat\mu=\overline X$，并记

    \[
    u=\sum_{i=1}^d(X_i-\overline X)^2.
    \]

    取边际方差的插件估计

    \[
    1+\widehat\tau^2=\frac{u}{d-3},
    \qquad d\geq4.
    \]

    则

    \[
    \widehat\theta_i^{\mathrm{EB}}
    =\overline X
    +\left(1-\frac{d-3}{u}\right)(X_i-\overline X).
    \]

    写成向量形式，

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

    这正是 $C=d-3$ 时的修正 James-Stein 估计量。
