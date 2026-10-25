---
tags:
  - 高等数理统计
  - 完备性
  - UMVUE
  - U-统计量
  - 信息不等式
  - 课后习题
---

# 高等数理统计第二次作业


## 习题 1：二项分布中的 UMVUE

**题目。** 设

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(p),
\qquad p\in(0,1),
\]

并令

\[
T=\sum_{i=1}^nX_i.
\]

(i) 当 $m\leq n$ 时，求 $p^m$ 的 UMVUE。

(ii) 当 $k\leq m\leq n$ 时，求 $P(X_1+\cdots+X_m=k)$ 的 UMVUE。

??? success "解答（点击展开）"

    联合概率质量函数为

    \[
    f_p(x_1,\ldots,x_n)
    =p^{\sum_{i=1}^n x_i}(1-p)^{n-\sum_{i=1}^n x_i}.
    \]

    由因子分解定理，$T=\sum_{i=1}^nX_i$ 是 $p$ 的充分统计量。又因为

    \[
    T\sim\operatorname{Binomial}(n,p)
    \]

    二项分布族是满的一参数指数族，所以 $T$ 完备。因此，由 Lehmann-Scheffé 定理，任何关于 $T$ 的无偏估计量都是相应参数函数的 UMVUE。

    **(i) $p^m$ 的 UMVUE**

    令

    \[
    W=X_1X_2\cdots X_m.
    \]

    则

    \[
    E_p(W)=p^m.
    \]

    对 $W$ 关于 $T$ 作 Rao-Blackwell 化。对 $t=0,\ldots,n$，

    \[
    E_p(W\mid T=t)
    =P(X_1=\cdots=X_m=1\mid T=t).
    \]

    条件于总成功次数为 $t$，成功位置在 $n$ 个位置中等可能，因此

    \[
    E_p(W\mid T=t)
    =\frac{\binom{n-m}{t-m}}{\binom nt}
    =\frac{\binom tm}{\binom nm}.
    \]

    约定组合数的下标不在允许范围内时其值为零。故

    \[
    \boxed{
    \widehat{p^m}_{\mathrm{UMVUE}}
    =\frac{\binom Tm}{\binom nm}
    }.
    \]

    **(ii) $P(X_1+\cdots+X_m=k)$ 的 UMVUE**

    令

    \[
    W=\mathbf 1_{\{X_1+\cdots+X_m=k\}}.
    \]

    则

    \[
    E_p(W)
    =P(X_1+\cdots+X_m=k)
    =\binom mkp^k(1-p)^{m-k}.
    \]

    条件于 $T=t$，前 $m$ 个观测中的成功次数服从超几何分布，所以

    \[
    E_p(W\mid T=t)
    =\frac{\binom mk\binom{n-m}{t-k}}{\binom nt}.
    \]

    因此

    \[
    \boxed{
    \widehat{P(X_1+\cdots+X_m=k)}_{\mathrm{UMVUE}}
    =\frac{\binom mk\binom{n-m}{T-k}}{\binom nT}
    }.
    \]

---

## 习题 2：变换参数的 Poisson UMVUE

**题目。** 设

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Pois}(\theta),
\]

令

\[
S=\sum_{i=1}^nX_i,
\]

并固定 $a>0$。求

\[
g(\theta)=e^{-a\theta}
\]

的 UMVUE。

??? success "解答（点击展开）"

    联合概率质量函数可写成

    \[
    f_\theta(x_1,\ldots,x_n)
    =\left(\prod_{i=1}^n\frac1{x_i!}\right)
    e^{-n\theta}\theta^{\sum_{i=1}^nx_i}.
    \]

    因此 $S$ 是充分统计量。又因为

    \[
    S\sim\operatorname{Pois}(n\theta),
    \]

    Poisson 分布族完备，所以 $S$ 是完备充分统计量。

    对任意实数 $q$，Poisson 分布的概率母函数给出

    \[
    E_\theta(q^S)=\exp\{n\theta(q-1)\}.
    \]

    取

    \[
    q=1-\frac an
    \]

    则

    \[
    E_\theta\left[\left(1-\frac an\right)^S\right]
    =e^{-a\theta}.
    \]

    由 Lehmann-Scheffé 定理，

    \[
    \boxed{
    \widehat{e^{-a\theta}}_{\mathrm{UMVUE}}
    =\left(1-\frac an\right)^S
    }.
    \]

    当 $a=n$ 且 $S=0$ 时，将 $0^0$ 解释为 $1$，这与概率母函数的幂级数表示一致。

---

## 习题 3：端点分布族中的最佳无偏估计

**题目。** 设 $X_1,\ldots,X_n$ 独立同分布，其密度为

\[
f(x\mid\theta)=c(\theta)m(x)\mathbf 1_{\{a<x<\theta\}},
\]

其中

\[
M(\theta)=\int_a^\theta m(x)\,dx=\frac1{c(\theta)}.
\]

令 $Y=X_{(n)}$。其密度为

\[
f_Y(y\mid\theta)
=nm(y)c(\theta)^nM(y)^{n-1}\mathbf 1_{\{a<y<\theta\}}.
\]

(i) 利用因子分解定理证明 $Y$ 充分，并直接根据完备性的定义证明 $Y$ 完备。

(ii) 设 $h$ 可微。写出 $\delta(Y)$ 是 $h(\theta)$ 无偏估计量的方程，对该方程关于 $\theta$ 求导，从而得到 $h(\theta)$ 的 UMVUE。

(iii) 对 $X_i\sim\operatorname{Unif}(0,\theta)$ 和 $h(\theta)=\theta^r$（$r>0$）写出结果，并说明为什么 $(n+r)Y^r/n$ 是一般公式的特殊情形。

??? success "解答（点击展开）"

    **(i) 充分性与完备性**

    联合密度为

    \[
    f_\theta(x_1,\ldots,x_n)
    =c(\theta)^n\prod_{i=1}^nm(x_i)
    \prod_{i=1}^n\mathbf 1_{\{a<x_i<\theta\}}.
    \]

    将其改写为

    \[
    f_\theta(x_1,\ldots,x_n)
    =\left(\prod_{i=1}^nm(x_i)\mathbf 1_{\{x_i>a\}}\right)
    c(\theta)^n\mathbf 1_{\{x_{(n)}<\theta\}}.
    \]

    由因子分解定理，$Y=X_{(n)}$ 是 $\theta$ 的充分统计量。

    为证明完备性，设可测函数 $q$ 满足

    \[
    E_\theta[q(Y)]=0
    \qquad\text{对所有 }\theta>a.
    \]

    利用 $Y$ 的密度，得到

    \[
    nc(\theta)^n
    \int_a^\theta q(y)m(y)M(y)^{n-1}\,dy=0.
    \]

    因为 $c(\theta)>0$，所以等价地，

    \[
    \int_a^\theta q(y)m(y)M(y)^{n-1}\,dy=0
    \qquad\text{对所有 }\theta>a.
    \]

    对 $\theta$ 求导可得

    \[
    q(\theta)m(\theta)M(\theta)^{n-1}=0
    \]

    几乎处处成立。在分布的支撑上，$m(\theta)>0$ 且 $M(\theta)>0$，故 $q(\theta)=0$ 几乎处处成立。因此 $Y$ 完备。

    **(ii) 一般 UMVUE 公式**

    若 $\delta(Y)$ 是 $h(\theta)$ 的无偏估计量，则

    \[
    h(\theta)
    =nc(\theta)^n
    \int_a^\theta\delta(y)m(y)M(y)^{n-1}\,dy.
    \]

    因为 $c(\theta)=1/M(\theta)$，所以

    \[
    h(\theta)M(\theta)^n
    =n\int_a^\theta\delta(y)m(y)M(y)^{n-1}\,dy.
    \]

    对两边关于 $\theta$ 求导，得到

    \[
    h'(\theta)M(\theta)^n
    +nh(\theta)M(\theta)^{n-1}m(\theta)
    =n\delta(\theta)m(\theta)M(\theta)^{n-1}.
    \]

    因此

    \[
    \delta(\theta)
    =h(\theta)+\frac{M(\theta)}{nm(\theta)}h'(\theta).
    \]

    用 $Y$ 代替 $\theta$，得到无偏估计量

    \[
    \delta(Y)
    =h(Y)+\frac{M(Y)}{nm(Y)}h'(Y).
    \]

    由于 $Y$ 完备且充分，Lehmann-Scheffé 定理说明

    \[
    \boxed{
    \widehat{h(\theta)}_{\mathrm{UMVUE}}
    =h(Y)+\frac{M(Y)}{nm(Y)}h'(Y)
    }.
    \]

    **(iii) 均匀分布的特殊情形**

    对 $X_i\sim\operatorname{Unif}(0,\theta)$，有

    \[
    a=0,
    \qquad
    m(x)=1,
    \qquad
    M(y)=\int_0^y1\,dx=y.
    \]

    当 $h(\theta)=\theta^r$ 时，

    \[
    h'(y)=ry^{r-1}.
    \]

    代入一般公式，

    \[
    \delta(Y)
    =Y^r+\frac Yn rY^{r-1}
    =\left(1+\frac rn\right)Y^r.
    \]

    所以

    \[
    \boxed{
    \widehat{\theta^r}_{\mathrm{UMVUE}}
    =\frac{n+r}{n}Y^r
    }.
    \]

---

## 习题 4：Bernoulli 样本均值与 Cramér-Rao 下界

**题目。** 设

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(p).
\]

证明 $\overline X$ 是 $p$ 的无偏估计量，并且达到 Cramér-Rao 下界。

??? success "解答（点击展开）"

    因为 $E_p(X_i)=p$，所以

    \[
    E_p(\overline X)
    =\frac1n\sum_{i=1}^nE_p(X_i)
    =p.
    \]

    因此 $\overline X$ 无偏。又由独立性，

    \[
    \operatorname{Var}_p(\overline X)
    =\frac1{n^2}\sum_{i=1}^n\operatorname{Var}_p(X_i)
    =\frac{p(1-p)}n.
    \]

    对数似然函数为

    \[
    \ell(p)
    =\left(\sum_{i=1}^nX_i\right)\log p
    +\left(n-\sum_{i=1}^nX_i\right)\log(1-p).
    \]

    二阶导数为

    \[
    \ell''(p)
    =-\frac{\sum_{i=1}^nX_i}{p^2}
    -\frac{n-\sum_{i=1}^nX_i}{(1-p)^2}.
    \]

    样本中的 Fisher 信息为

    \[
    I_n(p)
    =-E_p[\ell''(p)]
    =\frac{np}{p^2}+\frac{n(1-p)}{(1-p)^2}
    =\frac{n}{p(1-p)}.
    \]

    对 $p$ 的任意无偏估计量，Cramér-Rao 下界为

    \[
    \frac1{I_n(p)}=\frac{p(1-p)}n.
    \]

    由于

    \[
    \operatorname{Var}_p(\overline X)
    =\frac{p(1-p)}n
    =\frac1{I_n(p)},
    \]

    所以

    \[
    \boxed{\overline X\text{ 达到 Cramér-Rao 下界}}.
    \]

---

## 习题 5：矩乘积的一个具体 U-统计量

**题目。** 设 $X_1,\ldots,X_n$ 独立同分布，并且

\[
E(X_1^2)<\infty,
\qquad
E(X_1^{-2})<\infty.
\]

记

\[
\mu=E(X_1),
\qquad
\overline\mu=E(X_1^{-1}).
\]

求一个对 $\mu\overline\mu$ 无偏的 U-统计量。

??? success "解答（点击展开）"

    对两个独立观测，定义对称核函数

    \[
    h(x_1,x_2)
    =\frac12\left(\frac{x_1}{x_2}+\frac{x_2}{x_1}\right).
    \]

    题设的矩条件保证所需期望有限。由独立性，

    \[
    E[h(X_1,X_2)]
    =\frac12\left\{E(X_1)E(X_2^{-1})
    +E(X_2)E(X_1^{-1})\right\}.
    \]

    因而

    \[
    E[h(X_1,X_2)]=\mu\overline\mu.
    \]

    所以，与 $h$ 对应的二阶 U-统计量

    \[
    \boxed{
    U_n
    =\binom n2^{-1}
    \sum_{1\leq i<j\leq n}
    \frac12\left(\frac{X_i}{X_j}+\frac{X_j}{X_i}\right)
    }
    \]

    是 $\mu\overline\mu$ 的无偏估计量。

---

## 习题 6：Cramér-Rao 下界在非正则模型中失效的原因

**题目。** 设

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Unif}(0,\theta),
\qquad \theta>0,
\]

并令 $Y=X_{(n)}$。已知

\[
f_Y(y\mid\theta)
=\frac{ny^{n-1}}{\theta^n}\mathbf 1_{\{0<y<\theta\}}.
\]

(i) 对 $E_\theta(Y)$ 应用通常的积分与微分交换时，说明该交换为何失效。

(ii) 验证

\[
\widehat\theta=\frac{n+1}{n}Y
\]

无偏。利用习题 3 中端点分布族的完备充分性说明它是 UMVUE。计算其方差，并说明形式上的正则 Fisher 信息不等式为何不成立。

(iii) 指出通常 Cramér-Rao 证明中的哪一条正则性条件失效，并解释移动支撑集的作用。

??? success "解答（点击展开）"

    **(i) 积分与微分不能直接交换**

    直接计算，

    \[
    E_\theta(Y)
    =\int_0^\theta y\frac{ny^{n-1}}{\theta^n}\,dy
    =\frac{n}{n+1}\theta.
    \]

    因此

    \[
    \frac d{d\theta}E_\theta(Y)=\frac{n}{n+1}.
    \]

    如果忽略积分上限随 $\theta$ 变化，只对密度求导，就会得到

    \[
    \int_0^\theta y\frac\partial{\partial\theta}
    \left(\frac{ny^{n-1}}{\theta^n}\right)dy
    =-\frac{n^2}{n+1},
    \]

    这与 $n/(n+1)$ 不相等。缺少的正是边界项。由 Leibniz 法则，

    \[
    \frac d{d\theta}
    \int_0^\theta yf_Y(y\mid\theta)\,dy
    =\int_0^\theta y\frac\partial{\partial\theta}
    f_Y(y\mid\theta)\,dy
    +\theta f_Y(\theta\mid\theta).
    \]

    右端为

    \[
    -\frac{n^2}{n+1}+n
    =\frac{n}{n+1}.
    \]

    **(ii) UMVUE、方差及形式 CRLB 的失效**

    因为

    \[
    E_\theta(Y)=\frac{n}{n+1}\theta,
    \]

    所以

    \[
    E_\theta\left(\frac{n+1}{n}Y\right)=\theta.
    \]

    因此 $\widehat\theta=(n+1)Y/n$ 无偏。由习题 3，$Y$ 在该端点分布族中完备且充分，故 Lehmann-Scheffé 定理说明 $\widehat\theta$ 是 UMVUE。

    其次，

    \[
    E_\theta(Y^2)
    =\int_0^\theta y^2\frac{ny^{n-1}}{\theta^n}\,dy
    =\frac{n}{n+2}\theta^2.
    \]

    因而

    \[
    \operatorname{Var}_\theta(Y)
    =\frac{n}{n+2}\theta^2
    -\left(\frac{n}{n+1}\theta\right)^2
    =\frac{n}{(n+2)(n+1)^2}\theta^2.
    \]

    所以

    \[
    \boxed{
    \operatorname{Var}_\theta(\widehat\theta)
    =\frac{\theta^2}{n(n+2)}
    }.
    \]

    如果形式上对单个观测使用正则模型的 Fisher 信息计算，那么在支撑集内部，

    \[
    \frac\partial{\partial\theta}\log f(X\mid\theta)
    =-\frac1\theta.
    \]

    因此形式上的单个观测 Fisher 信息为

    \[
    I_1^{\mathrm{formal}}(\theta)
    =E_\theta\left[
    \left\{\frac\partial{\partial\theta}
    \log f(X\mid\theta)\right\}^2
    \right]
    =\frac1{\theta^2}.
    \]

    按照正则模型的信息可加性，形式上将有

    \[
    I_n^{\mathrm{formal}}(\theta)=\frac n{\theta^2},
    \]

    从而得到形式上的 Cramér-Rao 下界

    \[
    \frac1{I_n^{\mathrm{formal}}(\theta)}
    =\frac{\theta^2}{n}.
    \]

    但是

    \[
    \frac{\theta^2}{n(n+2)}<\frac{\theta^2}{n}.
    \]

    因此形式上的正则 Fisher 信息不等式被“违反”。这并不构成矛盾，因为该模型并非正则模型。

    **(iii) 失效的正则性条件**

    关键问题是密度的支撑集依赖于参数，使得积分与微分不能按通常方式交换。这里的支撑集为 $(0,\theta)$，它随 $\theta$ 移动。

    对单个观测，

    \[
    \int_0^\theta f(x\mid\theta)\,dx=1.
    \]

    但如果只对被积函数求导，则

    \[
    \int_0^\theta
    \frac\partial{\partial\theta}f(x\mid\theta)\,dx
    =\int_0^\theta\left(-\frac1{\theta^2}\right)dx
    =-\frac1\theta\neq0.
    \]

    缺失的边界项为

    \[
    f(\theta\mid\theta)=\frac1\theta,
    \]

    它恰好使总概率的导数恢复为零。等价地，通常的恒等式

    \[
    E_\theta\left[
    \frac\partial{\partial\theta}\log f(X\mid\theta)
    \right]=0
    \]

    在这里不成立。因此，标准 Cramér-Rao 证明不能应用。
