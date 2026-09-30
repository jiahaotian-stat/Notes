# 第三章：极限定理、随机收敛与无穷可分分布

本章从 Bernoulli 和的经典极限定理出发，系统介绍几乎处处收敛、依概率收敛与依分布收敛，并讨论连续映射定理、Slutsky 定理、Lévy 连续性定理、Borel-Cantelli 引理以及子列刻画。最后介绍无穷可分分布的定义、基本性质、典型例子、Lévy-Khintchine 表示及其与三角阵列极限的关系。

---

## 1. Bernoulli 和的经典极限定理

设

\[
S_n=\sum_{i=1}^n\xi_i,
\]

其中 $\xi_1,\ldots,\xi_n$ 独立同分布，且

\[
P(\xi_i=1)=p,
\qquad
P(\xi_i=0)=1-p.
\]

于是 $S_n\sim\operatorname{Binomial}(n,p)$。

### 1.1 Bernoulli 弱大数定律

!!! success "定理 1.1（Bernoulli 弱大数定律）"

    当 $n\to\infty$ 时，

    \[
    \frac{S_n}{n}\xrightarrow{P}p.
    \]

??? proof "定理 1.1 的证明（点击展开）"

    由于

    \[
    E\left(\frac{S_n}{n}\right)=p,
    \qquad
    \operatorname{Var}\left(\frac{S_n}{n}\right)
    =\frac{p(1-p)}{n},
    \]

    对任意 $\varepsilon>0$，由 Chebyshev 不等式，

    \[
    P\left(\left|\frac{S_n}{n}-p\right|>\varepsilon\right)
    \leq
    \frac{p(1-p)}{n\varepsilon^2}
    \longrightarrow0.
    \]

    因而 $S_n/n\xrightarrow{P}p$。$\square$

### 1.2 Bernoulli 中心极限定理

!!! success "定理 1.2（De Moivre-Laplace 中心极限定理）"

    若 $0<p<1$，则

    \[
    Z_n
    =\frac{S_n-np}{\sqrt{np(1-p)}}
    \xrightarrow{d}N(0,1).
    \]

记 $F_n$ 为 $Z_n$ 的分布函数，$\Phi$ 为标准正态分布函数。由于 $\Phi$ 连续，Pólya 定理进一步给出

\[
\Delta_n
=\sup_{x\in\mathbb{R}}|F_n(x)-\Phi(x)|
\longrightarrow0.
\]

因此，当 $n$ 很大时，二项分布的标准化尾概率可以用正态尾概率近似：

\[
P(Z_n>x)
=1-F_n(x)
\approx1-\Phi(x).
\]

### 1.3 Poisson 极限定理

对每个 $n$，设 $\xi_{n,1},\ldots,\xi_{n,n}$ 独立同分布，且

\[
P(\xi_{n,i}=1)=p_n,
\qquad
P(\xi_{n,i}=0)=1-p_n.
\]

令

\[
S_n=\sum_{i=1}^n\xi_{n,i}.
\]

!!! success "定理 1.3（Poisson 极限定理）"

    若

    \[
    np_n\longrightarrow\lambda>0,
    \]

    则

    \[
    S_n\xrightarrow{d}S,
    \qquad
    S\sim\operatorname{Poisson}(\lambda).
    \]

??? proof "定理 1.3 的证明（点击展开）"

    对任意固定的非负整数 $k$，

    \[
    P(S_n=k)
    =\binom{n}{k}p_n^k(1-p_n)^{n-k}.
    \]

    由于 $np_n\to\lambda$，有

    \[
    \binom{n}{k}p_n^k
    =\frac{n(n-1)\cdots(n-k+1)}{k!}p_n^k
    \longrightarrow\frac{\lambda^k}{k!},
    \]

    且

    \[
    (1-p_n)^{n-k}\longrightarrow e^{-\lambda}.
    \]

    因此

    \[
    P(S_n=k)
    \longrightarrow
    e^{-\lambda}\frac{\lambda^k}{k!}
    =P(S=k).
    \]

    右端对 $k\geq0$ 求和为 $1$，故 $S_n\xrightarrow{d}S$。$\square$

!!! note "适用范围"

    独立同分布假设可以推广。例如，大数定律和中心极限定理可以研究非独立序列、Markov 链及满足混合条件的过程；极值统计量则常收敛到 Gumbel 等极值分布。随机变量还可以取值于一般 Banach 空间，此时极限定理与空间的几何结构密切相关。

---

## 2. 随机变量序列的三种收敛

设 $X_n,X$ 定义在同一概率空间 $(\Omega,\mathcal{A},P)$ 上。

### 2.1 几乎处处收敛

!!! info "定义 2.1（几乎处处收敛）"

    若存在 $\Omega_0\in\mathcal{A}$，满足 $P(\Omega_0)=1$，并且对每个 $\omega\in\Omega_0$，

    \[
    X_n(\omega)\longrightarrow X(\omega),
    \]

    则称 $X_n$ 几乎处处收敛于 $X$，记为

    \[
    X_n\xrightarrow{a.s.}X.
    \]

等价地，对任意 $\varepsilon>0$，

\[
P\left(\limsup_{n\to\infty}
\{|X_n-X|>\varepsilon\}\right)=0.
\]

其中事件的上极限为

\[
\limsup_{n\to\infty}A_n
=\bigcap_{N=1}^{\infty}\bigcup_{n\geq N}A_n,
\]

表示事件 $A_n$ 发生无穷多次。

### 2.2 依概率收敛

!!! info "定义 2.2（依概率收敛）"

    若对任意 $\varepsilon>0$，

    \[
    P(|X_n-X|>\varepsilon)\longrightarrow0,
    \]

    则称 $X_n$ 依概率收敛于 $X$，记为

    \[
    X_n\xrightarrow{P}X.
    \]

Markov 不等式给出一个常用判别准则。若存在 $r>0$ 使得

\[
E|X_n-X|^r\longrightarrow0,
\]

则对每个 $\varepsilon>0$，

\[
P(|X_n-X|>\varepsilon)
\leq\frac{E|X_n-X|^r}{\varepsilon^r}
\longrightarrow0,
\]

即 $X_n\xrightarrow{P}X$。

### 2.3 依分布收敛

设 $F_n$ 与 $F$ 分别为 $X_n$ 与 $X$ 的分布函数，并记

\[
D_F=\{x\in\mathbb{R}:F\text{ 在 }x\text{ 处连续}\}.
\]

!!! info "定义 2.3（依分布收敛）"

    若对每个 $x\in D_F$，

    \[
    F_n(x)\longrightarrow F(x),
    \]

    则称 $X_n$ 依分布收敛于 $X$，记为

    \[
    X_n\xrightarrow{d}X.
    \]

依分布收敛只取决于边缘分布，不要求 $X_n$ 与 $X$ 定义在同一概率空间上。若

\[
X_n\xrightarrow{d}X
\quad\text{且}\quad
X_n\xrightarrow{d}Y,
\]

则 $X$ 与 $Y$ 同分布。

!!! example "例 2.4（为什么只要求连续点）"

    令 $X_n\equiv1/n$，$X\equiv0$。则 $X_n\to X$，但在 $x=0$ 处，

    \[
    F_n(0)=0,
    \qquad
    F(0)=1.
    \]

    由于 $0$ 是 $F$ 的不连续点，这并不妨碍 $X_n\xrightarrow{d}X$。

---

## 3. 依概率收敛与依分布收敛的运算规则

### 3.1 依概率收敛的代数运算

!!! success "定理 3.1（依概率收敛的代数性质）"

    若

    \[
    X_n\xrightarrow{P}X,
    \qquad
    Y_n\xrightarrow{P}Y,
    \]

    则

    \[
    X_n\pm Y_n\xrightarrow{P}X\pm Y,
    \]

    \[
    X_nY_n\xrightarrow{P}XY,
    \]

    并且当 $P(Y\neq0)=1$ 时，

    \[
    \frac{X_n}{Y_n}\xrightarrow{P}\frac{X}{Y}.
    \]

### 3.2 连续映射定理

!!! success "定理 3.2（连续映射定理）"

    若 $X_n\xrightarrow{P}X$，且 $f:\mathbb{R}\to\mathbb{R}$ 连续，则

    \[
    f(X_n)\xrightarrow{P}f(X).
    \]

    若 $X_n\xrightarrow{d}X$，并且 $f$ 的不连续点集合 $D_f$ 满足

    \[
    P(X\in D_f)=0,
    \]

    则

    \[
    f(X_n)\xrightarrow{d}f(X).
    \]

### 3.3 收敛类型定理与 Slutsky 定理

!!! success "定理 3.3（收敛类型定理）"

    设 $X_n\xrightarrow{d}X$，其中 $X$ 非退化。若

    \[
    a_nX_n+b_n\xrightarrow{d}Y,
    \]

    且 $Y$ 非退化，则存在常数 $a>0$ 与 $b\in\mathbb{R}$，使得

    \[
    a_n\longrightarrow a,
    \qquad
    b_n\longrightarrow b,
    \]

    并且

    \[
    Y\overset{d}{=}aX+b.
    \]

!!! success "定理 3.4（Slutsky 定理）"

    若

    \[
    X_n\xrightarrow{d}X,
    \qquad
    Y_n\xrightarrow{P}c,
    \]

    其中 $c$ 为常数，则

    \[
    X_n+Y_n\xrightarrow{d}X+c,
    \]

    \[
    X_nY_n\xrightarrow{d}cX,
    \]

    且当 $c\neq0$ 时，

    \[
    \frac{X_n}{Y_n}\xrightarrow{d}\frac{X}{c}.
    \]

!!! warning "仅有边缘依分布收敛并不足够"

    一般而言，$X_n\xrightarrow{d}X$ 与 $Y_n\xrightarrow{d}Y$ 并不能推出

    \[
    X_n+Y_n\xrightarrow{d}X+Y.
    \]

    若进一步有联合收敛

    \[
    (X_n,Y_n)\xrightarrow{d}(X,Y),
    \]

    则由连续映射定理，

    \[
    X_n+Y_n\xrightarrow{d}X+Y.
    \]

---

## 4. 特征函数与弱收敛

设 $X_n$ 与 $X$ 的特征函数分别为

\[
\varphi_n(t)=E(e^{itX_n}),
\qquad
\varphi(t)=E(e^{itX}).
\]

### 4.1 Lévy 连续性定理

!!! success "定理 4.1（Lévy 连续性定理）"

    下列两个命题等价：

    1. $X_n\xrightarrow{d}X$；
    2. 对每个 $t\in\mathbb{R}$，

    \[
    \varphi_n(t)\longrightarrow\varphi(t).
    \]

    更一般地，若 $\varphi_n(t)$ 逐点收敛到函数 $\varphi(t)$，且 $\varphi$ 在 $0$ 处连续，则 $\varphi$ 是某个概率分布的特征函数，并且相应分布弱收敛到该分布。

这个定理把分布函数的弱收敛问题转化为特征函数的逐点收敛问题，是证明中心极限定理和研究独立随机变量之和的基本工具。

### 4.2 收敛关系

三种收敛之间有

\[
X_n\xrightarrow{a.s.}X
\quad\Longrightarrow\quad
X_n\xrightarrow{P}X
\quad\Longrightarrow\quad
X_n\xrightarrow{d}X.
\]

一般而言，两个逆命题都不成立。但当极限为常数 $c$ 时，

\[
X_n\xrightarrow{d}c
\quad\Longleftrightarrow\quad
X_n\xrightarrow{P}c.
\]

---

## 5. 几乎处处收敛、Borel-Cantelli 引理与子列

### 5.1 Borel-Cantelli 判别法

!!! success "定理 5.1（完全收敛蕴含几乎处处收敛）"

    若对任意 $\varepsilon>0$，

    \[
    \sum_{n=1}^{\infty}
    P(|X_n-X|>\varepsilon)<\infty,
    \]

    则

    \[
    X_n\xrightarrow{a.s.}X.
    \]

??? proof "定理 5.1 的证明（点击展开）"

    固定 $\varepsilon>0$，令

    \[
    A_n^{(\varepsilon)}
    =\{|X_n-X|>\varepsilon\}.
    \]

    由第一 Borel-Cantelli 引理，

    \[
    \sum_{n=1}^{\infty}P(A_n^{(\varepsilon)})<\infty
    \quad\Longrightarrow\quad
    P(A_n^{(\varepsilon)}\ \text{i.o.})=0.
    \]

    对所有正有理数 $\varepsilon$ 取可数交，即得 $X_n\to X$ 几乎处处。$\square$

!!! example "例 5.2（Bernoulli 强大数定律的四阶矩证明）"

    对 Bernoulli 和 $S_n$，有

    \[
    E\left|\frac{S_n}{n}-p\right|^4
    =O(n^{-2}).
    \]

    因此由 Markov 不等式，

    \[
    P\left(\left|\frac{S_n}{n}-p\right|>\varepsilon\right)
    \leq
    \frac{E|S_n/n-p|^4}{\varepsilon^4}
    =O(n^{-2}).
    \]

    由于 $\sum_n n^{-2}<\infty$，定理 5.1 给出

    \[
    \frac{S_n}{n}\xrightarrow{a.s.}p.
    \]

### 5.2 依概率收敛的子列刻画

!!! success "定理 5.3（Riesz 子列原理）"

    $X_n\xrightarrow{P}X$ 当且仅当每个子列 $\{X_{n_k}\}$ 都存在进一步子列 $\{X_{n_{k_j}}\}$，使得

    \[
    X_{n_{k_j}}\xrightarrow{a.s.}X.
    \]

??? proof "充分方向的构造（点击展开）"

    假设 $X_n\xrightarrow{P}X$。从任意子列中递归选取 $n_{k_j}$，使得

    \[
    P\left(|X_{n_{k_j}}-X|>2^{-j}\right)<2^{-j}.
    \]

    因而

    \[
    \sum_{j=1}^{\infty}
    P\left(|X_{n_{k_j}}-X|>2^{-j}\right)<\infty.
    \]

    由 Borel-Cantelli 引理，$X_{n_{k_j}}\to X$ 几乎处处。逆方向可用反证法直接得到。$\square$

### 5.3 Skorohod 表示定理

!!! success "定理 5.4（Skorohod 表示定理）"

    若

    \[
    X_n\xrightarrow{d}X,
    \]

    则存在另一个概率空间及其上的随机变量 $X_n',X'$，满足

    \[
    X_n'\overset{d}{=}X_n,
    \qquad
    X'\overset{d}{=}X,
    \]

    并且

    \[
    X_n'\xrightarrow{a.s.}X'.
    \]

!!! warning "注意概率空间已经改变"

    Skorohod 表示定理并不说明原概率空间上的 $X_n$ 几乎处处收敛。它只说明可以在一个新的共同概率空间上构造具有相同边缘分布的版本，使这些版本几乎处处收敛。

---

## 6. 无穷可分分布

### 6.1 定义与等价刻画

!!! info "定义 6.1（无穷可分分布）"

    若对每个正整数 $n$，都存在独立同分布随机变量

    \[
    X_{n,1},\ldots,X_{n,n}
    \]

    使得

    \[
    X\overset{d}{=}
    X_{n,1}+\cdots+X_{n,n},
    \]

    则称 $X$ 的分布为无穷可分分布。

若 $F$ 是 $X$ 的分布函数，$F_n$ 是 $X_{n,1}$ 的分布函数，则定义等价于

\[
F=F_n^{*n},
\]

其中 $*$ 表示卷积。若 $\varphi$ 是 $X$ 的特征函数，$\varphi_n$ 是 $X_{n,1}$ 的特征函数，则等价地，

\[
\varphi(t)=\varphi_n(t)^n,
\qquad t\in\mathbb{R}.
\]

### 6.2 典型例子

!!! example "例 6.2（退化分布）"

    若 $X\equiv c$，则

    \[
    X=\frac{c}{n}+\cdots+\frac{c}{n},
    \]

    所以退化分布无穷可分。其特征函数满足

    \[
    e^{itc}
    =\left(e^{itc/n}\right)^n.
    \]

!!! example "例 6.3（正态分布）"

    若 $X\sim N(\mu,\sigma^2)$，则

    \[
    \varphi(t)
    =\exp\left(i\mu t-\frac{\sigma^2t^2}{2}\right)
    =\left[
    \exp\left(
    i\frac{\mu}{n}t-
    \frac{\sigma^2}{2n}t^2
    \right)
    \right]^n.
    \]

    因此可取

    \[
    X_{n,k}\sim N\left(\frac{\mu}{n},\frac{\sigma^2}{n}\right).
    \]

!!! example "例 6.4（Poisson 分布）"

    若 $X\sim\operatorname{Poisson}(\lambda)$，则

    \[
    \varphi(t)
    =\exp\{\lambda(e^{it}-1)\}
    =\left[
    \exp\left\{\frac{\lambda}{n}(e^{it}-1)\right\}
    \right]^n.
    \]

    因此可取

    \[
    X_{n,k}\sim\operatorname{Poisson}\left(\frac{\lambda}{n}\right).
    \]

!!! example "例 6.5（Cauchy 分布）"

    标准 Cauchy 分布的特征函数为

    \[
    \varphi(t)=e^{-|t|}
    =\left(e^{-|t|/n}\right)^n.
    \]

    因此标准 Cauchy 分布无穷可分。

### 6.3 非无穷可分的例子

Rademacher 分布的特征函数为

\[
\varphi(t)=\cos t,
\]

而 $\operatorname{Uniform}(-1,1)$ 的特征函数为

\[
\varphi(t)=\frac{\sin t}{t}.
\]

两者的特征函数都存在零点。下一节将证明，无穷可分分布的特征函数不能有零点，因此这两个分布都不是无穷可分的。

---

## 7. 无穷可分分布的基本性质

### 7.1 卷积与对称化

!!! success "命题 7.1（封闭性）"

    若 $f$ 与 $g$ 是无穷可分分布的特征函数，则以下函数仍是无穷可分分布的特征函数：

    \[
    f(t)g(t),
    \]

    \[
    f(-t)=\overline{f(t)},
    \]

    以及

    \[
    |f(t)|^2=f(t)f(-t).
    \]

这分别对应独立无穷可分随机变量之和、相反数以及独立同分布随机变量之差。

### 7.2 特征函数没有零点

!!! success "定理 7.2（无零点性质）"

    若 $f$ 是无穷可分分布的特征函数，则

    \[
    f(t)\neq0,
    \qquad \forall t\in\mathbb{R}.
    \]

??? proof "定理 7.2 的证明（点击展开）"

    对每个 $n$，存在特征函数 $f_n$，使得

    \[
    f(t)=f_n(t)^n.
    \]

    假设存在 $t_0$ 使 $f(t_0)=0$。由于 $f(0)=1$ 且 $f$ 连续，可取离原点最近的正零点 $t_0$。在某个 $\delta>0$ 上，$f$ 在 $[-\delta,\delta]$ 内没有零点。

    由

    \[
    |f_n(t)|=|f(t)|^{1/n},
    \]

    可知对固定 $t$ 且 $f(t)\neq0$，

    \[
    |f_n(t)|\longrightarrow1.
    \]

    对任意特征函数 $h$，有不等式

    \[
    0\leq1-|h(2t)|^2
    \leq4\{1-|h(t)|^2\}.
    \]

    从 $|f_n(t)|\to1$ 出发反复使用该不等式，可以把非零性从原点邻域逐步延伸到整个实轴，从而推出 $f(t_0)\neq0$，矛盾。$\square$

### 7.3 非退化有界分布不是无穷可分的

!!! success "定理 7.3（有界无穷可分分布必退化）"

    若 $X$ 无穷可分且存在 $M<\infty$ 使得

    \[
    |X|\leq M
    \quad\text{几乎处处},
    \]

    则 $X$ 必为退化随机变量。

??? proof "定理 7.3 的证明（点击展开）"

    对每个 $n$，写成

    \[
    X\overset{d}{=}
    X_{n,1}+\cdots+X_{n,n},
    \]

    其中 $X_{n,1},\ldots,X_{n,n}$ 独立同分布。由和几乎处处落在 $[-M,M]$ 内以及独立性，可以推出每个分量的本质振幅至多为 $2M/n$；平移各分量并保持总和平移不变后，可令

    \[
    |X_{n,k}|\leq\frac{M}{n}
    \quad\text{几乎处处}.
    \]

    因而

    \[
    \operatorname{Var}(X_{n,1})
    \leq E(X_{n,1}^2)
    \leq\frac{M^2}{n^2}.
    \]

    利用独立性，

    \[
    \operatorname{Var}(X)
    =n\operatorname{Var}(X_{n,1})
    \leq\frac{M^2}{n}.
    \]

    令 $n\to\infty$，得到 $\operatorname{Var}(X)=0$，故 $X$ 退化。$\square$

### 7.4 弱极限封闭性

!!! success "定理 7.4（无穷可分分布对弱收敛封闭）"

    若 $F_m$ 均为无穷可分分布，且

    \[
    F_m\xrightarrow{w}F,
    \]

    则 $F$ 仍为无穷可分分布。

??? proof "定理 7.4 的证明思路（点击展开）"

    设 $f_m$ 与 $f$ 分别为 $F_m$ 与 $F$ 的特征函数。由 Lévy 连续性定理，

    \[
    f_m(t)\longrightarrow f(t).
    \]

    固定 $n$。因为 $F_m$ 无穷可分，存在特征函数 $f_{n,m}$ 使得

    \[
    f_{n,m}(t)^n=f_m(t).
    \]

    由紧性选择一个子列，使 $f_{n,m}$ 沿该子列收敛到某个特征函数 $f_n$。取极限得到

    \[
    f_n(t)^n=f(t).
    \]

    由于这对每个 $n$ 都成立，$F$ 无穷可分。无零点性质保证连续的 $n$ 次根可以一致选取。$\square$

---

## 8. Lévy-Khintchine 表示

!!! success "定理 8.1（Lévy-Khintchine 公式）"

    特征函数 $\varphi$ 对应一个无穷可分分布，当且仅当存在 $\gamma\in\mathbb{R}$ 和 Lévy 测度 $\nu$，满足

    \[
    \nu(\{0\})=0,
    \qquad
    \int_{\mathbb{R}}(1\wedge x^2)\,\nu(dx)<\infty,
    \]

    以及某个 $\sigma^2\geq0$，使得

    \[
    \varphi(t)
    =\exp\left\{
    i\gamma t-\frac{\sigma^2t^2}{2}
    +\int_{\mathbb{R}}
    \left(e^{itx}-1-
    \frac{itx}{1+x^2}\right)
    \nu(dx)
    \right\}.
    \]

等价地，可以把 Gaussian 部分吸收到一个有界非降函数 $G$ 中，写成 Kolmogorov 典范形式

\[
\varphi(t)
=\exp\left\{
i\gamma t+
\int_{-\infty}^{\infty}
\left(e^{itx}-1-
\frac{itx}{1+x^2}\right)
\frac{1+x^2}{x^2}\,dG(x)
\right\}.
\]

其中 $G$ 可以取为左连续、有界、非降函数，并满足

\[
G(-\infty)=0,
\qquad
G(+\infty)<\infty.
\]

三项分别描述漂移、Gaussian 波动与跳跃结构。参数 $\gamma$、$\sigma^2$ 和 $\nu$ 在固定截断函数后唯一。

---

## 9. 无穷小三角阵列的极限

设对每个 $n$，随机变量

\[
\xi_{n,1},\ldots,\xi_{n,k_n}
\]

相互独立，并令

\[
S_n=\sum_{k=1}^{k_n}\xi_{n,k}.
\]

!!! info "定义 9.1（无穷小三角阵列）"

    若对任意 $\varepsilon>0$，

    \[
    \max_{1\leq k\leq k_n}
    P(|\xi_{n,k}|>\varepsilon)
    \longrightarrow0.
    \]

    则称该三角阵列为无穷小的。

!!! success "定理 9.2（无穷小阵列的极限必无穷可分）"

    若三角阵列无穷小，且

    \[
    S_n\xrightarrow{d}S,
    \]

    则 $S$ 的分布是无穷可分分布。

??? proof "定理 9.2 的证明思路（点击展开）"

    固定正整数 $m$。把每一行的独立变量依次分成 $m$ 个相邻块，使每一块对总特征指数的贡献尽量相等，并把每个块内的变量相加。无穷小条件保证单个变量的贡献趋于零，因此各块和在适当子列上具有相同的极限分布。

    记这些块和为 $Y_{n,1},\ldots,Y_{n,m}$。它们相互独立，且

    \[
    S_n=Y_{n,1}+\cdots+Y_{n,m}+o_P(1).
    \]

    沿子列取极限后，得到独立同分布随机变量 $Y_1,ldots,Y_m$，使得

    \[
    S\overset{d}{=}Y_1+\cdots+Y_m.
    \]

    因为 $m$ 任意，$S$ 的分布无穷可分。$\square$

!!! note "极限定理中的作用"

    无穷可分分布正是无穷小独立三角阵列可能出现的极限分布。正态分布、Poisson 分布和更一般的稳定分布都属于这一类。

---

## 10. 本章小结

- Bernoulli 和在不同标准化下分别产生大数定律、中心极限定理与 Poisson 极限。
- 几乎处处收敛蕴含依概率收敛，依概率收敛蕴含依分布收敛；反向一般不成立。
- 连续映射定理、Slutsky 定理与 Lévy 连续性定理是处理弱收敛的核心工具。
- Borel-Cantelli 引理把概率级数的可和性转化为几乎处处收敛；依概率收敛还可由几乎处处收敛子列刻画。
- 无穷可分分布可以被分解为任意多个独立同分布随机变量之和，其特征函数没有零点。
- 无穷可分分布在卷积与弱收敛下封闭，并由 Lévy-Khintchine 公式完全刻画。
- 无穷小独立三角阵列的任何弱极限必为无穷可分分布。
