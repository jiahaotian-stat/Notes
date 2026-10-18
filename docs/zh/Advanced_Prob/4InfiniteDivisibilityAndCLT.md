# 第四节：无穷可分分布、无穷小阵列与中心极限定理

本节首先完成无穷可分分布的 Lévy–Khintchine 表示，并以正态分布和 Poisson 分布为例说明其结构。随后引入无穷小三角阵列，解释无穷可分分布为什么恰好是这类阵列行和的可能极限。最后进入中心极限定理，分别介绍特征函数方法、Lindeberg 替换技巧与 Stein 方法。

---

## 1. 无穷可分分布

### 1.1 定义

!!! info "定义 1.1（无穷可分分布）"

    设随机变量 $X$ 的分布函数为 $F$，特征函数为 $f$。如果对每个正整数 $n$，都存在独立同分布随机变量

    \[
    X_{n,1},\ldots,X_{n,n},
    \]

    使得

    \[
    X\overset{d}{=}X_{n,1}+cdots+X_{n,n},
    \]

    则称 $F$ 为**无穷可分分布**，称 $f$ 为**无穷可分特征函数**，简记为 i.d.d.f. 与 i.d.c.f.

从特征函数的角度看，这等价于：对每个 $n\geq1$，存在某个特征函数 $f_n$，使得

\[
f(t)=\bigl(f_n(t)\bigr)^n,
\qquad t\in\mathbb{R}.
\]

!!! note "基本例子"

    正态分布、Poisson 分布以及复合 Poisson 分布都是无穷可分分布。

### 1.2 Lévy–Khintchine 表示

!!! success "定理 1.2（Lévy–Khintchine 典范表示）"

    特征函数 $f$ 是无穷可分特征函数，当且仅当存在常数 $\gamma\in\mathbb{R}$ 和有界、右连续、单调不减函数 $G$，满足 $G(-\infty)=0$，使得

    \[
    f(t)=\exp\{\psi(t)\},
    \]

    其中

    \[
    \psi(t)
    =i\gamma t
    +\int_{-\infty}^{\infty}
    \left(
    e^{itx}-1-\frac{itx}{1+x^2}
    \right)
    \frac{1+x^2}{x^2}\,dG(x).
    \]

    当 $x=0$ 时，积分核按连续延拓解释为

    \[
    \lim_{x\to0}
    \left(
    e^{itx}-1-\frac{itx}{1+x^2}
    \right)
    \frac{1+x^2}{x^2}
    =-\frac{t^2}{2}.
    \]

有序对 $(\gamma,G)$ 称为该无穷可分分布的一组**典范参数**。函数 $G$ 在原点的跳跃给出高斯部分；远离原点的增量给出跳跃部分。

??? proof "充分性的构造证明（点击展开）"

    设 $\psi$ 具有上述形式。线性项 $i\gamma t$ 是常数平移的特征指数。原点处若有跳跃

    \[
    \sigma^2=G(0+)-G(0-),
    \]

    则它对指数的贡献为

    \[
    -\frac{\sigma^2t^2}{2},
    \]

    这正是中心正态分布 $N(0,\sigma^2)$ 的特征指数。

    下面考虑不含原点的有限区间 $[a,b]$。取分割

    \[
    a=x_0<x_1<\cdots<x_m=b
    \]

    以及 $\xi_j\in(x_{j-1},x_j]$。相应的 Riemann–Stieltjes 和中的第 $j$ 项可写成

    \[
    \lambda_j\left(e^{it\xi_j}-1\right)-itc_j,
    \]

    其中

    \[
    \lambda_j
    =\frac{1+\xi_j^2}{\xi_j^2}
    \bigl(G(x_j)-G(x_{j-1})\bigr)\geq0,
    \]

    且

    \[
    c_j
    =\frac{\xi_j}{1+\xi_j^2}\lambda_j.
    \]

    因而该项的指数函数分解为

    \[
    \exp\left\{\lambda_j(e^{it\xi_j}-1)\right\}
    \exp\{-itc_j\}.
    \]

    第一因子是跳跃大小为 $\xi_j$、强度为 $\lambda_j$ 的复合 Poisson 分布的特征函数，第二因子是常数平移的特征函数，所以二者乘积仍是无穷可分特征函数。

    令分割网格趋于零，再令 $a\downarrow-\infty$、$b\uparrow\infty$。由特征函数的连续性定理，极限 $e^{\psi(t)}$ 仍是特征函数。并且对每个 $n$，

    \[
    e^{\psi(t)}
    =\left(e^{\psi(t)/n}\right)^n,
    \]

    而 $e^{\psi(t)/n}$ 仍由参数 $(\gamma/n,G/n)$ 给出，因此 $f$ 无穷可分。$\square$

??? proof "必要性的证明思路（教材补充，点击展开）"

    若 $f$ 无穷可分，则对每个 $n$ 可写成 $f=f_n^n$。把对应的 $n$ 个独立同分布变量视为一行，得到一个无穷小三角阵列。对每一行的截断均值、截断方差与大跳跃分布分别作紧性分析，并抽取收敛子列，可以得到极限参数 $\gamma$ 与 $G$。再对

    \[
    n\log f_n(t)
    \]

    作一致的小量展开，极限正是上述 $\psi(t)$。参数表示的唯一性保证不同子列得到同一个 $(\gamma,G)$，从而得到完整序列的 Lévy–Khintchine 表示。$\square$

### 1.3 正态分布

若 $X\sim N(\mu,\sigma^2)$，则

\[
f(t)=\exp\left\{i\mu t-\frac{\sigma^2t^2}{2}\right\}.
\]

在典范表示中，可以取

\[
\gamma=\mu,
\]

并令 $G$ 只在原点有大小为 $\sigma^2$ 的跳跃，即

\[
G(x)=
\begin{cases}
0, & x<0,\\
\sigma^2, & x\geq0.
\end{cases}
\]

因此正态分布的无穷可分性来自高斯部分。

### 1.4 Poisson 分布

若 $X\sim\operatorname{Pois}(\lambda)$，则

\[
f(t)=\exp\left\{\lambda(e^{it}-1)\right\}.
\]

可以取

\[
G(x)=
\begin{cases}
0, & x<1,\\
\lambda/2, & x\geq1,
\end{cases}
\qquad
\gamma=\frac{\lambda}{2}.
\]

事实上，$G$ 在 $1$ 处的跳跃给出

\[
2\cdot\frac{\lambda}{2}
\left(e^{it}-1-\frac{it}{2}\right),
\]

再与 $i\gamma t=i\lambda t/2$ 合并，便得到 $\lambda(e^{it}-1)$。

---

## 2. 无穷小三角阵列与极限分布

### 2.1 无穷小三角阵列

考虑三角阵列

\[
\{X_{n,k}:1\leq k\leq k_n,\ n\geq1\},
\]

并假设每一行中的随机变量相互独立。记第 $n$ 行的和为

\[
S_n=\sum_{k=1}^{k_n}X_{n,k}.
\]

!!! info "定义 2.1（无穷小阵列）"

    如果对任意 $\varepsilon>0$，都有

    \[
    \max_{1\leq k\leq k_n}
    P\bigl(|X_{n,k}|>\varepsilon\bigr)
    \longrightarrow0,
    \]

    则称该三角阵列为**无穷小的**。

该条件表示：随着 $n\to\infty$，一行中任何单独一项都不再对总和产生宏观影响。

### 2.2 极限必为无穷可分分布

!!! success "定理 2.2（无穷小阵列的极限）"

    若三角阵列每行独立且无穷小，并且

    \[
    S_n\overset{d}{\longrightarrow}S,
    \]

    则 $S$ 的分布是无穷可分的。

??? proof "证明思路（点击展开）"

    固定正整数 $m$。利用无穷小性，可以把第 $n$ 行按顺序分成 $m$ 个相互独立的块，使每一块都只贡献总特征指数的大约 $1/m$。设这些块和为

    \[
    Y_{n,1}^{(m)},\ldots,Y_{n,m}^{(m)}.
    \]

    通过紧性与对角子列论证，可以使每个块和沿同一子列收敛，并调整可忽略的中心化常数，使 $m$ 个极限块具有相同分布。于是

    \[
    S\overset{d}{=}Y_1^{(m)}+cdots+Y_m^{(m)},
    \]

    其中 $Y_1^{(m)},\ldots,Y_m^{(m)}$ 独立同分布。由于 $m$ 任意，$S$ 的分布无穷可分。$\square$

### 2.3 逆命题

!!! success "定理 2.3（无穷可分分布的阵列表示）"

    若 $S$ 具有无穷可分分布，则存在一个每行独立的无穷小三角阵列，使其行和依分布收敛到 $S$；事实上，可以让每一行的和都与 $S$ 同分布。

??? proof "证明（点击展开）"

    由无穷可分性，对每个 $n$，存在独立同分布的

    \[
    X_{n,1},\ldots,X_{n,n}
    \]

    使得

    \[
    S\overset{d}{=}\sum_{k=1}^nX_{n,k}.
    \]

    若 $f$ 是 $S$ 的特征函数，则可以选取第 $n$ 个卷积根 $f_n$，使

    \[
    \bigl(f_n(t)\bigr)^n=f(t).
    \]

    因为 $f_n(t)\to1$，Lévy 连续性定理给出

    \[
    X_{n,1}\overset{d}{\longrightarrow}0.
    \]

    极限为常数时，依分布收敛等价于依概率收敛，所以

    \[
    P\bigl(|X_{n,1}|>\varepsilon\bigr)\longrightarrow0.
    \]

    同一行内各变量同分布，故

    \[
    \max_{1\leq k\leq n}
    P\bigl(|X_{n,k}|>\varepsilon\bigr)
    =P\bigl(|X_{n,1}|>\varepsilon\bigr)
    \longrightarrow0.
    \]

    因此该阵列无穷小，而每一行的和都与 $S$ 同分布。$\square$

### 2.4 Bernoulli 和的三个极限

设 $\xi_1,\xi_2,\ldots$ 独立同分布，且

\[
P(\xi_i=1)=p,
\qquad
P(\xi_i=0)=q=1-p.
\]

记

\[
T_n=\sum_{i=1}^n\xi_i.
\]

#### 大数定律尺度

由弱大数定律，

\[
\frac{T_n}{n}\overset{P}{\longrightarrow}p.
\]

将 $X_{n,i}=\xi_i/n$ 看成三角阵列元素，则对任意固定 $\varepsilon>0$，当 $n>1/\varepsilon$ 时，

\[
P\bigl(|X_{n,i}|>\varepsilon\bigr)=0,
\]

所以该阵列无穷小。

#### 中心极限定理尺度

有

\[
\frac{T_n-np}{\sqrt{npq}}
=\sum_{i=1}^n\frac{\xi_i-p}{\sqrt{npq}}
\overset{d}{\longrightarrow}N(0,1).
\]

对应阵列元素

\[
X_{n,i}=\frac{\xi_i-p}{\sqrt{npq}}
\]

满足 $|X_{n,i}|\leq1/\sqrt{npq}\to0$，因而也是无穷小阵列。

#### Poisson 稀有事件尺度

设 $\xi_{n,1},\ldots,\xi_{n,n}$ 在第 $n$ 行独立同分布，且

\[
P(\xi_{n,k}=1)=p_n,
\qquad
P(\xi_{n,k}=0)=1-p_n.
\]

若

\[
np_n\longrightarrow\lambda\in(0,\infty),
\]

则 $p_n\to0$，阵列无穷小，并且

\[
\sum_{k=1}^n\xi_{n,k}
\overset{d}{\longrightarrow}
\operatorname{Pois}(\lambda).
\]

其特征函数证明为

\[
\mathbb{E}\left[
\exp\left\{
it\sum_{k=1}^{n}\xi_{n,k}
\right\}
\right]
=
\left(
1+p_n\left(e^{it}-1\right)
\right)^n
\longrightarrow
\exp\left\{
\lambda\left(e^{it}-1\right)
\right\}.
\]

---

## 3. Lévy 连续性定理

设 $Z_n$ 与 $Z$ 的特征函数分别为 $\varphi_n$ 与 $\varphi$。

!!! success "定理 3.1（Lévy 连续性定理）"

    若

    \[
    Z_n\overset{d}{\longrightarrow}Z,
    \]

    则对每个 $t\in\mathbb{R}$，

    \[
    \varphi_n(t)\longrightarrow\varphi(t).
    \]

    反之，若 $\varphi_n(t)$ 对每个 $t$ 收敛到函数 $\varphi(t)$，且 $\varphi$ 在 $0$ 处连续，则 $\varphi$ 是某个随机变量 $Z$ 的特征函数，并且

    \[
    Z_n\overset{d}{\longrightarrow}Z.
    \]

??? proof "正向结论（点击展开）"

    对固定 $t$，函数 $x\mapsto e^{itx}$ 有界且连续。由依分布收敛的 Portmanteau 定理，

    \[
    Ee^{itZ_n}\longrightarrow Ee^{itZ}.
    \]

    即 $\varphi_n(t)\to\varphi(t)$。$\square$

??? proof "逆向结论的证明思路（点击展开）"

    $\varphi$ 在原点连续，再结合 $\varphi_n(0)=1$，可推出 $\{Z_n\}$ 的分布族紧。由 Helly 选择定理，任意子列都有进一步子列依分布收敛到某个 $Z$。正向结论说明该子列的极限特征函数必须是 $\varphi$。特征函数的唯一性表明所有收敛子列具有相同极限分布，因此整个序列依分布收敛，且 $\varphi$ 正是极限分布的特征函数。$\square$

---

## 4. 独立同分布中心极限定理

设 $X_1,X_2,ldots$ 独立同分布，并满足

\[
EX_1=\mu,
\qquad
\operatorname{Var}(X_1)=\sigma^2\in(0,\infty).
\]

令

\[
Z_n=rac{\sum_{i=1}^n(X_i-\mu)}{\sigma\sqrt{n}}.
\]

!!! success "定理 4.1（Lindeberg–Lévy 中心极限定理）"

    在上述条件下，

    \[
    Z_n\overset{d}{\longrightarrow}N(0,1).
    \]

### 4.1 特征函数证明

??? proof "证明（点击展开）"

    令

    \[
    Y_i=\frac{X_i-\mu}{\sigma},
    \]

    则 $EY_i=0$、$EY_i^2=1$。记 $Y_1$ 的特征函数为 $g$。由二阶矩有限，$g$ 在原点具有展开

    \[
    g(u)=1-\frac{u^2}{2}+o(u^2),
    \qquad u\to0.
    \]

    利用独立性，$Z_n=n^{-1/2}\sum_{i=1}^nY_i$ 的特征函数为

    \[
    \varphi_n(t)
    =\left[g\left(\frac{t}{\sqrt n}\right)\right]^n.
    \]

    因此

    \[
    g\left(\frac{t}{\sqrt n}\right)
    =1-\frac{t^2}{2n}+o\left(\frac1n\right),
    \]

    从而

    \[
    \varphi_n(t)
    \longrightarrow e^{-t^2/2}.
    \]

    右端是标准正态分布的特征函数，故由 Lévy 连续性定理，

    \[
    Z_n\overset{d}{\longrightarrow}N(0,1).
    \]

    $\square$

### 4.2 Lindeberg 替换技巧

将变量标准化后，可设 $EX_i=0$、$EX_i^2=1$。再取独立同分布的 $Y_i\sim N(0,1)$，并假设两列变量彼此独立。对足够光滑的有界函数 $f$，定义

\[
W_{n,r}
=\frac{Y_1+\cdots+Y_{r-1}+X_{r+1}+\cdots+X_n}{\sqrt n}.
\]

则望远镜分解给出

\[
\begin{aligned}
&E f\left(\frac{X_1+\cdots+X_n}{\sqrt n}\right)
-E f\left(\frac{Y_1+\cdots+Y_n}{\sqrt n}\right)\\
&=\sum_{r=1}^n
E\left[
f\left(W_{n,r}+\frac{X_r}{\sqrt n}\right)
-f\left(W_{n,r}+\frac{Y_r}{\sqrt n}\right)
\right].
\end{aligned}
\]

??? proof "Lindeberg 替换论证（点击展开）"

    在点 $W_{n,r}$ 处作二阶 Taylor 展开：

    \[
    f(W_{n,r}+u)
    =f(W_{n,r})+uf'(W_{n,r})
    +\frac{u^2}{2}f''(W_{n,r})+R_{n,r}(u).
    \]

    由于 $X_r$、$Y_r$ 与 $W_{n,r}$ 独立，并且

    \[
    EX_r=EY_r=0,
    \qquad
    EX_r^2=EY_r^2=1,
    \]

    每个差值中的零阶、一次与二次项全部抵消，只剩余 Taylor 余项。若先对 $X_r$ 作截断，再利用 $E X_r^2<\infty$ 与 $f''$ 的一致连续性，则

    \[
    \sum_{r=1}^n
    E\left|R_{n,r}\left(\frac{X_r}{\sqrt n}\right)
    -R_{n,r}\left(\frac{Y_r}{\sqrt n}\right)\right|
    \longrightarrow0.
    \]

    因而

    \[
    E f\left(\frac{X_1+\cdots+X_n}{\sqrt n}\right)
    -E f\left(\frac{Y_1+\cdots+Y_n}{\sqrt n}\right)
    \longrightarrow0.
    \]

    而

    \[
    \frac{Y_1+\cdots+Y_n}{\sqrt n}\sim N(0,1),
    \]

    所以对一类决定弱收敛的光滑测试函数均有期望收敛，再由光滑逼近推广到所有 $f\in C_b(\mathbb R)$，得到中心极限定理。$\square$

### 4.3 教材中的 Lindeberg 条件

上述替换思想可推广到每行不必同分布的独立三角阵列。设 $EX_{n,k}=0$，并记

\[
B_n=\sum_{k=1}^{k_n}EX_{n,k}^2.
\]

!!! info "定义 4.2（Lindeberg 条件）"

    对任意 $\varepsilon>0$，若

    \[
    \frac1{B_n}
    \sum_{k=1}^{k_n}
    E\left[
    X_{n,k}^2
    \mathbf 1_{\{|X_{n,k}|>\varepsilon\sqrt{B_n}\}}
    \right]
    \longrightarrow0,
    \]

    则称该阵列满足 Lindeberg 条件。

!!! success "定理 4.3（Lindeberg–Feller 定理的充分性）"

    若每行变量相互独立、中心化，$B_n>0$，且满足 Lindeberg 条件，则

    \[
    \frac{\sum_{k=1}^{k_n}X_{n,k}}{\sqrt{B_n}}
    \overset{d}{\longrightarrow}N(0,1).
    \]

Lindeberg 条件排除了单个大跳跃主导总和的可能，并蕴含 Feller 无穷小条件

\[
\max_{1\leq k\leq k_n}
\frac{EX_{n,k}^2}{B_n}
\longrightarrow0.
\]

---

## 5. Stein 方法

### 5.1 正态分布的分部积分恒等式

设 $Z\sim N(0,1)$。对适当光滑且可积的函数 $f$，有

\[
E\bigl[Zf(Z)\bigr]=E\bigl[f'(Z)\bigr].
\]

??? proof "证明（点击展开）"

    记标准正态密度为

    \[
    \phi(z)=\frac1{\sqrt{2\pi}}e^{-z^2/2}.
    \]

    因为 $\phi'(z)=-z\phi(z)$，分部积分得到

    \[
    \begin{aligned}
    E[Zf(Z)]
    &=\int_{-\infty}^{\infty}zf(z)\phi(z)\,dz\\
    &=-\int_{-\infty}^{\infty}f(z)\,d\phi(z)\\
    &=\int_{-\infty}^{\infty}f'(z)\phi(z)\,dz\\
    &=E[f'(Z)].
    \end{aligned}
    \]

    边界项为零。$\square$

### 5.2 Stein 刻画

!!! success "定理 5.1（标准正态分布的 Stein 刻画）"

    随机变量 $W$ 服从 $N(0,1)$，当且仅当对一类足够丰富的绝对连续函数 $f$，

    \[
    E\bigl[f'(W)-Wf(W)\bigr]=0.
    \]

??? proof "逆向证明（点击展开）"

    固定 $z\in\mathbb R$，令

    \[
    h_z(x)=\mathbf 1_{(-\infty,z]}(x)-\Phi(z).
    \]

    考虑 Stein 方程

    \[
    f_z'(x)-xf_z(x)=h_z(x).
    \]

    它有一个有界绝对连续解。若 $W$ 满足 Stein 恒等式，则

    \[
    0=E\bigl[f_z'(W)-Wf_z(W)\bigr]
    =E h_z(W).
    \]

    因而

    \[
    P(W\leq z)=\Phi(z)
    \]

    对每个 $z$ 成立，所以 $W\sim N(0,1)$。$\square$

### 5.3 用 Stein 恒等式证明中心极限定理

再次设 $EX_i=0$、$EX_i^2=1$，并令

\[
Z_n=\frac1{\sqrt n}\sum_{i=1}^nX_i,
\qquad
Z_{n,i}=Z_n-\frac{X_i}{\sqrt n}.
\]

注意 $Z_{n,i}$ 与 $X_i$ 独立。对光滑函数 $f$，

\[
\begin{aligned}
E[Z_nf(Z_n)]
&=\frac1{\sqrt n}\sum_{i=1}^nE[X_if(Z_n)]\\
&=\frac1{\sqrt n}\sum_{i=1}^n
E\left[X_i\bigl(f(Z_n)-f(Z_{n,i})\bigr)\right].
\end{aligned}
\]

第二个等号利用了

\[
E[X_if(Z_{n,i})]=EX_i\,Ef(Z_{n,i})=0.
\]

对差值作一阶展开，

\[
f(Z_n)-f(Z_{n,i})
=\frac{X_i}{\sqrt n}f'(Z_{n,i})+R_{n,i}.
\]

于是

\[
E[Z_nf(Z_n)]
=\frac1n\sum_{i=1}^nE[f'(Z_{n,i})]
+\frac1{\sqrt n}\sum_{i=1}^nE[X_iR_{n,i}].
\]

由 $EX_i^2=1$、截断论证以及 $f'$ 的一致连续性，

\[
\frac1{\sqrt n}\sum_{i=1}^nE[X_iR_{n,i}]
\longrightarrow0,
\]

同时

\[
\frac1n\sum_{i=1}^nE[f'(Z_{n,i})]
-E[f'(Z_n)]
\longrightarrow0.
\]

因此

\[
E\bigl[f'(Z_n)-Z_nf(Z_n)\bigr]
\longrightarrow0.
\]

这说明 $Z_n$ 渐近满足标准正态分布的 Stein 恒等式。结合 Stein 方程与紧性，即得

\[
Z_n\overset{d}{\longrightarrow}N(0,1).
\]

---

## 6. 本节逻辑关系

!!! note "主线总结"

    **无穷可分分布：**可以对任意 $n$ 分解为 $n$ 个独立同分布随机变量之和。

    **Lévy–Khintchine 表示：**用平移、高斯部分与跳跃部分完整刻画无穷可分特征函数。

    **无穷小三角阵列：**其行和的任何弱极限都无穷可分；反之，每个无穷可分分布都能由这种行和得到。

    **Lévy 连续性定理：**把特征函数的逐点收敛转化为分布收敛。

    **中心极限定理：**可以由特征函数、Lindeberg 替换或 Stein 恒等式证明；三种方法分别突出乘积结构、矩匹配与分布刻画。
