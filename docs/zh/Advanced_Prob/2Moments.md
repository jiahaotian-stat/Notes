# 第二章：独立性、数字特征、特征函数与分布收敛

本章首先讨论随机变量的独立性及其在可测变换下的保持性，随后介绍期望、方差、矩与 Chebyshev 不等式。在此基础上，我们引入矩母函数与特征函数，并利用特征函数研究独立随机变量之和与分布的唯一性。最后介绍依概率收敛、依分布收敛、Bernoulli 和的极限定理以及中心极限定理。

---

## 1. 随机变量的独立性

### 1.1 两个随机变量的独立性

设 $X,Y$ 是概率空间 $(\Omega,\mathcal{A},P)$ 上的随机变量，联合分布函数为

\[
F_{X,Y}(x,y)=P(X\leq x,Y\leq y).
\]

!!! info "定义 1.1（随机变量的独立性）"

    如果对任意 Borel 集 $A,B\in\mathcal{B}(\mathbb{R})$，都有

    \[
    P(X\in A,Y\in B)
    =P(X\in A)P(Y\in B),
    \]

    则称随机变量 $X$ 与 $Y$ 相互独立。

只需在能够生成 Borel $\sigma$-代数的集合族上验证这个等式。因此，$X$ 与 $Y$ 独立当且仅当

\[
F_{X,Y}(x,y)=F_X(x)F_Y(y),
\qquad \forall x,y\in\mathbb{R}.
\]

对于离散随机变量，若 $X$ 的可能取值为 $x_i$，$Y$ 的可能取值为 $y_j$，则独立性等价于

\[
P(X=x_i,Y=y_j)
=P(X=x_i)P(Y=y_j),
\qquad \forall i,j.
\]

若联合分布具有密度，则独立性等价于

\[
f_{X,Y}(x,y)=f_X(x)f_Y(y)
\]

对 Lebesgue 几乎处处的 $(x,y)$ 成立。

### 1.2 独立性在可测变换下保持

!!! success "命题 1.2（独立随机变量的函数仍独立）"

    若 $X$ 与 $Y$ 独立，$f,g:\mathbb{R}\to\mathbb{R}$ 为 Borel 可测函数，则

    \[
    U=f(X),
    \qquad
    V=g(Y)
    \]

    也相互独立。

??? proof "命题 1.2 的证明（点击展开）"

    对任意 Borel 集 $A,B$，

    \[
    \begin{aligned}
    P(U\in A,V\in B)
    &=P\bigl(X\in f^{-1}(A),Y\in g^{-1}(B)\bigr)\\
    &=P\bigl(X\in f^{-1}(A)\bigr)
      P\bigl(Y\in g^{-1}(B)\bigr)\\
    &=P(U\in A)P(V\in B).
    \end{aligned}
    \]

    其中 $f^{-1}(A)$ 与 $g^{-1}(B)$ 都是 Borel 集。$\square$

### 1.3 多个随机变量的独立性

!!! info "定义 1.3（多个随机变量相互独立）"

    随机变量 $X_1,\ldots,X_m$ 称为相互独立，如果对任意 Borel 集 $A_1,\ldots,A_m$，都有

    \[
    P(X_1\in A_1,\ldots,X_m\in A_m)
    =\prod_{j=1}^mP(X_j\in A_j).
    \]

若联合密度存在，则相互独立等价于

\[
f_{X_1,\ldots,X_m}(x_1,\ldots,x_m)
=\prod_{j=1}^mf_{X_j}(x_j)
\]

几乎处处成立。

!!! success "命题 1.4（独立子组的函数）"

    设 $X_1,\ldots,X_m$ 相互独立。若

    \[
    I,J\subseteq\{1,\ldots,m\},
    \qquad I\cap J=\varnothing,
    \]

    且 $f$ 与 $g$ 分别是相应维数的 Borel 可测函数，则

    \[
    f\bigl((X_i)_{i\in I}\bigr)
    \quad\text{与}\quad
    g\bigl((X_j)_{j\in J}\bigr)
    \]

    相互独立。

---

## 2. 随机变量的数字特征

### 2.1 期望的测度论定义

设随机变量 $X$ 的分布为 $P_X=P\circ X^{-1}$。期望可以写成关于分布的积分：

\[
EX=\int_{\mathbb{R}}x\,P_X(dx)
=\int_{\mathbb{R}}x\,dF_X(x).
\]

最后一个积分可理解为 Lebesgue-Stieltjes 积分。它统一了离散与连续情形：

\[
EX=\sum_i x_iP(X=x_i)
\]

以及

\[
EX=\int_{-\infty}^{\infty}xf_X(x)\,dx.
\]

等价地，也可以直接在原概率空间上积分：

\[
EX=\int_\Omega X(\omega)\,P(d\omega).
\]

### 2.2 从简单随机变量到一般随机变量

若 $A\in\mathcal{A}$，则示性函数 $\mathbf{1}_A$ 的期望为

\[
E\mathbf{1}_A=P(A).
\]

若 $X$ 是简单随机变量

\[
X=\sum_{i=1}^na_i\mathbf{1}_{A_i},
\]

其中 $A_1,\ldots,A_n$ 两两不交，则定义

\[
EX=\sum_{i=1}^na_iP(A_i).
\]

对非负随机变量 $X\geq0$，可以取非负简单随机变量列 $X_n$，使得

\[
0\leq X_n\uparrow X,
\]

并定义

\[
EX=\lim_{n\to\infty}EX_n.
\]

对一般实值随机变量，定义正部与负部

\[
X^+=\max\{X,0\},
\qquad
X^-=\max\{-X,0\},
\]

则

\[
X=X^+-X^-,
\qquad
|X|=X^++X^-.
\]

当 $EX^+$ 与 $EX^-$ 不同时为无穷时，可以定义

\[
EX=EX^+-EX^-.
\]

特别地，$X$ 可积当且仅当

\[
E|X|=EX^++EX^-<\infty.
\]

此时 $EX$ 是有限实数。

### 2.3 期望的基本性质

!!! success "命题 2.1（期望的基本性质）"

    设相关期望存在，则有：

    (i) 若 $a\leq X\leq b$，则

    \[
    a\leq EX\leq b.
    \]

    (ii) 对常数 $a,b$，

    \[
    E(aX+b)=aEX+b.
    \]

    (iii) 若 $X,Y$ 可积，则

    \[
    E(X+Y)=EX+EY.
    \]

    (iv) 若 $X\leq Y$ 几乎处处，则

    \[
    EX\leq EY.
    \]

!!! success "定理 2.2（单调收敛定理）"

    若 $X_n\geq0$ 且

    \[
    X_n\uparrow X
    \quad\text{几乎处处},
    \]

    则

    \[
    EX_n\uparrow EX.
    \]

### 2.4 无意识统计学家法则

!!! success "定理 2.3（LOTUS）"

    设 $Y=g(X)$。若期望存在，则

    \[
    E\{g(X)\}
    =\int_{\mathbb{R}}g(x)\,P_X(dx)
    =\int_{\mathbb{R}}g(x)\,dF_X(x).
    \]

    在离散情形，

    \[
    E\{g(X)\}
    =\sum_i g(x_i)P(X=x_i),
    \]

    在具有密度的连续情形，

    \[
    E\{g(X)\}
    =\int_{-\infty}^{\infty}g(x)f_X(x)\,dx.
    \]

这个公式无需先求 $Y=g(X)$ 的分布。

### 2.5 方差、协方差与矩

!!! info "定义 2.4（方差）"

    若 $EX^2<\infty$，则 $X$ 的方差定义为

    \[
    \operatorname{Var}(X)
    =E(X-EX)^2
    =EX^2-(EX)^2.
    \]

!!! info "定义 2.5（矩与中心矩）"

    若期望存在，则 $EX^k$ 称为 $X$ 的 $k$ 阶原点矩，$E(X-EX)^k$ 称为 $X$ 的 $k$ 阶中心矩。

!!! info "定义 2.6（协方差）"

    若 $X,Y$ 二阶可积，则

    \[
    \operatorname{Cov}(X,Y)
    =E\{(X-EX)(Y-EY)\}.
    \]

由直接展开可得

\[
\operatorname{Var}(X+Y)
=\operatorname{Var}(X)
+\operatorname{Var}(Y)
+2\operatorname{Cov}(X,Y).
\]

若 $X_1,\ldots,X_n$ 两两不相关，则

\[
\operatorname{Var}\left(\sum_{i=1}^nX_i\right)
=\sum_{i=1}^n\operatorname{Var}(X_i).
\]

特别地，相互独立且二阶可积的随机变量两两不相关，因此上式成立。

### 2.6 Chebyshev 不等式

!!! success "定理 2.7（Chebyshev 不等式）"

    若 $EX^2<\infty$，则对任意 $\varepsilon>0$，

    \[
    P(|X-EX|\geq\varepsilon)
    \leq\frac{\operatorname{Var}(X)}{\varepsilon^2}.
    \]

??? proof "定理 2.7 的证明（点击展开）"

    在事件 $\{|X-EX|\geq\varepsilon\}$ 上，

    \[
    \mathbf{1}_{\{|X-EX|\geq\varepsilon\}}
    \leq\frac{|X-EX|^2}{\varepsilon^2}.
    \]

    两边取期望，得到

    \[
    \begin{aligned}
    P(|X-EX|\geq\varepsilon)
    &=E\mathbf{1}_{\{|X-EX|\geq\varepsilon\}}\\
    &\leq\frac{E|X-EX|^2}{\varepsilon^2}\\
    &=\frac{\operatorname{Var}(X)}{\varepsilon^2}.
    \end{aligned}
    \]

    证明完毕。$\square$

### 2.7 常见分布的矩

!!! example "例 2.8（Cauchy 分布）"

    若 $X$ 服从标准 Cauchy 分布，则

    \[
    f_X(x)=\frac{1}{\pi(1+x^2)}.
    \]

    因为 $E|X|=\infty$，所以 $EX$ 不存在。

!!! example "例 2.9（正态分布）"

    若 $X\sim N(0,\sigma^2)$，则对 $k\geq0$，

    \[
    EX^{2k+1}=0,
    \]

    \[
    EX^{2k}=(2k-1)!!\,\sigma^{2k}.
    \]

!!! example "例 2.10（Poisson 分布）"

    若 $X\sim\operatorname{Poisson}(\lambda)$，则下降阶乘矩满足

    \[
    E\{X(X-1)\cdots(X-m+1)\}
    =\lambda^m.
    \]

---

## 3. 矩母函数与特征函数

### 3.1 矩母函数

!!! info "定义 3.1（矩母函数）"

    随机变量 $X$ 的矩母函数定义为

    \[
    M_X(t)=Ee^{tX},
    \]

    其中 $t$ 取使该期望有限的实数。

若 $M_X(t)$ 在 $t=0$ 的某个邻域内有限，则它在该邻域内唯一确定 $X$ 的分布，并且在相应可微条件下

\[
M_X^{(k)}(0)=EX^k.
\]

!!! warning "矩并不总能唯一确定分布"

    仅仅知道所有阶矩相同，一般不能保证两个随机变量同分布。若矩母函数在原点邻域内存在，则可以利用矩母函数的唯一性得到分布唯一性。

### 3.2 特征函数的定义

!!! info "定义 3.2（特征函数）"

    随机变量 $X$ 的特征函数定义为

    \[
    \varphi_X(t)=Ee^{itX},
    \qquad t\in\mathbb{R}.
    \]

    等价地，

    \[
    \varphi_X(t)
    =E\cos(tX)+iE\sin(tX)
    =\int_{\mathbb{R}}e^{itx}\,P_X(dx).
    \]

与矩母函数不同，特征函数对任意实值随机变量、任意 $t\in\mathbb{R}$ 都存在，因为

\[
|e^{itX}|=1.
\]

### 3.3 常见分布的特征函数

!!! example "例 3.3（Poisson 分布）"

    若 $X\sim\operatorname{Poisson}(\lambda)$，则

    \[
    \varphi_X(t)
    =\exp\{\lambda(e^{it}-1)\}.
    \]

!!! example "例 3.4（正态分布）"

    若 $X\sim N(\mu,\sigma^2)$，则

    \[
    \varphi_X(t)
    =\exp\left(i\mu t-\frac{\sigma^2t^2}{2}\right).
    \]

!!! example "例 3.5（标准 Cauchy 分布）"

    若 $X$ 服从标准 Cauchy 分布，则

    \[
    \varphi_X(t)=e^{-|t|}.
    \]

### 3.4 特征函数的基本性质

!!! success "命题 3.6（特征函数的基本性质）"

    任意特征函数都满足：

    (i)

    \[
    \varphi_X(0)=1.
    \]

    (ii)

    \[
    |\varphi_X(t)|\leq1,
    \qquad \forall t\in\mathbb{R}.
    \]

    (iii)

    \[
    \varphi_X(-t)=\overline{\varphi_X(t)}.
    \]

    (iv) $\varphi_X$ 在 $\mathbb{R}$ 上一致连续。

??? proof "一致连续性的证明（点击展开）"

    对任意 $t,h\in\mathbb{R}$，

    \[
    \begin{aligned}
    |\varphi_X(t+h)-\varphi_X(t)|
    &=\left|E\left[e^{itX}(e^{ihX}-1)\right]\right|\\
    &\leq E|e^{ihX}-1|.
    \end{aligned}
    \]

    当 $h\to0$ 时，$|e^{ihX}-1|\to0$ 几乎处处，并且

    \[
    |e^{ihX}-1|\leq2.
    \]

    由控制收敛定理，

    \[
    E|e^{ihX}-1|\to0.
    \]

    右端与 $t$ 无关，所以 $\varphi_X$ 一致连续。$\square$

### 3.5 非负定性

!!! success "命题 3.7（Bochner 非负定性）"

    对任意 $n\geq1$、复数 $a_1,\ldots,a_n$ 和实数 $t_1,\ldots,t_n$，

    \[
    \sum_{j=1}^n\sum_{k=1}^n
    a_j\overline{a_k}\,
    \varphi_X(t_j-t_k)
    \geq0.
    \]

??? proof "命题 3.7 的证明（点击展开）"

    由特征函数定义，

    \[
    \begin{aligned}
    \sum_{j,k}a_j\overline{a_k}\varphi_X(t_j-t_k)
    &=E\left(\sum_ja_je^{it_jX}\right)
      \overline{\left(\sum_ka_ke^{it_kX}\right)}\\
    &=E\left|\sum_{j=1}^na_je^{it_jX}\right|^2
    \geq0.
    \end{aligned}
    \]

    证明完毕。$\square$

### 3.6 特征函数的导数与矩

!!! success "定理 3.8（特征函数的导数）"

    若 $E|X|^k<\infty$，则 $\varphi_X$ 至少 $k$ 次可微，并且

    \[
    \varphi_X^{(j)}(t)
    =E\left[(iX)^je^{itX}\right],
    \qquad j=1,\ldots,k.
    \]

    特别地，

    \[
    \varphi_X^{(j)}(0)=i^jEX^j.
    \]

因此，当相应矩存在时，特征函数在原点附近具有展开

\[
\varphi_X(t)
=\sum_{j=0}^k\frac{i^jEX^j}{j!}t^j+o(t^k).
\]

### 3.7 独立随机变量之和

!!! success "命题 3.9（独立和的特征函数）"

    若 $X$ 与 $Y$ 独立，则

    \[
    \varphi_{X+Y}(t)
    =\varphi_X(t)\varphi_Y(t).
    \]

    更一般地，若 $X_1,\ldots,X_n$ 相互独立，则

    \[
    \varphi_{\sum_{j=1}^nX_j}(t)
    =\prod_{j=1}^n\varphi_{X_j}(t).
    \]

??? proof "命题 3.9 的证明（点击展开）"

    由独立性，

    \[
    \begin{aligned}
    \varphi_{X+Y}(t)
    &=E e^{it(X+Y)}\\
    &=E(e^{itX}e^{itY})\\
    &=Ee^{itX}\,Ee^{itY}.
    \end{aligned}
    \]

    多个随机变量的情形由归纳法得到。$\square$

!!! example "例 3.10（独立正态变量之和）"

    若

    \[
    X_i\sim N(\mu_i,\sigma_i^2),
    \qquad i=1,\ldots,n,
    \]

    且 $X_1,\ldots,X_n$ 相互独立，则

    \[
    \sum_{i=1}^nX_i
    \sim N\left(\sum_{i=1}^n\mu_i,
    \sum_{i=1}^n\sigma_i^2\right).
    \]

    这是因为相应特征函数的乘积为

    \[
    \exp\left{
    it\sum_{i=1}^n\mu_i
    -\frac{t^2}{2}\sum_{i=1}^n\sigma_i^2
    \right},
    \]

    它正是上述正态分布的特征函数。

---

## 4. 特征函数与分布

### 4.1 唯一性定理

!!! success "定理 4.1（特征函数的唯一性）"

    特征函数唯一确定概率分布。也就是说，

    \[
    X\overset d=Y
    \quad\Longleftrightarrow\quad
    \varphi_X(t)=\varphi_Y(t),
    \qquad \forall t\in\mathbb{R}.
    \]

### 4.2 反演公式

!!! success "定理 4.2（Lévy 反演公式）"

    设 $F$ 是 $X$ 的分布函数，$x_1<x_2$，并且 $x_1,x_2$ 都是 $F$ 的连续点，则

    \[
    F(x_2)-F(x_1)
    =\lim_{T\to\infty}\frac{1}{2\pi}
    \int_{-T}^T
    \frac{e^{-itx_1}-e^{-itx_2}}{it}
    \varphi_X(t)\,dt.
    \]

该公式说明，可以从特征函数恢复分布函数在连续点之间的概率质量，从而得到唯一性定理。

!!! note "分布函数的不连续点"

    分布函数在点 $x$ 处的跳跃大小为

    \[
    F(x)-F(x-)=P(X=x).
    \]

    因此分布函数的不连续点至多可数。

### 4.3 密度的 Fourier 反演

在适当的可积条件下，若 $X$ 具有密度 $f_X$，则

\[
\varphi_X(t)
=\int_{-\infty}^{\infty}e^{itx}f_X(x)\,dx,
\]

并可由 Fourier 反演得到

\[
f_X(x)
=\frac{1}{2\pi}
\int_{-\infty}^{\infty}e^{-itx}\varphi_X(t)\,dt.
\]

!!! example "例 4.3（对称两点分布）"

    若

    \[
    P(X=1)=P(X=-1)=\frac12,
    \]

    则

    \[
    \varphi_X(t)
    =\frac12e^{it}+\frac12e^{-it}
    =\cos t.
    \]

---

## 5. 随机变量序列的收敛

### 5.1 依概率收敛

!!! info "定义 5.1（依概率收敛）"

    若对任意 $\varepsilon>0$，

    \[
    P(|X_n-X|>\varepsilon)\to0,
    \]

    则称 $X_n$ 依概率收敛到 $X$，记为

    \[
    X_n\xrightarrow{P}X.
    \]

### 5.2 依分布收敛

!!! info "定义 5.2（依分布收敛）"

    设 $F_n$ 与 $F$ 分别是 $X_n$ 与 $X$ 的分布函数。如果对 $F$ 的每个连续点 $x$，都有

    \[
    F_n(x)\to F(x),
    \]

    则称 $X_n$ 依分布收敛到 $X$，记为

    \[
    X_n\xrightarrow{d}X.
    \]

!!! warning "为什么只要求在极限分布的连续点收敛"

    令 $x_n=1/n$，并把 $X_n\equiv x_n$、$X\equiv0$ 看作退化随机变量。显然

    \[
    X_n\to X
    \]

    逐点成立，因而 $X_n\xrightarrow{d}X$。但是在 $x=0$ 处，

    \[
    F_n(0)=0,
    \qquad
    F(0)=1.
    \]

    所以分布函数在极限分布的不连续点未必收敛。

依概率收敛蕴含依分布收敛，但反之一般不成立。

### 5.3 Lévy 连续性定理

!!! success "定理 5.3（Lévy 连续性定理）"

    若 $X_n$ 与 $X$ 的特征函数分别为 $\varphi_n$ 和 $\varphi$，则

    \[
    X_n\xrightarrow{d}X
    \quad\Longleftrightarrow\quad
    \varphi_n(t)\to\varphi(t),
    \qquad \forall t\in\mathbb{R}.
    \]

更一般地，如果 $\varphi_n(t)$ 逐点收敛到函数 $\varphi(t)$，且 $\varphi$ 在 $0$ 处连续，则 $\varphi$ 是某个概率分布的特征函数，并且相应的分布弱收敛到该分布。

---

## 6. Bernoulli 和的极限定理

设事件 $A$ 的概率为 $p$，令

\[
\xi_i(\omega)=
\begin{cases}
1,&\omega_i\in A,\\
0,&\omega_i\notin A,
\end{cases}
\]

其中 $\xi_1,\ldots,\xi_n$ 相互独立，并记

\[
S_n=\sum_{i=1}^n\xi_i.
\]

则

\[
S_n\sim\operatorname{Binomial}(n,p),
\]

即

\[
P(S_n=k)=\binom nkp^k(1-p)^{n-k}.
\]

### 6.1 Bernoulli 弱大数定律

!!! success "定理 6.1（Bernoulli 弱大数定律）"

    对固定的 $0<p<1$，

    \[
    \frac{S_n}{n}\xrightarrow{P}p.
    \]

??? proof "定理 6.1 的证明（点击展开）"

    因为

    \[
    E\left(\frac{S_n}{n}\right)=p,
    \qquad
    \operatorname{Var}\left(\frac{S_n}{n}\right)
    =\frac{p(1-p)}{n},
    \]

    由 Chebyshev 不等式，对任意 $\varepsilon>0$，

    \[
    P\left(\left|\frac{S_n}{n}-p\right|>\varepsilon\right)
    \leq\frac{p(1-p)}{n\varepsilon^2}
    \to0.
    \]

    证明完毕。$\square$

### 6.2 Poisson 极限定理

!!! success "定理 6.2（Poisson 极限定理）"

    设

    \[
    S_n\sim\operatorname{Binomial}(n,p_n),
    \qquad np_n\to\lambda\in(0,\infty).
    \]

    则对每个固定的非负整数 $k$，

    \[
    P(S_n=k)
    \to e^{-\lambda}\frac{\lambda^k}{k!}.
    \]

    因此

    \[
    S_n\xrightarrow{d}\operatorname{Poisson}(\lambda).
    \]

### 6.3 De Moivre-Laplace 定理

令 $q=1-p$。

!!! success "定理 6.3（De Moivre-Laplace 中心极限定理）"

    若 $S_n\sim\operatorname{Binomial}(n,p)$，其中 $0<p<1$ 固定，则

    \[
    \frac{S_n-np}{\sqrt{npq}}
    \xrightarrow{d}N(0,1).
    \]

    等价地，对任意 $x\in\mathbb{R}$，

    \[
    P\left(
    \frac{S_n-np}{\sqrt{npq}}\leq x
    \right)
    \to
    \Phi(x)
    =\frac{1}{\sqrt{2\pi}}
    \int_{-\infty}^xe^{-t^2/2}\,dt.
    \]

该定理可以利用 Stirling 公式

\[
n!\sim\sqrt{2\pi n}\left(\frac ne\right)^n
\]

从二项分布概率质量函数出发证明。

---

## 7. 中心极限定理的特征函数证明

!!! success "定理 7.1（Lindeberg-Lévy 中心极限定理）"

    设 $X_1,X_2,\ldots$ 独立同分布，并且

    \[
    EX_1=\mu,
    \qquad
    \operatorname{Var}(X_1)=\sigma^2\in(0,\infty).
    \]

    令

    \[
    S_n=\sum_{i=1}^nX_i.
    \]

    则

    \[
    \frac{S_n-n\mu}{\sigma\sqrt n}
    \xrightarrow{d}N(0,1).
    \]

??? proof "特征函数证明（点击展开）"

    令

    \[
    Y_i=\frac{X_i-\mu}{\sigma},
    \]

    则 $EY_i=0$、$EY_i^2=1$。记 $Y_i$ 的特征函数为 $\varphi_Y$。由原点附近的二阶展开，

    \[
    \varphi_Y(t)
    =1-\frac{t^2}{2}+o(t^2),
    \qquad t\to0.
    \]

    由独立性，标准化和的特征函数为

    \[
    \varphi_{n}(t)
    =\left[\varphi_Y\left(\frac{t}{\sqrt n}\right)\right]^n.
    \]

    因此

    \[
    \begin{aligned}
    \varphi_n(t)
    &=\left[
    1-\frac{t^2}{2n}+o\left(\frac1n\right)
    \right]^n\\
    &\to e^{-t^2/2}.
    \end{aligned}
    \]

    函数 $e^{-t^2/2}$ 是标准正态分布的特征函数。由 Lévy 连续性定理，

    \[
    \frac{S_n-n\mu}{\sigma\sqrt n}
    \xrightarrow{d}N(0,1).
    \]

    证明完毕。$\square$

---

## 8. 本章小结

本章的主要结论如下：

1. 随机变量的独立性可以用联合分布函数、概率质量函数或联合密度的乘积形式刻画，并在分别进行可测变换后保持。
2. 期望可以统一表示为关于分布的 Lebesgue-Stieltjes 积分；可积性等价于 $E|X|<\infty$。
3. 方差、协方差与矩刻画随机变量的数字特征，Chebyshev 不等式把方差转化为尾概率上界。
4. 特征函数总是存在，并且唯一确定分布；独立随机变量之和的特征函数等于各自特征函数的乘积。
5. Lévy 连续性定理把特征函数的逐点收敛与依分布收敛联系起来。
6. Bernoulli 弱大数定律、Poisson 极限定理、De Moivre-Laplace 定理和 Lindeberg-Lévy 中心极限定理分别描述频率稳定、稀有事件和标准化和的极限行为。
