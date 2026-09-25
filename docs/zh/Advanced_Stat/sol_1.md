---
tags:
  - 高等数理统计
  - 统计推断
  - 课后习题
---

# 高等数理统计第一次作业


## 习题 1：无偏性与平方误差风险

**题目。** 设

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(\theta).
\]

证明 $\overline X$ 是 $\theta$ 的无偏估计量，并计算其在平方误差损失下的风险。

??? success "解答（点击展开）"

    由于

    \[
    E_\theta(X_i)=\theta,
    \]

    因此

    \[
    E_\theta(\overline X)
    =E_\theta\left(\frac{1}{n}\sum_{i=1}^nX_i\right)
    =\frac{1}{n}\sum_{i=1}^nE_\theta(X_i)
    =\theta.
    \]

    所以 $\overline X$ 是 $\theta$ 的无偏估计量。

    在平方误差损失下，

    \[
    L(\theta,a)=(a-\theta)^2,
    \]

    $\overline X$ 的风险为

    \[
    \begin{aligned}
    R(\theta,\overline X)
    &=E_\theta(\overline X-\theta)^2\\
    &=\operatorname{Var}_\theta(\overline X)
      +\{E_\theta(\overline X)-\theta\}^2\\
    &=\operatorname{Var}_\theta(\overline X).
    \end{aligned}
    \]

    由独立性，

    \[
    \begin{aligned}
    \operatorname{Var}_\theta(\overline X)
    &=\frac{1}{n^2}\sum_{i=1}^n
    \operatorname{Var}_\theta(X_i)\\
    &=\frac{1}{n^2}n\theta(1-\theta)\\
    &=\frac{\theta(1-\theta)}{n}.
    \end{aligned}
    \]

    因此，

    \[
    \boxed{R(\theta,\overline X)
    =\frac{\theta(1-\theta)}{n}}.
    \]

---

## 习题 2：Student (t) 分布

**题目。** 设 $X$ 服从自由度为 $p$ 的 Student $t$ 分布。

(a) 推导 $X$ 的均值和方差。

(b) 证明 $X^2$ 服从自由度分别为 $1$ 和 $p$ 的 $F$ 分布。

(c) 记 $f_X(x\mid p)$ 为 $X$ 的密度函数。证明对每个固定的 $x$，

\[
\lim_{p\to\infty}f_X(x\mid p)
=\frac{1}{\sqrt{2\pi}}e^{-x^2/2}.
\]

??? success "解答（点击展开）"

    使用标准表示

    \[
    X=\frac{Z}{\sqrt{W/p}},
    \qquad
    Z\sim N(0,1),
    \qquad
    W\sim\chi_p^2,
    \qquad
    Z\perp W.
    \]

    **(a) 均值与方差**

    当 $r<p/2$ 且 $W\sim\chi_p^2$ 时，

    \[
    E(W^{-r})
    =2^{-r}\frac{\Gamma(p/2-r)}{\Gamma(p/2)}.
    \]

    特别地，

    \[
    E(W^{-1/2})<\infty
    \iff p>1,
    \qquad
    E(W^{-1})<\infty
    \iff p>2.
    \]

    因此，当 $p>1$ 时，$X$ 可积，并且

    \[
    E(X)
    =\sqrt{p}\,E(Z)E(W^{-1/2})
    =0.
    \]

    所以，

    \[
    \boxed{E(X)=0\quad(p>1)}.
    \]

    当 $p\leq1$ 时，均值不存在。当 $p>2$ 时，

    \[
    \begin{aligned}
    E(X^2)
    &=pE(Z^2)E(W^{-1})\\
    &=p\cdot1\cdot
    \frac{1}{2}\frac{\Gamma(p/2-1)}{\Gamma(p/2)}\\
    &=p\cdot\frac{1}{p-2}\\
    &=\frac{p}{p-2}.
    \end{aligned}
    \]

    由于 $E(X)=0$，

    \[
    \boxed{\operatorname{Var}(X)
    =\frac{p}{p-2}\quad(p>2)}.
    \]

    当 $1<p\leq2$ 时，均值存在但方差为无穷；当 $p\leq1$ 时，由于均值不存在，通常意义下的方差也没有定义。

    **(b) (X^2) 的分布**

    由于 $Z^2\sim\chi_1^2$ 且 $Z^2\perp W$，

    \[
    X^2
    =\frac{Z^2}{W/p}
    =\frac{Z^2/1}{W/p}.
    \]

    根据 $F$ 分布的定义，

    \[
    \boxed{X^2\sim F_{1,p}}.
    \]

    **(c) 密度函数的极限**

    $t_p$ 分布的密度为

    \[
    f_X(x\mid p)
    =\frac{\Gamma\{(p+1)/2\}}
    {\sqrt{p\pi}\,\Gamma(p/2)}
    \left(1+\frac{x^2}{p}\right)^{-(p+1)/2}.
    \]

    记

    \[
    A_p
    =\frac{\Gamma\{(p+1)/2\}}
    {\sqrt{p\pi}\,\Gamma(p/2)},
    \qquad
    B_p
    =\left(1+\frac{x^2}{p}\right)^{-(p+1)/2}.
    \]

    令 $z=p/2$。由标准的 Gamma 函数比值渐近公式

    \[
    \frac{\Gamma(z+a)}{\Gamma(z+b)}
    \sim z^{a-b},
    \qquad z\to\infty,
    \]

    可得

    \[
    A_p
    \sim\frac{(p/2)^{1/2}}{\sqrt{p\pi}}
    =\frac{1}{\sqrt{2\pi}}.
    \]

    另一方面，

    \[
    \begin{aligned}
    \log B_p
    &=-\frac{p+1}{2}
    \log\left(1+\frac{x^2}{p}\right)\\
    &=-\frac{1}{2}\left(1+\frac{1}{p}\right)
    p\log\left(1+\frac{x^2}{p}\right)\\
    &\longrightarrow-\frac{x^2}{2},
    \end{aligned}
    \]

    这是因为

    \[
    p\log\left(1+\frac{x^2}{p}\right)
    \longrightarrow x^2.
    \]

    因此

    \[
    B_p\longrightarrow e^{-x^2/2}.
    \]

    所以对每个固定的 $x$，

    \[
    \boxed{\lim_{p\to\infty}f_X(x\mid p)
    =\frac{1}{\sqrt{2\pi}}e^{-x^2/2}}.
    \]

---

## 习题 3：样本方差恒等式及其与样本均值的协方差

**题目。** 设 $X_1,\ldots,X_n$ 为随机样本，$\overline X$ 与 $S^2$ 按通常方式定义。

(a) 证明

\[
S^2
=\frac{1}{2n(n-1)}
\sum_{i=1}^n\sum_{j=1}^n(X_i-X_j)^2.
\]

(b) 假设 $X_i$ 具有有限四阶矩。令

\[
\theta_1=E(X_i),
\qquad
\theta_j=E\{(X_i-\theta_1)^j\},
\qquad j=2,3,4.
\]

用 $\theta_1,\ldots,\theta_4$ 表示 $\operatorname{Cov}(\overline X,S^2)$。在什么条件下 $\operatorname{Cov}(\overline X,S^2)=0$？

??? success "解答（点击展开）"

    **(a) 成对差表示**

    展开二重求和可得

    \[
    \begin{aligned}
    \sum_{i=1}^n\sum_{j=1}^n(X_i-X_j)^2
    &=\sum_{i,j}(X_i^2-2X_iX_j+X_j^2)\\
    &=2n\sum_{i=1}^nX_i^2
      -2\left(\sum_{i=1}^nX_i\right)^2\\
    &=2n\left(\sum_{i=1}^nX_i^2-n\overline X^2\right)\\
    &=2n\sum_{i=1}^n(X_i-\overline X)^2.
    \end{aligned}
    \]

    由于

    \[
    S^2=\frac{1}{n-1}
    \sum_{i=1}^n(X_i-\overline X)^2,
    \]

    因此

    \[
    \boxed{S^2
    =\frac{1}{2n(n-1)}
    \sum_{i=1}^n\sum_{j=1}^n(X_i-X_j)^2}.
    \]

    **(b) 与样本均值的协方差**

    令

    \[
    Y_i=X_i-\theta_1,
    \qquad
    \overline Y=\overline X-\theta_1.
    \]

    则

    \[
    E(Y_i)=0,
    \qquad
    E(Y_i^j)=\theta_j,
    \qquad j=2,3,4,
    \]

    并且

    \[
    S^2
    =\frac{1}{n-1}
    \left(\sum_{i=1}^nY_i^2-n\overline Y^2\right).
    \]

    由于 $E(\overline Y)=0$，

    \[
    \operatorname{Cov}(\overline X,S^2)
    =\operatorname{Cov}(\overline Y,S^2)
    =E(\overline YS^2).
    \]

    首先，

    \[
    \begin{aligned}
    E\left(\overline Y\sum_{i=1}^nY_i^2\right)
    &=\frac{1}{n}
    \sum_{j=1}^n\sum_{i=1}^nE(Y_jY_i^2)\\
    &=\frac{1}{n}\sum_{i=1}^nE(Y_i^3)\\
    &=\theta_3,
    \end{aligned}
    \]

    其中所有 $i\neq j$ 的项都因独立性和中心化而为零。另外，

    \[
    \begin{aligned}
    E(\overline Y^3)
    &=\frac{1}{n^3}E\left(\sum_{i=1}^nY_i\right)^3\\
    &=\frac{1}{n^3}\sum_{i=1}^nE(Y_i^3)\\
    &=\frac{\theta_3}{n^2},
    \end{aligned}
    \]

    因为每个混合项都含有一个均值为零的因子。因此，

    \[
    \begin{aligned}
    \operatorname{Cov}(\overline X,S^2)
    &=\frac{1}{n-1}
    \left\{\theta_3-n\left(\frac{\theta_3}{n^2}\right)\right\}\\
    &=\frac{\theta_3}{n}.
    \end{aligned}
    \]

    所以，

    \[
    \boxed{\operatorname{Cov}(\overline X,S^2)
    =\frac{\theta_3}{n}}.
    \]

    从而，

    \[
    \operatorname{Cov}(\overline X,S^2)=0
    \iff\theta_3=0.
    \]

    特别地，任何关于其均值对称的分布都满足 $\theta_3=0$，因此 $\overline X$ 与 $S^2$ 不相关。对称性是充分条件，但不是必要条件。

---

## 习题 4：|X| 的充分性

**题目。** 设 $X$ 是来自 $N(0,\sigma^2)$ 总体的一个观测。$|X|$ 是否为 $\sigma^2$ 的充分统计量？

??? success "解答（点击展开）"

    密度函数为

    \[
    f_\sigma(x)
    =\frac{1}{\sqrt{2\pi}\sigma}
    \exp\left(-\frac{x^2}{2\sigma^2}\right).
    \]

    由于 $x^2=|x|^2$，

    \[
    f_\sigma(x)
    =\underbrace{
    \frac{1}{\sqrt{2\pi}\sigma}
    \exp\left(-\frac{|x|^2}{2\sigma^2}\right)
    }_{g_\sigma(|x|)}
    \underbrace{1}_{h(x)}.
    \]

    支撑集为 $\mathbb{R}$，且不依赖于 $\sigma^2$。因此，由 Fisher-Neyman 因子分解定理，

    \[
    \boxed{|X|\text{ 是 }\sigma^2\text{ 的充分统计量}}.
    \]

---

## 习题 5：指数族表示与最小充分性

**题目。** 设 $X_1,\ldots,X_n$ 的密度为

\[
f_{\mu,\sigma}(x)
=\exp\left\{-\left(\frac{x-\mu}{\sigma}\right)^4
-\xi(\mu,\sigma)\right\},
\qquad
\mu\in\mathbb{R},
\quad \sigma>0.
\]

证明其联合分布族是指数族，并找出一个最小充分统计量。

??? success "解答（点击展开）"

    首先展开

    \[
    -\left(\frac{x-\mu}{\sigma}\right)^4
    =-\frac{x^4}{\sigma^4}
    +\frac{4\mu x^3}{\sigma^4}
    -\frac{6\mu^2x^2}{\sigma^4}
    +\frac{4\mu^3x}{\sigma^4}
    -\frac{\mu^4}{\sigma^4}.
    \]

    因此联合密度为

    \[
    \begin{aligned}
    f_{\mu,\sigma}(\mathbf{x})
    &=\prod_{i=1}^nf_{\mu,\sigma}(x_i)\\
    &=\exp\Biggl\{
    -\frac{1}{\sigma^4}\sum_{i=1}^nx_i^4
    +\frac{4\mu}{\sigma^4}\sum_{i=1}^nx_i^3
    -\frac{6\mu^2}{\sigma^4}\sum_{i=1}^nx_i^2\\
    &\qquad
    +\frac{4\mu^3}{\sigma^4}\sum_{i=1}^nx_i
    -n\left(\frac{\mu^4}{\sigma^4}
    +\xi(\mu,\sigma)\right)
    \Biggr\}.
    \end{aligned}
    \]

    因而该分布族是一个曲指数族，其自然统计量向量为

    \[
    T(\mathbf{X})
    =\left(
    \sum_{i=1}^nX_i,
    \sum_{i=1}^nX_i^2,
    \sum_{i=1}^nX_i^3,
    \sum_{i=1}^nX_i^4
    \right).
    \]

    相应地，可以取自然参数向量

    \[
    \eta(\mu,\sigma)
    =\left(
    \frac{4\mu^3}{\sigma^4},
    -\frac{6\mu^2}{\sigma^4},
    \frac{4\mu}{\sigma^4},
    -\frac{1}{\sigma^4}
    \right).
    \]

    为证明最小充分性，使用似然比判别准则。对两个样本点 $\mathbf{x}$ 与 $\mathbf{y}$，

    \[
    \begin{aligned}
    \log\frac{f_{\mu,\sigma}(\mathbf{x})}
    {f_{\mu,\sigma}(\mathbf{y})}
    &=-\frac{1}{\sigma^4}
    \sum_{i=1}^n\{(x_i-\mu)^4-(y_i-\mu)^4\}\\
    &=-\frac{1}{\sigma^4}
    \left(D_4-4\mu D_3+6\mu^2D_2-4\mu^3D_1\right),
    \end{aligned}
    \]

    其中

    \[
    D_k=\sum_{i=1}^nx_i^k-
    \sum_{i=1}^ny_i^k,
    \qquad k=1,2,3,4.
    \]

    似然比与 $(\mu,\sigma)$ 无关，当且仅当

    \[
    D_4-4\mu D_3+6\mu^2D_2-4\mu^3D_1=0
    \]

    对每个 $\mu\in\mathbb{R}$ 都成立。一个多项式恒等于零，当且仅当其每个系数都为零，因此

    \[
    D_1=D_2=D_3=D_4=0.
    \]

    这等价于

    \[
    T(\mathbf{x})=T(\mathbf{y}).
    \]

    因此，由最小充分性的似然比刻画，

    \[
    \boxed{
    T(\mathbf{X})
    =\left(
    \sum_{i=1}^nX_i,
    \sum_{i=1}^nX_i^2,
    \sum_{i=1}^nX_i^3,
    \sum_{i=1}^nX_i^4
    \right)
    }
    \]

    是 $(\mu,\sigma)$ 的最小充分统计量。

---

## 习题 6：最小充分但不完备

**题目。** 对来自 $\operatorname{Uniform}(\theta,\theta+1)$ 的样本，证明最小充分统计量 $(X_{(1)},X_{(n)})$ 不完备。

??? success "解答（点击展开）"

    联合密度为

    \[
    \begin{aligned}
    f_\theta(\mathbf{x})
    &=\prod_{i=1}^n
    \mathbf{1}_{\{\theta<x_i<\theta+1\}}\\
    &=\mathbf{1}_{\{\theta<X_{(1)},\ X_{(n)}<\theta+1\}}\\
    &=\mathbf{1}_{\{X_{(n)}-1<\theta<X_{(1)}\}}.
    \end{aligned}
    \]

    因此，由 Fisher-Neyman 因子分解定理，

    \[
    T=(X_{(1)},X_{(n)})
    \]

    是充分统计量。

    为证明最小充分性，对样本点 $\mathbf{x}$ 定义

    \[
    I_\mathbf{x}=(x_{(n)}-1,x_{(1)}).
    \]

    在所有支撑集的并上，也就是对满足

    \[
    x_{(n)}-x_{(1)}\leq1,
    \]

    的样本点，有

    \[
    f_\theta(\mathbf{x})
    =\mathbf{1}_{\{\theta\in I_\mathbf{x}\}}.
    \]

    对两个样本点 $\mathbf{x}$ 与 $\mathbf{y}$，似然比与 $\theta$ 无关，当且仅当

    \[
    I_\mathbf{x}=I_\mathbf{y},
    \]

    这又等价于

    \[
    x_{(1)}=y_{(1)},
    \qquad
    x_{(n)}=y_{(n)}.
    \]

    因此，似然比判别准则表明 $(X_{(1)},X_{(n)})$ 是最小充分统计量。

    接下来证明它不完备。先设 $n\geq2$，考虑极差

    \[
    R=X_{(n)}-X_{(1)},
    \]

    它是 $T$ 的可测函数。写成

    \[
    X_i=\theta+U_i,
    \qquad
    U_i\overset{\mathrm{iid}}{\sim}\operatorname{Uniform}(0,1).
    \]

    则

    \[
    R=U_{(n)}-U_{(1)},
    \]

    所以其分布不依赖于 $\theta$。此外，

    \[
    E\{U_{(n)}\}=\frac{n}{n+1},
    \qquad
    E\{U_{(1)}\}=\frac{1}{n+1},
    \]

    从而

    \[
    E_\theta(R)=\frac{n-1}{n+1}.
    \]

    定义

    \[
    g(T)
    =X_{(n)}-X_{(1)}-\frac{n-1}{n+1}.
    \]

    则对每个 $\theta$，

    \[
    E_\theta\{g(T)\}=0,
    \]

    但由于极差不是退化随机变量，$g(T)$ 并不几乎处处等于零。因此 $T$ 不完备。综上，

    \[
    \boxed{(X_{(1)},X_{(n)})
    \text{ 是最小充分但不完备的统计量}}.
    \]

    当 $n=1$ 时，可以另取

    \[
    g(X_1)=\sin(2\pi X_1).
    \]

    来说明不完备性。事实上，

    \[
    E_\theta\{g(X_1)\}
    =\int_\theta^{\theta+1}\sin(2\pi x)\,dx
    =0
    \]

    对每个 $\theta$ 都成立，但 $g(X_1)$ 并不几乎处处等于零。
