# 第二章：矩阵范数、扰动理论与奇异值分解

本章首先继续讨论常用矩阵范数之间的关系，并介绍谱范数与 Frobenius 范数的正交不变性。随后利用 Neumann 级数研究逆矩阵在扰动下的变化。最后系统介绍奇异值分解（SVD），并讨论奇异值的基本性质、扰动估计以及低秩逼近。

---

## 1. 常用矩阵范数的性质

设 $A=(a_{ij})\in\mathbb{R}^{m\times n}$。本章主要使用以下四种矩阵范数：

\[
\|A\|_F
=\left(\sum_{i=1}^m\sum_{j=1}^n|a_{ij}|^2\right)^{1/2},
\]

\[
\|A\|_1
=\max_{1\leq j\leq n}\sum_{i=1}^m|a_{ij}|,
\]

\[
\|A\|_\infty
=\max_{1\leq i\leq m}\sum_{j=1}^n|a_{ij}|,
\]

以及诱导 $2$-范数（谱范数）

\[
\|A\|_2
=\max_{x\neq0}\frac{\|Ax\|_2}{\|x\|_2}.
\]

### 1.1 有限维空间中范数的等价性

!!! success "定理 1.1（有限维空间中范数等价）"

    在有限维线性空间 $\mathbb{R}^{m\times n}$ 上，任意两个范数都是等价的。也就是说，对任意两个矩阵范数 $\|\cdot\|_a$ 和 $\|\cdot\|_b$，存在仅依赖于 $m,n$ 的常数 $c,C>0$，使得

    \[
    c\|A\|_a\leq\|A\|_b\leq C\|A\|_a,
    \qquad A\in\mathbb{R}^{m\times n}.
    \]

因此，在有限维空间中，一个矩阵序列是否收敛并不依赖于具体选择的范数。不过，不同范数给出的数值大小和误差界可能不同。

### 1.2 谱范数与 Frobenius 范数

!!! success "定理 1.2（谱范数与 Frobenius 范数的关系）"

    对任意 $A\in\mathbb{R}^{m\times n}$，有

    \[
    \|A\|_2
    \leq\|A\|_F
    \leq\sqrt{\min\{m,n\}}\,\|A\|_2.
    \]

??? proof "定理 1.2 的证明（点击展开）"

    由于 $A^TA$ 是半正定矩阵，设其非零特征值为

    \[
    \lambda_1\geq\lambda_2\geq\cdots\geq\lambda_r>0,
    \qquad r=\operatorname{rank}(A).
    \]

    注意到

    \[
    \|A\|_2^2=\lambda_1,
    \qquad
    \|A\|_F^2=\operatorname{tr}(A^TA)=\sum_{i=1}^r\lambda_i.
    \]

    因此

    \[
    \lambda_1
    \leq\sum_{i=1}^r\lambda_i
    \leq r\lambda_1
    \leq\min\{m,n\}\lambda_1.
    \]

    两端开平方即得结论。$\square$

另一个直接证明上界

\[
\|A\|_2\leq\|A\|_F
\]

的方法是使用 Cauchy-Schwarz 不等式。对任意 $\|x\|_2=1$，

\[
\begin{aligned}
\|Ax\|_2^2
&=\sum_{i=1}^m\left(\sum_{j=1}^na_{ij}x_j\right)^2\\
&\leq\sum_{i=1}^m
\left(\sum_{j=1}^na_{ij}^2\right)
\left(\sum_{j=1}^nx_j^2\right)\\
&=\|A\|_F^2.
\end{aligned}
\]

对单位向量取上确界即可。

### 1.3 转置与常用范数

!!! success "命题 1.3（转置下的不变性与对偶关系）"

    对任意 $A\in\mathbb{R}^{m\times n}$，有

    \[
    \|A^T\|_F=\|A\|_F,
    \qquad
    \|A^T\|_2=\|A\|_2,
    \]

    以及

    \[
    \|A\|_1=\|A^T\|_\infty,
    \qquad
    \|A\|_\infty=\|A^T\|_1.
    \]

第一个等式直接来自 Frobenius 范数的定义；第二个等式来自 $A^TA$ 与 $AA^T$ 具有相同的非零特征值；后两个等式来自行和与列和在转置下的互换。

### 1.4 谱范数与 $1$-范数、$\infty$-范数

!!! success "定理 1.4（常用诱导范数之间的关系）"

    对任意 $A\in\mathbb{R}^{m\times n}$，有

    \[
    \frac{1}{\sqrt{n}}\|A\|_\infty
    \leq\|A\|_2
    \leq\sqrt{m}\,\|A\|_\infty,
    \]

    以及

    \[
    \frac{1}{\sqrt{m}}\|A\|_1
    \leq\|A\|_2
    \leq\sqrt{n}\,\|A\|_1.
    \]

??? proof "定理 1.4 的证明（点击展开）"

    对任意 $x\in\mathbb{R}^n$，由向量范数之间的关系可得

    \[
    \|Ax\|_2
    \leq\sqrt{m}\|Ax\|_\infty
    \leq\sqrt{m}\|A\|_\infty\|x\|_\infty
    \leq\sqrt{m}\|A\|_\infty\|x\|_2.
    \]

    因而

    \[
    \|A\|_2\leq\sqrt{m}\|A\|_\infty.
    \]

    另一方面，记 $e_i$ 为第 $i$ 个标准基向量。对任意 $i=1,\ldots,m$，

    \[
    \sum_{j=1}^na_{ij}^2
    =\|A^Te_i\|_2^2
    \leq\|A^T\|_2^2
    =\|A\|_2^2.
    \]

    再由 Cauchy-Schwarz 不等式，

    \[
    \left(\sum_{j=1}^n|a_{ij}|\right)^2
    \leq n\sum_{j=1}^na_{ij}^2
    \leq n\|A\|_2^2.
    \]

    对 $i$ 取最大值得到

    \[
    \|A\|_\infty\leq\sqrt{n}\|A\|_2.
    \]

    这就证明了第一组不等式。对 $A^T$ 应用第一组不等式，并使用

    \[
    \|A^T\|_2=\|A\|_2,
    \qquad
    \|A^T\|_\infty=\|A\|_1,
    \]

    即得第二组不等式。$\square$

### 1.5 谱范数的特征值刻画

!!! success "定理 1.5（谱范数的取值向量）"

    对任意 $A\in\mathbb{R}^{m\times n}$，存在单位向量 $z\in\mathbb{R}^n$，使得

    \[
    A^TAz=\|A\|_2^2z.
    \]

    换言之，$\|A\|_2^2$ 是 $A^TA$ 的最大特征值。

??? proof "定理 1.5 的证明（点击展开）"

    在单位球面上定义 Rayleigh 商

    \[
    q(x)=x^TA^TAx=\|Ax\|_2^2,
    \qquad \|x\|_2=1.
    \]

    单位球面是紧集，$q$ 是连续函数，因此最大值能够取到。另一方面，由实对称矩阵的谱定理，设 $A^TA$ 的特征值满足

    \[
    \lambda_1\geq\lambda_2\geq\cdots\geq\lambda_n\geq0,
    \]

    对应一组标准正交特征向量 $v_1,\ldots,v_n$。若

    \[
    x=\sum_{i=1}^n\alpha_iv_i,
    \qquad
    \sum_{i=1}^n\alpha_i^2=1,
    \]

    则

    \[
    x^TA^TAx
    =\sum_{i=1}^n\lambda_i\alpha_i^2
    \leq\lambda_1.
    \]

    当 $x=v_1$ 时等号成立，所以

    \[
    \|A\|_2^2=\lambda_1,
    \]

    并可取 $z=v_1$。$\square$

!!! success "推论 1.6"

    对任意 $A\in\mathbb{R}^{m\times n}$，

    \[
    \|A\|_2
    \leq\sqrt{\|A\|_1\|A\|_\infty}.
    \]

??? proof "推论 1.6 的证明（点击展开）"

    取非零向量 $z$ 使得

    \[
    A^TAz=\|A\|_2^2z.
    \]

    对两边取 $1$-范数，并利用诱导范数的相容性，得到

    \[
    \|A\|_2^2\|z\|_1
    =\|A^TAz\|_1
    \leq\|A^T\|_1\|A\|_1\|z\|_1.
    \]

    由于 $\|A^T\|_1=\|A\|_\infty$，约去 $\|z\|_1$ 后即得

    \[
    \|A\|_2^2\leq\|A\|_1\|A\|_\infty.
    \]

    两端开平方即可。$\square$

### 1.6 正交不变性

!!! success "定理 1.7（正交不变性）"

    设 $A\in\mathbb{R}^{m\times n}$，$Q\in\mathbb{R}^{m\times m}$ 和 $P\in\mathbb{R}^{n\times n}$ 均为正交矩阵，则

    \[
    \|QAP^T\|_2=\|A\|_2,
    \qquad
    \|QAP^T\|_F=\|A\|_F.
    \]

??? proof "定理 1.7 的证明（点击展开）"

    对谱范数，由正交矩阵保持 Euclidean 范数可得

    \[
    \begin{aligned}
    \|QAP^T\|_2
    &=\max_{\|x\|_2=1}\|QAP^Tx\|_2\\
    &=\max_{\|x\|_2=1}\|AP^Tx\|_2\\
    &=\max_{\|y\|_2=1}\|Ay\|_2
    =\|A\|_2,
    \end{aligned}
    \]

    其中令 $y=P^Tx$。

    对 Frobenius 范数，利用迹的循环不变性，

    \[
    \begin{aligned}
    \|QAP^T\|_F^2
    &=\operatorname{tr}\bigl(PA^TQ^TQAP^T\bigr)\\
    &=\operatorname{tr}\bigl(PA^TAP^T\bigr)\\
    &=\operatorname{tr}(A^TA)
    =\|A\|_F^2.
    \end{aligned}
    \]

    因而结论成立。$\square$

---

## 2. 矩阵扰动与逆矩阵

### 2.1 Neumann 级数

!!! success "引理 2.1（Neumann 级数）"

    设 $F\in\mathbb{R}^{n\times n}$，并且对某个相容矩阵范数 $\|\cdot\|_p$ 有

    \[
    \|F\|_p<1.
    \]

    则 $I-F$ 非奇异，并且

    \[
    (I-F)^{-1}=\sum_{k=0}^{\infty}F^k.
    \]

    此外，

    \[
    \|(I-F)^{-1}\|_p
    \leq\frac{1}{1-\|F\|_p}.
    \]

??? proof "引理 2.1 的证明（点击展开）"

    令

    \[
    S_N=\sum_{k=0}^NF^k.
    \]

    当 $M>N$ 时，由范数的三角不等式和相容性，

    \[
    \begin{aligned}
    \|S_M-S_N\|_p
    &\leq\sum_{k=N+1}^M\|F^k\|_p\\
    &\leq\sum_{k=N+1}^M\|F\|_p^k\\
    &\leq\frac{\|F\|_p^{N+1}}{1-\|F\|_p}.
    \end{aligned}
    \]

    因为 $\|F\|_p<1$，右端趋于 $0$，所以 $\{S_N\}$ 是 Cauchy 列。有限维赋范空间是完备的，故存在矩阵 $S$ 使得 $S_N\to S$。

    又因为

    \[
    S_N(I-F)=(I-F)S_N=I-F^{N+1},
    \]

    且 $F^{N+1}\to0$，令 $N\to\infty$ 得

    \[
    S(I-F)=(I-F)S=I.
    \]

    因此 $S=(I-F)^{-1}$。最后，

    \[
    \|(I-F)^{-1}\|_p
    \leq\sum_{k=0}^{\infty}\|F\|_p^k
    =\frac{1}{1-\|F\|_p}.
    \]

    证明完毕。$\square$

### 2.2 逆矩阵的扰动界

!!! success "定理 2.2（逆矩阵的扰动）"

    设 $A\in\mathbb{R}^{n\times n}$ 非奇异，$E\in\mathbb{R}^{n\times n}$，并令

    \[
    r=\|A^{-1}E\|_p<1.
    \]

    则 $A+E$ 非奇异，并且

    \[
    \|(A+E)^{-1}-A^{-1}\|_p
    \leq
    \frac{\|E\|_p\|A^{-1}\|_p^2}{1-r}.
    \]

??? proof "定理 2.2 的证明（点击展开）"

    令

    \[
    F=-A^{-1}E.
    \]

    则 $\|F\|_p=r<1$，并且

    \[
    A+E=A(I-F).
    \]

    由 Neumann 级数引理，$I-F$ 非奇异，因而 $A+E$ 非奇异，且

    \[
    (A+E)^{-1}=(I-F)^{-1}A^{-1}.
    \]

    利用逆矩阵之差恒等式，

    \[
    (A+E)^{-1}-A^{-1}
    =-A^{-1}E(A+E)^{-1}.
    \]

    又有

    \[
    \|(A+E)^{-1}\|_p
    \leq\|(I-F)^{-1}\|_p\|A^{-1}\|_p
    \leq\frac{\|A^{-1}\|_p}{1-r}.
    \]

    因此

    \[
    \begin{aligned}
    \|(A+E)^{-1}-A^{-1}\|_p
    &\leq\|A^{-1}\|_p\|E\|_p\|(A+E)^{-1}\|_p\\
    &\leq\frac{\|E\|_p\|A^{-1}\|_p^2}{1-r}.
    \end{aligned}
    \]

    证明完毕。$\square$

---

## 3. 奇异值分解

### 3.1 SVD 的存在性

!!! success "定理 3.1（实矩阵的奇异值分解）"

    对任意矩阵 $A\in\mathbb{R}^{m\times n}$，存在正交矩阵

    \[
    U\in\mathbb{R}^{m\times m},
    \qquad
    V\in\mathbb{R}^{n\times n},
    \]

    以及矩形对角矩阵 $\Sigma\in\mathbb{R}^{m\times n}$，使得

    \[
    A=U\Sigma V^T.
    \]

    其中 $\Sigma$ 的对角元满足

    \[
    \sigma_1\geq\sigma_2\geq\cdots\geq\sigma_p\geq0,
    \qquad
    p=\min\{m,n\}.
    \]

    $\sigma_i$ 称为 $A$ 的奇异值；$U$ 的列向量 $u_i$ 称为左奇异向量；$V$ 的列向量 $v_i$ 称为右奇异向量。

??? proof "定理 3.1 的证明（点击展开）"

    下面给出课堂笔记中的归纳证明。只需证明 $m\leq n$ 的情形；当 $m>n$ 时，可以对 $A^T$ 应用结论后再转置。

    当 $A=0$ 时结论显然成立。以下设 $A\neq0$，并对 $m$ 作归纳。

    当 $m=1$ 时，取

    \[
    v_1=\frac{A^T}{\|A\|_2},
    \]

    并将 $v_1$ 扩充为 $\mathbb{R}^n$ 的一组标准正交基。以这些向量为列组成正交矩阵 $V$，再取 $U=[1]$，便有

    \[
    U^TAV=(\|A\|_2,0,\ldots,0).
    \]

    现在设结论对行数小于 $m$ 的矩阵成立。由定理 1.5，存在单位向量 $x\in\mathbb{R}^n$，使得

    \[
    A^TAx=\|A\|_2^2x.
    \]

    令

    \[
    \sigma_1=\|A\|_2,
    \qquad
    y=\frac{Ax}{\sigma_1}.
    \]

    则 $\|y\|_2=1$ 且 $Ax=\sigma_1y$。分别把 $y$ 和 $x$ 扩充为 $\mathbb{R}^m$ 与 $\mathbb{R}^n$ 的标准正交基，并令

    \[
    U_0=(y,u_2,\ldots,u_m),
    \qquad
    V_0=(x,v_2,\ldots,v_n).
    \]

    于是

    \[
    U_0^TAV_0=
    \begin{pmatrix}
    \sigma_1&w^T\\
    0&B
    \end{pmatrix}.
    \]

    由于正交变换保持谱范数，左端的谱范数等于 $\sigma_1$。若 $w\neq0$，令

    \[
    s=\begin{pmatrix}\sigma_1\\w\end{pmatrix},
    \]

    则

    \[
    \left\|
    \begin{pmatrix}
    \sigma_1&w^T\\
    0&B
    \end{pmatrix}s
    \right\|_2
    \geq\sigma_1^2+\|w\|_2^2
    >\sigma_1\|s\|_2,
    \]

    这与该矩阵的谱范数等于 $\sigma_1$ 矛盾，因此 $w=0$。所以

    \[
    U_0^TAV_0=
    \begin{pmatrix}
    \sigma_1&0\\
    0&B
    \end{pmatrix}.
    \]

    对 $B\in\mathbb{R}^{(m-1)\times(n-1)}$ 应用归纳假设，再把所得正交矩阵嵌入分块对角矩阵，即可得到 $A$ 的奇异值分解。$\square$

!!! note "SVD 的非唯一性"

    奇异值本身是唯一确定的，但奇异向量不一定唯一。例如，若某个奇异值具有重数大于 $1$，则可以在相应的奇异子空间中作任意正交变换；即使奇异值互异，也可以同时改变一对左右奇异向量的符号。

### 3.2 左、右奇异向量的关系

由

\[
AV=U\Sigma,
\qquad
A^TU=V\Sigma^T,
\]

立即得到以下关系。

!!! success "推论 3.2（奇异向量方程）"

    对 $i=1,\ldots,p$，有

    \[
    Av_i=\sigma_i u_i,
    \qquad
    A^Tu_i=\sigma_i v_i.
    \]

进而，

\[
A^TAv_i=\sigma_i^2v_i,
\qquad
AA^Tu_i=\sigma_i^2u_i.
\]

因此，$A$ 的非零奇异值的平方恰好是 $A^TA$（也等价地是 $AA^T$）的非零特征值。

### 3.3 范数、秩与基本子空间

!!! success "推论 3.3（范数与奇异值）"

    对任意 $A\in\mathbb{R}^{m\times n}$，

    \[
    \|A\|_2=\sigma_1,
    \qquad
    \|A\|_F=\left(\sum_{i=1}^p\sigma_i^2\right)^{1/2}.
    \]

设 $r=\operatorname{rank}(A)$，则

\[
\sigma_1\geq\cdots\geq\sigma_r>0,
\qquad
\sigma_{r+1}=\cdots=\sigma_p=0,
\]

并且 SVD 可以写成秩一矩阵之和：

\[
A=\sum_{i=1}^r\sigma_i u_iv_i^T.
\]

!!! success "推论 3.4（四个基本子空间）"

    若 $\operatorname{rank}(A)=r$，则

    \[
    \operatorname{Range}(A)
    =\operatorname{span}\{u_1,\ldots,u_r\},
    \]

    \[
    \operatorname{Null}(A^T)
    =\operatorname{span}\{u_{r+1},\ldots,u_m\},
    \]

    \[
    \operatorname{Range}(A^T)
    =\operatorname{span}\{v_1,\ldots,v_r\},
    \]

    \[
    \operatorname{Null}(A)
    =\operatorname{span}\{v_{r+1},\ldots,v_n\}.
    \]

### 3.4 最小奇异值的刻画

!!! success "引理 3.5（最小奇异值）"

    设 $A\in\mathbb{R}^{m\times n}$，其中 $m\geq n$，且 $A$ 列满秩。则对任意 $x\in\mathbb{R}^n$，

    \[
    \|Ax\|_2\geq\sigma_n(A)\|x\|_2.
    \]

    此外，存在非零向量 $x_0$ 使等号成立。因此

    \[
    \sigma_n(A)
    =\min_{\|x\|_2=1}\|Ax\|_2.
    \]

??? proof "引理 3.5 的证明（点击展开）"

    设 $A=U\Sigma V^T$，并令 $y=V^Tx$。由于 $V$ 正交，$\|y\|_2=\|x\|_2$。于是

    \[
    \|Ax\|_2^2
    =\|\Sigma y\|_2^2
    =\sum_{i=1}^n\sigma_i^2y_i^2
    \geq\sigma_n^2\sum_{i=1}^ny_i^2
    =\sigma_n^2\|x\|_2^2.
    \]

    取 $x=v_n$ 时等号成立。$\square$

---

## 4. 奇异值的扰动性质

### 4.1 最大与最小奇异值的扰动界

!!! success "定理 4.1（奇异值扰动的基本界）"

    设 $A,E\in\mathbb{R}^{m\times n}$，则

    \[
    \sigma_{\max}(A+E)
    \leq\sigma_{\max}(A)+\|E\|_2.
    \]

    若进一步有 $m\geq n$，并且 $A+E$ 列满秩，则

    \[
    \sigma_{\min}(A+E)
    \geq\sigma_{\min}(A)-\|E\|_2.
    \]

??? proof "定理 4.1 的证明（点击展开）"

    第一式由谱范数的三角不等式直接得到：

    \[
    \sigma_{\max}(A+E)
    =\|A+E\|_2
    \leq\|A\|_2+\|E\|_2.
    \]

    对第二式，取单位向量 $x_0$，使得

    \[
    \|(A+E)x_0\|_2=\sigma_{\min}(A+E).
    \]

    于是

    \[
    \begin{aligned}
    \sigma_{\min}(A)
    &\leq\|Ax_0\|_2\\
    &\leq\|(A+E)x_0\|_2+\|Ex_0\|_2\\
    &\leq\sigma_{\min}(A+E)+\|E\|_2.
    \end{aligned}
    \]

    移项即得结论。$\square$

!!! example "推论 4.2（添加一列后的奇异值）"

    设 $A\in\mathbb{R}^{m\times n}$、$m>n$，并且 $z\in\mathbb{R}^m$。记 $[A\ z]$ 为在 $A$ 右侧添加列向量 $z$ 后得到的矩阵，则

    \[
    \sigma_{\max}([A\ z])\geq\sigma_{\max}(A),
    \]

    \[
    \sigma_{\min}([A\ z])\leq\sigma_{\min}(A).
    \]

这是因为 $[A\ z](x,0)^T=Ax$，而最大奇异值是在单位球面上取最大值，最小奇异值是在单位球面上取最小值。增加一个变量方向会扩大最大化问题的可行集合，同时也会扩大最小化问题的可行集合。

---

## 5. 最佳低秩逼近

设 $A\in\mathbb{R}^{m\times n}$ 的秩为 $r$，其紧奇异值分解为

\[
A=\sum_{i=1}^r\sigma_i u_iv_i^T,
\qquad
\sigma_1\geq\cdots\geq\sigma_r>0.
\]

对 $k<r$，定义截断 SVD

\[
A_k=\sum_{i=1}^k\sigma_i u_iv_i^T.
\]

显然 $\operatorname{rank}(A_k)=k$。

### 5.1 谱范数下的 Eckart-Young 定理

!!! success "定理 5.1（Eckart-Young 定理）"

    若 $k<r=\operatorname{rank}(A)$，则

    \[
    \min_{\operatorname{rank}(B)\leq k}\|A-B\|_2
    =\|A-A_k\|_2
    =\sigma_{k+1}.
    \]

    因为 $A_k$ 的秩恰好为 $k$，将约束写成 $\operatorname{rank}(B)=k$ 时也得到相同的最小值。

??? proof "定理 5.1 的证明（点击展开）"

    首先，由正交不变性，

    \[
    \begin{aligned}
    \|A-A_k\|_2
    &=\left\|U^T(A-A_k)V\right\|_2\\
    &=\left\|\operatorname{diag}
    (0,\ldots,0,\sigma_{k+1},\ldots,\sigma_r)\right\|_2\\
    &=\sigma_{k+1}.
    \end{aligned}
    \]

    下面证明任何秩不超过 $k$ 的矩阵都不能得到更小的误差。设 $\operatorname{rank}(B)\leq k$，则

    \[
    \dim\operatorname{Null}(B)\geq n-k.
    \]

    令

    \[
    S=\operatorname{span}\{v_1,\ldots,v_{k+1}\}.
    \]

    因为 $\dim S=k+1$，由维数公式，

    \[
    \dim\bigl(\operatorname{Null}(B)\cap S\bigr)
    \geq(n-k)+(k+1)-n=1.
    \]

    因此可以取单位向量 $z\in\operatorname{Null}(B)\cap S$。将其写为

    \[
    z=\sum_{i=1}^{k+1}(v_i^Tz)v_i.
    \]

    由于 $Bz=0$，

    \[
    \begin{aligned}
    \|A-B\|_2^2
    &\geq\|(A-B)z\|_2^2\\
    &=\|Az\|_2^2\\
    &=\sum_{i=1}^{k+1}\sigma_i^2(v_i^Tz)^2\\
    &\geq\sigma_{k+1}^2\sum_{i=1}^{k+1}(v_i^Tz)^2
    =\sigma_{k+1}^2.
    \end{aligned}
    \]

    所以 $\|A-B\|_2\geq\sigma_{k+1}$。结合 $A_k$ 达到该下界，定理得证。$\square$

### 5.2 满秩矩阵集合的开性与稠密性

Eckart-Young 定理说明，矩阵到低秩矩阵集合的谱范数距离由相应的奇异值刻画。

!!! success "命题 5.2（满秩矩阵集合是开且稠密的）"

    令 $p=\min\{m,n\}$。所有秩为 $p$ 的 $m\times n$ 实矩阵组成的集合在 $\mathbb{R}^{m\times n}$ 中既是开集，也是稠密集。

??? proof "命题 5.2 的说明（点击展开）"

    若 $A$ 满秩，则 $\sigma_p(A)>0$。当 $m\geq n$ 时直接应用奇异值扰动界；当 $m<n$ 时对转置矩阵应用同一结论。于是，只要

    \[
    \|E\|_2<\sigma_p(A),
    \]

    就有

    \[
    \sigma_p(A+E)
    \geq\sigma_p(A)-\|E\|_2>0.
    \]

    因而 $A+E$ 仍满秩，所以满秩矩阵集合是开集。

    另一方面，设 $A=U\Sigma V^T$ 不满秩。对任意 $\varepsilon>0$，把 $\Sigma$ 中所有为零的对角元替换为 $\varepsilon$，得到 $\Sigma_\varepsilon$，并令

    \[
    A_\varepsilon=U\Sigma_\varepsilon V^T.
    \]

    则 $A_\varepsilon$ 满秩，且

    \[
    \|A_\varepsilon-A\|_2=\varepsilon.
    \]

    因此任意矩阵都可以被满秩矩阵任意逼近，即满秩矩阵集合是稠密集。$\square$

### 5.3 奇异值的和不等式

!!! success "引理 5.3"

    设 $X,Y\in\mathbb{R}^{m\times n}$，且指标满足 $i+j-1\leq p=\min\{m,n\}$，则

    \[
    \sigma_i(X)+\sigma_j(Y)
    \geq\sigma_{i+j-1}(X+Y).
    \]

??? proof "引理 5.3 的证明（点击展开）"

    分别令 $X_{i-1}$ 和 $Y_{j-1}$ 为 $X$ 与 $Y$ 的截断 SVD。由 Eckart-Young 定理，

    \[
    \sigma_i(X)=\|X-X_{i-1}\|_2,
    \qquad
    \sigma_j(Y)=\|Y-Y_{j-1}\|_2.
    \]

    因而

    \[
    \begin{aligned}
    \sigma_i(X)+\sigma_j(Y)
    &\geq\|X+Y-(X_{i-1}+Y_{j-1})\|_2.
    \end{aligned}
    \]

    又因为

    \[
    \operatorname{rank}(X_{i-1}+Y_{j-1})
    \leq(i-1)+(j-1)=i+j-2,
    \]

    再次应用 Eckart-Young 定理，得到

    \[
    \|X+Y-(X_{i-1}+Y_{j-1})\|_2
    \geq\sigma_{i+j-1}(X+Y).
    \]

    合并两式即得结论。$\square$

### 5.4 Frobenius 范数下的最佳低秩逼近

!!! success "定理 5.4（Frobenius 范数下的低秩逼近）"

    若 $k<r=\operatorname{rank}(A)$，则

    \[
    \min_{\operatorname{rank}(B)\leq k}\|A-B\|_F
    =\|A-A_k\|_F
    =\left(\sum_{i=k+1}^r\sigma_i^2\right)^{1/2}.
    \]

??? proof "定理 5.4 的证明（点击展开）"

    首先，由 Frobenius 范数的正交不变性，

    \[
    \|A-A_k\|_F^2
    =\sum_{i=k+1}^r\sigma_i^2.
    \]

    下面任取满足 $\operatorname{rank}(B)\leq k$ 的矩阵 $B$。因为

    \[
    \sigma_{k+1}(B)=0,
    \]

    在引理 5.3 中取 $X=A-B$、$Y=B$ 和 $j=k+1$，可得

    \[
    \sigma_i(A-B)
    \geq\sigma_{i+k}(A),
    \qquad i=1,\ldots,r-k.
    \]

    因此

    \[
    \begin{aligned}
    \|A-B\|_F^2
    &=\sum_{i=1}^p\sigma_i^2(A-B)\\
    &\geq\sum_{i=1}^{r-k}\sigma_i^2(A-B)\\
    &\geq\sum_{i=1}^{r-k}\sigma_{i+k}^2(A)
    =\sum_{i=k+1}^r\sigma_i^2(A).
    \end{aligned}
    \]

    截断 SVD $A_k$ 达到这个下界，所以它是 Frobenius 范数下的最佳秩 $k$ 逼近。$\square$

---

## 6. 本章小结

本章的核心结论如下：

1. 在有限维矩阵空间中，所有范数等价；常用矩阵范数之间存在显式比较不等式。
2. 谱范数和 Frobenius 范数在左右正交变换下保持不变。
3. 当 $\|F\|<1$ 时，$(I-F)^{-1}$ 可以表示为 Neumann 级数；这一结论给出了逆矩阵的扰动界。
4. 任意实矩阵都存在奇异值分解 $A=U\Sigma V^T$，奇异值同时刻画矩阵的范数、秩和四个基本子空间。
5. 截断 SVD 在谱范数和 Frobenius 范数下都是最佳低秩逼近。
