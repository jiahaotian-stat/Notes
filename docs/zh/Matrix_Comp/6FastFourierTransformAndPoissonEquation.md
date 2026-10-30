# 第六章：快速傅里叶变换与 Poisson 方程

本章从离散傅里叶变换出发，推导基 2 快速傅里叶变换及蝶形运算，并讨论离散正弦变换、离散余弦变换、循环矩阵的对角化，以及 Poisson 方程有限差分离散后的快速求解。

---

## 1. 离散傅里叶变换

### 1.1 从连续变换到离散变换

连续 Fourier 变换的一种规范为

\[
\widehat{f}(\xi)=\frac{1}{2\pi}\int_{-\infty}^{\infty}f(x)e^{-\mathrm{i}\xi x}\,\mathrm{d}x.
\]

离散情形中，令

\[
\omega_n=e^{-2\pi\mathrm{i}/n}.
\]

!!! info "定义 1.1（DFT 矩阵）"

    $n$ 阶离散 Fourier 变换矩阵定义为

    \[
    F_n=
    \left[
    \omega_n^{jk}
    \right]_{j,k=0}^{n-1}.
    \]

    对 $x=(x_0,\ldots,x_{n-1})^{\mathsf T}\in\mathbb{C}^n$，其离散 Fourier 变换为

    \[
    y=F_nx,
    \]

    即

    \[
    y_k=\sum_{j=0}^{n-1}\omega_n^{kj}x_j,
    \qquad k=0,\ldots,n-1.
    \]

直接矩阵向量乘法需要 $O(n^2)$ 次运算，而 FFT 利用 $F_n$ 的递归结构将复杂度降为 $O(n\log n)$。

### 1.2 正交性与逆变换

若 $p\neq q$，则

\[
\sum_{k=0}^{n-1}\omega_n^{k(q-p)}
=\frac{1-\omega_n^{n(q-p)}}{1-\omega_n^{q-p}}
=0.
\]

当 $p=q$ 时，上式等于 $n$。因此

\[
F_n^{\mathsf H}F_n=nI.
\]

!!! success "定理 1.2（DFT 矩阵的逆）"

    \[
    F_n^{-1}=\frac{1}{n}F_n^{\mathsf H}
    =\frac{1}{n}\overline{F_n}.
    \]

    因而归一化矩阵

    \[
    \frac{1}{\sqrt n}F_n
    \]

    是酉矩阵。

逆变换也可以用相同的 FFT 结构在 $O(n\log n)$ 时间内完成。

---

## 2. 基 2 快速傅里叶变换

### 2.1 偶数项与奇数项分裂

设

\[
n=2m.
\]

将输入按偶数下标和奇数下标拆分为

\[
x^{(e)}=(x_0,x_2,\ldots,x_{2m-2})^{\mathsf T},
\]

\[
x^{(o)}=(x_1,x_3,\ldots,x_{2m-1})^{\mathsf T}.
\]

因为

\[
\omega_{2m}^2=\omega_m,
\]

所以对 $k=0,\ldots,m-1$，

\[
y_k
=\sum_{r=0}^{m-1}\omega_m^{kr}x_{2r}
+\omega_{2m}^k\sum_{r=0}^{m-1}\omega_m^{kr}x_{2r+1}.
\]

记

\[
y^{(e)}=F_mx^{(e)},
\qquad
y^{(o)}=F_mx^{(o)},
\]

则

\[
y_k=y_k^{(e)}+\omega_{2m}^ky_k^{(o)},
\]

\[
y_{k+m}=y_k^{(e)}-\omega_{2m}^ky_k^{(o)}.
\]

这两个式子构成一次蝶形运算。

### 2.2 矩阵分解

令 $P_{2m}$ 为把偶数下标分量排在前、奇数下标分量排在后的置换矩阵，并令

\[
D_m=\operatorname{diag}
\left(
1,\omega_{2m},\ldots,\omega_{2m}^{m-1}
\right).
\]

!!! success "定理 2.1（Cooley-Tukey 分解）"

    \[
    F_{2m}P_{2m}
    =
    \begin{bmatrix}
    I_m&D_m\\
    I_m&-D_m
    \end{bmatrix}
    \begin{bmatrix}
    F_m&0\\
    0&F_m
    \end{bmatrix}.
    \]

该分解把一个 $2m$ 点 DFT 化为两个 $m$ 点 DFT 加 $O(m)$ 次蝶形运算。

### 2.3 复杂度

设 $T(n)$ 为 $n$ 点 FFT 的运算量。对 $n=2^t$，

\[
T(n)=2T(n/2)+O(n).
\]

因此

\[
T(n)=O(n\log_2 n).
\]

!!! note "蝶形运算"

    已知 $a=y_k^{(e)}$、$b=y_k^{(o)}$ 和旋转因子 $w=\omega_{2m}^k$，一次蝶形计算输出

    \[
    a+wb
    \]

    与

    \[
    a-wb.
    \]

---

## 3. 离散正弦变换与离散余弦变换

### 3.1 离散正弦变换

!!! info "定义 3.1（DST-I）"

    对 $x=(x_1,\ldots,x_{m-1})^{\mathsf T}$，定义

    \[
    y_k=\sum_{j=1}^{m-1}\sin\left(\frac{kj\pi}{m}\right)x_j,
    \qquad k=1,\ldots,m-1.
    \]

    记 $y=S_{m-1}x$。

例如，当 $m=3$ 时，

\[
y_1=\sin\left(\frac{\pi}{3}\right)x_1
+\sin\left(\frac{2\pi}{3}\right)x_2,
\]

\[
y_2=\sin\left(\frac{2\pi}{3}\right)x_1
+\sin\left(\frac{4\pi}{3}\right)x_2.
\]

DST-I 可以通过对 $x$ 作奇延拓后调用 FFT 来计算，因此复杂度为 $O(m\log m)$。

### 3.2 离散余弦变换

!!! info "定义 3.2（DCT-I）"

    对 $x=(x_0,\ldots,x_m)^{\mathsf T}$，定义

    \[
    y_k=\frac{x_0}{2}
    +\sum_{j=1}^{m-1}\cos\left(\frac{kj\pi}{m}\right)x_j
    +\frac{(-1)^kx_m}{2},
    \qquad k=0,\ldots,m.
    \]

DCT-I 可以通过对 $x$ 作偶延拓后调用 FFT 来计算。

### 3.3 与 DFT 的关系

对 DST，构造长度 $2m$ 的奇延拓

\[
\widetilde{x}
=
(0,x_1,\ldots,x_{m-1},0,-x_{m-1},\ldots,-x_1)^{\mathsf T}.
\]

其 DFT 的虚部给出正弦变换。对 DCT，构造偶延拓

\[
\widehat{x}
=
(x_0,x_1,\ldots,x_{m-1},x_m,x_{m-1},\ldots,x_1)^{\mathsf T},
\]

其 DFT 的实部给出余弦变换。

!!! info "边界条件的对应"

    奇延拓自然对应齐次 Dirichlet 边界条件，因此 DST 常用于固定端点问题。

    偶延拓自然对应齐次 Neumann 边界条件，因此 DCT 常用于零法向导数问题。

---

## 4. 循环矩阵

### 4.1 循环移位矩阵

令

\[
D_n=
\begin{bmatrix}
0&0&\cdots&0&1\\
1&0&\cdots&0&0\\
0&1&\ddots&0&0\\
\vdots&\ddots&\ddots&\ddots&\vdots\\
0&0&\cdots&1&0
\end{bmatrix}.
\]

则 $D_n$ 表示循环移位，并满足

\[
D_n^n=I.
\]

!!! info "定义 4.1（循环矩阵）"

    循环矩阵可写成

    \[
    C(z)=\sum_{k=0}^{n-1}z_kD_n^k.
    \]

    它的每一列都是前一列的循环移位。

循环矩阵向量乘法对应离散循环卷积。

### 4.2 Fourier 对角化

!!! success "定理 4.2（循环矩阵的特征分解）"

    DFT 矩阵同时对角化所有 $n$ 阶循环矩阵。具体地，

    \[
    F_nD_nF_n^{-1}
    =\operatorname{diag}
    \left(
    1,\omega_n,\ldots,\omega_n^{n-1}
    \right).
    \]

    因而

    \[
    F_nC(z)F_n^{-1}
    =\operatorname{diag}(F_nz).
    \]

??? proof "证明（点击展开）"

    $F_n$ 的第 $j$ 个 Fourier 向量是 $D_n$ 的特征向量，对应特征值 $\omega_n^j$。因此 $F_n$ 对角化 $D_n$。由于 $C(z)$ 是 $D_n$ 的多项式，同一组特征向量也对角化 $C(z)$，其特征值为该多项式在各个 $n$ 次单位根处的取值，即 $F_nz$。 $\square$

若循环矩阵 $C(z)$ 非奇异，则方程

\[
C(z)x=b
\]

可通过两次 FFT 与一次逐分量相除在 $O(n\log n)$ 时间内求解。

---

## 5. 一维 Poisson 方程

考虑边值问题

\[
-u''(x)=f(x),
\qquad \alpha<x<\beta,
\]

以及 Dirichlet 边界条件

\[
u(\alpha)=u_\alpha,
\qquad
u(\beta)=u_\beta.
\]

取

\[
h=\frac{\beta-\alpha}{m},
\qquad
x_i=\alpha+ih,
\qquad i=0,\ldots,m.
\]

在内部网格点使用中心差分，得到

\[
\frac{-u_{i-1}+2u_i-u_{i+1}}{h^2}=f_i,
\qquad i=1,\ldots,m-1.
\]

令

\[
T_{m-1}^{(DD)}=
\begin{bmatrix}
2&-1&&&\\
-1&2&-1&&\\
&\ddots&\ddots&\ddots&\\
&&-1&2&-1\\
&&&-1&2
\end{bmatrix}.
\]

则离散系统为

\[
T_{m-1}^{(DD)}u=b,
\]

其中

\[
u=(u_1,\ldots,u_{m-1})^{\mathsf T},
\]

\[
b=h^2(f_1,\ldots,f_{m-1})^{\mathsf T}
+(u_\alpha,0,\ldots,0,u_\beta)^{\mathsf T}.
\]

!!! success "定理 5.1（Dirichlet 差分矩阵的正弦对角化）"

    $T_{m-1}^{(DD)}$ 的特征向量由 DST-I 给出，其特征值为

    \[
    \lambda_k=2-2\cos\left(\frac{k\pi}{m}\right)
    =4\sin^2\left(\frac{k\pi}{2m}\right),
    \qquad k=1,\ldots,m-1.
    \]

因此可先对 $b$ 作 DST，逐分量除以 $\lambda_k$，再作逆 DST，从而以 $O(m\log m)$ 复杂度求解。

---

## 6. 二维 Poisson 方程与 Kronecker 积

### 6.1 五点差分格式

设

\[
\Omega={(x,y):\alpha_x\leq x\leq\beta_x,
\ \alpha_y\leq y\leq\beta_y},
\]

并考虑

\[
-\Delta u=f,
\qquad
\Delta=\frac{\partial^2}{\partial x^2}
+\frac{\partial^2}{\partial y^2}.
\]

令

\[
h_x=\frac{\beta_x-\alpha_x}{m_x+1},
\qquad
h_y=\frac{\beta_y-\alpha_y}{m_y+1}.
\]

在内部网格点，五点差分格式为

\[
\frac{-u_{i-1,j}+2u_{i,j}-u_{i+1,j}}{h_x^2}
+\frac{-u_{i,j-1}+2u_{i,j}-u_{i,j+1}}{h_y^2}
=f_{i,j}.
\]

### 6.2 Kronecker 积表示

!!! info "定义 6.1（Kronecker 积）"

    若 $A=(a_{ij})\in\mathbb{R}^{m\times n}$，则

    \[
    A\otimes B=
    \begin{bmatrix}
    a_{11}B&\cdots&a_{1n}B\\
    \vdots&\ddots&\vdots\\
    a_{m1}B&\cdots&a_{mn}B
    \end{bmatrix}.
    \]

按列堆叠所有内部未知量后，二维离散系统可写成

\[
Au=b,
\]

其中

\[
A=\frac{1}{h_x^2}
\left(I_{m_y}\otimes T_{m_x}^{(DD)}\right)
+\frac{1}{h_y^2}
\left(T_{m_y}^{(DD)}\otimes I_{m_x}\right).
\]

这是两个一维离散算子的 Kronecker 和。

### 6.3 可分离对角化

若

\[
T_{m_x}^{(DD)}=S_x^{-1}\Lambda_xS_x,
\]

\[
T_{m_y}^{(DD)}=S_y^{-1}\Lambda_yS_y,
\]

则

\[
(S_y\otimes S_x)A(S_y^{-1}\otimes S_x^{-1})
=\frac{1}{h_x^2}(I\otimes\Lambda_x)
+\frac{1}{h_y^2}(\Lambda_y\otimes I),
\]

右端为对角矩阵。因此二维 Poisson 方程可以通过沿两个坐标方向分别执行 DST、逐点相除以及逆 DST 来快速求解。

---

## 7. 本章总结

!!! success "核心结论"

    1. DFT 矩阵满足 $F_n^{\mathsf H}F_n=nI$，其逆由共轭转置直接给出。

    2. 基 2 FFT 将一个 $n$ 点 DFT 递归分解为两个 $n/2$ 点 DFT，并把复杂度从 $O(n^2)$ 降为 $O(n\log n)$。

    3. DST 与 DCT 分别对应奇延拓与偶延拓，可借助 FFT 快速计算。

    4. DFT 同时对角化所有循环矩阵，循环卷积因此可以用 FFT 快速完成。

    5. 一维 Dirichlet Poisson 差分矩阵由 DST 对角化。

    6. 二维 Poisson 离散矩阵是两个一维算子的 Kronecker 和，可利用可分离变换快速求解。


