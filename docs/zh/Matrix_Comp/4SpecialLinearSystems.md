# 第四章：特殊线性系统与结构化分解

本章研究具有特殊结构的线性方程组。与一般稠密系统相比，这些矩阵往往不需要选主元，或者可以利用结构显著减少计算量。我们依次讨论对角占优矩阵、对称矩阵的 \(LDL^{\mathsf T}\) 分解、正定系统、Cholesky 分解、极分解以及带状线性系统。

---

## 1. 对角占优矩阵

### 1.1 行对角占优与列对角占优

设

\[
A=(a_{ij})\in\mathbb{R}^{n\times n}.
\]

!!! info "定义 1.1（对角占优）"

    如果对每个 \(i=1,\ldots,n\)，都有

    \[
    |a_{ii}|\geq\sum_{j\neq i}|a_{ij}|,
    \]

    则称 \(A\) 为行对角占优矩阵。

    如果对每个 \(j=1,\ldots,n\)，都有

    \[
    |a_{jj}|\geq\sum_{i\neq j}|a_{ij}|,
    \]

    则称 \(A\) 为列对角占优矩阵。

    如果相应不等式对每一行或每一列都严格成立，则分别称 \(A\) 严格行对角占优或严格列对角占优。

对角占优矩阵不一定可逆。例如，零矩阵既行对角占优又列对角占优，但它显然不可逆。严格对角占优则能保证可逆性。

### 1.2 严格对角占优矩阵的可逆性

!!! success "定理 1.2（严格对角占优矩阵可逆）"

    严格行对角占优矩阵与严格列对角占优矩阵都是非奇异矩阵。

??? proof "严格行对角占优情形的证明（点击展开）"

    反设 \(A\) 不可逆，则存在 \(x\neq0\) 使得

    \[
    Ax=0.
    \]

    取

    \[
    k\in\operatorname*{arg\,max}_{1\leq i\leq n}|x_i|.
    \]

    由于 \(x\neq0\)，有 \(|x_k|>0\)。由方程的第 \(k\) 行，

    \[
    a_{kk}x_k=-\sum_{j\neq k}a_{kj}x_j.
    \]

    因而

    \[
    \begin{aligned}
    |a_{kk}||x_k|
    &\leq\sum_{j\neq k}|a_{kj}||x_j|\\
    &\leq\left(\sum_{j\neq k}|a_{kj}|\right)|x_k|.
    \end{aligned}
    \]

    消去 \(|x_k|\) 后得到

    \[
    |a_{kk}|\leq\sum_{j\neq k}|a_{kj}|,
    \]

    这与严格行对角占优矛盾。因此 \(A\) 可逆。严格列对角占优情形可对 \(A^{\mathsf T}\) 应用同一结论。 \(\square\)

### 1.3 逆矩阵的范数估计

定义行占优因子

\[
\delta
=\min_{1\leq i\leq n}
\left(
|a_{ii}|-\sum_{j\neq i}|a_{ij}|
\right),
\]

以及列占优因子

\[
\mu
=\min_{1\leq j\leq n}
\left(
|a_{jj}|-\sum_{i\neq j}|a_{ij}|
\right).
\]

!!! success "定理 1.3（逆矩阵范数界）"

    若 \(A\) 严格行对角占优，则

    \[
    \|A^{-1}\|_{\infty}\leq\frac{1}{\delta}.
    \]

    若 \(A\) 严格列对角占优，则

    \[
    \|A^{-1}\|_1\leq\frac{1}{\mu}.
    \]

??? proof "无穷范数估计的证明（点击展开）"

    任取 \(x\in\mathbb{R}^n\)，并取 \(k\) 使得

    \[
    |x_k|=\|x\|_{\infty}.
    \]

    则

    \[
    \begin{aligned}
    \|Ax\|_{\infty}
    &\geq |(Ax)_k|\\
    &\geq |a_{kk}||x_k|-\sum_{j\neq k}|a_{kj}||x_j|\\
    &\geq
    \left(
    |a_{kk}|-\sum_{j\neq k}|a_{kj}|
    \right)\|x\|_{\infty}\\
    &\geq\delta\|x\|_{\infty}.
    \end{aligned}
    \]

    令 \(x=A^{-1}y\)，得到

    \[
    \|A^{-1}y\|_{\infty}
    \leq\frac{1}{\delta}\|y\|_{\infty}.
    \]

    对所有 \(y\neq0\) 取上确界，即得结论。列对角占优情形由

    \[
    \|A^{-1}\|_1
    =\|A^{-\mathsf T}\|_{\infty}
    \]

    以及 \(A^{\mathsf T}\) 的严格行对角占优性得到。 \(\square\)

### 1.4 列对角占优矩阵的 LU 分解

!!! success "定理 1.4（列对角占优与 LU 分解）"

    若 \(A\) 非奇异且列对角占优，则 \(A\) 存在唯一的 LU 分解

    \[
    A=LU,
    \]

    其中 \(L\) 为单位下三角矩阵，并且

    \[
    |l_{ij}|\leq1,
    \qquad i>j.
    \]

??? proof "定理 1.4 的证明思路（点击展开）"

    对矩阵阶数作归纳。将 \(A\) 写成

    \[
    A=
    \begin{bmatrix}
    \alpha & w^{\mathsf T}\\
    v & C
    \end{bmatrix}.
    \]

    非奇异性与列对角占优保证 \(\alpha\neq0\)，且第一步乘子满足

    \[
    \left|\frac{v_i}{\alpha}\right|\leq1.
    \]

    第一步消元后的 Schur 补为

    \[
    B=C-\frac{vw^{\mathsf T}}{\alpha}.
    \]

    可以逐列验证 \(B\) 仍然非奇异且列对角占优。由归纳假设，\(B\) 具有满足同样乘子界的 LU 分解。将第一步乘子与该分解组合，即得 \(A=LU\) 且所有 \(|l_{ij}|\leq1\)。唯一性来自单位下三角规范。 \(\square\)

!!! info "数值意义"

    这一结论说明，对非奇异列对角占优矩阵进行无主元高斯消去时，消元乘子天然得到控制，因此通常不需要额外的选主元操作。

---

## 2. 对称矩阵与 \(LDL^{\mathsf T}\) 分解

### 2.1 从 \(LDM^{\mathsf T}\) 到 \(LDL^{\mathsf T}\)

设非奇异矩阵 \(A\) 已有分解

\[
A=LDM^{\mathsf T},
\]

其中 \(L,M\) 为单位下三角矩阵，\(D\) 为非奇异对角矩阵。

!!! success "定理 2.1（对称矩阵的 \(LDL^{\mathsf T}\) 分解）"

    若 \(A\) 对称，则

    \[
    M=L,
    \]

    从而

    \[
    A=LDL^{\mathsf T}.
    \]

??? proof "定理 2.1 的证明（点击展开）"

    由 \(A=A^{\mathsf T}\)，有

    \[
    LDM^{\mathsf T}=MDL^{\mathsf T}.
    \]

    左乘 \(M^{-1}\)，右乘 \(L^{-\mathsf T}\)，得到

    \[
    M^{-1}LD=D(M^{-1}L)^{\mathsf T}.
    \]

    左边是下三角矩阵，右边是上三角矩阵，因此二者必须都是对角矩阵。又 \(M^{-1}L\) 是单位下三角矩阵，所以

    \[
    M^{-1}L=I.
    \]

    因而 \(M=L\)。 \(\square\)

### 2.2 计算与存储优势

对称矩阵只需存储一半元素。在高斯消去过程中，上三角部分可以由下三角部分通过对称性恢复，因此 \(LDL^{\mathsf T}\) 分解所需的计算量和存储量约为一般 LU 分解的一半。

若需要选主元，为保持对称性必须同时交换相应的行和列，得到

\[
PAP^{\mathsf T}=LDL^{\mathsf T}.
\]

!!! warning "对称性与选主元"

    只交换行会破坏矩阵的对称结构。对称不定系统通常采用对称选主元策略，并可能需要 \(1\times1\) 或 \(2\times2\) 的对角块。

---

## 3. 正定系统

### 3.1 正定与半正定

!!! info "定义 3.1（正定性）"

    若对任意 \(x\neq0\)，都有

    \[
    x^{\mathsf T}Ax>0,
    \]

    则称 \(A\) 为正定矩阵。

    若对任意 \(x\)，都有

    \[
    x^{\mathsf T}Ax\geq0,
    \]

    则称 \(A\) 为半正定矩阵。

这里不预先要求 \(A\) 对称。二次型只由 \(A\) 的对称部分决定，这一点将在后面说明。

### 3.2 主子矩阵保持正定性

!!! success "定理 3.2（正定矩阵的主子矩阵）"

    正定矩阵的每个主子矩阵仍然正定。

??? proof "定理 3.2 的证明（点击展开）"

    取指标

    \[
    1\leq i_1<\cdots<i_k\leq n,
    \]

    并令

    \[
    X=
    \begin{bmatrix}
    e_{i_1} & \cdots & e_{i_k}
    \end{bmatrix}.
    \]

    相应主子矩阵为

    \[
    B=X^{\mathsf T}AX.
    \]

    对任意 \(z\neq0\)，由于 \(X\) 满列秩，所以 \(Xz\neq0\)，从而

    \[
    z^{\mathsf T}Bz
    =(Xz)^{\mathsf T}A(Xz)>0.
    \]

    因此 \(B\) 正定。 \(\square\)

一个直接推论是：正定矩阵的所有顺序主子式均为正。

### 3.3 正定矩阵的 LU 分解

!!! success "定理 3.3（正定矩阵的 LU 分解）"

    若 \(A\in\mathbb{R}^{n\times n}\) 正定，则 \(A\) 存在唯一的 LU 分解

    \[
    A=LU,
    \]

    并且 \(U\) 的所有对角元均为正数。

??? proof "定理 3.3 的证明（点击展开）"

    由定理 3.2，\(A\) 的所有顺序主子矩阵都正定，因此所有顺序主子式均为正。由 LU 分解的存在唯一性定理，\(A\) 具有唯一的 LU 分解。

    又

    \[
    UL^{-\mathsf T}
    =L^{-1}AL^{-\mathsf T}.
    \]

    右端通过合同变换保持正定性，所以 \(UL^{-\mathsf T}\) 正定。由于 \(U\) 与 \(L^{-\mathsf T}\) 均为上三角矩阵，

    \[
    (UL^{-\mathsf T})_{ii}
    =u_{ii}>0.
    \]

    \(\square\)

### 3.4 对称部分与反对称部分

任意实矩阵都可以唯一写成

\[
A=T+S,
\]

其中

\[
T=\frac{A+A^{\mathsf T}}{2},
\qquad
S=\frac{A-A^{\mathsf T}}{2}.
\]

矩阵 \(T\) 对称，矩阵 \(S\) 反对称。

!!! success "命题 3.4（二次型只取决于对称部分）"

    对任意 \(x\in\mathbb{R}^n\)，都有

    \[
    x^{\mathsf T}Sx=0,
    \qquad
    x^{\mathsf T}Ax=x^{\mathsf T}Tx.
    \]

    因此 \(A\) 正定当且仅当其对称部分 \(T\) 正定。

??? proof "命题 3.4 的证明（点击展开）"

    因为 \(S^{\mathsf T}=-S\)，标量 \(x^{\mathsf T}Sx\) 满足

    \[
    x^{\mathsf T}Sx
    =(x^{\mathsf T}Sx)^{\mathsf T}
    =x^{\mathsf T}S^{\mathsf T}x
    =-x^{\mathsf T}Sx.
    \]

    所以 \(x^{\mathsf T}Sx=0\)，进而

    \[
    x^{\mathsf T}Ax
    =x^{\mathsf T}(T+S)x
    =x^{\mathsf T}Tx.
    \]

    \(\square\)

!!! note "稳定性提示"

    对称正定矩阵的反对称部分为零，因此其无主元分解通常具有很好的数值稳定性。对于一般非对称正定矩阵，稳定性还取决于反对称部分相对于对称部分的大小。

---

## 4. Cholesky 分解

### 4.1 分解的存在性与唯一性

!!! success "定理 4.1（Cholesky 分解）"

    若 \(A\in\mathbb{R}^{n\times n}\) 为对称正定矩阵，则存在唯一的下三角矩阵 \(G\)，使得

    \[
    A=GG^{\mathsf T},
    \]

    且 \(G\) 的所有对角元均为正。

??? proof "定理 4.1 的证明（点击展开）"

    对称正定矩阵具有唯一分解

    \[
    A=LDL^{\mathsf T},
    \]

    其中 \(L\) 为单位下三角矩阵，\(D\) 的对角元 \(d_i\) 均为正。令

    \[
    D^{1/2}
    =\operatorname{diag}(\sqrt{d_1},\ldots,\sqrt{d_n}),
    \qquad
    G=LD^{1/2}.
    \]

    则

    \[
    GG^{\mathsf T}
    =LD^{1/2}D^{1/2}L^{\mathsf T}
    =LDL^{\mathsf T}
    =A.
    \]

    正对角元的要求排除了符号不确定性；结合 \(LDL^{\mathsf T}\) 分解的唯一性，得到 \(G\) 的唯一性。 \(\square\)

### 4.2 Cholesky 算法

设 \(G=(g_{ij})\) 为下三角矩阵。由

\[
a_{ij}
=\sum_{s=1}^{\min(i,j)}g_{is}g_{js},
\]

可以逐列计算 \(G\)。

对角元满足

\[
g_{kk}
=\sqrt{
a_{kk}-\sum_{j=1}^{k-1}g_{kj}^2
},
\]

而第 \(k\) 列的非对角元满足

\[
g_{ik}
=\frac{
a_{ik}-\sum_{j=1}^{k-1}g_{ij}g_{kj}
}{g_{kk}},
\qquad
i=k+1,\ldots,n.
\]

!!! note "算法 4.2（逐列 Cholesky 分解）"

    对 \(k=1,\ldots,n\)，先计算

    \[
    g_{kk}
    =\sqrt{
    a_{kk}-\sum_{j=1}^{k-1}g_{kj}^2
    },
    \]

    再对 \(i=k+1,\ldots,n\) 计算

    \[
    g_{ik}
    =\frac{
    a_{ik}-\sum_{j=1}^{k-1}g_{ij}g_{kj}
    }{g_{kk}}.
    \]

Cholesky 分解约需

\[
\frac{1}{3}n^3+O(n^2)
\]

次浮点运算，约为一般 LU 分解的一半。

### 4.3 利用 Cholesky 分解求解方程组

由

\[
Ax=b
\quad\Longleftrightarrow\quad
GG^{\mathsf T}x=b,
\]

先通过向前代入求解

\[
Gy=b,
\]

再通过向后代入求解

\[
G^{\mathsf T}x=y.
\]

在标准浮点模型下，计算解 \(\widehat{x}\) 满足

\[
(A+E)\widehat{x}=b,
\qquad
\|E\|_2\leq c_nu\|A\|_2,
\]

其中 \(u\) 为单位舍入误差，\(c_n\) 是依赖于维数的温和常数。

!!! info "结论"

    对称正定系统应优先使用 Cholesky 分解：它不需要选主元，计算量和存储量更少，并且具有良好的向后稳定性。

### 4.4 对称正定平方根

!!! info "定义 4.3（矩阵平方根）"

    若矩阵 \(X\) 满足

    \[
    X^2=A,
    \]

    则称 \(X\) 为 \(A\) 的一个平方根。

矩阵平方根通常不唯一，但对称正定平方根是唯一的。

!!! success "定理 4.4（唯一的对称正定平方根）"

    若 \(A\) 为对称正定矩阵，则存在唯一的对称正定矩阵 \(A^{1/2}\)，使得

    \[
    (A^{1/2})^2=A.
    \]

由谱分解

\[
A=Q\Lambda Q^{\mathsf T},
\qquad
\Lambda=\operatorname{diag}(\lambda_1,\ldots,\lambda_n),
\]

其中所有 \(\lambda_i>0\)，可得

\[
A^{1/2}
=Q\Lambda^{1/2}Q^{\mathsf T}.
\]

---

## 5. 半正定矩阵与极分解

### 5.1 对称半正定矩阵的元素性质

!!! success "定理 5.1（半正定矩阵的元素界）"

    若 \(A\) 为对称半正定矩阵，则

    \[
    |a_{ij}|
    \leq\frac{a_{ii}+a_{jj}}{2},
    \]

    \[
    |a_{ij}|
    \leq\sqrt{a_{ii}a_{jj}},
    \]

    并且

    \[
    \max_{i,j}|a_{ij}|=\max_i a_{ii}.
    \]

    特别地，若 \(a_{ii}=0\)，则

    \[
    A[i,:]=A[:,i]=0.
    \]

??? proof "定理 5.1 的证明（点击展开）"

    对 \(x=e_i+e_j\) 与 \(x=e_i-e_j\) 使用 \(x^{\mathsf T}Ax\geq0\)，得到

    \[
    -\frac{a_{ii}+a_{jj}}{2}
    \leq a_{ij}
    \leq\frac{a_{ii}+a_{jj}}{2}.
    \]

    再对任意 \(\tau\in\mathbb{R}\) 取 \(x=\tau e_i+e_j\)，则

    \[
    a_{ii}\tau^2+2a_{ij}\tau+a_{jj}\geq0.
    \]

    该二次多项式的判别式非正，所以

    \[
    a_{ij}^2\leq a_{ii}a_{jj}.
    \]

    其余结论随即得到。 \(\square\)

### 5.2 极分解

!!! success "定理 5.2（极分解）"

    对任意 \(A\in\mathbb{R}^{n\times n}\)，存在正交矩阵 \(Q\) 与对称半正定矩阵 \(P\)，使得

    \[
    A=QP.
    \]

??? proof "定理 5.2 的证明（点击展开）"

    取 \(A\) 的奇异值分解

    \[
    A=U\Sigma V^{\mathsf T}.
    \]

    定义

    \[
    Q=UV^{\mathsf T},
    \qquad
    P=V\Sigma V^{\mathsf T}.
    \]

    则

    \[
    Q^{\mathsf T}Q=I,
    \qquad
    P=P^{\mathsf T}\succeq0,
    \qquad
    QP=A.
    \]

    \(\square\)

由于

\[
P^2=A^{\mathsf T}A,
\]

所以

\[
P=(A^{\mathsf T}A)^{1/2}.
\]

因此，极分解可以理解为先由 \(P\) 进行非负的伸缩，再由 \(Q\) 进行正交旋转或反射。

### 5.3 最佳正交逼近

!!! success "定理 5.3（最佳正交逼近）"

    若 \(A=QP\) 是极分解，则

    \[
    \|A-Q\|_F
    =\min_{R^{\mathsf T}R=I}\|A-R\|_F.
    \]

??? proof "定理 5.3 的证明思路（点击展开）"

    对任意正交矩阵 \(R\)，

    \[
    \|A-R\|_F^2
    =\operatorname{tr}(A^{\mathsf T}A)+n
    -2\operatorname{tr}(R^{\mathsf T}A).
    \]

    因此问题等价于最大化 \(\operatorname{tr}(R^{\mathsf T}A)\)。若

    \[
    A=U\Sigma V^{\mathsf T},
    \]

    令 \(Z=U^{\mathsf T}RV\)，则 \(Z\) 正交，并且

    \[
    \operatorname{tr}(R^{\mathsf T}A)
    =\operatorname{tr}(Z^{\mathsf T}\Sigma)
    \leq\sum_{i=1}^n\sigma_i.
    \]

    等号在 \(Z=I\)，即 \(R=UV^{\mathsf T}=Q\) 时取得。 \(\square\)

!!! example "正交 Procrustes 问题"

    对 \(A,B\in\mathbb{R}^{m\times n}\)，考虑

    \[
    \min_{Q^{\mathsf T}Q=I}\|A-BQ\|_F.
    \]

    该问题等价于

    \[
    \max_{Q^{\mathsf T}Q=I}
    \operatorname{tr}(Q^{\mathsf T}B^{\mathsf T}A).
    \]

    若

    \[
    B^{\mathsf T}A=U\Sigma V^{\mathsf T},
    \]

    则一个最优解为

    \[
    Q=UV^{\mathsf T}.
    \]

---

## 6. 带状线性系统

### 6.1 上带宽与下带宽

!!! info "定义 6.1（矩阵带宽）"

    若

    \[
    a_{ij}=0
    \qquad\text{whenever}\qquad
    j>i+q,
    \]

    则称 \(q\) 为 \(A\) 的上带宽，记为

    \[
    \operatorname{ub}(A)=q.
    \]

    若

    \[
    a_{ij}=0
    \qquad\text{whenever}\qquad
    i>j+p,
    \]

    则称 \(p\) 为 \(A\) 的下带宽，记为

    \[
    \operatorname{lb}(A)=p.
    \]

对角矩阵的上下带宽均为 \(0\)，三对角矩阵的上下带宽均为 \(1\)，一般稠密 \(n\times n\) 矩阵的上下带宽可达到 \(n-1\)。

### 6.2 带状三角方程组

若下三角矩阵 \(L\) 的下带宽为 \(p\)，则向前代入中每一行至多需要访问前 \(p\) 个未知量，因此复杂度为

\[
O(np).
\]

若上三角矩阵 \(U\) 的上带宽为 \(q\)，则向后代入的复杂度为

\[
O(nq).
\]

### 6.3 带状高斯消去

对上下带宽分别为 \(q\) 与 \(p\) 的矩阵进行无主元高斯消去，每一步只需更新带状区域，因此分解复杂度为

\[
O(npq).
\]

!!! success "定理 6.2（带状 LU 分解保持带宽）"

    假设 \(A\) 存在 LU 分解 \(A=LU\)，且

    \[
    \operatorname{ub}(A)=q,
    \qquad
    \operatorname{lb}(A)=p.
    \]

    则

    \[
    \operatorname{ub}(U)=q,
    \qquad
    \operatorname{lb}(L)=p.
    \]

??? proof "定理 6.2 的证明思路（点击展开）"

    对矩阵阶数作归纳。将 \(A\) 分块为

    \[
    A=
    \begin{bmatrix}
    \alpha & w^{\mathsf T}\\
    v & B
    \end{bmatrix}.
    \]

    第一列消元后的 Schur 补为

    \[
    B-\frac{vw^{\mathsf T}}{\alpha}.
    \]

    由于 \(v\) 只有前 \(p\) 个分量可能非零，\(w\) 只有前 \(q\) 个分量可能非零，秩一修正不会扩大相应带宽。对 Schur 补应用归纳假设，即得 \(L\) 与 \(U\) 继承原矩阵的下带宽和上带宽。 \(\square\)

### 6.4 选主元对带宽的影响

!!! success "定理 6.3（部分选主元后的上带宽）"

    若 \(A\) 非奇异，且

    \[
    \operatorname{ub}(A)=q,
    \qquad
\operatorname{lb}(A)=p.
    \]

    则采用按列部分选主元后，所得上三角因子满足

    \[
    \operatorname{ub}(U)\leq p+q.
    \]

行交换可能把下带区域中的元素移到主对角线以上，因此带宽可能从 \(q\) 增长到 \(p+q\)，但仍远小于一般稠密矩阵的带宽。

### 6.5 总计算复杂度

对于不需要选主元且保持带宽的系统，

\[
\text{LU 分解： }O(npq),
\]

\[
\text{求解 }Ly=b\text{： }O(np),
\]

\[
\text{求解 }Ux=y\text{： }O(nq).
\]

因此总复杂度为

\[
O\bigl(n(pq+p+q)\bigr).
\]

特别地，对于三对角矩阵，

\[
p=q=1,
\]

所以分解与求解总共只需要

\[
O(n)
\]

次运算。

!!! info "结构带来的收益"

    当 \(p,q\ll n\) 时，带状 LU、\(LDL^{\mathsf T}\) 和 Cholesky 分解在时间与存储上都远优于相应的稠密算法。

---

## 7. 本章总结

!!! success "核心结论"

    1. 严格对角占优保证矩阵可逆，并给出逆矩阵范数的直接上界。

    2. 非奇异列对角占优矩阵可以无主元地进行 LU 分解，且消元乘子的绝对值不超过 \(1\)。

    3. 对称矩阵的 \(LDM^{\mathsf T}\) 分解退化为 \(LDL^{\mathsf T}\) 分解。

    4. 对称正定矩阵具有唯一 Cholesky 分解 \(A=GG^{\mathsf T}\)，计算量约为一般 LU 分解的一半。

    5. 极分解 \(A=QP\) 把矩阵分成正交部分与对称半正定部分，并给出最佳正交逼近。

    6. 带状结构能够把分解复杂度从 \(O(n^3)\) 降为 \(O(npq)\)；三对角系统甚至可以在线性时间内求解。
