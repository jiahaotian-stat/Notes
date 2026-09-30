# 第三章：高斯消去与 LU 分解

本章讨论一般线性方程组的直接解法。我们先从三角方程组的向前代入与向后代入出发，再用高斯变换描述消元过程，并由此得到 LU 分解。最后讨论舍入误差、部分选主元以及带主元的 LU 分解。

---

## 1. 三角方程组

设

\[
A\in\mathbb{R}^{n\times n},
\qquad
b\in\mathbb{R}^n.
\]

我们希望求解线性方程组

\[
Ax=b.
\]

当 $A$ 为三角矩阵时，可以依次求出每个未知量，这构成一般消元法的基础。

### 1.1 下三角方程组与向前代入

设 $L=(l_{ij})$ 为非奇异下三角矩阵，则方程组 $Lx=b$ 可以写成

\[
\begin{bmatrix}
l_{11} & & & \\
l_{21} & l_{22} & & \\
\vdots & \vdots & \ddots & \\
l_{n1} & l_{n2} & \cdots & l_{nn}
\end{bmatrix}
\begin{bmatrix}
x_1\\x_2\\\vdots\\x_n
\end{bmatrix}
=
\begin{bmatrix}
b_1\\b_2\\\vdots\\b_n
\end{bmatrix}.
\]

从第一行开始逐行求解，得到向前代入公式

\[
x_1=\frac{b_1}{l_{11}},
\qquad
x_i=\frac{1}{l_{ii}}
\left(b_i-\sum_{j=1}^{i-1}l_{ij}x_j\right),
\quad i=2,\ldots,n.
\]

!!! note "算法 1.1（向前代入）"

    对 $i=1,\ldots,n$，依次计算

    \[
    x_i=\frac{b_i-l_{i,1:i-1}x_{1:i-1}}{l_{ii}}.
    \]

    若 $L$ 为单位下三角矩阵，则 $l_{ii}=1$，不需要除法。

### 1.2 上三角方程组与向后代入

设 $U=(u_{ij})$ 为非奇异上三角矩阵。从最后一行开始逐行向上求解，得到

\[
x_n=\frac{b_n}{u_{nn}},
\qquad
x_i=\frac{1}{u_{ii}}
\left(b_i-\sum_{j=i+1}^{n}u_{ij}x_j\right),
\quad i=n-1,\ldots,1.
\]

!!! note "算法 1.2（向后代入）"

    对 $i=n,n-1,\ldots,1$，依次计算

    \[
    x_i=\frac{b_i-u_{i,i+1:n}x_{i+1:n}}{u_{ii}}.
    \]

向前代入与向后代入都需要约 $n^2$ 次浮点运算，因此时间复杂度为 $O(n^2)$。

### 1.3 三角矩阵的代数性质

!!! success "命题 1.3（三角矩阵的封闭性）"

    上三角矩阵的乘积与逆仍为上三角矩阵；下三角矩阵的乘积与逆仍为下三角矩阵。

    若三角矩阵的对角元全为 $1$，则称其为单位三角矩阵。单位三角矩阵的乘积与逆仍为同类型的单位三角矩阵。

??? proof "命题 1.3 的证明（点击展开）"

    以上三角矩阵为例。若 $A$ 与 $B$ 都是上三角矩阵，则当 $i>j$ 时

    \[
    (AB)_{ij}=\sum_{k=1}^n a_{ik}b_{kj}=0,
    \]

    因为每一项至少有 $a_{ik}=0$ 或 $b_{kj}=0$。因此 $AB$ 仍为上三角矩阵。

    若 $A$ 非奇异，方程 $AX=I$ 的每一列都可由向后代入求得，并且所得 $X=A^{-1}$ 仍为上三角矩阵。下三角情形完全类似。若对角元全为 $1$，乘积与逆的对角元也全为 $1$。$\square$

### 1.4 三角求解的向后误差

设单位舍入误差为 $u$。浮点运算得到的向前代入解 $\widehat{x}$，可以看成某个邻近下三角系统的精确解：

\[
(L+\Delta L)\widehat{x}=b,
\qquad
|\Delta L|\leq \gamma_n|L|,
\]

其中

\[
\gamma_n=\frac{nu}{1-nu}=nu+O(u^2).
\]

类似地，向后代入满足

\[
(U+\Delta U)\widehat{x}=b,
\qquad
|\Delta U|\leq \gamma_n|U|.
\]

!!! info "解释"

    这说明三角求解是向后稳定的：计算结果恰好是一个与原问题非常接近的三角方程组的精确解。

---

## 2. 高斯消去

### 2.1 初等消元与高斯变换

设 $e_i$ 为第 $i$ 个标准基向量。左乘矩阵

\[
I+\lambda e_j e_i^{\mathsf T}
\]

等价于把矩阵的第 $i$ 行乘以 $\lambda$ 后加到第 $j$ 行。

在消去第 $1$ 列时，若 $a_{11}\neq0$，定义乘子

\[
\tau_i^{(1)}=\frac{a_{i1}}{a_{11}},
\qquad i=2,\ldots,n,
\]

以及向量

\[
\tau^{(1)}=
\begin{bmatrix}
0 & \tau_2^{(1)} & \cdots & \tau_n^{(1)}
\end{bmatrix}^{\mathsf T}.
\]

于是第一步高斯变换为

\[
M_1=I-\tau^{(1)}e_1^{\mathsf T},
\qquad
A^{(2)}=M_1A^{(1)},
\qquad
A^{(1)}=A.
\]

它把第 $1$ 列主元下方的元素全部消为零。

### 2.2 第 $k$ 步消元

假设前 $k-1$ 列已经消元。在第 $k$ 步，主元为 $a_{kk}^{(k)}$。若

\[
a_{kk}^{(k)}\neq0,
\]

定义

\[
\tau_i^{(k)}=
\frac{a_{ik}^{(k)}}{a_{kk}^{(k)}},
\quad i=k+1,\ldots,n,
\]

并令 $\tau_1^{(k)},\ldots,\tau_k^{(k)}$ 均为零。相应的高斯变换为

\[
M_k=I-\tau^{(k)}e_k^{\mathsf T},
\qquad
A^{(k+1)}=M_kA^{(k)}.
\]

完成 $n-1$ 步后，得到上三角矩阵

\[
M_{n-1}M_{n-2}\cdots M_1A=A^{(n)}=U.
\]

!!! warning "主元条件"

    不进行行交换的高斯消去要求每一步主元 $a_{kk}^{(k)}$ 都非零。矩阵非奇异本身并不能保证这一点。

### 2.3 高斯变换的性质

由于

\[
e_k^{\mathsf T}\tau^{(k)}=0,
\]

所以

\[
M_k^{-1}=I+\tau^{(k)}e_k^{\mathsf T}.
\]

事实上，

\[
\begin{aligned}
M_kM_k^{-1}
&=(I-\tau^{(k)}e_k^{\mathsf T})
  (I+\tau^{(k)}e_k^{\mathsf T})\\
&=I-\tau^{(k)}
  (e_k^{\mathsf T}\tau^{(k)})e_k^{\mathsf T}\\
&=I.
\end{aligned}
\]

当 $i<j$ 时，还有

\[
M_iM_j
=I-\tau^{(i)}e_i^{\mathsf T}
-\tau^{(j)}e_j^{\mathsf T},
\]

因为交叉项中 $e_i^{\mathsf T}\tau^{(j)}=0$。

### 2.4 消元更新与计算量

在第 $k$ 步，对 $i=k+1,ldots,n$，先计算乘子

\[
l_{ik}=\frac{a_{ik}^{(k)}}{a_{kk}^{(k)}},
\]

再更新余下元素

\[
a_{ij}^{(k+1)}
=a_{ij}^{(k)}-l_{ik}a_{kj}^{(k)},
\qquad
i,j=k+1,\ldots,n.
\]

高斯消去的主要计算量约为

\[
\frac{2}{3}n^3+O(n^2)
\]

次浮点运算，因此时间复杂度为 $O(n^3)$。

---

## 3. LU 分解

### 3.1 从高斯消去到 LU 分解

由

\[
M_{n-1}\cdots M_1A=U,
\]

可得

\[
A=M_1^{-1}M_2^{-1}\cdots M_{n-1}^{-1}U.
\]

定义

\[
L=M_1^{-1}M_2^{-1}\cdots M_{n-1}^{-1}.
\]

则 $L$ 为单位下三角矩阵，于是

\[
A=LU.
\]

进一步地，由高斯变换的结构可得

\[
L=I+\sum_{k=1}^{n-1}\tau^{(k)}e_k^{\mathsf T}.
\]

因此，$L$ 的严格下三角部分恰好保存消元过程中的乘子，而 $U$ 保存消元后的上三角部分。

### 3.2 LU 分解的存在性与唯一性

记 $A_k=A(1:k,1:k)$ 为 $A$ 的第 $k$ 个顺序主子矩阵。

!!! success "定理 3.1（LU 分解）"

    若

    \[
    \det(A_k)\neq0,
    \qquad k=1,\ldots,n,
    \]

    则存在唯一的单位下三角矩阵 $L$ 和上三角矩阵 $U$，使得

    \[
    A=LU.
    \]

    此外，

    \[
    \det(A)=\prod_{i=1}^n u_{ii}.
    \]

??? proof "定理 3.1 的证明（点击展开）"

    对矩阵阶数作归纳。$n=1$ 时结论显然成立。设结论对 $n-1$ 阶矩阵成立，并将 $A$ 分块为

    \[
    A=
    \begin{bmatrix}
    A_1 & v\\
    w^{\mathsf T} & \alpha
    \end{bmatrix}.
    \]

    由假设，$A_1$ 的所有顺序主子式均非零，所以存在唯一分解

    \[
    A_1=L_1U_1.
    \]

    于是

    \[
    A=
    \begin{bmatrix}
    L_1 & 0\\
    w^{\mathsf T}U_1^{-1} & 1
    \end{bmatrix}
    \begin{bmatrix}
    U_1 & L_1^{-1}v\\
    0 & \alpha-w^{\mathsf T}A_1^{-1}v
    \end{bmatrix}.
    \]

    这就证明了存在性。

    若还有

    \[
    A=L_1U_1=L_2U_2,
    \]

    则

    \[
    L_2^{-1}L_1=U_2U_1^{-1}.
    \]

    左边是单位下三角矩阵，右边是上三角矩阵，二者只能同时为单位矩阵。因此 $L_1=L_2$ 且 $U_1=U_2$，唯一性得证。

    最后，由 $\det(L)=1$，有

    \[
    \det(A)=\det(L)\det(U)=\prod_{i=1}^n u_{ii}.
    \]

    $\square$

!!! note "必要性说明"

    顺序主子式非零是无主元交换 LU 分解存在的常用充分条件。等价地，只要消元过程中所有主元均非零，就可以完成这种 LU 分解。

### 3.3 紧凑存储与计算算法

由于 $L$ 的对角元恒为 $1$，可以把 $L$ 的严格下三角部分与 $U$ 的上三角部分同时存放在原矩阵的空间中：

\[
\begin{bmatrix}
u_{11} & u_{12} & \cdots & u_{1n}\\
l_{21} & u_{22} & \cdots & u_{2n}\\
\vdots & \vdots & \ddots & \vdots\\
l_{n1} & l_{n2} & \cdots & u_{nn}
\end{bmatrix}.
\]

!!! note "算法 3.2（原位 LU 分解）"

    对 $k=1,\ldots,n-1$：

    先计算

    \[
    a_{k+1:n,k}\leftarrow
    \frac{a_{k+1:n,k}}{a_{kk}},
    \]

    再进行秩一更新

    \[
    A_{k+1:n,k+1:n}\leftarrow
    A_{k+1:n,k+1:n}
    -a_{k+1:n,k}a_{k,k+1:n}.
    \]

    算法结束后，严格下三角部分存储 $L$ 的乘子，上三角部分存储 $U$。

### 3.4 用 LU 分解求解方程组

一旦得到 $A=LU$，原方程组等价于

\[
LUx=b.
\]

引入中间变量 $y=Ux$，先通过向前代入求解

\[
Ly=b,
\]

再通过向后代入求解

\[
Ux=y.
\]

!!! info "多右端项情形"

    LU 分解只需进行一次，代价为 $O(n^3)$；此后每增加一个右端项，只需两次三角求解，代价为 $O(n^2)$。

### 3.5 $LDM^{\mathsf T}$ 分解

若 $A=LU$ 且 $U$ 的对角元均非零，令

\[
D=\operatorname{diag}(u_{11},\ldots,u_{nn}),
\qquad
M^{\mathsf T}=D^{-1}U.
\]

则 $M$ 为单位下三角矩阵，并且

\[
A=LDM^{\mathsf T}.
\]

!!! success "推论 3.3（$LDM^{\mathsf T}$ 分解）"

    在定理 3.1 的条件下，$A$ 存在唯一分解

    \[
    A=LDM^{\mathsf T},
    \]

    其中 $L,M$ 为单位下三角矩阵，$D$ 为非奇异对角矩阵。

---

## 4. 舍入误差与选主元

### 4.1 LU 分解的向后误差

设浮点运算得到 $\widehat{L}$ 和 $\widehat{U}$。在标准舍入模型下，可以写成

\[
\widehat{L}\widehat{U}=A+H,
\]

且有逐元素误差界

\[
|H|\leq
2(n-1)u
\left(|A|+|\widehat{L}||\widehat{U}|\right)
+O(u^2).
\]

利用 $\widehat{L}$、$\widehat{U}$ 和两次三角求解得到的解 $\widehat{x}$ 满足

\[
(A+E)\widehat{x}=b,
\]

其中

\[
|E|\leq
nu\left(2|A|+4|\widehat{L}||\widehat{U}|\right)
+O(u^2).
\]

!!! warning "小主元的危险"

    误差界依赖于 $|\widehat{L}||\widehat{U}|$。如果某一步主元很小，消元乘子可能很大，从而使 $L$ 和 $U$ 的元素显著增长，放大舍入误差。

例如，对

\[
A=
\begin{bmatrix}
0.001 & 1\\
1 & 0
\end{bmatrix},
\]

不交换行时得到

\[
L=
\begin{bmatrix}
1 & 0\\
1000 & 1
\end{bmatrix},
\qquad
U=
\begin{bmatrix}
0.001 & 1\\
0 & -1000
\end{bmatrix}.
\]

这说明非零但很小的主元同样可能造成数值不稳定。

### 4.2 对换矩阵

交换第 $i$ 行与第 $j$ 行的对换矩阵可写为

\[
P_{ij}
=I-(e_i-e_j)(e_i-e_j)^{\mathsf T}.
\]

它满足

\[
P_{ij}^{\mathsf T}=P_{ij}^{-1}=P_{ij},
\qquad
P_{ij}^2=I.
\]

左乘 $P_{ij}$ 等价于交换矩阵的第 $i$ 行与第 $j$ 行。

### 4.3 部分选主元

在第 $k$ 步，从当前第 $k$ 列的第 $k$ 行到第 $n$ 行中选择绝对值最大的元素：

\[
p_k=\operatorname*{arg\,max}_{k\leq i\leq n}
|a_{ik}^{(k)}|.
\]

随后交换第 $k$ 行与第 $p_k$ 行，再进行消元。于是所有乘子都满足

\[
|l_{ik}|\leq1,
\qquad i>k.
\]

!!! note "算法 4.1（带部分选主元的高斯消去）"

    对 $k=1,\ldots,n-1$：

    选择

    \[
    p_k=\operatorname*{arg\,max}_{k\leq i\leq n}|a_{ik}|;
    \]

    交换第 $k$ 行与第 $p_k$ 行；

    计算乘子并更新 Schur 补：

    \[
    l_{ik}=\frac{a_{ik}}{a_{kk}},
    \qquad
    a_{ij}\leftarrow a_{ij}-l_{ik}a_{kj}.
    \]

多次行交换的乘积构成置换矩阵 $P$，最终得到

\[
PA=LU.
\]

求解 $Ax=b$ 时，应依次求解

\[
Ly=Pb,
\qquad
Ux=y.
\]

### 4.4 复杂度与稳定性

部分选主元增加约 $O(n^2)$ 次比较，但主要消元运算仍为 $O(n^3)$。由于 $|l_{ik}|\leq1$，有

\[
\|L\|_{\infty}\leq n.
\]

这控制了消元乘子的增长。虽然 $U$ 在极端构造下仍可能出现较大的元素增长，但这种方法在实际计算中通常稳定，因此是求解一般稠密线性方程组的标准方法。

!!! info "完全选主元"

    完全选主元在剩余子矩阵中选择绝对值最大的元素，并同时交换相应的行与列。它能提供更强的增长控制，但搜索需要 $O(n^3)$ 次比较，因此实际中通常优先使用部分选主元。

---

## 5. 本章总结

!!! success "核心结论"

    1. 三角方程组可以用向前代入或向后代入在 $O(n^2)$ 时间内求解。

    2. 高斯变换把消元过程写成矩阵乘法，其逆矩阵中的向量正好记录消元乘子。

    3. 当消元主元均非零时，矩阵具有分解 $A=LU$；在单位下三角规范下，该分解唯一。

    4. LU 分解的主要计算量为 $O(n^3)$，但每个新右端项只需 $O(n^2)$ 的三角求解。

    5. 为避免零主元或过小主元，实际计算通常采用部分选主元，得到 $PA=LU$。
