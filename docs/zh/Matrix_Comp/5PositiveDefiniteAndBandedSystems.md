# 第五章：正定矩阵、矩阵分解与带状系统

本章讨论正定与半正定矩阵、Cholesky 分解、极分解和矩阵平方根，并进一步研究带状矩阵在 LU 分解、部分选主元及求逆过程中的结构性质。

---

## 1. 正定矩阵

### 1.1 正定性的定义

!!! info "定义 1.1（正定矩阵）"

    对矩阵 $A\in\mathbb{R}^{n\times n}$，若对任意 $x\neq0$ 都有

    \[
    x^{\mathsf T}Ax>0,
    \]

    则称 $A$ 为正定矩阵。

矩阵 $A$ 不必预先假设为对称矩阵。令

\[
T=\frac{A+A^{\mathsf T}}{2},
\qquad
S=\frac{A-A^{\mathsf T}}{2},
\]

则 $T$ 对称、$S$ 反对称，并且对任意 $x$，

\[
x^{\mathsf T}Sx=0.
\]

因此

\[
x^{\mathsf T}Ax=x^{\mathsf T}Tx.
\]

!!! success "定理 1.2（正定性由对称部分决定）"

    $A$ 正定当且仅当其对称部分

    \[
    T=\frac{A+A^{\mathsf T}}{2}
    \]

    正定。

### 1.2 正定矩阵的 LU 分解

正定矩阵的每个顺序主子矩阵仍然正定，因而所有顺序主子式均为正。

!!! success "定理 1.3（正定矩阵的 LU 分解）"

    若 $A\in\mathbb{R}^{n\times n}$ 正定，则 $A$ 存在唯一的 LU 分解

    \[
    A=LU,
    \]

    其中 $L$ 为单位下三角矩阵，并且 $U$ 的所有对角元均为正。

??? proof "证明（点击展开）"

    因为 $A$ 的所有顺序主子式均非零，所以无主元高斯消去可以顺利进行，从而存在唯一的 LU 分解。进一步，

    \[
    UL^{-\mathsf T}=L^{-1}AL^{-\mathsf T}.
    \]

    右端与 $A$ 合同，故仍为正定矩阵。左端是上三角矩阵，其第 $i$ 个对角元恰为 $u_{ii}$，因此

    \[
    u_{ii}>0,
    \qquad i=1,\ldots,n.
    \]

    $\square$

!!! example "例 1.4"

    考虑

    \[
    A=
    \begin{bmatrix}
    \varepsilon & m\\
    -m & \varepsilon
    \end{bmatrix},
    \qquad \varepsilon>0.
    \]

    其对称部分为 $\varepsilon I$，因此 $A$ 正定。其 LU 分解为

    \[
    A=
    \begin{bmatrix}
    1&0\\
    -m/\varepsilon&1
    \end{bmatrix}
    \begin{bmatrix}
    \varepsilon&m\\
    0&\varepsilon+m^2/\varepsilon
    \end{bmatrix}.
    \]

    当 $|m/\varepsilon|$ 很大时，消元乘子也很大。这说明一般非对称正定矩阵虽然存在无主元 LU 分解，但数值稳定性仍可能受到反对称部分的影响。

---

## 2. Cholesky 分解

!!! success "定理 2.1（Cholesky 分解）"

    若 $A\in\mathbb{R}^{n\times n}$ 为对称正定矩阵，则存在唯一的下三角矩阵 $G$，使得

    \[
    A=GG^{\mathsf T},
    \]

    且 $G$ 的所有对角元均为正。

??? proof "证明（点击展开）"

    对称正定矩阵具有唯一分解

    \[
    A=LDL^{\mathsf T},
    \]

    其中 $L$ 为单位下三角矩阵，$D$ 为正对角矩阵。令

    \[
    G=LD^{1/2},
    \]

    则

    \[
    GG^{\mathsf T}=LDL^{\mathsf T}=A.
    \]

    $G$ 的对角元均为正，而正对角元条件与 $LDL^{\mathsf T}$ 分解的唯一性共同保证了 $G$ 的唯一性。 $\square$

!!! note "计算公式"

    若 $G=(g_{ij})$，则逐列 Cholesky 算法为

    \[
    g_{kk}=\sqrt{a_{kk}-\sum_{j=1}^{k-1}g_{kj}^2},
    \]

    以及

    \[
    g_{ik}=\frac{a_{ik}-\sum_{j=1}^{k-1}g_{ij}g_{kj}}{g_{kk}},
    \qquad i=k+1,\ldots,n.
    \]

Cholesky 分解约需

\[
\frac{1}{3}n^3+O(n^2)
\]

次浮点运算，约为一般 LU 分解的一半。

---

## 3. 对称半正定矩阵

### 3.1 定义与谱性质

!!! info "定义 3.1（对称半正定矩阵）"

    若 $A=A^{\mathsf T}\in\mathbb{R}^{n\times n}$，且对任意 $x\in\mathbb{R}^n$ 都有

    \[
    x^{\mathsf T}Ax\geq0,
    \]

    则称 $A$ 为对称半正定矩阵，记作 $A\succeq0$。

若 $\operatorname{rank}(A)=r$，则其特征值可排列为

\[
\lambda_1(A)\geq\cdots\geq\lambda_r(A)>0,
\]

以及

\[
\lambda_{r+1}(A)=\cdots=\lambda_n(A)=0.
\]

### 3.2 元素界

!!! success "定理 3.2（半正定矩阵的元素性质）"

    若 $A=(a_{ij})\succeq0$，则

    \[
    |a_{ij}|\leq\frac{a_{ii}+a_{jj}}{2},
    \]

    \[
    |a_{ij}|\leq\sqrt{a_{ii}a_{jj}},
    \]

    并且

    \[
    \max_{i,j}|a_{ij}|=\max_i a_{ii}.
    \]

    若 $a_{ii}=0$，则第 $i$ 行与第 $i$ 列全为零。

??? proof "证明（点击展开）"

    分别取 $x=e_i+e_j$ 与 $x=e_i-e_j$，由 $x^{\mathsf T}Ax\geq0$ 得

    \[
    a_{ii}+2a_{ij}+a_{jj}\geq0,
    \]

    \[
    a_{ii}-2a_{ij}+a_{jj}\geq0.
    \]

    因而

    \[
    |a_{ij}|\leq\frac{a_{ii}+a_{jj}}{2}.
    \]

    再取 $x=te_i+e_j$。对任意 $t\in\mathbb{R}$，

    \[
    a_{ii}t^2+2a_{ij}t+a_{jj}\geq0.
    \]

    该二次多项式的判别式非正，所以

    \[
    a_{ij}^2\leq a_{ii}a_{jj}.
    \]

    其余结论由此立即得到。 $\square$

---

## 4. 极分解与矩阵平方根

### 4.1 极分解

设 $A\in\mathbb{R}^{m\times n}$ 的秩为 $r$，其紧奇异值分解为

\[
A=U_1\Sigma_1V_1^{\mathsf T},
\]

其中

\[
\Sigma_1=\operatorname{diag}(\sigma_1,\ldots,\sigma_r),
\qquad \sigma_1\geq\cdots\geq\sigma_r>0.
\]

!!! success "定理 4.1（极分解）"

    $A$ 可以写成

    \[
    A=ZP,
    \]

    其中 $P\succeq0$，并且 $Z$ 在 $A$ 的有效子空间上保持长度。

    对方阵 $A\in\mathbb{R}^{n\times n}$，可取

    \[
    Z=UV^{\mathsf T},
    \qquad
    P=V\Sigma V^{\mathsf T},
    \]

    从而 $Z$ 正交且 $P\succeq0$。若 $A$ 非奇异，则 $P\succ0$。

由奇异值分解直接得到

\[
P^2=A^{\mathsf T}A,
\]

因此

\[
P=(A^{\mathsf T}A)^{1/2}.
\]

### 4.2 对称正定平方根

!!! success "定理 4.2（唯一的对称正定平方根）"

    若 $A$ 为对称正定矩阵，则存在唯一的对称正定矩阵 $X$，使得

    \[
    X^2=A.
    \]

由谱分解

\[
A=Q\Lambda Q^{\mathsf T}
\]

可构造

\[
A^{1/2}=Q\Lambda^{1/2}Q^{\mathsf T}.
\]

其唯一性可由正定平方根与 $A$ 具有相同的特征向量、且每个正特征值只有一个正平方根得到。

---

## 5. 带状矩阵与 LU 分解

### 5.1 带宽

!!! info "定义 5.1（上下带宽）"

    若当 $j>i+q$ 时总有 $a_{ij}=0$，则称 $A$ 的上带宽为 $q$。

    若当 $i>j+p$ 时总有 $a_{ij}=0$，则称 $A$ 的下带宽为 $p$。

### 5.2 无主元分解保持带宽

!!! success "定理 5.2（带状 LU 分解）"

    设 $A=LU$，且 $A$ 的上带宽与下带宽分别为 $q$ 和 $p$。则 $U$ 的上带宽为 $q$，$L$ 的下带宽为 $p$。

??? proof "证明（点击展开）"

    对矩阵阶数作归纳。写成分块形式

    \[
    A=
    \begin{bmatrix}
    \alpha&w^{\mathsf T}\\
    v&B
    \end{bmatrix}.
    \]

    第一步消元后的 Schur 补为

    \[
    B-\frac{vw^{\mathsf T}}{\alpha}.
    \]

    向量 $v$ 至多在前 $p$ 个位置非零，$w$ 至多在前 $q$ 个位置非零，所以秩一修正不会扩大原有带宽。对 Schur 补应用归纳假设即可。 $\square$

无主元带状 LU 分解的计算复杂度为

\[
O(npq),
\]

向前与向后代入的复杂度分别为 $O(np)$ 和 $O(nq)$。

### 5.3 部分选主元后的带宽

!!! success "定理 5.3（部分选主元产生的上带宽）"

    设 $A$ 非奇异，上带宽为 $q$，下带宽为 $p$。若使用按列部分选主元得到

    \[
    PA=LU,
    \]

    则

    \[
    \operatorname{ub}(U)\leq p+q.
    \]

??? proof "证明思路（点击展开）"

    设

    \[
    P=P_{n-1}\cdots P_1,
    \]

    并写

    \[
    P^{\mathsf T}=
    \begin{bmatrix}
    e_{s_1}&\cdots&e_{s_n}
    \end{bmatrix}.
    \]

    在第 $i$ 步，主元只能从第 $i$ 行到第 $i+p$ 行之间选取。否则，前 $i+p$ 行的前 $i$ 列会形成奇异的领先块，从而与 $A$ 非奇异矛盾。因此

    \[
    s_i\leq i+p.
    \]

    原矩阵第 $s_i$ 行在第 $s_i+q$ 列之后为零，所以置换后第 $i$ 行在第 $i+p+q$ 列之后为零。消元不会在更右侧产生填充，故

    \[
    \operatorname{ub}(U)\leq p+q.
    \]

    $\square$

---

## 6. 带状矩阵逆的块秩性质

一般而言，非奇异带状矩阵的逆是稠密矩阵，但其远离对角线的块具有低秩结构。

!!! success "定理 6.1（逆矩阵非对角块的秩界）"

    设非奇异带状矩阵分块为

    \[
    A=
    \begin{bmatrix}
    A_{11}&A_{12}\\
    A_{21}&A_{22}
    \end{bmatrix},
    \]

    其下带宽为 $p$，上带宽为 $q$。若

    \[
    A^{-1}=X=
    \begin{bmatrix}
    X_{11}&X_{12}\\
    X_{21}&X_{22}
    \end{bmatrix},
    \]

    则

    \[
    \operatorname{rank}(X_{21})\leq p,
    \]

    \[
    \operatorname{rank}(X_{12})\leq q.
    \]

??? proof "证明（点击展开）"

    先假设 $A_{22}$ 非奇异。由 $AX=I$ 的左下分块方程，

    \[
    A_{21}X_{11}+A_{22}X_{21}=0.
    \]

    因此

    \[
    X_{21}=-A_{22}^{-1}A_{21}X_{11}.
    \]

    由于 $A$ 的下带宽为 $p$，跨越分块边界的 $A_{21}$ 至多有 $p$ 个非零方向，所以

    \[
    \operatorname{rank}(A_{21})\leq p.
    \]

    从而

    \[
    \operatorname{rank}(X_{21})\leq p.
    \]

    同理由右上分块方程可得

    \[
    \operatorname{rank}(X_{12})\leq q.
    \]

    若 $A_{11}$ 或 $A_{22}$ 奇异，可用同带宽的任意小扰动使相应对角块非奇异，再令扰动趋于零。秩不超过给定整数的矩阵集合是闭集，因此极限仍满足相同秩界。 $\square$

!!! info "结构意义"

    带状矩阵的逆虽然通常不再带状，但其非对角块是低秩的。这一性质是半可分矩阵、分层矩阵以及快速直接法的重要基础。

---

## 7. 本章总结

!!! success "核心结论"

    1. 矩阵的正定性只取决于其对称部分。

    2. 正定矩阵存在唯一 LU 分解，且 $U$ 的对角元为正；对称正定矩阵存在唯一 Cholesky 分解。

    3. 对称半正定矩阵满足 $|a_{ij}|\leq\sqrt{a_{ii}a_{jj}}$。

    4. 极分解将矩阵写成等距变换与对称半正定伸缩的乘积。

    5. 无主元 LU 分解保持带宽；部分选主元后 $U$ 的上带宽至多增至 $p+q$。

    6. 带状矩阵的逆通常稠密，但其远离对角线的块具有由带宽控制的低秩结构。


