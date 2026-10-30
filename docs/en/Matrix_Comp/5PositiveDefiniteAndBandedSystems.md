# Chapter 5: Positive Definite Matrices, Factorizations, and Banded Systems

This chapter studies positive definite and positive semidefinite matrices, Cholesky factorization, polar decomposition, and matrix square roots. It then develops the structural properties of banded matrices under LU factorization, partial pivoting, and inversion.

---

## 1. Positive Definite Matrices

### 1.1 Definition of Positive Definiteness

!!! info "Definition 1.1 (Positive Definite Matrix)"

    A matrix $A\in\mathbb{R}^{n\times n}$ is called positive definite if

    \[
    x^{\mathsf T}Ax>0,
    \]

    for every $x\neq0$.

The matrix $A$ is not assumed to be symmetric. Let

\[
T=\frac{A+A^{\mathsf T}}{2},
\qquad
S=\frac{A-A^{\mathsf T}}{2},
\]

Then $T$ is symmetric, $S$ is skew-symmetric, and

\[
x^{\mathsf T}Sx=0.
\]

for every $x$. Hence

\[
x^{\mathsf T}Ax=x^{\mathsf T}Tx.
\]

!!! success "Theorem 1.2 (Positive Definiteness Is Determined by the Symmetric Part)"

    $A$ is positive definite if and only if its symmetric part

    \[
    T=\frac{A+A^{\mathsf T}}{2}
    \]

    is positive definite.

### 1.2 LU Factorization of a Positive Definite Matrix

Every leading principal submatrix of a positive definite matrix is positive definite. Consequently, all leading principal minors are positive.

!!! success "Theorem 1.3 (LU Factorization of a Positive Definite Matrix)"

    If $A\in\mathbb{R}^{n\times n}$ is positive definite, then it has a unique LU factorization

    \[
    A=LU,
    \]

    where $L$ is unit lower triangular and every diagonal entry of $U$ is positive.

??? proof "Proof (click to expand)"

    Since every leading principal minor of $A$ is nonzero, Gaussian elimination without pivoting can be completed, giving a unique LU factorization. Moreover,

    \[
    UL^{-\mathsf T}=L^{-1}AL^{-\mathsf T}.
    \]

    The right-hand side is congruent to $A$, so it remains positive definite. The left-hand side is upper triangular and its $i$th diagonal entry is $u_{ii}$. Therefore

    \[
    u_{ii}>0,
    \qquad i=1,\ldots,n.
    \]

    $\square$

!!! example "Example 1.4"

    Consider

    \[
    A=
    \begin{bmatrix}
    \varepsilon & m\\
    -m & \varepsilon
    \end{bmatrix},
    \qquad \varepsilon>0.
    \]

    Its symmetric part is $\varepsilon I$, so $A$ is positive definite. Its LU factorization is

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

    When $|m/\varepsilon|$ is large, the elimination multiplier is also large. Thus, although a nonsymmetric positive definite matrix admits an LU factorization without pivoting, its numerical stability may still depend on the size of its skew-symmetric part.

---

## 2. Cholesky Factorization

!!! success "Theorem 2.1 (Cholesky Factorization)"

    If $A\in\mathbb{R}^{n\times n}$ is symmetric positive definite, then there exists a unique lower triangular matrix $G$ such that

    \[
    A=GG^{\mathsf T},
    \]

    and all diagonal entries of $G$ are positive.

??? proof "Proof (click to expand)"

    A symmetric positive definite matrix has the unique factorization

    \[
    A=LDL^{\mathsf T},
    \]

    where $L$ is unit lower triangular and $D$ is positive diagonal. Set

    \[
    G=LD^{1/2},
    \]

    Then

    \[
    GG^{\mathsf T}=LDL^{\mathsf T}=A.
    \]

    The diagonal entries of $G$ are positive. This sign convention, together with the uniqueness of the $LDL^{\mathsf T}$ factorization, gives uniqueness of $G$. $\square$

!!! note "Computational Formulas"

    If $G=(g_{ij})$, the column-oriented Cholesky algorithm uses

    \[
    g_{kk}=\sqrt{a_{kk}-\sum_{j=1}^{k-1}g_{kj}^2},
    \]

    followed by

    \[
    g_{ik}=\frac{a_{ik}-\sum_{j=1}^{k-1}g_{ij}g_{kj}}{g_{kk}},
    \qquad i=k+1,\ldots,n.
    \]

The Cholesky factorization requires approximately

\[
\frac{1}{3}n^3+O(n^2)
\]

floating-point operations, about half the cost of a general LU factorization.

---

## 3. Symmetric Positive Semidefinite Matrices

### 3.1 Definition and Spectral Properties

!!! info "Definition 3.1 (Symmetric Positive Semidefinite Matrix)"

    A symmetric matrix $A=A^{\mathsf T}\in\mathbb{R}^{n\times n}$ is called positive semidefinite, written $A\succeq0$, if

    \[
    x^{\mathsf T}Ax\geq0,
    \]

    for every $x\in\mathbb{R}^n$.

If $\operatorname{rank}(A)=r$, its eigenvalues can be ordered so that

\[
\lambda_1(A)\geq\cdots\geq\lambda_r(A)>0,
\]

and

\[
\lambda_{r+1}(A)=\cdots=\lambda_n(A)=0.
\]

### 3.2 Entrywise Bounds

!!! success "Theorem 3.2 (Entrywise Properties of a Positive Semidefinite Matrix)"

    If $A=(a_{ij})\succeq0$, then

    \[
    |a_{ij}|\leq\frac{a_{ii}+a_{jj}}{2},
    \]

    \[
    |a_{ij}|\leq\sqrt{a_{ii}a_{jj}},
    \]

    and

    \[
    \max_{i,j}|a_{ij}|=\max_i a_{ii}.
    \]

    If $a_{ii}=0$, then the entire $i$th row and $i$th column are zero.

??? proof "Proof (click to expand)"

    Apply $x^{\mathsf T}Ax\geq0$ to $x=e_i+e_j$ and $x=e_i-e_j$. This gives

    \[
    a_{ii}+2a_{ij}+a_{jj}\geq0,
    \]

    \[
    a_{ii}-2a_{ij}+a_{jj}\geq0.
    \]

    Therefore

    \[
    |a_{ij}|\leq\frac{a_{ii}+a_{jj}}{2}.
    \]

    Next take $x=te_i+e_j$. For every $t\in\mathbb{R}$,

    \[
    a_{ii}t^2+2a_{ij}t+a_{jj}\geq0.
    \]

    The discriminant of this quadratic polynomial is nonpositive, so

    \[
    a_{ij}^2\leq a_{ii}a_{jj}.
    \]

    The remaining claims follow immediately. $\square$

---

## 4. Polar Decomposition and Matrix Square Roots

### 4.1 Polar Decomposition

Let $A\in\mathbb{R}^{m\times n}$ have rank $r$, with compact singular value decomposition

\[
A=U_1\Sigma_1V_1^{\mathsf T},
\]

where

\[
\Sigma_1=\operatorname{diag}(\sigma_1,\ldots,\sigma_r),
\qquad \sigma_1\geq\cdots\geq\sigma_r>0.
\]

!!! success "Theorem 4.1 (Polar Decomposition)"

    The matrix $A$ can be written as

    \[
    A=ZP,
    \]

    where $P\succeq0$ and $Z$ preserves lengths on the effective subspace of $A$.

    For a square matrix $A\in\mathbb{R}^{n\times n}$, one may take

    \[
    Z=UV^{\mathsf T},
    \qquad
    P=V\Sigma V^{\mathsf T},
    \]

    Then $Z$ is orthogonal and $P\succeq0$. If $A$ is nonsingular, then $P\succ0$.

The singular value decomposition gives

\[
P^2=A^{\mathsf T}A,
\]

so

\[
P=(A^{\mathsf T}A)^{1/2}.
\]

### 4.2 Symmetric Positive Definite Square Root

!!! success "Theorem 4.2 (Unique Symmetric Positive Definite Square Root)"

    If $A$ is symmetric positive definite, then there exists a unique symmetric positive definite matrix $X$ such that

    \[
    X^2=A.
    \]

From the spectral decomposition

\[
A=Q\Lambda Q^{\mathsf T}
\]

we construct

\[
A^{1/2}=Q\Lambda^{1/2}Q^{\mathsf T}.
\]

Uniqueness follows because a positive definite square root shares eigenvectors with $A$, and every positive eigenvalue has a unique positive square root.

---

## 5. Banded Matrices and LU Factorization

### 5.1 Bandwidth

!!! info "Definition 5.1 (Upper and Lower Bandwidths)"

    If $a_{ij}=0$ whenever $j>i+q$, then $A$ has upper bandwidth $q$.

    If $a_{ij}=0$ whenever $i>j+p$, then $A$ has lower bandwidth $p$.

### 5.2 Factorization without Pivoting Preserves Bandwidth

!!! success "Theorem 5.2 (Banded LU Factorization)"

    Suppose $A=LU$, and $A$ has upper bandwidth $q$ and lower bandwidth $p$. Then $U$ has upper bandwidth $q$, while $L$ has lower bandwidth $p$.

??? proof "Proof (click to expand)"

    Proceed by induction on the matrix order. Partition

    \[
    A=
    \begin{bmatrix}
    \alpha&w^{\mathsf T}\\
    v&B
    \end{bmatrix}.
    \]

    After the first elimination step, the Schur complement is

    \[
    B-\frac{vw^{\mathsf T}}{\alpha}.
    \]

    At most the first $p$ entries of $v$ and the first $q$ entries of $w$ can be nonzero. Hence the rank-one correction cannot enlarge the original bandwidths. The result follows by applying the induction hypothesis to the Schur complement. $\square$

The complexity of banded LU factorization without pivoting is

\[
O(npq),
\]

Forward and backward substitution cost $O(np)$ and $O(nq)$, respectively.

### 5.3 Bandwidth after Partial Pivoting

!!! success "Theorem 5.3 (Upper Bandwidth after Partial Pivoting)"

    Let $A$ be nonsingular with upper bandwidth $q$ and lower bandwidth $p$. If columnwise partial pivoting produces

    \[
    PA=LU,
    \]

    then

    \[
    \operatorname{ub}(U)\leq p+q.
    \]

??? proof "Proof sketch (click to expand)"

    Write

    \[
    P=P_{n-1}\cdots P_1,
    \]

    and

    \[
    P^{\mathsf T}=
    \begin{bmatrix}
    e_{s_1}&\cdots&e_{s_n}
    \end{bmatrix}.
    \]

    At step $i$, the pivot can only be selected from rows $i$ through $i+p$. Otherwise the leading block formed by the first $i+p$ rows and first $i$ columns would be singular, contradicting nonsingularity of $A$. Thus

    \[
    s_i\leq i+p.
    \]

    Row $s_i$ of the original matrix vanishes after column $s_i+q$. Therefore, after permutation, row $i$ vanishes after column $i+p+q$. Elimination creates no fill farther to the right, so

    \[
    \operatorname{ub}(U)\leq p+q.
    \]

    $\square$

---

## 6. Block-Rank Properties of the Inverse of a Banded Matrix

The inverse of a nonsingular banded matrix is generally dense, but its off-diagonal blocks have low-rank structure.

!!! success "Theorem 6.1 (Rank Bounds for Off-Diagonal Blocks of the Inverse)"

    Partition a nonsingular banded matrix as

    \[
    A=
    \begin{bmatrix}
    A_{11}&A_{12}\\
    A_{21}&A_{22}
    \end{bmatrix},
    \]

    where the lower bandwidth is $p$ and the upper bandwidth is $q$. If

    \[
    A^{-1}=X=
    \begin{bmatrix}
    X_{11}&X_{12}\\
    X_{21}&X_{22}
    \end{bmatrix},
    \]

    then

    \[
    \operatorname{rank}(X_{21})\leq p,
    \]

    \[
    \operatorname{rank}(X_{12})\leq q.
    \]

??? proof "Proof (click to expand)"

    First suppose that $A_{22}$ is nonsingular. The lower-left block equation in $AX=I$ is

    \[
    A_{21}X_{11}+A_{22}X_{21}=0.
    \]

    Hence

    \[
    X_{21}=-A_{22}^{-1}A_{21}X_{11}.
    \]

    Since $A$ has lower bandwidth $p$, the block $A_{21}$, which crosses the partition, has at most $p$ independent nonzero directions. Therefore

    \[
    \operatorname{rank}(A_{21})\leq p.
    \]

    and consequently

    \[
    \operatorname{rank}(X_{21})\leq p.
    \]

    The upper-right block equation similarly gives

    \[
    \operatorname{rank}(X_{12})\leq q.
    \]

    If $A_{11}$ or $A_{22}$ is singular, perturb the corresponding diagonal block by an arbitrarily small matrix of the same bandwidth so that it becomes nonsingular, and then let the perturbation tend to zero. The set of matrices of rank at most a fixed integer is closed, so the same rank bounds hold in the limit. $\square$

!!! info "Structural Meaning"

    Although the inverse of a banded matrix is generally not banded, its off-diagonal blocks are low rank. This observation underlies semiseparable matrices, hierarchical matrices, and fast direct solvers.

---

## 7. Chapter Summary

!!! success "Main Conclusions"

    1. Positive definiteness is determined entirely by the symmetric part of a matrix.

    2. A positive definite matrix has a unique LU factorization with positive pivots; a symmetric positive definite matrix has a unique Cholesky factorization.

    3. A symmetric positive semidefinite matrix satisfies $|a_{ij}|\leq\sqrt{a_{ii}a_{jj}}$.

    4. Polar decomposition separates a matrix into an isometric transformation and a symmetric positive semidefinite scaling.

    5. LU factorization without pivoting preserves bandwidth; with partial pivoting, the upper bandwidth of $U$ is at most $p+q$.

    6. The inverse of a banded matrix is usually dense, but its off-diagonal blocks have ranks controlled by the bandwidths.


