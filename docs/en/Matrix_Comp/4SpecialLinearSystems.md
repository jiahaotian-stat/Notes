# Chapter 4: Special Linear Systems and Structured Factorizations

This chapter studies linear systems with special structure. Compared with general dense systems, these matrices often require no pivoting or permit substantial reductions in computational cost. We discuss diagonally dominant matrices, the \(LDL^{\mathsf T}\) factorization of symmetric matrices, positive definite systems, Cholesky factorization, polar decomposition, and banded linear systems.

---

## 1. Diagonally Dominant Matrices

### 1.1 Row and Column Diagonal Dominance

Let

\[
A=(a_{ij})\in\mathbb{R}^{n\times n}.
\]

!!! info "Definition 1.1 (Diagonal Dominance)"

    If, for every \(i=1,\ldots,n\),

    \[
    |a_{ii}|\geq\sum_{j\neq i}|a_{ij}|,
    \]

    then \(A\) is called row diagonally dominant.

    If, for every \(j=1,\ldots,n\),

    \[
    |a_{jj}|\geq\sum_{i\neq j}|a_{ij}|,
    \]

    then \(A\) is called column diagonally dominant.

    If the corresponding inequality is strict for every row or every column, then \(A\) is called strictly row diagonally dominant or strictly column diagonally dominant, respectively.

A diagonally dominant matrix need not be invertible. For example, the zero matrix is both row and column diagonally dominant but is clearly singular. Strict diagonal dominance, however, guarantees invertibility.

### 1.2 Invertibility of Strictly Diagonally Dominant Matrices

!!! success "Theorem 1.2 (Invertibility under Strict Diagonal Dominance)"

    Every strictly row diagonally dominant matrix and every strictly column diagonally dominant matrix is nonsingular.

??? proof "Proof for the strictly row diagonally dominant case (click to expand)"

    Suppose, to the contrary, that \(A\) is singular. Then there exists \(x\neq0\) such that

    \[
    Ax=0.
    \]

    Choose

    \[
    k\in\operatorname*{arg\,max}_{1\leq i\leq n}|x_i|.
    \]

    Since \(x\neq0\), we have \(|x_k|>0\). The \(k\)th equation gives

    \[
    a_{kk}x_k=-\sum_{j\neq k}a_{kj}x_j.
    \]

    Hence

    \[
    \begin{aligned}
    |a_{kk}||x_k|
    &\leq\sum_{j\neq k}|a_{kj}||x_j|\\
    &\leq\left(\sum_{j\neq k}|a_{kj}|\right)|x_k|.
    \end{aligned}
    \]

    Cancelling \(|x_k|\) yields

    \[
    |a_{kk}|\leq\sum_{j\neq k}|a_{kj}|,
    \]

    contradicting strict row diagonal dominance. Thus \(A\) is nonsingular. The strictly column diagonally dominant case follows by applying the same result to \(A^{\mathsf T}\). \(\square\)

### 1.3 Norm Bounds for the Inverse

Define the row dominance factor

\[
\delta
=\min_{1\leq i\leq n}
\left(
|a_{ii}|-\sum_{j\neq i}|a_{ij}|
\right),
\]

and the column dominance factor

\[
\mu
=\min_{1\leq j\leq n}
\left(
|a_{jj}|-\sum_{i\neq j}|a_{ij}|
\right).
\]

!!! success "Theorem 1.3 (Inverse Norm Bounds)"

    If \(A\) is strictly row diagonally dominant, then

    \[
    \|A^{-1}\|_{\infty}\leq\frac{1}{\delta}.
    \]

    If \(A\) is strictly column diagonally dominant, then

    \[
    \|A^{-1}\|_1\leq\frac{1}{\mu}.
    \]

??? proof "Proof of the infinity-norm bound (click to expand)"

    For any \(x\in\mathbb{R}^n\), choose \(k\) such that

    \[
    |x_k|=\|x\|_{\infty}.
    \]

    Then

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

    Setting \(x=A^{-1}y\) gives

    \[
    \|A^{-1}y\|_{\infty}
    \leq\frac{1}{\delta}\|y\|_{\infty}.
    \]

    Taking the supremum over all \(y\neq0\) proves the claim. The column diagonally dominant case follows from

    \[
    \|A^{-1}\|_1
    =\|A^{-\mathsf T}\|_{\infty}
    \]

    and the strict row diagonal dominance of \(A^{\mathsf T}\). \(\square\)

### 1.4 LU Factorization of Column Diagonally Dominant Matrices

!!! success "Theorem 1.4 (Column Diagonal Dominance and LU)"

    If \(A\) is nonsingular and column diagonally dominant, then it has a unique LU factorization

    \[
    A=LU,
    \]

    where \(L\) is unit lower triangular and

    \[
    |l_{ij}|\leq1,
    \qquad i>j.
    \]

??? proof "Proof sketch for Theorem 1.4 (click to expand)"

    We use induction on the order of the matrix. Partition \(A\) as

    \[
    A=
    \begin{bmatrix}
    \alpha & w^{\mathsf T}\\
    v & C
    \end{bmatrix}.
    \]

    Nonsingularity and column diagonal dominance guarantee \(\alpha\neq0\), and the first-step multipliers satisfy

    \[
    \left|\frac{v_i}{\alpha}\right|\leq1.
    \]

    The Schur complement after the first elimination step is

    \[
    B=C-\frac{vw^{\mathsf T}}{\alpha}.
    \]

    A columnwise verification shows that \(B\) remains nonsingular and column diagonally dominant. By the induction hypothesis, \(B\) has an LU factorization satisfying the same multiplier bound. Combining the first-step multipliers with this factorization gives \(A=LU\) with every \(|l_{ij}|\leq1\). Uniqueness follows from the unit lower triangular normalization. \(\square\)

!!! info "Numerical Meaning"

    Thus, when Gaussian elimination without pivoting is applied to a nonsingular column diagonally dominant matrix, the elimination multipliers are automatically controlled, so additional pivoting is usually unnecessary.

---

## 2. Symmetric Matrices and the \(LDL^{\mathsf T}\) Factorization

### 2.1 From \(LDM^{\mathsf T}\) to \(LDL^{\mathsf T}\)

Suppose a nonsingular matrix \(A\) has a factorization

\[
A=LDM^{\mathsf T},
\]

where \(L\) and \(M\) are unit lower triangular and \(D\) is nonsingular and diagonal.

!!! success "Theorem 2.1 (The \(LDL^{\mathsf T}\) Factorization of a Symmetric Matrix)"

    If \(A\) is symmetric, then

    \[
    M=L,
    \]

    and therefore

    \[
    A=LDL^{\mathsf T}.
    \]

??? proof "Proof of Theorem 2.1 (click to expand)"

    Since \(A=A^{\mathsf T}\),

    \[
    LDM^{\mathsf T}=MDL^{\mathsf T}.
    \]

    Multiplying on the left by \(M^{-1}\) and on the right by \(L^{-\mathsf T}\) gives

    \[
    M^{-1}LD=D(M^{-1}L)^{\mathsf T}.
    \]

    The left-hand side is lower triangular and the right-hand side is upper triangular, so both must be diagonal. Moreover, \(M^{-1}L\) is unit lower triangular, and hence

    \[
    M^{-1}L=I.
    \]

    Therefore \(M=L\). \(\square\)

### 2.2 Computational and Storage Advantages

Only half of a symmetric matrix needs to be stored. During elimination, the upper triangular part can be recovered from the lower triangular part by symmetry, so an \(LDL^{\mathsf T}\) factorization requires roughly half the work and storage of a general LU factorization.

If pivoting is required, the corresponding rows and columns must be interchanged simultaneously to preserve symmetry, giving

\[
PAP^{\mathsf T}=LDL^{\mathsf T}.
\]

!!! warning "Symmetry and Pivoting"

    Interchanging rows alone destroys the symmetric structure. Symmetric indefinite systems usually require a symmetric pivoting strategy and may use \(1\times1\) or \(2\times2\) diagonal blocks.

---

## 3. Positive Definite Systems

### 3.1 Positive Definiteness and Positive Semidefiniteness

!!! info "Definition 3.1 (Definiteness)"

    If, for every \(x\neq0\),

    \[
    x^{\mathsf T}Ax>0,
    \]

    then \(A\) is called positive definite.

    If, for every \(x\),

    \[
    x^{\mathsf T}Ax\geq0,
    \]

    then \(A\) is called positive semidefinite.

Symmetry is not initially assumed here. A quadratic form depends only on the symmetric part of \(A\), as shown below.

### 3.2 Principal Submatrices Preserve Positive Definiteness

!!! success "Theorem 3.2 (Principal Submatrices of a Positive Definite Matrix)"

    Every principal submatrix of a positive definite matrix is positive definite.

??? proof "Proof of Theorem 3.2 (click to expand)"

    Choose indices

    \[
    1\leq i_1<\cdots<i_k\leq n,
    \]

    and let

    \[
    X=
    \begin{bmatrix}
    e_{i_1} & \cdots & e_{i_k}
    \end{bmatrix}.
    \]

    The corresponding principal submatrix is

    \[
    B=X^{\mathsf T}AX.
    \]

    For every \(z\neq0\), the full column rank of \(X\) implies \(Xz\neq0\), and hence

    \[
    z^{\mathsf T}Bz
    =(Xz)^{\mathsf T}A(Xz)>0.
    \]

    Thus \(B\) is positive definite. \(\square\)

A direct consequence is that all leading principal minors of a positive definite matrix are positive.

### 3.3 LU Factorization of a Positive Definite Matrix

!!! success "Theorem 3.3 (LU Factorization of a Positive Definite Matrix)"

    If \(A\in\mathbb{R}^{n\times n}\) is positive definite, then it has a unique LU factorization

    \[
    A=LU,
    \]

    and every diagonal entry of \(U\) is positive.

??? proof "Proof of Theorem 3.3 (click to expand)"

    By Theorem 3.2, every leading principal submatrix of \(A\) is positive definite, so every leading principal minor is positive. The existence and uniqueness theorem for LU factorization therefore gives a unique LU factorization of \(A\).

    Moreover,

    \[
    UL^{-\mathsf T}
    =L^{-1}AL^{-\mathsf T}.
    \]

    The right-hand side remains positive definite under the congruence transformation, so \(UL^{-\mathsf T}\) is positive definite. Since both \(U\) and \(L^{-\mathsf T}\) are upper triangular,

    \[
    (UL^{-\mathsf T})_{ii}
    =u_{ii}>0.
    \]

    \(\square\)

### 3.4 Symmetric and Skew-Symmetric Parts

Every real matrix has the unique decomposition

\[
A=T+S,
\]

where

\[
T=\frac{A+A^{\mathsf T}}{2},
\qquad
S=\frac{A-A^{\mathsf T}}{2}.
\]

The matrix \(T\) is symmetric, while \(S\) is skew-symmetric.

!!! success "Proposition 3.4 (A Quadratic Form Depends Only on the Symmetric Part)"

    For every \(x\in\mathbb{R}^n\),

    \[
    x^{\mathsf T}Sx=0,
    \qquad
    x^{\mathsf T}Ax=x^{\mathsf T}Tx.
    \]

    Therefore \(A\) is positive definite if and only if its symmetric part \(T\) is positive definite.

??? proof "Proof of Proposition 3.4 (click to expand)"

    Since \(S^{\mathsf T}=-S\), the scalar \(x^{\mathsf T}Sx\) satisfies

    \[
    x^{\mathsf T}Sx
    =(x^{\mathsf T}Sx)^{\mathsf T}
    =x^{\mathsf T}S^{\mathsf T}x
    =-x^{\mathsf T}Sx.
    \]

    Thus \(x^{\mathsf T}Sx=0\), and consequently

    \[
    x^{\mathsf T}Ax
    =x^{\mathsf T}(T+S)x
    =x^{\mathsf T}Tx.
    \]

    \(\square\)

!!! note "Stability Remark"

    The skew-symmetric part of a symmetric positive definite matrix vanishes, so its factorization without pivoting usually has excellent numerical stability. For a general nonsymmetric positive definite matrix, stability also depends on the size of the skew-symmetric part relative to the symmetric part.

---

## 4. Cholesky Factorization

### 4.1 Existence and Uniqueness

!!! success "Theorem 4.1 (Cholesky Factorization)"

    If \(A\in\mathbb{R}^{n\times n}\) is symmetric positive definite, then there exists a unique lower triangular matrix \(G\) such that

    \[
    A=GG^{\mathsf T},
    \]

    and every diagonal entry of \(G\) is positive.

??? proof "Proof of Theorem 4.1 (click to expand)"

    A symmetric positive definite matrix has a unique factorization

    \[
    A=LDL^{\mathsf T},
    \]

    where \(L\) is unit lower triangular and the diagonal entries \(d_i\) of \(D\) are positive. Define

    \[
    D^{1/2}
    =\operatorname{diag}(\sqrt{d_1},\ldots,\sqrt{d_n}),
    \qquad
    G=LD^{1/2}.
    \]

    Then

    \[
    GG^{\mathsf T}
    =LD^{1/2}D^{1/2}L^{\mathsf T}
    =LDL^{\mathsf T}
    =A.
    \]

    Requiring positive diagonal entries removes the sign ambiguity. Together with uniqueness of the \(LDL^{\mathsf T}\) factorization, this proves uniqueness of \(G\). \(\square\)

### 4.2 The Cholesky Algorithm

Let \(G=(g_{ij})\) be lower triangular. From

\[
a_{ij}
=\sum_{s=1}^{\min(i,j)}g_{is}g_{js},
\]

we can compute \(G\) one column at a time.

The diagonal entry is

\[
g_{kk}
=\sqrt{
a_{kk}-\sum_{j=1}^{k-1}g_{kj}^2
},
\]

and the off-diagonal entries in column \(k\) are

\[
g_{ik}
=\frac{
a_{ik}-\sum_{j=1}^{k-1}g_{ij}g_{kj}
}{g_{kk}},
\qquad
i=k+1,\ldots,n.
\]

!!! note "Algorithm 4.2 (Columnwise Cholesky Factorization)"

    For \(k=1,\ldots,n\), first compute

    \[
    g_{kk}
    =\sqrt{
    a_{kk}-\sum_{j=1}^{k-1}g_{kj}^2
    },
    \]

    and then, for \(i=k+1,\ldots,n\), compute

    \[
    g_{ik}
    =\frac{
    a_{ik}-\sum_{j=1}^{k-1}g_{ij}g_{kj}
    }{g_{kk}}.
    \]

The Cholesky factorization requires approximately

\[
\frac{1}{3}n^3+O(n^2)
\]

floating-point operations, roughly half the work of a general LU factorization.

### 4.3 Solving Linear Systems with Cholesky

From

\[
Ax=b
\quad\Longleftrightarrow\quad
GG^{\mathsf T}x=b,
\]

first solve by forward substitution

\[
Gy=b,
\]

and then solve by backward substitution

\[
G^{\mathsf T}x=y.
\]

Under the standard floating-point model, the computed solution \(\widehat{x}\) satisfies

\[
(A+E)\widehat{x}=b,
\qquad
\|E\|_2\leq c_nu\|A\|_2,
\]

where \(u\) is the unit roundoff and \(c_n\) is a modest dimension-dependent constant.

!!! info "Conclusion"

    Cholesky factorization should be preferred for symmetric positive definite systems: it requires no pivoting, uses less computation and storage, and has good backward stability.

### 4.4 The Symmetric Positive Definite Square Root

!!! info "Definition 4.3 (Matrix Square Root)"

    If a matrix \(X\) satisfies

    \[
    X^2=A,
    \]

    then \(X\) is called a square root of \(A\).

Matrix square roots are generally not unique, but the symmetric positive definite square root is unique.

!!! success "Theorem 4.4 (Unique Symmetric Positive Definite Square Root)"

    If \(A\) is symmetric positive definite, then there exists a unique symmetric positive definite matrix \(A^{1/2}\) such that

    \[
    (A^{1/2})^2=A.
    \]

From the spectral decomposition

\[
A=Q\Lambda Q^{\mathsf T},
\qquad
\Lambda=\operatorname{diag}(\lambda_1,\ldots,\lambda_n),
\]

where every \(\lambda_i>0\), we obtain

\[
A^{1/2}
=Q\Lambda^{1/2}Q^{\mathsf T}.
\]

---

## 5. Positive Semidefinite Matrices and Polar Decomposition

### 5.1 Entrywise Properties of Symmetric Positive Semidefinite Matrices

!!! success "Theorem 5.1 (Entrywise Bounds for a Positive Semidefinite Matrix)"

    If \(A\) is symmetric positive semidefinite, then

    \[
    |a_{ij}|
    \leq\frac{a_{ii}+a_{jj}}{2},
    \]

    \[
    |a_{ij}|
    \leq\sqrt{a_{ii}a_{jj}},
    \]

    and

    \[
    \max_{i,j}|a_{ij}|=\max_i a_{ii}.
    \]

    In particular, if \(a_{ii}=0\), then

    \[
    A[i,:]=A[:,i]=0.
    \]

??? proof "Proof of Theorem 5.1 (click to expand)"

    Applying \(x^{\mathsf T}Ax\geq0\) to \(x=e_i+e_j\) and \(x=e_i-e_j\) gives

    \[
    -\frac{a_{ii}+a_{jj}}{2}
    \leq a_{ij}
    \leq\frac{a_{ii}+a_{jj}}{2}.
    \]

    Next, for arbitrary \(\tau\in\mathbb{R}\), take \(x=\tau e_i+e_j\). Then

    \[
    a_{ii}\tau^2+2a_{ij}\tau+a_{jj}\geq0.
    \]

    The discriminant of this quadratic polynomial must be nonpositive, so

    \[
    a_{ij}^2\leq a_{ii}a_{jj}.
    \]

    The remaining statements follow immediately. \(\square\)

### 5.2 Polar Decomposition

!!! success "Theorem 5.2 (Polar Decomposition)"

    For every \(A\in\mathbb{R}^{n\times n}\), there exist an orthogonal matrix \(Q\) and a symmetric positive semidefinite matrix \(P\) such that

    \[
    A=QP.
    \]

??? proof "Proof of Theorem 5.2 (click to expand)"

    Take a singular value decomposition

    \[
    A=U\Sigma V^{\mathsf T}.
    \]

    Define

    \[
    Q=UV^{\mathsf T},
    \qquad
    P=V\Sigma V^{\mathsf T}.
    \]

    Then

    \[
    Q^{\mathsf T}Q=I,
    \qquad
    P=P^{\mathsf T}\succeq0,
    \qquad
    QP=A.
    \]

    \(\square\)

Since

\[
P^2=A^{\mathsf T}A,
\]

we have

\[
P=(A^{\mathsf T}A)^{1/2}.
\]

Thus the polar decomposition can be interpreted as a nonnegative scaling by \(P\), followed by an orthogonal rotation or reflection by \(Q\).

### 5.3 Best Orthogonal Approximation

!!! success "Theorem 5.3 (Best Orthogonal Approximation)"

    If \(A=QP\) is a polar decomposition, then

    \[
    \|A-Q\|_F
    =\min_{R^{\mathsf T}R=I}\|A-R\|_F.
    \]

??? proof "Proof sketch for Theorem 5.3 (click to expand)"

    For any orthogonal matrix \(R\),

    \[
    \|A-R\|_F^2
    =\operatorname{tr}(A^{\mathsf T}A)+n
    -2\operatorname{tr}(R^{\mathsf T}A).
    \]

    Thus the problem is equivalent to maximizing \(\operatorname{tr}(R^{\mathsf T}A)\). If

    \[
    A=U\Sigma V^{\mathsf T},
    \]

    let \(Z=U^{\mathsf T}RV\). Then \(Z\) is orthogonal, and

    \[
    \operatorname{tr}(R^{\mathsf T}A)
    =\operatorname{tr}(Z^{\mathsf T}\Sigma)
    \leq\sum_{i=1}^n\sigma_i.
    \]

    Equality holds when \(Z=I\), or equivalently when \(R=UV^{\mathsf T}=Q\). \(\square\)

!!! example "The Orthogonal Procrustes Problem"

    For \(A,B\in\mathbb{R}^{m\times n}\), consider

    \[
    \min_{Q^{\mathsf T}Q=I}\|A-BQ\|_F.
    \]

    This problem is equivalent to

    \[
    \max_{Q^{\mathsf T}Q=I}
    \operatorname{tr}(Q^{\mathsf T}B^{\mathsf T}A).
    \]

    If

    \[
    B^{\mathsf T}A=U\Sigma V^{\mathsf T},
    \]

    then one optimal solution is

    \[
    Q=UV^{\mathsf T}.
    \]

---

## 6. Banded Linear Systems

### 6.1 Upper and Lower Bandwidths

!!! info "Definition 6.1 (Matrix Bandwidth)"

    If

    \[
    a_{ij}=0
    \qquad\text{whenever}\qquad
    j>i+q,
    \]

    then \(q\) is called the upper bandwidth of \(A\), denoted by

    \[
    \operatorname{ub}(A)=q.
    \]

    If

    \[
    a_{ij}=0
    \qquad\text{whenever}\qquad
    i>j+p,
    \]

    then \(p\) is called the lower bandwidth of \(A\), denoted by

    \[
    \operatorname{lb}(A)=p.
    \]

A diagonal matrix has upper and lower bandwidth \(0\), a tridiagonal matrix has upper and lower bandwidth \(1\), and a general dense \(n\times n\) matrix can have both bandwidths as large as \(n-1\).

### 6.2 Banded Triangular Systems

If a lower triangular matrix \(L\) has lower bandwidth \(p\), each row of forward substitution accesses at most the preceding \(p\) unknowns, so the complexity is

\[
O(np).
\]

If an upper triangular matrix \(U\) has upper bandwidth \(q\), the complexity of backward substitution is

\[
O(nq).
\]

### 6.3 Banded Gaussian Elimination

For a matrix with upper bandwidth \(q\) and lower bandwidth \(p\), Gaussian elimination without pivoting updates only entries in the band, so the factorization complexity is

\[
O(npq).
\]

!!! success "Theorem 6.2 (Banded LU Factorization Preserves Bandwidth)"

    Suppose \(A\) has an LU factorization \(A=LU\), and

    \[
    \operatorname{ub}(A)=q,
    \qquad
    \operatorname{lb}(A)=p.
    \]

    Then

    \[
    \operatorname{ub}(U)=q,
    \qquad
    \operatorname{lb}(L)=p.
    \]

??? proof "Proof sketch for Theorem 6.2 (click to expand)"

    We use induction on the matrix order. Partition \(A\) as

    \[
    A=
    \begin{bmatrix}
    \alpha & w^{\mathsf T}\\
    v & B
    \end{bmatrix}.
    \]

    The Schur complement after eliminating the first column is

    \[
    B-\frac{vw^{\mathsf T}}{\alpha}.
    \]

    Since only the first \(p\) entries of \(v\) and the first \(q\) entries of \(w\) can be nonzero, the rank-one correction does not enlarge the corresponding bandwidths. Applying the induction hypothesis to the Schur complement shows that \(L\) and \(U\) inherit the lower and upper bandwidths of the original matrix. \(\square\)

### 6.4 The Effect of Pivoting on Bandwidth

!!! success "Theorem 6.3 (Upper Bandwidth after Partial Pivoting)"

    Suppose \(A\) is nonsingular and

    \[
    \operatorname{ub}(A)=q,
    \qquad
    \operatorname{lb}(A)=p.
    \]

    After columnwise partial pivoting, the resulting upper triangular factor satisfies

    \[
    \operatorname{ub}(U)\leq p+q.
    \]

Row interchanges may move entries from the lower band above the main diagonal, so the upper bandwidth can grow from \(q\) to \(p+q\), but it remains much smaller than the bandwidth of a general dense matrix.

### 6.5 Total Computational Complexity

For a system that requires no pivoting and preserves its bandwidth,

\[
\text{LU factorization: }O(npq),
\]

\[
\text{solve }Ly=b\text{: }O(np),
\]

\[
\text{solve }Ux=y\text{: }O(nq).
\]

Therefore, the total complexity is

\[
O\bigl(n(pq+p+q)\bigr).
\]

In particular, for a tridiagonal matrix,

\[
p=q=1,
\]

so factorization and solution together require only

\[
O(n)
\]

operations.

!!! info "Benefits of Structure"

    When \(p,q\ll n\), banded LU, \(LDL^{\mathsf T}\), and Cholesky factorizations are much more efficient in both time and storage than their dense counterparts.

---

## 7. Chapter Summary

!!! success "Key Conclusions"

    1. Strict diagonal dominance guarantees invertibility and gives a direct upper bound for the norm of the inverse.

    2. A nonsingular column diagonally dominant matrix admits LU factorization without pivoting, with elimination multipliers bounded in absolute value by \(1\).

    3. The \(LDM^{\mathsf T}\) factorization of a symmetric matrix reduces to an \(LDL^{\mathsf T}\) factorization.

    4. A symmetric positive definite matrix has a unique Cholesky factorization \(A=GG^{\mathsf T}\), requiring roughly half the work of a general LU factorization.

    5. The polar decomposition \(A=QP\) separates a matrix into orthogonal and symmetric positive semidefinite factors and yields the best orthogonal approximation.

    6. Banded structure reduces factorization complexity from \(O(n^3)\) to \(O(npq)\); a tridiagonal system can even be solved in linear time.
