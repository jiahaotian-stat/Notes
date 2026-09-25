# Chapter 2: Matrix Norms, Perturbation Theory, and Singular Value Decomposition

This chapter first continues the discussion of relationships among common matrix norms and introduces the orthogonal invariance of the spectral norm and the Frobenius norm. Then it uses the Neumann series to study how inverse matrices change under perturbations. Finally, it systematically introduces the singular value decomposition (SVD) and discusses basic properties of singular values, perturbation estimates, and low-rank approximation.

---

## 1. Properties of Common Matrix Norms

Let $A=(a_{ij})\in\mathbb{R}^{m\times n}$. This chapter mainly uses the following four matrix norms:

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

and the induced $2$-norm (spectral norm)

\[
\|A\|_2
=\max_{x\neq0}\frac{\|Ax\|_2}{\|x\|_2}.
\]

### 1.1 Equivalence of Norms in Finite-Dimensional Spaces

!!! success "Theorem 1.1 (Equivalence of Norms in Finite-Dimensional Spaces)"

    On the finite-dimensional linear space $\mathbb{R}^{m\times n}$, any two norms are equivalent. That is, for any two matrix norms $\|\cdot\|_a$ and $\|\cdot\|_b$, there exist constants $c,C>0$ depending only on $m,n$, such that

    \[
    c\|A\|_a\leq\|A\|_b\leq C\|A\|_a,
    \qquad A\in\mathbb{R}^{m\times n}.
    \]

Therefore, in a finite-dimensional space, whether a matrix sequence converges does not depend on the specific norm chosen. However, different norms may give different numerical magnitudes and error bounds.

### 1.2 Spectral Norm and Frobenius Norm

!!! success "Theorem 1.2 (Relationship Between the Spectral Norm and the Frobenius Norm)"

    For any $A\in\mathbb{R}^{m\times n}$,

    \[
    \|A\|_2
    \leq\|A\|_F
    \leq\sqrt{\min\{m,n\}}\,\|A\|_2.
    \]

??? proof "Proof of Theorem 1.2 (click to expand)"

    Since $A^TA$ is positive semidefinite, let its nonzero eigenvalues be

    \[
    \lambda_1\geq\lambda_2\geq\cdots\geq\lambda_r>0,
    \qquad r=\operatorname{rank}(A).
    \]

    Note that

    \[
    \|A\|_2^2=\lambda_1,
    \qquad
    \|A\|_F^2=\operatorname{tr}(A^TA)=\sum_{i=1}^r\lambda_i.
    \]

    Therefore

    \[
    \lambda_1
    \leq\sum_{i=1}^r\lambda_i
    \leq r\lambda_1
    \leq\min\{m,n\}\lambda_1.
    \]

    Taking square roots of both sides gives the conclusion. $\square$

Another direct way to prove the upper bound

\[
\|A\|_2\leq\|A\|_F
\]

is to use the Cauchy-Schwarz inequality. For any $\|x\|_2=1$,

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

Taking the supremum over unit vectors gives the result.

### 1.3 Transpose and Common Norms

!!! success "Proposition 1.3 (Invariance Under Transposition and Duality Relations)"

    For any $A\in\mathbb{R}^{m\times n}$,

    \[
    \|A^T\|_F=\|A\|_F,
    \qquad
    \|A^T\|_2=\|A\|_2,
    \]

    and

    \[
    \|A\|_1=\|A^T\|_\infty,
    \qquad
    \|A\|_\infty=\|A^T\|_1.
    \]

The first equality follows directly from the definition of the Frobenius norm; the second equality follows from the fact that $A^TA$ and $AA^T$ have the same nonzero eigenvalues; the last two equalities follow from the interchange of row sums and column sums under transposition.

### 1.4 Spectral Norm and the 1-Norm, ∞-Norm

!!! success "Theorem 1.4 (Relationships Among Common Induced Norms)"

    For any $A\in\mathbb{R}^{m\times n}$,

    \[
    \frac{1}{\sqrt{n}}\|A\|_\infty
    \leq\|A\|_2
    \leq\sqrt{m}\,\|A\|_\infty,
    \]

    and

    \[
    \frac{1}{\sqrt{m}}\|A\|_1
    \leq\|A\|_2
    \leq\sqrt{n}\,\|A\|_1.
    \]

??? proof "Proof of Theorem 1.4 (click to expand)"

    For any $x\in\mathbb{R}^n$, from the relationships among vector norms,

    \[
    \|Ax\|_2
    \leq\sqrt{m}\|Ax\|_\infty
    \leq\sqrt{m}\|A\|_\infty\|x\|_\infty
    \leq\sqrt{m}\|A\|_\infty\|x\|_2.
    \]

    Hence

    \[
    \|A\|_2\leq\sqrt{m}\|A\|_\infty.
    \]

    On the other hand, let $e_i$ be the $i$-th standard basis vector. For any $i=1,\ldots,m$,

    \[
    \sum_{j=1}^na_{ij}^2
    =\|A^Te_i\|_2^2
    \leq\|A^T\|_2^2
    =\|A\|_2^2.
    \]

    Then, by the Cauchy-Schwarz inequality,

    \[
    \left(\sum_{j=1}^n|a_{ij}|\right)^2
    \leq n\sum_{j=1}^na_{ij}^2
    \leq n\|A\|_2^2.
    \]

    Taking the maximum over $i$ gives

    \[
    \|A\|_\infty\leq\sqrt{n}\|A\|_2.
    \]

    This proves the first set of inequalities. Applying the first set of inequalities to $A^T$ and using

    \[
    \|A^T\|_2=\|A\|_2,
    \qquad
    \|A^T\|_\infty=\|A\|_1,
    \]

    yields the second set of inequalities. $\square$

### 1.5 Eigenvalue Characterization of the Spectral Norm

!!! success "Theorem 1.5 (Achieving Vector for the Spectral Norm)"

    For any $A\in\mathbb{R}^{m\times n}$, there exists a unit vector $z\in\mathbb{R}^n$ such that

    \[
    A^TAz=\|A\|_2^2z.
    \]

    In other words, $\|A\|_2^2$ is the largest eigenvalue of $A^TA$.

??? proof "Proof of Theorem 1.5 (click to expand)"

    Define the Rayleigh quotient on the unit sphere:

    \[
    q(x)=x^TA^TAx=\|Ax\|_2^2,
    \qquad \|x\|_2=1.
    \]

    The unit sphere is compact and $q$ is continuous, so the maximum is attained. On the other hand, by the spectral theorem for real symmetric matrices, let the eigenvalues of $A^TA$ satisfy

    \[
    \lambda_1\geq\lambda_2\geq\cdots\geq\lambda_n\geq0,
    \]

    with corresponding orthonormal eigenvectors $v_1,\ldots,v_n$. If

    \[
    x=\sum_{i=1}^n\alpha_iv_i,
    \qquad
    \sum_{i=1}^n\alpha_i^2=1,
    \]

    then

    \[
    x^TA^TAx
    =\sum_{i=1}^n\lambda_i\alpha_i^2
    \leq\lambda_1.
    \]

    Equality holds when $x=v_1$, so

    \[
    \|A\|_2^2=\lambda_1,
    \]

    and one may take $z=v_1$. $\square$

!!! success "Corollary 1.6"

    For any $A\in\mathbb{R}^{m\times n}$,

    \[
    \|A\|_2
    \leq\sqrt{\|A\|_1\|A\|_\infty}.
    \]

??? proof "Proof of Corollary 1.6 (click to expand)"

    Take a nonzero vector $z$ such that

    \[
    A^TAz=\|A\|_2^2z.
    \]

    Taking the $1$-norm on both sides and using the compatibility of induced norms gives

    \[
    \|A\|_2^2\|z\|_1
    =\|A^TAz\|_1
    \leq\|A^T\|_1\|A\|_1\|z\|_1.
    \]

    Since $\|A^T\|_1=\|A\|_\infty$, canceling $\|z\|_1$ gives

    \[
    \|A\|_2^2\leq\|A\|_1\|A\|_\infty.
    \]

    Taking square roots of both sides gives the result. $\square$

### 1.6 Orthogonal Invariance

!!! success "Theorem 1.7 (Orthogonal Invariance)"

    Let $A\in\mathbb{R}^{m\times n}$, and let $Q\in\mathbb{R}^{m\times m}$ and $P\in\mathbb{R}^{n\times n}$ be orthogonal matrices. Then

    \[
    \|QAP^T\|_2=\|A\|_2,
    \qquad
    \|QAP^T\|_F=\|A\|_F.
    \]

??? proof "Proof of Theorem 1.7 (click to expand)"

    For the spectral norm, since orthogonal matrices preserve the Euclidean norm,

    \[
    \begin{aligned}
    \|QAP^T\|_2
    &=\max_{\|x\|_2=1}\|QAP^Tx\|_2\\
    &=\max_{\|x\|_2=1}\|AP^Tx\|_2\\
    &=\max_{\|y\|_2=1}\|Ay\|_2
    =\|A\|_2,
    \end{aligned}
    \]

    where we set $y=P^Tx$.

    For the Frobenius norm, using the cyclic invariance of the trace,

    \[
    \begin{aligned}
    \|QAP^T\|_F^2
    &=\operatorname{tr}\bigl(PA^TQ^TQAP^T\bigr)\\
    &=\operatorname{tr}\bigl(PA^TAP^T\bigr)\\
    &=\operatorname{tr}(A^TA)
    =\|A\|_F^2.
    \end{aligned}
    \]

    Hence the conclusion holds. $\square$

---

## 2. Matrix Perturbation and Inverse Matrices

### 2.1 Neumann Series

!!! success "Lemma 2.1 (Neumann Series)"

    Let $F\in\mathbb{R}^{n\times n}$, and suppose that for some compatible matrix norm $\|\cdot\|_p$,

    \[
    \|F\|_p<1.
    \]

    Then $I-F$ is nonsingular, and

    \[
    (I-F)^{-1}=\sum_{k=0}^{\infty}F^k.
    \]

    Moreover,

    \[
    \|(I-F)^{-1}\|_p
    \leq\frac{1}{1-\|F\|_p}.
    \]

??? proof "Proof of Lemma 2.1 (click to expand)"

    Let

    \[
    S_N=\sum_{k=0}^NF^k.
    \]

    When $M>N$, by the triangle inequality and compatibility of the norm,

    \[
    \begin{aligned}
    \|S_M-S_N\|_p
    &\leq\sum_{k=N+1}^M\|F^k\|_p\\
    &\leq\sum_{k=N+1}^M\|F\|_p^k\\
    &\leq\frac{\|F\|_p^{N+1}}{1-\|F\|_p}.
    \end{aligned}
    \]

    Since $\|F\|_p<1$, the right-hand side tends to $0$, so $\{S_N\}$ is a Cauchy sequence. A finite-dimensional normed space is complete, so there exists a matrix $S$ such that $S_N\to S$.

    Also, since

    \[
    S_N(I-F)=(I-F)S_N=I-F^{N+1},
    \]

    and $F^{N+1}\to0$, letting $N\to\infty$ gives

    \[
    S(I-F)=(I-F)S=I.
    \]

    Therefore $S=(I-F)^{-1}$. Finally,

    \[
    \|(I-F)^{-1}\|_p
    \leq\sum_{k=0}^{\infty}\|F\|_p^k
    =\frac{1}{1-\|F\|_p}.
    \]

    This completes the proof. $\square$

### 2.2 Perturbation Bounds for the Inverse Matrix

!!! success "Theorem 2.2 (Perturbation of the Inverse Matrix)"

    Let $A\in\mathbb{R}^{n\times n}$ be nonsingular, let $E\in\mathbb{R}^{n\times n}$, and let

    \[
    r=\|A^{-1}E\|_p<1.
    \]

    Then $A+E$ is nonsingular, and

    \[
    \|(A+E)^{-1}-A^{-1}\|_p
    \leq
    \frac{\|E\|_p\|A^{-1}\|_p^2}{1-r}.
    \]

??? proof "Proof of Theorem 2.2 (click to expand)"

    Let

    \[
    F=-A^{-1}E.
    \]

    Then $\|F\|_p=r<1$, and

    \[
    A+E=A(I-F).
    \]

    By the Neumann series lemma, $I-F$ is nonsingular, hence $A+E$ is nonsingular, and

    \[
    (A+E)^{-1}=(I-F)^{-1}A^{-1}.
    \]

    Using the identity for the difference of inverse matrices,

    \[
    (A+E)^{-1}-A^{-1}
    =-A^{-1}E(A+E)^{-1}.
    \]

    Also,

    \[
    \|(A+E)^{-1}\|_p
    \leq\|(I-F)^{-1}\|_p\|A^{-1}\|_p
    \leq\frac{\|A^{-1}\|_p}{1-r}.
    \]

    Therefore

    \[
    \begin{aligned}
    \|(A+E)^{-1}-A^{-1}\|_p
    &\leq\|A^{-1}\|_p\|E\|_p\|(A+E)^{-1}\|_p\\
    &\leq\frac{\|E\|_p\|A^{-1}\|_p^2}{1-r}.
    \end{aligned}
    \]

    This completes the proof. $\square$

---

## 3. Singular Value Decomposition

### 3.1 Existence of the SVD

!!! success "Theorem 3.1 (Singular Value Decomposition of a Real Matrix)"

    For any matrix $A\in\mathbb{R}^{m\times n}$, there exist orthogonal matrices

    \[
    U\in\mathbb{R}^{m\times m},
    \qquad
    V\in\mathbb{R}^{n\times n},
    \]

    and a rectangular diagonal matrix $\Sigma\in\mathbb{R}^{m\times n}$, such that

    \[
    A=U\Sigma V^T.
    \]

    The diagonal entries of $\Sigma$ satisfy

    \[
    \sigma_1\geq\sigma_2\geq\cdots\geq\sigma_p\geq0,
    \qquad
    p=\min\{m,n\}.
    \]

    The $\sigma_i$ are called the singular values of $A$; the columns $u_i$ of $U$ are called left singular vectors; the columns $v_i$ of $V$ are called right singular vectors.

??? proof "Proof of Theorem 3.1 (click to expand)"

    The following gives the inductive proof in the lecture notes. It suffices to prove the case $m\leq n$; when $m>n$, one can apply the result to $A^T$ and then transpose.

    When $A=0$ the conclusion is obvious. In what follows assume $A\neq0$, and perform induction on $m$.

    When $m=1$, take

    \[
    v_1=\frac{A^T}{\|A\|_2},
    \]

    and extend $v_1$ to an orthonormal basis of $\mathbb{R}^n$. Form an orthogonal matrix $V$ with these vectors as columns, and take $U=[1]$; then

    \[
    U^TAV=(\|A\|_2,0,\ldots,0).
    \]

    Now suppose the conclusion holds for matrices with fewer than $m$ rows. By Theorem 1.5, there exists a unit vector $x\in\mathbb{R}^n$ such that

    \[
    A^TAx=\|A\|_2^2x.
    \]

    Let

    \[
    \sigma_1=\|A\|_2,
    \qquad
    y=\frac{Ax}{\sigma_1}.
    \]

    Then $\|y\|_2=1$ and $Ax=\sigma_1y$. Extend $y$ and $x$ to orthonormal bases of $\mathbb{R}^m$ and $\mathbb{R}^n$, respectively, and let

    \[
    U_0=(y,u_2,\ldots,u_m),
    \qquad
    V_0=(x,v_2,\ldots,v_n).
    \]

    Then

    \[
    U_0^TAV_0=
    \begin{pmatrix}
    \sigma_1&w^T\\
    0&B
    \end{pmatrix}.
    \]

    Since orthogonal transformations preserve the spectral norm, the spectral norm of the left-hand side equals $\sigma_1$. If $w\neq0$, let

    \[
    s=\begin{pmatrix}\sigma_1\\w\end{pmatrix},
    \]

    then

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

    which contradicts the fact that the spectral norm of this matrix equals $\sigma_1$. Therefore $w=0$. Hence

    \[
    U_0^TAV_0=
    \begin{pmatrix}
    \sigma_1&0\\
    0&B
    \end{pmatrix}.
    \]

    Apply the induction hypothesis to $B\in\mathbb{R}^{(m-1)\times(n-1)}$, then embed the resulting orthogonal matrices into block diagonal matrices to obtain the singular value decomposition of $A$. $\square$

!!! note "Non-uniqueness of the SVD"

    The singular values themselves are uniquely determined, but the singular vectors are not necessarily unique. For example, if a singular value has multiplicity greater than $1$, one can perform any orthogonal transformation in the corresponding singular subspace; even if the singular values are distinct, one can simultaneously change the signs of a pair of left and right singular vectors.

### 3.2 Relationship Between Left and Right Singular Vectors

From

\[
AV=U\Sigma,
\qquad
A^TU=V\Sigma^T,
\]

we immediately obtain the following relationships.

!!! success "Corollary 3.2 (Singular Vector Equations)"

    For $i=1,\ldots,p$,

    \[
    Av_i=\sigma_i u_i,
    \qquad
    A^Tu_i=\sigma_i v_i.
    \]

Furthermore,

\[
A^TAv_i=\sigma_i^2v_i,
\qquad
AA^Tu_i=\sigma_i^2u_i.
\]

Therefore, the squares of the nonzero singular values of $A$ are exactly the nonzero eigenvalues of $A^TA$ (and equivalently of $AA^T$).

### 3.3 Norms, Rank, and Fundamental Subspaces

!!! success "Corollary 3.3 (Norms and Singular Values)"

    For any $A\in\mathbb{R}^{m\times n}$,

    \[
    \|A\|_2=\sigma_1,
    \qquad
    \|A\|_F=\left(\sum_{i=1}^p\sigma_i^2\right)^{1/2}.
    \]

Let $r=\operatorname{rank}(A)$. Then

\[
\sigma_1\geq\cdots\geq\sigma_r>0,
\qquad
\sigma_{r+1}=\cdots=\sigma_p=0,
\]

and the SVD can be written as a sum of rank-one matrices:

\[
A=\sum_{i=1}^r\sigma_i u_iv_i^T.
\]

!!! success "Corollary 3.4 (Four Fundamental Subspaces)"

    If $\operatorname{rank}(A)=r$, then

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

### 3.4 Characterization of the Smallest Singular Value

!!! success "Lemma 3.5 (Smallest Singular Value)"

    Let $A\in\mathbb{R}^{m\times n}$, where $m\geq n$, and suppose $A$ has full column rank. Then for any $x\in\mathbb{R}^n$,

    \[
    \|Ax\|_2\geq\sigma_n(A)\|x\|_2.
    \]

    Moreover, there exists a nonzero vector $x_0$ for which equality holds. Therefore

    \[
    \sigma_n(A)
    =\min_{\|x\|_2=1}\|Ax\|_2.
    \]

??? proof "Proof of Lemma 3.5 (click to expand)"

    Let $A=U\Sigma V^T$, and let $y=V^Tx$. Since $V$ is orthogonal, $\|y\|_2=\|x\|_2$. Therefore

    \[
    \|Ax\|_2^2
    =\|\Sigma y\|_2^2
    =\sum_{i=1}^n\sigma_i^2y_i^2
    \geq\sigma_n^2\sum_{i=1}^ny_i^2
    =\sigma_n^2\|x\|_2^2.
    \]

    Equality holds when $x=v_n$. $\square$

---

## 4. Perturbation Properties of Singular Values

### 4.1 Perturbation Bounds for the Largest and Smallest Singular Values

!!! success "Theorem 4.1 (Basic Bounds for Perturbation of Singular Values)"

    Let $A,E\in\mathbb{R}^{m\times n}$. Then

    \[
    \sigma_{\max}(A+E)
    \leq\sigma_{\max}(A)+\|E\|_2.
    \]

    If in addition $m\geq n$, and $A+E$ has full column rank, then

    \[
    \sigma_{\min}(A+E)
    \geq\sigma_{\min}(A)-\|E\|_2.
    \]

??? proof "Proof of Theorem 4.1 (click to expand)"

    The first inequality follows directly from the triangle inequality for the spectral norm:

    \[
    \sigma_{\max}(A+E)
    =\|A+E\|_2
    \leq\|A\|_2+\|E\|_2.
    \]

    For the second inequality, take a unit vector $x_0$ such that

    \[
    \|(A+E)x_0\|_2=\sigma_{\min}(A+E).
    \]

    Then

    \[
    \begin{aligned}
    \sigma_{\min}(A)
    &\leq\|Ax_0\|_2\\
    &\leq\|(A+E)x_0\|_2+\|Ex_0\|_2\\
    &\leq\sigma_{\min}(A+E)+\|E\|_2.
    \end{aligned}
    \]

    Rearranging yields the conclusion. $\square$

!!! example "Corollary 4.2 (Singular Values After Adding a Column)"

    Let $A\in\mathbb{R}^{m\times n}$, $m>n$, and let $z\in\mathbb{R}^m$. Denote by $[A\ z]$ the matrix obtained by appending the column vector $z$ to the right of $A$. Then

    \[
    \sigma_{\max}([A\ z])\geq\sigma_{\max}(A),
    \]

    \[
    \sigma_{\min}([A\ z])\leq\sigma_{\min}(A).
    \]

This is because $[A\ z](x,0)^T=Ax$, while the largest singular value is the maximum over the unit sphere and the smallest singular value is the minimum over the unit sphere. Adding a variable direction enlarges the feasible set of the maximization problem and also enlarges the feasible set of the minimization problem.

---

## 5. Best Low-Rank Approximation

Let $A\in\mathbb{R}^{m\times n}$ have rank $r$, and let its compact singular value decomposition be

\[
A=\sum_{i=1}^r\sigma_i u_iv_i^T,
\qquad
\sigma_1\geq\cdots\geq\sigma_r>0.
\]

For $k<r$, define the truncated SVD

\[
A_k=\sum_{i=1}^k\sigma_i u_iv_i^T.
\]

Clearly $\operatorname{rank}(A_k)=k$.

### 5.1 Eckart-Young Theorem in the Spectral Norm

!!! success "Theorem 5.1 (Eckart-Young Theorem)"

    If $k<r=\operatorname{rank}(A)$, then

    \[
    \min_{\operatorname{rank}(B)\leq k}\|A-B\|_2
    =\|A-A_k\|_2
    =\sigma_{k+1}.
    \]

    Since $A_k$ has rank exactly $k$, the same minimum is obtained when the constraint is written as $\operatorname{rank}(B)=k$.

??? proof "Proof of Theorem 5.1 (click to expand)"

    First, by orthogonal invariance,

    \[
    \begin{aligned}
    \|A-A_k\|_2
    &=\left\|U^T(A-A_k)V\right\|_2\\
    &=\left\|\operatorname{diag}
    (0,\ldots,0,\sigma_{k+1},\ldots,\sigma_r)\right\|_2\\
    &=\sigma_{k+1}.
    \end{aligned}
    \]

    Next we prove that no matrix of rank at most $k$ can yield a smaller error. Let $\operatorname{rank}(B)\leq k$. Then

    \[
    \dim\operatorname{Null}(B)\geq n-k.
    \]

    Let

    \[
    S=\operatorname{span}\{v_1,\ldots,v_{k+1}\}.
    \]

    Since $\dim S=k+1$, by the dimension formula,

    \[
    \dim\bigl(\operatorname{Null}(B)\cap S\bigr)
    \geq(n-k)+(k+1)-n=1.
    \]

    Therefore we can choose a unit vector $z\in\operatorname{Null}(B)\cap S$. Write it as

    \[
    z=\sum_{i=1}^{k+1}(v_i^Tz)v_i.
    \]

    Since $Bz=0$,

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

    Hence $\|A-B\|_2\geq\sigma_{k+1}$. Combining this with the fact that $A_k$ attains this lower bound proves the theorem. $\square$

### 5.2 Openness and Density of the Set of Full-Rank Matrices

The Eckart-Young theorem shows that the spectral norm distance from a matrix to the set of low-rank matrices is characterized by the corresponding singular values.

!!! success "Proposition 5.2 (The Set of Full-Rank Matrices Is Open and Dense)"

    Let $p=\min\{m,n\}$. The set of all $m\times n$ real matrices of rank $p$ is both open and dense in $\mathbb{R}^{m\times n}$.

??? proof "Explanation of Proposition 5.2 (click to expand)"

    If $A$ has full rank, then $\sigma_p(A)>0$. When $m\geq n$, apply the singular value perturbation bound directly; when $m<n$, apply the same conclusion to the transpose. Thus, as long as

    \[
    \|E\|_2<\sigma_p(A),
    \]

    we have

    \[
    \sigma_p(A+E)
    \geq\sigma_p(A)-\|E\|_2>0.
    \]

    Hence $A+E$ still has full rank, so the set of full-rank matrices is open.

    On the other hand, suppose $A=U\Sigma V^T$ does not have full rank. For any $\varepsilon>0$, replace all zero diagonal entries of $\Sigma$ with $\varepsilon$ to obtain $\Sigma_\varepsilon$, and let

    \[
    A_\varepsilon=U\Sigma_\varepsilon V^T.
    \]

    Then $A_\varepsilon$ has full rank, and

    \[
    \|A_\varepsilon-A\|_2=\varepsilon.
    \]

    Therefore every matrix can be approximated arbitrarily well by full-rank matrices, that is, the set of full-rank matrices is dense. $\square$

### 5.3 Sum Inequality for Singular Values

!!! success "Lemma 5.3"

    Let $X,Y\in\mathbb{R}^{m\times n}$, and suppose the indices satisfy $i+j-1\leq p=\min\{m,n\}$. Then

    \[
    \sigma_i(X)+\sigma_j(Y)
    \geq\sigma_{i+j-1}(X+Y).
    \]

??? proof "Proof of Lemma 5.3 (click to expand)"

    Let $X_{i-1}$ and $Y_{j-1}$ be the truncated SVDs of $X$ and $Y$, respectively. By the Eckart-Young theorem,

    \[
    \sigma_i(X)=\|X-X_{i-1}\|_2,
    \qquad
    \sigma_j(Y)=\|Y-Y_{j-1}\|_2.
    \]

    Thus

    \[
    \begin{aligned}
    \sigma_i(X)+\sigma_j(Y)
    &\geq\|X+Y-(X_{i-1}+Y_{j-1})\|_2.
    \end{aligned}
    \]

    Also, since

    \[
    \operatorname{rank}(X_{i-1}+Y_{j-1})
    \leq(i-1)+(j-1)=i+j-2,
    \]

    applying the Eckart-Young theorem again gives

    \[
    \|X+Y-(X_{i-1}+Y_{j-1})\|_2
    \geq\sigma_{i+j-1}(X+Y).
    \]

    Combining the two inequalities yields the conclusion. $\square$

### 5.4 Best Low-Rank Approximation in the Frobenius Norm

!!! success "Theorem 5.4 (Low-Rank Approximation in the Frobenius Norm)"

    If $k<r=\operatorname{rank}(A)$, then

    \[
    \min_{\operatorname{rank}(B)\leq k}\|A-B\|_F
    =\|A-A_k\|_F
    =\left(\sum_{i=k+1}^r\sigma_i^2\right)^{1/2}.
    \]

??? proof "Proof of Theorem 5.4 (click to expand)"

    First, by the orthogonal invariance of the Frobenius norm,

    \[
    \|A-A_k\|_F^2
    =\sum_{i=k+1}^r\sigma_i^2.
    \]

    Now take any matrix $B$ satisfying $\operatorname{rank}(B)\leq k$. Since

    \[
    \sigma_{k+1}(B)=0,
    \]

    taking $X=A-B$, $Y=B$, and $j=k+1$ in Lemma 5.3 gives

    \[
    \sigma_i(A-B)
    \geq\sigma_{i+k}(A),
    \qquad i=1,\ldots,r-k.
    \]

    Therefore

    \[
    \begin{aligned}
    \|A-B\|_F^2
    &=\sum_{i=1}^p\sigma_i^2(A-B)\\
    &\geq\sum_{i=1}^{r-k}\sigma_i^2(A-B)\\
    &\geq\sum_{i=1}^{r-k}\sigma_{i+k}^2(A)
    =\sum_{i=k+1}^r\sigma_i^2(A).
    \end{aligned}
    \]

    The truncated SVD $A_k$ attains this lower bound, so it is the best rank $k$ approximation in the Frobenius norm. $\square$

---

## 6. Summary of This Chapter

The core conclusions of this chapter are as follows:

1. In a finite-dimensional matrix space, all norms are equivalent; there are explicit comparison inequalities among common matrix norms.
2. The spectral norm and Frobenius norm are invariant under left and right orthogonal transformations.
3. When $\|F\|<1$, $(I-F)^{-1}$ can be represented as a Neumann series; this conclusion gives perturbation bounds for inverse matrices.
4. Every real matrix has a singular value decomposition $A=U\Sigma V^T$; the singular values simultaneously characterize the matrix's norm, rank, and four fundamental subspaces.
5. The truncated SVD is the best low-rank approximation in both the spectral norm and the Frobenius norm.