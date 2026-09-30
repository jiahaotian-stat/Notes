# Chapter 3: Gaussian Elimination and LU Factorization

This chapter studies direct methods for general linear systems. We begin with forward and backward substitution for triangular systems, describe elimination through Gauss transformations, and then derive the LU factorization. Finally, we discuss rounding errors, partial pivoting, and LU factorization with pivoting.

---

## 1. Triangular Systems

Let

\[
A\in\mathbb{R}^{n\times n},
\qquad
b\in\mathbb{R}^n.
\]

We wish to solve the linear system

\[
Ax=b.
\]

When $A$ is triangular, the unknowns can be determined successively. This is the foundation of general elimination methods.

### 1.1 Lower Triangular Systems and Forward Substitution

Let $L=(l_{ij})$ be a nonsingular lower triangular matrix. The system $Lx=b$ can be written as

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

Solving row by row from the first row gives the forward-substitution formulas

\[
x_1=\frac{b_1}{l_{11}},
\qquad
x_i=\frac{1}{l_{ii}}
\left(b_i-\sum_{j=1}^{i-1}l_{ij}x_j\right),
\quad i=2,\ldots,n.
\]

!!! note "Algorithm 1.1 (Forward Substitution)"

    For $i=1,\ldots,n$, compute successively

    \[
    x_i=\frac{b_i-l_{i,1:i-1}x_{1:i-1}}{l_{ii}}.
    \]

    If $L$ is unit lower triangular, then $l_{ii}=1$, so no division is needed.

### 1.2 Upper Triangular Systems and Backward Substitution

Let $U=(u_{ij})$ be a nonsingular upper triangular matrix. Solving upward from the last row gives

\[
x_n=\frac{b_n}{u_{nn}},
\qquad
x_i=\frac{1}{u_{ii}}
\left(b_i-\sum_{j=i+1}^{n}u_{ij}x_j\right),
\quad i=n-1,\ldots,1.
\]

!!! note "Algorithm 1.2 (Backward Substitution)"

    For $i=n,n-1,\ldots,1$, compute successively

    \[
    x_i=\frac{b_i-u_{i,i+1:n}x_{i+1:n}}{u_{ii}}.
    \]

Both forward and backward substitution require approximately $n^2$ floating-point operations and therefore have complexity $O(n^2)$.

### 1.3 Algebraic Properties of Triangular Matrices

!!! success "Proposition 1.3 (Closure of Triangular Matrices)"

    The product and inverse of upper triangular matrices are upper triangular; the product and inverse of lower triangular matrices are lower triangular.

    A triangular matrix whose diagonal entries are all $1$ is called a unit triangular matrix. The product and inverse of unit triangular matrices remain unit triangular matrices of the same type.

??? proof "Proof of Proposition 1.3 (click to expand)"

    Consider upper triangular matrices. If both $A$ and $B$ are upper triangular, then for $i>j$,

    \[
    (AB)_{ij}=\sum_{k=1}^n a_{ik}b_{kj}=0,
    \]

    because every term has either $a_{ik}=0$ or $b_{kj}=0$. Hence $AB$ is upper triangular.

    If $A$ is nonsingular, each column of $AX=I$ can be found by backward substitution, and the resulting matrix $X=A^{-1}$ is upper triangular. The lower triangular case is analogous. If every diagonal entry equals $1$, the diagonal entries of both the product and inverse also equal $1$. $\square$

### 1.4 Backward Error of Triangular Solves

Let the unit roundoff be $u$. The floating-point solution $\widehat{x}$ produced by forward substitution can be interpreted as the exact solution of a nearby lower triangular system:

\[
(L+\Delta L)\widehat{x}=b,
\qquad
|\Delta L|\leq \gamma_n|L|,
\]

where

\[
\gamma_n=\frac{nu}{1-nu}=nu+O(u^2).
\]

Similarly, backward substitution satisfies

\[
(U+\Delta U)\widehat{x}=b,
\qquad
|\Delta U|\leq \gamma_n|U|.
\]

!!! info "Interpretation"

    Thus triangular solves are backward stable: the computed result is the exact solution of a triangular system that is very close to the original one.

---

## 2. Gaussian Elimination

### 2.1 Elementary Elimination and Gauss Transformations

Let $e_i$ denote the $i$th standard basis vector. Left multiplication by

\[
I+\lambda e_j e_i^{\mathsf T}
\]

is equivalent to adding $\lambda$ times row $i$ to row $j$.

To eliminate the first column, suppose $a_{11}\neq0$ and define the multipliers

\[
\tau_i^{(1)}=\frac{a_{i1}}{a_{11}},
\qquad i=2,\ldots,n,
\]

and the vector

\[
\tau^{(1)}=
\begin{bmatrix}
0 & \tau_2^{(1)} & \cdots & \tau_n^{(1)}
\end{bmatrix}^{\mathsf T}.
\]

The first Gauss transformation is then

\[
M_1=I-\tau^{(1)}e_1^{\mathsf T},
\qquad
A^{(2)}=M_1A^{(1)},
\qquad
A^{(1)}=A.
\]

It annihilates every entry below the pivot in the first column.

### 2.2 The $k$th Elimination Step

Assume that the first $k-1$ columns have already been eliminated. At step $k$, the pivot is $a_{kk}^{(k)}$. If

\[
a_{kk}^{(k)}\neq0,
\]

define

\[
\tau_i^{(k)}=
\frac{a_{ik}^{(k)}}{a_{kk}^{(k)}},
\quad i=k+1,\ldots,n,
\]

and set $\tau_1^{(k)},\ldots,\tau_k^{(k)}$ equal to zero. The corresponding Gauss transformation is

\[
M_k=I-\tau^{(k)}e_k^{\mathsf T},
\qquad
A^{(k+1)}=M_kA^{(k)}.
\]

After $n-1$ steps, we obtain the upper triangular matrix

\[
M_{n-1}M_{n-2}\cdots M_1A=A^{(n)}=U.
\]

!!! warning "Pivot Condition"

    Gaussian elimination without row interchanges requires every pivot $a_{kk}^{(k)}$ to be nonzero. Nonsingularity of the matrix alone does not guarantee this condition.

### 2.3 Properties of Gauss Transformations

Since

\[
e_k^{\mathsf T}\tau^{(k)}=0,
\]

we have

\[
M_k^{-1}=I+\tau^{(k)}e_k^{\mathsf T}.
\]

Indeed,

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

When $i<j$, we also have

\[
M_iM_j
=I-\tau^{(i)}e_i^{\mathsf T}
-\tau^{(j)}e_j^{\mathsf T},
\]

because $e_i^{\mathsf T}\tau^{(j)}=0$ in the cross term.

### 2.4 Elimination Updates and Operation Count

At step $k$, for $i=k+1,\ldots,n$, first compute the multiplier

\[
l_{ik}=\frac{a_{ik}^{(k)}}{a_{kk}^{(k)}},
\]

and then update the remaining entries by

\[
a_{ij}^{(k+1)}
=a_{ij}^{(k)}-l_{ik}a_{kj}^{(k)},
\qquad
i,j=k+1,\ldots,n.
\]

The leading operation count of Gaussian elimination is approximately

\[
\frac{2}{3}n^3+O(n^2)
\]

floating-point operations, so its complexity is $O(n^3)$.

---

## 3. LU Factorization

### 3.1 From Gaussian Elimination to LU Factorization

From

\[
M_{n-1}\cdots M_1A=U,
\]

we obtain

\[
A=M_1^{-1}M_2^{-1}\cdots M_{n-1}^{-1}U.
\]

Define

\[
L=M_1^{-1}M_2^{-1}\cdots M_{n-1}^{-1}.
\]

Then $L$ is unit lower triangular, and hence

\[
A=LU.
\]

Moreover, the structure of the Gauss transformations gives

\[
L=I+\sum_{k=1}^{n-1}\tau^{(k)}e_k^{\mathsf T}.
\]

Thus the strict lower triangular part of $L$ stores exactly the elimination multipliers, while $U$ stores the upper triangular result of elimination.

### 3.2 Existence and Uniqueness of the LU Factorization

Let $A_k=A(1:k,1:k)$ denote the $k$th leading principal submatrix of $A$.

!!! success "Theorem 3.1 (LU Factorization)"

    If

    \[
    \det(A_k)\neq0,
    \qquad k=1,\ldots,n,
    \]

    then there exist unique matrices $L$ and $U$, with $L$ unit lower triangular and $U$ upper triangular, such that

    \[
    A=LU.
    \]

    Moreover,

    \[
    \det(A)=\prod_{i=1}^n u_{ii}.
    \]

??? proof "Proof of Theorem 3.1 (click to expand)"

    We use induction on the order of the matrix. The result is immediate for $n=1$. Assume it holds for matrices of order $n-1$, and partition $A$ as

    \[
    A=
    \begin{bmatrix}
    A_1 & v\\
    w^{\mathsf T} & \alpha
    \end{bmatrix}.
    \]

    By assumption, all leading principal minors of $A_1$ are nonzero, so there is a unique factorization

    \[
    A_1=L_1U_1.
    \]

    Therefore,

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

    This proves existence.

    If also

    \[
    A=L_1U_1=L_2U_2,
    \]

    then

    \[
    L_2^{-1}L_1=U_2U_1^{-1}.
    \]

    The left-hand side is unit lower triangular, whereas the right-hand side is upper triangular. They can be equal only if both are the identity. Therefore $L_1=L_2$ and $U_1=U_2$, proving uniqueness.

    Finally, since $\det(L)=1$,

    \[
    \det(A)=\det(L)\det(U)=\prod_{i=1}^n u_{ii}.
    \]

    $\square$

!!! note "Remark on the Condition"

    Nonzero leading principal minors are a standard sufficient condition for the existence of an LU factorization without pivoting. Equivalently, this factorization can be completed whenever every pivot generated during elimination is nonzero.

### 3.3 Compact Storage and the Computational Algorithm

Because every diagonal entry of $L$ equals $1$, the strict lower triangular part of $L$ and the upper triangular part of $U$ can be stored together in the space originally occupied by $A$:

\[
\begin{bmatrix}
u_{11} & u_{12} & \cdots & u_{1n}\\
l_{21} & u_{22} & \cdots & u_{2n}\\
\vdots & \vdots & \ddots & \vdots\\
l_{n1} & l_{n2} & \cdots & u_{nn}
\end{bmatrix}.
\]

!!! note "Algorithm 3.2 (In-Place LU Factorization)"

    For $k=1,\ldots,n-1$:

    First compute

    \[
    a_{k+1:n,k}\leftarrow
    \frac{a_{k+1:n,k}}{a_{kk}},
    \]

    and then perform the rank-one update

    \[
    A_{k+1:n,k+1:n}\leftarrow
    A_{k+1:n,k+1:n}
    -a_{k+1:n,k}a_{k,k+1:n}.
    \]

    At termination, the strict lower triangular part stores the multipliers of $L$, and the upper triangular part stores $U$.

### 3.4 Solving Linear Systems with LU

Once $A=LU$ has been computed, the original system is equivalent to

\[
LUx=b.
\]

Introducing the intermediate variable $y=Ux$, first solve by forward substitution

\[
Ly=b,
\]

and then solve by backward substitution

\[
Ux=y.
\]

!!! info "Multiple Right-Hand Sides"

    The LU factorization is computed only once at a cost of $O(n^3)$. Each additional right-hand side then requires only two triangular solves, with cost $O(n^2)$.

### 3.5 The $LDM^{\mathsf T}$ Factorization

If $A=LU$ and every diagonal entry of $U$ is nonzero, define

\[
D=\operatorname{diag}(u_{11},\ldots,u_{nn}),
\qquad
M^{\mathsf T}=D^{-1}U.
\]

Then $M$ is unit lower triangular and

\[
A=LDM^{\mathsf T}.
\]

!!! success "Corollary 3.3 (The $LDM^{\mathsf T}$ Factorization)"

    Under the assumptions of Theorem 3.1, $A$ has a unique factorization

    \[
    A=LDM^{\mathsf T},
    \]

    where $L$ and $M$ are unit lower triangular and $D$ is nonsingular and diagonal.

---

## 4. Rounding Errors and Pivoting

### 4.1 Backward Error of LU Factorization

Let $\widehat{L}$ and $\widehat{U}$ be the computed floating-point factors. Under the standard rounding model, we can write

\[
\widehat{L}\widehat{U}=A+H,
\]

with the componentwise error bound

\[
|H|\leq
2(n-1)u
\left(|A|+|\widehat{L}||\widehat{U}|\right)
+O(u^2).
\]

The solution $\widehat{x}$ obtained from $\widehat{L}$, $\widehat{U}$, and two triangular solves satisfies

\[
(A+E)\widehat{x}=b,
\]

where

\[
|E|\leq
nu\left(2|A|+4|\widehat{L}||\widehat{U}|\right)
+O(u^2).
\]

!!! warning "The Danger of a Small Pivot"

    The error bound depends on $|\widehat{L}||\widehat{U}|$. If a pivot is very small, the elimination multipliers may become large, causing substantial element growth in $L$ and $U$ and amplifying rounding errors.

For example, for

\[
A=
\begin{bmatrix}
0.001 & 1\\
1 & 0
\end{bmatrix},
\]

elimination without a row interchange gives

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

Thus a nonzero but very small pivot can still cause numerical instability.

### 4.2 Interchange Matrices

The interchange matrix that swaps rows $i$ and $j$ can be written as

\[
P_{ij}
=I-(e_i-e_j)(e_i-e_j)^{\mathsf T}.
\]

It satisfies

\[
P_{ij}^{\mathsf T}=P_{ij}^{-1}=P_{ij},
\qquad
P_{ij}^2=I.
\]

Left multiplication by $P_{ij}$ is equivalent to interchanging rows $i$ and $j$ of a matrix.

### 4.3 Partial Pivoting

At step $k$, choose the entry with the largest absolute value among rows $k$ through $n$ of the current $k$th column:

\[
p_k=\operatorname*{arg\,max}_{k\leq i\leq n}
|a_{ik}^{(k)}|.
\]

Then interchange rows $k$ and $p_k$ before performing elimination. Consequently, every multiplier satisfies

\[
|l_{ik}|\leq1,
\qquad i>k.
\]

!!! note "Algorithm 4.1 (Gaussian Elimination with Partial Pivoting)"

    For $k=1,\ldots,n-1$:

    Choose

    \[
    p_k=\operatorname*{arg\,max}_{k\leq i\leq n}|a_{ik}|;
    \]

    interchange rows $k$ and $p_k$;

    compute the multipliers and update the Schur complement:

    \[
    l_{ik}=\frac{a_{ik}}{a_{kk}},
    \qquad
    a_{ij}\leftarrow a_{ij}-l_{ik}a_{kj}.
    \]

The product of all row interchanges forms a permutation matrix $P$, and the final factorization is

\[
PA=LU.
\]

To solve $Ax=b$, solve successively

\[
Ly=Pb,
\qquad
Ux=y.
\]

### 4.4 Complexity and Stability

Partial pivoting adds approximately $O(n^2)$ comparisons, while the main elimination work remains $O(n^3)$. Since $|l_{ik}|\leq1$, we have

\[
\|L\|_{\infty}\leq n.
\]

This controls the growth of the elimination multipliers. Although $U$ can still exhibit large element growth for specially constructed extreme examples, the method is usually stable in practice and is therefore the standard method for solving general dense linear systems.

!!! info "Complete Pivoting"

    Complete pivoting selects the entry of largest absolute value in the entire remaining submatrix and interchanges both the corresponding rows and columns. It provides stronger growth control, but the search requires $O(n^3)$ comparisons, so partial pivoting is usually preferred in practice.

---

## 5. Chapter Summary

!!! success "Key Conclusions"

    1. A triangular system can be solved by forward or backward substitution in $O(n^2)$ time.

    2. Gauss transformations express elimination as matrix multiplication, and the vectors in their inverses record exactly the elimination multipliers.

    3. If all elimination pivots are nonzero, then $A=LU$; with the unit lower triangular normalization, this factorization is unique.

    4. LU factorization costs $O(n^3)$, but each new right-hand side requires only $O(n^2)$ work in triangular solves.

    5. To avoid zero or excessively small pivots, practical computations usually use partial pivoting and obtain $PA=LU$.
