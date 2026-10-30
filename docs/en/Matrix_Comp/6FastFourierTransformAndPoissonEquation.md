# Chapter 6: Fast Fourier Transform and the Poisson Equation

This chapter begins with the discrete Fourier transform, derives the radix-2 fast Fourier transform and butterfly operation, and then discusses discrete sine and cosine transforms, diagonalization of circulant matrices, and fast solution of finite-difference discretizations of the Poisson equation.

---

## 1. Discrete Fourier Transform

### 1.1 From the Continuous Transform to the Discrete Transform

One convention for the continuous Fourier transform is

\[
\widehat{f}(\xi)=\frac{1}{2\pi}\int_{-\infty}^{\infty}f(x)e^{-\mathrm{i}\xi x}\,\mathrm{d}x.
\]

For the discrete case, let

\[
\omega_n=e^{-2\pi\mathrm{i}/n}.
\]

!!! info "Definition 1.1 (DFT Matrix)"

    The discrete Fourier transform matrix of order $n$ is

    \[
    F_n=
    \left[
    \omega_n^{jk}
    \right]_{j,k=0}^{n-1}.
    \]

    For $x=(x_0,\ldots,x_{n-1})^{\mathsf T}\in\mathbb{C}^n$, its discrete Fourier transform is

    \[
    y=F_nx,
    \]

    that is,

    \[
    y_k=\sum_{j=0}^{n-1}\omega_n^{kj}x_j,
    \qquad k=0,\ldots,n-1.
    \]

Direct matrix-vector multiplication costs $O(n^2)$ operations, whereas the FFT exploits the recursive structure of $F_n$ to reduce the cost to $O(n\log n)$.

### 1.2 Orthogonality and the Inverse Transform

If $p\neq q$, then

\[
\sum_{k=0}^{n-1}\omega_n^{k(q-p)}
=\frac{1-\omega_n^{n(q-p)}}{1-\omega_n^{q-p}}
=0.
\]

When $p=q$, the sum equals $n$. Hence

\[
F_n^{\mathsf H}F_n=nI.
\]

!!! success "Theorem 1.2 (Inverse of the DFT Matrix)"

    \[
    F_n^{-1}=\frac{1}{n}F_n^{\mathsf H}
    =\frac{1}{n}\overline{F_n}.
    \]

    Therefore the normalized matrix

    \[
    \frac{1}{\sqrt n}F_n
    \]

    is unitary.

The inverse transform can be computed with the same FFT structure in $O(n\log n)$ time.

---

## 2. Radix-2 Fast Fourier Transform

### 2.1 Even-Odd Splitting

Let

\[
n=2m.
\]

Split the input into its even-indexed and odd-indexed parts:

\[
x^{(e)}=(x_0,x_2,\ldots,x_{2m-2})^{\mathsf T},
\]

\[
x^{(o)}=(x_1,x_3,\ldots,x_{2m-1})^{\mathsf T}.
\]

Since

\[
\omega_{2m}^2=\omega_m,
\]

for $k=0,\ldots,m-1$,

\[
y_k
=\sum_{r=0}^{m-1}\omega_m^{kr}x_{2r}
+\omega_{2m}^k\sum_{r=0}^{m-1}\omega_m^{kr}x_{2r+1}.
\]

Define

\[
y^{(e)}=F_mx^{(e)},
\qquad
y^{(o)}=F_mx^{(o)},
\]

Then

\[
y_k=y_k^{(e)}+\omega_{2m}^ky_k^{(o)},
\]

\[
y_{k+m}=y_k^{(e)}-\omega_{2m}^ky_k^{(o)}.
\]

These two formulas constitute one butterfly operation.

### 2.2 Matrix Factorization

Let $P_{2m}$ be the permutation matrix that places even-indexed entries before odd-indexed entries, and define

\[
D_m=\operatorname{diag}
\left(
1,\omega_{2m},\ldots,\omega_{2m}^{m-1}
\right).
\]

!!! success "Theorem 2.1 (Cooley-Tukey Factorization)"

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

This factorization reduces one $2m$-point DFT to two $m$-point DFTs plus $O(m)$ butterfly operations.

### 2.3 Complexity

Let $T(n)$ denote the operation count of an $n$-point FFT. For $n=2^t$,

\[
T(n)=2T(n/2)+O(n).
\]

Therefore

\[
T(n)=O(n\log_2 n).
\]

!!! note "Butterfly Operation"

    Given $a=y_k^{(e)}$, $b=y_k^{(o)}$, and the twiddle factor $w=\omega_{2m}^k$, one butterfly produces

    \[
    a+wb
    \]

    and

    \[
    a-wb.
    \]

---

## 3. Discrete Sine and Cosine Transforms

### 3.1 Discrete Sine Transform

!!! info "Definition 3.1 $DST-I$"

    For $x=(x_1,\ldots,x_{m-1})^{\mathsf T}$, define

    \[
    y_k=\sum_{j=1}^{m-1}\sin\left(\frac{kj\pi}{m}\right)x_j,
    \qquad k=1,\ldots,m-1.
    \]

    We write $y=S_{m-1}x$.

For example, when $m=3$,

\[
y_1=\sin\left(\frac{\pi}{3}\right)x_1
+\sin\left(\frac{2\pi}{3}\right)x_2,
\]

\[
y_2=\sin\left(\frac{2\pi}{3}\right)x_1
+\sin\left(\frac{4\pi}{3}\right)x_2.
\]

DST-I can be computed by taking an odd extension of $x$ and applying an FFT, so its complexity is $O(m\log m)$.

### 3.2 Discrete Cosine Transform

!!! info "Definition 3.2 $DCT-I$"

    For $x=(x_0,\ldots,x_m)^{\mathsf T}$, define

    \[
    y_k=\frac{x_0}{2}
    +\sum_{j=1}^{m-1}\cos\left(\frac{kj\pi}{m}\right)x_j
    +\frac{(-1)^kx_m}{2},
    \qquad k=0,\ldots,m.
    \]

DCT-I can be computed by taking an even extension of $x$ and applying an FFT.

### 3.3 Relationship with the DFT

For the DST, construct the odd extension of length $2m$:

\[
\widetilde{x}
=
(0,x_1,\ldots,x_{m-1},0,-x_{m-1},\ldots,-x_1)^{\mathsf T}.
\]

The imaginary part of its DFT gives the sine transform. For the DCT, construct the even extension

\[
\widehat{x}
=
(x_0,x_1,\ldots,x_{m-1},x_m,x_{m-1},\ldots,x_1)^{\mathsf T},
\]

The real part of its DFT gives the cosine transform.

!!! info "Connection with Boundary Conditions"

    An odd extension naturally corresponds to homogeneous Dirichlet boundary conditions, so the DST is commonly used for fixed-endpoint problems.

    An even extension naturally corresponds to homogeneous Neumann boundary conditions, so the DCT is commonly used for zero-normal-derivative problems.

---

## 4. Circulant Matrices

### 4.1 Cyclic Shift Matrix

Let

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

Then $D_n$ represents a cyclic shift and satisfies

\[
D_n^n=I.
\]

!!! info "Definition 4.1 (Circulant Matrix)"

    A circulant matrix can be written as

    \[
    C(z)=\sum_{k=0}^{n-1}z_kD_n^k.
    \]

    Every column is a cyclic shift of the preceding column.

Multiplication by a circulant matrix corresponds to discrete circular convolution.

### 4.2 Fourier Diagonalization

!!! success "Theorem 4.2 (Spectral Decomposition of a Circulant Matrix)"

    The DFT matrix simultaneously diagonalizes all circulant matrices of order $n$. In particular,

    \[
    F_nD_nF_n^{-1}
    =\operatorname{diag}
    \left(
    1,\omega_n,\ldots,\omega_n^{n-1}
    \right).
    \]

    Therefore

    \[
    F_nC(z)F_n^{-1}
    =\operatorname{diag}(F_nz).
    \]

??? proof "Proof (click to expand)"

    The $j$th Fourier vector is an eigenvector of $D_n$ with eigenvalue $\omega_n^j$. Hence $F_n$ diagonalizes $D_n$. Since $C(z)$ is a polynomial in $D_n$, the same eigenvectors diagonalize $C(z)$. Its eigenvalues are the values of that polynomial at the $n$th roots of unity, namely the entries of $F_nz$. $\square$

If $C(z)$ is nonsingular, the system

\[
C(z)x=b
\]

can be solved in $O(n\log n)$ time using two FFTs and one componentwise division.

---

## 5. One-Dimensional Poisson Equation

Consider the boundary-value problem

\[
-u''(x)=f(x),
\qquad \alpha<x<\beta,
\]

with Dirichlet boundary conditions

\[
u(\alpha)=u_\alpha,
\qquad
u(\beta)=u_\beta.
\]

Let

\[
h=\frac{\beta-\alpha}{m},
\qquad
x_i=\alpha+ih,
\qquad i=0,\ldots,m.
\]

Using the centered difference formula at the interior grid points gives

\[
\frac{-u_{i-1}+2u_i-u_{i+1}}{h^2}=f_i,
\qquad i=1,\ldots,m-1.
\]

Define

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

The discrete system is

\[
T_{m-1}^{(DD)}u=b,
\]

where

\[
u=(u_1,\ldots,u_{m-1})^{\mathsf T},
\]

\[
b=h^2(f_1,\ldots,f_{m-1})^{\mathsf T}
+(u_\alpha,0,\ldots,0,u_\beta)^{\mathsf T}.
\]

!!! success "Theorem 5.1 (Sine Diagonalization of the Dirichlet Difference Matrix)"

    The eigenvectors of $T_{m-1}^{(DD)}$ are given by DST-I, and its eigenvalues are

    \[
    \lambda_k=2-2\cos\left(\frac{k\pi}{m}\right)
    =4\sin^2\left(\frac{k\pi}{2m}\right),
    \qquad k=1,\ldots,m-1.
    \]

Thus one may apply a DST to $b$, divide componentwise by $\lambda_k$, and apply the inverse DST. The total complexity is $O(m\log m)$.

---

## 6. Two-Dimensional Poisson Equation and the Kronecker Product

### 6.1 Five-Point Difference Scheme

Let

\[
\Omega={(x,y):\alpha_x\leq x\leq\beta_x,
\ \alpha_y\leq y\leq\beta_y},
\]

and consider

\[
-\Delta u=f,
\qquad
\Delta=\frac{\partial^2}{\partial x^2}
+\frac{\partial^2}{\partial y^2}.
\]

Set

\[
h_x=\frac{\beta_x-\alpha_x}{m_x+1},
\qquad
h_y=\frac{\beta_y-\alpha_y}{m_y+1}.
\]

At an interior grid point, the five-point difference scheme is

\[
\frac{-u_{i-1,j}+2u_{i,j}-u_{i+1,j}}{h_x^2}
+\frac{-u_{i,j-1}+2u_{i,j}-u_{i,j+1}}{h_y^2}
=f_{i,j}.
\]

### 6.2 Kronecker Product Representation

!!! info "Definition 6.1 (Kronecker Product)"

    If $A=(a_{ij})\in\mathbb{R}^{m\times n}$, then

    \[
    A\otimes B=
    \begin{bmatrix}
    a_{11}B&\cdots&a_{1n}B\\
    \vdots&\ddots&\vdots\\
    a_{m1}B&\cdots&a_{mn}B
    \end{bmatrix}.
    \]

After stacking the interior unknowns column by column, the two-dimensional discrete system becomes

\[
Au=b,
\]

where

\[
A=\frac{1}{h_x^2}
\left(I_{m_y}\otimes T_{m_x}^{(DD)}\right)
+\frac{1}{h_y^2}
\left(T_{m_y}^{(DD)}\otimes I_{m_x}\right).
\]

This is the Kronecker sum of two one-dimensional discrete operators.

### 6.3 Separable Diagonalization

If

\[
T_{m_x}^{(DD)}=S_x^{-1}\Lambda_xS_x,
\]

\[
T_{m_y}^{(DD)}=S_y^{-1}\Lambda_yS_y,
\]

then

\[
(S_y\otimes S_x)A(S_y^{-1}\otimes S_x^{-1})
=\frac{1}{h_x^2}(I\otimes\Lambda_x)
+\frac{1}{h_y^2}(\Lambda_y\otimes I),
\]

and the right-hand side is diagonal. Hence the two-dimensional Poisson equation can be solved rapidly by applying DSTs in both coordinate directions, dividing pointwise, and then applying inverse DSTs.

---

## 7. Chapter Summary

!!! success "Main Conclusions"

    1. The DFT matrix satisfies $F_n^{\mathsf H}F_n=nI$, and its inverse follows directly from its conjugate transpose.

    2. The radix-2 FFT recursively splits an $n$-point DFT into two $n/2$-point DFTs, reducing the cost from $O(n^2)$ to $O(n\log n)$.

    3. The DST and DCT correspond to odd and even extensions, respectively, and can be computed using FFTs.

    4. The DFT simultaneously diagonalizes all circulant matrices, so circular convolution can be performed efficiently by FFT.

    5. The one-dimensional Dirichlet Poisson difference matrix is diagonalized by the DST.

    6. The two-dimensional Poisson discretization is a Kronecker sum of one-dimensional operators and can be solved using separable transforms.


