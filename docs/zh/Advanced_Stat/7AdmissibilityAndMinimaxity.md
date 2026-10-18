# 第七章：可容许性、MRU 与 Minimax

本章讨论统计决策理论中的四类最优性准则：Bayes、可容许性、最小风险无偏估计（MRU）与 minimax。随后给出 Bayes 规则推出可容许性的三条途径、Blyth 方法、最不利先验，并通过 Bernoulli 与正态均值问题比较这些准则之间的区别。

---

## 1. 四种最优性准则

设频率学派风险函数为

\[
R(\theta,\delta)
=E_\theta L\{\theta,\delta(X)\}.
\]

对于先验概率分布 $\pi$，定义 Bayes 风险

\[
r(\pi,\delta)
=\int_\Theta R(\theta,\delta)\,\pi(d\theta),
\qquad
B(\pi)=\inf_\delta r(\pi,\delta).
\]

对于一个非零但可能非正常的测度 $\nu$，同一个积分记作 $r(\nu,\delta)$。若某个规则在有限积分风险的规则中使该积分最小，则称其为关于 $\nu$ 的广义 Bayes 规则。

| 准则 | 定义 | 所防范的问题 |
| --- | --- | --- |
| Bayes | 最小化先验平均风险 | 所选先验下的平均表现 |
| 可容许 | 不被任何规则一致支配 | 在每个参数点都可以避免的改进 |
| MRU | 在无偏规则中使风险最小 | 无偏类中的方差 |
| Minimax | 最小化 $\sup_\theta R(\theta,\delta)$ | 最坏参数值 |

这些准则彼此有关，但不能互相混用。特别地，MRU 是相对于无偏规则这个受限类别定义的；可容许性与 minimax 通常需要同所有决策规则比较。

---

## 2. Bayes 规则何时可容许

如果

\[
R(\theta,\delta_1)
\leq R(\theta,\delta_0)
\qquad
\text{对每个 }\theta
\]

成立，并且至少在一个参数点严格不等，则称 $\delta_1$ 支配 $\delta_0$。若不存在支配某规则的其他规则，则称该规则可容许。

!!! success "定理 2.1（从 Bayes 性到可容许性的三条途径）"

    设 $\delta_\nu$ 关于非零测度 $\nu$ 最小化 $r(\nu,\delta)$。下列任一组条件都足以推出 $\delta_\nu$ 可容许。

    **(i) 连续性与满支撑：**每个风险函数都连续；$\Theta$ 的每个非空开子集都有正的 $\nu$-测度；并且至少有一个规则的积分风险有限。

    **(ii) 强唯一性：**$r(\nu,\delta_\nu)<\infty$，并且 $\delta_\nu$ 在如下决策论意义下是唯一 Bayes 规则：任意其他最小化规则 $\delta$ 都满足

    \[
    P_\theta\{\delta(X)=\delta_\nu(X)\}=1
    \qquad
    \text{对每个 }\theta.
    \]

    **(iii) 可数参数空间：**$\Theta$ 可数，并且

    \[
    0<\nu\{\theta\}<\infty
    \qquad
    \text{对每个 }\theta,
    \]

    同时 $r(\nu,\delta_\nu)<\infty$。

    在条件 (i) 中，只要积分风险定义良好，$\nu$ 可以是非正常的 $\sigma$-有限测度。

??? proof "定理 2.1 的证明（点击展开）"

    **情形 (i)。**最小化性质及有限风险竞争者的存在性推出

    \[
    r(\nu,\delta_\nu)<\infty.
    \]

    假设 $\delta$ 支配 $\delta_\nu$，并在 $\theta_0$ 处严格改进。连续函数

    \[
    g(\theta)
    =R(\theta,\delta_\nu)-R(\theta,\delta)
    \]

    非负且在 $\theta_0$ 处为正。因此，在 $\theta_0$ 的某个非空开邻域 $U$ 上存在 $\epsilon>0$，使得

    \[
    g(\theta)\geq\epsilon.
    \]

    满支撑条件给出 $\nu(U)>0$，从而

    \[
    r(\nu,\delta_\nu)-r(\nu,\delta)
    =\int_\Theta g(\theta)\,\nu(d\theta)>0,
    \]

    这与 $\delta_\nu$ 的最优性矛盾。

    **情形 (ii)。**支配关系给出

    \[
    r(\nu,\delta)\leq r(\nu,\delta_\nu),
    \]

    而最优性给出反向不等式，因此 $\delta$ 也是 Bayes 规则。强唯一性说明二者对每个 $\theta$ 都 $P_\theta$-几乎处处相等，因而不可能存在严格的风险改进。

    **情形 (iii)。**若在 $\theta_0$ 处有严格改进，则可数和

    \[
    r(\nu,\delta_\nu)-r(\nu,\delta)
    =\sum_{\theta\in\Theta}
    \{R(\theta,\delta_\nu)-R(\theta,\delta)\}\nu\{\theta\}
    \]

    中对应 $\theta_0$ 的项为正，其余各项非负。这再次与 Bayes 最优性矛盾。 $\square$

!!! note "点质量先验不一定危险"

    满支撑只是一个充分条件，并非必要条件。若能验证定理 2.1(ii) 中的强唯一性，则点质量先验下的 Bayes 规则仍可能是可容许的。

---

## 3. Bernoulli 抽样与非正常先验

设

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}
\operatorname{Bernoulli}(\theta),
\qquad
0\leq\theta\leq1,
\qquad
T=\sum_{i=1}^nX_i,
\]

并采用平方误差损失。考虑非正常测度

\[
\nu(d\theta)
=\frac{d\theta}{\theta(1-\theta)}.
\tag{1}
\]

当 $1\leq t\leq n-1$ 时，形式后验为

\[
\theta\mid T=t
\sim\operatorname{Beta}(t,n-t),
\qquad
E(\theta\mid T=t)=\frac{t}{n}.
\tag{2}
\]

但是，当 $t=0$ 或 $t=n$ 时，形式后验是非正常的。因此，如果不单独处理两个端点，直接说“后验均值为 $T/n$”并不完整。

!!! success "命题 3.1（严格的广义 Bayes 计算）"

    在先验测度 (1) 下，具有最小有限积分风险的唯一规则为

    \[
    \delta_0(T)=\frac{T}{n}.
    \]

    该规则可容许，并且在平方误差损失下也是 MRU 估计量。

??? proof "命题 3.1 的证明（点击展开）"

    对规则 $a_t=\delta(t)$，由 Tonelli 定理，

    \[
    r(\nu,\delta)
    =\sum_{t=0}^n
    \binom{n}{t}
    \int_0^1
    (a_t-\theta)^2
    \theta^{t-1}(1-\theta)^{n-t-1}\,d\theta.
    \tag{3}
    \]

    当 $1\leq t\leq n-1$ 时，第 $t$ 项由有限 Beta 测度的均值 $t/n$ 唯一最小化。

    当 $t=0$ 时，除非 $a_0=0$，否则 (3) 中的积分在零点附近发散；当 $a_0=0$ 时积分有限。类似地，$t=n$ 项有限当且仅当 $a_n=1$。因此，包括两个端点在内，$a_t=t/n$ 是唯一的有限风险最小化规则。

    它的频率学派风险与积分风险分别为

    \[
    R(\theta,T/n)
    =\frac{\theta(1-\theta)}{n},
    \qquad
    r(\nu,T/n)
    =\int_0^1\frac{1}{n}\,d\theta
    =\frac{1}{n}.
    \]

    风险函数连续，并且测度 (1) 对每个开子区间都赋予正测度。由定理 2.1(i)，$T/n$ 可容许。

    最后，$T/n$ 无偏，而 $T$ 是二项分布族的完全充分统计量。由 Lehmann-Scheffe 定理，$T/n$ 是 UMVU；等价地，它在平方误差损失下是 MRU。 $\square$

!!! warning "经验"

    非正常先验不是概率分布，形式后验也未必在每个样本点都存在。当出现例外样本时，广义 Bayes 结论应通过积分风险进行验证。

---

## 4. 正态均值：仿射规则 $c\overline{X}+d$

设

\[
X_1,\ldots,X_n
\overset{\mathrm{iid}}{\sim}
N(\theta,\sigma^2),
\qquad
v=\frac{\sigma^2}{n},
\qquad
Y=\overline{X}\sim N(\theta,v),
\]

其中 $\sigma^2$ 已知。对于

\[
\delta_{c,d}(Y)=cY+d,
\]

风险函数为

\[
R(\theta,\delta_{c,d})
=c^2v+\{(c-1)\theta+d\}^2.
\tag{4}
\]

!!! success "命题 4.1（仿射规则的完全分类）"

    在所有仿射规则中，并且事实上相对于任意竞争规则，下列规则可容许：

    \[
    0\leq c<1,
    \qquad
    d\in\mathbb{R},
    \qquad
    \text{以及 }(c,d)=(1,0).
    \]

    其他所有仿射规则均不可容许。

??? proof "命题 4.1 的证明（点击展开）"

    **当 $0<c<1$ 时。**选择先验

    \[
    \theta\sim N(\mu,\tau^2),
    \qquad
    \tau^2=\frac{cv}{1-c},
    \qquad
    \mu=\frac{d}{1-c}.
    \]

    后验均值为

    \[
    \frac{\tau^2}{\tau^2+v}Y
    +\frac{v}{\tau^2+v}\mu
    =cY+d.
    \]

    这是具有全实轴正密度的正常先验下的唯一 Bayes 规则，且风险连续有限，因此由定理 2.1 可容许。

    **当 $c=0$ 时。**常数规则 $d$ 是点质量先验 $\theta=d$ 下的 Bayes 规则。任何其他 Bayes 规则在 $\theta=d$ 处风险为零，因此与 $d$ 在 $P_d$-几乎处处意义下相等。所有一维正态分布都有正密度并且彼此绝对连续，所以该规则对每个 $\theta$ 都与 $d$ 在 $P_\theta$-几乎处处意义下相等。强唯一性推出可容许性。

    **当 $c>1$ 时。**由 (4)，

    \[
    R(\theta,cY+d)
    \geq c^2v
    >v
    =R(\theta,Y),
    \]

    因此 $Y$ 支配 $cY+d$。

    **当 $c<0$ 时。**令

    \[
    m=\frac{d}{1-c}.
    \]

    则

    \[
    R(\theta,cY+d)
    =c^2v+(1-c)^2(\theta-m)^2
    >(\theta-m)^2
    =R(\theta,m),
    \]

    所以常数规则 $m$ 支配 $cY+d$。

    **当 $c=1$ 且 $d\neq0$ 时。**

    \[
    R(\theta,Y+d)
    =v+d^2
    >v,
    \]

    因而 $Y$ 支配 $Y+d$。

    最后还需证明 $(c,d)=(1,0)$ 的情形。规则 $Y$ 是 Lebesgue 测度下的广义 Bayes 规则，但积分风险为无穷，因此不能直接使用定理 2.1(i)。它的可容许性将在第 5 节利用 Blyth 方法证明。 $\square$

统计量 $Y=\overline{X}$ 同时也是 MRU：它无偏，并且是关于 $\theta$ 的完全充分统计量，所以 Lehmann-Scheffe 定理适用。但这一事实本身尚不能证明它可容许或 minimax。

---

## 5. Blyth 方法：由先验序列证明可容许性

当目标规则是某个非正常测度下的广义 Bayes 规则，而普通 Bayes 风险论证失效时，可以使用 Blyth 方法。

!!! success "定理 5.1（Blyth 方法）"

    设 $\mathcal{D}_c$ 是风险函数连续的完全类，即 $\mathcal{D}_c$ 外的每个规则都被 $\mathcal{D}_c$ 中的某个规则弱支配。设 $\delta_0\in\mathcal{D}_c$。

    假设对每个 $j$，存在非零测度 $\nu_j$，使得 $r(\nu_j,\delta_0)<\infty$，并且对每个 $\theta_0$ 与 $\eta>0$，

    \[
    0<\nu_j\{B(\theta_0,\eta)\}<\infty,
    \]

    以及

    \[
    \frac{
    r(\nu_j,\delta_0)
    -\inf_{\delta\in\mathcal{D}_c}r(\nu_j,\delta)
    }{
    \nu_j\{B(\theta_0,\eta)\}
    }
    \longrightarrow0.
    \tag{5}
    \]

    则 $\delta_0$ 可容许。

??? proof "定理 5.1 的证明（点击展开）"

    假设 $\delta_1\in\mathcal{D}_c$ 支配 $\delta_0$，并在 $\theta_0$ 处严格改进。由连续性，存在 $\eta>0$ 与 $\epsilon>0$，使得

    \[
    R(\theta,\delta_0)-R(\theta,\delta_1)
    \geq\epsilon
    \qquad
\text{在集合 }B(\theta_0,\eta)\text{ 上}.
    \]

    因而

    \[
    \begin{aligned}
    r(\nu_j,\delta_0)
    -\inf_{\delta\in\mathcal{D}_c}r(\nu_j,\delta)
    &\geq
    r(\nu_j,\delta_0)-r(\nu_j,\delta_1)\\
    &\geq
    \epsilon\,\nu_j\{B(\theta_0,\eta)\}.
    \end{aligned}
    \]

    所以 (5) 中的比值对每个 $j$ 都至少为 $\epsilon$，矛盾。完全类性质则排除了来自 $\mathcal{D}_c$ 之外的支配规则。 $\square$

### 5.1 应用：一维正态均值的可容许性

回到

\[
Y\sim N(\theta,v),
\qquad
\delta_0(Y)=Y.
\]

令

\[
\nu_j=N(0,j).
\]

相应的 Bayes 规则和 Bayes 风险为

\[
\delta_j(Y)=\frac{j}{j+v}Y,
\qquad
B(\nu_j)=\frac{vj}{v+j}.
\]

由于 $r(\nu_j,Y)=v$，Bayes 风险差为

\[
r(\nu_j,Y)-B(\nu_j)
=\frac{v^2}{v+j}
=O(j^{-1}).
\tag{6}
\]

对于固定区间

\[
B(\theta_0,\eta)
=(\theta_0-\eta,\theta_0+\eta),
\]

有

\[
\begin{aligned}
\nu_j\{B(\theta_0,\eta)\}
&=
\Phi\left(\frac{\theta_0+\eta}{\sqrt{j}}\right)
-\Phi\left(\frac{\theta_0-\eta}{\sqrt{j}}\right)\\
&\sim\frac{2\eta\phi(0)}{\sqrt{j}}.
\end{aligned}
\tag{7}
\]

将 (6) 除以 (7)，得到

\[
O(j^{-1/2})\longrightarrow0.
\]

在该一维正态问题的标准完全类正则条件下，由定理 5.1，

\[
Y=\overline{X}
\]

可容许。

!!! note "Blyth 论证的核心"

    关键在于比较收敛速度：Bayes 风险差以 $j^{-1}$ 消失，比任意固定邻域的先验质量 $j^{-1/2}$ 消失得更快。只说正常先验“趋近于平坦先验”并不足以完成证明。

---

## 6. Minimax 性与最不利先验

定义

\[
M(\delta)
=\sup_{\theta\in\Theta}R(\theta,\delta),
\qquad
V=\inf_\delta M(\delta).
\]

对每个正常先验 $\pi$ 和任意规则 $\delta$，

\[
B(\pi)
\leq r(\pi,\delta)
\leq M(\delta).
\tag{8}
\]

因此，每个 Bayes 风险都是 minimax 风险的下界。

!!! success "定理 6.1（证明 minimax 性的两种 Bayes 方法）"

    **(i)** 若 $\delta_\pi$ 是正常先验 $\pi$ 下的 Bayes 规则，并且

    \[
    B(\pi)
    =r(\pi,\delta_\pi)
    =M(\delta_\pi)
    <\infty,
    \tag{9}
    \]

    则 $\delta_\pi$ 是 minimax 规则，且 $\pi$ 是最不利先验。特别地，具有常风险的正常 Bayes 规则是 minimax 的。

    **(ii)** 若 $\{\pi_j\}$ 是正常先验序列，并且

    \[
    B(\pi_j)
    \longrightarrow M(\delta_0)
    <\infty,
    \tag{10}
    \]

    则 $\delta_0$ 是 minimax 规则，且 $\{\pi_j\}$ 是最不利先验序列。

??? proof "定理 6.1 的证明（点击展开）"

    在 (9) 下，对任意 $\delta$，由 (8)，

    \[
    M(\delta)
    \geq B(\pi)
    =M(\delta_\pi),
    \]

    这证明了 minimax 性。

    此外，对任意先验 $\pi'$，

    \[
    B(\pi')
    \leq V
    =M(\delta_\pi)
    =B(\pi),
    \]

    因而 $\pi$ 最不利。

    在 (10) 下，(8) 对每个 $j$ 与每个 $\delta$ 给出

    \[
    B(\pi_j)\leq M(\delta).
    \]

    令 $j\to\infty$，得到

    \[
    M(\delta_0)\leq M(\delta),
    \]

    从而证明 $\delta_0$ minimax，并同时证明 $\{\pi_j\}$ 是最不利先验序列。 $\square$

---

## 7. 二项分布例子：MRU 不一定 Minimax

设

\[
T\sim\operatorname{Bin}(n,\theta),
\]

在平方误差损失下估计 $\theta$。MRU 估计量 $T/n$ 的风险为

\[
R(\theta,T/n)
=\frac{\theta(1-\theta)}{n},
\qquad
M(T/n)=\frac{1}{4n}.
\]

由命题 3.1，它可容许，但它不是 minimax 规则。

取正常先验

\[
\theta\sim\operatorname{Beta}(a,a).
\]

Bayes 规则为

\[
\delta_a(T)
=\frac{T+a}{n+2a}.
\]

对任意 $a>0$，

\[
R(\theta,\delta_a)
=\frac{
n\theta(1-\theta)+a^2(1-2\theta)^2
}{
(n+2a)^2
}.
\tag{11}
\]

选择

\[
a=\frac{\sqrt{n}}{2}.
\]

此时 $a^2=n/4$，而 (11) 的分子为

\[
n\left\{
\theta(1-\theta)
+\frac{1}{4}(1-2\theta)^2
\right\}
=\frac{n}{4}.
\]

所以

\[
R(\theta,\delta_a)
=\frac{1}{4(\sqrt{n}+1)^2}
\qquad
\text{对每个 }\theta.
\]

该规则既是正常 Bayes 规则，又具有常风险。由定理 6.1(i)，它是 minimax 规则。而且

\[
\frac{1}{4(\sqrt{n}+1)^2}
<\frac{1}{4n}
=M(T/n).
\]

!!! info "解释"

    MRU 规则拒绝引入偏差，并在 $\theta=0,1$ 处具有零风险，但代价是在 $\theta=1/2$ 附近承担更大的风险。Minimax 规则允许在边界附近产生偏差，以压平并降低最坏情形风险。两个规则都是可容许的。

---

## 8. 正态例子：相同先验，不同结论

对于

\[
Y\sim N(\theta,v),
\]

规则 $Y$ 的风险恒为 $v$。取

\[
\pi_j=N(0,j).
\]

则

\[
B(\pi_j)
=\frac{vj}{v+j}
\longrightarrow v
=M(Y).
\]

由定理 6.1(ii)，$Y$ 是 minimax 规则。

同一个先验序列承担了两种不同作用：

- 对于 minimax 性，比较绝对 Bayes 风险 $B(\pi_j)$ 与最大风险 $M(Y)$；

- 对于可容许性，比较 Bayes 风险差 $r(\pi_j,Y)-B(\pi_j)$ 与每个邻域的局部先验质量。

Minimax 性并不推出可容许性。对于维数至少为三的多元正态均值，通常估计量是 minimax 的，却被 James-Stein 收缩估计量支配。该结论属于第五章，此处不再重复。

---

## 9. 需要记住的关系

!!! success "命题 9.1（常风险与唯一性）"

    **(i)** 具有有限常风险的可容许规则是 minimax 的。

    **(ii)** 唯一的 minimax 规则是可容许的。

    **(iii)** 在严格凸损失和凸行动空间下，具有常风险的可容许规则是唯一的 minimax 规则，忽略几乎处处相等的差异。

??? proof "命题 9.1 的证明（点击展开）"

    设

    \[
    R(\theta,\delta_0)=C.
    \]

    如果 $\delta_0$ 不是 minimax 规则，则存在 $\delta_1$ 使得

    \[
    M(\delta_1)<C.
    \]

    因而 $\delta_1$ 支配 $\delta_0$，与可容许性矛盾。这证明了 (i)。

    如果唯一的 minimax 规则不可容许，则支配它的规则具有不更大的最大风险，因而也是 minimax 规则，与唯一性矛盾。这证明了 (ii)。

    对于 (iii)，设 $\delta_1$ 是另一个 minimax 规则。由行动空间凸性，中点规则

    \[
    \delta_m
    =\frac{\delta_0+\delta_1}{2}
    \]

    仍然允许。由损失函数的严格凸性，

    \[
    R(\theta,\delta_m)
    \leq
    \frac{1}{2}R(\theta,\delta_0)
    +\frac{1}{2}R(\theta,\delta_1)
    \leq C,
    \]

    并且在两个规则以正的 $P_\theta$-概率不同的地方严格不等。如果它们对某个 $\theta$ 不同，则 $\delta_m$ 支配 $\delta_0$，与可容许性矛盾。因此 minimax 规则在几乎处处相等的意义下唯一。 $\square$

---

## 10. 本章总结

!!! success "核心结论"

    1. Bayes、可容许、MRU 与 minimax 是不同的最优性概念，比较的规则类别与风险目标不同。

    2. 连续性与满支撑、强唯一性或可数参数空间都可以把 Bayes 最优性提升为可容许性。

    3. 非正常先验需要通过积分风险严格处理，尤其要检查形式后验失效的样本点。

    4. Blyth 方法通过比较 Bayes 风险差与局部先验质量的收敛速度证明广义 Bayes 规则的可容许性。

    5. Bayes 风险为 minimax 风险提供下界；常风险 Bayes 规则或最不利先验序列可用于证明 minimax 性。

    6. MRU 不必是 minimax，minimax 也不必可容许；附加的常风险、唯一性与凸性条件可以建立更强关系。
