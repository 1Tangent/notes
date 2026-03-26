## 先验知识

### 1.1 概率论相关知识



### 1.2 马尔可夫过程



### 1.3 矩阵相关
#### 1.3.1 向量

#### 5.0.1 矩阵指数
矩阵指数可由线性微分方程的解的形式导出。这里我们先考察一个一阶线性微分方程
$$
\dot x(t)= ax(t)
$$
其中，$x(t) \in \mathbb{R}$，$a \in R$，且初始条件为$x(0) = x_0$，那么这个方程的解为
$$
x(t) = e^{at}x_0
$$
然后我们考察一阶齐次线性微分方程组
$$
\left\{
\begin{array}{c}
\frac{d}{dt} x_{1}(t)=a_{11} x_{1}(t)+a_{12} x_{2}(t)+\cdots+a_{1n} x_{n}(t) \\
\frac{d}{dt} x_{2}(t)=a_{21} x_{1}(t)+a_{22} x_{2}(t)+\cdots+a_{2n} x_{n}(t) \\
\cdots \\
\frac{d}{dt} x_{n}(t)=a_{n1} x_{1}(t)+a_{n2} x_{2}(t)+\cdots+a_{nn} x_{n}(t)
\end{array}\right.
$$
我们用矩阵写出来，即
$$
\frac{d}{dt}x(t) =  Ax(t)
$$
其中 $x(t) = \begin{bmatrix} x_1(t) \\ x_2(t) \\ \cdots \\ x_n(t)\end{bmatrix}$，$A=(a_{ij})_{n\times n}$
可以发现，这个形式与前面的一阶线性微分方程形式类似，那我们也希望预期的方程组的解也有类似的形式，即$e^{(\cdots)}$。从形式上来看，若$x(t) = e^{At}x(0)$，则有
$$
\frac{d}{dt}x(t) = \frac{d}{dt}(e^{At})x(0) = Ae^{At}x(0) = Ax(t)
$$
可以发现$x(t) = e^{At}x(0)$在形式上确实是方程的解，那我们现在试图找到这样的矩阵函数。

首先我们回忆对指数函数$e^x$在级数上的定义：

对于$x \in \mathbb{R}$，指数函数$e^x$可以用收敛的级数刻画，即
$$
e^x = \sum_{n=0}^{+\infty}\frac{x^n}{n!}
$$
那么类似地，对于矩阵而言
$$
e^A = \sum_{n=0}^{+\infty}\frac{A^n}{n!}
$$
让我们验证一下：
$$
\frac{d}{dt}e^{At} = \frac{d}{dt} (\sum_{n=0}^{+\infty}\frac{(At)^n}{n!}) = \sum_{n=0}^{+\infty}\frac{d}{dt}\frac{(At)^n}{n!} = \sum_{n=1}^{+\infty}\frac{A^nt^{n-1}}{(n-1)!} = Ae^{At}
$$
即可以验证对于方程$\frac{d}{dt}x(t) =  Ax(t)$，有形如$x(t) = e^{At}x(0)$的解。

下面证明$e^A$是良好定义的，即级数$\underset{n=0}{\overset{+\infty}{\sum}}\frac{A^n}{n!}$收敛。

任取$i,j(1\le i,j \le n)$，记矩阵$\frac{A^n}{n!}$的第$i$行，第$j$列的元素为$b_n$，即
$e^A$的第$i$行，第$j$列的元素为$\sum_{n=0}^{+\infty}b_n$，所以只需证明该级数收敛。
取矩阵最大值范数$\|A\|=\underset{1 \leq i, j \leq n}{\max} \left|a_{i j}\right|$，则有
$$
|b_n| \le \| \frac{A^n}{n!} \| \le \frac{\|A\|^n}{n!}
$$
$n\rightarrow \infty$时，$\frac{\|A\|^n}{n!} \rightarrow 0$，故$b \rightarrow 0,n\rightarrow\infty$
$$
\sum_{n=0}^{+\infty}|b_n| \le \sum_{n=0}^{+\infty}\| \frac{A^n}{n!} \| \le \sum_{n=0}^{+\infty}\frac{\|A\|^n}{n!} = e^{\|A\|} < +\infty 
$$
因此，级数$\underset{n=0}{\overset{+\infty}{\sum}}|b_n|$收敛，原命题得证。
下面介绍矩阵指数的一些性质：

1. 若$A$的特征值为$\lambda_1,\lambda_2\cdots\lambda_n$，则$e^A$的特征值为$e^{\lambda_1},e^{\lambda_2},\cdots,e^{\lambda_n}$
  
   可以得到$e^A$的行列式为：
   $$
   \det(e^A) = e^{\lambda_1}e^{\lambda_2}\cdots e^{\lambda_n} = e^{trA}
   $$
   显然可以得到$e^A$可逆。
2. 若$AB=BA$，则$e^Ae^B=e^{A+B}$
  
   证明如下：
   $$
   e^{A+B} = \sum_{n=0}^{+\infty}\frac{(A+B)^n}{n!} = 		 
   \sum_{n=0}^{+\infty}\sum_{k=0}^{n}\frac{C^k_nA^kB^{n-k}}{n!}
   = \sum_{n=0}^{+\infty}\sum_{k=0}^{n}\frac{A^kB^{n-k}}{k!(n-k)!}
   $$
   $$
   \begin{aligned}
   e^Ae^B&= (\sum_{n=0}^{+\infty}\frac{A^n}{n!})(\sum_{n=0}^{+\infty}\frac{B^n}{n!})\\
   &=(I+A+\frac{A^2}{2!}+\cdots+\frac{A^m}{m!}+\cdots)
   (I+B+\frac{B^2}{2!}+\cdots+\frac{B^m}{m!}+\cdots) \\
   &=I+(A+B)+(\frac{A^2}{2!} + AB + BA + \frac{B^2}{2!}) + \cdots \\
   &+ (\frac{A^m}{m!}+\frac{A^{m-1}}{(m-1)!}\cdot B + \cdots + 
   \frac{A^{m-k}}{(m-k)!}\cdot \frac{B^k}{k!} + \cdots + \frac{B^n}{n!}) \\
   &= \sum_{n=0}^{+\infty}\sum_{k=0}^{n}\frac{A^{n-k}B^{k}}{k!(n-k)!} \\
   &= \sum_{n=0}^{+\infty}\sum_{k=0}^{n}\frac{A^kB^{n-k}}{k!(n-k)!}
   \end{aligned}
   $$
   （如果对上述推导在理解上有困难，可以类比在高中学的二项式定理）

   即
   $$e^Ae^B=e^{A+B}$$
3. 对于矩阵指数函数$e^{At}$，可以将其展开成无穷级数的形式：
   $$
   e^{At}=I+At+\frac{(At)^2}{2!}+\frac{(At)^3}{3!}+\cdots
   $$
4. 若$A$的$Jordan$矩阵为$J$，则令$A=PJP^{-1}$，可得
   $$
   P^{-1}e^{At}P = e^{Jt} 
   $$
   注：若$A$的特征值为$\lambda_1,\lambda_2\cdots\lambda_n$，特征向量为 
   $\boldsymbol{v}_1,\boldsymbol{v}_2\cdots \boldsymbol{v}_n$，则 
	$J =  \begin{bmatrix} \lambda_1 & 0 & \cdots & 0 \\ 0 & \lambda_2 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & \lambda_n\end{bmatrix}$， 
	$P = \begin{bmatrix} \boldsymbol{v}_1,\boldsymbol{v}_2\cdots \boldsymbol{v}_n \end{bmatrix}$
#### 5.0.2 矩阵求导

##### 5.0.2.1 函数与标量、向量、矩阵

我们先约定符号如下表所示：
|函数&输入|标量|向量|矩阵|
|:---|:---:|:---:|:---:|
|标量|$f(x)$|$f(\boldsymbol{x})$|$f(X)$|
|向量|$\boldsymbol{f}(x)$|$\boldsymbol{f}(\boldsymbol{x})$|$\boldsymbol{f}(X)$|
|矩阵|$F(x)$|$F(\boldsymbol{x})$|$F(X)$|

##### 5.0.2.2 矩阵求导的本质

我们在数学分析学过对于一个多元函数$f(x_1, x_2, x_3)$，我们可以求$f$对$x_1,x_2,x_3$的偏导得到
$\frac{\partial f}{\partial x_1},\frac{\partial f}{\partial x_2},\frac{\partial f}{\partial x_3}$。把它写成向量的形式，令$\boldsymbol{x}=\begin{bmatrix} x_1 & x_2 & x_3 \end{bmatrix}^T$，则有
$$\frac{\partial f(\boldsymbol{x})}{\partial \boldsymbol{x}} = 
\begin{bmatrix} 
	\frac{\partial f}{\partial x_1} &
	\frac{\partial f}{\partial x_2} &
	\frac{\partial f}{\partial x_3} 
\end{bmatrix}^T$$

推广到$n$维，可以定义标量方程$f(\boldsymbol{x})$对向量$\boldsymbol{x}$的偏导为
$$
\frac{\partial f(\boldsymbol{x})}{\partial \boldsymbol{x}} = 
\begin{bmatrix} 
	\frac{\partial f}{\partial x_1} \\
	\frac{\partial f}{\partial x_2} \\
	\vdots \\
	\frac{\partial f}{\partial x_n} 
\end{bmatrix}
$$
这种写法的维度和列向量相同，又被称为**分母布局**。

当我们把写法改成$\begin{bmatrix}\frac{\partial f}{\partial x_1} & \frac{\partial f}{\partial x_2} &\cdots & \frac{\partial f}{\partial x_n} \end{bmatrix}$ 的时候，这时的维度和分子（方程$f(\boldsymbol{x})$）相同，这种写法被称为**分子布局**。

现在我们推广到更高阶的情况，假定我们的函数中有$m$个$f$，输入有$n$个元素，那么我们对每个$f_1,\cdots, f_m$对$x_1,\cdots,x_n$求偏导，总共有$m\times n$个结果。由此可以得出，矩阵求导就是每个$f$对$x$中每个元素求偏导。

#### 5.0.2.3 矩阵求导结果的布局和雅可比$(Jacobian)$矩阵

这里我们只给出$f(\boldsymbol{x})$和$\boldsymbol{f}(\boldsymbol{x})$的结果形式，更多表示请查看[这里](https://zhuanlan.zhihu.com/p/263777564)

1. 向量变元的标量函数
   2. 行向量偏导形式
      $$D_xf(\boldsymbol{x})=\frac{\partial f(\boldsymbol{x})}{\partial \boldsymbol{x}} = \begin{bmatrix} \frac{\partial f}{\partial x_1} & \frac{\partial f}{\partial x_2} & \cdots & \frac{\partial f}{\partial x_n} \end{bmatrix}$$
   3. 梯度向量（列向量偏导）形式
      $$\nabla_xf(\boldsymbol{x})=\frac{\partial f(\boldsymbol{x})}{\partial \boldsymbol{x}} = \begin{bmatrix} \frac{\partial f}{\partial x_1} & \frac{\partial f}{\partial x_2} & \cdots & \frac{\partial f}{\partial x_n} \end{bmatrix}^T$$
	   这两种情况互为转置。
4. 向量变元的向量函数$\boldsymbol{f}_{m \times 1}(\boldsymbol{x}_{n\times 1})$
  
   我们使用行向量偏导形式，可以得到
   $$
   \begin{aligned}
   \frac{\partial \boldsymbol{f}(\boldsymbol{x})}{\partial \boldsymbol{x}} 
   &= \begin{bmatrix} \frac{\partial \boldsymbol{f}}{\partial x_1} & \frac{\partial \boldsymbol{f}}{\partial x_2} & \cdots & \frac{\partial \boldsymbol{f}}{\partial x_n} \end{bmatrix} \\ 
   &=\begin{bmatrix}\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n} \\ 
   \frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n} \\ 
   \vdots & \vdots & \ddots & \vdots \\ 
   \frac{\partial f_m}{\partial x_1} & \frac{\partial f_m}{\partial x_2} & \cdots & \frac{\partial f_m}{\partial x_n}
   \end{bmatrix}
   \end{aligned}
   $$
   上述矩阵又被称为**雅可比矩阵$(Jacobian$ $Matrix)$**

#### 5.0.2.4 矩阵求导的一些基本性质

1. 向量变元的标量函数
   + $\frac{\partial c}{\partial\boldsymbol{x}} = \boldsymbol{0}_{n\times1}$
   + $\frac{\partial(c_1f(\boldsymbol{x})+c_2f(\boldsymbol{x}))}{\partial\boldsymbol{x}} = c_1 \frac{\partial f}{\partial\boldsymbol{x}} + c_2 \frac{\partial g}{\partial\boldsymbol{x}}$
   + $\frac{\partial(f(\boldsymbol{x})g(\boldsymbol{x}))}{\partial\boldsymbol{x}} =\frac{\partial f}{\partial\boldsymbol{x}}g(\boldsymbol{x}) + f(\boldsymbol{x})\frac{\partial g}{\partial\boldsymbol{x}}$
   + $\frac{\partial(\frac{f(\boldsymbol{x})}{g(\boldsymbol{x})})}{\partial\boldsymbol{x}} =\frac{1}{g^2(\boldsymbol{x})}(\frac{\partial f}{\partial\boldsymbol{x}}g(\boldsymbol{x}) - f(\boldsymbol{x})\frac{\partial g}{\partial\boldsymbol{x}})$
2. 向量变元的一些常用公式
   + 对于两个$n\times1$向量$\boldsymbol{x}, \boldsymbol{a}$，有$\frac{\partial(\boldsymbol{x}^T\boldsymbol{a})}{\partial\boldsymbol{x}} = \frac{\partial(\boldsymbol{a}^T\boldsymbol{x})}{\partial\boldsymbol{x}} = \boldsymbol{a}$
   + $\frac{\partial(\boldsymbol{x}^T\boldsymbol{x})}{\partial\boldsymbol{x}} =2\boldsymbol{x}$
   + 对于$A_{n\times n}$和$\boldsymbol{x}_{n\times1}$，有$\frac{\partial(A\boldsymbol{x})}{\partial\boldsymbol{x}} =A^T$，$\frac{\partial(\boldsymbol{x}^TA\boldsymbol{x})}{\partial\boldsymbol{x}} =A\boldsymbol{x}+A^T\boldsymbol{x}$，$\frac{\partial^2(\boldsymbol{x}^TA\boldsymbol{x})}{\partial\boldsymbol{x}^2} =A+A^T$
   + 对于三个$n\times1$向量$\boldsymbol{x},\boldsymbol{a},\boldsymbol{b}$，有$\frac{\partial(\boldsymbol{a}^T\boldsymbol{x}\boldsymbol{x}^T\boldsymbol{b})}{\partial\boldsymbol{x}} =\boldsymbol{a}\boldsymbol{b}^T\boldsymbol{x}+\boldsymbol{b}\boldsymbol{a}^T\boldsymbol{x}$
3. 矩阵变元的一些常用公式
   + 对于向量$\boldsymbol{a}_{m\times1},\boldsymbol{b}_{n\times1}$，矩阵$X_{m\times n}$，有
   $\frac{\partial(\boldsymbol{a}^TX\boldsymbol{b})}{\partial X}=\boldsymbol{a}\boldsymbol{b}^T$
   + 对于向量$\boldsymbol{a}_{m\times1},\boldsymbol{b}_{n\times1}$，矩阵$X_{m\times n}$，有
   $\frac{\partial(\boldsymbol{b}^TX^T\boldsymbol{a})}{\partial X}=\boldsymbol{a}\boldsymbol{b}^T$
   + 对于向量$\boldsymbol{a}_{m\times1},\boldsymbol{b}_{m\times1}$，矩阵$X_{m\times n}$，有
   $\frac{\partial(\boldsymbol{a}^TXX^T\boldsymbol{b})}{\partial X}=\boldsymbol{a}\boldsymbol{b}^TX + \boldsymbol{b}\boldsymbol{a}^TX$
   + 对于向量$\boldsymbol{a}_{n\times1},\boldsymbol{b}_{n\times1}$，矩阵$X_{m\times n}$，有
   $\frac{\partial(\boldsymbol{a}^TX^TX\boldsymbol{b})}{\partial X}=X\boldsymbol{b}\boldsymbol{a}^T + X\boldsymbol{a}\boldsymbol{b}^T$

以上公式推导过程在[这里](https://zhuanlan.zhihu.com/p/273729929)，这里建议先自己试着推导一下再看链接，也可以在`Eigen`库里面进行验证。

### 1.4 李群与李代数
#### 5.1.1 群、$\mathrm{SO}(3)$、$\mathrm{SE}(3)$

让我们先回忆一点离散数学的内容。群是一种集合加上一种运算的代数结构，记作$G=(A,\cdot)$，运算$\cdot$应该满足以下条件：

1. 封闭性： $\forall a_1,a_2 \in A,a_1 \cdot a_2 \in A$
2. 结合律：$\forall a_1,a_2,a_3 \in A,(a_1 \cdot a_2) \cdot a_3 = a_1 \cdot (a_2 \cdot a_3)$
3. 幺元：$\exists a_0 \in A, \forall a \in A, a \cdot a_0 = a_0 \cdot a = a$
4. 逆：$\forall a \in A, \exists a^{-1} \in A, a \cdot a^{-1} = a_0$

我们容易发现，由所有旋转矩阵构成的集合和矩阵乘法构成群，所有变换矩阵构成的集合和矩阵乘法构成群。这里我们对这两个群进行定义：
1. 三维特殊正交群 $\mathrm{SO}(3)$ 指旋转矩阵群，这里由于$\det(R) = 1$，所以称为特殊正交群
   $$
   \mathrm{SO}(3) = \{R \in \mathbb{R}^{3 \times 3} \mid RR^T=I,\det(R)=1\}
   $$
   这里我们也可以推广到$n$维，即
   $$
   \mathrm{SO}(n) = \{R \in \mathbb{R}^{n \times n} \mid RR^T=I,\det(R)=1\}
   $$

2. 三维特殊欧式群 $\mathrm{SE}(3)$ 指变换矩阵群
   $$
   \mathrm{SE}(3) = \{T = \begin{bmatrix} R & t \\ \boldsymbol{0}^T & 1
   \end{bmatrix} \in \mathbb{R}^{4 \times 4} \mid R \in \mathrm{SO}(3), t \in \mathbb{R}^3\}
   $$

**李群**是指具有连续性质（光滑）的群。我们可以直观想象一个可以在空间中连续运动的刚体，容易发现他对应的旋转矩阵和变换矩阵都是连续的，因而$\mathrm{SO}(n)$和$\mathrm{SE}(n)$都是李群。在这里我们主要讨论$\mathrm{SO}(3),\mathrm{SE}(3)$来进行相机的姿态估计。

#### 5.1.2 矩阵指数和旋转矩阵的关系

为了对李代数的理解更加深入，我们先寻找一下矩阵指数和旋转矩阵之间的关系。

![向量旋转](.\assets\向量旋转.png)

如上图所示，三维向量$\boldsymbol{p}(0)$绕着单位旋转轴$\hat\omega(\hat\omega\in\mathbb{R}^3, \|\hat\omega\|=1)$旋转$\theta$角度，得到三维向量$\boldsymbol{p}(\theta)$。该旋转运动也可以视作是三维向量$\boldsymbol{p}(0)$以$1\,rad/s$的角速度绕着单位旋转轴$\omega$从时间$t=0$运动到时间$t=\theta$得到三维向量$\boldsymbol{p}(\theta)$。

我们用$p(t)$表示旋转路径向量端点的位置，用$\dot p(t)$表示该点的瞬时速度值，则有
$$\dot p=\hat\omega \times p=\hat\omega^{\wedge}p$$
解这个方程，可得
$$p(t)=e^{\hat\omega^{\wedge}t}p(0)$$
由于这里$\theta = 1t$,可得
$$p(t)=e^{\hat\omega^{\wedge}\theta}p(0)$$
由$(\hat\omega^{\wedge})^3=-\hat\omega^{\wedge}$，对矩阵指数展开后，可得
$$
\begin{aligned}
e^{\hat\omega^{\wedge}\theta} 
&= I + 
(\theta - \frac{\theta^3}{3!}+\frac{\theta^5}{5!}-\cdots)\hat\omega^{\wedge} + 
(\frac{\theta^2}{2!} - \frac{\theta^4}{4!}+\frac{\theta^6}{6!}-\cdots)(\hat\omega^{\wedge})^2 \\
&= I + \sin\theta\hat\omega^{\wedge} + (1-\cos\theta)(\hat\omega^{\wedge})^2
\end{aligned}
$$
因此，三维旋转可以用矩阵指数的形式来表示。

#### 5.1.4 李代数的引出与定义

根据上文的描述，我们发现，三维旋转矩阵$R$可以用三维向量$\hat\omega\theta\in\mathbb{R}^3$的反对称矩阵的矩阵指数进行表示。如果称$\hat\omega\theta$为$R$对应的指数坐标，那么它的反对称矩阵$\hat\omega^\wedge\theta$称为$R$的矩阵对数，也称为$\mathrm{SO}(3)$对应的李代数$\mathfrak{so}(3)$。

可以发现，李代数描述了李群的局部性质，一般李代数的定义如下：

李代数由一个集合$\mathbb{V}$，一个数域$\mathbb{F}$和一个二元运算（**李括号**）$[,]$组成，如果满足以下性质，则称为李代数$\mathbb{g}=(\mathbb{V,F,[,]})$：
1. 封闭性 $\forall X,Y\in\mathbb{V},[X,Y]\in\mathbb{V}$
2. 双线性 $\forall X,Y,Z\in\mathbb{V},a,b\in\mathbb{F},$
  
   $[aX+bY,Z]=a[X,Z]+b[Y,Z],[Z,aX+bY]=a[Z,X]+b[Z,Y]$
3. 自反性 $\forall X\in\mathbb{V},[X,X]=0$
4. 雅可比等价 $\forall X,Y,Z\in\mathbb{V},[X,[Y,Z]]+[Z,[X,Y]]+[Y,[Z,X]]=0$
#### 5.1.5 $\mathfrak{so}(3)$
前文已经对$\mathrm{SO}(3)$和$\mathfrak{so}(3)$之间的关系进行推导，这里我们给出定义
$$
\mathfrak{so}(3)=\{\phi\in\mathbb{R^3},\Phi=\phi^\wedge\in\mathbb{R^{3\times3}}\}
$$
此时的李括号为
$$
[\phi_1,\phi_2]=(\Phi_1\Phi_2-\Phi_2\Phi_1)^\vee
$$
关系表达式由指数映射给定，即
$$
R=\exp(\phi^\wedge)
$$
#### 5.1.6 $\mathfrak{se}(3)$
与$\mathfrak{so}(3)$类似，$\mathfrak{se}(3)$在$\mathbb{R^6}$空间中，定义如下
$$
\mathfrak{se}(3)=\{
\xi=\begin{bmatrix}\rho\\\phi\end{bmatrix}\in\mathbb{R^6},
\rho\in\mathbb{R^3},\phi\in\mathfrak{so}(3),
\xi^\wedge=
\begin{bmatrix}
	\phi^\wedge & \rho \\
	\boldsymbol{0}^T & 0
\end{bmatrix}\in\mathbb{R^{4\times4}}\}
$$
六维向量$\xi$的前三维$\rho$是平移，但是有别于平移矩阵（后面会提），后三维$\phi$是旋转，需要注意$\wedge$符号在这里不表示反对称矩阵，仅是代表从向量到矩阵的变换符号。这里我们给出$\mathfrak{se}(3)$的李括号
$$
[\xi_1,\xi_2]=(\xi_1^\wedge\xi_2^\wedge-\xi_2^\wedge\xi_1^\wedge)^\vee
$$
### 5.2 指数和对数映射
#### 5.1.2 $\mathrm{SO}(3)$指数映射
前面我们已经推导过，对于$\phi=\theta\hat\omega$，有
$$
\exp(\phi^{\wedge}) = I + \sin\theta\hat\omega^{\wedge} + (1-\cos\theta)(\hat\omega^{\wedge})^2
$$
容易发现，对于前文的$\hat\omega$，有
$$
(\hat\omega^{\wedge})^2=\hat\omega\hat\omega^T-I
$$
代入可得
$$
\exp(\theta\hat\omega^\wedge)=
\cos\theta I+(1-\cos\theta)\hat\omega\hat\omega^T+\sin\theta\hat\omega^\wedge
$$
这里我们发现这个式子和$Rodrigues$ $Formula$一致，表明$\mathfrak{so}(3)$就是由全部旋转向量组成的空间。

类似地，我们可以构造一个对数映射
$$
\phi = \ln(R)^\vee=(\sum_{n=0}^\infty\frac{(-1)^n}{n+1}(R-I)^{(n+1)})^\vee
$$
当然我们没有必要用这样的方法找旋转矩阵对应的李代数，可以使用矩阵的迹的性质分别求解角度和转轴，这一点我们在前面有介绍过。

为保证一一对应的关系，我们把旋转角固定在$\pm\pi$之间。
#### 5.2.3 $\mathrm{SE}(3)$指数映射
根据矩阵指数的定义，可得
$$
\begin{aligned}
\exp(\xi^\wedge)&=\exp(\begin{bmatrix}
	\phi^\wedge & \rho \\ \boldsymbol{0}^T & 0
\end{bmatrix}) \\
&=I + \begin{bmatrix}
	\phi^\wedge & \rho \\ \boldsymbol{0}^T & 0
\end{bmatrix} + \frac{1}{2!}
\begin{bmatrix}
	(\phi^\wedge)^2 & (\phi^\wedge)\rho \\ \boldsymbol{0}^T & 0
\end{bmatrix} + \cdots \\
&=\begin{bmatrix}
\sum_{n=0}^{\infty}\frac{1}{n!}(\phi^\wedge)^n & \sum_{n=0}^{\infty}\frac{1}{(n+1)!}(\phi^\wedge)^n\rho \\
	\boldsymbol{0}^T & 1
\end{bmatrix} \\
&= \begin{bmatrix}
	R & J\rho \\ \boldsymbol{0}^T & 1
\end{bmatrix} = T
\end{aligned}
$$
其中
$$
\begin{aligned}
J&=\sum_{n=0}^{\infty}\frac{1}{(n+1)!}(\phi^\wedge)^n \\
&= I + \frac{1}{2!}\theta\boldsymbol{a}^\wedge+\frac{1}{3!}\theta^2(\boldsymbol{a}^\wedge)^2+\cdots \\
&=I+\frac{1}{\theta}(\frac{1}{2!}\theta^2 - \frac{1}{4!}\theta^4+\frac{1}{6!}\theta^6-\cdots)(\boldsymbol{a}^\wedge)+\frac{1}{\theta}(\frac{1}{3!}\theta^3 - \frac{1}{5!}\theta^5+\frac{1}{7!}\theta^7-\cdots)(\boldsymbol{a}^\wedge)^2 \\
&=I+\frac{1}{\theta}(1-\cos\theta)(\boldsymbol{a}^\wedge)+\frac{\theta-\sin\theta}{\theta}(\boldsymbol{a}^\wedge)^2 \\
&=\frac{\sin\theta}{\theta}I+(1-\frac{\sin\theta}{\theta})\boldsymbol{a}\boldsymbol{a}^T+\frac{1-\cos\theta}{\theta}\boldsymbol{a}^\wedge
\end{aligned}
$$
由此可以看出，平移向量$\boldsymbol{t}$满足$\boldsymbol{t}=J\rho$

与$\mathrm{SO}(3)$类似，我们也可以给出对数映射，这里我们可以根据上面的公式，使用如下方法进行定义：
$$
\begin{bmatrix}
	R & \boldsymbol{t} \\ \boldsymbol{0}^T & 1
\end{bmatrix} \in \mathrm{SE}(3) \\
\log(\begin{bmatrix}
	R & \boldsymbol{t} \\ \boldsymbol{0}^T & 1
\end{bmatrix}) = \begin{bmatrix}
	\log(R) & J^{-1}\boldsymbol{t} \\ \boldsymbol{0}^T & 0
\end{bmatrix} \\
\theta = \arccos\frac{tr(R)-1}{2} \\
J^{-1}=I-\frac{1}{2}\Phi+\frac{1}{\theta^2}(1-\frac{\theta \sin\theta}{2(1-\cos\theta)})\Phi^2 \\
\Phi=\phi^\wedge=\theta\boldsymbol{a}^\wedge,\|\boldsymbol{a}\|=1
$$
如果感觉上述推导不是很理解可以看[这里](https://zhuanlan.zhihu.com/p/88771394)

下图展示了$\mathrm{SO}(3)$和$\mathfrak{so}(3)$，$\mathrm{SE}(3)$和$\mathfrak{se}(3)$之间的相互关系

![指数映射和对数映射](.\assets\指数映射和对数映射.png)

### 5.3 李代数求导与扰动模型

#### 5.3.1 BCH近似

还是以$\mathrm{SO}(3)$为例子，当我们在$\mathrm{SO}(3)$完成两个矩阵乘法的时候，李代数$\mathfrak{so}(3)$发生了什么改变？当两个李代数$\mathfrak{so}(3)$做李代数加法的时候在$\mathrm{SO}(3)$上面是否是对应两个矩阵的乘积？即
$$
\exp(\phi_1^\wedge)\exp(\phi_2^\wedge) \stackrel{?}{=}\exp((\phi_1+\phi_2)^\wedge)
$$
容易验证，这个等式是不恒成立的。描述两个李代数的指数映射的完整形式由$Baker-Campbell-Hausdorff(BCH)$公式给出，即
$$
\ln(\exp(A)\exp(B))=A+B+\frac{1}{2}[A,B]+\frac{1}{12}[A,[A,B]]-\frac{1}{12}[B,[A,B]]+\cdots
$$
上述公式的推导可以参考[这篇](https://zhuanlan.zhihu.com/p/610054972)文章，这里我们考虑线性近似表达：

对于$\mathrm{SO}(3)$，有
$$
\ln(\exp(\phi_1^\wedge)\exp(\phi_2^\wedge))^\vee \approx \begin{cases}
	J_l(\phi_2)^{-1}\phi_1+\phi_2 & \phi_1 \approx 0 \\
	J_r(\phi_1)^{-1}\phi_2+\phi_1 & \phi_2 \approx 0
\end{cases} \\
J_l = \frac{\sin\theta}{\theta}I+(1-\frac{\sin\theta}{\theta})\boldsymbol{a}\boldsymbol{a}^T+\frac{1-\cos\theta}{\theta}\boldsymbol{a}^\wedge \\
J_l^{-1} = \frac{\theta}{2}\cot\frac{\theta}{2}I+(1-\frac{\theta}{2}\cot\frac{\theta}{2})\boldsymbol{a}\boldsymbol{a}^T-\frac{\theta}{2}\boldsymbol{a}^\wedge \\
J_r(\phi)=J_l(-\phi)
$$

于是这里就有左乘近似和右乘近似两种情况，后文如果没有特别提及，都使用左乘近似。

这时我们重新叙述$BCH$近似的意义。假设对某个旋转$R$，对应的李代数为$\phi$，我们给他左乘一个微小旋转$\Delta R$，得到的结果是$\Delta R\cdot R$。在李代数上，根据$BCH$近似，结果为$J_l^{-1}(\phi)\Delta\phi+\phi$。合并，可得
$$
\exp(\Delta\phi^\wedge)\exp(\phi^\wedge)=\exp((J_l^{-1}(\phi)\Delta\phi+\phi)^\wedge)
$$

在$\mathfrak{so}(3)$进行加法，可以近似为$\mathrm{SO}(3)$上带左右雅可比的乘法
$$
\exp((\phi+\Delta\phi)^\wedge)=\exp((J_1\Delta\phi)^\wedge)=\exp(\phi^\wedge)\exp((J_r\Delta\phi)^\wedge)
$$
对于$\mathrm{SE}(3)$，也有类似的近似手段：
$$
\exp(\Delta\xi^\wedge)\exp(\xi^\wedge)=\exp((\mathcal{J}_l^{-1}(\xi)\Delta\xi+\xi)^\wedge) \\
\exp(\xi^\wedge)\exp(\Delta\xi^\wedge)=\exp((\mathcal{J}_r^{-1}(\xi)\Delta\xi+\xi)^\wedge)
$$
由于在计算的时候并没有用到$\mathcal{J}_l$和$\mathcal{J}_r$，故这里不写出表达式。

#### 5.3.2 $\mathrm{SO}(3)$中李代数求导
我们想要估计一个相机的位置和姿态，首先就会自然想到要使用$\mathrm{SO}(3)$的旋转矩阵和$\mathrm{SE}(3)$的变换矩阵进行描述。那么要想具体去描述位置和姿态的变化就想到要建立与位姿相关的函数，然后借助这个函数的导数来调整对于位姿的估计值。然而，$\mathrm{SO}(3),\mathrm{SE}(3)$上面并没有良好定义的加法，$T$自身又有比较多的约束；同时$\mathfrak{so}(3),\mathfrak{se}(3)$李代数由向量组成，具有良好的加法运算，因此有以下两种求解方法：
1. 用李代数表示姿态，然后根据李代数加法对李代数求导。
2. 对李群左乘或右乘微小扰动，然后对扰动求导。
#### 5.3.3 李代数求导
假设我们对空间点$\boldsymbol{p}$进行旋转，得到$R\boldsymbol{p}$。现在想要计算旋转之后的坐标相对于旋转的导数，先非正式地记为
$$ \frac{\partial(R\boldsymbol{p})}{\partial R} $$
由于$\mathrm{SO}(3)$没有加法，所以这个导数不能按照导数的定义进行计算。设$R$对应的李代数是$\phi$，则转而计算
$$ \frac{\partial \exp(\phi^\wedge)\boldsymbol{p}}{\partial\phi} $$
根据导数定义，有
$$
\begin{aligned}
\frac{\partial \exp(\phi^\wedge)\boldsymbol{p}}{\partial\phi} &= 
\lim_{\delta\phi \to 0}\frac{\exp((\phi+\delta\phi)^\wedge)\boldsymbol{p}-\exp(\phi^\wedge)\boldsymbol{p}}{\delta\phi} \\ &=
\lim_{\delta\phi \to 0}\frac{\exp((J_l\delta\phi)^\wedge)\exp(\phi^\wedge)\boldsymbol{p}-\exp(\phi^\wedge)\boldsymbol{p}}{\delta\phi} \\ &=
\lim_{\delta\phi \to 0}\frac{(I+(J_l\delta\phi)^\wedge)\exp(\phi^\wedge)\boldsymbol{p}-\exp(\phi^\wedge)\boldsymbol{p}}{\delta\phi} \\ &=
\lim_{\delta\phi \to 0}\frac{(J_l\delta\phi)^\wedge \exp(\phi^\wedge)\boldsymbol{p}}{\delta\phi} \\ &=
\lim_{\delta\phi \to 0}\frac{-(\exp(\phi)^\wedge\boldsymbol{p})^\wedge J_l\delta\phi}{\delta\phi} = -(R\boldsymbol{p})^\wedge J_l
\end{aligned}
$$
下面对上述推导做简要解释：第二行使用前面提到的$BCH$近似，第三行为泰勒展开舍去高阶项的近似（取极限可以写等号），四五行将反对称符号看作叉积，交换之后变号，就可以得到旋转后的点相对于李代数的导数
$$\frac{\partial \exp(\phi^\wedge)\boldsymbol{p}}{\partial\phi}=-(R\boldsymbol{p})^\wedge J_l$$
由于我们不希望计算形式比较复杂的$J_l$，所以后面将讲解更简单的导数计算方式。
#### 5.3.4 左乘扰动模型
我们对$R$进行一次扰动$\Delta R$，看结果相对于扰动的变化率，这个扰动可以乘在左边或者右边，值得注意的是，**左扰动和右扰动对应的结果会有差异**。令左扰动$\Delta R$对应的李代数为$\varphi$，对$\varphi$求导，即
$$
\begin{aligned}
\frac{\partial R\boldsymbol{p}}{\partial\varphi} &=
\lim_{\varphi\to0}\frac{\exp(\varphi^\wedge)\exp(\phi^\wedge)\boldsymbol{p}-\exp(\phi^\wedge)\boldsymbol{p}}{\delta\phi} \\ &=
\lim_{\varphi\to0}\frac{(I+\varphi^\wedge)\exp(\phi^\wedge)\boldsymbol{p}-\exp(\phi^\wedge)\boldsymbol{p}}{\delta\phi} \\ &=
\lim_{\varphi\to0}\frac{\varphi^\wedge R\boldsymbol{p}}{\varphi} =
\lim_{\varphi\to0}\frac{-(R\boldsymbol{p})^\wedge\varphi}{\varphi} = -(R\boldsymbol{p})^\wedge
\end{aligned}
$$
可以发现，相比于直接对李代数求导不需要计算$J_l$，这使得扰动模型更加实用。
#### 5.3.5 $\mathrm{SE}(3)$中李代数求导（左扰动模型）
假设某空间点$\boldsymbol{p}$经过一次变换$T$（对应的李代数为$\xi$），得到$T\boldsymbol{p}$。给$T$左乘一个扰动$\Delta T=\exp(\delta\xi^\wedge)$，设扰动项李代数$\delta\xi=[\delta\rho,\delta\phi]^T$，则有
$$
\begin{aligned}
\frac{\partial(T\boldsymbol{p})}{\partial(\delta\xi)}&=	\lim_{\delta\xi\to0}\frac{\exp(\delta\xi^\wedge)\exp(\xi^\wedge)\boldsymbol{p}-\exp(\xi^\wedge)\boldsymbol{p}}{\delta\xi} \\ &=
\lim_{\delta\xi\to0}\frac{(I+\delta\xi^\wedge)\exp(\xi^\wedge)\boldsymbol{p}-\exp(\xi^\wedge)\boldsymbol{p}}{\delta\xi} \\ &=
\lim_{\delta\xi\to0}\frac{\delta\xi^\wedge \exp(\xi^\wedge)\boldsymbol{p}}{\delta\xi} \\ &=
\lim_{\delta\xi\to0}\frac{\begin{bmatrix}\delta\phi^\wedge & \delta\rho \\ \boldsymbol{0}^T & 0\end{bmatrix}\begin{bmatrix}R\boldsymbol{p}+\boldsymbol{t} \\ 1\end{bmatrix}}{\delta\xi} \\ &=
\lim_{\delta\xi\to0}\frac{\begin{bmatrix}\delta\phi^\wedge(R\boldsymbol{p}+\boldsymbol{t})+\delta\rho \\ \boldsymbol{0}^T\end{bmatrix}}{[\delta\rho, \delta\phi]^T} =
\begin{bmatrix}
I & -(R\boldsymbol{p}+\boldsymbol{t})^\wedge \\ \boldsymbol{0}^T & 0^T
\end{bmatrix} \stackrel{def}{=} (T\boldsymbol{p})^\odot
\end{aligned}
$$
我们把最终结果定义为一个算符$^\odot$，可以将一个齐次坐标的空间点变换成一个$4\times6$的矩阵。假设$\boldsymbol{a,b,x,y}$都是列向量，那么这个符号有以下规则：
$$
\frac{d\begin{bmatrix}\boldsymbol{a} \\ \boldsymbol{b}\end{bmatrix}}{d\begin{bmatrix}\boldsymbol{x} \\ \boldsymbol{y}\end{bmatrix}} = 
\begin{bmatrix}
	\frac{d\boldsymbol{a}}{d\boldsymbol{x}} & \frac{d\boldsymbol{a}}{d\boldsymbol{y}} \\
	\frac{d\boldsymbol{b}}{d\boldsymbol{x}} & \frac{d\boldsymbol{b}}{d\boldsymbol{y}}
\end{bmatrix}
$$

对于使用左扰动模型和右扰动模型对$\mathrm{SE}(3)$进行求偏导之间的区别，请参考下面这两篇文章[1](https://blog.csdn.net/qq_45060471/article/details/127705516)，[2](https://www.zhihu.com/question/454486535)（主要是不想写了:sleeping:）


## KF
考查一维线性系统
$$
\begin{align}
&\dot{x} = Ax+Bu+w\\
&z = Hx+v
\end{align}\tag{2.1}
$$
其中$w\sim N(0,Q), v\sim N(0,R), cov(w,v)=0$，离散化，可得
$$
\begin{align}
& x_{k+1} = Ax_{k}+Bu_k+w\\
& z_{k+1}=Hx_{k+1}+v
\end{align}\tag{2.2}
$$
这里$w$称为**过程噪声**,$v$称为**测量噪声**，$z$称为**观测值**
假定初始状态$x_0$是已知均值和方差的随机向量，希望进行最优估计，这里首先进行定义
**最佳估值**$\hat{x}_k=E[x_k|z_0,\cdots,z_k]=E[x_k|Z^k]$
其中$Z^k=\{z_0,\cdots,z_k\}$满足$E[\{x_k-\hat{x}_k\}\{x_k-\hat{x}_k\}^T|Z^k]$最小
**先验估计**$\hat{x}^-_k=\hat{x}_{k+1/k}=E[x_{k+1}|z_0,\cdots,z_k]=E[x_{k+1}|Z^k]$
**后验估计**$\hat{x}^+_k=\hat{x}_{k/k}=E[x_k|z_0,\cdots,z_k]=E[x_k|Z^k]$
这样理解起来显然不太方便，于是可以确定以下关系式：
$$
\hat{x}^-_{k+1}=A\hat{x}^+_k+Bu_k\tag{2.3}
$$
通俗来说，就是先验估计是通过上一次的后验估计和输入通过线性系统得到的计算值。
这时可以确定目标如下：
对于线性状态估计
$$
\hat{x}^+_k=\hat{x}^-_k+K(z_k-H\hat{x}^-_k)\tag{2.4}
$$
使其均方误差最小，即
$$
{\min}\limits_{K_k}\,E[(e^+_k)^2]
$$
其中
$$
e^+_k=x_k-\hat{x}^+_k
$$
（说人话就是希望最优估计值（后验估计）接近实际值）
通过上述方法对线性状态估计进行表示的原因可以理解为通过增益矩阵$K_k$来描述先验与后验之间的误差和测量值与观测值之间的误差的关系。
类似地，有$e^-_k=x_k-\hat{x}_k^-$，这里首先观察$e^-_k$和$e^+_k$之间的关系，可以发现
$$
\begin{align}
 \hat{x}^+_k&=(I-KH)\hat{x}^-_k+Kz_k\\
 &=(I-KH)\hat{x}^-_k+KHx_k+Kv_k\\
 \Rightarrow \hat{x}^+_k-x_k&=(I-KH)(\hat{x}^-_k-x_k)+Kv_k\\
 \Rightarrow e^+_k&=(I-KH)e^-_k-Kv_k \tag{2.5}
\end{align}
$$
引入**先验协方差**$P_k^-=E[(e^-_k)^2]$，**后验协方差**$P^+_k=E[(e^+_k)^2]$，容易发现$P_k=P_k^T$，并注意到$e_k$与$v_k$独立，基于公式2.5，可得
$$
\begin{align}
P^+_k&=E[e^+_ke^{+T}_k]\\
&=E[((I-KH)e^-_k-Kv_k)((I-KH)e^-_k-Kv_k)^T]\\
&=E[(I-KH)e^-_ke^{-T}_k(I-KH)^T]-0-0+E[Kv_kv_k^TK^T]\\
&=(I-KH)P^-_k(I-KH)^T+KRK^T\\
&=P^-_k-KHP^-_k-P^-_kH^TK^T+K(HP^-_kH^T+R)K^T\tag{2.6}
\end{align}
$$
接下来让$P^+_k$对$K$求偏导，即
$$
\begin{align}
&\frac{d\,tr(P^+_k)}{dK}=0 \\
\Leftrightarrow\,\, & 0-(HP_k^-)^T-(HP^{-T}_k)^T+K(HP_k^-H^T+R+(HP^-_kH^T+R)^T)=0 \\
\Rightarrow\,\, & P^-_kH^T=K(R+HP^-_kH^T)\\
\Rightarrow\,\,& K=\frac{P^-_kH^T}{R+HP_k^-H^T}\tag{2.7}
\end{align}
$$
类似地，针对$P^-_k$，通过公式2.3，可以推导得到
$$
\begin{align}
&\hat{x}^-_{k+1}-x_{k+1}=A(\hat{x}^+_k-x^+_k)-w \\
\Leftrightarrow\,\, & e^-_{k+1}=Ae^+_k+w\tag{2.8}
\end{align}
$$
从而有
$$
\begin{align}
P^-_k&=E[(Ae^+_{k-1}+w)(Ae^+_{k-1}+w)^T]\\
&=E(Ae^T_{k-1}e^{+T}_{k-1}A^T)-0-0+Q\\
&=AP^-_{k-1}A^T+Q\tag{2.9}
\end{align}
$$
这样就得到了KF整个流程：
+ 初始化$P,Q,R$
+ 预测步骤：
	1. 状态更新：$x_{k+1}=Ax_k+Bu_k$
	2. 更新协方差矩阵：$P_{k+1}=AP_{k}Q^T+Q$
+ 修正步骤：
	3. 计算卡尔曼增益矩阵：$K=P_kH^T(HP_kH^T+R)^{-1}$
	4. 更新估计值：$x_k=x_k+K(z_k-Hx_k)$
	5. 更新协方差矩阵：$P_{k+1}=(I-KH)P_k-P_kH^TK^T+K(HP_kH^T+R)K^T$
	   通常在更新协方差矩阵这步选择简化形式：$P_{k+1}=(I-KH)P_k$
	
 


## EKF


## ESKF



## IEKF



## IESKF



## MSCKF

