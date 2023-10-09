\section{Brachistochrone}
Consider a particle subject to a uniform gravitationnal field of strength $g$ oriented in the negative $y$ direction of a planar $Oxy$ cartesian coordinate system. The particle starts at point $A$ of coordinates $x_A, y_A$ and ends at point $B$ of coordinates $x_B, y_B$. We are interested in the curve linking $A$ and $B$ minimizing the travel time of the particle. 
%\input{figures/P2_1.pdf}
\subsection{Travel time functionnal}
The first step in constructing the travel time functionnal is the parametrisation of its curve argument. If $x_A = x_B$, the trajectory has to be a straight vertical line. Indeed to minimize the travel time, gravity has to be exploited as much as possible and and disformation of the vertical path will reduce components of the gravitationnal acceleration along the trajectory (increasing the travel time). We therefore restrict the analysis to $x_B > x_A$.


Suppose the trajectory is parametrised by $x$ so that it takes the form of a function $y(x)$ (this supposes that $x$ is always increasing from $x_A$ to $x_B$). The length element at $x$ is $\text{d} l = \sqrt{\dot{x}^2 + \dot{y}(x)^2} \text{d} x = \sqrt{1 + \dot{y}(x)^2} \text{d} x$ with $\dot{()} = \frac{d}{dx}$. By a mass independant analogue of the conservation of energy $0 = \frac{1}{2} v^2-g y$ (with "potential" energy zero at $y_B$) where $v$ is the velocity, we have $v = \sqrt{2gy}$. For motion along an element $\text{d}l$ of the curve, the time elapsed is $\text{d}t = \text{d}l/\sqrt{2gy}$ leading to te following travel time functionnal 
\begin{align*}
    T[y(x)] = \int_{x_A}^{x_B} \text{d} x \dfrac{\sqrt{1 + \dot{y}(x)^2}}{\sqrt{2gy}}. 
\end{align*}



%\input{figures/P1_1.pdf}
\subsection{Conserved quantity}
Since the integrand (effective lagrangian $L = \sqrt{1 + \dot{y}(x)^2}/{\sqrt{2gy}}$) of the previous functionnal is independant of $x$, we have the following invariant quantity
\begin{align*}
    Q = \dfrac{\partial L}{\partial \dot{y}} \dot{y} - L = \dfrac{2\dot{y}(x)^2}{\sqrt{1 + \dot{y}(x)^2}\sqrt{2gy}} - \dfrac{\sqrt{1 + \dot{y}(x)^2}}{\sqrt{2gy}}.
\end{align*}
We can now write 
\begin{align*}
    Q^2 = \left(\dfrac{2\dot{y}(x)^2}{\sqrt{1 + \dot{y}(x)^2}\sqrt{2gy}} - \dfrac{\sqrt{1 + \dot{y}(x)^2}}{\sqrt{2gy}}\right)^2 = \dfrac{1}{2gy} \left(\dfrac{4\dot{y}(x)^4}{1 + \dot{y}(x)^2} - 1 + \dot{y}(x)^2-4 \dot{y}(x)^2\right)
\end{align*}

%gy Q^2(1 + \dot{x}(y)^2) = 2\dot{x}(y)^2 &\implies \dot{x}(y) = \sqrt{\dfrac{gy Q^2}{2-gy Q^2}}\\ 
%&\implies \dot{x}(y) = \int \text{d} y \sqrt{\dfrac{gy Q^2}{2-gy Q^2}}


\begin{align*}
    Q =\dfrac{2\sqrt{2g y_B}}{\sqrt{1 + 2g y_B}\sqrt{2gy_B}}=\dfrac{2}{\sqrt{1 + 2g y_B}}
\end{align*}
so that the radius of the cycloid motion is 
\begin{align*}
    a = \dfrac{1 + 2g y_B}{2g}.
\end{align*}
Finally, using the expression for $y$, 
\begin{align*}
    y_B = \dfrac{1 + 2g y_B}{4g}(1 + \cos(\theta_B)) \iff \cos(\theta_B) = \dfrac{4gy_B - 1 - 2g y_B}{1 + 2g y_B} = \dfrac{2gy_B - 1}{2g y_B + 1}.
\end{align*} 

% The new paramatrisation provided by integration makes y retake its values

We want to relate the value of $Q$ to $y_B$, $x_B$. By "conservation of the energy" $0 = \frac{1}{2}\dot{x}(y_B)^2 - g y_B \implies \dot{x}(y_B) = \pm \sqrt{2g y_B}$. Because $x_B > x_A$, the particle will reach $x_B$ with positive velocity ($\dot{x}(y_B) = \sqrt{2g y_B}$). At $x_B, y_B$, the quantity $Q$ can be expressed as 

Going further, setting $y_A = 0$ we have an initial angle $\theta_A = -\pi$ (setting the expression for $y$ to $0$) and $C =  0$ implies a reference $x_A = \dfrac{\pi}{gQ^2}$.


Consider the case $E=0$ and $\omega = \omega_0$ where the particle's history is a circular motion at fixed height. The energy $E_0$ is then the sum of Kinetic $\frac12 m (\rho \omega)^2$ and potential $m g \alpha \rho^2 = \frac12 m \omega_0^2 \rho^2  = \frac12 m \omega^2 \rho^2$ so the missing half of motor energy is due to maintaining the particle's height in this case. 

where we used the fact $\mathbf{u} \times \mathbf{n}/u$ is a unit vector in the direction orthogonal to $\mathbf{u}$ and orthogonal to $\mathbf{n}$. Since $\mathbf{u}$, $\mathbf{n}$ are orthogonal to $\mathbf{v}'_\perp$, $\mathbf{u} \times \mathbf{n}/u$ must be parellel to $\mathbf{v}'_\perp$. More precisely, the right hand rule for the vector product ensures that $|\mathbf{v}_\perp|\mathbf{u} \times (\mathbf{u} \times \mathbf{m})/u^2 = |\mathbf{v}_\perp'|\mathbf{m} = \mathbf{v}'_\perp$where we used the fact $\mathbf{u} \times \mathbf{n}/u$ is a unit vector in the direction orthogonal to $\mathbf{u}$ and orthogonal to $\mathbf{n}$. Since $\mathbf{u}$, $\mathbf{n}$ are orthogonal to $\mathbf{v}'_\perp$, $\mathbf{u} \times \mathbf{n}/u$ must be parellel to $\mathbf{v}'_\perp$. More precisely, the right hand rule for the vector product ensures that $|\mathbf{v}_\perp|\mathbf{u} \times (\mathbf{u} \times \mathbf{m})/u^2 = |\mathbf{v}_\perp'|\mathbf{m} = \mathbf{v}'_\perp$


\begin{align*}
    \dfrac{\mathbf{u}}{u} \times \mathbf{v}' = \dfrac{\mathbf{u} \times (\mathbf{v}'_\parallel + \mathbf{v}'_\perp)}{u} = \dfrac{\mathbf{u}}{u} \times \mathbf{v}'_\perp.
\end{align*}
Taking the vector product with $\mathbf{u}$ again leads to 

= \dfrac{\frac{\delta x_\parallel}{\delta t} \frac{\mathbf{u}}{u} - \mathbf{u}'}{1 - u \frac{\delta x_\parallel}{\delta t}/c^2}

\dfrac{\textbf{v}' + \mathbf{u}}{1 + \mathbf{u} \cdot \mathbf{v}'/c^2} + (\gamma(u)^{-1}-1)\dfrac{\textbf{v}'_\perp}{1 + \textbf{u}' \cdot \textbf{v}'/c^2}


\mathbf{E}\times\mathbf{B} = (a_0 \mathbf{k} - \mathbf{a} k_0)  \cos\left(k_\mu x^\mu + \pi/2\right) \times \left(\mathbf{k} \times \mathbf{a}\cos\left(k_\mu x^\mu + \pi/2\right) \right) = (- \mathbf{a} k_0)   \times \left(\mathbf{k} \times \mathbf{a}\right) \cos^2\left(k_\mu x^\mu + \pi/2\right)

\mathbf{S} = -\dot{\mathbf{A}} \times (\mathbf{n} \times \dot{\mathbf{A}})  = (\dot{\mathbf{A}} \cdot \dot{\mathbf{A}}) \mathbf{n} - (\dot{\mathbf{A}} \cdot \mathbf{n}) \mathbf{\dot{\mathbf{A}}} =

\dfrac{1}{|\mathbf{r}| - \mathbf{r}'\cdot\dfrac{\mathbf{r}}{|\mathbf{r}|} + O(|\mathbf{r}'|^2/|\mathbf{r}|^2)} =

\dfrac{1}{|\mathbf{r}|} \int \text{d}^3r' \left[\nabla \cdot (\mathbf{r'} j_i(t', \mathbf{r}'))- \mathbf{r'} \cdot \nabla j_i(t', \mathbf{r}')\right]

= -\dfrac{1}{|\mathbf{r}|^2}\dfrac{\partial \dot{\mathbf{d}}}{\partial t} \times (\underbrace{\nabla \times \dot{\mathbf{d}}}_0) -\dfrac{1}{|\mathbf{r}|}\dfrac{\partial \dot{\mathbf{d}}}{\partial t} \times (\nabla \dfrac{1}{|\mathbf{r}|} \times \dot{\mathbf{d}})  = \dfrac{1}{|\mathbf{r}|^4}\dfrac{\partial \dot{\mathbf{d}}}{\partial t} \times ( \mathbf{r} \times \dot{\mathbf{d}}) = \dfrac{1}{|\mathbf{r}|^4}\left(\left[\dfrac{\partial \dot{\mathbf{d}}}{\partial t} \cdot \dot{\mathbf{d}}\right] \mathbf{r} - \left[\dfrac{\partial \dot{\mathbf{d}}}{\partial t} \cdot  \mathbf{r} \right] \dot{\mathbf{d}}\right)

\dfrac{1}{|\mathbf{r}|^}\left(\left[\dfrac{\partial \dot{\mathbf{d}}}{\partial t} \cdot \dot{\mathbf{d}}\right] \mathbf{r} - \left[\dfrac{\partial \dot{\mathbf{d}}}{\partial t} \cdot  \mathbf{r} \right] \dot{\mathbf{d}}\right)

= \dfrac{1}{|\mathbf{r}|^4}\dfrac{\partial \dot{\mathbf{d}}}{\partial t} \times ( \mathbf{r} \times \dot{\mathbf{d}}) = \dfrac{1}{|\mathbf{r}|^4}\left(\left[\dfrac{\partial \dot{\mathbf{d}}}{\partial t} \cdot \dot{\mathbf{d}}\right] \mathbf{r} - \left[\dfrac{\partial \dot{\mathbf{d}}}{\partial t} \cdot  \mathbf{r} \right] \dot{\mathbf{d}}\right) 


-\dfrac{1}{|\mathbf{r}|^2}\nabla_{\dot{\mathbf{d}}}\left[\dfrac{\partial \dot{\mathbf{d}}}{\partial t} \cdot \dot{\mathbf{d}}\right] + \dfrac{1}{|\mathbf{r}|^2}\dfrac{\partial \dot{\mathbf{d}}}{\partial t} \cdot \nabla \times (\nabla \times \dot{\mathbf{d}})  + \left[\sim \dfrac{1}{|\mathbf{r}|^3}\right] 

\begin{align*}
    \mathbf{S} = -\dfrac{|\mathbf{\ddot{d}}|}{|\mathbf{r}|^2} \mathbf{\ddot{d}} \times \mathbf{u}\cos(\theta) + \left[\sim \dfrac{1}{|\mathbf{r}|^3}\right] = \dfrac{|\mathbf{\ddot{d}}|^2}{|\mathbf{r}|^2}  \dfrac{\mathbf{r}}{|\mathbf{r}|} + \left[\sim \dfrac{1}{|\mathbf{r}|^3}\right].
\end{align*}
with $\mathbf{u} = (\mathbf{r} \times \mathbf{\ddot{d}})/|\mathbf{r} \times \mathbf{\ddot{d}}|$  and $\mathbf{\ddot{d}} \times \mathbf{u} =  $