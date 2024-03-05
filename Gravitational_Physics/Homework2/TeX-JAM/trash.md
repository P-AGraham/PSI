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


To evaluate the second integral, we notice that the ground sate of the harmonic oscillator is an even function and has an odd first $x'$ derivative making the integrand globally odd. The bounds of integration go from $-\infty$ to $\infty$ (symmetric interval around $x' = 0$) and this integral vanishes.

\begin{align*}
A &= C\left(\prod_{n = 1}^{N-1}\int_{-\infty}^{\infty} \text{d} q_n \right) e^{\frac{i}{\hbar}\sum_{n=0}^{N-1} \frac12 m \left(\frac{q_{n+1}-q_{n}}{\Delta t}\right)^2  \Delta t} \\
&= C\left(\prod_{n = 1}^{N-1}\int_{-\infty}^{\infty} \text{d} q_n e^{\frac{i}{\hbar} \frac12 m \left(\frac{q_{n+1}-q_{n}}{\Delta t}\right)^2 \Delta t}\right) \times \int_{-\infty}^{\infty} \text{d} q_1 e^{\frac{i}{\hbar} \frac12 m \left(\frac{q_{1}-q_{0}}{\Delta t}\right)^2 \Delta t}\\
&= C\left(\prod_{n = 1}^{N-1}\int_{q_{n+1}+\infty}^{q_{n+1}-\infty} -\text{d} u_n e^{\frac{im}{2\Delta t \hbar}  u_n^2}\right), \quad u_n = q_{n+1}-q_{n}\\
&= C\left(\prod_{n = 1}^{N-1}\int_{-\infty}^{\infty} \text{d} u_n e^{\frac{i m}{2\Delta t \hbar} u_n^2}\right)\\
&= C\left(\prod_{n = 1}^{N-1} \left(\dfrac{2\pi \Delta t \hbar}{i m}\right)^{\frac12}\right)
\end{align*}


\begin{align*}
    &\int_{-\infty}^{\infty} \text{d} q_{N-1} \exp\left(\frac{im}{2\hbar} \left(\frac{q_{N}-q_{N-1}}{n_1\Delta t}\right)^2  n_1\Delta t + \frac{im}{2\hbar} \left(\frac{q_{N-1}-q_{N-2}}{\Delta t}\right)^2  \Delta t\right)\\
    =& \int_{-\infty}^{\infty} \text{d} q_{N-1} \exp\left(\frac{im}{\hbar n_1\Delta t} \left(\left(\frac{q_{N-1}}{n_1}\right)^2 - \left(\frac{q_{N}}{n_1^2} - q_{N-2}\right)q_{N-1}\right) \right)\exp\left(\frac{im}{2\hbar\Delta t} \left(q_{N-2}^2 + \left(\dfrac{q_{N}}{n_1}\right)^2\right)\right)\\
    =& \int_{-\infty}^{\infty} \text{d} q_{N-1} \exp\left(\frac{im}{\hbar\Delta t} \left(\left(\frac{q_{N-1}}{n_1}\right)^2 - \left(\frac{q_N}{n_1^2} - q_{N-2}\right)q_{N-1} +  \dfrac{1}{4}\left(\frac{q_{N}}{n_1^2} - q_{N-2}\right)^2\right) \right)\exp\left(\frac{im}{2\hbar\Delta t} \left(q_{N-2}^2 + \left(\frac{q_{N}}{n_1^2}\right)^2 - \dfrac{1}{2}\left(\dfrac{q_{N}}{n_1^2} - q_{N-2}\right)^2\right)\right)\\
    =& \left(\dfrac{\hbar\pi \Delta t}{m i}\right)^{1/2}\exp\left(\frac{m i}{4\hbar\Delta t} \left(q_{N-2}^2 +  \left(\frac{q_N}{n_1}\right)^2 - 2q_{N}q_{N-2}\right)\right) = \left(\dfrac{\hbar\pi \Delta t}{m i}\right)^{1/2}\exp\left(\frac{m i}{2\hbar} \left(\dfrac{q_{N}-q_{N-2}}{2\Delta t}\right)^2 2\Delta t\right) 
\end{align*}


\left(\prod_{n = 1}^{N-2}\int_{-\infty}^{\infty} \text{d} q_n \right) \exp\left(\frac{i}{\hbar}\sum_{n=0}^{N-3} \frac12 m \left(\frac{q_{n+1}-q_{n}}{\Delta t}\right)^2  \Delta t\right)  \left(\dfrac{\hbar\pi \Delta t r}{m i(r+1)}\right)^{1/2}\exp\left(\frac{m i}{2\hbar} \left(\dfrac{q_{M}-q_{M-2}}{(r+1)\Delta t}\right)^2 (r+1)\Delta t\right) 


&=C\left(\prod_{n = 1}^{N-r-1}\int_{-\infty}^{\infty} \text{d} q_n \right) \exp\left(\frac{i}{\hbar}\sum_{n=0}^{N-3} \frac12 m \left(\frac{q_{n+1}-q_{n}}{\Delta t}\right)^2  \Delta t\right)  \left(\dfrac{\hbar\pi \Delta t r}{m i(r+1)}\right)^{1/2}\exp\left(\frac{m i}{2\hbar} \left(\dfrac{q_{M}-q_{M-2}}{(r+1)\Delta t}\right)^2 (r+1)\Delta t\right)

+
      \begin{tikzpicture}[baseline=-\the\dimexpr\fontdimen22\textfont2\relax]
        \begin{feynman}
          % External legs
          \vertex[label={left:$\mathbf{p}$}] (a) at (-1, 0);
          \vertex[label={below right:$\mathbf{k}_1$}] (b) at (1, -1);
          \vertex[label={above right:$\mathbf{k}_2$}] (c) at (1, 1);
          

          \vertex (e) at (0.3, 0.3);
          \vertex (f) at (0.6, 0.6);

          \fill (a) circle (2pt);
          \fill (b) circle (2pt);
          \fill (c) circle (2pt);
          \fill (e) circle (2pt);
          \fill (f) circle (2pt);
    
          \diagram* {
            (a) -- [dashed] (v),
            (v) -- (b),
            (v) -- (e),
            (e) -- [dashed, half left] (f) -- [half left] (e), 
            (f) -- (c),
          };
        \end{feynman}
      \end{tikzpicture}
      +
      \begin{tikzpicture}[baseline=-\the\dimexpr\fontdimen22\textfont2\relax]
        \begin{feynman}
          % External legs
          \vertex[label={left:$\mathbf{p}$}] (a) at (-1, 0);
          \vertex[label={below right:$\mathbf{k}_1$}] (b) at (1, -1);
          \vertex[label={above right:$\mathbf{k}_2$}] (c) at (1, 1);
          

          \vertex (e) at (0.3, -0.3);
          \vertex (f) at (0.6, -0.6);

          \fill (a) circle (2pt);
          \fill (b) circle (2pt);
          \fill (c) circle (2pt);
          \fill (e) circle (2pt);
          \fill (f) circle (2pt);
    
          \diagram* {
            (a) -- [dashed] (v),
            (v) -- (c),
            (v) -- (e),
            (e) -- [dashed, half left] (f) -- [half left] (e), 
            (f) -- (b),
          };
        \end{feynman}
      \end{tikzpicture}


\Lambda_{\mu}{}^{\rho}  \omega_{\rho\sigma}  (\Lambda^{-1})^{\sigma}{}_{\nu} = 


We can also write $\epsilon_{ij}{}^k J^{ij} =  J^k$ where $-\epsilon_{ij}{}^k = \epsilon_{ijk}$ is the levi-civita tensor with one up index in our mostly minus signature.

\epsilon_{ij}{}^{k} \epsilon_{lm}{}^{n} [J^{ij}, J^{lm}] &= -i \epsilon_{ij}{}^{k} \epsilon_{lm}{}^{n} \left(\eta^{m j} J^{l i} -\eta^{m i} J^{lj} + \eta^{jl} J^{im} - \eta^{il} J^{jm}\right)\\
    &= -i \left(- \epsilon_{i}{}^{mk} \epsilon_{lm}{}^{n} J^{l i} + \epsilon^{m}{}_{j}{}^{k} \epsilon_{lm}{}^{n} J^{lj} - \epsilon_{i}{}^{lk} \epsilon_{lm}{}^{n} J^{im} + \epsilon^{l}{}_{j}{}^{k} \epsilon_{lm}{}^{n} J^{jm}\right)\\
    &= -i \left(- \epsilon_{i}{}^{km} \epsilon_{l}{}^{nm} J^{l i} + \epsilon^{k}{}_{j}{}^{m} \epsilon_{l}{}^{n}{}_m J^{lj} - \epsilon_{i}{}^{lk} \epsilon_{lm}{}^{n} J^{im} + \epsilon^{l}{}_{j}{}^{k} \epsilon_{lm}{}^{n} J^{jm}\right)


The only non-vanishing one-loop contribution contracts fields sharing the same type and the effective interaction term is of the form $\bar{\psi}_a \psi^a \bar{\psi}_a \psi^a$ (with no summation implied). This implies that the incomming and outgoing fields are contracted with the same color forcing them to share this color. Therefore, the one-loop correction only consitutes a self-interaction correction. 

Supposing there exists a change of variables from $z, z'$ to $y, y'$  bringing $\delta(z-z') (i\slashed{\partial}_{z} - \sigma)$ to $f(y) \delta(y-y')$ (a diagonalizing change of variables), the determinant can be expanded as 

\\
   &\sim \frac{1}{(4\pi)^2}\left[\frac{1}{u}-\frac{m_R^2}{2u^2}\right]_{m_R^2}^{(\Lambda/p)^2 + m_R^2} = \frac{1}{(4\pi)^2}\left[\frac{1}{(\Lambda/p)^2 + m_R^2}-\frac{m_R^2}{2((\Lambda/p)^2 + m_R^2)^2}\right]_{m_R^2}^{(\Lambda/p)^2 + m_R^2} 

   Since the derivative in $m_R$ is finite, it is $\sim 0$ from the point of view of UV divergence and the integral divergence behaves as a constant (the same constant as for $m_R = 0$ found in (a)).

   Decomposing $B_1$ in a term $B_{1, 0}$ (independant of $m_R$) and a term $m_R^2 B_{1, 1}$ (proportionnal to $m_R^2$), we get $B_{1, 0} = -\frac{g_R}{2} \frac{1}{(4\pi)^2} \Lambda^2$ and $B_{1, 1} = \frac{g_R}{2} \frac{1}{(4\pi)^2}\ln\left(\frac{\Lambda^2}{m_R^2}\right)$


Expanding $\underline{\omega}^b$, we have
  \begin{align*}
    g_{\mu \nu} \underline{d} x^\mu \underline{d} x^\nu = \eta_{ab} c_\mu^a c_\nu^b \underline{d} x^\mu \underline{d} x^\nu &\iff g_{\mu \nu} = \eta_{ab}  c_\mu^a c_\nu^b \quad \text{(linear independance on the basis of two-forms)}\\
    &\iff 
    \begin{cases}
      -c_\mu^0 c_\nu^0 = g_{\mu\nu} \impliedby c_\mu^0 = (\sqrt{-g_{00}}, 0, 0, 0) = (1, 0, 0, 0) \\
      c_\mu^1 c_\nu^1 = g_{\mu\nu} \impliedby c_\mu^1 = (0, \sqrt{-g_{11}}, 0, 0) = (1, 0, 0, 0) \\
      c_\mu^2 c_\nu^2 = g_{\mu\nu} \impliedby c_\mu^1 = (0, \sqrt{-g_{11}}, 0, 0) = (1, 0, 0, 0) 
    \end{cases}
  \end{align*}



Decomposing the connection one-form as $\underline{\theta}^a{ }_b = c^a_{db} \underline{d}\omega^d$, we have 
  \begin{align*}
    \begin{cases}
      0\\
     \frac{a'(t)}{a(t)} \underline{\omega}^0 \wedge \underline{\omega}^1 + \frac{1}{a(t)r}\sqrt{1-k r^2}\underline{\omega}^3 \wedge \underline{\omega}^1\\
     \frac{a'(t)}{a(t)} \underline{\omega}^0 \wedge \underline{\omega}^2 +  \frac{1}{a(t)r}\sqrt{1-k r^2} \underline{\omega}^3 \wedge \underline{\omega}^2 + \frac{1}{a(t)r}\cot \theta \underline{\omega}^1 \wedge \underline{ \omega}^2\\
     \frac{a'(t)}{a(t)}\underline{\omega}^0 \wedge \underline{\omega}^3
    \end{cases}
    =
    \begin{cases}
      \underline{\theta}^0{ }_b \wedge \underline{\omega}^b \\
      \underline{\theta}^1{ }_b \wedge \underline{\omega}^b \\
      \underline{\theta}^2{ }_b \wedge \underline{\omega}^b \\
      \underline{\theta}^3{ }_b \wedge \underline{\omega}^b 
    \end{cases}
    \iff 
    \begin{cases}
      \underline{\theta}^0{ }_b \propto \omega^b \\
      \underline{\theta}^1{ }_0 = -\frac{a'(t)}{a(t)},\ \underline{\theta}^1{ }_3 = - \frac{1}{a(t)r}\sqrt{1-k r^2} \\
      \underline{\theta}^2{ }_b \wedge \underline{\omega}^b \\
      \underline{\theta}^3{ }_b \wedge \underline{\omega}^b 
    \end{cases}
  \end{align*}


  Reversing the order of the wedge products used to read the connection one-form, we get
  \begin{align*}
    \begin{cases}
    \underline{\theta}^1{ }_0 = \frac{a'(t)}{a(t)}\underline{\omega}^1 + [\cdots] \underline{\omega}^0,\quad \underline{\theta}^1{ }_3 = \frac{1}{a(t)r}\sqrt{1-k r^2}\underline{\omega}^1 + [\cdots] \underline{\omega}^3\\
    \underline{\theta}^2{ }_0 = \frac{a'(t)}{a(t)}\underline{\omega}^2 + [\cdots]\underline{\omega}^0,\quad \underline{\theta}^2{ }_3 = \frac{1}{a(t)r}\sqrt{1-k r^2}\underline{\omega}^2 + [\cdots]\underline{\omega}^3,\quad \underline{\theta}^2{ }_1 = \frac{1}{a(t)r}\cot \theta \underline{\omega}^2 + [\cdots]\underline{\omega}^1\\
    \underline{\theta}^3{ }_0 = \frac{a'(t)}{a(t)}\underline{\omega}^0 + [\cdots]\underline{\omega}^3, \quad \underline{\theta}^3{ }_0 = [\cdots]\underline{\omega}^3, \quad  \underline{\theta}^3{ }_0 = [\cdots]\underline{\omega}^3
    \end{cases}
    
  \end{align*} 
    \begin{align*}
    \eta_{ca}\underline{\theta}^c{ }_b + \eta_{ca} \underline{\theta}_b{ }^c = 0 \iff  \eta_{ca}\eta^{ad}\underline{\theta}^c{ }_b + \eta_{ca}\eta^{ad} \underline{\theta}_b{ }^c = \underline{\theta}^d{ }_b + \underline{\theta}_b{ }^d = 0.
  \end{align*}


   then we can write $\omega = g(x)$ with $g(x) \in C^{\infty}(\mathbb{R})$. This time $\text{d}\omega = 0$ imposes the restriction $\partial_x g(x) \text{d}x = 0 \iff g(x) = c \in \mathbb{R} \ \forall x \in \mathbb{R}$.  


   As a first step towards this evaluation, $U\subset\mathbb{S}^2$ is mapped to $V\subset\mathbb{R}^2$ through the spherical chart map $\varphi : p\in U \mapsto (x, y)$ where integration can be performed under the rules of real analysis. To bring the integrand on $U$ to an integrand on $V$, we use the pullback of the diffeomorphism defined by $\varphi^{-1}: (x, y) \to (\sin(x) \cos(y), \sin(x)\sin(y), \cos(x))$. We have that the pullback of $F^{(2)}$ under this map is
  \begin{align*}
    \varphi^{-1}_*F^{(2)} = \sin(\theta) (\varphi^{-1}_*\text{d}\theta) \wedge (\varphi^{-1}_*\text{d}\phi) = 
  \end{align*} 

 $\phi = \arctan_2(v_\pm, u_\pm)$ and $\theta = \arctan_2(\pm \frac{1 - u^2_\pm + v_\pm^2}{1 + u^2_\pm + v_\pm^2}, \sqrt{u^2_\pm + v_\pm^2})$ since our sphere has radius 1 allowing us to get 
  \begin{align*}
    u^2_\pm + v_\pm^2 = (1-z^2)/(1\mp z)^2 = (1 \pm z)/(1 \mp z) \implies  z = \pm \frac{1 - u^2_\pm + v_\pm^2}{1 + u^2_\pm + v_\pm^2}.
  \end{align*}
  From these expressions, we have the following results 
  \begin{align*}
    \sin\theta = \sqrt{1-z^2} = \frac{((1 + u^2_\pm + v_\pm^2)^2 -(1 - u^2_\pm + v_\pm^2)^2)^{1/2}}{1 + u^2_\pm + v_\pm^2} = \frac{2\sqrt{u^2_\pm + v_\pm^2}}{1 + u^2_\pm + v_\pm^2},\quad \text{d}\phi =  \dfrac{-v_\pm \text{d}u_\pm + u_\pm\text{d}v_\pm}{u_\pm^2 + v_\pm^2}
  \end{align*}
  Wolfram alpha was used for the last one
  \item[(g)]

  &=
    \begin{cases}
      -Q (1 + z) \left(\frac{(1-u^2-v^2)((1+u^2+v^2) - (1-u^2-v^2))}{(1+u^2+v^2)^2}  +  4\frac{v^2 + u^2}{(1 + u^2 + v^2)^2}\right), \quad (+)\\
      -Q (1 - z) \left(\frac{(1-u^2-v^2)((1+u^2+v^2) + (1-u^2-v^2))}{(1+u^2+v^2)^2}  - 4\frac{v^2 + u^2}{(1 + u^2 + v^2)^2}\right), \quad (-)
    \end{cases}\\


To express $F^{(2)}$ in these new coordinates, we notice that $x = u_\pm (1\mp z)$ and $y = v_\pm (1\mp z)$. Since our sphere has radius 1, we have 
  \begin{align*}
    u^2_\pm + v_\pm^2 = (1-z^2)/(1\mp z)^2 = (1 \pm z)/(1 \mp z) \implies u^2_\pm + v_\pm^2 \mp z(u^2_\pm + v_\pm^2)  = 1 \pm z\implies  z = \pm \frac{1 - u^2_\pm - v_\pm^2}{1 + u^2_\pm + v_\pm^2}
  \end{align*}
  leading to $\text{d}x = (1\mp z) \text{d}u_\pm \mp u_\pm \text{d}z$, $\text{d}y = (1\mp z) \text{d}v_\pm \mp v_\pm \text{d}z$ and $\text{d}z = A \text{d}u_\pm + B \text{d}v_\pm$
  where 
  \begin{align*}
    A = -\pm\frac{2u_\pm (1 + u^2_\pm + v_\pm^2)}{(1 + u^2_\pm + v_\pm^2)^2}-\pm\frac{2u_\pm(1 - u^2_\pm - v_\pm^2)}{(1 + u^2_\pm + v_\pm^2)^2} = \mp\frac{4u_\pm}{(1 + u^2_\pm + v_\pm^2)^2}, \quad B = \mp\frac{4v_\pm}{(1 + u^2_\pm + v_\pm^2)^2}.
  \end{align*}
  We can also relate the two-form frame fields in cartesian coordinates to the $\text{d}u \wedge \text{d}v$ frame field as (omitting $\pm$ on $u,v$ symbols from now on)
  \begin{align*}
    \text{d}x \wedge \text{d}y &= ((1\mp z) \text{d}u \mp u\text{d}z)\wedge((1\mp z) \text{d}v \mp v \text{d}z)\\ &= (1\mp z)^2 \text{d}u \wedge \text{d}v \mp (1\mp z)(Bv + Au)\text{d}u \wedge \text{d}v\\ &= (1 \mp z)^2 \text{d}u \wedge \text{d}v + 4(1\mp z)\frac{v^2 + u^2}{(1 + u^2 + v^2)^2}\text{d}u \wedge \text{d}v\\
    \text{d}y \wedge \text{d}z &= (1\mp z) A\text{d}v \wedge \text{d}u = \pm(1\mp z)\frac{4u}{(1 + u^2 + v^2)^2}\text{d}u \wedge \text{d}v\\
    \text{d}z \wedge \text{d}x &= (1\mp z) B\text{d}v \wedge \text{d}u = \pm(1\mp z)\frac{4v}{(1 + u^2 + v^2)^2}\text{d}u \wedge \text{d}v
  \end{align*}
  With these expressions we are ready to express $F^{(3)}$ in the $\text{d}u$ and $\text{d}v$ frame field (we omit the $\pm$ on $u, v$ in what follows) as 
  \begin{align*}
    F^{(3)} &=  Q \frac{1}{r^3}\ \left(-y\text{d}z \wedge \text{d}x - x\text{d}y \wedge \text{d}z - z\text{d}x \wedge \text{d} y\right)\\
    &= -Q \left(z(1 \mp z)^2  + 4z(1\mp z)\frac{v^2 + u^2}{(1 + u^2 + v^2)^2} \pm(1\mp z)^2\frac{4u^2 + 4v^2}{(1 + u^2 + v^2)^2}\right)\\
    &=-Q \left(z(1 \mp z)^2  + 4(z \pm (1\mp z))(1\mp z)\frac{v^2 + u^2}{(1 + u^2 + v^2)^2}\right)\\
    &=-Q (1 \mp z) \left(\frac{(1-u^2-v^2)((1+u^2+v^2) \mp (1-u^2-v^2))}{(1+u^2+v^2)^2}  \pm  4\frac{v^2 + u^2}{(1 + u^2 + v^2)^2}\right)\\    
    &=
    \begin{cases}
      -Q \frac{2}{1+u^2 + v^2} \left(-\frac{2(u^2+v^2)^2}{(1+u^2+v^2)^2}  +  6\frac{v^2 + u^2}{(1 + u^2 + v^2)^2}\right), \quad (+)\\
      -Q \frac{2}{1+u^2 + v^2} \left(\frac{2}{(1+u^2+v^2)^2}  - 6\frac{v^2 + u^2}{(1 + u^2 + v^2)^2}\right), \quad (-)
    \end{cases}
  \end{align*}


  We now want to apply the chart transition map from the coordinates on $U_+$ to coordinates on $U_-$ to express $A_+$ in terms of $A_-$. To do so, we consider a fixed set of coordinates $x, y, z$ and pick a compare their $u, v$ coordinates. We have the relations
  \begin{align*}
    u_{+} = \frac{u_- (1+z)}{1-z} = \frac{u_- \left(1-\frac{1 - u_-^2 - v_-^2}{1 + u_-^2 + v_-^2}\right)}{1+\frac{1 - u_-^2 - v_-^2}{1 + u_-^2 + v_-^2}} = \frac{u_- \left(1-\frac{1 - u_-^2 - v_-^2}{1 + u_-^2 + v_-^2}\right)}{1+\frac{1 - u_-^2 - v_-^2}{1 + u_-^2 + v_-^2}} 
  \end{align*}


   $\star F_{\mu\nu} \text{d}x^{\mu} \wedge \text{d}x^{\nu}/2 = (\partial_\mu A_\nu - \partial_\nu A_\mu) \text{d}x^{\mu} \wedge \text{d}x^{\nu}/2 $


   &\star \text{d}\phi\wedge \text{d}\theta = \sqrt{(r^2 + n^2)^2 \sin^2 \theta}\varepsilon_{\mu \nu tr} g^{\mu\phi}g^{\nu\theta}\text{d}t \wedge \text{d}r = (r^2 + n^2)\sin \theta \varepsilon_{\phi \theta t r} g^{\phi \phi}g^{\theta \theta}\text{d}\phi \wedge \text{d}\theta = \frac{1}{(r^2 + n^2) \sin \theta} \text{d}t \wedge \text{d}r,\\
    %
    &\star \text{d}\phi \wedge \text{d}r = \sqrt{(r^2 + n^2)^2 \sin^2 \theta} \varepsilon_{\mu \nu \phi r} g^{\mu t}g^{\nu \theta} \text{d}t \wedge \text{d}\theta = (r^2 + n^2) \sin {\theta} \varepsilon_{t \theta \phi r} g^{t t}g^{\theta \theta} \text{d}t \wedge \text{d}\theta = -\left(\frac{4 A_\sigma^2}{(r^2 + n^2)\sin \theta} - \frac{\sin\theta}{f}\right)  \text{d}t \wedge \text{d}\theta,\\
    %
    &\star \text{d}t \wedge \text{d}\theta = \sqrt{(r^2 + n^2)^2 \sin^2 \theta} \varepsilon_{\mu \nu t\theta} g^{\mu \phi}g^{\nu r} \text{d}\phi \wedge \text{d}r = (r^2 + n^2) \sin \theta \varepsilon_{\phi r t\theta} g^{\phi \phi}g^{r r} \text{d}\phi \wedge \text{d}r = \frac{f}{\sin\theta} \text{d}\phi \wedge \text{d}r,\\
    %
    &

    We can now simplify further calculation by taking the limit $r \to \infty$ to get 
  \begin{align*}
    &\star \text{d}t\wedge \text{d}r = 0,\\
    &\star \text{d}\phi \wedge \text{d}r = - \sin\theta \text{d}t \wedge \text{d}\theta,\\
    &\star \text{d}t \wedge \text{d}\theta = -\frac{1}{\sin\theta} \text{d}\phi \wedge \text{d}r,\\
    &\star \text{d}\phi \wedge \text{d}\theta = \left(\frac{4 A_\sigma^2}{\sin \theta} + (r^2 + n^2) \sin \theta\right)\text{d}t \wedge \text{d}r
  \end{align*}