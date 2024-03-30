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


  Suppose $I_-(\Omega_0) = 0$, then we can write $I_-(\Omega)  = f(\Omega) (\Omega-\Omega_0)$ which is equivalent


  \frac{\partial \hat{L}}{\partial v^j \partial v^i} (\hat{q}, \hat{v})  \text{d}v^{j}_{\mathbf{q}, \mathbf{v}} \wedge \text{d}q^{i}_{\mathbf{q}, \mathbf{v}}


  To simplify further calculations we make the product states more explicit with the binomial expansion 
  \begin{align*}
    \ket{\psi(t)}= \sum_{n=0}^{2^L - 1} (\cos(t/2))^{L-\alpha_n} (\sin(t/2))^{\alpha_n}\ket{n}
  \end{align*}
  where $n$ is the natural number having binary expansion equal to the sequence of $1$ and $0$ (each term in the expansion of $(\cos(t/2)\ket{0} + \sin(t/2)\ket{1})^{\otimes L}$ is obtained by making $L$ binary choices and every number from $0$ to $2^{L}-1$ is reached for one and only one choice). The quantity $\alpha_n$ represents the number $1$ in the binary expansion of $n$. 


  \begin{align*}
    \left(- J \sum_{i=0}^{L-1} Z_i Z_{i+1}\right)\ket{n} = -J(L - 2JD_n) \ket{n}
  \end{align*}
  where $D_n$ is the number of domain walls in $\ket{n}$ (number of sites with different neighbours). The action of each $Z_i Z_{i+1}$ for a given $n$ will contribute $-J$ if the neighbours $i$ and $i+1$ are the same and $+J$ if they are different. If they were all the same we would have a total of $-JL$. From this total we can obtain the general energy by adding by adding $2JD_n$ (the factor of two splits in a contribution canceling the overshoot of alignement energy cost and a contribution for the actual anialigned energy cost). From this result, we get 
  \begin{align*}
    \bra{\psi(t)} \left(- J \sum_{i=0}^{L-1} Z_i Z_{i+1}\right) \ket{\psi(t)} = \sum_{n=0}^{2^L - 1} -J(L - 2JD_n) (\cos(t/2))^{2(L-\alpha_n)} (\sin(t/2))^{2\alpha_n}
  \end{align*}


    We check that for $t=0$, the only non zero term in the previous result is realized when the power of $\sin$ vanishes (preventing $\sin$ to cancel the entire expression) which implies that $\alpha_n = 0 \implies n = 0$ (no domain walls) leading to an average energy $-JL$ as expected. 

    This property is acheived by taking $M(0)$ and $M(1)$ to anticommute and square to identity. A valid choice is $M(0) = Z$ and $M(1)= X$ (pauli matrices for $D=2$). Using anticommutation we can regroup all the $X$ and $Z$ of a product at the cost of a global $\pm 1$ factor ($-1$ for an odd number of transpositions and $+1$ for an even number of transpositions). Then, we can use the fact $X, Z$ square to one to reduce the $X, Z$ blocks of the chain a $X, Z, XZ, 1$   

    M(0) = 
    \begin{pmatrix}
      1 & 1 & 0 \\
      0 & 0 & 1 \\
      0 & 0 & 0
    \end{pmatrix}, 
    M(0)^2 = 
    \begin{pmatrix}
      1 & 1 & 0 \\
      0 & 0 & 1 \\
      0 & 0 & 0
    \end{pmatrix}
    \begin{pmatrix}
      1 & 1 & 0 \\
      0 & 0 & 1 \\
      0 & 0 & 0
    \end{pmatrix}
    = 
    \begin{pmatrix}
      1 & 1 & 0 \\
      0 & 0 & 1 \\
      0 & 0 & 0
    \end{pmatrix}


    $[M(0), M(1)] = 0$ $M(0)^3 = 0$,  $M(0)^2 \neq 0$. Following these constraints, matrix $M(1)$ is taken to be the identiy matrix up to some normalisation factor $a \in \mathbb{C}$. Since identity commutes with $M(0)$ we can regroup all $M(0)$ matrices. For powers higher than two, the entire chain of matrix multiplication will vanish. We also want single powers of $M(0)$ to lead to vanishing contribution so we have $Tr(M(0) 1^{L-1}) = Tr(M(0)) = 0$. However, $Tr(M(0)^2)$ must be non-zero to produce the $\frac{J}{4g}$ contribution. A matrix satisfying all these constraints exists for $D=3$ and reads 


    Through a local unitary mapping (this unitary might affect entanglement within each component of the bipartition, but does not entangle/disentangle the components), we can send $\ket{u_i}$, $\ket{v_i}$ to states in our original basis $\ket{S_i}$, $\ket{S_i'}$. The state we obtain is $\ket{\psi} = \sum_{i} \alpha_i \ket{S_i}\ket{S_i'}$.



For a spacial conformal transformation, we have $\xi^\mu = x^2 b^\mu -2 x^\mu x_\lambda b^\lambda$ paramatrized by the transaltion vector $b^\lambda$ around $\infty$. For this vector, we have 
  \begin{align*}
    R_{\mu}^{\nu} (x) &= \delta_\nu^\mu + \delta_\nu^\mu \partial_\sigma (x_\lambda x^\lambda b^\sigma)  -2 \delta_\nu^\mu \partial_\sigma  (x^\sigma x_\lambda b^\lambda) + 2\partial_\mu  x^\nu x_\lambda b^\lambda - \partial_\mu  (x^\lambda x_\lambda b^\nu) +  O(\xi^2)\\
    &= \delta_\nu^\mu + 2\delta_\nu^\mu x_\sigma b^\sigma -2 (3) \delta_\nu^\mu x_\lambda b^\lambda -2 \delta_\nu^\mu x^\sigma  \eta_{\rho \lambda}  \delta^\rho_\sigma b^\lambda + 2  \delta^\nu_\mu x_\lambda b^\lambda + 2 x^\nu \eta_{\rho\lambda} \delta^{\rho}_{\mu} b^{\lambda} - \partial_\mu  (x^\lambda x_\lambda b^\nu) +  O(\xi^2) \\
    &= \delta_\nu^\mu + 2\delta_\nu^\mu x_\sigma b^\sigma -2 (3) \delta_\nu^\mu x_\lambda b^\lambda -2 \delta_\nu^\mu x_\sigma b^\sigma + 2 \delta^\nu_\mu x_\lambda b^\lambda + 2 x^\nu b_{\mu} - \partial_\mu  (x^\lambda x_\lambda b^\nu) +  O(\xi^2)  \\
    &= \delta_\nu^\mu +  O(\xi^2)
  \end{align*}
  As expected for a special conformal transformation.


  The inverse transformation (given at first order in $\xi$ by ${x}^{\mu} = f^{-1}(\tilde{x}) = \tilde{x}^{\mu} - \xi^{\mu}(\tilde{x})$) reads 
  \begin{align*}
    F_{\mu\nu}(x)
    &= \tilde{F}_{\mu \nu}(x) - \tilde{F}_{\mu \nu}(x) \frac{\Delta}{D}\partial_\lambda \xi^\lambda(x) - \tilde{A}_{(\nu}(x) \frac{\Delta}{D}\partial_{\mu)} \partial_\lambda \xi^\lambda(x)  - (\partial_{(\mu} \tilde{A}_\lambda(x))M^{\lambda}_{\nu)} - \xi^{\lambda}(\tilde{x}) \partial_\lambda F_{\mu\nu}(\tilde{x})
  \end{align*}


  which reads 
  \begin{align*}
   \begin{pmatrix}
      1 & 1\\
      -1 & 1\\
   \end{pmatrix}^{1/2}
    = \left(\begin{pmatrix}
      2 & 0\\
      0 & 0\\
   \end{pmatrix}
   \begin{pmatrix}
    1 & 1\\
    1 & -1\\
 \end{pmatrix}
   \right)

  \end{align*}
  where we used the symmetry of the block to deduce its eigenvactor have components $(1, -1)$ and $(1, 1)$ with respective eigenvalues $0$ and $2$. 



  As before, we simplify the calculation by applying the $\sigma_x$ local operation to subsystem $A$ to work with $\rho = \ket{\Psi^{+}_{AB}}\bra{\Psi^{+}_{AB}}$.


  (1_1 \otimes \ket{0}\bra{0}_2 + X_1 \otimes \ket{1}\bra{1}_2)(\ket{0}\bra{0}_1 \otimes 1_2  + \ket{1}_1\bra{1}_1 \otimes X_2)(1_1 \otimes \ket{0}\bra{0}_2 + X_1 \otimes \ket{1}\bra{1}_2)\\
  &= (1_1 \otimes \ket{0}\bra{0}_2 + X_1 \otimes \ket{1}\bra{1}_2)(\ket{0}\bra{0}_1 \otimes 1_2  + \ket{1}_1\bra{1}_1 \otimes X_2)
  
  \ket{0}\bra{0}_1 \otimes \ket{0}\bra{0}_2  + \ket{1}_1\bra{1}_1 \otimes \ket{0}\bra{1}_2 + \ket{0}\bra{0}_1 \otimes 1_2  + \ket{1}_1\bra{1}_1 \otimes X_2

  
  
  
  (1_1 \otimes \ket{0}\bra{0}_2 + X_1 \otimes \ket{1}\bra{1}_2)




  \begin{pmatrix}
      1 & 0 & 0 & 0\\
      0 &  e^{i \pi/4} & 0 & 0\\
      0 & 0 & 1 & 0\\
      0 & 0 & 0 & e^{i \pi/4}
    \end{pmatrix}
    \begin{pmatrix}
      1 & 0 & 0 & 0\\
      0 & 0 & 0 & 1\\
      0 & 0 & 1 & 0\\
      0 & 1 & 0 & 0
    \end{pmatrix}


    In what follows, we evaluate the matrix expressions representing a quantum circuit unitary acting on a sequence of qubit input. We work in the computationnal basis $\{\ket{0}, \ket{1}\}$ and use the notation $X, Y, Z$ for the Pauli gates in this basis. 
\begin{enumerate}
  \item[(a)] First we consider the conjugation of a \verb|CNOT| by two \verb|CNOT| with control and target qubit reversed:
  \begin{align*}
  \begin{quantikz}[baseline={([yshift=-.5ex]current bounding box.center)}]
    \lstick{$1$} & \targ{}  & \ctrl{1}  & \targ{}  & \qw\\
    \lstick{$2$} & \ctrl{-1}  & \targ{}  & \ctrl{-1} & \qw
  \end{quantikz}
  &= 
  \begin{pmatrix}
    1 & 0 & 0 & 0\\
    0 & 0 & 0 & 1\\
    0 & 0 & 1 & 0\\
    0 & 1 & 0 & 0
  \end{pmatrix}
  \begin{pmatrix}
    1 & 0 & 0 & 0\\
    0 & 1 & 0 & 0\\
    0 & 0 & 0 & 1\\
    0 & 0 & 1 & 0
  \end{pmatrix}
  \begin{pmatrix}
    1 & 0 & 0 & 0\\
    0 & 0 & 0 & 1\\
    0 & 0 & 1 & 0\\
    0 & 1 & 0 & 0
  \end{pmatrix}
  =
  \begin{pmatrix}
    1 & 0 & 0 & 0\\
    0 & 0 & 1 & 0\\
    0 & 1 & 0 & 0\\
    0 & 0 & 0 & 1
  \end{pmatrix}
  \end{align*}
  which exchanges the qubits ($\ket{00} \to \ket{00}, \ket{01} \to \ket{10},  \ket{10} \to \ket{01}, \ket{11} \to \ket{11}$) and constitutes a \verb|SWAP| gate. The matrix expression for the reversed \verb|CNOT| was obtained by writting its action on the computationnal basis wich reads $\ket{00} \to \ket{00}, \ket{01} \to \ket{11},  \ket{10} \to \ket{10}, \ket{11} \to \ket{01}$. 
  \item[(b)] Then we calculate the matrix expression of the entanglement generating circuit 
  \begin{align*}
    \begin{quantikz}[baseline={([yshift=-.5ex]current bounding box.center)}]
      \lstick{$1$} &      &   & \targ{}  & \qw\\
      \lstick{$2$} & \gate{H}  & \gate{R_{\pi/4}}  & \ctrl{-1} & \qw
    \end{quantikz} &= \left( 1_1 \otimes \frac{1}{\sqrt{2}} (X_2 + Z_2) \right) R_{\pi/4, 2} (1_1 \otimes \ket{0}\bra{0}_2 + X_1 \otimes \ket{1}\bra{1}_2)\\
    &= \left( 1_1 \otimes \frac{1}{\sqrt{2}} (X_2 + Z_2) \right) (1_1 \otimes \ket{0}\bra{0}_2 + X_1 \otimes e^{i\pi/4} \ket{1}\bra{1}_2)\\
    &= \frac{1}{\sqrt{2}} \left(1_1 \otimes (\ket{0}_2 + \ket{1}_2)\bra{0}_2 + X_1 \otimes e^{i\pi/4} (\ket{0}_2 - \ket{1}_2)\bra{1}_2\right)\\
    &= \frac{1}{\sqrt{2}} \left(\begin{pmatrix}
      1 & 0\\ 
      0 & 1 
    \end{pmatrix} \otimes 
    \begin{pmatrix}
      1 & 0 \\
      1 & 0 
    \end{pmatrix}
       + \begin{pmatrix}
        0 & 1\\
        1 & 0
      \end{pmatrix} \otimes  \begin{pmatrix}
        0 & e^{i\pi/4} \\
        0 & -e^{i\pi/4}
      \end{pmatrix}\right)\\
      &= \frac{1}{\sqrt{2}} \begin{pmatrix}
        1 & 0 & 0 & 0 \\ 
        1 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0 \\
        0 & 0 & 1 & 0 
      \end{pmatrix} +
      \frac{1}{\sqrt{2}}
      \begin{pmatrix}
        0 & 0 & 0 & e^{i\pi/4} \\ 
        0 & 0 & 0 & -e^{i\pi/4} \\
        0 & e^{i\pi/4} & 0 & 0 \\
        0 & -e^{i\pi/4} & 0 & 0 
      \end{pmatrix} = 
      \frac{1}{\sqrt{2}}
      \begin{pmatrix}
        1 & 0 & 0 & e^{i\pi/4} \\ 
        1 & 0 & 0 & -e^{i\pi/4} \\
        0 & e^{i\pi/4} &  1 & 0 \\
        0 & -e^{i\pi/4} & 1 & 0 
      \end{pmatrix}.
  \end{align*}
  If we set the phases to $1$, we recover the Bell state mapping $\ket{00} \to (\ket{00} + \ket{11})/\sqrt{2}$, $\ket{01} \to (\ket{00} - \ket{11})/\sqrt{2}$, $\ket{10} \to (\ket{01} + \ket{10})/\sqrt{2}$ and $\ket{11} \to (-\ket{01} + \ket{10})/\sqrt{2}$. 
  \item[(c)]  Finally, we calculate the matrix expression associated with a three qubit circuit as follows:
  \begin{align*}
    \begin{quantikz}[baseline={([yshift=-.5ex]current bounding box.center)}]
      \lstick{$1$} &  \gate{H}         & \ctrl{1} &          & \ctrl{1} & \qw\\
      \lstick{$2$} &                   & \targ{1} & \gate{H} & \ctrl{1} & \qw\\
      \lstick{$3$} &  \gate{R_{\pi/4}} &          &          & \targ{1} & \qw
    \end{quantikz}
    &= H_1 (\ket{0}\bra{0}_1 \otimes 1_2 + \ket{1}\bra{1}_1 \otimes X_2) H_2 R_{\pi/4, 3}\verb!TOFFOLI!\\
    &= \frac{1}{\sqrt{2}}((\ket{0} + \ket{1})\bra{0}_1 \otimes H_2 +  (\ket{0} - \ket{1})\bra{1}_1 \otimes X_2 H_2)R_{\pi/4, 3}\verb!TOFFOLI!\\
    &= \frac{1}{2}\left(
      \begin{pmatrix}
        1 & 0 \\
        1 & 0
      \end{pmatrix}
       \otimes \begin{pmatrix}
        1 & 1 \\
        1 & -1
      \end{pmatrix} + \begin{pmatrix}
        0 & 1 \\
        0 & -1
      \end{pmatrix} \otimes \begin{pmatrix}
        0 & 1 \\
        1 & 0
      \end{pmatrix}\begin{pmatrix}
        1 & 1 \\
        1 & -1
      \end{pmatrix}\right)R_{\pi/4, 3}\verb!TOFFOLI!\\
      & = \frac{1}{2}\left(
      \begin{pmatrix}
        1 & 1 & 0 & 0\\
        1 & - 1 & 0 & 0\\
        1 & 1 & 0 & 0\\
        1 & - 1 & 0 & 0
      \end{pmatrix}
      +
      \begin{pmatrix}
        0 & 0 & 1 & -1\\
        0 & 0 & 1 & 1\\
        0 & 0 & -1 & 1\\
        0 & 0 & -1 & -1\\
      \end{pmatrix}\right)R_{\pi/4, 3}\verb!TOFFOLI!\\ &= \frac{1}{2}
        \begin{pmatrix}
          1 & 1 & 1 & -1\\
          1 & - 1 & 1 & 1\\
          1 & 1 & -1 & 1\\
          \textcolor{blue}{1} & \textcolor{blue}{-1} & \textcolor{blue}{-1} & \textcolor{blue}{-1}
        \end{pmatrix}
       R_{\pi/4, 3}\verb!TOFFOLI!
  \end{align*}
Where the blue elements contribute to the $\ket{11}$ component of the output of the matrix. This component is assocated with the only state on qubits $1, 2$ that activates the \verb|TOFFOLI| gate on qubit $3$. We can explicit the effect of the \verb|TOFFOLI| gate as follows 
\begin{align*}
  \frac{1}{2}
        \begin{pmatrix}
          1 & 1 & 1 & -1\\
          1 & - 1 & 1 & 1\\
          1 & 1 & -1 & 1\\
          1 & - 1 & -1 & \textcolor{blue}{-1}
        \end{pmatrix}
       R_{\pi/4, 3}\verb!TOFFOLI!
\end{align*}
\end{enumerate}


We can explicit the effect of the \verb|TOFFOLI| gate as follows 
\begin{align*}
  \verb!TOFFOLI! \ R_{\pi/4, 3}\frac{1}{2}
  \begin{pmatrix}
    1 & 1 & 1 & 1\\
    1 & - 1 & 1 & -1\\
    1 & 1 & -1 & -1\\
    -1 & 1 & 1 & -1
  \end{pmatrix}
  =
  .
\end{align*}


 \begin{quantikz}[baseline={([yshift=-.5ex]current bounding box.center)}]
    \lstick{$a$} & \ctrl{1}  & \ctrl{1} & \\
    \lstick{$b$} & \ctrl{1} & \targ{}  & \\
    \lstick{$c$} & \targ{}   &          &
  \end{quantikz}