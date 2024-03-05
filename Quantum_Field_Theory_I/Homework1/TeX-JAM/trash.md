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