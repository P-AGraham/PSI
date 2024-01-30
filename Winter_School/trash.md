We are interested in compact trajectories which make both the real and imaginary parts of $I_-$ vanish at a specific frequency $\Omega$. A basic trajectory choice is given by a circle with varying angular velocity. Chosing to align the wave vector $k_\mu = (\omega_k, -k, 0, 0)$ with the the $x^1$ axis and expressing the circular trajectory as $x(\tau) = (aN(\tau), aM(\tau)\cos(\theta(\tau)), aM(\tau)\sin(\theta(\tau)))$ (normalizing factor $N(\tau)$)
the argument of the exponential becomes 
$$
 \Omega \tau + \omega_k x^0(\tau) - k x^1(\tau) = \Omega \tau + \omega_k N(\tau) - a k M(\tau) \cos(\theta(\tau)). 
$$
To fix $N(\tau)$, we require that
$$
\begin{align*}
1 = \dot{x}(\tau)^2 &= a^2\dot{N}(\tau)^2 - a^2(\dot{M}(\tau) \sin(\theta(\tau)) + M(\tau) \dot{\theta}(\tau )\cos(\theta(\tau)))^2 - a^2(\dot{M}(\tau) \cos(\theta(\tau)) - M(\tau) \dot{\theta}(\tau )\sin(\theta(\tau)))^2\\
&= a^2\dot{N}(\tau)^2 - a^2\dot{M}(\tau)^2 - (aM(\tau) \dot{\theta}(\tau))^2
\end{align*}
$$
Then take, $N= \sinh(f(\tau)), M = \cosh(f(\tau))$ to get
$$
\begin{align*}
1 
&= a^2\dot{f}(\tau)^2 - (aM(\tau) \dot{\theta}(\tau))^2
\end{align*}
$$
A natural choice of $\dot{\theta}(\tau)$ is then $\omega/\cosh(f(\tau))$ associated to $f(\tau) = \tau \sqrt{\omega+1/a} \equiv \bar{\omega} \tau$. Combined with $a = \sqrt{2}$, this choice implies that the argument of the exponential will involve a term 
$$
a k \cosh(\bar{\omega} \tau) \dfrac{1}{\sqrt{1-\bar{\omega}^2\tau^2}}
$$