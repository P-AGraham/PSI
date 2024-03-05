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

We express the state of the system as $\ket{\psi} = \sum_{s'} c_{s;n} e^{i\Omega [s] t} \ket{s;n}$ ($[s]$ is the number of $e$ minus the number of $g$ over $2$ in $s$). For simplicity we specialize to a single detector in state $\ket{g; 1}$ to get 
\begin{align*}
    \frac{\text{d}}{\text{d}\tau} e^{i\Omega t/2} c_{g;1} &= \bra{g;1}\sum_{j = 0}^{N-1} (e^{i\Omega t}\sigma^+_j + e^{-i\Omega t} \sigma^-_j)(e^{-i x^\mu_j k_\mu} a_k + e^{i x_j^\mu k_\mu} a_k^{\dagger})\ket{g;1} = 0\\
    \frac{\text{d}}{\text{d}\tau} e^{-i\Omega t/2} c_{e;0} &= \bra{e;0}\sum_{j = 0}^{N-1} (e^{i\Omega t}\sigma^+_j + e^{-i\Omega t} \sigma^-_j)(e^{-i x^\mu_j k_\mu} a_k + e^{i x_j^\mu k_\mu} a_k^{\dagger})\ket{g;1} = 
\end{align*}

\begin{align*}
    \frac{\text{d}}{\text{d}\tau} e^{i\Omega[s] t/2} c_{s;n} =  
\end{align*}


We express the state of the system as $\ket{\psi} = \sum_{s'} c_{s;n} e^{i\Omega [s] t} \ket{s;n}$ ($[s]$ is the number of $e$ minus the number of $g$ over $2$ in $s$). For simplicity we specialize to a single detector in state $\ket{g; 1}$ to get 
\begin{align*}
    \frac{\text{d}}{\text{d}\tau} e^{i\Omega t/2} c_{g;n} &= \bra{g;n}\sum_{m}\sum_{j = 0}^{N-1} e^{-i\Omega t}\sigma^-_j(e^{-i x^\mu_j k_\mu} a_k + e^{i x_j^\mu k_\mu} a_k^{\dagger})c_{e;m}\ket{e;m}\\
    &= \bra{g;n}\sum_{m}\sum_{j = 0}^{N-1} e^{-i\Omega t}\sigma^-_j(e^{-i x^\mu_j k_\mu} \sqrt{m}c_{e;m}\ket{g;m-1} + e^{i x_j^\mu k_\mu}   \sqrt{m + 1}c_{e;m}\ket{e;m+1})\\
 &= e^{-i\Omega t-i x^\mu_j k_\mu} \sqrt{n+1}c_{e;n+1} + e^{-i\Omega t + i x_j^\mu k_\mu} \sqrt{n-1}c_{e;n-1}\\
    \frac{\text{d}}{\text{d}\tau} e^{-i\Omega t/2} c_{e;n} &= \bra{e;n}\sum_{j = 0}^{N-1} e^{i\Omega t}\sigma^+_j (e^{-i x^\mu_j k_\mu} a_k + e^{i x_j^\mu k_\mu} a_k^{\dagger})\ket{g;1} = 
\end{align*}


def fourier_wrap(f, start, end, delta_tau):

    # see https://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy

    tau_grid=np.arange(start, end, delta_tau)

    f_eval = np.array([f(tau) for tau in tau_grid])

    #f_inter = interp1d(tau_grid, f_eval)

    
    g_eval=np.fft.fft(f_eval)
    Omega_grid = np.fft.fftfreq(f_eval.size) * 2*np.pi/delta_tau

    g_eval *= delta_tau*np.exp(-1j * Omega_grid * delta_tau)

    return g_eval[:-100], Omega_grid[:-100]