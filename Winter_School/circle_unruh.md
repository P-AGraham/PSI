# Acceleration induced transparency for a circular trajectory

We consider an UDW detector moving through as scalar field $\phi$ representing photons with with wavevector $k_{\mu}$. The peicewise uniform circular trajectory $x(\tau)$ of the detctor is paramatrized by proper time $\tau$ as 
$$
\begin{align*}
x(\tau) = 
\begin{cases}
\gamma_{\omega_1 a} (\tau,\quad a\cos(\omega_1\tau),\quad a\sin(\omega_1 \tau),\quad 0)\quad\tau \in (-T, -T/2)\\
\gamma_{\omega_2 a} (\tau,\quad a\cos(\omega_2\tau - \omega_1 T/2),\quad a\sin(\omega_2 \tau - \omega_1 T/2),\quad 0)\quad\tau \in (-T/2, +T/2)\\
\gamma_{\omega_3 a} (\tau,\quad a\cos(\omega_3\tau-\omega_1 T/2+\omega_2 T/2),\quad a\sin(\omega_3\tau-\omega_1 T/2+\omega_2 T/2),\quad 0)\quad\tau \in (+T/2, T)
\end{cases}
\end{align*}
$$
where $a$ is the radius and $\omega_1, \omega_2, \omega_3$ are angular frequencies in different pieces of the trajectory seperated by two infinite acceleration shifts at $\tau = \pm T/2$. The $\gamma = (1-\omega_{i}^2a^2)^{-1/2}$ factors ensure $\dot {x}_\mu \dot {x}^\mu = 1$. The resonant excitation amplitude of the detector movin on this trajectory is given by the integral
$$
I_{-}(\Omega, \mathbf{k})=\int \mathrm{d} \tau \ e^{i \Omega \tau - i k_\mu x^\mu(\tau)}.
$$
Chosing to align the wave vector $k_\mu = (k, k, 0, 0)$ with the the $x^1$ axis, the argument of the exponential becomes 
$$
 \Omega \tau + \omega_k x^0(\tau) - k x^1(\tau) = (\Omega + k \gamma_{\omega_i a}) \tau - a k \cos(\omega_i \tau + \delta_i). 
$$
where $\delta_i = \sum_j^i \tau_j \omega_j$ represents the phase shifts making the trajectory contiuous. 
To simplify, we denote $\Omega_i = \Omega + k \gamma_{\omega_i a}$ and use the Jaccobi-Anger expansion to get 
$$
\begin{align*}
I_{-}(\Omega, \mathbf{k}) &= \sum_i \int_i \mathrm{d} \tau \ e^{i \Omega_i \tau} e^{-i a k \cos(\omega_i \tau+\delta_i)}\\
&=\sum_i \int_i \mathrm{d} \tau \ e^{i \Omega_i \tau} \sum_{n=-\infty}^{\infty} i^n J_n(-ka) e^{in \omega_i \tau + i n\delta_i}\\
&=  \sum_{n=-\infty}^{\infty} i^n  J_n(-ka) \sum_i  \int_{\tau_i}^{\tau_{i+1}} \mathrm{d} \tau \  e^{i \Omega_i \tau + i n \omega_i \tau + i n\delta_i}\\
&= \sum_{n=-\infty}^{\infty} i^n  J_n(-ka) \sum_i e^{i n\delta_i}\frac{1}{i (\Omega_i +  n \omega_i)}\left( e^{i (\Omega_i +  n \omega_i) \tau_{i+1}} - e^{i (\Omega_i +  n \omega_i) \tau_{i}}\right) 
\end{align*}
$$
For our specific peicewise inertial trajectory, we have 
$$
\begin{align*}
\sum_i e^{i n\delta_i}\frac{1}{i (\Omega_i +  n \omega_i)}\left( e^{i (\Omega_i +  n \omega_i) \tau_{i+1}} - e^{i (\Omega_i +  n \omega_i) \tau_{i}}\right)  &= e^{i n\delta_1}\frac{1}{i (\Omega_1 +  n \omega_1)}\left(e^{-i (\Omega_1 +  n \omega_1)T/2} - e^{-i (\Omega_1 +  n \omega_1) T}\right) \\
&+ e^{i n\delta_2}\frac{1}{i (\Omega_2 +  n \omega_2)}\left( e^{i (\Omega_2 +  n \omega_2) T/2} - e^{-i (\Omega_2 +  n \omega_2) T/2}\right)\\
&+ e^{i n\delta_3}\frac{1}{i (\Omega_3 +  n \omega_3)}\left( e^{i (\Omega_3 +  n \omega_3) T} - e^{i (\Omega_3 +  n \omega_3) T/2}\right)\\
\end{align*}
$$
Lets now treat the case where all frequencies are equal. We have the simplified result 
$$
\begin{align*}
I_{-}(\Omega, \mathbf{k}) &= \sum_{n=-\infty}^{\infty} i^{n}  J_n(-ka) \frac{1}{i(\Omega_3 +  n \omega_3)}\left(e^{i (\Omega_1 +  n \omega_1) T} - e^{-i (\Omega_1 +  n \omega_1) T}\right)
\\
&= \sum_{n=-\infty}^{\infty} i^{n}  J_n(-ka) \frac{\sin((\Omega + k \gamma_{\omega_i a} + n \omega_1)T)}{(\Omega + k \gamma_{\omega_i a} +  n \omega_1)}
\end{align*}
$$
