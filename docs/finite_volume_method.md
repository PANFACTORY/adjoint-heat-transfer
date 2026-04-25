# Finite Volume Method

The Finite Volume Method (FVM) is one of the partial differential equation solvers, and in this manuscript, we explain its application to a thermal diffusion problem as an example.

The governing equations are as follows:
$$
\begin{align}
    &\frac{\partial u}{\partial t}=-\frac{\partial q_\alpha}{\partial x_\alpha}+q_v && \text{in }D \\
    &q_\alpha=-k\frac{\partial u}{\partial x_\alpha} && \text{in }D \\
    &u=\bar{u} && \text{on }\Gamma_\text{D} \\
    &n_\alpha q_\alpha=-\bar{q} && \text{on }\Gamma_\text{N} \\
    &u\left(t_0\right)=u_0 && \text{in }D
\end{align}
$$ where $u$ and $q_\alpha$ are a temperature and heat flux, respectively. $q_v$, $k$, $n_\alpha$, $\bar{u}$, $\bar{q}$ and $u_0$ are a volumetric heat source, thermal diffusivity, an outward normal vactor on a boundary, temperature on the Dirichlet boundary $\Gamma_\text{D}$, heat flux on the Neumann boundary $\Gamma_\text{N}$, and the initial temperature, respectively. It should be noted that $\bar{q}$ is defined as an inward heat flux.

Eq. (1) is integrated over a control volume $V$, which is a subdomain of the analysis domain $D$, and Gauss's theorem is applied as follows:
$$
\begin{align}
    \int_V\frac{\partial u}{\partial t}dV=&-\int_V\frac{\partial q_\alpha}{\partial x_\alpha}dV+\int_Vq_vdV \notag\\
    =&-\int_{\partial V}n_\alpha q_\alpha dS+\int_Vq_vdV
\end{align}
$$ Based on the FVM, a charactaric temperature is defined at the center of $V$; therefore, Eq. (6) is descretized as follows:
$$
\begin{align}
    \therefore\:&\frac{u^{t+1}-u^t}{\Delta t}\Delta x\Delta y=\left(q_x\right)_\text{W}\Delta y-\left(q_x\right)_\text{E}\Delta y+\left(q_y\right)_\text{S}\Delta x-\left(q_y\right)_\text{N}\Delta x+q_v\Delta x\Delta y \notag\\
    \therefore\:&u^{t+1}=u^t-\Delta t\left\{\frac{\left(q_x\right)_\text{E}-\left(q_x\right)_\text{W}}{\Delta x}+\frac{\left(q_y\right)_\text{N}-\left(q_y\right)_\text{S}}{\Delta y}-q_v\right\}
\end{align}
$$

The heat flux is calculated as follows:
$$
\begin{align}
    &q_x^L=-k^L\frac{u^f-u^L}{\Delta x/2} \notag\\
    &q_x^R=-k^R\frac{u^R-u^f}{\Delta x/2} \notag\\
    &q_x^L=q_x^R \notag\\
    \therefore\:&-k^L\frac{u^f-u^L}{\Delta x/2}=-k^R\frac{u^R-u^f}{\Delta x/2} \notag\\
    \therefore\:&u^f=\frac{k^Lu^L+k^Ru^R}{k^L+k^R} \notag\\
    \therefore\:&q_x=q_x^L=q_x^R=-\frac{2k^Lk^R}{k^L+k^R}\frac{u^R-u^L}{\Delta x}
\end{align}
$$ $q_y$ is caluculated in the same manner as $q_x$.

For the Dirichlet boundary condition at $x_\text{min}$, $q_x$ is given as follows:
$$
\begin{equation}
    q_x=-k^c\frac{u^c-\bar{u}}{\Delta x/2}
\end{equation}
$$ Here, $k^c$ and $u^c$ are the thermal diffusivity and the cell-centered temperature for a control volume at $x_\text{min}$. On the other boundaries ($x_\text{max}$, $y_\text{min}$, and $y_\text{max}$), the same approach is applied. 

For Neumann boundary condition at $x_\text{min}$, $q_x$ is given as follows:
$$
\begin{align}
    -q_x&=-\bar{q} \notag\\
    \therefore\;q_x&=\bar{q}
\end{align}
$$ On the other boundaries ($x_\text{max}$, $y_\text{min}$ and $y_\text{max}$), the same approach is applied. 

By solving Eqs. (7)-(9) as a time marching scheme with the initial condition Eq. (5), we can obtain the temperature at time step $t$.