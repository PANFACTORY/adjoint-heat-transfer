# Formulation

## Governing equation

The governing equations are as follows:
$$
\begin{align}
    &\frac{\partial u}{\partial t}=-\frac{\partial q_\alpha}{\partial x_\alpha}+q_v && \text{in }D \\
    &q_\alpha=-k\frac{\partial u}{\partial x_\alpha} && \text{in }D \\
    &u=\bar{u} && \text{on }\Gamma_\text{D} \\
    &n_\alpha q_\alpha=\bar{q} && \text{on }\Gamma_\text{N} \\
    &u\left(t_0\right)=u_0 && \text{in }D
\end{align}
$$
Based on the Finite Volume Method (FVM), they are discretized as follows:
$$
\begin{align}
    &\begin{aligned}
        \int_V\frac{\partial u}{\partial t}dV=&-\int_V\frac{\partial q_\alpha}{\partial x_\alpha}dV+\int_Vq_vdV \\
        =&-\int_Sn_\alpha q_\alpha dS+\int_Vq_vdV
    \end{aligned} \\
    \therefore\:&\frac{u^{t+1}-u^t}{\Delta t}\Delta x\Delta y=\left(q_x\right)_\text{W}\Delta y-\left(q_x\right)_\text{E}\Delta y+\left(q_y\right)_\text{S}\Delta x-\left(q_y\right)_\text{N}\Delta x+q_v\Delta x\Delta y \\
    \therefore\:&u^{t+1}=u^t-\Delta t\left\{\frac{\left(q_x\right)_\text{E}-\left(q_x\right)_\text{W}}{\Delta x}+\frac{\left(q_y\right)_\text{N}-\left(q_y\right)_\text{S}}{\Delta y}-q_v\right\}
\end{align}
$$

For flux:
$$
\begin{align}
    &q_x^L=-k^L\frac{u^f-u^L}{\Delta x/2} \\
    &q_x^R=-k^R\frac{u^R-u^f}{\Delta x/2} \\
    &q_x^L=q_x^R \\
    \therefore\:&-k^L\frac{u^f-u^L}{\Delta x/2}=-k^R\frac{u^R-u^f}{\Delta x/2} \\
    \therefore\:&u^f=\frac{k^Lu^L+k^Ru^R}{k^L+k^R} \\
    \therefore\:&q_x=q_x^L=q_x^R=-\frac{2k^Lk^R}{k^L+k^R}\frac{u^R-u^L}{\Delta x}
\end{align}
$$

For Dirichlet boundary condition on x min:
$$
\begin{equation}
    qx=-k^c\frac{u^c-u^b}{\Delta x/2}
\end{equation}
$$

## Adjoint equation

The lagrangian $L$ is defined as follows:
$$
\begin{align}
    &L=J+R_u+R_q+R_{db}+R_{nb} \\
    &R_u=\int_\mathcal{I}\int_\mathcal{O}\tilde{u}\left(\frac{\partial u}{\partial t}+\frac{\partial q_\alpha}{\partial x_\alpha}-q_v\right)d\Omega dt \\
    &R_q=\int_\mathcal{I}\int_\mathcal{O}\tilde{q}_\alpha\left(q_\alpha+k\frac{\partial u}{\partial x_\alpha}\right)d\Omega dt \\
    &R_{db}=\int_\mathcal{I}\int_{\Gamma_\text{D}}\tilde{\lambda}\left(u-\bar{u}\right)d\Gamma dt \\
    &R_{nb}=\int_\mathcal{I}\int_{\Gamma_\text{N}}\tilde{\lambda}\left(n_\alpha q_\alpha-\bar{q}\right)d\Gamma dt
\end{align}
$$
Here,
$$
\begin{align}
    \delta R_u=&\int_\mathcal{I}\int_\mathcal{O}\tilde{u}\left(\frac{\partial\delta u}{\partial t}+\frac{\partial\delta q_\alpha}{\partial x_\alpha}\right)d\Omega dt \notag\\
    =&\int_\mathcal{O}\left[\tilde{u}\delta u\right]_{t_0}^{t_1}d\Omega+\int_\mathcal{I}\int_{\partial\mathcal{O}}n_\alpha\tilde{u}\delta q_\alpha d\Gamma dt+\int_\mathcal{I}\int_\mathcal{O}\left(-\frac{\partial\tilde{u}}{\partial t}\delta u-\frac{\partial\tilde{u}}{\partial x_\alpha}\delta q_\alpha\right)d\Omega dt \\
    \delta R_q=&\int_\mathcal{I}\int_\mathcal{O}\tilde{q}_\alpha\left(\delta q_\alpha+\delta k\frac{\partial u}{\partial x_\alpha}+k\frac{\partial\delta u}{\partial x_\alpha}\right)d\Omega dt \notag\\
    =&\int_\mathcal{I}\int_{\partial\mathcal{O}}n_\alpha\tilde{q}_\alpha k\delta ud\Gamma dt+\int_\mathcal{I}\int_\mathcal{O}\left(\tilde{q}_\alpha\delta q_\alpha+\tilde{q}_\alpha\delta k\frac{\partial u}{\partial x_\alpha}-\frac{\partial\left(\tilde{q}_\alpha k\right)}{\partial x_\alpha}\delta u\right)d\Omega dt \\
    \delta R_{db}=&\int_\mathcal{I}\int_{\Gamma_\text{D}}\tilde{\lambda}\delta ud\Gamma dt \\
    \delta R_{nb}=&\int_\mathcal{I}\int_{\Gamma_\text{N}}\tilde{\lambda}n_\alpha\delta q_\alpha dt
\end{align}
$$
Therefore,
$$
\begin{align}
    &-\frac{\partial\tilde{u}}{\partial t}=\frac{\partial\left(\tilde{q}_\alpha k\right)}{\partial x_\alpha}-\frac{\partial J}{\partial u} &&\text{in }D \\
    &\tilde{q}_\alpha=\frac{\partial\tilde{u}}{\partial x_\alpha}-\frac{\partial J}{\partial q_\alpha} &&\text{in }D
\end{align}
$$
On Dirichlet boundary,
$$
\begin{align}
    &\tilde{\lambda}+n_\alpha\tilde{q}_\alpha k=0 \\
    &n_\alpha\tilde{u}+\frac{\partial J}{\partial q_\alpha}=0
\end{align}
$$
And on Neumann boundary,
$$
\begin{align}
    &n_\alpha\tilde{\lambda}+n_\alpha\tilde{u}=0 \\
    &n_\alpha\tilde{q}_\alpha k+\frac{\partial J}{\partial u}=0
\end{align}
$$

Same as the forward problem, the adjoint equations are discretized as follows:
$$
\begin{align}
    &\begin{aligned}
    -\int_V\frac{\partial\tilde{u}}{\partial t}dV=&\int_V\frac{\partial\left(\tilde{q}_\alpha k\right)}{\partial x_\alpha}dV-\int_V\frac{\partial J}{\partial u}dV \\
    =&\int_Sn_\alpha\tilde{q}_\alpha kdS-\int_V\frac{\partial J}{\partial u}dV
    \end{aligned} \\
    \therefore\:&-\frac{\tilde{u}^t-\tilde{u}^{t-1}}{\Delta t}\Delta x\Delta y=-\left(\tilde{q}_xk\right)_\text{W}\Delta y+\left(\tilde{q}_xk\right)_\text{E}\Delta y-\left(\tilde{q}_yk\right)_\text{S}\Delta x+\left(\tilde{q}_yk\right)_\text{N}\Delta x-\frac{\partial J}{\partial u}\Delta x\Delta y \\
    \therefore\:&\tilde{u}^{t-1}=\tilde{u}^t-\Delta t\left\{\frac{\left(-\tilde{q}_xk\right)_\text{E}-\left(-\tilde{q}_xk\right)_\text{W}}{\Delta x}+\frac{\left(-\tilde{q}_yk\right)_\text{N}-\left(-\tilde{q}_yk\right)_\text{S}}{\Delta y}+\frac{\partial J}{\partial u}\right\}
\end{align}
$$

For flux:
$$
\begin{align}
    &-k^L\tilde{q}_x^L=-k^L\frac{\tilde{u}^f-\tilde{u}^L}{\Delta x/2} \\
    &-k^R\tilde{q}_x^R=-k^R\frac{\tilde{u}^R-\tilde{u}^f}{\Delta x/2} \\
    &-k^L\tilde{q}_x^L=-k^R\tilde{q}_x^R \\
    \therefore\:&-k^L\frac{\tilde{u}^f-\tilde{u}^L}{\Delta x/2}=-k^R\frac{\tilde{u}^R-\tilde{u}^f}{\Delta x/2} \\
    \therefore\:&\tilde{u}^f=\frac{k^Lu^L+k^Ru^R}{k^L+k^R} \\
    \therefore\:&-k\tilde{q}_x=-k^L\tilde{q}_x^L=-k^R\tilde{q}_x^R=-\frac{2k^Lk^R}{k^L+k^R}\frac{\tilde{u}^R-\tilde{u}^L}{\Delta x}
\end{align}
$$

For Dirichlet boundary condition on x min:
$$
\begin{equation}
    -k\tilde{q}_x=-k^c\frac{\tilde{u}^c-\tilde{u}^b}{\Delta x/2}
\end{equation}
$$

## Design sensitivity

$$
\begin{equation}
    \langle J^\prime,\delta\gamma\rangle=\int_\mathcal{I}\int_\mathcal{O}\tilde{q}_\alpha\delta k\frac{\partial u}{\partial x_\alpha}d\Omega dt
\end{equation}
$$

## Objective functional

$$
\begin{equation}
    J=\int_{\Gamma_\text{q}}u^bd\Gamma
\end{equation}
$$

$$
\begin{align}
    &\bar{q}_y=-k^c\frac{u^c-u^b}{\Delta x/2} \\
    \therefore\;&u^b=u^c+\frac{\Delta x}{2}\frac{\bar{q}_y}{k^c} \\
    \therefore\;&\delta u^b=\delta u^c-\frac{\Delta x}{2}\frac{\bar{q}_y}{k^{c^2}}\delta k^c
\end{align}
$$