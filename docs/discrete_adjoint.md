# Discrete adjoint variable method

The adjoint variable method is one of the sensitivity analysis methods, and its feature is efficiency for a huge number of design variables.
In this method, it is only required to obtain the design sensitivity for all design variables with the same order of computational cost as one forward analysis. 
There are two types of the adjoint variable method: (a) continuous form (first optimize, then discretize), and (b) discrete form (first discretize, then optimize).
In this manuscript, we explain (b) the discrete form. 

## Formulation

We can interpret the FVM as a matrix operation that maps $u^t$ to $u^{t+1}$ using $k$, $q_\alpha$, $\bar{u}$, $\bar{q}$ and $u_0$. This is represented as follows:
$$
\begin{align}
    u^0&=u_0 \\
    q^t_\alpha&=\left(\begin{array}{c|c|c|c}M_{\alpha,u^{t-1}}&M_{\alpha,\bar{u}}&M_{\alpha,\bar{q}}&M_{\alpha,k}\end{array}\right)\left(\begin{array}{l}u^{t-1}\\\bar{u}\\\bar{q}\\k\end{array}\right) \\
    u^t&=\left(\begin{array}{c|c}N_{u^{t-1}}&N_{q^t_\alpha}\end{array}\right)\left(\begin{array}{l}u^{t-1}\\q^t_\alpha\end{array}\right)
\end{align}
$$ The residual is given by
$$
\begin{align}
    R=&\tilde{u}^0\left(u^0-u_0\right) \notag\\
    &+\sum_{t=1}^{N_\text{t}}\tilde{q}_\alpha^t\left\{q^t_\alpha-\left(\begin{array}{c|c|c|c}M^\alpha_{u^{t-1}}&M^\alpha_{\bar{u}}&M^\alpha_{\bar{q}}&M^\alpha_k\end{array}\right)\left(\begin{array}{l}u^{t-1}\\\bar{u}\\\bar{q}\\k\end{array}\right)\right\} \notag\\
    &+\sum_{t=1}^{N_\text{t}}\tilde{u}^t\left\{u^t-\left(\begin{array}{c|c}N_{u^{t-1}}&N_{q^t_\alpha}\end{array}\right)\left(\begin{array}{l}u^{t-1}\\q^t_\alpha\end{array}\right)\right\}
\end{align}
$$ Therefore, the variation of $R$ is given by
$$
\begin{align}
    \delta R=&\tilde{u}^0\delta u^0 \notag\\
    &+\sum_{t=1}^{N_\text{t}}\tilde{q}_\alpha^t\left\{\delta q^t_\alpha-\left(\begin{array}{c|c}M_{\alpha,u^{t-1}}&M_{\alpha,k}\end{array}\right)\left(\begin{array}{l}\delta u^{t-1}\\\delta k\end{array}\right)\right\} \notag\\
    &+\sum_{t=1}^{N_\text{t}}\tilde{u}^t\left\{\delta u^t-\left(\begin{array}{c|c}N_{u^{t-1}}&N_{q^t_\alpha}\end{array}\right)\left(\begin{array}{l}\delta u^{t-1}\\\delta q^t_\alpha\end{array}\right)\right\} \notag\\
    =&\tilde{u}^{N_\text{t}}\delta u^{N_\text{t}} \notag\\
    &+\sum_{t=1}^{N_\text{t}-1}\left\{\tilde{u}^{t-1}-\left(\begin{array}{c|c}M_{\alpha,u^{t-1}}&N_{u^{t-1}}\end{array}\right)^\text{T}\left(\begin{array}{c}\tilde{q}^t_\alpha\\\tilde{u}^t\end{array}\right)\right\}\delta u^{t-1} \notag\\
    &+\sum_{t=1}^{N_\text{t}}\left\{\tilde{q}_\alpha^t-N_{q_\alpha^t}^\text{T}\tilde{u}^t\right\}\delta q_\alpha^t \notag\\
    &-\sum_{t=1}^{N_\text{t}}\tilde{q}^t_\alpha M_{\alpha,k}^T\delta k
\end{align}
$$ where the superscript $T$ denotes the transpose. The design sensitivity of a general functional $J$ is given by
$$
\begin{equation}
    J^\prime=-\sum_{t=1}^{N_\text{t}}\tilde{q}^t_\alpha M_k^{\alpha^T}\delta k
\end{equation}
$$ and the corresponding adjoint equations are
$$
\begin{align}
    &\tilde{u}^{N_\text{t}}=-\frac{\partial J}{\partial u^{N_\text{t}}} \\
    &\tilde{u}^{t-1}=\left(\begin{array}{c|c}M_{\alpha,u^{t-1}}&N_{u^{t-1}}\end{array}\right)^\text{T}\left(\begin{array}{c}\tilde{q}^t_\alpha\\\tilde{u}^t\end{array}\right)-\frac{\partial J}{\partial u^{t-1}} \\
    &\tilde{q}_\alpha^t=N_{q_\alpha^t}^\text{T}\tilde{u}^t-\frac{\partial J}{\partial q_\alpha^t}
\end{align}
$$ The adjoint equations are solved sequentially backward in time, whereas the governing equations are computed sequentially forward in time.

## Implementation

Here, the procedure to derive the implementation of adjoint equations from governing equations are explained. This consists three steps:

- Governing equation -> Adjoint equation
- Scatter form -> Gather form
- Random access -> numpy style

### Cell-centered value

#### Governing equation -> Adjoint equation

```python
for i in range(nx):
    for j in range(ny):
        aqx[i + 1, j] += -dt / dx * au_t[i, j]
        aqx[i, j] += dt / dx * au_t[i, j]
```

#### Scatter form -> Gather form

```python
for i in range(1, nx + 1):
    for j in range(ny):
        aqx[i, j] += -dt / dx * au_t[i - 1, j]
        
for i in range(nx):
    for j in range(ny):
        aqx[i, j] += dt / dx * au_t[i, j]  
```

If necessary, they can be combined into a single for loop sets.

```python
for i in range(nx + 1):
    for j in range(ny):
        if i > 0:
            aqx[i, j] += -dt / dx * au_t[i - 1, j]
        if i < nx:
            aqx[i, j] += dt / dx * au_t[i, j]
```

#### Random access -> numpy style

```python
aqx[1 : nx + 1, :] += -dt / dx * au_t
aqx[:nx, :] += dt / dx * au_t
```

### Y flux

#### Governing equation -> Adjoint equation

```python
for i in range(nx):
    for j in range(1, ny):
        k_face = 2 * k[i, j - 1] * k[i, j] / (k[i, j - 1] + k[i, j])
        # qy[i, j] = -k_face * (u_t[i, j] - u_t[i, j - 1]) / dy
        au_tm1[i, j - 1] += k_face / dy * aqy[i, j]
```

#### Scatter form -> Gather form

```python
for i in range(nx):
    for j in range(ny - 1):
        k_face = 2 * k[i, j] * k[i, j + 1] / (k[i, j] + k[i, j + 1])
        # qy[i, j] = -k_face * (u_t[i, j] - u_t[i, j - 1]) / dy
        au_tm1[i, j] += k_face / dy * aqy[i, j + 1]
```

If necessary, they can be combined into a single for loop sets.

```python
for i in range(nx):
    for j in range(ny):
        if j > 0:
            k_face = 2 * k[i, j - 1] * k[i, j] / (k[i, j - 1] + k[i, j])
            au_tm1[i, j] += -k_face / dy * aqy[i, j]
        if j < ny - 1:
            k_face = 2 * k[i, j] * k[i, j + 1] / (k[i, j] + k[i, j + 1])
            au_tm1[i, j] += k_face / dy * aqy[i, j + 1]
```

#### Random access -> numpy style

```python
k_face = 2 * k[:, : ny - 1] * k[:, 1:ny] / (k[:, : ny - 1] + k[:, 1:ny])
au_tm1[:, 1:ny] += -k_face / dy * aqy[:, 1:ny]
au_tm1[:, : ny - 1] += k_face / dy * aqy[:, 1:ny]
```

### Thermal diffusivity

#### Governing equation -> Adjoint equation

```python
for i in range(nx):
    for j in range(1, ny):
        ak_face = -(u_tm1[i, j] - u_tm1[i, j - 1]) / dy * aqy[i, j]
        ak[i, j - 1] += 2 * (k[i, j] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face
        ak[i, j] += 2 * (k[i, j - 1] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face
```

#### Scatter form -> Gather form

```python
for i in range(nx):
    for j in range(1, ny):
        ak_face = -(u_tm1[i, j] - u_tm1[i, j - 1]) / dy * aqy[i, j]
        ak[i, j - 1] += 2 * (k[i, j] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face
        
for i in range(nx):
    for j in range(1, ny):
        ak_face = -(u_tm1[i, j] - u_tm1[i, j - 1]) / dy * aqy[i, j]
        ak[i, j] += 2 * (k[i, j - 1] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face
```

If necessary, they can be combined into a single for loop sets.

```python
for i in range(nx):
    for j in range(ny):
        if j < ny - 1:
            ak_face = -(u_tm1[i, j + 1] - u_tm1[i, j]) / dy * aqy[i, j + 1]
            ak[i, j] += 2 * (k[i, j + 1] / (k[i, j] + k[i, j + 1])) ** 2 * ak_face
        if j > 0:
            ak_face = -(u_tm1[i, j] - u_tm1[i, j - 1]) / dy * aqy[i, j]
            ak[i, j] += 2 * (k[i, j - 1] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face
```

#### Random access -> numpy style

```python
ak_face = -(u_tm1[:, 1:ny] - u_tm1[:, : ny - 1]) / dy * aqy[:, 1:ny]
ak[:, 1:ny] += 2 * (k[:, : ny - 1] / (k[:, : ny - 1] + k[:, 1:ny])) ** 2 * ak_face
ak[:, : ny - 1] += 2 * (k[:, 1:ny] / (k[:, : ny - 1] + k[:, 1:ny])) ** 2 * ak_face
```
