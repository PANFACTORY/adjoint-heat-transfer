# Cell-centered value

## Scatter form -> Gather form

```python
for i in range(nx):
    for j in range(ny):
        aqx[i + 1, j] += -dt / dx * au_tp1[i, j]
        aqx[i, j] += dt / dx * au_tp1[i, j]
```

```python
for i in range(1, nx + 1):
    for j in range(ny):
        aqx[i, j] += -dt / dx * au_tp1[i - 1, j]
        
for i in range(nx):
    for j in range(ny):
        aqx[i, j] += dt / dx * au_tp1[i, j]  
```

```python
for i in range(nx + 1):
    for j in range(ny):
        if i > 0:
            aqx[i, j] += -dt / dx * au_tp1[i - 1, j]
        if i < nx:
            aqx[i, j] += dt / dx * au_tp1[i, j]
```

## Random access -> numpy style

# Compute y flux

## Scatter form -> Gather form

```python
for i in range(nx):
    for j in range(1, ny):
        k_face = 2 * k[i, j - 1] * k[i, j] / (k[i, j - 1] + k[i, j])
        # qy[i, j] = -k_face * (u_t[i, j] - u_t[i, j - 1]) / dy
        au_t[i, j - 1] += k_face / dy * aqy[i, j]
```

```python
for i in range(nx):
    for j in range(ny - 1):
        k_face = 2 * k[i, j] * k[i, j + 1] / (k[i, j] + k[i, j + 1])
        # qy[i, j] = -k_face * (u_t[i, j] - u_t[i, j - 1]) / dy
        au_t[i, j] += k_face / dy * aqy[i, j + 1]
```

```python
for i in range(nx):
    for j in range(ny):
        if j > 0:
            k_face = 2 * k[i, j - 1] * k[i, j] / (k[i, j - 1] + k[i, j])
            au_t[i, j] += -k_face / dy * aqy[i, j]
        if j < ny - 1:
            k_face = 2 * k[i, j] * k[i, j + 1] / (k[i, j] + k[i, j + 1])
            au_t[i, j] += k_face / dy * aqy[i, j + 1]
```

```python
for i in range(nx):
    # boundaries
    au_t[i, ny - 1] += k[i, ny - 1] / (0.5 * dy) * aqy[i, ny]

    # inner faces
    for j in range(ny):
        if j > 0:
            k_face = 2 * k[i, j - 1] * k[i, j] / (k[i, j - 1] + k[i, j])
            au_t[i, j] += -k_face / dy * aqy[i, j]
        if j < ny - 1:
            k_face = 2 * k[i, j] * k[i, j + 1] / (k[i, j] + k[i, j + 1])
            au_t[i, j] += k_face / dy * aqy[i, j + 1]
```

## Random access -> numpy style

# Compute thermal diffusivity

## Scatter form -> Gather form

```python
for i in range(nx):
    for j in range(1, ny):
        ak_face = -(u_t[i, j] - u_t[i, j - 1]) / dy * aqy[i, j]
        ak[i, j - 1] += 2 * (k[i, j] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face
        ak[i, j] += 2 * (k[i, j - 1] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face
```

```python
for i in range(nx):
    for j in range(1, ny):
        ak_face = -(u_t[i, j] - u_t[i, j - 1]) / dy * aqy[i, j]
        ak[i, j - 1] += 2 * (k[i, j] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face
        
for i in range(nx):
    for j in range(1, ny):
        ak_face = -(u_t[i, j] - u_t[i, j - 1]) / dy * aqy[i, j]
        ak[i, j] += 2 * (k[i, j - 1] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face
```

```python
for i in range(nx):
    for j in range(ny):
        if j < ny - 1:
            ak_face = -(u_t[i, j + 1] - u_t[i, j]) / dy * aqy[i, j + 1]
            ak[i, j] += 2 * (k[i, j + 1] / (k[i, j] + k[i, j + 1])) ** 2 * ak_face
        if j > 0:
            ak_face = -(u_t[i, j] - u_t[i, j - 1]) / dy * aqy[i, j]
            ak[i, j] += 2 * (k[i, j - 1] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face
```