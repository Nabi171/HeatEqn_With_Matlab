import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 1.0          # length of rod
k = 1.0          # thermal diffusivity
N = 40           # number of cosine modes

# Spatial grid
x = np.linspace(0, L, 400)

# Time values
t_values = [0, 0.01, 0.05, 0.2]

# Initial condition (change as needed)
def u0(x):
    return np.cos(np.pi*x/L) + 0.5*np.cos(3*np.pi*x/L) + 1

# Compute Fourier coefficients
A0 = (1/L) * np.trapz(u0(x), x)
A = np.zeros(N)

for n in range(1, N+1):
    integrand = u0(x) * np.cos(n*np.pi*x/L)
    A[n-1] = (2/L) * np.trapz(integrand, x)

# Plot solution at different times
plt.figure(figsize=(8,5))

for t in t_values:
    u = A0 * np.ones_like(x)
    for n in range(1, N+1):
        u += A[n-1] * np.cos(n*np.pi*x/L) * np.exp(-k*(n*np.pi/L)**2 * t)
    plt.plot(x, u, label=f"t = {t}")

plt.xlabel("x")
plt.ylabel("u(x,t)")
plt.title("1D Heat Equation (Neumann Boundary Condition)")
plt.legend()
plt.grid(True)
plt.show()
