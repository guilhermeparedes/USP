import numpy as np


def euler(t_ini, y_ini, z_ini, h, tf):
    t, y, z = t_ini, y_ini, z_ini
    result = [(t, y, z)]

    for _ in range(tf):
        z_next = z + h * g(t, y, z)
        y_next = y + h * z
        y, z, t = y_next, z_next, t + h
        result.append((t, y, z))

    return np.array(result)


def g(t, y, z):
    return z + y - t ** 3 - 3 * t ** 2 + 7 * t + 1


def rk4(t_ini, y_ini, z_ini, h, n):
    t, y, z = t_ini, y_ini, z_ini
    result = [(t, y, z)]

    for _ in range(n):
        k1y = h * z
        k1z = h * g(t, y, z)

        k2y = h * (z + k1z / 2)
        k2z = h * g(t + h / 2, y + k1y / 2, z + k1z / 2)

        k3y = h * (z + k2z / 2)
        k3z = h * g(t + h / 2, y + k2y / 2, z + k2z / 2)

        k4y = h * (z + k3z)
        k4z = h * g(t + h, y + k3y, z + k3z)

        y_next = y + (k1y + 2 * k2y + 2 * k3y + k4y) / 6
        z_next = z + (k1z + 2 * k2z + 2 * k3z + k4z) / 6

        y, z, t = y_next, z_next, t + h
        result.append((t, y, z))

    return np.array(result)


t0 = 0
y0 = 0
z0 = -1
step = 0.01
N = 500

euler_result = euler(t0, y0, z0, step, N)
rk4_result = rk4(t0, y0, z0, step, N)

t_values = np.linspace(0, 5)
exact_values = exact_solution(t_values)

print(f"Método de Euler: y(5) = {euler_result[-1, 1]:.8f}, y'(5) = {euler_result[-1, 2]:.8f}")
print(f"Método de RK4: y(5) = {rk4_result[-1, 1]:.8f}, y'(5) = {rk4_result[-1, 2]:.8f}")
print(f"Solução exata: y(5) = {exact_solution(5)}, y'(5) = {exact_derivative(5)}")

erro = abs(rk4_result[-1, 2] - 74)

print(erro)

'''

import matplotlib.pyplot as plt


# Plotando as soluções
plt.plot(euler_result[:, 0], euler_result[:, 1], label='Euler - y(t)')
plt.plot(rk4_result[:, 0], rk4_result[:, 1], label='RK4 - y(t)')
plt.plot(t_values, exact_values, label='Exata - y(t)', linestyle='--')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.legend()
plt.title('Comparação de soluções: Euler, RK4 e Exata')
plt.grid(True)
plt.show()
'''