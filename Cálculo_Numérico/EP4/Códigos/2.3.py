import numpy as np
import matplotlib.pyplot as plt

omega = 1


def f_duffing_forced(t, x, v, F, omega):
    dxdt = v
    dvdt = 0.5 * x * (1 - 4 * x ** 2) + F * np.cos(omega * t) - 0.25 * v
    return dxdt, dvdt


def rk4_step(t, x, v, F, omega, h):
    k1x, k1v = f_duffing_forced(t, x, v, F, omega)
    k2x, k2v = f_duffing_forced(t + h / 2, x + k1x * h / 2, v + k1v * h / 2, F, omega)
    k3x, k3v = f_duffing_forced(t + h / 2, x + k2x * h / 2, v + k2v * h / 2, F, omega)
    k4x, k4v = f_duffing_forced(t + h, x + k3x * h, v + k3v * h, F, omega)

    x_new = x + (h / 6) * (k1x + 2 * k2x + 2 * k3x + k4x)
    v_new = v + (h / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)

    return x_new, v_new


def bifurcation_diagram(omega=1.0, F=0.26, h=0.001 * 2 * np.pi / omega, trans_steps=200000, per_steps=1000,
                        per_count=3000):
    t = 0.0
    x = -0.5
    v = 0.5

    x_values = []
    v_values = []

    # Evoluindo o transiente
    for i in range(trans_steps):
        x, v = rk4_step(t, x, v, F, omega, h)
        t += h

    # Coletando pontos após o transiente
    for i in range(per_count):
        for j in range(per_steps):
            x, v = rk4_step(t, x, v, F, omega, h)
            t += h
        # Armazenando o valor de x e v para cada período
        x_values.append(x)
        v_values.append(v)

    # Convertendo para arrays numpy para facilitar o plot
    x_vals = np.array(x_values)
    v_vals = np.array(v_values)

    # Plotando o gráfico no plano de fase (x(t) vs. v(t))
    plt.figure(figsize=(8, 6))
    plt.scatter(x_vals, v_vals, color='black', s=0.1)
    plt.xlabel("x(t)")
    plt.ylabel("v(t)")
    plt.show()


# Chamada da função para F = 0.26
bifurcation_diagram()
