'''
import numpy as np
import matplotlib.pyplot as plt


omega = 1

# Definindo a equação de Duffing com amortecimento e força externa
def f_duffing_forced(t, x, v, F, omega):
    dxdt = v
    dvdt = 0.5 * x * (1 - 4 * x ** 2) + F * np.cos(omega * t) - 0.25 * v
    return dxdt, dvdt

# Método Runge-Kutta de 4ª ordem
def rk4_step(t, x, v, F, omega, h):
    k1x, k1v = f_duffing_forced(t, x, v, F, omega)
    k2x, k2v = f_duffing_forced(t + h / 2, x + k1x * h / 2, v + k1v * h / 2, F, omega)
    k3x, k3v = f_duffing_forced(t + h / 2, x + k2x * h / 2, v + k2v * h / 2, F, omega)
    k4x, k4v = f_duffing_forced(t + h, x + k3x * h, v + k3v * h, F, omega)

    x_new = x + (h / 6) * (k1x + 2 * k2x + 2 * k3x + k4x)
    v_new = v + (h / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)

    return x_new, v_new


# Função principal para calcular o diagrama de bifurcação
def bifurcation_diagram(omega=1.0, max_F=0.35, delta_F=0.00025, h=0.001 * 2 * np.pi / omega, trans_steps=200000, per_steps=1000, per_count=100):
    # Faixa de F de 0 até max_F com incremento delta_F
    F_values = np.arange(0.1, max_F + delta_F, delta_F)
    x_vals = []

    # Loop sobre os valores de F
    for F in F_values:
        t = 0.0
        x = -0.5  # Condição inicial para x
        v = 0.5   # Condição inicial para v

        # Evoluindo o transiente
        for i in range(trans_steps):
            x, v = rk4_step(t, x, v, F, omega, h)
            t += h

        # Coletando pontos após o transiente
        for i in range(per_count):
            for j in range(per_steps):
                x, v = rk4_step(t, x, v, F, omega, h)
                t += h
            print(F, x)
            # Armazenando o valor de x para cada período
            x_vals.append(x)

    # Convertendo para um array numpy para facilitar o plot
    x_vals = np.array(x_vals)

    # Plotando o gráfico de bifurcação
    plt.figure(figsize=(8, 6))
    for i, F in enumerate(F_values):
        plt.scatter([F] * per_count, x_vals[i * per_count: (i + 1) * per_count], color='black', s=0.1)
    plt.title("Diagrama de Bifurcação")
    plt.xlabel("F ")
    plt.ylabel("x")
    plt.show()

# Chamada da função
bifurcation_diagram()
'''

import numpy as np
import matplotlib.pyplot as plt

# Definindo os parâmetros fora das funções
omega = 1
max_F = 0.11
delta_F = 0.00025
h = 0.001 * 2 * np.pi / omega  # Passo temporal
trans_steps = 200000
per_steps = 1000
per_count = 100
x0 = -0.5  # Condição inicial para x
v0 = 0.5   # Condição inicial para v
t0 = 0.0   # Tempo inicial

# Potencial poço duplo com força externa
def duffing_forced(t, x, v, forcing=None, frequency=None):
    return v, 0.5 * x * (1 - 4 * x ** 2) + forcing * np.cos(frequency * t) - 0.25 * v

# Método de Runge-Kutta de 4ª ordem
def rk4(t_ini, x_ini, v_ini, h, n, f, forcing=None, frequency=None):
    t = t_ini
    x = x_ini
    v = v_ini
    result = [(t, x, v)]

    for _ in range(n):
        k1x, k1v = f(t, x, v, forcing, frequency)
        k2x, k2v = f(t + h / 2, x + h * k1x / 2, v + h * k1v / 2, forcing, frequency)
        k3x, k3v = f(t + h / 2, x + h * k2x / 2, v + h * k2v / 2, forcing, frequency)
        k4x, k4v = f(t + h, x + h * k3x, v + h * k3v, forcing, frequency)

        x_next = x + (h / 6) * (k1x + 2 * k2x + 2 * k3x + k4x)
        v_next = v + (h / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)

        x, v, t = x_next, v_next, t + h
        result.append((t, x, v))

    return np.array(result)

# Função para remover o transiente (os primeiros 500 pontos)
def remove_transient(data):
    transient_cutoff = int(0.5 * len(data))
    return data[transient_cutoff:]

# Função principal para calcular o diagrama de bifurcação
def bifurcation_diagram():
    # Faixa de F de 0 até max_F com incremento delta_F
    F_values = np.arange(0.1, max_F + delta_F, delta_F)
    x_vals = []

    # Loop sobre os valores de F
    for F in F_values:
        t = t0
        x = x0  # Condição inicial para x
        v = v0  # Condição inicial para v

        # Evoluindo o transiente
        result = rk4(t, x, v, h, trans_steps, duffing_forced, forcing=F, frequency=omega)

        # Coletando pontos após o transiente
        for i in range(per_count):
            for j in range(per_steps):
                result = rk4(t, x, v, h, 1, duffing_forced, forcing=F, frequency=omega)
                t += h
            print(F, x)
            # Armazenando o valor de x para cada período
            x_vals.append(result[-1, 1])

    # Convertendo para um array numpy para facilitar o plot
    x_vals = np.array(x_vals)

    # Plotando o gráfico de bifurcação
    plt.figure(figsize=(8, 6))
    for i, F in enumerate(F_values):
        plt.scatter([F] * per_count, x_vals[i * per_count: (i + 1) * per_count], color='black', s=0.1)
    plt.title("Diagrama de Bifurcação")
    plt.xlabel("F ")
    plt.ylabel("x")
    plt.show()

# Chamada da função
bifurcation_diagram()

