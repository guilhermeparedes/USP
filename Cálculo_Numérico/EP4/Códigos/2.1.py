import numpy as np
import matplotlib.pyplot as plt


# Potencial poco duplo
def duffing_potential(t, x, v, damping=None, forcing=None, frequency=None):
    return v, 0.5 * x * (1 - 4 * x ** 2)


# Potencial poço duplo com amortecimento
def duffing_damped(t, x, v, damping, forcing=None, frequency=None):
    return v, 0.5 * x * (1 - 4 * x ** 2) - 2 * damping * v


# Potencial poço duplo com forca externa
def duffing_forced(t, x, v, damping=None, forcing=None, frequency=None):
    return v, 0.5 * x * (1 - 4 * x ** 2) + forcing * np.cos(frequency * t) - 0.25 * v


def rk4(t_ini, x_ini, v_ini, h, n, f, damping=None, forcing=None, frequency=None):
    t = t_ini
    x = x_ini
    v = v_ini
    result = [(t, x, v)]

    for _ in range(n):
        k1x, k1v = f(t, x, v, damping, forcing, frequency)
        k2x, k2v = f(t + h / 2, x + h * k1x / 2, v + h * k1v / 2, damping, forcing, frequency)
        k3x, k3v = f(t + h / 2, x + h * k2x / 2, v + h * k2v / 2, damping, forcing, frequency)
        k4x, k4v = f(t + h, x + h * k3x, v + h * k3v, damping, forcing, frequency)

        x_next = x + (h / 6) * (k1x + 2 * k2x + 2 * k3x + k4x)
        v_next = v + (h / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)

        x, v, t = x_next, v_next, t + h
        result.append((t, x, v))

    return np.array(result)


t0 = 0
x0 = -0.5
v0_values = [0.1, 0.25, 0.5]
gamma_values = [0.25 / 2, 0.8 / 2]
F_values = [0.11, 0.115, 0.14, 0.35]
omega = 1
step = 0.01
N = 10000


# Função para remover o transiente (os primeiros 500 dos pontos)
def remove_transient(data):
    transient_cutoff = int(0.5 * len(data))
    return data[transient_cutoff:]


plt.figure(figsize=(8, 6))

# a) Potencial poço duplo
for v0 in v0_values:
    results = rk4(t0, x0, v0, step, N, duffing_potential)
    plt.plot(results[:, 1], results[:, 2], label=f'v0 = {v0}')

plt.xlabel('x(t)')
plt.ylabel('v(t)')
plt.title('Espaço de Fase - Potencial Poço Duplo')
plt.legend()
plt.show()

plt.figure(figsize=(8, 6))

# b) Amortecimento
for gamma in gamma_values:
    results = rk4(t0, x0, 0.5, step, N, duffing_damped, damping=gamma)
    plt.plot(results[:, 1], results[:, 2], label=f'$\gamma$ = {gamma}')

plt.xlabel('x(t)')
plt.ylabel('v(t)')
plt.title('Espaço de Fase - Potencial com Amortecimento')
plt.legend()
plt.show()

plt.figure(figsize=(8, 6))

# c) Força externa
for F in F_values:
    results = rk4(t0, x0, 0.5, step, N, duffing_forced, forcing=F, frequency=omega)
    result_no_transient = remove_transient(results)
    plt.plot(result_no_transient[:, 1], result_no_transient[:, 2], label=f'F = {F}')

plt.xlabel('x(t)')
plt.ylabel('v(t)')
plt.title('Espaço de Fase - Potencial com Força Externa')
plt.legend()
plt.show()



'''
Caso (a): Poço duplo sem amortecimento
Atrator: Ciclo limite.
O sistema oscila indefinidamente dentro de uma trajetória fechada no espaço de fase.
Isso ocorre porque não há perda de energia (sem amortecimento), então o sistema se mantém em movimento periódico ao longo do ciclo.
Caso (b): Poço duplo com amortecimento
Atrator: Ponto fixo.
O amortecimento remove energia do sistema, fazendo com que as oscilações diminuam gradualmente até o sistema atingir o equilíbrio no ponto de menor energia do potencial (o atrator estável).
Dependendo da intensidade do amortecimento (
𝛾
γ), o sistema pode oscilar antes de se estabilizar (amortecimento fraco) ou simplesmente decair sem oscilar (amortecimento forte).
Caso (c): Poço duplo com força externa
Atrator: Ciclo limite.
A força externa aplicada periodicamente impede que o sistema alcance o equilíbrio. Em vez disso, o sistema oscila continuamente em resposta à força externa.
No espaço de fase, isso se traduz em uma trajetória fechada estável que corresponde a oscilações periódicas no sistema, mas essas oscilações agora dependem da frequência ($\omega$) e amplitude (F) da força externa.
'''
