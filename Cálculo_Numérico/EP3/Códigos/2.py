import numpy as np
import matplotlib.pyplot as plt

# Constantes
l = 1.0
g = 9.81

# Função para calcular o valor da integral usando o método dos trapézios
def trapezoidal_integration(a, b, n, k):
    h = (b - a) / n
    integral = 0.5 * (integrando(a, k) + integrando(b, k))

    for i in range(1, n):
        integral += integrando(a + i * h, k)

    return integral * h

# Função para a integral
def integrando(phi, k):
    term = 1 - k ** 2 * np.sin(phi) ** 2
    if np.isclose(term, 0):     # Evitar divisão por zero
        return 0
    return 1 / np.sqrt(term)


# Função principal para calcular T/TGalileu
def razao_periodos(num_theta_valores):
    theta_valores = np.linspace(0, np.pi, num_theta_valores)
    T_razao = []

    for theta0 in theta_valores:
        k_quadrado = (1 - np.cos(theta0)) / 2
        integral_value = trapezoidal_integration(0, np.pi / 2, 1000, k_quadrado)
        T = 4 * np.sqrt(l / g) * integral_value
        T_Galileu = 2 * np.pi * np.sqrt(l / g)
        T_razao.append(T / T_Galileu)

    return theta_valores, T_razao


# Parâmetros
num_theta_valores = 200
theta_valores, T_razao = razao_periodos(num_theta_valores)

# Tabela com 10 valores de θ0
print("-" * 37)
print(f"| {'θ0 (rad)':^15} | {'T/T_Galileu':^15} |")
print("-" * 37)
# 10 valores igualmente espaçados utilizando fatiamento [inicio:fim:passo]
for theta, razao in zip(theta_valores[::num_theta_valores // 10],
                        T_razao[::num_theta_valores // 10]):
    print(f"| {theta:^15.8f} | {razao:^15.8f} |")
print("-" * 37)

plt.figure(figsize=(10, 6))
plt.plot(theta_valores, T_razao, label='$T/T_{Galileu}$', color='blue')
plt.xlabel('$θ_0$ [rad]')
plt.ylabel('$T/T_{Galileu}$')
plt.axhline(1, color='black', linewidth=0.5, linestyle='--')
plt.legend()
plt.show()
