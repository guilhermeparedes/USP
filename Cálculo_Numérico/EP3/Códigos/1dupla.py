'''
Cálculo da integral com precisão simples usando Simpson e análise de erros
'''

import numpy as np
import matplotlib.pyplot as plt


# Função a ser integrada
def f(x):
    return np.float64(6 - 6 * x ** 5)


I_analitico = np.float64(5)


# Função para o método de Simpson
def simpson_integral(a, b, N):
    # Passo de integração em precisão simples
    h = np.float64((b - a) / N)
    x_values = [np.float64(a + i * h) for i in range(N + 1)]
    y_values = [f(x) for x in x_values]

    # Cálculo da integral usando método de Simpson
    integral = np.float64(0)
    for i in range(1, N, 2):  # Termos ímpares
        integral += np.float64(4) * y_values[i]
    for i in range(2, N - 1, 2):  # Termos pares
        integral += np.float64(2) * y_values[i]

    # Acrescenta os extremos
    integral += y_values[0] + y_values[N]
    integral *= h / np.float64(3)
    return integral

# Parâmetros
a, b = np.float64(0), np.float64(1)
p_values = np.arange(1, 26, dtype=int)

# Tabela de resultados e cálculo dos erros
print("-" * 55)
print(f"|{'p':^10}|{'N':^10}|{'I_num':^15}|{'Erro':^15}|")
print("-" * 55)

# Cálculo da integral e análise de erro ponto a ponto
params = []
for p in p_values:
    N = 2 ** p
    I_num = simpson_integral(a, b, N)

    erro = np.float64(abs(I_num - I_analitico))

    params.append((p, N, I_num, erro))
    print(f"|{p:^10}|{N:^10}|{I_num:^15.8f}|{erro:^15.8f}|")

print("-" * 55)

'''
# Dados para o gráfico
p_values_simple = [p for p, N, I_num, erro in params if erro > 0]
erro_values_simple = [np.log2(erro) for p, N, I_num, erro in params if erro > 0]

# Divisão em duas regiões para análise visual
p_values_region1 = p_values_simple[:11]
erro_values_region1 = erro_values_simple[:11]
p_values_region2 = p_values_simple[12:18]
erro_values_region2 = erro_values_simple[12:18]

# Ajuste linear para a primeira região
alpha1, logA1 = np.polyfit(p_values_region1, erro_values_region1, 1)

# Ajuste linear para a segunda região
alpha2, logA2 = np.polyfit(p_values_region2, erro_values_region2, 1)

# Gráficos dos erros
plt.figure(figsize=(10, 6))
plt.scatter(p_values_simple, erro_values_simple, label='erro (dupla)', marker='o', color='blue')

# Plotando os ajustes lineares
plt.plot(p_values_region1, alpha1 * np.array(p_values_region1) + logA1, color='orange', label=f'ajuste linear alpha = {-alpha1:.4f}')
plt.plot(p_values_region2, alpha2 * np.array(p_values_region2) + logA2, color='red', label=f'ajuste linear alpha = {-alpha2:.4f}')

# Configurações do gráfico
plt.xlabel('p')
plt.ylabel('$\log_{2}(\mathrm{erro})$')
plt.axhline(0, color='black', linewidth=0.5, linestyle='--')
plt.legend()
plt.show()
'''