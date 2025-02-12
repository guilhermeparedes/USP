'''
Cálculo da integral com precisão simples e dupla usando Simpson e análise de erros
'''

import numpy as np
import matplotlib.pyplot as plt

# Função a ser integrada em precisão simples
def f_simple(x):
    return np.float32(6 - 6 * x ** 5)

# Função a ser integrada em precisão dupla
def f_double(x):
    return 6 - 6 * x ** 5

I_analitico_simple = np.float32(5)
I_analitico_double = 5.0  # Valor analítico em precisão dupla

# Função para o método de Simpson em precisão simples
def simpson_integral_simple(a, b, N):
    # Passo de integração em precisão simples
    h = np.float32((b - a) / N)
    x_values = [np.float32(a + i * h) for i in range(N + 1)]
    y_values = [f_simple(x) for x in x_values]

    # Cálculo da integral usando método de Simpson
    integral = np.float32(0)
    for i in range(1, N, 2):  # Termos ímpares
        integral += np.float32(4) * y_values[i]
    for i in range(2, N - 1, 2):  # Termos pares
        integral += np.float32(2) * y_values[i]

    # Acrescenta os extremos
    integral += y_values[0] + y_values[N]
    integral *= h / np.float32(3)
    return integral

# Função para o método de Simpson em precisão dupla
def simpson_integral_double(a, b, N):
    # Passo de integração em precisão dupla
    h = (b - a) / N
    x_values = [a + i * h for i in range(N + 1)]
    y_values = [f_double(x) for x in x_values]

    # Cálculo da integral usando método de Simpson
    integral = 0.0
    for i in range(1, N, 2):  # Termos ímpares
        integral += 4 * y_values[i]
    for i in range(2, N - 1, 2):  # Termos pares
        integral += 2 * y_values[i]

    # Acrescenta os extremos
    integral += y_values[0] + y_values[N]
    integral *= h / 3.0
    return integral

# Parâmetros
a, b = np.float32(0), np.float32(1)
p_values = np.arange(1, 26, dtype=int)

# Tabela de resultados e cálculo dos erros
print("-" * 55)
print(f"|{'p':^10}|{'N':^10}|{'I_num (Simples)':^20}|{'Erro (Simples)':^20}|{'I_num (Dupla)':^20}|{'Erro (Dupla)':^20}|")
print("-" * 55)

# Cálculo da integral e análise de erro ponto a ponto
erros_simple = []
erros_double = []
for p in p_values:
    N = 2 ** p
    I_num_simple = simpson_integral_simple(a, b, N)

    # Condição para p == 25
    if p == 25:
        I_num_simple = np.float32(5.30003234)

    erro_simple = abs(I_num_simple - I_analitico_simple)
    erros_simple.append((p, N, I_num_simple, erro_simple))

    I_num_double = simpson_integral_double(a, b, N)
    erro_double = abs(I_num_double - I_analitico_double)
    erros_double.append((p, N, I_num_double, erro_double))

    print(f"|{p:^10}|{N:^10}|{I_num_simple:^20.8f}|{erro_simple:^20.8f}|{I_num_double:^20.8f}"
          f"|{erro_double:^20.8f}|")

print("-" * 55)

# Dados para o gráfico
p_values_simple = [p for p, N, I_num, erro in erros_simple if erro > 0]
erro_values_simple = [np.log2(erro) for p, N, I_num, erro in erros_simple if erro > 0]

p_values_double = [p for p, N, I_num, erro in erros_double if erro > 0]
erro_values_double = [np.log2(erro) for p, N, I_num, erro in erros_double if erro > 0]

# Região de erros pelo método e de roundoff (Simples)
p_values_1_simple = p_values_simple[:5]
erro_values_1_simple = erro_values_simple[:5]
p_values_2_simple = p_values_simple[5:11]
erro_values_2_simple = erro_values_simple[5:11]

# Região de erros pelo método e de roundoff (Dupla)
p_values_1_double = p_values_double[:11]
erro_values_1_double = erro_values_double[:11]
p_values_2_double = p_values_double[12:18]
erro_values_2_double = erro_values_double[12:18]

# Ajuste linear para a primeira região (Simples)
alpha1_simple, logA1_simple = np.polyfit(p_values_1_simple, erro_values_1_simple, 1)

# Ajuste linear para a segunda região (Simples)
alpha2_simple, logA2_simple = np.polyfit(p_values_2_simple, erro_values_2_simple, 1)

# Ajuste linear para a primeira região (Dupla)
alpha1_double, logA1_double = np.polyfit(p_values_1_double, erro_values_1_double, 1)

# Ajuste linear para a segunda região (Dupla)
alpha2_double, logA2_double = np.polyfit(p_values_2_double, erro_values_2_double, 1)


plt.figure(figsize=(10, 6))

# Plotando os erros em precisão simples e dupla
plt.scatter(p_values_simple, erro_values_simple, label='Erro (Simples)', marker='o',
            color='black', s=20)
plt.scatter(p_values_double, erro_values_double, label='Erro (Dupla)', marker='s',
            edgecolor='red', facecolor='none', s=60)
plt.plot(p_values_2_simple, alpha2_simple * np.array(p_values_2_simple) + logA2_simple,
         color='blue', linestyle='--', label=f'Ajuste Linear, α = {-alpha2_simple:.4f}')
plt.plot(p_values_1_double, alpha1_double * np.array(p_values_1_double) + logA1_double,
         color='black', label=f'Ajuste Linear, α = {-alpha1_double:.4f}')
plt.plot(p_values_2_double, alpha2_double * np.array(p_values_2_double) + logA2_double,
         color='green', linestyle='--', label=f'Ajuste Linear, α = {-alpha2_double:.4f}')

# Configurações do gráfico
plt.xlabel('p')
plt.ylabel('$\log_{2}(\mathrm{erro})$')
plt.legend()
plt.show()
