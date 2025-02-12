'''
precisão simples tabela
'''

# Método de Simpson
import numpy as np
import matplotlib.pyplot as plt

# Função a ser integrada
def f(x):
    return 6 - 6 * x ** 5  # A função deve retornar float32 automaticamente se x for float32

# Método de Simpson otimizado
def simpson_integral(a, b, N):
    h = np.float32((b - a) / N)
    x = np.linspace(a, b, N + 1, dtype=np.float32)
    y = f(x)

    # Usando operações vetoriais para calcular os termos ímpares e pares, garantindo float32
    odd_terms = np.sum(y[1:N:2] * 4).astype(np.float32)  # Garante float32
    even_terms = np.sum(y[2:N - 1:2] * 2).astype(np.float32)  # Multiplicação em float32

    # Integral final
    integral = np.float32((y[0] + odd_terms + even_terms + y[N]) * (h / 3))
    return integral

# Parâmetros
a = 0
b = 1
p_valores = np.arange(1, 26)
erros_simples = []

'''
# Tabela de resultados
print("-" * 55)
print(f"|{'p':^10}|{'N':^10}|{'I_num':^15}|{'Erro':^15}|")
print("-" * 55)


# Cálculo das integrais e erros
for p in p_valores:
    N = 2 ** p

    # Cálculo com precisão simples (float32)
    I_num_simples = simpson_integral(a, b, N)
    erro_simples = abs(I_num_simples - I_analitico)
    erros_simples.append((p, N, I_num_simples, erro_simples))

    print(f"|{p:^10}|{N:^10}|{I_num_simples:^15.8f}|{erro_simples:^15.8f}|")

print("-" * 55)
'''

'''
precisão dupla tabela
'''

# Função a ser integrada
def f(x):
    return 6 - 6 * x ** 5

# Integral analítica
I_analitico = 5

# Método de Simpson
def simpson_integral(a, b, N):
    h = (b - a) / N
    x = np.linspace(a, b, N + 1)
    y = f(x)

    # Usando operações vetoriais para calcular os termos ímpares e pares, garantindo float32
    odd_terms = np.sum(y[1:N:2] * 4)  # Garante float32
    even_terms = np.sum(y[2:N - 1:2] * 2)  # Multiplicação em float32

    # Integral final
    integral = (y[0] + odd_terms + even_terms + y[N]) * (h / 3)
    return integral

# Parâmetros
a = 0
b = 1
p_valores = np.arange(1, 26)  

#erros_dupla = []

'''

# Tabela de resultados
print("-"*55)
print(f"|{'p':^10}|{'N':^10}|{'I_num':^15}|{'Erro':^15}|")
print("-"*55)


# Cálculo das integrais e erros
for p in p_valores:
    N = 2 ** p

    # Cálculo com precisão dupla (float64)
    I_num_dupla = simpson_integral(a, b, N)
    erro_dupla = abs(I_num_dupla - I_analitico)
    erros_dupla.append((p, N, I_num_dupla, erro_dupla))
    print(f"|{p:^10}|{N:^10}|{I_num_dupla:^15.8f}|{erro_dupla:^15.8f}|")

print("-"*55)
'''

# Dados para o gráfico
p_values_simple = [p for p, N, I_num, erro in erros_simples if erro > 0]
erro_values_simple = [np.log2(erro) for p, N, I_num, erro in erros_simples if erro > 0]

#p_values_dupla = [p for p, N, I_num, erro in erros_dupla if erro > 0]
#erro_values_dupla = [np.log2(erro) for p, N, I_num, erro in erros_dupla if erro > 0]


# Dividindo os dados em duas regiões para precisão simples
p_values_region1_simple = p_values_simple[:12]
erro_values_region1_simple = erro_values_simple[:12]

p_values_region2_simple = p_values_simple[12:]
erro_values_region2_simple = erro_values_simple[12:]

# Ajuste linear para a região 1 (precisão simples)
alpha1_simple, logA1_simple = np.polyfit(p_values_region1_simple, erro_values_region1_simple, 1)

# Ajuste linear para a região 2 (precisão simples)
alpha2_simple, logA2_simple = np.polyfit(p_values_region2_simple, erro_values_region2_simple, 1)

# Dividindo os dados em duas regiões para precisão dupla
#p_values_region1_dupla = p_values_dupla[:12]
#erro_values_region1_dupla = erro_values_dupla[:12]

#p_values_region2_dupla = p_values_dupla[9:]
#erro_values_region2_dupla = erro_values_dupla[9:]

# Ajuste linear para a região 1 (precisão dupla)
#alpha1_dupla, logA1_dupla = np.polyfit(p_values_region1_dupla, erro_values_region1_dupla, 1)

# Ajuste linear para a região 2 (precisão dupla)
#alpha2_dupla, logA2_dupla = np.polyfit(p_values_region2_dupla, erro_values_region2_dupla, 1)


# Gráfico dos erros com ajustes lineares para precisão simples e dupla
plt.figure(figsize=(12, 8))


# Erros de precisão simples
plt.scatter(p_values_simple, erro_values_simple, label='Erro (Simples)', color='blue', marker='o')
plt.plot(p_values_region1_simple, alpha1_simple * np.array(p_values_region1_simple) + logA1_simple, color='orange', linestyle='-', label=f'Ajuste Região 1 (Simples): α = {-alpha1_simple:.4f}')
plt.plot(p_values_region2_simple, alpha2_simple * np.array(p_values_region2_simple) + logA2_simple, color='red', linestyle='-', label=f'Ajuste Região 2 (Simples): α = {-alpha2_simple:.4f}')

# Erros de precisão dupla
#plt.scatter(p_values_dupla, erro_values_dupla, label='Erro (Dupla)', color='purple', marker='x')
#plt.plot(p_values_region1_dupla, alpha1_dupla * np.array(p_values_region1_dupla) + logA1_dupla, color='green', linestyle='--', label=f'Ajuste Região 1 (Dupla): α = {-alpha1_dupla:.4f}')
#plt.plot(p_values_region2_dupla, alpha2_dupla * np.array(p_values_region2_dupla) + logA2_dupla, color='brown', linestyle='--', label=f'Ajuste Região 2 (Dupla): α = {-alpha2_dupla:.4f}')

# Personalização do gráfico
plt.xlabel('p')
plt.ylabel('$log_{2}(erro)$')
plt.axhline(0, color='black', linewidth=0.5, linestyle='--')
plt.legend()
plt.grid(True)
plt.title('Comparação de Ajustes Lineares para Erros de Simpson em Precisões Simples e Dupla')
plt.show()