'''
Método de Jacobi
'''

def jacobi(A, b, x0, epsilon=1e-3):
    n = len(A)
    x = x0[:]
    dados_iteracao = []

    while True:
        x_novo = [0.0] * n
        for i in range(n):
            s1 = sum(A[i][j] * x[j] for j in range(i))  # Soma dos elementos anteriores
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, n))  # Soma dos elementos seguintes
            x_novo[i] = (b[i] - s1 - s2) / A[i][i]  # Atualização de x_novo[i]

        # Critério de parada baseado na diferença máxima
        max_diff = max(abs(x_novo[i] - x[i]) for i in range(n))

        dados_iteracao.append((len(dados_iteracao) + 1, x_novo[0], x_novo[1],
                               x_novo[2], max_diff))  # Armazena os dados da iteração

        if max_diff < epsilon:
            break

        x = x_novo[:]

    return x_novo, dados_iteracao


# Matriz A e vetor b fornecidos
A = [[11.9, 0.0, 1.8], [0.0, 5.3, -1.8], [1.0, -1.0, -1.0]]
b = [15.0, 3.1, 0.0]

# Vetor inicial
x0 = [1.0, 1.0, 1.0]

# Resolução do sistema
solucao, dados_iteracao = jacobi(A, b, x0)

# Imprimir a tabela de iterações
print("\nTabela de Iterações:")
print(f"\n|{'Iteração':^10}|{'I1':^10}|{'I2':^10}|{'I3':^10}|{'Erro':^10}|")
print("-" * 56)
for dado in dados_iteracao:
    print(f"|{dado[0]:^10}|{dado[1]:^10.7f}|{dado[2]:^10.7f}|{dado[3]:^10.7f}|"
          f"{dado[4]:^10.7f}|")

# Imprimir a solução
print("\nSolução:")
print(f"I1 = {solucao[0]:.7f}")
print(f"I2 = {solucao[1]:.7f}")
print(f"I3 = {solucao[2]:.7f}")

