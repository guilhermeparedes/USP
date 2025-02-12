'''
Método de Gauss-Seidel
'''

def gauss_seidel(A, b, x0, epsilon=1e-3):
    n = len(A)  # Dimensão da matriz
    x = x0[:]  # Vetor de aproximação inicial
    dados_iteracao = []  # Para armazenar os dados das iterações

    while True:
        x_ant = x[:] # Loop infinito até que o critério de parada seja atendido
        for i in range(n):
            s1 = sum(A[i][j] * x[j] for j in range(i))  # Usa os valores já atualizados
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, n))  # Usa os valores antigos
            x[i] = (b[i] - s1 - s2) / A[i][i]  # Atualiza x[i]

        # Cálculo do critério de parada baseado na diferença máxima
        max_diff = max(abs(x[i] - x_ant[i]) for i in range(n))  # Calcula a diferença máxima
        dados_iteracao.append((len(dados_iteracao) + 1, x[0], x[1], x[2], max_diff))  # Armazena os dados da iteração

        if max_diff < epsilon:  # Verifica o critério de parada
            break

         # Atualiza o vetor original para o próximo ciclo

    return x, dados_iteracao

# Matriz A e vetor b fornecidos
A = [[11.9, 0.0, 1.8], [0.0, 5.3, -1.8], [1.0, -1.0, -1.0]]
b = [15.0, 3.1, 0.0]

# Vetor inicial
x0 = [1.0, 1.0, 1.0]

# Resolução do sistema
solucao, dados_iteracao = gauss_seidel(A, b, x0)

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