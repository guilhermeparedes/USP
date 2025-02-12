def jacobi(A, b, x0, epsilon=1e-3, max_iteracao=100):
    """
    Resolve o sistema linear Ax = b usando o método de Jacobi.
    Sem numpy
    :param A: Matriz dos coeficientes
    :param b: Vetor dos termos constantes
    :param x0: Vetor inicial
    :param epsilon: Critério de parada (tolerância)
    :param max_iteracao: Número máximo de iterações
    :return: Solução aproximada e tabela de iterações
    """

    n = len(A)
    x = x0[:]
    dados_iteracao = []

    for k in range(max_iteracao):
        x_novo = [0.0] * n
        for i in range(n):
            # Soma dos termos conhecidos do lado esquerdo da equação
            s1 = sum(A[i][j] * x[j] for j in range(i))
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            # Cálculo do novo valor de x[i] usando valores antigos de x
            x_novo[i] = (b[i] - s1 - s2) / A[i][i]

        # Cálculo do critério de parada
        max_diff = max(abs(x_novo[i] - x[i]) for i in range(n))
        dados_iteracao.append((k + 1, x_novo[0], x_novo[1], x_novo[2], max_diff))

        if max_diff < epsilon:
            break

        # Atualiza x com os valores da iteração atual
        x = x_novo[:]

    return x_novo, dados_iteracao


def gauss_seidel(A, b, x0, epsilon=1e-3, max_iteracao=100):
    """
    Resolve o sistema linear Ax = b usando o método de Gauss-Seidel.
    Sem numpy
    :param A: Matriz dos coeficientes
    :param b: Vetor dos termos constantes
    :param x0: Vetor inicial
    :param epsilon: Critério de parada (tolerância)
    :param max_iteracao: Número máximo de iterações
    :return: Solução aproximada e tabela de iterações
    """

    n = len(A)
    x = x0[:]
    dados_iteracao = []

    for k in range(max_iteracao):
        x_novo = x[:]  # Copia os valores atuais de x
        for i in range(n):
            # Soma dos termos conhecidos do lado esquerdo da equação
            s1 = sum(A[i][j] * x_novo[j] for j in range(i))
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            # Atualiza o valor de x[i] com o valor mais recente disponível
            x_novo[i] = (b[i] - s1 - s2) / A[i][i]

        # Cálculo do critério de parada
        max_diff = max(abs(x_novo[i] - x[i]) for i in range(n))
        dados_iteracao.append((k + 1, x_novo[0], x_novo[1], x_novo[2], max_diff))

        if max_diff < epsilon:
            break

        # Atualiza x com os valores mais recentes
        x = x_novo[:]

    return x_novo, dados_iteracao


# Matriz A e vetor b fornecidos
A = [[11.9, 0.0, 1.8], [0.0, 5.3, -1.8], [1.0, -1.0, -1.0]]
b = [15.0, 3.1, 0.0]

# Vetor inicial
x0 = [1.0, 1.0, 1.0]

# Resolução do sistema usando o Método de Jacobi
solucao_jacobi, dados_iteracao_jacobi = jacobi(A, b, x0)
print("\nMétodo de Jacobi")
print("\nTabela de Iterações:")
print(f"\n|{'Iteração':^10}|{'I1':^10}|{'I2':^10}|{'I3':^10}|{'Erro':^10}|")
print("-"*56)
for dado in dados_iteracao_jacobi:
    print(f"|{dado[0]:^10}|{dado[1]:^10.6f}|{dado[2]:^10.6f}|{dado[3]:^10.6f}|{dado[4]:^10.6f}|")

print("\nSolução:")
print(f"I1 = {solucao_jacobi[0]:.6f}")
print(f"I2 = {solucao_jacobi[1]:.6f}")
print(f"I3 = {solucao_jacobi[2]:.6f}")

# Resolução do sistema usando o Método de Gauss-Seidel
solucao_gauss_seidel, dados_iteracao_gauss_seidel = gauss_seidel(A, b, x0)
print("\nMétodo de Gauss-Seidel")
print("\nTabela de Iterações:")
print(f"\n|{'Iteração':^10}|{'I1':^10}|{'I2':^10}|{'I3':^10}|{'Erro':^10}|")
print("-"*56)
for dado in dados_iteracao_gauss_seidel:
    print(f"|{dado[0]:^10}|{dado[1]:^10.6f}|{dado[2]:^10.6f}|{dado[3]:^10.6f}|{dado[4]:^10.6f}|")

print("\nSolução:")
print(f"I1 = {solucao_gauss_seidel[0]:.6f}")
print(f"I2 = {solucao_gauss_seidel[1]:.6f}")
print(f"I3 = {solucao_gauss_seidel[2]:.6f}")