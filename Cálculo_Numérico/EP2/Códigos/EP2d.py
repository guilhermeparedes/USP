'''
Matriz Jocobi
'''


def compute_matriz_jacobi(A):
    """
    Calcula a matriz de Jacobi J para a matriz A fornecida.
    :param A: Matriz dos coeficientes (lista de listas)
    :return: Matriz de Jacobi J (lista de listas)
    """
    n = len(A)  # Número de linhas/colunas da matriz A

    # Inicializar as matrizes D, L e U
    D = [[0.0] * n for _ in range(n)]
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    # Preencher a matriz D com os elementos da diagonal de A
    for i in range(n):
        D[i][i] = A[i][i]

    # Preencher as matrizes L e U com os elementos de A
    for i in range(n):
        for j in range(n):
            if i > j:
                L[i][j] = A[i][j]
            elif i < j:
                U[i][j] = A[i][j]

    # Calcular a matriz de Jacobi J = (L + U) / D
    J = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(n):
            if D[i][i] != 0:
                if i != j:
                    J[i][j] = (-1)*(L[i][j] + U[i][j]) / D[i][i]
                else:
                    J[i][j] = 0.0  # A diagonal de J deve ser zero

    return J


# Matriz A fornecida
A = [[11.9, 0.0, 1.8],[0.0, 5.3, -1.8],[1.0, -1.0, -1.0]]

# Calcular a matriz de Jacobi J
J = compute_matriz_jacobi(A)

# Imprimir a matriz J
print("Matriz de Jacobi J:")
for row in J:
    print("".join(f"{elem:^10.4f}" for elem in row))
print()

'''
os autovalores sao {0. + 0.700631 I, 0. - 0.700631 I, 0. + 0. I}.
o raio de convergencia é \rho =0.700631.
fazendo que \rho^k ~ 10^{-precisão}. nesse caso precisão = 3. fazendo as contas obtem-se p ~3.09 
'''