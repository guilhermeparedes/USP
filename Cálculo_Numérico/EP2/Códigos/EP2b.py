'''
Método de eliminação de Gauss
'''

def eliminacao_gauss(A, b):
    n = len(b)

    # Criar a matriz aumentada
    Ab = [row + [b[i]] for i, row in enumerate(A)]
    print("Matriz aumentada inicial:")

    for row in Ab:
        print(" ".join(f"{elem:8.1f}" for elem in row))
    print()

    # Etapa 1: Eliminação para forma escalonada com pivotamento parcial
    for i in range(n):
        # Encontrar o pivô e trocar linhas se necessário
        max_row = i
        for k in range(i + 1, n):
            if abs(Ab[k][i]) > abs(Ab[max_row][i]):
                max_row = k
        if i != max_row:
            Ab[i], Ab[max_row] = Ab[max_row], Ab[i]
            print(f"\nTrocando linha {i + 1} com linha {max_row + 1}:")
            for row in Ab:
                print(" ".join(f"{elem:8.1f}" for elem in row))
            print()

        # Zerar elementos abaixo do pivô
        for j in range(i + 1, n):
            if Ab[j][i] != 0:  # Verifique se o elemento abaixo do pivô não é zero
                factor = Ab[j][i] / Ab[i][i]
                print(f"\nEliminando elemento na linha {j + 1}, coluna {i + 1} com fator {factor:.7f}:")
                for k in range(i, n + 1):
                    Ab[j][k] -= factor * Ab[i][k]
                for row in Ab:
                    print(" ".join(f"{elem:8.1f}" for elem in row))
                print()

    # Etapa 2: Substituição reversa para resolver as variáveis
    x = [0] * n
    for i in reversed(range(n)):
        x[i] = Ab[i][n] / Ab[i][i]
        print(f"Calculando x{i + 1} = {Ab[i][n]:.7f} / {Ab[i][i]:.7f} = {x[i]:.7f}")
        for j in range(i + 1, n):
            x[i] -= Ab[i][j] * x[j] / Ab[i][i]
            print(f"Ajustando x{i + 1} com x{j + 1} = {x[j]:.7f}")
        print()

    return x


A = [
    [0.0, 5.3, -1.8],
    [11.9, 0.0, 1.8],
    [1, -1, -1]
]

b = [3.1, 15.0, 0]

# Resolve o sistema
solucao = eliminacao_gauss(A, b)

print("Solução:")
print(f"I1 = {solucao[0]:.7f}")
print(f"I2 = {solucao[1]:.7f}")
print(f"I3 = {solucao[2]:.7f}")
