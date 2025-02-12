'''
Método de secante
'''

import numpy as np

def secante(f, x0, x1, e=10 ** -5):

    iter_count = 0
    print("-" * 60)
    print(f"|{'n':^10}|{'xn':^15}|{'f(xn)':^15}|{'Precisao':^15}|")
    print("-" * 60)

    while abs(f(x1)) > e:  # Enquanto o valor absoluto de f(x1) for maior que a tolerância

        f_x0 = f(x0)
        f_x1 = f(x1)

        if f_x1 == f_x0:
            raise ZeroDivisionError("A diferença f(x_n) - f(x_{n-1}) é zero. O método da secante não pode continuar.")

        # Impressão dos valores atuais
        print(f"|{iter_count + 1:^10}|{x1:^15.8f}|{f_x1:^15.8f}|{abs(f_x1):^15.8f}|")

        # Atualização do ponto x2
        x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)

        # Atualizando os pontos
        x0, x1 = x1, x2

        iter_count += 1

    # Impressão final dos valores quando a precisão desejada é alcançada
    print(f"|{iter_count + 1:^10}|{x1:^15.8f}|{f(x1):^15.8f}|{abs(f(x1)):^15.8f}|")
    print("-" * 60)

    return x1


# Definindo a função x^(3/4) - cos(x^2) e sua derivada
def minha_funcao(x):
    return -(14.4 / x ** 2) + (1.38e3 / 0.328) * np.exp(-(x / 0.328))


# Valores iniciais
x0 = 0.5
x1 = 1.5

# Calculando a raiz
raiz = secante(minha_funcao, x0, x1)
print(f"\nA raiz aproximada é: {raiz:.8f}")