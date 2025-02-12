'''
Método de bisseção
'''

import numpy as np

def bissecao(f, a, b, e=10 ** -8):
    # Verificar se a funcao muda de sinal no intervalo [a, b]
    if f(a) * f(b) >= 0:
        raise ValueError("A funcao deve ter sinais opostos em a e b.")

    iter_count = 0  # conta o numero de iteracoes
    print("-" * 108)
    print(f"|{'Iteracao':^10}|{'x1':^15}|{'x2':^15}|{'xm':^15}|{'f(x1)':^15}|{'f(xm)':^15}|{'Precisao':^15}|")
    print("-" * 108)
    while (b - a) / 2.0 > e:
        m = (a + b) / 2.0  # Ponto medio
        f_a = f(a)
        f_m = f(m)

        # Impressao dos valores atuais
        print(f"|{iter_count + 1:^10}|{a:^15.8f}|{b:^15.8f}|{m:^15.8f}|{f_a:^15.8f}|{f_m:^15.8f}|{abs(b - a):^15.8f}|")

        if f_m == 0.0:
            return m  # A raiz exata foi encontrada
        elif f_a * f_m < 0:
            b = m  # A raiz esta no intervalo [a, c]
        else:
            a = m  # A raiz esta no intervalo [c, b]
        iter_count += 1

    # Impressao final dos valores quando a precisao desejada e alcancada
    print(
        f"|{iter_count + 1:^10}|{a:^15.8f}|{b:^15.8f}|{(a + b) / 2.0:^15.8f}|{f(a):^15.8f}|{f((a + b) / 2.0):^15.8f}|{abs(b - a):^15.8f}|")
    print("-" * 108)

    return (a + b) / 2.0  # Retorna a raiz aproximada


# Definindo a funcao x^(3/4) - cos(x^2)
def minha_funcao(x):
    return x ** (3 / 4) - np.cos(x ** 2)


# Definindo intervalo
a = 0.5
b = 1.0

# Calculando a raiz
raiz = bissecao(minha_funcao, a, b)
print(f"\nA raiz aproximada e: {raiz:.8f}")