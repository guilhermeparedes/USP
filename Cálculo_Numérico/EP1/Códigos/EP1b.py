'''
Método de Newton
'''
'''
import numpy as np

def newton_raphson(f, df, x0, e):
    iter_count = 0
    print("-" * 76)
    print(f"|{'n':^10}|{'xn':^15}|{'f(xn)':^15}|{'df(xn)':^15}|{'Precisao':^15}|")
    print("-" * 76)

    x = x0
    while abs(f(x)) > e:  # Sendo abs o valor absoluto
        f_x = f(x)
        df_x = df(x)

        if df_x == 0:
            raise ZeroDivisionError("A derivada e zero. O metodo de Newton-Raphson nao pode continuar.")

        # Impressao dos valores atuais
        print(f"|{iter_count + 1:^10}|{x:^15.8f}|{f_x:^15.8f}|{df_x:^15.8f}|{abs(f_x):^15.8f}|")

        x = x - f_x / df_x
        iter_count += 1

    # Impressao final dos valores quando a precisao desejada e alcancada
    print(
        f"|{iter_count + 1:^10}|{x:^15.8f}|{f(x):^15.8f}|{df(x):^15.8f}|{abs(f(x)):^15.8f}|")
    print("-" * 76)

    return x


# Definindo a funcao x^(3/4) - cos(x^2) e sua derivada
def minha_funcao(x):
    return x ** (3 / 4) - np.cos(x ** 2)


# Definindo derivada da funcao
def derivada_funcao(x):
    return (3 / 4) * x ** (-1 / 4) + 2 * x * np.sin(x ** 2)


# Valor inicial
x0 = 1
e = 10 ** -8
# Calculando a raiz
raiz = newton_raphson(minha_funcao, derivada_funcao, x0, e)
print(f"\nA raiz aproximada e: {raiz:.8f}")
'''

import numpy as np

def newton_raphson(f, df, x0, e):
    iter_count = 0
    print("-" * 76)
    print(f"|{'n':^10}|{'xn':^15}|{'f(xn)':^15}|{'df(xn)':^15}|{'Precisao':^15}|")
    print("-" * 76)

    x = x0
    while True:
        f_x = f(x)
        df_x = df(x)
        x_ant = x  # Guarda o valor atual de x

        if df_x == 0:
            raise ZeroDivisionError("A derivada é zero. O método de Newton-Raphson não pode continuar.")

        # Atualiza x
        x = x_ant - f_x / df_x

        # Cálculo da precisão
        precisao = abs((x - x_ant) / x)


        # Verifica se a precisão está dentro do limite
        if precisao < e:
            break

        iter_count += 1

        print(f"|{iter_count:^10}|{x:^15.8f}|{f_x:^15.8f}|{df_x:^15.8f}|{precisao:^15.8f}|")




    return x


# Definindo a função sin(x) - e^(-x)
def minha_funcao(x):
    return np.sin(x) - np.exp(-x)

# Definindo a derivada da função
def derivada_funcao(x):
    return np.cos(x) + np.exp(-x)

# Valor inicial
x0 = 1
e = 10 ** -5
# Calculando a raiz
raiz = newton_raphson(minha_funcao, derivada_funcao, x0, e)
print(f"\nA raiz aproximada é: {raiz:.8f}")