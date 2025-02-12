'''
3 a)
'''
# Parâmetros do gerador linear congruencial
a = 1103515245
c = 12345
m = 2147483647
seed = 12693754  # Semente inicial

# Função para gerar próximo número aleatório
def next_random(seed):
    new_seed = (a * seed + c) % m
    return new_seed / m, new_seed  # Retorna o valor normalizado e a nova semente


'''
3 b)
'''

# Função da curva
def f(x):
    return 1 - x**2

# Função para estimar a área com 100 pontos usando o método Monte Carlo com o LCG
def estimate_area_lcg(n_points, seed):
    count_inside = 0
    for _ in range(n_points):
        x_random, seed = next_random(seed)
        y_random, seed = next_random(seed)
        if y_random < f(x_random):
            count_inside += 1
    return count_inside / n_points, seed

# Parâmetros
n_points = 100

# Estimativa inicial da área
area_estimate, seed = estimate_area_lcg(n_points, seed)
print(f"Estimativa inicial da área sob a curva com 100 pontos: {area_estimate}")


'''
3 c)
'''

N_t_values = [2 ** i for i in range(1, 18)]
results = []

print("-" * 50)
print(f"|{'Nt':^10}|{'Im':^15}|{'σ':^10}|{'σm':^10}|")
print("-" * 50)

# Cálculos para diferentes N_t
for N_t in N_t_values:
    estimates = []
    current_seed = seed
    for _ in range(N_t):
        area, current_seed = estimate_area_lcg(n_points, current_seed)
        estimates.append(area)

    # Cálculo de Im (média das integrais estimadas)
    I_m = sum(estimates) / N_t

    # Cálculo de σ^2 (variância) manualmente
    sigma = sum((I - I_m) ** 2 for I in estimates) / (N_t - 1)

    # Cálculo de σm (desvio padrão da média)
    sigma_m = (sigma / N_t) ** 0.5

    # Armazenar os resultados
    results.append((N_t, I_m, sigma, sigma_m))

    print(f"|{N_t:^10}|{I_m:^15.6f}|{sigma:^10.6f}|{sigma_m:^10.6f}|")


