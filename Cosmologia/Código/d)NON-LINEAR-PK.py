
'''
D) NON-LINEAR TEST
'''

import camb
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid, quad
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import sici

########################################################################################################################

# Constants
H0 = 67.5
Om_m = 0.23
ombh2 = 0.022
omch2 = Om_m * (H0 / 100) ** 2 - ombh2
Om_r = 8e-5
Om_DE = 1 - (Om_m + Om_r)
w = -1
rho_crit_0 = 2.775 * (H0 / 100) ** 2 * 1e11
rho_m0 = Om_m * rho_crit_0
delta_c = 1.686


########################################################################################################################

# Function to get the power spectrum P(k) by CAMB parameters
def P():
    pars = camb.set_params(H0=H0, ombh2=ombh2, omch2=omch2, omk=0, TCMB=2.7255)
    pars.set_matter_power(redshifts=[0], kmax=1500)
    pars.set_dark_energy(w=w)
    results = camb.get_results(pars)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-3, maxkh=100, npoints=100)
    return kh, z, pk




# Window function W(k, R)
def W(kh, R):
    return (3 / (kh * R) ** 2) * (np.sin(kh * R) / (kh * R) - np.cos(kh * R))


# Hubble function H(z)
def H(z):
    return H0 * np.sqrt(Om_m * (1 + z) ** 3 + Om_r * (1 + z) ** 4 + Om_DE * (1 + z) ** (3 * (1 + w)))


# Growth factor normalized D(z)
def G(z):
    integrand = lambda z_prime: (1 + z_prime) / H(z_prime) ** 3
    integral, _ = quad(integrand, z, np.inf)
    integral_0, _ = quad(integrand, 0, np.inf)
    return H(z) / H0 * integral / integral_0


# Function to calculate sigma(z, R)
def sigma_zR(z, R, kh, pk):
    G_z = G(z)
    integrand = (kh ** 2) * (W(kh, R) ** 2) * pk[0, :]
    integral_value = trapezoid(integrand, kh)
    return G_z * np.sqrt(integral_value / (2 * np.pi ** 2))


# Combined function: calculates sigma, its derivatives, and interpolates
def calculate_interpolation(R_values, kh, pk):
    M_vir = (4 / 3) * np.pi * R_values ** 3 * rho_m0

    # Calculate sigma for each R and z value
    sigma_z = np.array([[sigma_zR(z, R, kh, pk) for R in R_values]])

    # Interpolated sigma and derivatives
    splines_sigma = []
    splines_dln_sigma_inv_dln_M = []

    for sigma in sigma_z:
        # Interpolated sigma
        spline_sigma = InterpolatedUnivariateSpline(M_vir, sigma, k=1)
        splines_sigma.append(spline_sigma)

        # Calculate derivative d ln σ^{-1} / d ln M
        d_sigma_dM = np.gradient(sigma, M_vir)
        dln_sigma_inv_dln_M = -M_vir / sigma * d_sigma_dM

        # Interpolated d ln σ^{-1} / d ln M
        spline_dln_sigma_inv_dln_M = InterpolatedUnivariateSpline(M_vir, dln_sigma_inv_dln_M, k=1)
        splines_dln_sigma_inv_dln_M.append(spline_dln_sigma_inv_dln_M)

    return M_vir, sigma_z, splines_sigma, splines_dln_sigma_inv_dln_M


# Função Tinker g(sigma)
def g_func(sigma):
    # Tinker et al. (2008) parameters
    B, d, e, f, g = 0.494, 2.30, 0.93, 0.48, 1.403

    sigma = np.asarray(sigma)
    return B * ((sigma / e) ** (-d) + sigma**(-f)) * np.exp(-g / (sigma ** 2))

########################################################################################################################

# Função de bias b_z(M)
def b_z(nu):
    y = np.log10(400)
    A = 1 + 0.24 * y * np.exp(-(4 / y) ** 4)
    a = 0.44 * y - 0.88
    B = 0.183
    b = 1.5
    C = 0.019 + 0.107 * y + 0.19 * np.exp(-(4 / y) ** 4)
    c = 2.4

    return 1 - A * nu ** a / (nu ** a + delta_c ** a) + B * nu ** b + C * nu ** c

########################################################################################################################

# Funtion to find M_star
def find_M_star(sigma_crit, kh, pk, z, R_values):
    # Inicializar sigma(R) com base em R_values e P(k) (usando a função sigma_zR já definida)
    sigma_R = np.array([sigma_zR(z, R, kh, pk) for R in R_values])

    # Encontrar o valor de R para o qual sigma(R) atinge sigma_crit
    R_crit = R_values[np.argmin(np.abs(sigma_R - sigma_crit))]

    # Calcular M_star a partir de R_crit, usando a fórmula dependente de R (baseada na densidade virializada)
    M_star = (4 / 3) * np.pi * R_crit ** 3 * rho_m0

    print(f"The value of M corresponding to sigma(z=0, M) ≈ 1.686"
          f" is approximately {M_star:.2e} M_sun/h.")
    return M_star


def rho(rho_s, c, r_values, r_vir):
    return rho_s / (c * r_values / r_vir) * (1 / (1 + c * r_values / r_vir) ** 2)


def E2(zet):
    return Om_m * (1 + zet)**3 + Om_DE * (1 + zet)**(3 * (1 + w))


def rho_cri(zet):
    return rho_crit_0 * E2(zet)


def w_m(zet):
    return Om_m * (1 + zet) ** 3 / E2(zet)


def Delta_c(zet):
    x = w_m(zet) - 1
    return 18 * np.pi ** 2 + 82 * x - 39 * x ** 2


def r_vir(zet, M_vir):
    return ((3 * M_vir) / (4 * np.pi * Delta_c(zet) * rho_cri(zet)))**(1/3)


def c(M_vir, M_star, zet):
    return 9 / (1 + zet) * (M_vir / M_star) ** (-0.13)


def rho_s(M_vir, M_star, zet):
    c_value = c(M_vir, M_star, zet)
    return (M_vir / (4 * np.pi * r_vir(zet, M_vir) ** 3)) * (c_value ** 3 / (np.log(1 + c_value) - c_value / (1 + c_value)))


# Função para calcular u(k|Mvir)
def u_k_Mvir(k, M_vir, r_vir, rho_s, c):
    rs = r_vir / c  # Raio de escala
    prefactor = 4 * np.pi * rho_s * rs ** 3 / M_vir  # Normalização

    u_values = []  # Lista para armazenar os valores de u(k, M_vir)

    for k in kh:
        # Cálculo das integrais seno e cosseno
        Si_1, Ci_1 = sici((1 + c) * k * rs)
        Si_2, Ci_2 = sici(k * rs)

        # Termos analíticos para u(k)
        term1 = np.sin(k * rs) * (Si_1 - Si_2)
        term2 = -np.sin(c * k * rs) / ((1 + c) * k * rs)
        term3 = np.cos(k * rs) * (Ci_1 - Ci_2)

        # Calcular u(k) para o valor atual de k
        u_k = prefactor * (term1 + term2 + term3)
        u_values.append(u_k)

    return u_values

########################################################################################################################

# Non-linear Halofit power spectrum by CAMB
def P_hf():
    pars = camb.set_params(H0=H0, ombh2=ombh2, omch2=omch2, omk=0, TCMB=2.7255)
    pars.set_matter_power(redshifts=[0], kmax=1500)
    pars.set_dark_energy(w=w)
    results = camb.get_results(pars)
    pars.NonLinearModel.set_params(halofit_version='takahashi')
    kh_taka, z_taka, pk_taka = results.get_nonlinear_matter_power_spectrum(params=pars)
    return kh_taka, z_taka, pk_taka


# Integrals for Pk terms
def pk_1term(M_vir, splines_sigma, splines_dln_sigma_inv_dln_M, u_values):
    sigma_values = splines_sigma[0](M_vir)
    dln_sigma_inv_dln_M_values = splines_dln_sigma_inv_dln_M[0](M_vir)

    dn_dlnM = g_func(sigma_values) * rho_m0 / M_vir * dln_sigma_inv_dln_M_values

    integrand = (M_vir / rho_m0)**2 * dn_dlnM * (np.abs(u_values)) ** 2

    pk_1 = trapezoid(integrand, np.log(M_vir))

    return pk_1

def pk_2term(M_vir, splines_sigma, splines_dln_sigma_inv_dln_M, u_values, pk):
    sigma_values = splines_sigma[0](M_vir)
    dln_sigma_inv_dln_M_values = splines_dln_sigma_inv_dln_M[0](M_vir)

    dn_dlnM = g_func(sigma_values) * rho_m0 / M_vir * dln_sigma_inv_dln_M_values

    nu = delta_c / sigma_values
    bias = b_z(nu)

    integrand_1 = (M_vir / rho_m0) * bias * dn_dlnM
    integrand_2 = (M_vir / rho_m0) * bias * u_values * dn_dlnM

    # Normalizing for small values of k with u(k,M) = 1
    A = trapezoid(integrand_1, np.log(M_vir))
    pk_2 = trapezoid(integrand_2, np.log(M_vir))

    return (pk_2 / A) ** 2 * pk


# Call for linear and non-linear power spectrum function
kh, zs, pk = P()
kh_nlin, z_nlin, pk_nlin = P_hf()

R_values = np.logspace(-2, 8, 1000)
z = 0
r_values = np.logspace(-6, 5, 1000)
M_star = find_M_star(1.686, kh, pk, z, R_values)

pk_1term_values = np.zeros_like(kh)
pk_2term_values = np.zeros_like(kh)
p_tot_values = np.zeros_like(kh)

M_vir, sigma_z, splines_sigma, splines_dln_sigma_inv_dln_M = calculate_interpolation(R_values, kh, pk)
rho_s_value = rho_s(M_vir, M_star, z)
r_vir_value = r_vir(z, M_vir)
c_value = c(M_vir, M_star, z)

u_values = u_k_Mvir(kh, M_vir, r_vir_value, rho_s_value, c_value)

for i, u in enumerate(u_values):
    pk_1term_values[i] = pk_1term(M_vir, splines_sigma, splines_dln_sigma_inv_dln_M, u)
    pk_2term_values[i] = pk_2term(M_vir, splines_sigma, splines_dln_sigma_inv_dln_M, u, pk[0, i])
    p_tot_values[i] = pk_1term_values[i] + pk_2term_values[i]

f, axs = plt.subplots(2, 1, figsize=(8, 6),
                       gridspec_kw={'height_ratios': [6, 2], 'hspace': 0.5},
                       sharex=True)

# First plot: Halo Model, linear and non-linear Power Spectrum
axs[0].loglog(kh, pk[0, :], label="Linear P(k)", color='red')
axs[0].loglog(kh_nlin, pk_nlin[0, :], label="Halofit", color='purple')
axs[0].loglog(kh, pk_1term_values, label="One-halo term", color='royalblue', linestyle='--')
axs[0].loglog(kh, pk_2term_values, label="Two-halo term", color='forestgreen', linestyle='--')
axs[0].loglog(kh, p_tot_values, label="Halo Model P(k)", color='black')
axs[0].grid(visible=True, color='gray', linestyle='-.', linewidth=0.5)
axs[0].xaxis.set_tick_params(direction='in',  which='both', top=True, bottom=True)
axs[0].yaxis.set_tick_params(direction='in',  which='both', left=True, right=True)
axs[0].set_xlim(1e-3, 1e1)
axs[0].set_ylim(1e-1, 1e5)
axs[0].set_ylabel(r"$P(k) \, [\mathrm{Mpc}^3/h^3]$")
axs[0].set_title('Halo Model Power Spectrum')
axs[0].legend()


from scipy.interpolate import interp1d

# Interpolation values for the second plot
interpolated_pk_nlin = interp1d(kh_nlin, pk_nlin[0, :], kind='linear', bounds_error=False, fill_value="extrapolate")
pk_nlin_interpolated = interpolated_pk_nlin(kh)
difference = (p_tot_values - pk_nlin_interpolated) / pk_nlin_interpolated

# Second plot for relative difference between Halofit and Halo Model
axs[1].semilogx(kh, difference, color='black')
axs[1].set_ylim(-0.5, 0.5)
axs[1].grid(visible=True, color='gray', linestyle='-.', linewidth=0.5)
axs[1].xaxis.set_tick_params(direction='in',  which='both', top=True, bottom=True)
axs[1].yaxis.set_tick_params(direction='in',  which='both', left=True, right=True)
axs[1].set_xlabel(r"$k \, [h \, \mathrm{Mpc}^{-1}]$", fontsize=12)
axs[1].set_ylabel(r"$[P(k)^{model} - P(k)^{fit}] \  / \ P(k)^{fit}$", fontsize=12)
axs[1].set_title("Difference Between Halo Model and Halofit Spectra")

plt.tight_layout()
plt.show()