'''
Halo-Profile
'''

import camb
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid, quad
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import sici

# Constants
H0 = 67.5
ombh2 = 0.022
omch2 = 0.122
Om_m = (ombh2 + omch2) / (H0 / 100) ** 2
Om_r = 8e-5
Om_DE = 1 - (Om_m + Om_r)
w = -1
rho_crit_0 = 2.775e11 * (H0 / 100) ** 2
rho_m0 = Om_m * rho_crit_0


# Function to get the power spectrum P(k) by CAMB parameters
def P():
    pars = camb.set_params(H0=H0, ombh2=ombh2, omch2=omch2, omk=0, TCMB=2.7255)
    pars.set_matter_power(redshifts=[0], kmax=10000)
    pars.set_dark_energy(w=w)
    results = camb.get_results(pars)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-3, maxkh=10000, npoints=1000)
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
def calculate_interpolation(R_values, kh, pk, z_values):
    M_vir = (4 / 3) * np.pi * R_values ** 3 * rho_m0

    # Calculate sigma for each R and z value
    sigma_z = np.array([[sigma_zR(z, R, kh, pk) for R in R_values] for z in z_values])

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




'''
Halo Profile
'''

def find_M_star(sigma_crit, kh, pk, z, R_values):
    sigma_R = np.array([sigma_zR(z, R, kh, pk) for R in R_values])
    R_crit = R_values[np.argmin(np.abs(sigma_R - sigma_crit))]
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
    return ((3 * M_vir) / (4 * np.pi * Delta_c(zet) * rho_cri(zet))) ** (1/3)


def c(M_vir, M_star, zet):
    return 9 / (1 + zet) * (M_vir / M_star) ** (-0.13)


def rho_s(M_vir, M_star, zet):
    c_value = c(M_vir, M_star, zet)
    return (M_vir / (4 * np.pi * r_vir(zet, M_vir) ** 3)) * (c_value **3 / (np.log(1 + c_value) - c_value / (1 + c_value)))


# Call for power spectrum function
kh, zs, pk = P()

R_values = np.logspace(-2, 8, 1000)
z_values = [0, 1]
r_values = np.logspace(-6, 2, 1000)
M_vir_values1 = [1e14, 1e15]

f = plt.figure(figsize=(6, 6))
ax = f.add_subplot(111)
colors_1 = ['royalblue', 'red', 'forestgreen', 'cyan']  # Lista de cores

for i, z in enumerate(z_values):
    for j, M_vir in enumerate(M_vir_values1):
        M_star = find_M_star(1.686, kh, pk, z, R_values)
        rho_s_value = rho_s(M_vir, M_star, z)
        r_vir_value = r_vir(z, M_vir)
        c_value = c(M_vir, M_star, z)
        rho_values = rho(rho_s_value, c_value, r_values, r_vir_value)
        color = colors_1[i * len(M_vir_values1) + j]
        plt.loglog(r_values, rho_values, label=f'M = {M_vir:.0e}, z = {z}', color=color)


plt.grid(visible=True, color='gray', linestyle='-.', linewidth=0.5)
ax.xaxis.set_tick_params(direction='in',  which='both', top=True, bottom=True)
ax.yaxis.set_tick_params(direction='in',  which='both', left=True, right=True)

plt.xlim(1e-2, 1e2)
plt.ylim(1e6, 1e17)

plt.xlabel(r'$r$', fontsize=12)
plt.ylabel(r'$\rho(r|M_{vir},z)$', fontsize=12)
plt.legend(loc='best')
plt.title('Halo Profile - Real Space')
plt.show()


M_vir_values2 = [1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16]


# Function to calculate u(k|M_vir) numerically
def u_numeric(k, M_vir, r_vir, rho_s, c):
    rs = r_vir / c
    prefactor = 4 * np.pi * rho_s * rs ** 3 / M_vir

    def integrand(r):
        rho_r = rho(rho_s, c, r, r_vir)
        return np.sin(k * r) / (k * r) * rho_r * r ** 2

    integral_value, _ = quad(integrand, 0, r_vir)

    return prefactor * integral_value

k_values = kh
z = 0
colors_2 = ['darkred', 'royalblue', 'forestgreen', 'gold', 'indigo', 'black', 'magenta', 'teal', 'orchid', 'darkviolet']
f = plt.figure(figsize=(6, 6))
ax = f.add_subplot(111)

for i, M_vir in enumerate(M_vir_values2):
    M_star = find_M_star(1.686, kh, pk, z, R_values)
    rho_s_value = rho_s(M_vir, M_star, z)
    r_vir_value = r_vir(z, M_vir)
    c_value = c(M_vir, M_star, z)

    u_values_num = []

    for k in k_values:
        u_values_num.append(u_numeric(k, M_vir, r_vir_value, rho_s_value, c_value))

    plt.loglog(k_values, u_values_num, label=f'M = {M_vir:.0e}', color=colors_2[i])


plt.grid(visible=True, color='gray', linestyle='-.', linewidth=0.5)
ax.xaxis.set_tick_params(direction='in',  which='both', top=True, bottom=True)
ax.yaxis.set_tick_params(direction='in',  which='both', left=True, right=True)
plt.xlim(1e-2, 1e4)
plt.ylim(1e-3, 1e1)
plt.xlabel('$k$ [h/Mpc]', fontsize=12)
plt.ylabel(r'$u[k|M_{vir}]$', fontsize=12)
plt.title('Halo Profile - Fourier Space ')
plt.legend(loc='best')
plt.show()

