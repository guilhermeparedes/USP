'''
Halo-Bias
'''

import camb
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid, quad
from scipy.interpolate import InterpolatedUnivariateSpline

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

# Function to get the power spectrum P(k) by CAMB parameters
def P():
    pars = camb.set_params(H0=H0, ombh2=ombh2, omch2=omch2, omk=0, TCMB=2.7255)
    pars.set_matter_power(redshifts=[0],  kmax=100)
    pars.set_dark_energy(w=w)
    results = camb.get_results(pars)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=10, npoints=1000)
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
    B = 0.494
    d = 2.30
    e = 0.93
    f = 0.48
    g = 1.403

    sigma = delta_c / np.asarray(sigma)
    return B * ((sigma / e) ** (-d) + sigma**(-f)) * np.exp(-g / (sigma ** 2))



# Call for power spectrum function
kh, zs, pk = P()

R_values = np.logspace(-2, 3, 1000)
z_values = [0, 1]

'''
Halo Bias
'''

# The Bias Function
def b_z(nu):
    y = np.log10(400)
    A = 1 + 0.24 * y * np.exp(-(4 / y) ** 4)
    a = 0.44 * y - 0.88
    B = 0.183
    b = 1.5
    C = 0.019 + 0.107 * y + 0.19 * np.exp(-(4 / y) ** 4)
    c = 2.4

    return 1 - A * nu ** a / (nu ** a + delta_c ** a) + B * nu ** b + C * nu ** c



f = plt.figure(figsize=(6, 6))
ax = f.add_subplot(111)

# Loop principal
for z in z_values:
    M_vir, sigma_z, splines_sigma, splines_dln_sigma_inv_dln_M = calculate_interpolation(R_values, kh, pk, [z])
    sigma_at_z = sigma_z[0]
    nu = delta_c / sigma_at_z
    bias_values = b_z(nu)
    if z == 0:
        plt.plot(M_vir, bias_values, label=f'z = {z}', color='blue', linestyle='-')
    else:
        plt.plot(M_vir, bias_values, label=f'z = {z}', color='red', linestyle='--')


plt.xscale('log')
plt.yscale('linear')
plt.grid(visible=True, color='gray', linestyle='-.', linewidth=0.5)
ax.xaxis.set_tick_params(direction='in',  which='both', top=True, bottom=True)
ax.yaxis.set_tick_params(direction='in',  which='both', left=True, right=True)
plt.xlim(1e12, 1e16)
plt.ylim(0, 30)
plt.xlabel(r'$M \ [h^{-1}M_\odot]$', fontsize=12)
plt.ylabel(r'$b(z, M)$ ', fontsize=12)
plt.title(r'Halo Bias')
plt.legend()
plt.show()
