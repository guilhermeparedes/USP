'''
Halo-Mass Function
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
def P_L():
    pars = camb.set_params(H0=H0, ombh2=ombh2, omch2=omch2, omk=0, TCMB=2.7255)
    pars.set_matter_power(redshifts=[0], kmax=100)
    pars.set_dark_energy(w=w)
    results = camb.get_results(pars)
    k_lin, z_lin, pk_lin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=100, npoints=1000)
    return k_lin, z_lin, pk_lin


# Window function W(k, R)
def W(k, R):
    return (3 / (k * R) ** 2) * (np.sin(k * R) / (k * R) - np.cos(k * R))


# Hubble function H(z)
def H(z):
    return H0 * np.sqrt(Om_m * (1 + z) ** 3 + Om_r * (1 + z) ** 4 + Om_DE * (1 + z) ** (3 * (1 + w)))


# Growth factor G(z) normalized by D(z) / D(z=0)
def G(z):
    integrand = lambda z_prime: (1 + z_prime) / H(z_prime) ** 3
    integral, _ = quad(integrand, z, np.inf)
    integral_0, _ = quad(integrand, 0, np.inf)
    return H(z) / H0 * integral / integral_0


# Function to calculate sigma(z, R)
def sigma_zR(zet, R, k, p_l):
    G_z = G(zet)
    integrand = (k ** 2) * (W(k, R) ** 2) * p_l[0, :]
    integral_value = trapezoid(integrand, k)
    return G_z * np.sqrt(integral_value / (2 * np.pi ** 2))


# Combined function: calculates sigma, its derivatives, and interpolates
def calculate_interpolation(R_value, k, p_l, zet_values):
    M_vir = (4 / 3) * np.pi * R_value ** 3 * rho_m0

    # Calculate sigma for each R and z value
    sigma_z = np.array([[sigma_zR(zet, R_value, k, p_l) for R_value in R_values] for zet in zet_values])

    # Interpolated sigma and derivatives
    splines_sigma = []
    splines_dln_sigma_inv_dln_M = []

    for sigma in sigma_z:
        # Interpolated sigma
        spline_sigma = InterpolatedUnivariateSpline(M_vir, sigma, k=3)
        splines_sigma.append(spline_sigma)

        # Calculate derivative d ln σ^{-1} / d ln M
        d_sigma_dM = np.gradient(sigma, M_vir)
        dln_sigma_inv_dln_M = -M_vir / sigma * d_sigma_dM

        # Interpolated d ln σ^{-1} / d ln M
        spline_dln_sigma_inv_dln_M = InterpolatedUnivariateSpline(M_vir, dln_sigma_inv_dln_M, k=3)
        splines_dln_sigma_inv_dln_M.append(spline_dln_sigma_inv_dln_M)


    return M_vir, sigma_z, splines_sigma, splines_dln_sigma_inv_dln_M


# Function g(sigma) by Tinker et al. 2008
def g_func(sigma):
    B, d, e, f, g = 0.494, 2.30, 0.93, 0.48, 1.403

    sigma = np.asarray(sigma)
    return B * ((sigma / e) ** (-d) + sigma**(-f)) * np.exp(-g / (sigma ** 2))



# Call for power spectrum function
k_value, z_value, pk_value = P_L()

R_values = np.logspace(-2, 8, 1000)
z_values = [0, 1]

# Compute dn/dlnM for each z value
dn_dlnM_curves = []

# Loop principal
for z in z_values:
    M_values, sigma_z, splines_sigma, splines_dln_sigma_inv_dln_M = \
        calculate_interpolation(R_values, k_value, pk_value, [z])
    sigma_at_z = np.asarray(sigma_z[0])

    dn_dlnM = g_func(sigma_at_z) * rho_m0 / M_values * \
                splines_dln_sigma_inv_dln_M[0](M_values)

    dn_dlnM_curves.append((M_values, dn_dlnM))

# Plot dn/dlnM versus M for z = 0 and z = 1
f = plt.figure(figsize=(6, 6))
ax = f.add_subplot(111)

for z, (M_values, dn_dlnM_g) in zip(z_values, dn_dlnM_curves):
    if z == 0:
        plt.loglog(M_values, dn_dlnM_g, label=f'z = {z}', color='blue', linestyle='-')
    else:
        plt.loglog(M_values, dn_dlnM_g, label=f'z = {z}', color='red', linestyle='--')



plt.grid(visible=True, color='gray', linestyle='-.', linewidth=0.5)
ax.xaxis.set_tick_params(direction='in',  which='both', top=True, bottom=True)
ax.yaxis.set_tick_params(direction='in',  which='both', left=True, right=True)
plt.xlim(1e12, 1e16)
plt.ylim(1e-10, 1e-2)
plt.xlabel(r'$M \ [h^{-1}M_\odot]$', fontsize=12)
plt.ylabel(r'$dn/d\ln M \  [h^3 \mathrm{Mpc}^{-3}]$', fontsize=12)
plt.title('Halo Mass Function')
plt.legend()
plt.show()
