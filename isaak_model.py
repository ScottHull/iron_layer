import matplotlib.pyplot as plt
import numpy as np

plt.style.use("dark_background")

"""
Taken from Isaak & Anderson 2002
https://www.sciencedirect.com/science/article/pii/S0921452602018586?casa_token=TQWBXXUQvmIAAAAA:l3QXlcBb7MNegpJPPgdz6urpRgOB6pVUXhOKEPqA7RosbRG1P7xxHQZpx0HH2hiePZXPp67p_g
"""

def vibrational_gruneisen(V, a=97.3, b=2996, c=0.33):
    """
    Taken from Anderson et al. 2001.
    https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2000GL008544
    The vibrational Gruneisen parameter from a fitted Debye model.
    gamma = (d theta / dV)_T
    where theta = a + b e^(-cV)
    where theta is the Debye temperature.
    :param V:
    :return:
    """
    return V * ((b * c) / ((a * np.exp(c * V)) + b))

def debye_temperature(V, a=97.3, b=2996, c=0.33):
    """
    Taken from Anderson et al. 2001.
    https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2000GL008544
    theta = a + b e^(-cV)
    where theta is the Debye temperature.
    :param V:
    :return:
    """
    return a + (b * np.exp(-c * V))

def birch_murnaghan_3rd_order(K0, V0, V, K0_prime):
    """
    Third-order birtch-Murnaghan isothermal EoS, P(V).
    Anderson et al. 2001, with fitted values of:
    K0 = 155.8 GPa
    K0_prime = 5.81
    V0 = 6.73 cm^3 mol^-1
    :param K0:
    :param V0:
    :param V:
    :param K0_prime:
    :return: P(V)
    """
    term1 = 3 * K0 / 2
    term2 = (V0 / V) ** (7 / 3) - (V0 / V) ** (5 / 3)
    term3 = 1 + (3 / 4) * (K0_prime - 4) * ((V0 / V) ** (2 / 3) - 1)
    return term1 * term2 * term3

def birch_murnaghan_3rd_order_bulk_modulus(K0, K0_prime, P):
    """
    Assume bulk modulus varies with pressure: K = K0 + K'_0 P
    https://mcbrennan.github.io/BMderivation.pdf
    :param K0:
    :param K0_prime:
    :param P: pressure
    :return:
    """
    return K0 + (K0_prime * P)

volumes = np.arange(4, 7 + 0.1, 0.1)
fig, axs = plt.subplots(1, 3, figsize=(16, 9), sharex='all', gridspec_kw={"wspace": 0.20})
ax1, ax2, ax3 = axs.flatten()
ax1.plot(
    volumes,
    [vibrational_gruneisen(V) for V in volumes],
    linewidth=2.0
)
ax2.plot(
    volumes,
    [debye_temperature(V) for V in volumes],
    linewidth=2.0
)
ax3.plot(
    volumes,
    [birch_murnaghan_3rd_order(V=V, V0=6.73, K0=155.8, K0_prime=5.81) for V in volumes]
)
ax1.set_title("Vibrational Gruneisen Parameter")
ax1.set_ylabel("Gamma")
ax2.set_title("Fitted Debye Temperature")
ax2.set_ylabel("Theta")
ax3.set_title("3rd Order Birch-Murnaghan EoS")
ax3.set_ylabel("Pressure (GPa)")
for ax in axs.flatten():
    ax.set_xlabel("Volume (cm^3 mol^-1)")
    ax.grid(alpha=0.4)
plt.show()



def P_th(deltaT):
    a = 12.1 * 10 ** -3  # GPa / K
    b = 7.8 * 10 ** -7  # GPa / K^2
    return (a * deltaT) + (0.5 * b * (deltaT ** 2))

def P_0(V, T=300):
    return 0

def T_P_EoS(deltaT):
    """
    P (V,T) = P_0 (V, 300 K) + P_TH (T)
    :param deltaT:
    :return:
    """
    return P_0(0) + P_th(deltaT)

# volume dependence of parameters of HPC iron at 300K
test_data_table1 = {
    "V": 5.90,  # cm/mol
    "P": 30,  # GPa
    "K_T": 319,  # GPa,
    "alpha": 3.88,  # 10^-5 K^-1
    "alpha * K_T": 12.4,  # 10^-3 GPa/K
    "debye_T": 524,  # K
    "gruniesen": 1.72,
}
