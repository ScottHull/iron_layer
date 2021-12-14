import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("dark_background")

def birch_murnaghan_3rd_order(V, V0=6.73, K0=155.8, K0_prime=5.81):
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

def P_th(deltaT):
    a = 12.1 * 10 ** -3  # GPa / K
    b = 7.8 * 10 ** -7  # GPa / K^2
    return (a * deltaT) + (0.5 * b * (deltaT ** 2))

def pressure(V, deltaT):
    """
    From Anderson et al. 2003
    https://reader.elsevier.com/reader/sd/pii/S0921452602018586?token=79DF393A37E61505739170DECD664BC38CDC50155A39DCBD0C05FE3B0A43813A5204D17AEABEAB9914D8706505BB487B&originRegion=us-east-1&originCreation=20211214012248
    The P-V-T EoS used for finding thermal expansion.
    P(V, T) = P_0 (V, 3000) + P_th (T)
    :return:
    """
    return birch_murnaghan_3rd_order(V) + P_th(deltaT)

def expanded_volume(V1, assumed_alpha, T1, T2):
    deltaV = assumed_alpha * V1 * (T2 - T1)  # integrated thermal expansion equation for volume
    return V1 + deltaV  # by definition, V1 = V1 + deltaV


def get_v_from_isotherm(T, target_P):
    """
    Returns the volume of a target P from an isotherm.
    :param T:
    :param target_P:
    :return:
    """
    volumes = list(np.arange(3.9, 7 + 0.1, 0.1))
    isobar = [pressure(V=V, deltaT=(T - 300)) for V in volumes]
    min_val = 1e10
    target_V = None
    for index, i in enumerate(isobar):
        test = abs(i - target_P)
        if test < min_val:
            min_val = test
            target_V = volumes[index]
    return target_V

alphas = np.arange(1e-6, 5e-5, 1e-6)
temperatures = np.arange(1000, 8000, 1000)
pres = []
v2s = []
target_P = 181.7
v1 = 0.68
T1 = 300
found_alphas = []
for T2 in temperatures:
    target_v = get_v_from_isotherm(T=T2, target_P=target_P)
    candidate_alpha = 1e99
    candidate_v2 = 1e99
    min_diff = 1e99
    for alpha in alphas:
        v2 = expanded_volume(V1=v1, assumed_alpha=alpha, T1=T1, T2=T2)
        test = abs(v2 - target_v)
        if test < min_diff:
            candidate_alpha = alpha
            candidate_v2 = v2
            min_diff = test
    found_alphas.append(candidate_alpha * 10 ** 5)
    v2s.append(candidate_v2)


fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    temperatures,
    found_alphas,
    linewidth=2.0
)

plt.show()
