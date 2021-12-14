import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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

def V2(V1, assumed_alpha, T1, T2):
    deltaV = assumed_alpha * V1 * (T2 - T1)  # integrated thermal expansion equation for volume
    return V1 + deltaV  # by definition, V1 = V1 + deltaV

print(pressure())

# df = pd.read_csv("anderson_iron_300k.xlsx")
# temps = [300] + list(np.arange(1000, 7000 + 1000, 1000))
# for row in df.index:
#     for index, T1 in enumerate(temps):
#         if index + 1 < len(temps):
#             isobar_pressure = df["P"][row]
#             V_1 = df["V"][row]
#             T2 = temps[index + 1]
#

