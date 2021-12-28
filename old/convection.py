from math import sqrt

class InviscidConvection:
    """
    The inviscid convection model from Prandlt Mixing Length Theory.
    Taken From Stevenson's Notes
    http://web.gps.caltech.edu/classes/ge131/notes2016/Ch14.pdf
    """

    """
    g * alpha * beta * L is an "effective gravity"
    """

    def convection_velocity(self, g, alpha, beta, length):
        """
        v = g * alpha * beta * l^2
        v = ([m / s^2] * [1 / K] * [K / m] * [m^2])^1/2 = ([1 / s^2] * [1 / m] * [m^2])^1/2 = [m / s]
        :param g:
        :param alpha:
        :param beta:
        :param length:
        :return:
        """
        return sqrt(g * alpha * beta * (length ** 2))

    def convective_heat_flux(self, density, spec_heat, velocity, beta, length):
        return density * spec_heat * velocity * beta * length

    # ------------ The two functions below are simplifications for an adiabatic flow ------------

    def velocity_adiabatic(self, length, convective_heat_flux, density, char_scale_height):
        """
        The characteristic temperature scale height for an adiabatic fluid is
        H_T = T / -(dT/dz)_ad = c_p / (alpha * g), so we can make the following simplification.
        :param length:
        :param convective_heat_flux:
        :param density:
        :param char_scale_height:
        :return:
        """
        return 0.1 * ((length * convective_heat_flux) / (density * char_scale_height)) ** (1 / 3)

    def prod_beta_L(self, alpha, velocity, g, length):
        """
        The characteristic temperature scale height for an adiabatic fluid is
        H_T = T / -(dT/dz)_ad = c_p / (alpha * g), so we can make the following simplification.
        :param alpha:
        :param velocity:
        :param g:
        :param length:
        :return:
        """
        term1 = (10 ** 2) / alpha
        term2 = (velocity / sqrt(g * length)) ** 2
        return term1 * term2
