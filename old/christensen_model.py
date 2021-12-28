from math import pi, sqrt

def magnetic_field_strength(grav_cmb, buoyancy_flux, core_radius, layer_radius, density_core=7500, alpha_core=3 * 10 ** -5,
                            spec_heat_core=840):
    """
    From Reese and Solomotov 2010, modified from Christensen & Aubert 2006
    :param grav_cmb:
    :param buoyancy_flux:
    :param core_radius:
    :param layer_radius:
    :param density_core:
    :param alpha_core: Core thermal expansivity
    :param spec_heat_core: Specific heat of the core
    :return:
    """
    delta_r = layer_radius - core_radius
    mu = 4 * pi * 10 * -7  # permeability of free space
    return 0.9 * sqrt(mu) * (density_core ** (1 / 6)) * (((alpha_core * grav_cmb * buoyancy_flux * delta_r) / (
                4 * pi * spec_heat_core * (core_radius ** 2))) ** (1 / 3))





