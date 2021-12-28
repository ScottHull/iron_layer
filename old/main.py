from christensen_model import magnetic_field_strength
from convection import InviscidConvection
from heat_gradient import GradientFromSPH

ic = InviscidConvection()
g = GradientFromSPH(path="/Users/scotthull/Desktop/2900.csv", radius_lim=3400 * 1000, tag=3)
mean_dT_dz = g.mean_dT_dz_from_avg(num_samples=100)

density_core = 7500
alpha_core = 3 * 10 ** -5
spec_heat_core = 840
grav = 3.8
length = 150 * 1000

convective_velocity = ic.convection_velocity(g=grav, alpha=alpha_core, beta=mean_dT_dz, length=length)
flux = ic.convective_heat_flux(density=density_core, spec_heat=spec_heat_core, velocity=convective_velocity, beta=mean_dT_dz, length=length)
print(mean_dT_dz, convective_velocity, flux)

