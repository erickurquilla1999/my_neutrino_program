include("constants.jl")

# grid parameters
number_of_cells=1
cell_x_lenght=1e10 # m

particles_number_density=0.0 # m^(-3)
perturbation_amplitud=0.0

background_electron_number_density=1e32 # m^(-3)
background_antielectron_number_density=0.0 # m^(-3)

#mixing and CP phase angles numerical values in degrees             
theta_12 = 39.0  * pi/180
theta_13 = 29.0  * pi/180
theta_23 = 27.0  * pi/180
delta_cp = 220.0 * pi/180
   
#mass numerical values in eV/c^2               
mass_1=0.01 * eV_to_J / c2 # kg 
mass_2=0.03 * eV_to_J / c2 # kg 
mass_3=0.02 * eV_to_J / c2 # kg

#energy numerical value in eV              
neutrino_energy=50.0e+6 * eV_to_J # J

#number of neutrino flavor to be simulated
number_of_flavors=2

#simulation parameters
simulation_initial_time=0.0 # s
time_step=0.0005169584630893724/100 # s
number_of_steps=200





