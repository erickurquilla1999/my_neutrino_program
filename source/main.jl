include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/grid_functions.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/particles_init.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/input_parameters.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/evolve.jl")

# compute cells centers in base of imput parameters
cells__center=cells_center() 	
# cells__center ---> 1: x center cell position

#create the particles with the initial condition
initial_particles_data=particles_initial_condition(cells__center)
# initial_particles_data ---> 1: x particle position
# initial_particles_data ---> 2: x particle direction
# initial_particles_data ---> 3: rho neutrinos
# initial_particles_data ---> 4: rho bar antineutrinos 
# initial_particles_data ---> 5: number of neutrinos
# initial_particles_data ---> 6: number of antineutrinos

# evolve de system of neutrinos
evolve_particles(simulation_initial_time,initial_particles_data)

