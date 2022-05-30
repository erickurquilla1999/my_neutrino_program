include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/grid_functions.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/particles_init.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/input_parameters.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/evolve.jl")

# compute cells centers in base of imput parameters
cells__center=cells_center() 	
# cells__center ---> 1: x center cell position
# cells__center ---> 2: y center cell position
# cells__center ---> 3: z center cell position

#create the particles with the initial condition
initial_particles_data=particles_initial_condition(cells__center[1],cells__center[2],cells__center[3])
# initial_particles_data ---> 1: x particle position
# initial_particles_data ---> 2: y particle position
# initial_particles_data ---> 3: z particle position
# initial_particles_data ---> 4: x particle direction
# initial_particles_data ---> 5: y particle direction
# initial_particles_data ---> 6: z particle direction
# initial_particles_data ---> 7: rho neutrinos
# initial_particles_data ---> 8: rho bar antineutrinos 
# initial_particles_data ---> 9: number of neutrinos
# initial_particles_data ---> 10: number of antineutrinos

evolve_particles(simulation_initial_time,initial_particles_data)


