include("grid_functions.jl")
include("particles_init.jl")
include("input_parameters.jl")
include("evolve.jl")
include("evolve.jl")

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

