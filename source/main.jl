include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/grid_functions.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/particles_init.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/input_parameters.jl")



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

@time begin

particles__number_cell=particles_number_cell(initial_particles_data[1],initial_particles_data[2],initial_particles_data[3])
# particles__number_cell ---> array of the number of cell the particles bellows

shape_func=shape_function(initial_particles_data[1],initial_particles_data[2],initial_particles_data[3],cells__center[1],cells__center[2],cells__center[3],particles__number_cell)
# shape_func ---> 1: w previous cell
# shape_func ---> 2: w actual cell
# shape_func ---> 3: w nest cell

grid_info_deposited=deposition_from_particles_to_grid(initial_particles_data[4],initial_particles_data[5],initial_particles_data[6],initial_particles_data[7],initial_particles_data[8],initial_particles_data[9],initial_particles_data[10],shape_func[1],shape_func[2],shape_func[3],particles__number_cell)
# grid_info_deposited ---> 1: neutrinos number density 
# grid_info_deposited ---> 2: antineutrinos number density
# grid_info_deposited ---> 3: x neutrinos number density flux
# grid_info_deposited ---> 4: y neutrinos number density flux
# grid_info_deposited ---> 5: z neutrinos number density flux
# grid_info_deposited ---> 6: x antineutrinos number density flux
# grid_info_deposited ---> 7: y antineutrinos number density flux
# grid_info_deposited ---> 8: z antineutrinos number density flux
# grid_info_deposited ---> 9: lepton number density
# grid_info_deposited ---> 10: x lepton number density flux
# grid_info_deposited ---> 11: y lepton number density flux
# grid_info_deposited ---> 12: z lepton number density flux

particle_interpolation_data=interpolation_from_grid_to_particles(shape_func[1],shape_func[2],shape_func[3],particles__number_cell,grid_info_deposited) 
# particle_interpolation_data ---> 1: neutrinos number density 
# particle_interpolation_data ---> 2: antineutrinos number density
# particle_interpolation_data ---> 3: x neutrinos number density flux
# particle_interpolation_data ---> 4: y neutrinos number density flux
# particle_interpolation_data ---> 5: z neutrinos number density flux
# particle_interpolation_data ---> 6: x antineutrinos number density flux
# particle_interpolation_data ---> 7: y antineutrinos number density flux
# particle_interpolation_data ---> 8: z antineutrinos number density flux
# particle_interpolation_data ---> 9: lepton number density
# particle_interpolation_data ---> 10: x lepton number density flux
# particle_interpolation_data ---> 11: y lepton number density flux
# particle_interpolation_data ---> 12: z lepton number density flux

end	




