include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/input_parameters.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/constants.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/numerical_methods.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/deposition_and_interpolation.jl")

function compute_hamiltonians(x_dir,interpotation_data)

	println("Computing hamiltonians")							

	# particle_interpolation_data ---> 1: neutrinos number density 
	# particle_interpolation_data ---> 2: antineutrinos number density
	# particle_interpolation_data ---> 3: x neutrinos number density flux
	# particle_interpolation_data ---> 4: x antineutrinos number density flux
	
	N          = interpotation_data[1]
	N_bar      = interpotation_data[2]
	x_flux     = interpotation_data[3]
	x_flux_bar = interpotation_data[4]	

	H_vacuum=fill(zeros(Complex{Float64},(3,3)),length(x_dir))
	H_matter=fill(zeros(Complex{Float64},(3,3)),length(x_dir))
	H_neutrino=fill(zeros(Complex{Float64},(3,3)),length(x_dir))

	H_vacuum_bar=fill(zeros(Complex{Float64},(3,3)),length(x_dir))
	H_matter_bar=fill(zeros(Complex{Float64},(3,3)),length(x_dir))
	H_neutrino_bar=fill(zeros(Complex{Float64},(3,3)),length(x_dir))

	s12=sin(theta_12)
	c12=cos(theta_12)
	s13=sin(theta_13)
	c13=cos(theta_13)
	s23=sin(theta_23)
	c23=cos(theta_23)

	U=[c12*c13 s12*c13 s13*e^(-im*delta_cp);-s12*c23-c12*s13*s23*e^(im*delta_cp) c12*c23-s12*s13*s23*e^(im*delta_cp) c13*s23;s12*s23-c12*s13*c23*e^(im*delta_cp) -c12*s23-s12*s13*c23*e^(im*delta_cp) c13*c23]	
	
	U_dagger=conj(transpose(U))
	
	mass=eV_to_J*[mass_1^2/(2*neutrino_energy)+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im mass_2^2/(2*neutrino_energy)+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im mass_3^2/(2*neutrino_energy)+0.0im]
	
	ham_vacuum=U*mass*U_dagger
	
	for i in eachindex(x_dir)
	
		H_vacuum[i]=ham_vacuum
		H_vacuum_bar[i]=conj(ham_vacuum)			
	
		H_matter[i]     = sqrt(2)*Gf*(hbar*c)^3*[(background_electron_number_density-background_antielectron_number_density) 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]	
		H_matter_bar[i] = -conj(H_matter[i])
						
		H_neutrino[i]     = sqrt(2)*Gf*(hbar*c)^3*((N[i]-conj(N_bar[i]))-x_dir[i]*(x_flux[i]-conj(x_flux_bar[i])))
		H_neutrino_bar[i] = -conj(H_neutrino[i])
		
	end
	
	return H_vacuum,H_matter,H_neutrino,H_vacuum_bar,H_matter_bar,H_neutrino_bar

end























function compute_next_step_particle_data(hamiltonians,initial_data,time)

	println("Computing next step particles varibles")							

	# initial_particles_data ---> 1: x particle position
	# initial_particles_data ---> 2: x particle direction
	# initial_particles_data ---> 3: rho neutrinos
	# initial_particles_data ---> 4: rho bar antineutrinos
	# initial_particles_data ---> 5: number of neutrinos
	# initial_particles_data ---> 6: number of antineutrinos

	x=initial_data[1]
	x_dir=initial_data[2]
	rho=initial_data[3]
	rho_bar=initial_data[4]

	H_vacuum=hamiltonians[1]
	H_matter=hamiltonians[2]
	H_neutrino=hamiltonians[3]
	H_vacuum_bar=hamiltonians[4]
	H_matter_bar=hamiltonians[5]
	H_neutrino_bar=hamiltonians[6]

	x_particle_position=fill(0.0,length(x_dir))
	particles_rho=fill(zeros(Complex{Float64},(3,3)),length(x_dir)) 
	particles_rho_bar=fill(zeros(Complex{Float64},(3,3)),length(x_dir))

	for i in eachindex(initial_data[1])

		function position_x_dot(t,x)
			return c*x_dir[i]
		end		
		function rho_dot(t,rho_)
			#return (-im/hbar)*((H_vacuum[i]+H_matter[i]+H_neutrino[i])*rho_-rho_*(H_vacuum[i]+H_matter[i]+H_neutrino[i]))
			#return (-im/hbar)*(H_vacuum[i]*rho_-rho_*H_vacuum[i])
			#return (-im/hbar)*((H_vacuum[i]+H_matter[i])*rho_-rho_*(H_vacuum[i]+H_matter[i]))
			return (-im/hbar)*(H_neutrino[i]*rho_-rho_*H_neutrino[i])
		end
		function rho_bar_dot(t,rho_bar_)
			#return (-im/hbar)*((H_vacuum_bar[i]+H_matter_bar[i]+H_neutrino_bar[i])*rho_bar_-rho_bar_*(H_vacuum_bar[i]+H_matter_bar[i]+H_neutrino_bar[i]))
			return (-im/hbar)*(H_vacuum_bar[i]*rho_bar_-rho_bar_*H_vacuum_bar[i])
			#return (-im/hbar)*((H_vacuum_bar[i]+H_matter_bar[i])*rho_bar_-rho_bar_*(H_vacuum_bar[i]+H_matter_bar[i]))
			return (-im/hbar)*(H_neutrino_bar[i]*rho_bar_-rho_bar_*H_neutrino_bar[i])
		end

		x_particle_position[i] = rk4(position_x_dot,time,x[i],time_step)
		particles_rho[i]       = rk4(rho_dot,time,rho[i],time_step)
		particles_rho_bar[i]   = rk4(rho_bar_dot,time,rho_bar[i],time_step)

		# periodic boundary conditions for particle position position

		if x_particle_position[i]<0.0
			x_particle_position[i]=cell_x_lenght*number_of_cells+x_particle_position[i]
		elseif x_particle_position[i]>=cell_x_lenght*number_of_cells
			x_particle_position[i]=x_particle_position[i]-cell_x_lenght*number_of_cells
		end
	end
	return x_particle_position,x_dir,particles_rho,particles_rho_bar,initial_data[5],initial_data[6]
end

































function evolve_particles(simulation_time,initial_particles_data)

	# initial_particles_data ---> 1: x particle position
	# initial_particles_data ---> 2: x particle direction
	# initial_particles_data ---> 3: rho neutrinos
	# initial_particles_data ---> 4: rho bar antineutrinos 
	# initial_particles_data ---> 5: number of neutrinos
	# initial_particles_data ---> 6: number of antineutrinos

	#writting particles initial condition
	write_info("step_0.txt",initial_particles_data,simulation_time)
	
	println(".......................................")							
	println("Start particle evolution")							
	println(".......................................")							
		
	for a in collect(1:number_of_steps)
				
		print("time: ")
		print(simulation_time)
		print(" step: ")
		println(a)
		
		particles__number_cell=particles_number_cell(initial_particles_data[1])
		# particles__number_cell ---> 1: array of the number of cell the particles bellows
		
		shape_func=shape_function(initial_particles_data[1],cells__center,particles__number_cell)
		# shape_func ---> 1: w previous cell
		# shape_func ---> 2: w actual cell
		# shape_func ---> 3: w nest cell

		grid_info_deposited=deposition_from_particles_to_grid(initial_particles_data,shape_func,particles__number_cell)
		# grid_info_deposited ---> 1: neutrinos number density 
		# grid_info_deposited ---> 2: antineutrinos number density
		# grid_info_deposited ---> 3: x neutrinos number density flux
		# grid_info_deposited ---> 4: x antineutrinos number density flux

		particle_interpolation_data=interpolation_from_grid_to_particles(shape_func,particles__number_cell,grid_info_deposited) 
		# particle_interpolation_data ---> 1: neutrinos number density 
		# particle_interpolation_data ---> 2: antineutrinos number density
		# particle_interpolation_data ---> 3: x neutrinos number density flux
		# particle_interpolation_data ---> 4: x antineutrinos number density flux
		
		hamiltonians=compute_hamiltonians(initial_particles_data[2],particle_interpolation_data)
		# hamiltonians ---> 1: neutrinos vacuum  
		# hamiltonians ---> 2: neutrinos-matter  
		# hamiltonians ---> 3: neutrinos-neutrinos  
		# hamiltonians ---> 4: antineutrinos vacuum  
		# hamiltonians ---> 5: antineutrinos-matter  
		# hamiltonians ---> 6: antineutrinos-antineutrinos  
						
		initial_particles_data=compute_next_step_particle_data(hamiltonians,initial_particles_data,simulation_time)
		
		simulation_time+=time_step

		write_info("step_$a.txt",initial_particles_data,simulation_time)
				
		println(".......................................")							
	end

end

