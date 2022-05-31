include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/input_parameters.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/constants.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/numerical_methods.jl")









function compute_hamiltonians(x_dir,y_dir,z_dir,interpotation_data)

	println("Computing hamiltonians")							

	H_vacuum=Array{Array{Complex{Float64},2},1}(undef,0)
	H_matter=Array{Array{Complex{Float64},2},1}(undef,0)
	H_neutrino=Array{Array{Complex{Float64},2},1}(undef,0)

	H_vacuum_bar=Array{Array{Complex{Float64},2},1}(undef,0)
	H_matter_bar=Array{Array{Complex{Float64},2},1}(undef,0)
	H_neutrino_bar=Array{Array{Complex{Float64},2},1}(undef,0)

	s12=sin(theta_12)
	c12=cos(theta_12)
	s13=sin(theta_13)
	c13=cos(theta_13)
	s23=sin(theta_23)
	c23=cos(theta_23)
	sdcp=sin(delta_cp)
	cdcp=cos(delta_cp)

	U=[c12*c13+0.0im s12*c13+0.0im s13*cdcp-s13*sdcp*im;-s12*c23-c12*s13*s23*cdcp-c12*s13*s23*sdcp*im c12*c23-s12*s13*s23*cdcp-s12*s13*s23*sdcp*im c13*s23+0.0im;s12*s23-c12*s13*c23*cdcp-c12*s13*c23*sdcp*im -c12*s23-s12*s13*c23*cdcp-s12*s13*c23*sdcp*im c13*c23+0.0im]
	U_dagger=conj(transpose(U))
	mass=eV_to_J*[mass_1^2/neutrino_energy+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im mass_2^2/neutrino_energy+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im mass_3^2/neutrino_energy+0.0im]
	ham_vacuum=U*mass*U_dagger
	
	for i in eachindex(interpotation_data[1])
	
		push!(H_vacuum,ham_vacuum)
		push!(H_vacuum_bar,conj(ham_vacuum))				
	
		ham_neutrinos=sqrt(2)*Gf*(hbar*c)^3*((interpotation_data[1][i]-conj(interpotation_data[2][i]))-x_dir[i]*(interpotation_data[3][i]-conj(interpotation_data[6][i]))-y_dir[i]*(interpotation_data[4][i]-conj(interpotation_data[7][i]))-z_dir[i]*(interpotation_data[5][i]-conj(interpotation_data[8][i])))
		
		push!(H_matter,[ham_neutrinos[1,1] 0.0+0.0im 0.0+0.0im;0.0+0.0im ham_neutrinos[2,2] 0.0+0.0im;0.0+0.0im 0.0+0.0im ham_neutrinos[3,3]])
		push!(H_matter_bar,-conj(H_matter[i]))
						
		push!(H_neutrino,ham_neutrinos)
		push!(H_neutrino_bar,-conj(ham_neutrinos))
	end

	return H_vacuum,H_matter,H_neutrino,H_vacuum_bar,H_matter_bar,H_neutrino_bar

end























function compute_next_step_particle_data(hamiltonians,initial_data,time)

	println("Computing next step particles varibles")							

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

	x=initial_data[1]
	y=initial_data[2]
	z=initial_data[3]
	rho=initial_data[7]
	rho_bar=initial_data[8]

	H_vacuum=hamiltonians[1]
	H_matter=hamiltonians[2]
	H_neutrino=hamiltonians[3]
	H_vacuum_bar=hamiltonians[4]
	H_matter_bar=hamiltonians[5]
	H_neutrino_bar=hamiltonians[6]

	x_particle_position=Array{Float64,1}(undef,0)
	y_particle_position=Array{Float64,1}(undef,0)
	z_particle_position=Array{Float64,1}(undef,0)
	particles_rho=Array{Array{Complex{Float64},2},1}(undef,0)  
	particles_rho_bar=Array{Array{Complex{Float64},2},1}(undef,0)

	for i in eachindex(initial_data[1])
		
		function position_x_dot(t,x)
			return c*initial_data[4][i]
		end
		function position_y_dot(t,y)
			return c*initial_data[5][i]
		end
		function position_z_dot(t,z)
			return c*initial_data[6][i]
		end		
		function rho_dot(t,rho_)
			return (-im/hbar)*((H_vacuum[i]+H_matter[i]+H_neutrino[i])*rho_-rho_*(H_vacuum[i]+H_matter[i]+H_neutrino[i]))
		end
		function rho_bar_dot(t,rho_bar)
			return (-im/hbar)*((H_vacuum_bar[i]+H_matter_bar[i]+H_neutrino_bar[i])*rho_bar-rho_bar*(H_vacuum_bar[i]+H_matter_bar[i]+H_neutrino_bar[i]))
		end

		push!(x_particle_position,rk4(position_x_dot,time,x[i],time_step))
		push!(y_particle_position,rk4(position_y_dot,time,y[i],time_step))	
		push!(z_particle_position,rk4(position_z_dot,time,z[i],time_step))
		push!(particles_rho,rk4(rho_dot,time,rho[i],time_step))	
		push!(particles_rho_bar,rk4(rho_bar_dot,time,rho_bar[i],time_step))	

		# periodic boundary conditions for particle position position

		if x_particle_position[i]<0.0
			x_particle_position[i]=cell_x_lenght*number_of_cells+x_particle_position[i]
		elseif x_particle_position[i]>=cell_x_lenght*number_of_cells
			x_particle_position[i]=x_particle_position[i]-cell_x_lenght*number_of_cells
		end
				
		if y_particle_position[i]<0.0
			y_particle_position[i]=cell_y_lenght+y_particle_position[i]
		elseif y_particle_position[i]>=cell_y_lenght
			y_particle_position[i]=y_particle_position[i]-cell_y_lenght
		end

		if z_particle_position[i]<0.0
			z_particle_position[i]=cell_z_lenght+z_particle_position[i]
		elseif z_particle_position[i]>=cell_z_lenght
			z_particle_position[i]=z_particle_position[i]-cell_z_lenght
		end
	
	end
	return x_particle_position,y_particle_position,z_particle_position,initial_data[4],initial_data[5],initial_data[6],particles_rho,particles_rho_bar,initial_data[9],initial_data[10]
end


















function evolve_particles(simulation_time,initial_particles_data)

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

	#writting particles initial condition
	write_info("step_0.txt",initial_particles_data,simulation_time)
	
	println(".......................................")							
	println("Start particle evolution")							
	println(".......................................")							
		

		
	for a in collect(1:number_of_steps)
	
		write_info("step_$a.txt",initial_particles_data,simulation_time)
			
		print("time: ")
		print(simulation_time)
		print(" step: ")
		println(a)
		
		particles__number_cell=particles_number_cell(initial_particles_data[1],initial_particles_data[2],initial_particles_data[3])
		# particles__number_cell ---> 1: array of the number of cell the particles bellows
		
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

		particle_interpolation_data=interpolation_from_grid_to_particles(shape_func[1],shape_func[2],shape_func[3],particles__number_cell,grid_info_deposited) 
		# particle_interpolation_data ---> 1: neutrinos number density 
		# particle_interpolation_data ---> 2: antineutrinos number density
		# particle_interpolation_data ---> 3: x neutrinos number density flux
		# particle_interpolation_data ---> 4: y neutrinos number density flux
		# particle_interpolation_data ---> 5: z neutrinos number density flux
		# particle_interpolation_data ---> 6: x antineutrinos number density flux
		# particle_interpolation_data ---> 7: y antineutrinos number density flux
		# particle_interpolation_data ---> 8: z antineutrinos number density flux
		
		hamiltonians=compute_hamiltonians(initial_particles_data[4],initial_particles_data[5],initial_particles_data[6],particle_interpolation_data)
		# hamiltonians ---> 1: neutrinos vacuum  
		# hamiltonians ---> 2: neutrinos-matter  
		# hamiltonians ---> 3: neutrinos-neutrinos  
		# hamiltonians ---> 4: antineutrinos vacuum  
		# hamiltonians ---> 5: antineutrinos-matter  
		# hamiltonians ---> 6: antineutrinos-antineutrinos  
						
		initial_particles_data=compute_next_step_particle_data(hamiltonians,initial_particles_data,simulation_time)
		
		simulation_time+=time_step
		
		println(".......................................")							
	end

end

