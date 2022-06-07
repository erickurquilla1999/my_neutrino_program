include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/input_parameters.jl")

# compute cells centers in base of imput parameters
# two ghost zone at the begining and at the end are included

function cells_center()	
	
	println("Computing cell center")	
	
	x_grid_center=fill(0.0,number_of_cells+2)
	
	x_grid_center[1]=-cell_x_lenght/2 #first ghost zone
	x_grid_center[number_of_cells+2]=cell_x_lenght*number_of_cells+cell_x_lenght/2 #last ghost zone

	for i in collect(2:(number_of_cells+1))
		x_grid_center[i]=cell_x_lenght*(i-2)+cell_x_lenght/2
	end

	return x_grid_center

end









# find the number cell a particle bellow
# return an array with the number of cell the particles below

function particles_number_cell(x_par)

	println("Computing number cell the particles below")	

	particle_num_cell=fill(0,length(x_par)) 
	
	for i in eachindex(x_par)

		dec,int=modf((x_par[i]/cell_x_lenght)+2.0)

		if (0.0<=dec<=1.0)&&(2.0<=int<=convert(Float64,(number_of_cells+1.0)))
			particle_num_cell[i]=convert(Int128,int)
		else
			println("error trying to find the particle, the particle start the simulation out of the cells or the time_step is too big")
			exit()	
		end
	end
	return particle_num_cell		
end












# compute the shape function of the particles

function shape_function(x_par,x_cell_center,num_cell)

	println("Computing particles shape function")		

	sf_w    = fill(0.0,length(x_par)) 	
	sf_w_p1 = fill(0.0,length(x_par)) 	
	sf_w_m1 = fill(0.0,length(x_par)) 		

	function shape_function_form(distance)
		sf_a=-1.0/(3*cell_x_lenght/2)^2
		return(sf_a*distance^2+1)		
	end
	
	for i in eachindex(x_par)

		w_cell         = shape_function_form(x_cell_center[num_cell[i]]-x_par[i])
		w_cell_plus_1  = shape_function_form(x_cell_center[num_cell[i]+1]-x_par[i])
		w_cell_minus_1 = shape_function_form(x_cell_center[num_cell[i]-1]-x_par[i])
		
		if w_cell<0.0
			w_cell=0.0
		end
		
		if w_cell_plus_1<0.0
			w_cell_plus_1=0.0
		end
		
		if w_cell_minus_1<0.0
			w_cell_minus_1=0.0
		end
		
		sf_w_m1[i] = w_cell_minus_1/(w_cell+w_cell_plus_1+w_cell_minus_1)		
		sf_w[i]    = w_cell/(w_cell+w_cell_plus_1+w_cell_minus_1)
		sf_w_p1[i] = w_cell_plus_1/(w_cell+w_cell_plus_1+w_cell_minus_1)
	end
	
	return sf_w_m1,sf_w,sf_w_p1
end



























function deposition_from_particles_to_grid(particle_info,shape_func,par_num_cell)
	
	println("Depositing information from particle to grid")							
	
	x_dir=particle_info[2]
	rho=particle_info[3]
	rho_bar=particle_info[4]
	N_neutrinos=particle_info[5]
	N_antineutrinos=particle_info[6]
	wm1=shape_func[1]
	w=shape_func[2]
	wp1=shape_func[3]
	
	nd_grid=fill(zeros(Complex{Float64},(3,3)),number_of_cells+2) 
	nd_bar_grid=fill(zeros(Complex{Float64},(3,3)),number_of_cells+2)
	
	x_flux_grid=fill(zeros(Complex{Float64},(3,3)),number_of_cells+2)
	x_flux_bar_grid=fill(zeros(Complex{Float64},(3,3)),number_of_cells+2)

	for i in eachindex(par_num_cell)

		nd_grid[par_num_cell[i]-1] +=N_neutrinos[i]*wm1[i]*rho[i]
		nd_grid[par_num_cell[i]]   +=N_neutrinos[i]*w[i]*rho[i]
		nd_grid[par_num_cell[i]+1] +=N_neutrinos[i]*wp1[i]*rho[i]
		
		nd_bar_grid[par_num_cell[i]-1] +=N_antineutrinos[i]*wm1[i]*rho_bar[i]
		nd_bar_grid[par_num_cell[i]]   +=N_antineutrinos[i]*w[i]*rho_bar[i]
		nd_bar_grid[par_num_cell[i]+1] +=N_antineutrinos[i]*wp1[i]*rho_bar[i]

		x_flux_grid[par_num_cell[i]-1] +=N_neutrinos[i]*wm1[i]*rho[i]*x_dir[i]
		x_flux_grid[par_num_cell[i]]   +=N_neutrinos[i]*w[i]*rho[i]*x_dir[i]
		x_flux_grid[par_num_cell[i]+1] +=N_neutrinos[i]*wp1[i]*rho[i]*x_dir[i]
						
		x_flux_bar_grid[par_num_cell[i]-1]     +=N_antineutrinos[i]*wm1[i]*rho_bar[i]*x_dir[i]
		x_flux_bar_grid[par_num_cell[i]]       +=N_antineutrinos[i]*w[i]*rho_bar[i]*x_dir[i]
		x_flux_bar_grid[par_num_cell[i]+1]     +=N_antineutrinos[i]*wp1[i]*rho_bar[i]*x_dir[i]
		
	end
	
	#adding ghost zones to valid zones

	nd_grid[1]         +=nd_grid[number_of_cells+1]
	nd_bar_grid[1]     +=nd_bar_grid[number_of_cells+1]
	x_flux_grid[1]     +=x_flux_grid[number_of_cells+1]
	x_flux_bar_grid[1] +=x_flux_bar_grid[number_of_cells+1]
	
	nd_grid[number_of_cells+2]         +=nd_grid[2]
	nd_bar_grid[number_of_cells+2]     +=nd_bar_grid[2]
	x_flux_grid[number_of_cells+2]     +=x_flux_grid[2]
	x_flux_bar_grid[number_of_cells+2] +=x_flux_bar_grid[2]
	
	cell_vol=cell_x_lenght^3
	
	nd_grid=nd_grid/cell_vol		
	nd_bar_grid=nd_bar_grid/cell_vol	
	x_flux_grid=x_flux_grid/cell_vol
	x_flux_bar_grid=x_flux_bar_grid/cell_vol

	return nd_grid,nd_bar_grid,x_flux_grid,x_flux_bar_grid

end






















function interpolation_from_grid_to_particles(shape_func,par_num_cell,data_from_grid) 

	println("Depositing information from grid to particles")							
	
	wm1=shape_func[1]
	w=shape_func[2]
	wp1=shape_func[3]
	
	num_den_grid=data_from_grid[1]
	num_den_bar_grid=data_from_grid[2]
	
	x_flux_grid=data_from_grid[3]
	x_flux_bar_grid=data_from_grid[4]

	number_den_particles     = fill(zeros(Complex{Float64},(3,3)),number_of_cells*2)
	number_den_bar_particles = fill(zeros(Complex{Float64},(3,3)),number_of_cells*2)
	
	x_flux_particles         = fill(zeros(Complex{Float64},(3,3)),number_of_cells*2)
	x_flux_bar_particles     = fill(zeros(Complex{Float64},(3,3)),number_of_cells*2)
	
	for i in eachindex(par_num_cell)
	
		number_den_particles[i]+=wm1[i]*num_den_grid[par_num_cell[i]-1]
		number_den_particles[i]+=w[i]*num_den_grid[par_num_cell[i]]
		number_den_particles[i]+=wp1[i]*num_den_grid[par_num_cell[i]+1]
		
		number_den_bar_particles[i]+=wm1[i]*num_den_bar_grid[par_num_cell[i]-1]
		number_den_bar_particles[i]+=w[i]*num_den_bar_grid[par_num_cell[i]]
		number_den_bar_particles[i]+=wp1[i]*num_den_bar_grid[par_num_cell[i]+1]
		
		x_flux_particles[i]+=wm1[i]*x_flux_grid[par_num_cell[i]-1]
		x_flux_particles[i]+=w[i]*x_flux_grid[par_num_cell[i]]
		x_flux_particles[i]+=wp1[i]*x_flux_grid[par_num_cell[i]+1]
		
		x_flux_bar_particles[i]+=wm1[i]*x_flux_bar_grid[par_num_cell[i]-1]
		x_flux_bar_particles[i]+=w[i]*x_flux_bar_grid[par_num_cell[i]]
		x_flux_bar_particles[i]+=wp1[i]*x_flux_bar_grid[par_num_cell[i]+1]
	
	end
	
		println(number_den_particles[4][1,1])		
		println(number_den_bar_particles[4][1,1])
		println(x_flux_particles[4][1,1])
		println(x_flux_bar_particles[4][1,1])		

	return number_den_particles,number_den_bar_particles,x_flux_particles,x_flux_bar_particles

end
