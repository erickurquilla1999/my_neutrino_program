include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/input_parameters.jl")

# compute cells centers in base of imput parameters

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
				
		if (0.0<=dec<=1.0)&&(1.0<=int<=convert(Float64,number_of_cells))
			particle_num_cell[i]=convert(Int128,int)
		else
			println("error trying to find the particle, the time_step is too big")
			exit()	
		end
	end
	return particle_num_cell		
end


3










# compute the shape function of the particles

function shape_function(x_par,x_cell_center,num_cell)

	println("Computing particles shape function")		

	sf_w    = fill(0.0,length(x_par)) 	
	sf_w_p1 = fill(0.0,length(x_par)) 	
	sf_w_m1 = fill(0.0,length(x_par)) 		

	function shape_function_form(distance)
		sf_a=-1.0/cell_x_lenght^2
		sf_c=1.0
		return(sf_a*distance^2+sf_c)		
	end
	
	for i in eachindex(x_par)

		#distance from particle to the center of its grid cell
		d_cell=x_cell_center[num_cell[i]]-x_par[i]

		#distance from particle to the center of the next grid cell
		d_cell_p1=x_cell_center[num_cell[i]+1]-x_par[i]

		#distance from particle to the center of the last grid cell		
		d_cell_m1=x_cell_center[num_cell[i]-1]-x_par[i]

	
		
		w_cell         = shape_function_form(d_cell)
		w_cell_plus_1  = shape_function_form(d_cell)
		w_cell_minus_1 = shape_function_form(d_cell)
		
		if w_cell<0.0
			w_cell=0.0
		end
		
		if w_cell_plus_1<0.0
			w_cell_plus_1=0.0
		end
		
		if w_cell_minus_1<0.0
			w_cell_minus_1=0.0
		end
		
		push!(sf_w_m1,w_cell_minus_1/(w_cell+w_cell_plus_1+w_cell_minus_1))		
		push!(sf_w,w_cell/(w_cell+w_cell_plus_1+w_cell_minus_1))
		push!(sf_w_p1,w_cell_plus_1/(w_cell+w_cell_plus_1+w_cell_minus_1))		
	end
	
	return sf_w_m1,sf_w,sf_w_p1
end




























function deposition_from_particles_to_grid(x_dir,y_dir,z_dir,rho,rho_bar,N_neutrinos,N_antineutrinos,wm1,w,wp1,par_num_cell)

	println("Depositing information from particle to grid")							
	
	nd_grid=Array{Array{Complex{Float64},2},1}(undef,number_of_cells)
	nd_bar_grid=Array{Array{Complex{Float64},2},1}(undef,number_of_cells)
	
	x_flux_grid=Array{Array{Complex{Float64},2},1}(undef,number_of_cells)
	y_flux_grid=Array{Array{Complex{Float64},2},1}(undef,number_of_cells)
	z_flux_grid=Array{Array{Complex{Float64},2},1}(undef,number_of_cells)
	
	x_flux_bar_grid=Array{Array{Complex{Float64},2},1}(undef,number_of_cells)
	y_flux_bar_grid=Array{Array{Complex{Float64},2},1}(undef,number_of_cells)
	z_flux_bar_grid=Array{Array{Complex{Float64},2},1}(undef,number_of_cells)
	
	for i in eachindex(nd_grid)
		nd_grid[i]     =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
		nd_bar_grid[i] =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
		
		x_flux_grid[i]     =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
		y_flux_grid[i]	   =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
		z_flux_grid[i]	   =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
		
		x_flux_bar_grid[i] =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
		y_flux_bar_grid[i] =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
		z_flux_bar_grid[i] =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]	
	end
	
	for i in eachindex(par_num_cell)
		if par_num_cell[i]==1
			nd_grid[number_of_cells] +=N_neutrinos[i]*wm1[i]*rho[i]
			nd_grid[1]               +=N_neutrinos[i]*w[i]*rho[i]
			nd_grid[2]       	+=N_neutrinos[i]*wp1[i]*rho[i]
			
			nd_bar_grid[number_of_cells] +=N_antineutrinos[i]*wm1[i]*rho_bar[i]
			nd_bar_grid[1]       	    +=N_antineutrinos[i]*w[i]*rho_bar[i]
			nd_bar_grid[2]       	    +=N_antineutrinos[i]*wp1[i]*rho_bar[i]

			x_flux_grid[number_of_cells] +=N_neutrinos[i]*wm1[i]*rho[i]*x_dir[i]
			x_flux_grid[1]               +=N_neutrinos[i]*w[i]*rho[i]*x_dir[i]
			x_flux_grid[2]       	     +=N_neutrinos[i]*wp1[i]*rho[i]*x_dir[i]
			
			y_flux_grid[number_of_cells] +=N_neutrinos[i]*wm1[i]*rho[i]*y_dir[i]
			y_flux_grid[1]               +=N_neutrinos[i]*w[i]*rho[i]*y_dir[i]
			y_flux_grid[2]       	     +=N_neutrinos[i]*wp1[i]*rho[i]*y_dir[i]
			
			z_flux_grid[number_of_cells] +=N_neutrinos[i]*wm1[i]*rho[i]*z_dir[i]
			z_flux_grid[1]               +=N_neutrinos[i]*w[i]*rho[i]*z_dir[i]
			z_flux_grid[2]       	     +=N_neutrinos[i]*wp1[i]*rho[i]*z_dir[i]
									
			x_flux_bar_grid[number_of_cells]     +=N_antineutrinos[i]*wm1[i]*rho_bar[i]*x_dir[i]
			x_flux_bar_grid[1]       	     +=N_antineutrinos[i]*w[i]*rho_bar[i]*x_dir[i]
			x_flux_bar_grid[2]       	     +=N_antineutrinos[i]*wp1[i]*rho_bar[i]*x_dir[i]
			
			y_flux_bar_grid[number_of_cells]     +=N_antineutrinos[i]*wm1[i]*rho_bar[i]*y_dir[i]
			y_flux_bar_grid[1]       	     +=N_antineutrinos[i]*w[i]*rho_bar[i]*y_dir[i]
			y_flux_bar_grid[2]       	     +=N_antineutrinos[i]*wp1[i]*rho_bar[i]*y_dir[i]
			
			z_flux_bar_grid[number_of_cells]     +=N_antineutrinos[i]*wm1[i]*rho_bar[i]*z_dir[i]
			z_flux_bar_grid[1]       	     +=N_antineutrinos[i]*w[i]*rho_bar[i]*z_dir[i]
			z_flux_bar_grid[2]       	     +=N_antineutrinos[i]*wp1[i]*rho_bar[i]*z_dir[i]
			
		elseif par_num_cell[i]==number_of_cells
			
			nd_grid[par_num_cell[i]-1] +=N_neutrinos[i]*wm1[i]*rho[i]
			nd_grid[par_num_cell[i]]   +=N_neutrinos[i]*w[i]*rho[i]
			nd_grid[1]                 +=N_neutrinos[i]*wp1[i]*rho[i]
			
			nd_bar_grid[par_num_cell[i]-1] +=N_antineutrinos[i]*wm1[i]*rho_bar[i]
			nd_bar_grid[par_num_cell[i]]   +=N_antineutrinos[i]*w[i]*rho_bar[i]
			nd_bar_grid[1]                 +=N_antineutrinos[i]*wp1[i]*rho_bar[i]

			x_flux_grid[par_num_cell[i]-1] +=N_neutrinos[i]*wm1[i]*rho[i]*x_dir[i]
			x_flux_grid[par_num_cell[i]]   +=N_neutrinos[i]*w[i]*rho[i]*x_dir[i]
			x_flux_grid[1]       	       +=N_neutrinos[i]*wp1[i]*rho[i]*x_dir[i]
			
			y_flux_grid[par_num_cell[i]-1] +=N_neutrinos[i]*wm1[i]*rho[i]*y_dir[i]
			y_flux_grid[par_num_cell[i]]   +=N_neutrinos[i]*w[i]*rho[i]*y_dir[i]
			y_flux_grid[1]       	       +=N_neutrinos[i]*wp1[i]*rho[i]*y_dir[i]
			
			z_flux_grid[par_num_cell[i]-1] +=N_neutrinos[i]*wm1[i]*rho[i]*z_dir[i]
			z_flux_grid[par_num_cell[i]]   +=N_neutrinos[i]*w[i]*rho[i]*z_dir[i]
			z_flux_grid[1]       	       +=N_neutrinos[i]*wp1[i]*rho[i]*z_dir[i]
									
			x_flux_bar_grid[par_num_cell[i]-1]   +=N_antineutrinos[i]*wm1[i]*rho_bar[i]*x_dir[i]
			x_flux_bar_grid[par_num_cell[i]]     +=N_antineutrinos[i]*w[i]*rho_bar[i]*x_dir[i]
			x_flux_bar_grid[1]       	     +=N_antineutrinos[i]*wp1[i]*rho_bar[i]*x_dir[i]
			
			y_flux_bar_grid[par_num_cell[i]-1]   +=N_antineutrinos[i]*wm1[i]*rho_bar[i]*y_dir[i]
			y_flux_bar_grid[par_num_cell[i]]     +=N_antineutrinos[i]*w[i]*rho_bar[i]*y_dir[i]
			y_flux_bar_grid[1]       	     +=N_antineutrinos[i]*wp1[i]*rho_bar[i]*y_dir[i]
			
			z_flux_bar_grid[par_num_cell[i]-1]   +=N_antineutrinos[i]*wm1[i]*rho_bar[i]*z_dir[i]
			z_flux_bar_grid[par_num_cell[i]]     +=N_antineutrinos[i]*w[i]*rho_bar[i]*z_dir[i]
			z_flux_bar_grid[1]       	     +=N_antineutrinos[i]*wp1[i]*rho_bar[i]*z_dir[i]
		else
			nd_grid[par_num_cell[i]-1] +=N_neutrinos[i]*wm1[i]*rho[i]
			nd_grid[par_num_cell[i]]   +=N_neutrinos[i]*w[i]*rho[i]
			nd_grid[par_num_cell[i]+1] +=N_neutrinos[i]*wp1[i]*rho[i]
			
			nd_bar_grid[par_num_cell[i]-1] +=N_antineutrinos[i]*wm1[i]*rho_bar[i]
			nd_bar_grid[par_num_cell[i]]   +=N_antineutrinos[i]*w[i]*rho_bar[i]
			nd_bar_grid[par_num_cell[i]+1] +=N_antineutrinos[i]*wp1[i]*rho_bar[i]

			x_flux_grid[par_num_cell[i]-1] +=N_neutrinos[i]*wm1[i]*rho[i]*x_dir[i]
			x_flux_grid[par_num_cell[i]]   +=N_neutrinos[i]*w[i]*rho[i]*x_dir[i]
			x_flux_grid[par_num_cell[i]+1] +=N_neutrinos[i]*wp1[i]*rho[i]*x_dir[i]
			
			y_flux_grid[par_num_cell[i]-1] +=N_neutrinos[i]*wm1[i]*rho[i]*y_dir[i]
			y_flux_grid[par_num_cell[i]]   +=N_neutrinos[i]*w[i]*rho[i]*y_dir[i]
			y_flux_grid[par_num_cell[i]+1] +=N_neutrinos[i]*wp1[i]*rho[i]*y_dir[i]
			
			z_flux_grid[par_num_cell[i]-1] +=N_neutrinos[i]*wm1[i]*rho[i]*z_dir[i]
			z_flux_grid[par_num_cell[i]]   +=N_neutrinos[i]*w[i]*rho[i]*z_dir[i]
			z_flux_grid[par_num_cell[i]+1] +=N_neutrinos[i]*wp1[i]*rho[i]*z_dir[i]
									
			x_flux_bar_grid[par_num_cell[i]-1]     +=N_antineutrinos[i]*wm1[i]*rho_bar[i]*x_dir[i]
			x_flux_bar_grid[par_num_cell[i]]       +=N_antineutrinos[i]*w[i]*rho_bar[i]*x_dir[i]
			x_flux_bar_grid[par_num_cell[i]+1]     +=N_antineutrinos[i]*wp1[i]*rho_bar[i]*x_dir[i]
			
			y_flux_bar_grid[par_num_cell[i]-1]     +=N_antineutrinos[i]*wm1[i]*rho_bar[i]*y_dir[i]
			y_flux_bar_grid[par_num_cell[i]]       +=N_antineutrinos[i]*w[i]*rho_bar[i]*y_dir[i]
			y_flux_bar_grid[par_num_cell[i]+1]     +=N_antineutrinos[i]*wp1[i]*rho_bar[i]*y_dir[i]
			
			z_flux_bar_grid[par_num_cell[i]-1]     +=N_antineutrinos[i]*wm1[i]*rho_bar[i]*z_dir[i]
			z_flux_bar_grid[par_num_cell[i]]       +=N_antineutrinos[i]*w[i]*rho_bar[i]*z_dir[i]
			z_flux_bar_grid[par_num_cell[i]+1]     +=N_antineutrinos[i]*wp1[i]*rho_bar[i]*z_dir[i]		

			z_flux_bar_grid[par_num_cell[i]-1]   +=N_antineutrinos[i]*wm1[i]*rho_bar[i]*z_dir[i]
			z_flux_bar_grid[par_num_cell[i]]     +=N_antineutrinos[i]*w[i]*rho_bar[i]*z_dir[i]
			z_flux_bar_grid[1]       	     +=N_antineutrinos[i]*wp1[i]*rho_bar[i]*z_dir[i]
		end
	end
	
	cell_vol=cell_x_lenght*cell_y_lenght*cell_z_lenght
	
	nd_grid=nd_grid/cell_vol		
	nd_bar_grid=nd_bar_grid/cell_vol	
	x_flux_grid=x_flux_grid/cell_vol
	y_flux_grid=y_flux_grid/cell_vol
	z_flux_grid=z_flux_grid/cell_vol
	x_flux_bar_grid=x_flux_bar_grid/cell_vol
	y_flux_bar_grid=y_flux_bar_grid/cell_vol
	z_flux_bar_grid=z_flux_bar_grid/cell_vol
	
	return nd_grid,nd_bar_grid,x_flux_grid,y_flux_grid,z_flux_grid,x_flux_bar_grid,y_flux_bar_grid,z_flux_bar_grid
end






















function interpolation_from_grid_to_particles(wm1,w,wp1,par_num_cell,data_from_grid) 

	println("Depositing information from grid to particles")							

	num_den_grid=data_from_grid[1]
	num_den_bar_grid=data_from_grid[2]
	
	x_flux_grid=data_from_grid[3]
	y_flux_grid=data_from_grid[4]
	z_flux_grid=data_from_grid[5]
	
	x_flux_bar_grid=data_from_grid[6]
	y_flux_bar_grid=data_from_grid[7]
	z_flux_bar_grid=data_from_grid[8]

	number_den_particles    =Array{Array{Complex{Float64},2},1}(undef,2*number_of_cells)
	number_den_bar_particles=Array{Array{Complex{Float64},2},1}(undef,2*number_of_cells)
	
	x_flux_particles        =Array{Array{Complex{Float64},2},1}(undef,2*number_of_cells)
	y_flux_particles        =Array{Array{Complex{Float64},2},1}(undef,2*number_of_cells)
	z_flux_particles        =Array{Array{Complex{Float64},2},1}(undef,2*number_of_cells)
	
	x_flux_bar_particles    =Array{Array{Complex{Float64},2},1}(undef,2*number_of_cells)
	y_flux_bar_particles    =Array{Array{Complex{Float64},2},1}(undef,2*number_of_cells)
	z_flux_bar_particles    =Array{Array{Complex{Float64},2},1}(undef,2*number_of_cells)
	
	for i in eachindex(number_den_particles)	
	
		number_den_particles[i]    =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]	
		number_den_bar_particles[i]=[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]	
		
		x_flux_particles[i]        =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
		y_flux_particles[i]	   =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
		z_flux_particles[i]	   =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
		
		x_flux_bar_particles[i]    =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
		y_flux_bar_particles[i]    =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
		z_flux_bar_particles[i]    =[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
	
	end


	for i in eachindex(par_num_cell)
		
		if par_num_cell[i]==1
			number_den_particles[i]+=wm1[i]*num_den_grid[number_of_cells]
			number_den_particles[i]+=w[i]*num_den_grid[par_num_cell[i]]
			number_den_particles[i]+=wp1[i]*num_den_grid[par_num_cell[i]+1]
			
			number_den_bar_particles[i]+=wm1[i]*num_den_bar_grid[number_of_cells]
			number_den_bar_particles[i]+=w[i]*num_den_bar_grid[par_num_cell[i]]
			number_den_bar_particles[i]+=wp1[i]*num_den_bar_grid[par_num_cell[i]+1]

			x_flux_particles[i]+=wm1[i]*x_flux_grid[number_of_cells]
			x_flux_particles[i]+=w[i]*x_flux_grid[par_num_cell[i]]
			x_flux_particles[i]+=wp1[i]*x_flux_grid[par_num_cell[i]+1]

			y_flux_particles[i]+=wm1[i]*y_flux_grid[number_of_cells]
			y_flux_particles[i]+=w[i]*y_flux_grid[par_num_cell[i]]
			y_flux_particles[i]+=wp1[i]*y_flux_grid[par_num_cell[i]+1]

			z_flux_particles[i]+=wm1[i]*z_flux_grid[number_of_cells]
			z_flux_particles[i]+=w[i]*z_flux_grid[par_num_cell[i]]
			z_flux_particles[i]+=wp1[i]*z_flux_grid[par_num_cell[i]+1]
			
			x_flux_bar_particles[i]+=wm1[i]*x_flux_bar_grid[number_of_cells]
			x_flux_bar_particles[i]+=w[i]*x_flux_bar_grid[par_num_cell[i]]
			x_flux_bar_particles[i]+=wp1[i]*x_flux_bar_grid[par_num_cell[i]+1]
			
			y_flux_bar_particles[i]+=wm1[i]*y_flux_bar_grid[number_of_cells]
			y_flux_bar_particles[i]+=w[i]*y_flux_bar_grid[par_num_cell[i]]
			y_flux_bar_particles[i]+=wp1[i]*y_flux_bar_grid[par_num_cell[i]+1]

			z_flux_bar_particles[i]+=wm1[i]*z_flux_bar_grid[number_of_cells]
			z_flux_bar_particles[i]+=w[i]*z_flux_bar_grid[par_num_cell[i]]
			z_flux_bar_particles[i]+=wp1[i]*z_flux_bar_grid[par_num_cell[i]+1]
																	
		elseif par_num_cell[i]==number_of_cells
			number_den_particles[i]+=wm1[i]*num_den_grid[par_num_cell[i]-1]
			number_den_particles[i]+=w[i]*num_den_grid[par_num_cell[i]]
			number_den_particles[i]+=wp1[i]*num_den_grid[1]
			
			number_den_bar_particles[i]+=wm1[i]*num_den_bar_grid[par_num_cell[i]-1]
			number_den_bar_particles[i]+=w[i]*num_den_bar_grid[par_num_cell[i]]
			number_den_bar_particles[i]+=wp1[i]*num_den_bar_grid[1]

			x_flux_particles[i]+=wm1[i]*x_flux_grid[par_num_cell[i]-1]
			x_flux_particles[i]+=w[i]*x_flux_grid[par_num_cell[i]]
			x_flux_particles[i]+=wp1[i]*x_flux_grid[1]

			y_flux_particles[i]+=wm1[i]*y_flux_grid[par_num_cell[i]-1]
			y_flux_particles[i]+=w[i]*y_flux_grid[par_num_cell[i]]
			y_flux_particles[i]+=wp1[i]*y_flux_grid[1]

			z_flux_particles[i]+=wm1[i]*z_flux_grid[par_num_cell[i]-1]
			z_flux_particles[i]+=w[i]*z_flux_grid[par_num_cell[i]]
			z_flux_particles[i]+=wp1[i]*z_flux_grid[1]
			
			x_flux_bar_particles[i]+=wm1[i]*x_flux_bar_grid[par_num_cell[i]-1]
			x_flux_bar_particles[i]+=w[i]*x_flux_bar_grid[par_num_cell[i]]
			x_flux_bar_particles[i]+=wp1[i]*x_flux_bar_grid[1]
			
			y_flux_bar_particles[i]+=wm1[i]*y_flux_bar_grid[par_num_cell[i]-1]
			y_flux_bar_particles[i]+=w[i]*y_flux_bar_grid[par_num_cell[i]]
			y_flux_bar_particles[i]+=wp1[i]*y_flux_bar_grid[1]

			z_flux_bar_particles[i]+=wm1[i]*z_flux_bar_grid[par_num_cell[i]-1]
			z_flux_bar_particles[i]+=w[i]*z_flux_bar_grid[par_num_cell[i]]
			z_flux_bar_particles[i]+=wp1[i]*z_flux_bar_grid[1]			 
		else
			number_den_particles[i]+=wm1[i]*num_den_grid[par_num_cell[i]-1]
			number_den_particles[i]+=w[i]*num_den_grid[par_num_cell[i]]
			number_den_particles[i]+=wp1[i]*num_den_grid[par_num_cell[i]+1]
			
			number_den_bar_particles[i]+=wm1[i]*num_den_bar_grid[par_num_cell[i]-1]
			number_den_bar_particles[i]+=w[i]*num_den_bar_grid[par_num_cell[i]]
			number_den_bar_particles[i]+=wp1[i]*num_den_bar_grid[par_num_cell[i]+1]
			
			x_flux_particles[i]+=wm1[i]*x_flux_grid[par_num_cell[i]-1]
			x_flux_particles[i]+=w[i]*x_flux_grid[par_num_cell[i]]
			x_flux_particles[i]+=wp1[i]*x_flux_grid[par_num_cell[i]+1]

			y_flux_particles[i]+=wm1[i]*y_flux_grid[par_num_cell[i]-1]
			y_flux_particles[i]+=w[i]*y_flux_grid[par_num_cell[i]]
			y_flux_particles[i]+=wp1[i]*y_flux_grid[par_num_cell[i]+1]

			z_flux_particles[i]+=wm1[i]*z_flux_grid[par_num_cell[i]-1]
			z_flux_particles[i]+=w[i]*z_flux_grid[par_num_cell[i]]
			z_flux_particles[i]+=wp1[i]*z_flux_grid[par_num_cell[i]+1]
			
			x_flux_bar_particles[i]+=wm1[i]*x_flux_bar_grid[par_num_cell[i]-1]
			x_flux_bar_particles[i]+=w[i]*x_flux_bar_grid[par_num_cell[i]]
			x_flux_bar_particles[i]+=wp1[i]*x_flux_bar_grid[par_num_cell[i]+1]
			
			y_flux_bar_particles[i]+=wm1[i]*y_flux_bar_grid[par_num_cell[i]-1]
			y_flux_bar_particles[i]+=w[i]*y_flux_bar_grid[par_num_cell[i]]
			y_flux_bar_particles[i]+=wp1[i]*y_flux_bar_grid[par_num_cell[i]+1]

			z_flux_bar_particles[i]+=wm1[i]*z_flux_bar_grid[par_num_cell[i]-1]
			z_flux_bar_particles[i]+=w[i]*z_flux_bar_grid[par_num_cell[i]]
			z_flux_bar_particles[i]+=wp1[i]*z_flux_bar_grid[par_num_cell[i]+1]
			
		end
	end
	
	return number_den_particles,number_den_bar_particles,x_flux_particles,y_flux_particles,z_flux_particles,x_flux_bar_particles,y_flux_bar_particles,z_flux_bar_particles
end
















