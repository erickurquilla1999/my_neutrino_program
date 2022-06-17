include("input_parameters.jl")

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

