include("input_parameters.jl")


function rk4(y_i_dot,t_i,y_i,h)

		k1=y_i_dot(t_i,y_i)
		k2=y_i_dot(t_i+0.5*h,y_i+0.5*k1*h)
		k3=y_i_dot(t_i+0.5*h,y_i+0.5*k2*h)
		k4=y_i_dot(t_i+h,y_i+k3*h)

		y_i_p1=y_i+(1/6)*h*(k1+2*k2+2*k3+k4)
	
		return y_i_p1
end
