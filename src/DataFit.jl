function fit_SAIM(file_path::String, file_name::String, optic::SAIMOptics, angles::AbstractArray, 
	init_params::AbstractArray, lower_bounds::AbstractArray, upper_bounds::AbstractArray; glb::Bool=false, 
	step::Float64=40.0, disp::Bool=false, color::Int64=1)
    
    angles_rad = angles*2*pi/360  	
    
    if color == 1
        if glb	#Run a global grid search when true
            fit_1c_global(file_path, file_name, optic, angles_rad, init_params, lower_bounds, 
            upper_bounds, step, disp)
        else	#Conduct a local search
            fit_1c_local(file_path, file_name, optic, angles_rad, init_params, lower_bounds, 
            upper_bounds, disp)
        end
    elseif color == 2
        if glb	#Run a global grid search when true
            fit_2c_global(file_path, file_name, optic, angles_rad, init_params, lower_bounds, 
            upper_bounds, step, disp)
            println("2c global")
        else	#Conduct a local search
            fit_2c_local(file_path, file_name, optic, angles_rad, init_params, lower_bounds, 
            upper_bounds, disp)
            println("2c local")
        end
    elseif color == 3
        #=
        if glb	#Run a global grid search when true
            fit_3c_global(file_path, file_name, optic1, optic2, optic3, angles_rad, init_params, 
            lower_bounds, upper_bounds, step, disp)
        else	#Conduct a local search
            fit_3c_local(file_path, file_name, optic1, optic2, optic3, angles_rad, init_params, 
            lower_bounds, upper_bounds, disp)
        end
        =#
    else
        println("Only 1, 2, and 3 colors supported")
    end
end
