function fit_SAIM(file_path::String, file_name::String, optic::SAIMOptics, angles::AbstractArray, 
	init_params::AbstractArray, lower_bounds::AbstractArray, upper_bounds::AbstractArray; glb::Bool=false, 
	step::Float64=40.0, show::Bool=false, color::Int64=1)
    
    if color == 1
        if glb	#Run a global grid search when true
            fit_1c_global(file_path, file_name, optic, angles, init_params, lower_bounds, 
            upper_bounds, step, show)
        else	#Conduct a local search
            fit_1c_local(file_path, file_name, optic, angles, init_params, lower_bounds, 
            upper_bounds, show)
        end
    elseif color == 2
        if glb	#Run a global grid search when true
            fit_2c_global(file_path, file_name, optic, angles, init_params, lower_bounds, 
            upper_bounds, step, show)
        else	#Conduct a local search
            fit_2c_local(file_path, file_name, optic, angles, init_params, lower_bounds, 
            upper_bounds, show)
        end
    else
        println("Only 1 and 2 color SAIM is currently supported.  Three color SAIM is under development")
    end
end
