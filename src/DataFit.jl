function fit_SAIM(file_path::String, file_name::String, optic1::SAIMOptics, angles::AbstractArray, 
	init_params::AbstractArray, lower_bounds::AbstractArray, upper_bounds::AbstractArray; glb::Bool=false, 
	step::Float64=40.0, disp::Bool=false)
    
    angles_rad = angles*2*pi/360  	
    
    if glb	#Run a global grid search when true
        fit_1c_global(file_path, file_name, optic1, angles_rad, init_params, lower_bounds, 
        upper_bounds, step, disp)
    else	#Conduct a local search
        fit_1c_local(file_path, file_name, optic1, angles_rad, init_params, lower_bounds, 
        upper_bounds, disp)
    end
end

function fit_SAIM(file_path::String, file_name::String, optic1::SAIMOptics, optic2::SAIMOptics, 
    angles::AbstractArray, init_params::AbstractArray, lower_bounds::AbstractArray, 
    upper_bounds::AbstractArray; glb::Bool=false, step::Float64=40.0, disp::Bool=false)
    
    angles_rad = angles*2*pi/360  	
    
    if glb	#Run a global grid search when true
        fit_2c_global(file_path, file_name, optic1, optic2, angles_rad, init_params, lower_bounds, 
        upper_bounds, step, disp)
        println("2c global")
    else	#Conduct a local search
        fit_2c_local(file_path, file_name, optic1, optic2, angles_rad, init_params, lower_bounds, 
        upper_bounds, disp)
        println("2c local")
    end
end

function fit_SAIM(file_path::String, file_name::String, optic1::SAIMOptics, optic2::SAIMOptics, 
    optic3::SAIMOptics, angles::AbstractArray, init_params::AbstractArray, lower_bounds::AbstractArray, 
    upper_bounds::AbstractArray; glb::Bool=false, step::Float64=40.0, disp::Bool=false)
    
    angles_rad = angles*2*pi/360  	
    
    if glb	#Run a global grid search when true
        fit_3c_global(file_path, file_name, optic1, optic2, optic3, angles_rad, init_params, 
        lower_bounds, upper_bounds, step, disp)
    else	#Conduct a local search
        fit_3c_local(file_path, file_name, optic1, optic2, optic3, angles_rad, init_params, 
        lower_bounds, upper_bounds, disp)
    end
end



