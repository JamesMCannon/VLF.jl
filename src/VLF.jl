module VLF

    #Import needed modules
    using MAT
    #Export all function names
    export read_multiple_mat_files
    export read_data

    function read_multiple_mat_files(folder_path::AbstractString, 
        

        include_pattern::AbstractString)
        mat_files = filter(f -> isfile(joinpath(folder_path, f)), readdir(folder_path))
        data = []
    
        for file in mat_files
            if occursin(".mat", file)
                file_path = joinpath(folder_path, file)
                
                # Check if the file should be included based on the common part of the file name
                include_file = occursin(include_pattern, file)
                
                if include_file
                    mat_contents = MAT.matopen(file_path)
                    push!(data, mat_contents)
                end
            end
        end
    
            return data
    end

    function read_data(files::Vector) #Get all information needed form each data file, push into array of arrays
		start_year = Vector{Float64}()
		start_month = Vector{Float64}()
		start_day = Vector{Float64}()
		start_hour = Vector{Float64}()
		start_minute = Vector{Float64}()
		start_second = Vector{Float64}()
		data = []
		for i in 1:length(files)
            try
			    curr_day = read(files[i],"start_day")
                curr_min = read(files[i],"start_minute")
                curr_hour = read(files[i],"start_hour")
                curr_sec = read(files[i],"start_second")
                curr_month = read(files[i],"start_month")
                curr_year = read(files[i],"start_year")
                append!(start_day, curr_day)
			    append!(start_minute, curr_min)
			    append!(start_hour, curr_hour)
			    append!(start_second, curr_sec)
			    append!(start_month, curr_month)
			    append!(start_year, curr_year)
			    push!(data, read(files[i],"data"))
            catch x
                println("Unable to read file:", files[i])
            end
		end
		return start_year,start_month,start_day,start_hour,start_minute,start_second,data
end

end
