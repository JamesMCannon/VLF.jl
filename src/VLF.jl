module VLF

    #Import needed modules
    using MAT
    using Dates
    using PyPlot

    #Export all function names
    export read_multiple_mat_files
    export read_data
    export start_time
    export merge_data
    export label_maker
    export plot_data
    export unwrap

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

function start_time(time_data::Tuple) #Get time series for each data set, push into array of arrays
	times = Any[]
	date_times = []
	for i in 1:length(time_data[3])
		push!(times, collect(Dates.value.(Time(time_data[4][i],time_data[5][i],time_data[6][i]))/1e9:Dates.value.(Time(time_data[4][i],time_data[5][i],time_data[6][i])+Dates.Second(length(time_data[7][i])-1))/1e9)/3600)
		push!(date_times, Dates.value(Date(time_data[1][i],time_data[2][i],time_data[3][i])))
	end
	return times,date_times
end

function merge_data(data_times::Tuple,time_data::Tuple,un_wrap::Bool) 
    #Collect each day of data into their own array, then push array of full day of data into final_data. 
    #Same with times for each day
    final_data = Any[]
    final_time = Any[]
    years = []
    months = []
    copy_data = copy(time_data[7])
    copy_start = unique(data_times[2]) #Unique Date-Time objects - covers the span
    duration = copy(data_times[2]) #Date-Time objects associated with each input file
    copy_time = copy(data_times[1])
    day_data = Any[]
    date_data = Any[]
    for i in 1:eachindex(copy_start) #cycle through unique days
        empty!(day_data) #clear current day data array
        empty!(date_data)#clear current day time array
        for j in eachindex(duration) #cycle through unique days... should this be 'duration' instead?
            #if j > length(copy_data)
            #    break
            #end
            if copy_start[i]==duration[j]
                append!(day_data, copy_data[j])
                append!(date_data, copy_time[j])
            end
        end
        if length(day_data) > 0
            if un_wrap == true
                copy_day_data = deepcopy(day_data)
                copy_date_data = deepcopy(date_data)
                push!(final_time, copy_date_data)
                push!(final_data, unwrap(copy_day_data))
                push!(years, year(copy_start[i]))
                push!(months, month(copy_start[i]))
                #push!(years, time_data[1][i])
                #push!(months, time_data[2][i])
            else
                copy_day_data = deepcopy(day_data)
                copy_date_data = deepcopy(date_data)
                push!(final_time, copy_date_data)
                push!(final_data, copy_day_data)
                push!(years, year(copy_start[i]))
                push!(months, month(copy_start[i]))
                #push!(years, time_data[1][i])
                #push!(months, time_data[2][i])
            end
        end
    end
    return final_data,final_time,years,months
end

function label_maker(merged_data::Tuple,time_data::Tuple)
    line_labels = []
    unique_days = []
    for n in 1:length(time_data[3])-1
        if time_data[3][n+1]!=time_data[3][n]
            append!(unique_days, time_data[3][n])
        end
    end
    if time_data[3][end]!=unique_days[end]
        append!(unique_days, time_data[3][end])
    end
    for i in 1:length(merged_data[3])
        push!(line_labels, Date(merged_data[3][i],merged_data[4][i],unique_days[i]))
    end
    return line_labels
end

function plot_data(axis_data1::Tuple, line_labels, xlabel1, plot_title)
    n = length(axis_data1[1])

    fig, ax = subplots()
    cmap = matplotlib.cm.get_cmap("jet")

    for i in 1:n
        line_color = cmap((i - 1) / (n - 1))
        ax.plot(axis_data1[2][i], axis_data1[1][i], color=line_color, 				label=line_labels[i], linewidth=0.5)
    end
    ax.grid(true)

    norm = matplotlib.colors.Normalize(vmin=0, vmax=n - 1)
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Date")
    xticks(0:2:24)
    ylabel(xlabel1)
    #ylim((-180,180)) #set y-limits
    xlabel("UTC Time (Hrs)")
    title(plot_title)

    # Set custom tick labels on the colorbar
    tick_positions = range(0, stop=n-1, length=n)
    cbar.set_ticks(tick_positions)
    cbar.set_ticklabels(line_labels)

    show()
end

function unwrap(phase::Vector{Any})
    unwrapped_phase = similar(phase)
    unwrapped_phase[1] = phase[1]

    for i = 2:length(phase)
        diff = phase[i] - phase[i - 1]
        if diff > 180
            phase[i] = phase[i] - 360 * ceil((diff - 180) / 360)
        elseif diff < -180
            phase[i] = phase[i] + 360 * ceil((-diff - 180) / 360)
        else
            phase[i] = phase[i]
        end
    end

    return phase
end

end
