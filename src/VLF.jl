module VLF

    #Import needed modules
    using MAT
    using Dates
    using PyPlot
    using Parameters

    #Export all function names; can use VLF.foo() for any function, exported or not
    export read_multiple_mat_files
    export read_data
    export start_time
    export merge_data
    export label_maker
    export plot_data
    export unwrap
    export calibrate_NB
    export combine_2ch
    export plot_day
    export plot_multi_site

    function read_multiple_mat_files(folder_path::AbstractString, include_pattern::AbstractString)
        mat_files = filter(f -> isfile(joinpath(folder_path, f)), readdir(folder_path))
        data = []
    
        for file in mat_files
            if occursin(".mat", file)
                file_path = joinpath(folder_path, file)
                
                # Check if the file should be included based on the common part of the file name
                include_file = occursin(include_pattern, file)
                
                if include_file
                    mat_contents = matopen(file_path)
                    push!(data, mat_contents)
                    close(file_path)
                end
            end
        end
    
        return data
    end

    function read_data(files::Vector,first_date::String="2000-01-01",last_date::String="2500-01-01") #Get all information needed form each data file, push into array of arrays
		start_year = Vector{Float64}()
		start_month = Vector{Float64}()
        start_day = Vector{Float64}()
        start_hour = Vector{Float64}()
        start_minute = Vector{Float64}()
        start_second = Vector{Float64}()
        data = Any[]
        Fs = Vector{Float64}()
        Fc = Vector{Float64}()
        adc_channel_number = Vector{Float64}()
        for i in eachindex(files)
            try 
                if Dates.value(Date(read(files[i],"start_year")[1],read(files[i],"start_month")[1],read(files[i],"start_day")[1])) >= Dates.value(Date(first_date)) && Dates.value(Date(read(files[i],"start_year")[1],read(files[i],"start_month")[1],read(files[i],"start_day")[1])) <= Dates.value(Date(last_date))

                    #Need the [1] behind each "read" function because "read" returns an array instead of an Int64
				
                    append!(start_day, read(files[i],"start_day"))
                    append!(start_minute, read(files[i],"start_minute"))
                    append!(start_hour, read(files[i],"start_hour"))
                    append!(start_second, read(files[i],"start_second"))
                    append!(start_month, read(files[i],"start_month"))
                    append!(start_year, read(files[i],"start_year"))
                    push!(data, read(files[i],"data"))
                    append!(Fs, read(files[i],"Fs"))
                    append!(Fc, read(files[i],"Fc"))
                    append!(adc_channel_number, read(files[i],"adc_channel_number"))

                    close(files[i])
                end
            catch x
                println("Unable to read file: ",files[i])
            end
        end
        return start_year,start_month,start_day,start_hour,start_minute,start_second,data,Fs,Fc,adc_channel_number
	end

    function calibrate_NB(raw_data::Tuple; cal_file::String="default",cal_num::Float64=-1.0)
        #=  
            raw_data: Narrowband data returned by the read_data() fcn. Expected to be a tuple with 10 dimensions 
            cal_file: the path to the relevant calibration file for the receiver who's data is contained in raw_data - if applicable
            cal_num: The cal_num coefficient for the receiver channel corresponding to the data in raw_data

            calibrate_NB calibrates the Narrowband data contained in raw_data using the frequency response contained in cal_file or specified in cal_num

            authors: James M Cannon
            date of last modification: 02/09/24
        =#

        cal_structure = deepcopy(raw_data) #copy the raw data structure to not lose timing information on return

        if cal_file=="default" && cal_num!=-1.0
            #use cal_num
            cal_factor = cal_num
        elseif cal_num==-1.0 && cal_file!="default"
            mat_contents = matopen(cal_file)
            adc_channel_number = raw_data[10][1] #identify the channel number from the data file which specifies ns/ew
            Fc = raw_data[9][1] #center frequency of the transmitter tracked
            if adc_channel_number == 0 #0 => NS, 1=> EW
                cal_curve = read(mat_contents,"CalibrationNumberNS")
            elseif adc_channel_number == 1
                 cal_curve = read(mat_contents,"CalibrationNumberEW")
            else
                error("Unexpected Channel Number: ", adc_channel_number)
            end
    
            freqs = abs.(cal_curve[:,1])
            response = abs.(cal_curve[:,2])
    
            Fc_idx = findmin(abs.(freqs.*1000 .-Fc))[2] #find the index of the entry closest to the transmitter frequency used in the data
    
            cal_factor = response[Fc_idx]
            #CalibrationNumber, contained in response[], is a coefficient to convert from DAQ counts to pT
            close(mat_contents)
        else
            error("Calibration information not provided")
        end 
        for i in 1:length(raw_data[1])
            calibrated_data = raw_data[7][i] .* cal_factor
    
            cal_structure[7][i] = calibrated_data
            #replace the uncalibrated data in cal_structure[7][] with the calibrated data  
        end
        close(cal_file)
        return cal_structure
    
    end

    function start_time(time_data::Tuple)
        times = Any[]
        date_times = []
        for i in 1:length(time_data[3])
            push!(times, collect(Dates.value.(Time(time_data[4][i],time_data[5][i],time_data[6][i]))/1e9:1/time_data[8][i]:Dates.value.(Time(time_data[4][i],time_data[5][i],time_data[6][i])+Dates.Millisecond(length(time_data[7][i])*1000/time_data[8][i]-1))/1e9)/3600)

            #for times, the increments are in milliseconds. We multiply the length of the data array by 1000/frequency to trick the dataset into thinking it is being sampled at 1Hz (This gives us the correct time duration for the data set). The time array has steps of 1/Frequency so there are enough "data points" every second.
			
            push!(date_times, Dates.value(Date(time_data[1][i],time_data[2][i],time_data[3][i]))) 
            #Same date_times variable, nothing changed
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
        for i in eachindex(copy_start) #cycle through unique days
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

    function combine_2ch(cal_data1::Tuple, cal_data2::Tuple)
        #=
            cal_data1: calibrated channel 1 data of the 2 channels to be combined
            cal_data2: calibrated channel 2 data of the 2 channels to be combined
        
            This function adds in quadrature 2 channels of calibrated data together
        
            authors: James M Cannon
            date of last modification: 02/09/24
        =#
        
           Combined_Data = deepcopy(cal_data1)
           #copy the timing data contained in the rest of the data tuple to preserve when the tuple is returned
        
           for i = 1:length(cal_data1[1])

                #check the lengths to be merged before merging
                #more robust solution desired for situations where data is off by only 1 second
                if length(cal_data1[1][i]) == length(cal_data2[1][i])
                    Combined_Data[1][i] = sqrt.(cal_data1[1][i].^2 .+ cal_data2[1][i].^2)
                    #add in quadrature (c = sqrt(a^2+b^2))
                end
           end
           return Combined_Data
        
        end

    function label_maker(merged_data::Tuple,time_data::Tuple,number)
        line_labels = []
        unique_days = []
            for n in 1:length(time_data[3])-1
                if time_data[3][n+1]!=time_data[3][n]
                    append!(unique_days, time_data[3][n])
                end
            end
            n=length(time_data[3])
            m = length(unique_days)
            if time_data[3][n]!=unique_days[m]
                append!(unique_days, time_data[3][n])
            end
            for i in 1:length(merged_data[3])
                push!(line_labels, Date(merged_data[3][i],merged_data[4][i],unique_days[i]))
            end
            final_line_labels = []
            k = round(length(line_labels)/(number-1))
            l = round((Dates.value(last(line_labels))-Dates.value(line_labels[1]))/(number-1))
            j = collect(1:number)
            h = 0
            for i in eachindex(line_labels)
                if i == 1 || i == length(line_labels)
                    push!(final_line_labels, line_labels[i])
                    h+=1
                elseif (i == j[h]*k && i+5 < length(line_labels))
                    push!(final_line_labels, Date(Dates.UTD(Dates.value(line_labels[1])+j[h]*l)))
                    h+=1
                else
                    push!(final_line_labels, "")
                end
    
            end
        return final_line_labels,line_labels
    end

    function plot_data(axis_data1::Tuple, line_labels::Tuple, xlabel1, plot_title, xlim_fcn=[0 24], xticks_fcn=(0:4:24))
	    n = length(axis_data1[1])
        base = Dates.value(line_labels[2][1]) #full value of first day
	    endpoint = Dates.value(last(line_labels[2])) #full value of last day
	
	    fig, ax = subplots()
	    cmap = matplotlib.cm.get_cmap("jet")
		  tick_positions = Vector{Float64}()
		  final_line_labels=[]
	
	    for i in 1:n
			#Line_color can be a number from 0 to 1, and divides the current days "value" by the last days "value". So the first day will have a cmap index of 0, and the last day will have a cmap index of 1
	        line_color = cmap((Dates.value(line_labels[2][i])-base) / (endpoint-base))

			#Same principal as the cmap index, but the ending tick of the colorbar corresponds to how many data points are represented by the colorbar. So, multiply by (n-1) so the index goes from 0 to (# of days plotted) instead of 0 to 1.
			    if line_labels[1][i] != ""
                    append!(tick_positions, ((Dates.value(line_labels[1][i])-base) / (endpoint-base))*(n-1))
                    push!(final_line_labels, line_labels[1][i])
                    ax.plot(axis_data1[2][i], axis_data1[1][i], color=line_color, 				label=line_labels[1][i], linewidth=0.5)
			    else
                    ax.plot(axis_data1[2][i],axis_data1[1][i], color=line_color, label="", linewidth=0.5)
			    end
	    end
	    ax.grid(true)
	
	    norm = matplotlib.colors.Normalize(vmin=0, vmax=n - 1)
	    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
	    sm.set_array([])
	    cbar = fig.colorbar(sm, ax=ax)
	    cbar.set_label("Date")
        xticks(xticks_fcn)
        #xticks(0:1:24)
        ylabel(xlabel1)
        xlabel("Time (Hrs)")
        title(plot_title)
        xlim(xlim_fcn)
	    #set custom tick labels on the colorbar
	    cbar.set_ticks(tick_positions)
	    cbar.set_ticklabels(final_line_labels)
	
	    show()
	end

    @with_kw struct Plot_Options
        #=
            xlabel: x-axis label 
            ylabel: y-axis label
            title: plot title
            xlim: x-axis limits
            xticks: x-axis tick marker locations
            primary_color: primary color used in the plot
            grid: which gridlines to show (minor, major, both)

            This structure contains common plotting options to be reused by various plot functions

            authors: James M Cannon
            date of last modification: 12/14/23
        =#

        xlabel = "Time (Hrs)"
        ylabel
        title
        xlim = [0, 24]
        xticks = (0:4:24)
        ylim = []
        yticks
        primary_color = "#3CA0FA"
        grid = "both"
    end

    function plot_day(axis_data_fcn::Vector; line_label="",opts::Plot_Options)
        #=
            axis_data_fcn: vector of [1][:] data and [2][:] timestamps
            line_label: desired label for line (default of "")
            opts: Plot_Options structure containing various options for modifying a plot

            This function plots a single day of narrowband data
        
            authors: James M Cannon
            date of last modification: 12/14/23
        =#

        plt.plot(axis_data_fcn[2][:],axis_data_fcn[1][:],color=opts.primary_color, label=line_label,linewidth=0.5)
        plt.grid(which = opts.grid)
        plt.xticks(opts.xticks)
        plt.ylabel(opts.ylabel)
        plt.xlabel(opts.xlabel)
        plt.title(opts.title)
        plt.xlim(opts.xlim)
        plt.show()
    end

    function plot_multi_site(axis_data_fcn::Vector; line_labels="",opts::Plot_Options)
        #=
            axis_data_fcn: vector of [n][1][:] data and [n][2][:] timestamps where n is the number of sites
            line_labels: desired labels for each of n lines (default of "")
            opts: Plot_Options structure containing various options for modifying a plot

            This function plots a single day of narrowband data across multiple sites
        
            authors: James M Cannon
            date of last modification: 12/15/23
        =#
        count=1
        for i in axis_data_fcn
            plt.plot(i[2][:],i[1][:],linewidth=0.5,color=opts.primary_color[count])
            count=count+1
        end
        plt.grid(which = opts.grid)
        plt.xticks(opts.xticks)
        plt.ylabel(opts.ylabel)
        plt.xlabel(opts.xlabel)
        plt.title(opts.title)
        plt.xlim(opts.xlim)
        if !isempty(opts.ylim)
            plt.ylim(opts.ylim)
            plt.yticks(opts.yticks)
        end
        if !isempty(line_labels)
            plt.legend(line_labels)
        end
        plt.show()
    end

    function unwrap(phase::Vector{Any})
        unwrapped_phase = copy(phase)
        unwrapped_phase[1] = phase[1]

        for i in 2:length(phase) #eachindex() doesn't work here and I don't know why
            diff = phase[i] - phase[i - 1]
            if diff > 180 #if jump goes from neg to pos
                phase[i] = phase[i] - 360 * ceil((diff - 180) / 360)
                #phase[i] = -180 - (180-phase[i])
            elseif diff < -180 #if jump goes from pos to neg
                phase[i] = phase[i] + 360 * ceil((-diff - 180) / 360)
                #phase[i] = 180 - (-180-phase[i])
            else
                phase[i] = phase[i]
            end
        end

        return phase
    end
    
end