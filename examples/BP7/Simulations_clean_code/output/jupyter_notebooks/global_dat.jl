using CSV
using DataFrames
using Statistics

output_dir = "../202503041618"
file_lists = readdir(output_dir)
station_files = filter(f -> occursin("fltst_strk", f), file_lists)

output_M =  station_files
Mvw_index = [1,2,4,6,7,8,10,11,12]
output_Mvw = station_files[Mvw_index]


file_name = file_lists[1]
file_path = string(output_dir,"/",file_name)
content = CSV.File(file_path; skipto=19, header=0) |> DataFrame

# slip_rate_2_log = content[:,4]
# slip_rate_2 = 10 .^ slip_rate_log 
# slip_rate_3_log = content[:,5]
# slip_rate_3 = 10 .^ slip_rate_3_log
# slip_rate = hypot.(slip_rate_2, slip_rate_3)


## Starting to obtain max_slip_rate
df = DataFrame()
for file_name in station_files
    file_path = string(output_dir,"/",file_name)
    content = CSV.File(file_path; skipto=19, header=0) |> DataFrame
    
    slip_rate_2_log = content[:,4]
    slip_rate_2 = 10 .^ slip_rate_log 
    slip_rate_3_log = content[:,5]
    slip_rate_3 = 10 .^ slip_rate_3_log
    slip_rate = hypot.(slip_rate_2, slip_rate_3) 

    header = first(split(file_name,"."))
    df[!, header] = slip_rate
end



max_output_M = maximum(eachcol(df))
mean_output_M =  mean(eachcol(df))
mean_output_Mvw = mean(eachcol(df[:,Mvw_index]))


maximum_slip_rate = log.(10, max_output_M)

moment_rate = mean_output_M .* π * 300^2 * 32.04 * 10^9


moment_rate_vw = mean_output_M .* π * 200^2 * 32.04 * 10^9