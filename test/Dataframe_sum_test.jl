using CSV, XLSX, DataFrames, Interpolations, Plots, LinearAlgebra, Dates

stage_seperate_time = 2.21; # 분리 시점 [s], 0.01s 최소 단위로 설정
second_fire_time = 1.0;  # 분리 시점으로 부터 2단 연소 시작 시간 [s], 0.01s 최소 단위로 설정
first_file_name = "db/data.csv" # 1st stage Thrust CSV file name
second_file_name = "db/thrust.csv" # 2nd stage Thrust CSV file name


thrust = CSV.read("db/data.csv", DataFrame);

push!(thrust, [stage_seperate_time + second_fire_time-0.01 ,0.0])
thrust2 = CSV.read(second_file_name, DataFrame);
thrust2[:,1] .= thrust2[:,1] .+ (stage_seperate_time + second_fire_time);
thrust = vcat(thrust, thrust2);

print(thrust);


