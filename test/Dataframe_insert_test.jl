using CSV, XLSX, DataFrames, Interpolations, Plots, LinearAlgebra, Dates



thrust = CSV.read("db/thrust.csv", DataFrame);

println(thrust);


stage_seperate_time = 2.0
second_fire_time = 1.0

push!(thrust, [0.0,0.0])    
push!(thrust, [0.0,0.0])

change_index = 0
for i in reverse(eachindex(thrust[:, 1]))
    
    if thrust[i-2, 1] >= stage_seperate_time
        thrust[i,:] = thrust[i-2, :]
    else
        print(i-1,"\n", i,"\n");
        change_index = i-1
        break;
    end

end

thrust[change_index, 1] = stage_seperate_time
thrust[change_index, 2] = stage_seperate_time+second_fire_time
thrust[change_index::end,1] += second_fire_time

println(thrust);