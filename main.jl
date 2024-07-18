# ============================ SPACE Y ====================================
#  6 Degree of Freedom Flight Simulation Program (ver 1.0)
#  Release Date : 2024.06.08
#  Script by : Sehwan An
#  Last Edited: 2024.06.12, Hanjun Oh
#  Function : 6 Degree of Freedom Flight Simulation Program
#             (Flat Earth, Atmosphere, No sloshing and bending)
# =========================================================================

# =========================== SPACE Y ====================================
#  6 Degree of Freedom Flight Simulation Program (ver 2.0)
#  Release Date : 2024.06.29
#  Updated by : Hyunsung Baek
#  Last Edited: 2024.06.29 Hyunsung Baek
#  Function : 6 Degree of Freedom Flight Simulation Program
#             (Flat Earth, Atmosphere, No sloshing and bending)
# =========================================================================

@time begin
# Import
using CSV, XLSX, DataFrames, Interpolations, Plots, LinearAlgebra, Dates, FilePathsBase
include("set/func.jl")
include("set/rot.jl")

# Setup
ThrustMode = "interp";  # interp : Thrust CSV file interpolation / calc : Thrust equation formula
WindMode = "on";       # on : with wind / off : without wind (if you want to use wind data, you must change "func.jl")
LaunchAngle = 85;
Phase_Change = "1stage"; # auto : 자동으로 계산, manual : 수동 단 분리 시점과 연소 시점 지정, none : 단 분리 없음(1단 로켓 분석)
                        # auto의 경우 z축 속도 0이 되는 시점을 분리 시점으로, 연소 시점은 지정(미구현)
                        # manual의 경우 stage_seperate_time 변수 선언 필요
                        # 1stage: 1단 로켓 추락까지 추적
                        # none
# ================================================================================================ #


# Parameters. 시뮬레이션 전에 맞춰서 수정 필요

tb1 = 2.18;                    # 1st stage Burn Time [s]
tb2 = 4.233;                    # 2nd stage Burn Time [s] 4.233 

stage_seperate_time = tb1+0.50;   # separate point [s], 0.01s 단위로 설정
second_fire_time = 2.8;          # 2nd stage ignition time from separation point [s], 0.01s 단위로 설정
#2.8
#const value
g    = 9.81;             # Gravity accel. [m/sec^2]
rho  = 1.229;            # Air density [kg/m^3]

# 1st+2nd stage info
m0   = 7.194;           # Mass [kg] 1st + 2nd stage
Isp0  = 104.9;               # Specific impulse [sec], 1st stage
Ixx0  = 6.091e-3;         # Moment of Inertia [kgm^2]
Iyy0  = 4.388e-2;         # Moment of Inertia [kgm^2]
Izz0  = 4.388e-2;         # Moment of Inertia [kgm^2]
diam0 = 0.104;            # Body diameter [m]
Cdb0  = 0.338;             # Body drag coefficient
Body_Area0 = 0.1879      # Body Area [m^2]  옆에서 본 면적

# Parachute info
Cdpu = 0.69;             # Parachute drag coefficient u방향
Cdpv = 0.69;             # Parachute drag coefficient v방향
Cdpw = 0.69;             # Parachute drag coefficient w방향 (Cdpv랑 같음)
Apu   = 0.518;            # Parachute Area [m^2]   위에서 본 면적  1단: 0.518 2단: 1.44
Aps   = 0.324;            # Parachute Area [m^2]   옆에서 본 면적   1단: 0.324 2단: 0.9


# 2nd stage info
m2   = 4.656;             # Mass [kg] 2nd stage
Isp2  = 101.54;               # Specific impulse [sec]
Ixx2  = 6.091e-3;         # Moment of Inertia [kgm^2]
Iyy2  = 4.388e-2;         # Moment of Inertia [kgm^2]
Izz2  = 4.388e-2;         # Moment of Inertia [kgm^2]
diam2 = 0.104;            # Body diameter [m]
Cdb2  = 0.338;             # Body drag coefficient
Body_Area2 = 0.12151;     # Body Area [m^2]  옆에서 본 면적


# 불러올 파일/ ThrustMode가 interp일 때만 사용
first_file_name = "db/tms_1s.csv"; # 1st stage Thrust CSV file name
second_file_name = "db/tms_2s.csv"; # 2nd stage Thrust CSV file name

# Save folder name
nt = "res/tms_angle_85_0715_4_wind_1stage/";




# ================================================================================================ #
# 아웃풋을 변경하고 싶지 않은 이상 건드리지 말 것.


# Simulation
(idxf, idxb, t, tb, x, ang, Ft, Fa, M) = simulation(ThrustMode, WindMode, LaunchAngle, Phase_Change, second_fire_time,
                                                    g, rho, tb1, tb2,
                                                    m0, Isp0, Ixx0, Iyy0, Izz0, diam0, Cdb0, Body_Area0,
                                                    Cdpu, Cdpv, Cdpw,Apu, Aps,
                                                    m2, Isp2, Ixx2, Iyy2, Izz2, diam2, Cdb2, Body_Area2,
                                                    first_file_name, second_file_name);



# Plot
fx = plot(t,x[1,:], title="Pos(ENU) E", xlabel="Time[s]",ylabel="x[m]", legend = false, dpi=300); 
fy = plot(t,x[2,:], title="Pos(ENU) N", xlabel="Time[s]",ylabel="y[m]", legend = false, dpi=300); 
fz = plot(t,x[3,:], title="Pos(ENU) U", xlabel="Time[s]",ylabel="z[m]", legend = false, dpi=300); 
fu = plot(t,x[4,:], title="Vel(Body) u", xlabel="Time[s]",ylabel="u[m/s]", legend = false, dpi=300); 
fv = plot(t,x[5,:], title="Vel(Body) v", xlabel="Time[s]",ylabel="v[m/s]", legend = false, dpi=300); 
fw = plot(t,x[6,:], title="Vel(Body) w", xlabel="Time[s]",ylabel="w[m/s]", legend = false, dpi=300);
fph = plot(t,ang[1,:], title="Angle(Euler) Phi", xlabel="Time[s]",ylabel="Phi[deg]", legend = false, dpi=300); 
fth = plot(t,ang[2,:], title="Angle(Euler) Theta", xlabel="Time[s]",ylabel="Theta[deg]", legend = false, dpi=300); 
fps = plot(t,ang[3,:], title="Angle(Euler) Psi", xlabel="Time[s]",ylabel="Psi[deg]", legend = false, dpi=300); 
fFt = plot(tb,Ft[1,1:idxb], title="Thrust Force", xlabel="Time[s]",ylabel="Force[N]", legend = false, dpi=300); 
fFau = plot(t,Fa[1,:], title="Aerodynamic Force", xlabel="Time[s]",ylabel="Force[N]", legend = false, dpi=300); 
fFav = plot(t,Fa[2,:], title="Aerodynamic Force", xlabel="Time[s]",ylabel="Force[N]", legend = false, dpi=300);
fFaw = plot(t,Fa[3,:], title="Aerodynamic Force", xlabel="Time[s]",ylabel="Force[N]", legend = false, dpi=300);
fM = plot(tb,M[1,1:idxb], title="Mass", xlabel="Time[s]",ylabel="Mass[kg]", legend = false, dpi=300); 
pos = plot3d(x[1,:],x[2,:],x[3,:], title="Pos(ENU)", xlabel="x[m]",ylabel="y[m]",zlabel="z[m]", legend = false, dpi=300);

# @gif for i in 1:1500
#     step!(attractor)
#     push!(pos, attractor.x, attractor.y, attractor.z)
# end every 10
# Save

# 저장할 폴더 없으면 생성

if !isdir(nt)
    mkpath(nt);
end

# 이미지 저장
savefig(fx,nt*"/x.png");
savefig(fy,nt*"/y.png");
savefig(fz,nt*"/z.png");
savefig(fu,nt*"/u.png");
savefig(fv,nt*"/v.png");
savefig(fw,nt*"/w.png");
savefig(fph,nt*"/phi.png");
savefig(fth,nt*"/theta.png");
savefig(fps,nt*"/psi.png");
savefig(fFt,nt*"/Ft.png");
savefig(fFau,nt*"/Fa_u.png");
savefig(fFav,nt*"/Fa_v.png");
savefig(fFaw,nt*"/Fa_w.png");
savefig(fM,nt*"/m.png");
savefig(pos,nt*"/pos.png");
data = DataFrames.DataFrame(time=t, E=vec(x[1,:]), N=vec(x[2,:]), U=vec(x[3,:])
                            ,u=vec(x[4,:]), v=vec(x[5,:]), w=vec(x[6,:])
                            ,phi=vec(ang[1,:]), theta=vec(ang[2,:]), psi=vec(ang[3,:])
                            ,Thrust=vec(Ft[1,:]), Aerodynamic_u=vec(Fa[1,:]), Aerodynamic_v=vec(Fa[2,:]), Aerodynamic_w=vec(Fa[3,:]), Mass=vec(M));
XLSX.writetable(nt*"/simdata.xlsx","Data"=>data,overwrite=true);

print("로켓 떨어진 거리: ", (x[1, end].^2 + x[2, end].^2)^(0.5), "m\n");
print("로켓 낙하 좌표: ", x[1, end], "m, ", x[2, end], "m\n");
print("로켓 최고점", maximum(x[3,:]), "m\n");
print("비행시간", t[end], "s\n");
end