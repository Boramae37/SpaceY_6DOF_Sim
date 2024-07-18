# ============================ SPACE Y ====================================
#  6 Degree of Freedom Flight Simulation Program (ver 1.0)
#  Release Date : 2024.06.08
#  Script by : Sehwan An
#  Last Edited: 2024.06.12, Hanjun Oh
#  Function : 6 Degree of Freedom Flight Simulation Program
#             (Flat Earth, Atmosphere, No sloshing and bending)
# =========================================================================

@time begin
# Import
using CSV, XLSX, DataFrames, Interpolations, Plots, LinearAlgebra, Dates
include("set/func.jl")
include("set/rot.jl")

# Setup
ThrustMode = "interp"; # interp : Thrust CSV file interpolation / calc : Thrust equation formula
WindMode = "off"; # on : with wind / off : without wind (if you want to use wind data, you must change "func.jl")
LaunchAngle = 90;


# ================================================================================================ #
# 아웃풋을 변경하고 싶지 않은 이상 건드리지 말 것.

# Simulation
(idxf, idxb, t, tb, x, ang, Ft, Fa, M) = simulation(ThrustMode, WindMode, LaunchAngle);


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
fFa = plot(t,Fa[1,:], title="Aerodynamic Force", xlabel="Time[s]",ylabel="Force[N]", legend = false, dpi=300); 
fM = plot(tb,M[1,1:idxb], title="Mass", xlabel="Time[s]",ylabel="Mass[kg]", legend = false, dpi=300); 


# Save
nt = "res/tms_2024_06_27_thrust/"    # 저장폴더의 경로 설정
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
savefig(fFa,nt*"/Fa.png");
savefig(fM,nt*"/m.png");
data = DataFrames.DataFrame(time=t, E=vec(x[1,:]), N=vec(x[2,:]), U=vec(x[3,:])
                            ,u=vec(x[4,:]), v=vec(x[5,:]), w=vec(x[6,:])
                            ,phi=vec(ang[1,:]), theta=vec(ang[2,:]), psi=vec(ang[3,:])
                            ,Thrust=vec(Ft[1,:]), Aerodynamic=vec(Fa[1,:]), Mass=vec(M));
XLSX.writetable(nt*"/simdata.xlsx","Data"=>data,overwrite=true);

end