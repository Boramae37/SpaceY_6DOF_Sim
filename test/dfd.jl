using Pkg, Plots, GLMakie, LinearAlgebra, Statistics

# t 값의 범위를 설정합니다.
t = 0:0.1:10π

# x, y, z 값 계산
x = t .* cos.(t)
y = t .* sin.(t)
z = t

# 3차원 곡선을 그립니다.
sd = plot3d(x, y, z, label="Parametric Curve", xlabel="X", ylabel="Y", zlabel="Z", title="3D Parametric Curve", linewidth=2)
sd