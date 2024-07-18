# RK4 Integral
function solver(f,h,x,Ft, Fa, I, V, g, m)
    k1 = h * f(x, Ft, Fa, I, V, g, m);
    k2 = h * f(x+1/2*k1, Ft, Fa, I, V, g, m); 
    k3 = h * f(x+1/2*k2, Ft, Fa, I, V, g, m);
    k4 = h * f(x+k3, Ft, Fa, I, V, g, m);
    xn = x + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    return xn
end

# Solve the Equation
function ode(x, Ft, Fa, I, V, g, m) 
    # variables
    u = x[4]; 
    v = x[5];
    w = x[6];
    q1 = x[7];
    q2 = x[8];
    q3 = x[9];
    q4 = x[10];
    p = x[11];
    q = x[12];
    r = x[13];

    # dot terms
    pos_dot = b2n(x[7:10]) * V; # ENU 좌표계로 변환
    vel_dot = -[0 -r q; r 0 -p; -q p 0] * [u; v; w] + n2b(x[7:10])*[0; 0; g] + (1/m)*Ft + (1/m)*Fa;
    ang_dot = 0.5 * [0 -p -q -r;
                     p  0  r -q;
                     q -r  0  p;
                     r  q -p  0]*[q1; q2; q3; q4];
    agr_dot = I\(-[0 -r q;r 0 -p;-q p 0]*I*[p; q; r]);

    # dot matrix
    dx = [pos_dot;
          vel_dot;
          ang_dot;
          agr_dot];

    return dx
end


# Simulation
function simulation(ThrustMode, WindMode, LaunchAngle, Phase_Change, second_fire_time,
                    g, rho, tb1, tb2,
                    m0, Isp0, Ixx0, Iyy0, Izz0, diam0, Cdb0, Body_Area0,
                    Cdpu, Cdpv, Cdpw, Apu, Aps,
                    m2, Isp2, Ixx2, Iyy2, Izz2, diam2, Cdb2, Body_Area2,
                    first_file_name, second_file_name)

    r2d = 180/pi;
    d2r = pi/180;
    dt = 0.01;                              # Sampling Time [s]
    tb = tb1;                               # tb = Burn Time [s] 2단: tb1 + tb2 + second_fire_time+0.1, 1stage모드: tb11
    Tb = 0:dt:tb;                           # Burn Time [s]  vector
    idxb = trunc(Int64,tb/dt+1);            # Burn Time [s] index 개수
    tf = 200;                                # Simulation Stop Time [s], 낙하시간 보다 길게
    Tt = 0:dt:(tb+tf);                      # Total Time [s] vector 
    idxf = trunc(Int64,(tb+tf)/dt+1);       # Time [s] index 개수

    Mode = 1;         # Do not change.

    m = m0;           # Total Mass [kg]
    Isp = Isp0;       # Specific impulse [sec]
    Ixx = Ixx0;       # Moment of Inertia [kgm^2]
    Iyy = Iyy0;       # Moment of Inertia [kgm^2]
    Izz = Izz0;       # Moment of Inertia [kgm^2]
    diam = diam0;     # Body diameter [m]
    Cdb = Cdb0;       # Body drag coefficient
    Body_Area = Body_Area0; # Body Area [m^2] 옆에서 본 면적


    # State variables
    x = zeros(13,trunc(Int64,idxf));
    Vw = zeros(3,trunc(Int64,idxf)); # 풍속 벡터
    Vr = zeros(3,trunc(Int64,idxf)); # 로켓 속도 벡터
    Vabs = zeros(1,trunc(Int64,idxf)); # 로켓 속도 크기
    ang = zeros(3,trunc(Int64,idxf));
    Ft = zeros(3,trunc(Int64,idxf));
    Fa = zeros(3,trunc(Int64,idxf));
    M = zeros(1,trunc(Int64,idxf));
    I = zeros(3,3,trunc(Int64,idxf));

    

    # Thrust Setup
    if ThrustMode == "interp";
        thrust = CSV.read(first_file_name, DataFrame);

        if Phase_Change == "manual"
            push!(thrust, [stage_seperate_time + second_fire_time-0.01 ,0.0])
            thrust2 = CSV.read(second_file_name, DataFrame);
            thrust2[:,1] .= thrust2[:,1] .+ (stage_seperate_time + second_fire_time);
            thrust = vcat(thrust, thrust2);
        end



        thrust = linear_interpolation(thrust[:,1],thrust[:,2],extrapolation_bc=Periodic()); # interpolation 함수  thrust의 1열을 x, 2열을 y로

    elseif ThrustMode == "calc";
        # 엔진 타들어가는 속도 계산해서 로켓방정식 수정.
        mf = 0.400;     # 연료 질량 [kg]
        mfr = mf/tb;    # 연료 감소량 [kg/s]
    end

    

    # Initial conditions
    x[:,1] = [0; 0; 0; 0; 0; 0;  0;  0;  0;  0; 0; 0; 0];
    #        [x  y  z  u  v  w  q1  q2  q3  q4  p  q  r]
    ang[:,1] = [0; LaunchAngle*d2r; 0];
    x[7:10,1] = e2q(ang[:,1]);
    I[1,1,:] .= Ixx;
    I[2,2,:] .= Iyy;
    I[3,3,:] .= Izz;


    # Controller (TBD)

    # Simulation
    for i in 1:1:idxf

        # phase2 단 분리 시점에, 질량, 관성모멘트, 지름, 항력계수 변경
        if Phase_Change == "manual"
            if i == stage_seperate_time*100
                #print(i,"\n");
                m = m2;
                Isp = Isp2;
                I[1,1,:] .= Ixx2;
                I[2,2,:] .= Iyy2;
                I[3,3,:] .= Izz2;
                diam = diam2;
                Cdb = Cdb2;
                Body_Area = Body_Area2;
            end
        
        elseif Phase_Change == "1stage"
            if i == stage_seperate_time*100
                m = m-m2;
                Body_Area = Body_Area0 - Body_Area2;
            end
        end


        M[1,i] = m;
        
        # Wind 
        # if WindMode == "on"
        #     Vw[:,i] = [0; -1; 0];  # 풍속 벡터
        #     if Mode == 1
        #         print("x",x[4:6,i],"\n");
        #         Vr[1,i] = x[4,i] + Vw[1,i] + 2*m*Vw[1,i]/(Cdb*(i/100)*Vw[1,i]-2*m);
        #         Vr[2,i] = x[5,i] + Vw[2,i] + 2*m*Vw[2,i]/(Cdb*(i/100)*Vw[2,i]-2*m);
        #         Vr[3,i] = x[6,i] + Vw[3,i] + 2*m*Vw[3,i]/(Cdb*(i/100)*Vw[3,i]-2*m);
                
        #     elseif Mode == 2
        #         print("x",x[4:6,i],"\n");
        #         Vr[1,i] = x[4,i] + Vw[1,i] + 2*m*Vw[1,i]/(Cdpu*(i/100)*Vw[1,i]-2*m);
        #         Vr[2,i] = x[5,i] + Vw[2,i] + 2*m*Vw[2,i]/(Cdpv*(i/100)*Vw[2,i]-2*m);
        #         Vr[3,i] = x[6,i] + Vw[3,i] + 2*m*Vw[3,i]/(Cdpw*(i/100)*Vw[3,i]-2*m);
        #     end
        # elseif WindMode == "off"
        #     Vw[:,i] = [0; 0; 0];
        # end
        #Vr[:,i] = x[4:6,i]-Vw[:,i]; # 속도 벡터(로켓 속도 - 풍속)
        #Vr[:,i] = Vr[:,i] + x[4:6,i];
        #print(i,Vr[:,i],"\n");
        #print("wind",Vw[:,i],"\n");
        #print("x",x[4:6,i],"\n");
        
        
        
        ang[:,i] = q2e(x[7:10,i]).*r2d;


        # Thrust Force
        if ThrustMode == "calc"
            if i < idxb
                Ft[1,i] = mfr*Isp*g;
            else
                Ft[1,i] = 0;
            end
        elseif ThrustMode == "interp"
            if i < idxb
                Ft[1,i] = thrust(Tt[i]);
            else
                Ft[1,i] = 0;
            end
        end

        # 로켓 속도
        Vr[:,i] = x[4:6,i];
        Vabs[i] = norm(Vr[:,i]);    # 속도 크기
        
        # Aero Force
       # 로켓 속도와 반대 방향으로 항력 발생 
        if WindMode == "off"
            if Mode == 1
                Fa[1,i] = -sign(Vr[1,i])*0.5 * Cdb * (pi*(diam^2)/4) * rho * (Vr[1,i].^2);
                Fa[2,i] = -sign(Vr[2,i])*0.5 * Cdb * Body_Area * rho * (Vr[2,i].^2);
                Fa[3,i] = -sign(Vr[3,i])*0.5 * Cdb * Body_Area * rho * (Vr[3,i].^2);
            elseif Mode == 2
                Fa[1,i] = -sign(Vr[1,i])*0.5 * Cdpu * Apu * rho * (Vr[1,i].^2);
                Fa[2,i] = -sign(Vr[2,i])*0.5 * Cdpv * Aps * rho * (Vr[2,i].^2);
                Fa[3,i] = -sign(Vr[3,i])*0.5 * Cdpw * Aps * rho * (Vr[3,i].^2);
            end
        elseif WindMode == "on"
            Vw[:,i] = [0; 0; 0]; # - 역풍/ [u, v, w]
            if Mode == 1
                Fa[1,i] = -sign(Vr[1,i]-Vw[1,i])*0.5 * Cdb * (pi*(diam^2)/4) * rho * ((Vr[1,i]-Vw[1,i]).^2);
                Fa[2,i] = -sign(Vr[2,i]-Vw[2,i])*0.5 * Cdb * Body_Area * rho * ((Vr[2,i]-Vw[2,i]).^2);
                Fa[3,i] = -sign(Vr[3,i]-Vw[3,i])*0.5 * Cdb * Body_Area * rho * ((Vr[3,i]-Vw[3,i]).^2);
                
            elseif Mode == 2
                Fa[1,i] = -sign(Vr[1,i]-Vw[1,i])*0.5 * Cdpu * Apu * rho * ((Vr[1,i]-Vw[1,i]).^2);
                Fa[2,i] = -sign(Vr[2,i]-Vw[2,i])*0.5 * Cdpv * Aps * rho * ((Vr[2,i]-Vw[2,i]).^2);
                Fa[3,i] = -sign(Vr[3,i]-Vw[3,i])*0.5 * Cdpw * Aps * rho * ((Vr[3,i]-Vw[3,i]).^2);
            end
        end

        
        # if i == idxf
        #     print(i);
        #     break;
        # end

        # Solve
        x[:,i+1] = solver(ode,dt,x[:,i],Ft[:,i],Fa[:,i],I[:,:,i],Vr[:,i],g,m);
        x[3,i] = -x[3,i];
        m = m - (Ft[1,i]/(Isp*g))*dt;

        # Parachute condition
        if i > idxb
            if abs(x[3,i+1])-abs(x[3,i]) < 0
                Mode = 2;
            end
            #if x[4,i]  < -40.0 #최대 고도 2초후 숫자 조정 필요
            #    Mode = 2;
            #end
        end
        
        # Ground condition
        if x[3,i] < 0
            x[3,i] = 0;
        end

        # Stopping Criteria
        if i > idxb
            if x[3,i] < 0.1
                idxf = i;
                print("최대 속력 ", maximum(Vabs[:]),"m\n");
                print("종단 속력: ", Vabs[i],"m/s\n");
                break;
            end
        end
        
        
            
    end
    
    return idxf, idxb, Tt[1:idxf], Tb, x[:,1:idxf], ang[:,1:idxf], Ft[:,1:idxf], Fa[:,1:idxf], M[:,1:idxf];
end