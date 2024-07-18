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
    pos_dot = b2n(x[7:10]) * V;
    vel_dot = -[0 -r q;r 0 -p;-q p 0] * [u; v; w] + n2b(x[7:10])*[0; 0; g] + (1/m)*Ft + (1/m)*Fa;
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
function simulation(ThrustMode,WindMode,LaunchAngle)
    r2d = 180/pi;   #degree to radian  
    d2r = pi/180;   #radian to degree   
    dt = 0.01;              # Sampling Time [s]
    tb = 2.21;                 # Burn Time [s]
    Tb = 0:dt:tb;   
    idxb = trunc(Int64,tb/dt+1);
    tf = 60;                # Simulation Stop Time [s]
    Tt = 0:dt:(tb+tf);
    idxf = trunc(Int64,(tb+tf)/dt+1); 

    Mode = 1;   # Do not change.

    # Parameters
    m    = 4.4932;           # Mass [kg]
    Isp  = 89;               # Specific impulse [sec]
    Ixx  = 6.091e-3;         # Moment of Inertia [kgm^2]
    Iyy  = 4.388e-2;         # Moment of Inertia [kgm^2]
    Izz  = 4.388e-2;         # Moment of Inertia [kgm^2]
    g    = 9.81;             # Gravity accel. [m/sec^2]
    rho  = 1.229;            # Air density [kg/m^3]
    diam = 0.104;            # Body diameter [m]
    Cdb  = 0.39;             # Body drag coefficient
    Cdp  = 0.69;             # Parachute drag coefficient
    Ap   = 1.063;            # Parachute Area [m^2]

    # State variables
    x = zeros(13,trunc(Int64,idxf)); # [x  y  z  u  v  w  q1  q2  q3  q4  p  q  r]
    Vw = zeros(3,trunc(Int64,idxf));  # Wind velocity
    Vr = zeros(3,trunc(Int64,idxf)); # Relative velocity
    Vabs = zeros(1,trunc(Int64,idxf)); # Absolute velocity 
    ang = zeros(3,trunc(Int64,idxf)); # [roll pitch yaw] angle
    Ft = zeros(3,trunc(Int64,idxf));    # Thrust Force
    Fa = zeros(3,trunc(Int64,idxf));   # Aero Force
    M = zeros(1,trunc(Int64,idxf));   # Mass
    I = zeros(3,3,trunc(Int64,idxf)); # Moment of Inertia

    # Thrust Setup
    if ThrustMode == "interp";
        thrust = CSV.read("db/data.csv",DataFrame);
        thrust = linear_interpolation(thrust[:,1],thrust[:,2],extrapolation_bc=Periodic());

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
        M[1,i] = m;

        # Wind
        if WindMode == "on"
            Vw[:,i] = [0; 3; 3];
        elseif WindMode == "off"
            Vw[:,i] = [0; 0; 0];
        end
        Vr[:,i] = x[4:6,i]-Vw[:,i];
        Vabs[i] = norm(Vr[:,i]);
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

        # Aero Force
        Fa[1,i] = -0.5 * Cdb * (pi*diam^2/4) * rho * Vabs[i].^2;
        if Mode == 2
            Fa[1,i] = 0.5 * Cdp * Ap * rho * Vabs[i].^2;
        end

        # Solve
        x[:,i+1] = solver(ode,dt,x[:,i],Ft[:,i],Fa[:,i],I[:,:,i],Vr[:,i],g,m);
        x[3,i] = -x[3,i];
        m = m - (Ft[1,i]/(Isp*g))*dt;

        # Parachute condition
        if i > idxb
            if abs(x[3,i+1])-abs(x[3,i]) < 0
                Mode = 2;
            end
        end

        # Ground condition
        if x[3,i] < 0
            x[3,i] = 0;
        end

        # Stopping Criteria
        if i > idxb
            if x[3,i] < 0.5
                idxf = i;
                break;
            end
        end
    
    end
    
    return idxf, idxb, Tt[1:idxf], Tb, x[:,1:idxf], ang[:,1:idxf], Ft[:,1:idxf], Fa[:,1:idxf], M[:,1:idxf];


end