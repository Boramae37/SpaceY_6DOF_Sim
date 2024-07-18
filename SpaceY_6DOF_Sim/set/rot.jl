function b2n(quat)
    # Body 좌표계에서 NED 좌표계로 회전변환하는 행렬 (쿼터니안)
    q1 = quat[1];
    q2 = quat[2];
    q3 = quat[3];
    q4 = quat[4];

    m11 = 1-2*(q3*q3+q4*q4);    
    m12 = 2*(q2*q3+q1*q4);    
    m13 = 2*(q2*q4-q1*q3);
    m21 = 2*(q2*q3-q1*q4);      
    m22 = 1-2*(q2*q2+q4*q4);  
    m23 = 2*(q3*q4+q1*q2);
    m31 = 2*(q2*q4+q1*q3);      
    m32 = 2*(q3*q4-q1*q2);    
    m33 = 1-2*(q2*q2+q3*q3);

    return [m11 m21 m31 ;
            m12 m22 m32 ;
            m13 m23 m33];
end

function n2b(quat)
    q1 = quat[1];
    q2 = quat[2];
    q3 = quat[3];
    q4 = quat[4];

    m11 = 1-2*(q3*q3+q4*q4);    
    m12 = 2*(q2*q3+q1*q4);    
    m13 = 2*(q2*q4-q1*q3);
    m21 = 2*(q2*q3-q1*q4);      
    m22 = 1-2*(q2*q2+q4*q4);  
    m23 = 2*(q3*q4+q1*q2);
    m31 = 2*(q2*q4+q1*q3);      
    m32 = 2*(q3*q4-q1*q2);    
    m33 = 1-2*(q2*q2+q3*q3);

    return [m11 m12 m13 ;
            m21 m22 m23 ;
            m31 m32 m33];
end

function e2q(ang)
    phi = ang[1];
    theta = ang[2];
    psi = ang[3];

    q1 = cos(phi/2)*cos(theta/2)*cos(psi/2)-sin(phi/2)*sin(theta/2)*sin(psi/2);
    q2 = cos(phi/2)*sin(theta/2)*sin(psi/2)+sin(phi/2)*cos(theta/2)*cos(psi/2);
    q3 = cos(phi/2)*sin(theta/2)*cos(psi/2)+sin(phi/2)*cos(theta/2)*sin(psi/2);
    q4 = cos(phi/2)*cos(theta/2)*sin(psi/2)-sin(phi/2)*sin(theta/2)*cos(psi/2);

    return [q1; q2; q3; q4];
end

function q2e(quat)
    q1 = quat[1];
    q2 = quat[2];
    q3 = quat[3];
    q4 = quat[4];

    m11 = 1-2*(q3*q3+q4*q4);    
    m12 = 2*(q2*q3+q1*q4);    
    m13 = 2*(q2*q4-q1*q3);
    m22 = 1-2*(q2*q2+q4*q4);  
    m32 = 2*(q3*q4-q1*q2);    
    
    phi = atan(-m32,m22);
    theta = atan(-m13,m11);
    psi = asin(m12);
    
    return [phi; theta; psi];
end
