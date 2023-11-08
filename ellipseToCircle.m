function deltaV = ellipseToCircle(e, a, m, theta)
    % Constants

    G = 6.6743e-11; %Nm^2/kg^2

    mu = m*G; %km^3/s^2

    p = a*(1-e^2);
    
    r1mag = p/(1+e*cos(theta));

    E = real(acos((a-r1mag)/(a*e)));

    tmintp = sqrt(a^3/mu)*(E - e*sin(E));

    v1minMag = sqrt(2*(mu/r1mag-mu/(2*a)));

    r1 = [r1mag, 0];

    gamma1min = acos(sqrt(mu*p)/(r1mag*v1minMag));

    v1min = v1minMag*[sin(gamma1min) cos(gamma1min)];

    v1plusMag = sqrt(mu/r1mag);

    gamma1plus = acos(sqrt(mu*r1mag)/(r1mag*v1plusMag));

    v1plus = v1plusMag*[sin(gamma1plus) cos(gamma1plus)];

    deltaV = sqrt(v1minMag^2 + v1plusMag^2 - 2*v1plusMag*v1minMag*cos(gamma1plus-gamma1min));

    deltaVvec = v1plus - v1min;

    alpha = real(acos(dot(deltaVvec, v1min)/(deltaV*v1minMag)));
        
    deltaV = struct('E', rad2deg(E), "WaitTime", tmintp, 'r1mag', r1mag, 'r1', r1, 'v1min', v1min, 'gamma1min', rad2deg(gamma1min), 'v1plus', v1plus, 'gamma1plus', rad2deg(gamma1plus), 'deltaV', deltaV, 'deltaVvec', deltaVvec, 'alpha', rad2deg(alpha));
end