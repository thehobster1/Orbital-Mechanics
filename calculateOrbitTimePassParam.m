function orbitTimePass = calculateOrbitTimePassParam (m, a, e, t, rVec, vVec, r1, theta1)
    
    % Constants

    G = 6.6743e-11; %Nm^2/kg^2

    mu = m*G; %km^3/s^2

    tol = 0.001;

    %n = sqrt(mu/a^3);

    p = a*(1 - e^2);

    %M
    %disp(rem(t, 2*pi()*sqrt(a^3/mu)));

    M = sqrt(mu/a^3)*rem(t, 2*pi()*sqrt(a^3/mu));
    %disp(M*180/pi());
    
    %E

    E = 1;
    E0 = 0;

    j = 0;

    while (abs(E - E0)/E > tol)
        E0 = E;
        E = E0 - (E0 - e*sin(E0) - M)/(1 - e*cos(E0));
        j = j+1;
    end

    E = E*180/pi();

    % True anomaly
    theta = 2*atand(sqrt((1+e)/(1-e))*tand(E/2));

    t2 = sqrt(a^3/mu)*(E*pi()/180 - e*sind(E));

    %r
    r = a*(1 - e*cosd(E));

    %f and g functions
    f = 1 - r/p*(1 - cosd(theta - theta1));
    g = r*r1/sqrt(mu*p)*sind(theta - theta1);
    fdot = (dot(rVec, vVec)/(p*r1)*(1 - cosd(theta - theta1)) - 1/r1*sqrt(mu/p)*sin(theta - theta1));
    gdot = 1 - r1/p*(1 - cosd(theta - theta1));

    %rVec
    b = sqrt(a^2*(1 - e^2));
    r2Vec = f*rVec + g*vVec;
    r2VecP = [a*(cosd(E) - e); b*sind(E)];

    
    %v
    v = sqrt(2*mu/r - mu/a);
    v2Vec = fdot*rVec + gdot*vVec;
    %v2VecP = [-a^2*n/r*sind(E); a*b*n/r*cosd(E)];

    % Eccentricity vector
    %e2Vec = ((v^2 - mu/r)*r2Vec - dot(r2Vec, v2Vec)*v2Vec)/mu;

    
    
    %gamma
    %gamma = asind(dot(r2Vec, v2Vec)/(r*v));
    gamma = atan2d(v2Vec(3), sqrt(v2Vec(1)^2 + v2Vec(2)^2));

    orbitTimePass = struct('Theta', theta, 'r2', r, 'v2', v, 'E2', E,'Gamma2', gamma, 't2mintp', t2, 'f', f, 'g', g, 'fdot', fdot, 'gdot', gdot, 'r2VecP', r2VecP, 'j', j, 'r2Vec', r2Vec, 'v2Vec', v2Vec);

end