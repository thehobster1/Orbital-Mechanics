function [a, p, e, TOF, rd, vd, gammad, thetaStard, ra, va, gammaa, thetaStara, delVd, alphad, delVa, alphaa,...
    viVec, vfVec, vdVec, vaVec] = lambertArcs(ri, rf, m, eta)

    % Constants
    G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
    mu = G*m;
    
    c = sqrt(ri^2 + rf^2 -2*ri*rf*cosd(eta));

    a = (ri + rf + c)/4;
    s = 2*a;

    spaceTriangleAlpha = 2*asin(sqrt(s/(2*a)));
    spaceTriangleBeta = 2*asin(sqrt((s - c)/(2*a)));

    p = 4*a*(s - ri)*(s - rf)/(c^2)*sin((spaceTriangleAlpha+spaceTriangleBeta)/2)^2;
    
    e = sqrt(1 - p/a);

    %Depart and Arival Radii
    rd = ri;
    ra = rf;

    %Depart and Arival Velocity
    vd = sqrt(2*mu*(1/rd - 1/(2*a)));
    va = sqrt(2*mu*(1/ra - 1/(2*a)));

    %Gamma (Assume gammaa is positive);
    gammad = acosd(sqrt(mu*p)/(rd*vd));
    gammaa = acosd(sqrt(mu*p)/(ra*va));

    if eta > 180
        gammad = -1*gammad;
        gammaa = -1*gammaa;
    end



    %ThetaStar
    thetaStard = acosd((p/ri - 1)/e);
    Vr = vd*sind(gammad);
    if Vr < 0
        thetaStard = -1*thetaStard;
    end
    thetaStara = thetaStard + eta;

    %TOF
    if eta <= 180
        TOF = sqrt(a^3/mu)*(spaceTriangleAlpha - spaceTriangleBeta - (sin(spaceTriangleAlpha) - sin(spaceTriangleBeta)));
    else
        TOF = sqrt(a^3/mu)*(2*pi() - (spaceTriangleAlpha - sin(spaceTriangleAlpha)) + (spaceTriangleBeta - sin(spaceTriangleBeta)));
    end

    %Delta V
    vi = sqrt(mu/ri);
    vf = sqrt(mu/rf);

    viVec = [1, 0]*vi;
    vfVec = [1, 0]*vf;

    vdVec = vd*[cosd(gammad), sind(gammad)];
    vaVec = va*[cosd(gammaa), sind(gammaa)];

    delVdVec = vdVec - viVec;
    delVaVec = vfVec - vaVec;

    delVd = norm(delVdVec);
    delVa = norm(delVaVec);

    %Alphad
    alphad = asind(sind(abs(gammad))*vd/delVd);
    betad = asind(sind(abs(gammad))*vi/delVd);

    if 180 - (alphad + betad + gammad) <= 0.01
        alphad = 180 - alphad;
    elseif 180 - (alphad + 180 - betad + gammad) <= 0.01
        alphad = 180 - alphad;
    end
    
    %Alphaa
    alphaa = asind(sind(abs(gammaa))*va/delVa);
    betaa = asind(sind(abs(gammaa))*vf/delVa);

    if 180 - (alphaa + betaa + gammaa) <= 0.01
        alphaa = 180 - alphaa;
    elseif 180 - (alphaa + 180 - betaa + gammaa) <= 0.01
        alphaa = 180 - alphaa;
    end

end