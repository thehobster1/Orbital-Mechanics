function orbital_parameters = calculate_orbital_parameters(r_vec, v_vec, m)
    % Constants
    G = 6.6743e-11; %Nm^2/kg^2
    mu = m*G; %km^3/s^2

    % Magnitudes
    r = norm(r_vec);
    v = norm(v_vec);

    % Angular momentum
    h_vec = cross(r_vec, v_vec);
    h = norm(h_vec);

    % Eccentricity vector
    %e_vec = ((v^2 - mu/r)*r_vec - dot(r_vec, v_vec)*v_vec)/mu;
    e_vec = cross(v_vec, h_vec) / mu - r_vec / r;
    eccentricity = norm(e_vec);

    % Semi-major axis
    a = h^2 / (mu * (1 - eccentricity^2));

    %Energy
    eta = -mu/(2*a);

    % Inclination
    inclination = acosd(h_vec(3)/h);

    % Longitude of ascending node
    n_vec = cross([0; 0; 1], h_vec);
    n = norm(n_vec);
    if n_vec(1) >= 0
        RAAN = acosd(n_vec(1)/n);
    else
        RAAN = 360 - acosd(n_vec(1)/n);
    end

    % Argument of periapsis
    if e_vec(3) >= 0
        argument_of_periapsis = acosd(dot(n_vec, e_vec)/(n*eccentricity));
    else
        argument_of_periapsis = 360 - acosd(dot(n_vec, e_vec)/(n*eccentricity));
    end

    %E
    E = acos((a - r)/(a*eccentricity));

    %M
    M = E - eccentricity*sin(E);
    
    %(t - tp)
    tmintp = sqrt(a^3/mu)*(E - eccentricity*sin(E));

    period = 2*pi()*sqrt(a^3/mu);

    p = a*(1-eccentricity^2);
    %Gamma
    gamma = acosd(sqrt(mu*p)/(r*v));
    rp = p/(1+eccentricity);

% True anomaly
    %true_anomaly = acosd(dot(e_vec, r_vec)/(eccentricity*r));
    true_anomaly = acosd((p/r-1)/eccentricity);

    theta = argument_of_periapsis + true_anomaly;

    % Orbital parameters structure
    orbital_parameters = struct("SemiMajorAxis", a, 'Eccentricity', eccentricity,  ...
        'Inclination', inclination, 'ArgumentOfPeriapsis', argument_of_periapsis, 'RAAN', RAAN, "r_m", r, "v_ms", v, 'Gamma', gamma, 'ThetaStar', true_anomaly, 'M', rad2deg(M), 'E', rad2deg(E), 'tmintp', tmintp, 'eta', eta, 'Period', period, 'rp', rp, 'theta', theta);
end
