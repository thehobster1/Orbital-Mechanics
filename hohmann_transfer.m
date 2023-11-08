function [time_of_flight, delta_V_total, phaseAngle, synodicPeriod] = hohmann_transfer(initial_radius, final_radius, m)
    % Constants
    G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
    mu = G*m;

    % Calculate semi-major axes of the transfer and target orbits
    a_transfer = (initial_radius + final_radius) / 2;
    %a_target = (initial_radius + final_radius) / 2;

    % Calculate velocities at each orbit
    v_initial = sqrt(mu / initial_radius);
    v_transfer1 = sqrt(mu * (2 / initial_radius - 1 / a_transfer));
    v_transfer2 = sqrt(mu * (2 / final_radius - 1 / a_transfer));
    v_final = sqrt(mu / final_radius);

    % Calculate delta Vs
    delta_V1 = v_transfer1 - v_initial;
    delta_V2 = v_final - v_transfer2;
    delta_V_total = delta_V1 + delta_V2;

    % Calculate the time of flight for the transfer
    time_of_flight = pi * sqrt((a_transfer^3) / mu);

    % Calculate the phase angle
    phaseAngle = rad2deg(pi()-time_of_flight*sqrt(mu/(final_radius)^3));

    % Synodic Period
    etaI = sqrt(mu/(initial_radius)^3);
    etaF = sqrt(mu/(final_radius)^3);
    synodicPeriod = 2*pi()/(etaI - etaF);
end