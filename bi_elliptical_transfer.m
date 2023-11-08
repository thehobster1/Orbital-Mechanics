function [time_of_flight, delta_V_total, v_transfer1_2] = bi_elliptical_transfer(initial_radius, intermediate_radius, final_radius, m)
    % Constants
    G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
    mu = m*G;

    % Calculate semi-major axes of the transfer and target orbits
    a_transfer1 = (initial_radius + intermediate_radius) / 2;
    a_transfer2 = (intermediate_radius + final_radius) / 2;

    % Calculate velocities at each orbit
    v_initial = sqrt(mu / initial_radius);
    v_transfer1_1 = sqrt(mu * (2 / initial_radius - 1 / a_transfer1));
    v_transfer1_2 = sqrt(mu * (2 / intermediate_radius - 1 / a_transfer1));
    v_transfer2_1 = sqrt(mu * (2 / intermediate_radius - 1 / a_transfer2));
    v_transfer2_2 = sqrt(mu * (2 / final_radius - 1 / a_transfer2));
    v_final = sqrt(mu / final_radius);

    % Calculate delta Vs
    delta_V1 = v_transfer1_1 - v_initial + v_transfer2_2 - v_final;
    disp(v_transfer1_1 - v_initial);
    disp(v_transfer2_1 - v_transfer1_2);
    delta_V2 = v_transfer1_2 - v_transfer2_1;
    delta_V_total = delta_V1 + delta_V2;

    % Calculate the time of flight for the transfer
    time_of_flight = pi * sqrt((a_transfer1^3) / mu) + pi * sqrt((a_transfer2^3) / mu);
end