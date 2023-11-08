function [inclination, RAAN] = transferPlane(r1, r2)
    % r1 and r2 are 3x1 position vectors in Cartesian coordinates

    % Calculate the angular momentum vector (h)
    h = cross(r1, r2);

    % Calculate the magnitude of h
    h_magnitude = norm(h);

    % Calculate the inclination (i)
    inclination = acos(h(3) / h_magnitude);

    % Calculate the unit vector along the line of nodes (N)
    N = cross([0; 0; 1], h);
    N_magnitude = norm(N);

    % Calculate the Right Ascension of Ascending Node (RAAN)
    if N(2) >= 0
        RAAN = acos(N(1) / N_magnitude);
    else
        RAAN = 2 * pi - acos(N(1) / N_magnitude);
    end

    % Convert inclination and RAAN from radians to degrees
    inclination = rad2deg(inclination);
    RAAN = rad2deg(RAAN);
end