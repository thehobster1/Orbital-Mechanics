function [deltaV_escape, deltaV_insert, TOF] = orbitalTransferGravityFields(earthOrbit, jupiterOrbit)
    % Constants
    mu_earth = 398600;  % Earth gravitational parameter (km^3/s^2)
    mu_jupiter = 126686511;  % Jupiter gravitational parameter (km^3/s^2)
    mu_sun = 132712440018;  % Sun gravitational parameter (km^3/s^2)

    % Extract orbital parameters
    r1 = earthOrbit(1);  % periapsis radius around Earth (km)
    e1 = earthOrbit(2);  % eccentricity around Earth
    r2 = jupiterOrbit(1);  % periapsis radius around Jupiter (km)
    e2 = jupiterOrbit(2);  % eccentricity around Jupiter

    % Calculate velocities at departure and arrival points
    a1 = r1 / (1 - e1);  % semi-major axis for Earth orbit (km)
    a2 = r2 / (1 - e2);  % semi-major axis for Jupiter orbit (km)

    % Calculate deltaV to escape Earth's orbit
    v_escape = sqrt(mu_earth / r1);  % velocity for escape (km/s)
    deltaV_escape = abs(v_escape - sqrt(mu_sun / r1 + mu_earth / a1));

    % Calculate deltaV to be placed in orbit around Jupiter
    v_insert = sqrt(mu_jupiter / r2);  % velocity for insertion (km/s)
    deltaV_insert = abs(v_insert - sqrt(mu_sun / r2 + mu_jupiter / a2));

    % Calculate time of flight (TOF)
    TOF = pi * sqrt((a1^3 + a2^3) / (8 * mu_sun));

    fprintf('DeltaV to escape Earth orbit: %.2f km/s\n', deltaV_escape);
    fprintf('DeltaV to insert into Jupiter orbit: %.2f km/s\n', deltaV_insert);
    fprintf('Total DeltaV: %.2f km/s\n', deltaV_escape + deltaV_insert);
    fprintf('Time of Flight: %.2f days\n', TOF / (24 * 3600));
end
