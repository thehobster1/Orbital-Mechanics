function gamma = flightPathAngle (R, V, e, theta)
    gamma = atand(e*sind(theta)/(1+e*cosd(theta)));
    %{
    vr = dot(R, V);

    if vr < 0
        gamma = -1*gamma;
    end
    %}
end