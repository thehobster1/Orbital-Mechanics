clear all
close all

re = 6378100;
me = 5.97*10^24;
G = 6.6743*10^-11;
mu = G*me;
A = 200000;

a = 1.4*re;

X0 = [0; 0; 0];

Xdot0 = [0.05; -0.01; 0];

phase = 0.6;

[XTF, XdotTF, x, y] = forceFree(X0, Xdot0, a, mu, 50*60);

function [XTF, XdotTF, x, y] = forceFree(X0, Xdot0, a, mu, TOF)

    n = sqrt(mu/a^3);

    tf = TOF;
    
    A = [4-3*cos(n*tf) 0 0; 6*(sin(n*tf)-n*tf) 1 0; 0 0 cos(n*tf)];
    
    B = [1/n*sin(n*tf) 2/n*(1-cos(n*tf)) 0; -2/n*(1-cos(n*tf)) 4/n*sin(n*tf)-3*tf 0; 0 0 1/n*sin(n*tf)];

    dA = [3*n*sin(n*tf) 0 0; 6*(n*cos(n*tf)-n) 0 0; 0 0 -n*sin(n*tf)];

    dB = [cos(n*tf) 2*sin(n*tf) 0; -2*sin(n*tf) 4*cos(n*tf)-3 0; 0 0 cos(n*tf)];

    %Xdot0 = inv(B)*(Xtf - A*X0);

    XTF = A*X0 + B*Xdot0;



    XdotTF = dA*X0 + dB*Xdot0;

    [x, y] = plotEuler(n, TOF, X0, Xdot0);

end

function [x, y] = plotEuler(n, tf, X0, Xdot0)

    t = 0:1:tf;

    %{

    Aplot = [(4-3*cos(n*t)) 0 0; 6*(sin(n*t)-n*t) 1 0; 0 0 cos(n*t)];

    Bplot = [(1/n*sin(n*t)) (2/n*(1-cos(n*t))) 0];
    Cplot = [(-2/n*(1-cos(n*t))) (4/n*sin(n*t)-3*t) 0];
    Dplot = [0 0 (1/n*sin(n*t))];
    %}
    %{
    x1 = [(4-3*cos(n.*t)) 0 0; 6*(sin(n.*t)-n.*t) 1 0; 0 0 cos(n.*t)].*X0;
    x2 = [(1/n*sin(n.*t)) (2/n*(1-cos(n.*t)))  0; (-2/n*(1-cos(n.*t))) (4/n*sin(n.*t)-3.*t) 0; 0 0 (1/n*sin(n.*t))].*Xdot0;
    %}
    
    x = Xdot0(1)/n*sin(n.*t) - (2*Xdot0(2)/n + 3*X0(1))*cos(n.*t) + 2*Xdot0(2)/n + 4*X0(1);

    y = 2*Xdot0(1)/n*cos(n.*t) - (4*Xdot0(2)/n + 6*X0(1))*sin(n.*t) + X0(2) - 2*Xdot0(1)/n -(3*Xdot0(2) + 6*n*X0(1)).*t;

    xmid = Xdot0(1)/n*sin(n*50*60) - (2*Xdot0(2)/n + 3*X0(1))*cos(n*50*60) + 2*Xdot0(2)/n + 4*X0(1);

    ymid = 2*Xdot0(1)/n*cos(n*50*60) - (4*Xdot0(2)/n + 6*X0(1))*sin(n*50*60) + X0(2) - 2*Xdot0(1)/n -(3*Xdot0(2) + 6*n*X0(1))*50*60;

    figure(1);
    hold on
    grid on
    axis equal
    plot(y,x);
    plot(xmid, ymid, 'o');
    xlabel('y (m)');
    ylabel('x (m)');

end