function [T, H] = nonlinear_simulation(systemFun,t0,tf,x0,us,d,p)

[T,X] = ode15s(systemFun,[t0 tf],x0,[],us,d,p);

% help variables
[nT,nX] = size(X);
rho = p(12,1); % Density of water [g/cm3]
A = p(5:8,1)'; % Tank cross sectional areas [cm2]

% Compute the measured variables
H = zeros(nT,nX);
for i=1:nT
    H(i,:) = (X(i,:)-x0')./(rho*A);
end

step = 4;
Y = zeros((tf-t0)/step, 4);
x = x0;

ys = p(17:20);
i = 1;
for t=t0:step:tf
    xdot = systemFun(t,x,us,d,p);
    % v = Rvv*randn(4,1);
    y = FourTankSystemSensor(x+xdot*step, p);
    Y(i, :) = y-ys;
    % x = y.*(rho*A');
    x = x + xdot*step;
    i = i+1;
end
H = Y;
T = t0:step:tf;

end