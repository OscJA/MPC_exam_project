function [T, Y] = noisyNonlinearSimulation(systemFun,step,x0,t0,tf,us,p,Rvv,Rdd,mud)

T = t0:step:tf;

rho = p(12,1); % Density of water [g/cm3]
A = p(5:8,1)'; % Tank cross sectional areas [cm2]
Y = zeros((tf-t0)/step, 4);
x = x0;
ys = x0./(rho*A');
p2 = p;
p2 = [p2; step];
p2(13,1) = step;

i = 1;
for t=t0:step:tf
    xdot = systemFun(t,x,us,Rdd,mud,p2);
    v = Rvv*randn(4,1);
    y = FourTankSystemSensor(x+xdot*step, p2) + v;
    Y(i, :) = y-ys;
    % x = y.*(rho*A');
    x = x + xdot*step;
    i = i+1;
end

end