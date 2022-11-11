function [T, X] = stochasticNonlinearSimulation(systemFun,step,x0,t0,tf,us,p,Rvv)

T = t0:step:tf;

rho = p(12,1); % Density of water [g/cm3]
A = p(5:8,1)'; % Tank cross sectional areas [cm2]
X = zeros(length((tf-t0)/step), 4);
x = x0;
p2 = p;
p2 = [p2; step];
p2(13,1) = step;

i = 1;
for t=t0:step:tf
    xdot = systemFun(t,x,us,p2);
    v = Rvv*randn(4,1);
    y = FourTankSystemSensor(x+xdot*step, p2) + v;
    x = y./(rho*A');
    X(i, :) = x;
    i = i+1;
end

end