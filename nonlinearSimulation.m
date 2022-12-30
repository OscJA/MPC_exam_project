function [T, H] = nonlinearSimulation(systemFun,t0,tf,x0,us,d,p)
%% 
% Perform a nonlinear simulation of systemFun from t0 to tf, using
% explicicit Euler
% Author: Oscar Juul Andersen, s194316
%%

% Initialize the variables
step = 4;
Y = zeros((tf-t0)/step, 4);
x = x0;

% Extract the steady state of the heights
ys = p(17:20);
i = 2;

% Perform the explicit Euler algorithm
for t=t0:step:tf-step
    xdot = systemFun(t,x,us,d,p);
    y = FourTankSystemSensor(x+xdot*step, p);
    Y(i, :) = y-ys;
    x = x + xdot*step;
    i = i+1;
end
H = Y;
T = t0:step:tf;

end