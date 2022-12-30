function [T, Y] = stochasticSimulation(x0,t0,tf,step,u,p,Rvv,Qd_chol,Ad,Bd,Cd)
%% 
% Perform a nonlinear simulation of systemFun from t0 to tf, using
% explicicit Euler
% Author: Oscar Juul Andersen, s194316
%%

% Initialize and extract variables
T = t0:step:tf;

rho = p(12,1); % Density of water [g/cm3]
A = p(5:8,1)'; % Tank cross sectional areas [cm2]
Y = zeros((tf-t0)/step, 4);
ys = x0(1:4)./(rho*A');
us = p(21:22);
x = zeros(6,1);

% Simulate the system
i = 1;
for t=t0:step:tf
    w = Qd_chol*randn(6, 1);
    x = Ad*x + Bd*(u-us);
    v = Rvv*randn(4,1);
    y = Cd*x + v;
    Y(i, :) = y;
    i = i+1;
end

end