function [T, H] = stepResponseSimulationNoisy(step,flowNo,T,systemFun,t0,tf,xs,u,p,Rvv,Rdd,mud)
%% 
% From the heights and timestep, find the best fitting transfer function
% parameters by using Brute force
% Author: Oscar Juul Andersen, s194316
%%

% Determine the size of the step change
stepSize = u(flowNo)*(step-1);

% Do the step change in the flow
u(flowNo) = u(flowNo)*step;

% Simulate using the systemFun in the interval [t0, tf]
[T, H] = nonlinearSimulationNoisy(systemFun,T,xs,t0,tf,u,p,Rvv,Rdd,mud);

% Calculate the normalized step change
H = H./stepSize;

end