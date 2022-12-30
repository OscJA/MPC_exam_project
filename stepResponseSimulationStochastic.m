function [T, H] = stepResponseSimulationStochastic(step,flowNo,T,t0,tf,xs,u,p,Rvv,Qd_chol,Ad,Bd,Cd)
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
[T, H] = stochasticSimulation(xs,t0,tf,T,u,p,Rvv,Qd_chol,Ad,Bd,Cd);

% Calculate the normalized step change
H = H./stepSize;

end