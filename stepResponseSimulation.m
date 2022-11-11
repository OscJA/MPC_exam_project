function [T, H] = stepResponseSimulation(step,flowNo,systemFun,t0,tf,xs,u,d,p)

stepSize = u(flowNo)*(step-1);
u(flowNo) = u(flowNo)*step;
[T, H] = nonlinear_simulation(systemFun,t0,tf,xs,u,d,p);
H = H./stepSize;

end