function [T, H] = stochasticStepResponseSimulation(step,flowNo,T,t0,tf,xs,u,p,Rvv,Qd_chol,Ad,Bd,Cd)

stepSize = u(flowNo)*(step-1);
u(flowNo) = u(flowNo)*step;
[T, H] = stochasticSimulation(xs,t0,tf,T,u,p,Rvv,Qd_chol,Ad,Bd,Cd);
H = H./stepSize;

end