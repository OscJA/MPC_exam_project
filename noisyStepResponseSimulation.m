function [T, H] = noisyStepResponseSimulation(step,flowNo,T,systemFun,t0,tf,xs,u,p,Rvv,Rdd,mud)

stepSize = u(flowNo)*(step-1);
u(flowNo) = u(flowNo)*step;
[T, H] = noisyNonlinearSimulation(systemFun,T,xs,t0,tf,u,p,Rvv,Rdd,mud);
H = H./stepSize;

end