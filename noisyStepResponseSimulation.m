function noisyStepResponseSimulation(step,flowNo,T,systemFun,t0,tf,xs,u,p,Rvv,Rdd,mud,title_)

stepSize = u(flowNo)*(step-1);
u(flowNo) = u(flowNo)*step;
[T, H] = noisyNonlinearSimulation(systemFun,T,xs,t0,tf,u,p,Rvv,Rdd,mud);
plotTankHeights(T, H./stepSize, title_);

end