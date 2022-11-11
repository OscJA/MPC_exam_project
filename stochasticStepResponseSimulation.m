function stochasticStepResponseSimulation(step,flowNo,T,systemFun,t0,tf,xs,u,d,p,Rvv,title_)

stepSize = u(flowNo)*(step-1);
u(flowNo) = u(flowNo)*step;
[T, H] = stochasticNonlinearSimulation(systemFun,T,xs,t0,tf,u,d,p,Rvv);
plotTankHeights(T, H./stepSize, title_);

end