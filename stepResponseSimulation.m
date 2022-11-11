function stepResponseSimulation(step,flowNo,systemFun,t0,tf,xs,u,d,p,title_)

stepSize = u(flowNo)*(step-1);
u(flowNo) = u(flowNo)*step;
[T, H] = nonlinear_simulation(systemFun,t0,tf,xs,u,d,p);
plotTankHeights(T, H./stepSize, title_);

end