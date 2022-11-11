function [T, H] = nonlinear_simulation(systemFun,t0,tf,x0,us,d,p)

[T,X] = ode15s(systemFun,[t0 tf],x0,[],us,d,p);

% help variables
[nT,nX] = size(X);
rho = p(12,1); % Density of water [g/cm3]
A = p(5:8,1)'; % Tank cross sectional areas [cm2]

% Compute the measured variables
H = zeros(nT,nX);
for i=1:nT
    H(i,:) = (X(i,:)-x0')./(rho*A);
end

end