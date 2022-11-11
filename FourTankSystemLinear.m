function xdot = FourTankSystemLinear(t,x,u,p)
% FOURTANKSYSTEM Model dx/dt = f(t,x,u,p) for 4-tank System
%
% --------------------------------------------------------------
% Parameters
% --------------------------------------------------------------

m = x; % Mass of liquid in each tank [g]
F = u; % Flow rates in pumps [cm3/s]
ap = p(1:4,1); % Pipe cross sectional areas [cm2]
At = p(5:8,1); % Tank cross sectional areas [cm2]
gam = p(9:10,1); % Valve positions [-]
g = p(11,1); % Acceleration of gravity [cm/s2]
rho = p(12,1); % Density of water [g/cm3]
xs = p(13:16,1); % Steady state vector
ys = p(17:20,1);
dt = p(21,1);
us = p(22:23,1);

% --------------------------------------------------------------
% Linearization
% --------------------------------------------------------------
hs = ys;
T = (At./ap).*sqrt(hs*4/(2*g));
A = [-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
B = [rho*gam(1) 0;0 rho*gam(2); 0 rho*(1-gam(2)); rho*(1-gam(1)) 0];
C = diag(1./(rho*At));
Cz = C(1:2,:);

M = expm([A, B; zeros(2,6)]*dt);
Ad = M(1:4, 1:4);
Bd = M(1:4, 5:6);


xdot = Ad*(x-xs) + Bd*(u-us);
xdot = xdot + xs;

end