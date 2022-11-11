clear;
close all;

%% Set up the system
a1 = 1.2272; %[cm2] Area of outlet pipe 1
a2 = 1.2272; %[cm2] Area of outlet pipe 2
a3 = 1.2272; %[cm2] Area of outlet pipe 3
a4 = 1.2272; %[cm2] Area of outlet pipe 4
A1 = 380.1327; %[cm2] Cross sectional area of tank 1
A2 = 380.1327; %[cm2] Cross sectional area of tank 2
A3 = 380.1327; %[cm2] Cross sectional area of tank 3
A4 = 380.1327; %[cm2] Cross sectional area of tank 4
gamma1 = 0.58; % Flow distribution constant. Valve 1
gamma2 = 0.68; % Flow distribution constant. Valve 2
g = 981; %[cm/s2] The acceleration of gravity
rho = 1.00; %[g/cm3] Density of water
p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];

us = [300; 300]; % [cm3/s] Steady state flow
d = [50; 50];
xs0 = [5000; 5000; 5000; 5000]; % [g] Initial guess on xs
xs = fsolve(@FourTankSystemWrap,xs0,[],us,d,p); % Steady state masses
ys = FourTankSystemSensor(xs,p); % Steady states heights
p = [p; xs; ys; us; d];

%% Solve the deterministic problem over time
t0 = 0;
tf = 15*60;

% 10% F1 response
stepResponseSimulation(1.1,1,@FourTankSystem,t0,tf,xs,us,d,p,"Simulation of 10% step response of F1")

% 10% F2 response
stepResponseSimulation(1.1,2,@FourTankSystem,t0,tf,xs,us,d,p,"Simulation of 10% step response of F2")

% 25% F1 response
stepResponseSimulation(1.25,1,@FourTankSystem,t0,tf,xs,us,d,p,"Simulation of 25% step response of F1")

% 25% F2 response
stepResponseSimulation(1.25,2,@FourTankSystem,t0,tf,xs,us,d,p,"Simulation of 25% step response of F2")

% 50% F1 response
stepResponseSimulation(1.5,1,@FourTankSystem,t0,tf,xs,us,d,p,"Simulation of 50% step response of F1")

% 50% F2 response
stepResponseSimulation(1.5,2,@FourTankSystem,t0,tf,xs,us,d,p,"Simulation of 50% step response of F2")

%% PROBLEM 4
xdot = FourTankSystemLinear(0,xs,us,p);
