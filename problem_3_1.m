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

% F1 RESPONSES
% 10% F1 response
[T10, H10] = stepResponseSimulation(1.1,1,@FourTankSystem,t0,tf,xs,us,d,p);

% 25% F1 response
[T25, H25] = stepResponseSimulation(1.25,1,@FourTankSystem,t0,tf,xs,us,d,p);

% 50% F1 response
[T50, H50] = stepResponseSimulation(1.5,1,@FourTankSystem,t0,tf,xs,us,d,p);

fig = figure;

subplot(2,2,3);
hold on;
plot(T10, H10(:,1));
plot(T25, H25(:,1));
plot(T50, H50(:,1));
title('Tank 1')
hold off;

subplot(2,2,4);
hold on;
plot(T10, H10(:,2));
plot(T25, H25(:,2));
plot(T50, H50(:,2));
title('Tank 2')
hold off;
legend('10% step', '25% step', '50% step', 'Location', 'SouthEast');

subplot(2,2,1);
hold on;
plot(T10, H10(:,3));
plot(T25, H25(:,3));
plot(T50, H50(:,3));
title('Tank 3')
hold off;

subplot(2,2,2);
hold on;
plot(T10, H10(:,4));
plot(T25, H25(:,4));
plot(T50, H50(:,4));
title('Tank 4')
hold off;
sgtitle(fig, "Step responses to changes in F1 flow");
saveas(fig, '../Exam project/Figures/deterministic_f1.png')

% F2 RESPONSES
% 10% F2 response
[T10, H10] = stepResponseSimulation(1.1,2,@FourTankSystem,t0,tf,xs,us,d,p);

% 25% F2 response
[T25, H25] = stepResponseSimulation(1.25,2,@FourTankSystem,t0,tf,xs,us,d,p);

% 50% F2 response
[T50, H50] = stepResponseSimulation(1.5,2,@FourTankSystem,t0,tf,xs,us,d,p);

fig = figure;
subplot(2,2,3);
hold on;
plot(T10, H10(:,1));
plot(T25, H25(:,1));
plot(T50, H50(:,1));
title('Tank 1')
hold off;

subplot(2,2,4);
hold on;
plot(T10, H10(:,2));
plot(T25, H25(:,2));
plot(T50, H50(:,2));
title('Tank 2')
hold off;
legend('10% step', '25% step', '50% step', 'Location', 'SouthEast');

subplot(2,2,1);
hold on;
plot(T10, H10(:,3));
plot(T25, H25(:,3));
plot(T50, H50(:,3));
title('Tank 3')
hold off;

subplot(2,2,2);
hold on;
plot(T10, H10(:,4));
plot(T25, H25(:,4));
plot(T50, H50(:,4));
title('Tank 4')
hold off;
sgtitle(fig, "Step responses to changes in F2 flow");
saveas(fig, '../Exam project/Figures/deterministic_f2.png')

% figure;
% plot(0:1:900, 0.1769*(1-exp(-(0:1:900)/100)))
% hold on;
% plot(T10, H10(:,2));
% hold off;

%% PROBLEM 4
xdot = FourTankSystemLinear(0,xs,us,p);
