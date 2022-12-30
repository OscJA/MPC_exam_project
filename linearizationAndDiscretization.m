clear;
close all;
addpath("../Realization/")

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

%% Pack the variables for further calculations
At = [A1; A2; A3; A4];
ap = [a1; a2; a3; a4];
gam = [gamma1; gamma2];
dt = 4; % Time step

%% Linearization
hs = ys;
T = sqrt(2*At.*xs)./(ap.*sqrt(rho*g));
A = [-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
B = [rho*gam(1) 0;0 rho*gam(2); 0 rho*(1-gam(2)); rho*(1-gam(1)) 0];
E = [
    0, 0;
    0, 0;
    rho, 0;
    0, rho
    ];
C = diag(1./(rho*At));
Cz = C(1:2,:);


%% Discretize
M = expm([A, B; zeros(2,6)]*dt);
Ad1 = M(1:4, 1:4);
Bd1 = M(1:4, 5:6);

M = expm([A, B, E; zeros(4,8)]*dt);
Ad = M(1:4, 1:4);
Bd = M(1:4, 5:6);
Ed = M(1:4, 7:8);
Cd = C;


%% Markov parameters
N = 250;
H = zeros(4, 2, N+1);
H(:, :, 1) = zeros(4, 2);
Obs_matrix = Cd;
for i=1:N
    Obs_matrix = Obs_matrix*Ad;
    H(:, :, i+1) = Obs_matrix*Bd;
end
Dd = zeros(4, 2);
H = mimodss2dimpulse(Ad,Bd,Cd,Dd,N);

fig = figure;
plot(reshape(H(1,1,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{11}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/linear_H11.png')

fig = figure;
plot(reshape(H(1,2,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{12}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/linear_H12.png')

fig = figure;
plot(reshape(H(2,1,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{21}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/linear_H21.png')

fig = figure;
plot(reshape(H(2,2,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{22}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/linear_H22.png')

fig = figure;
plot(reshape(H(3,1,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{31}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/linear_H31.png')

fig = figure;
plot(reshape(H(3,2,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{32}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/linear_H32.png')

fig = figure;
plot(reshape(H(4,1,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{41}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/linear_H41.png')

fig = figure;
plot(reshape(H(4,2,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{42}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/linear_H42.png')

