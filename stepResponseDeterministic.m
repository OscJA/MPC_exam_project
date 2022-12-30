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

%% Solve the deterministic problem over time
t0 = 0;
tf = 15*60;
Ts = 4;

%% Calculate and plot F1 responses
% 10% F1 response
[T10_1, H10_1] = stepResponseSimulation(1.1,1,@FourTankSystem,t0,tf,xs,us,d,p);

% 25% F1 response
[T25_1, H25_1] = stepResponseSimulation(1.25,1,@FourTankSystem,t0,tf,xs,us,d,p);

% 50% F1 response
[T50_1, H50_1] = stepResponseSimulation(1.5,1,@FourTankSystem,t0,tf,xs,us,d,p);

fig = figure;

subplot(2,2,3);
hold on;
plot(T10_1, H10_1(:,1));
plot(T25_1, H25_1(:,1));
plot(T50_1, H50_1(:,1));
title('Tank 1')
hold off;

subplot(2,2,4);
hold on;
plot(T10_1, H10_1(:,2));
plot(T25_1, H25_1(:,2));
plot(T50_1, H50_1(:,2));
title('Tank 2')
hold off;
legend('10% step', '25% step', '50% step', 'Location', 'SouthEast');

subplot(2,2,1);
hold on;
plot(T10_1, H10_1(:,3));
plot(T25_1, H25_1(:,3));
plot(T50_1, H50_1(:,3));
title('Tank 3')
hold off;

subplot(2,2,2);
hold on;
plot(T10_1, H10_1(:,4));
plot(T25_1, H25_1(:,4));
plot(T50_1, H50_1(:,4));
title('Tank 4')
hold off;
sgtitle(fig, "Step responses to changes in F1 flow");
saveas(fig, '../Exam project/Figures/deterministic_f1.png')

%% Calculate and plot F2 responses
% 10% F2 response
[T10_2, H10_2] = stepResponseSimulation(1.1,2,@FourTankSystem,t0,tf,xs,us,d,p);

% 25% F2 response
[T25_2, H25_2] = stepResponseSimulation(1.25,2,@FourTankSystem,t0,tf,xs,us,d,p);

% 50% F2 response
[T50_2, H50_2] = stepResponseSimulation(1.5,2,@FourTankSystem,t0,tf,xs,us,d,p);

fig = figure;
subplot(2,2,3);
hold on;
plot(T10_2, H10_2(:,1));
plot(T25_2, H25_2(:,1));
plot(T50_2, H50_2(:,1));
title('Tank 1')
hold off;

subplot(2,2,4);
hold on;
plot(T10_2, H10_2(:,2));
plot(T25_2, H25_2(:,2));
plot(T50_2, H50_2(:,2));
title('Tank 2')
hold off;
legend('10% step', '25% step', '50% step', 'Location', 'SouthEast');

subplot(2,2,1);
hold on;
plot(T10_2, H10_2(:,3));
plot(T25_2, H25_2(:,3));
plot(T50_2, H50_2(:,3));
title('Tank 3')
hold off;

subplot(2,2,2);
hold on;
plot(T10_2, H10_2(:,4));
plot(T25_2, H25_2(:,4));
plot(T50_2, H50_2(:,4));
title('Tank 4')
hold off;
sgtitle(fig, "Step responses to changes in F2 flow");
saveas(fig, '../Exam project/Figures/deterministic_f2.png')

%% Plot step responses in tank 1 and 2, for changes in F1 and F2 in one plot
fig = figure;
subplot(2,2,1);
hold on;
plot(T10_1, H10_1(:,1));
plot(T25_1, H25_1(:,1));
plot(T50_1, H50_1(:,1));
title('Flow 1 to tank 1');
hold off;

subplot(2,2,2);
hold on;
plot(T10_1, H10_1(:,2));
plot(T25_1, H25_1(:,2));
plot(T50_1, H50_1(:,2));
title('Flow 2 to tank 1');
hold off;

subplot(2,2,3);
hold on;
plot(T10_2, H10_2(:,1));
plot(T25_2, H25_2(:,1));
plot(T50_2, H50_2(:,1));
title('Flow 1 to tank 2');
hold off;

subplot(2,2,4);
hold on;
plot(T10_2, H10_2(:,2));
plot(T25_2, H25_2(:,2));
plot(T50_2, H50_2(:,2));
title('Flow 2 to tank 2');
legend('10% step', '25% step', '50% step', 'Location', 'SouthEast');
hold off;
saveas(fig, '../Exam project/Figures/deterministic_all_in_one.png')
%% Estimate the gains, poles and zeros

% G21
[r_21, den_21, num_21, s_21, ts_21] = findTransferParameters(H10_2(:,1), Ts);

figure;
plot(ts_21, s_21);
hold on;
plot(T10_2, H10_2(:,1));
hold off;
title('G21');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');
K21 = num_21/den_21(1);
a = 1;
b = den_21(2)/den_21(1);
c = den_21(3)/den_21(1);
tau1_21 = (-b-sqrt(b^2-4*a*c))/(2*a);
tau2_21 = (-b+sqrt(b^2-4*a*c))/(2*a);

% G12
[r_12, den_12, num_12, s_12, ts_12] = findTransferParameters(H10_1(:,2), Ts);

figure;
plot(ts_12, s_12);
hold on;
plot(T10_1, H10_1(:,2));
hold off;
title('G12');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');
K12 = num_12/den_12(1);
a = 1;
b = den_12(2)/den_12(1);
c = den_12(3)/den_12(1);
tau1_12 = (-b-sqrt(b^2-4*a*c))/(2*a);
tau2_12 = (-b+sqrt(b^2-4*a*c))/(2*a);


% G11
[r_11, den_11, num_11, s_11, ts_11] = findTransferParameters(H10_1(:,1), Ts);

figure;
plot(ts_11, s_11);
hold on;
plot(T10_1, H10_1(:,1));
hold off;
title('G11');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');
K11 = num_11/den_11(2);
tau11 = -den_11(3)/den_11(2);


% G22
[r_22, den_22, num_22, s_22, ts_22] = findTransferParameters(H10_2(:,2), Ts);

figure;
plot(ts_22, s_22);
hold on;
plot(T10_2, H10_2(:,2));
hold off;
title('G22');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');
K22 = num_22/den_22(2);
tau22 = -den_22(3)/den_22(2);


%% Find linear model from param estimates
tol = 1e-5;

[Ad11,Bd11,Cd11,Dd11,sH11] = mimoctf2dss({num_11},{den_11(2:3)},1,Ts,tf,tol);

[Ad22,Bd22,Cd22,Dd22,sH22] = mimoctf2dss({num_22},{den_22(2:3)},1,Ts,tf,tol);


% Find the last params
[r_13, den_13, num_13, s_13, ts_13] = findTransferParameters(H10_1(:,3), Ts);
[r_14, den_14, num_14, s_14, ts_14] = findTransferParameters(H10_1(:,4), Ts);

[r_23, den_23, num_23, s_23, ts_23] = findTransferParameters(H10_2(:,3), Ts);
[r_24, den_24, num_24, s_24, ts_24] = findTransferParameters(H10_2(:,4), Ts);

%% Calculate the impulse responses (Markov parameters)
numerators = [{num_11}, {num_12}, {num_13}, {num_14}; {num_21}, {num_22}, {num_23}, {num_24}]';
denominators = [{den_11}, {den_12}, {den_13}, {den_14}; {den_21}, {den_22}, {den_23}, {den_24}]';
t_delays = zeros(2,4)';

[Ad,Bd,Cd,Dd,sH] = mimoctf2dss(numerators,denominators,t_delays,Ts,tf,tol);

% Plots the step response for all tanks and flow changes
N = tf/Ts;
H = mimodss2dimpulse(Ad,Bd,Cd,Dd,N);
fig = figure;
plot(reshape(H(1,1,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{11}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/deterministic_H11.png')

fig = figure;
plot(reshape(H(1,2,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{12}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/deterministic_H12.png')

fig = figure;
plot(reshape(H(2,1,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{21}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/deterministic_H21.png')

fig = figure;
plot(reshape(H(2,2,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{22}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/deterministic_H22.png')

fig = figure;
plot(reshape(H(3,1,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{31}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/deterministic_H31.png')

fig = figure;
plot(reshape(H(3,2,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{32}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/deterministic_H32.png')

fig = figure;
plot(reshape(H(4,1,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{41}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/deterministic_H41.png')

fig = figure;
plot(reshape(H(4,2,:), 1, []), 'LineWidth', 2)
set(gca, "FontSize", 14)
title('$H_{42}$', 'FontSize', 20)
xlabel('Time step', 'FontSize', 16)
saveas(fig, '../Exam project/Figures/deterministic_H42.png')


%% Calculate the linearization of the system
[Ad11,Bd11,Cd11,Dd11,l11] = sisoctf2dss(num_11,den_11,0,Ts);
[Ad12,Bd12,Cd12,Dd12,l12] = sisoctf2dss(num_12,den_12,0,Ts);
[Ad13,Bd13,Cd13,Dd13,l13] = sisoctf2dss(num_13,den_13,0,Ts);
[Ad14,Bd14,Cd14,Dd14,l14] = sisoctf2dss(num_14,den_14,0,Ts);
[Ad21,Bd21,Cd21,Dd21,l21] = sisoctf2dss(num_21,den_21,0,Ts);
[Ad22,Bd22,Cd22,Dd22,l22] = sisoctf2dss(num_22,den_22,0,Ts);
[Ad23,Bd23,Cd23,Dd23,l23] = sisoctf2dss(num_23,den_23,0,Ts);
[Ad24,Bd24,Cd24,Dd24,l24] = sisoctf2dss(num_24,den_24,0,Ts);

A = [Ad11(1,1), Ad12(2,1); Ad21(2,1), Ad22(1,1)];


