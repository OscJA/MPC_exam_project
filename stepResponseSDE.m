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
p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho; xs; ys; us; d];

%% Solve the stochastic problem
t0 = 0;
tf = 15*60;
Ts = 4;
Rvv = 0.05*eye(4);
Rdd = 5*eye(2);
mud = d;
Sigma = 2*[zeros(4,2); eye(2)];

%% Linearize the system
hs = ys; At = [A1; A2; A3; A4]; ap = [a1; a2; a3; a4]; gam = [gamma1; gamma2];
T = sqrt(2*At.*xs)./(ap.*sqrt(rho*g));
A = [-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
B = [rho*gam(1) 0;0 rho*gam(2); 0 rho*(1-gam(2)); rho*(1-gam(1)) 0];
C = diag(1./(rho*At));
Cz = C(1:2,:);

% Include the disturbance variables in the matrices
A = [[A; zeros(2,4)], [zeros(2,2); [rho, 0; 0, rho]; zeros(2,2)]];
B = [B; zeros(2,2)];
C = [C, zeros(4,2)];

M = expm([A, B; zeros(2,8)]*Ts);
Ad = M(1:6, 1:6);
Bd = M(1:6, 7:8);
Cd = C;

RES = expm([-A, Sigma*Sigma'; zeros(size(A)), A']*Ts);
Qd = Ad*RES(1:6, 7:end);
Qd_chol = chol(Qd);

xs = [xs; zeros(2,1)];


%% Calculate and plot F1 responses
% 10% F1 response
[T10_1, H10_1] = stepResponseSimulationStochastic(1.1,1,Ts,t0,tf,xs,us,p,Rvv,Qd_chol,Ad,Bd,Cd);

% 25% F1 response
[T25_1, H25_1] = stepResponseSimulationStochastic(1.25,1,Ts,t0,tf,xs,us,p,Rvv,Qd_chol,Ad,Bd,Cd);

% 50% F1 response
[T50_1, H50_1] = stepResponseSimulationStochastic(1.5,1,Ts,t0,tf,xs,us,p,Rvv,Qd_chol,Ad,Bd,Cd);

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
saveas(fig, '../Exam project/Figures/SDEstochastic_f1.png')


%% Calculate and plot F2 responses
% 10% F2 response
[T10_2, H10_2] = stepResponseSimulationStochastic(1.1,2,Ts,t0,tf,xs,us,p,Rvv,Qd_chol,Ad,Bd,Cd);

% 25% F1 response
[T25_2, H25_2] = stepResponseSimulationStochastic(1.25,2,Ts,t0,tf,xs,us,p,Rvv,Qd_chol,Ad,Bd,Cd);

% 50% F1 response
[T50_2, H50_2] = stepResponseSimulationStochastic(1.5,2,Ts,t0,tf,xs,us,p,Rvv,Qd_chol,Ad,Bd,Cd);

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
saveas(fig, '../Exam project/Figures/SDEstochastic_f2.png')


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
saveas(fig, '../Exam project/Figures/SDEstochastic_all_in_one.png')

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

