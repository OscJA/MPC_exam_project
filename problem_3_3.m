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
Sigma = 0.5*[zeros(4,2); eye(2)];

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


%% F1
% 10% F1 response
%[T10_1, H10_1] = noisyStepResponseSimulation(1.1,1,Ts,@noisyFourTankSystem,t0,tf,xs,us,p,Rvv,Rdd,mud);
[T10_1, H10_1] = stochasticStepResponseSimulation(1.1,1,Ts,t0,tf,xs,us,p,Rvv,Qd_chol,Ad,Bd,Cd);

% 25% F1 response
% [T25_1, H25_1] = noisyStepResponseSimulation(1.25,1,Ts,@noisyFourTankSystem,t0,tf,xs,us,p,Rvv,Rdd,mud);
[T25_1, H25_1] = stochasticStepResponseSimulation(1.25,1,Ts,t0,tf,xs,us,p,Rvv,Qd_chol,Ad,Bd,Cd);

% 50% F1 response
% [T50_1, H50_1] = noisyStepResponseSimulation(1.5,1,Ts,@noisyFourTankSystem,t0,tf,xs,us,p,Rvv,Rdd,mud);
[T50_1, H50_1] = stochasticStepResponseSimulation(1.5,1,Ts,t0,tf,xs,us,p,Rvv,Qd_chol,Ad,Bd,Cd);

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


%% F2

% 10% F1 response
% [T10_2, H10_2] = stochasticStepResponseSimulation(1.1,2,Ts,@noisyFourTankSystem,t0,tf,xs,us,p,Rvv,Rdd,mud);
[T10_2, H10_2] = stochasticStepResponseSimulation(1.1,2,Ts,t0,tf,xs,us,p,Rvv,Qd_chol,Ad,Bd,Cd);

% 25% F1 response
% [T25_2, H25_2] = stochasticStepResponseSimulation(1.25,2,Ts,@noisyFourTankSystem,t0,tf,xs,us,p,Rvv,Rdd,mud);
[T25_2, H25_2] = stochasticStepResponseSimulation(1.25,2,Ts,t0,tf,xs,us,p,Rvv,Qd_chol,Ad,Bd,Cd);

% 50% F1 response
% [T50_2, H50_2] = stochasticStepResponseSimulation(1.5,2,Ts,@noisyFourTankSystem,t0,tf,xs,us,p,Rvv,Rdd,mud);
[T50_2, H50_2] = stochasticStepResponseSimulation(1.5,2,Ts,t0,tf,xs,us,p,Rvv,Qd_chol,Ad,Bd,Cd);

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


%% All in one plot
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

%% Approx params

%% G21
[r_min, den_opt, num_min, s_opt, ts_22] = find_transfer_params(H10_2(:,1), Ts);

figure;
plot(ts_22, s_opt);
hold on;
plot(T10_2, H10_2(:,1));
hold off;
title('G21');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');
K21 = num_min/den_opt(1);
a = 1;
b = den_opt(2)/den_opt(1);
c = den_opt(3)/den_opt(1);
tau1_21 = (-b-sqrt(b^2-4*a*c))/(2*a);
tau2_21 = (-b+sqrt(b^2-4*a*c))/(2*a);

%% G21
[r_min, den_opt, num_min, s_opt, ts_22] = find_transfer_params(H10_1(:,2), Ts);

figure;
plot(ts_22, s_opt);
hold on;
plot(T10_1, H10_1(:,2));
hold off;
title('G12');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');
K12 = num_min/den_opt(1);
a = 1;
b = den_opt(2)/den_opt(1);
c = den_opt(3)/den_opt(1);
tau1_12 = (-b-sqrt(b^2-4*a*c))/(2*a);
tau2_12 = (-b+sqrt(b^2-4*a*c))/(2*a);


%% G11
[r_min, den_opt, num_min, s_opt, ts_22] = find_transfer_params(H10_1(:,1), Ts);

figure;
plot(ts_22, s_opt);
hold on;
plot(T10_1, H10_1(:,1));
hold off;
title('G11');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');
K11 = num_min/den_opt(2);
tau11 = -den_opt(3)/den_opt(2);


%% G22
[r_22, den_22, num_22, s_22, ts_22] = find_transfer_params(H10_2(:,2), Ts);

figure;
plot(ts_22, s_22);
hold on;
plot(T10_2, H10_2(:,2));
hold off;
title('G22');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');
K22 = num_22/den_22(2);
tau22 = -den_22(3)/den_22(2);

%% Load the params to latex

param_names = ["K"; "\\tau_1"; "\\tau_2"];
trans_fun = ["G_{11}"; "G_{12}"; "G_{21}"; "G_{22}"];
K = num2str([K11; K12; K21; K22]);
tau_1s = num2str([tau11; tau1_12; tau1_21; tau22]);
tau_2s = [""; tau2_12; tau2_21; ""];

% T = table(trans_fun,K,tau_1s,tau_2s);
% table2latex(T, '../Exam project/Tables/T.tex'); % params_sim.tex
% Ttex = table2latex(T, []); % params_sim.tex

Mat = [string(K11), string(K12), string(K21), string(K22);
    string(tau11), string(tau1_12), string(tau1_21), string(tau22);
    "", string(tau2_12), string(tau2_21), ""];

T2L(param_names, trans_fun, Mat, '../Exam project/Tables/params_sim_SDE.tex');
