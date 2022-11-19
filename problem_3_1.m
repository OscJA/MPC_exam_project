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

% F1 RESPONSES
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

%% F2 RESPONSES
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
saveas(fig, '../Exam project/Figures/deterministic_all_in_one.png')
%% Approx params

r_min = inf;
den_opt = [0,0,0];
num_min = 0;
a1DEL = 10;
a2DEL = 100;
a1LB = 0;
a1UB = 1000;
a2LB = 0;
a2UB = 10000;
s_opt = [];
ts_opt = [];

num = max(H10_2(:, 1));
clc;
for LLL=1:20
    if LLL~=1
        a1LB = den_opt(2)-a1DEL;
        a1UB = den_opt(2)+a1DEL;
        a1DEL = a1DEL/10;

        a2LB = den_opt(1)-a2DEL;
        a2UB = den_opt(1)+a2DEL;
        a2DEL = a2DEL/10;
    end
    for alpha1=a1LB:a1DEL:a1UB
        for alpha2=a2LB:a2DEL:a2UB
            if alpha2 == 2735
                f = 2;
            end
            den = [alpha2, alpha1, 1];
            [s, ts] = sisoctf2dstep(num, den, 0, 4, size(H10_2,1)-1);
            r_new = sum((H10_2(:,1)-s(1:end)).^2);
            if r_new < r_min
                r_min = r_new;
                den_opt = den;
                num_min = num;
                s_opt = s;
                ts_opt = ts;
            end
        end
    end
    disp(LLL);
end
figure;
plot(ts_opt, s_opt);
hold on;
plot(T10_2, H10_2(:,1));
hold off;
title('G21');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');


%%

% G11
K11 = 0.129058;
tau11 = 90;

figure;
plot(0:1:900, K11*(1-exp(-(0:1:900)/tau11)));
hold on;
plot(T10_1, H10_1(:,1));
hold off;
title('G11');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');


% G22
K22 = 0.176909;
tau22 = 95;

figure;
plot(0:1:900, K22*(1-exp(-(0:1:900)/tau22)));
hold on;
plot(T10_2, H10_2(:,2));
hold off;
title('G22');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');

% G12
K12 = 0.109817;
tau1_12 = 120; tau1_12 = 104;
tau2_12 = 35;tau2_12 = 45;
c1 = tau1_12/(tau1_12-tau2_12);
c2 = tau2_12/(tau2_12-tau1_12);

figure;
ys = K12*(1-c1*exp(-(0:1:900)/(tau1_12))-c2*exp(-(0:1:900)/(tau2_12)));
ys = max([ys; zeros(1,length(ys))]); % Avoid zeros in the first variable
plot(0:1:900, ys);
hold on;
plot(T10_1, H10_1(:,2));
hold off;
title('G12');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');

% G21
K21 = 0.0703878;
tau1_21 = 85;
tau2_21 = 35;
c1 = tau1_21/(tau1_21-tau2_21);
c2 = tau2_21/(tau2_21-tau1_21);

figure;
plot(0:1:900, K21*(1-c1*exp(-(0:1:900)/(tau1_21))-c2*exp(-(0:1:900)/(tau2_21))));
hold on;
plot(T10_2, H10_2(:,1));
hold off;
title('G21');
legend('Transfer estimate', 'Simulation', 'Location', 'SouthEast');

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

T2L(param_names, trans_fun, Mat, '../Exam project/Tables/params_sim.tex');


%% Find linear model from param estimates
tol = 1e-4;

% [Ad,Bd,Cd,Dd,sH] = mimoctf2dss(num,den,tau,Ts,Nmax,tol);


