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
dt = 4; % Time step (I think)

%% Linearization
hs = ys;
T = sqrt(2*At.*xs)./(ap.*sqrt(rho*g));
% T = (At./ap).*sqrt(hs*4/(2*g));
A = [-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
B = [rho*gam(1) 0;0 rho*gam(2); 0 rho*(1-gam(2)); rho*(1-gam(1)) 0];
C = diag(1./(rho*At));
Cz = C(1:2,:);

M = expm([A, B; zeros(2,6)]*dt);
Ad = M(1:4, 1:4);
Bd = M(1:4, 5:6);

%% Print the poles and gains

param_names = ["K"; "\\tau_1"; "\\tau_2"];
trans_fun = ["G_{11}"; "G_{12}"; "G_{21}"; "G_{22}"];
K = num2str([gamma1/A1; (1-gamma2)/(A2*T(3)); (1-gamma1)/(A2*T(4)); gamma2/A2]);
tau_1s = num2str([-1/T(1); -1/T(3); -1/T(4); -1/T(2)]);
tau_2s = [""; -1/T(1); -1/T(2); ""];

% T = table(trans_fun,K,tau_1s,tau_2s);
% table2latex(T, '../Exam project/Tables/T.tex'); % params_sim.tex
% Ttex = table2latex(T, []); % params_sim.tex
precision = 4;

% Mat = [
%     string(round(gamma1/A1, precision)),...
%     string(round((1-gamma2)/(A2*T(3)), precision)),...
%     string(round((1-gamma1)/(A2*T(4)), precision)),...
%     string(round(gamma2/A2, precision));
%     string(round(-1/T(1), precision)),...
%     string(round(-1/T(3), precision)),...
%     string(round(-1/T(4), precision)),...
%     string(round(-1/T(2), precision));
%     "", ...
%     string(round(-1/T(1), precision)),...
%     string(round(-1/T(2), precision)),...
%     ""];

Mat = [
    string(gamma1/A1),...
    string((1-gamma2)/(A2*T(3))),...
    string((1-gamma1)/(A2*T(4))),...
    string(gamma2/A2);
    string(-1/T(1)),...
    string(-1/T(3)),...
    string(-1/T(4)),...
    string(-1/T(2));
    "", ...
    string(-1/T(1)),...
    string(-1/T(2)),...
    ""];

T2L(param_names, trans_fun, Mat, '../Exam project/Tables/params_anal.tex');

%% Discretize
M = expm([A, B; zeros(2,6)]*dt);
Ad1 = M(1:4, 1:4);
Bd1 = M(1:4, 5:6);

E = [
    0, 0;
    0, 0;
    rho, 0;
    0, rho
    ];

M = expm([A, B, E; zeros(4,8)]*dt);
Ad = M(1:4, 1:4);
Bd = M(1:4, 5:6);
Ed = M(1:4, 7:8);
Cd = C;

%% Markov parameters
N = 250;
H = zeros(4, 2, N+1);
H(:, :, 1) = Ed;
% Obs_matrix = zeros(4, 2, N+1);
Obs_matrix = Cd;
for i=1:N
    Obs_matrix = Obs_matrix*Ad;
    H(:, :, i+1) = Obs_matrix*Bd;
end

figure;
plot(reshape(H(1,1,:), 1, []))
title('Markov Parameter 1-1')

figure;
plot(reshape(H(1,2,:), 1, []))
title('Markov Parameter 1-2')

figure;
plot(reshape(H(4,1,:), 1, []))
title('Markov Parameter 4-1')

figure;
plot(reshape(H(4,2,:), 1, []))
title('Markov Parameter 4-2')

