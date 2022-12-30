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
t0 = 0;
Ts = dt;
Tstep = Ts;
tf = 45*60;

%% Linearization and discretization
Sigma = [zeros(4,2); eye(2)];

hs = ys;
T = sqrt(2*At.*xs)./(ap.*sqrt(rho*g));
A = [-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
B = [rho*gam(1) 0;0 rho*gam(2); 0 rho*(1-gam(2)); rho*(1-gam(1)) 0];
C = diag(1./(rho*At));
Cz = C(1:2,:);

% Include the disturbance variables in the matrices
A = [[A; zeros(2,4)], [zeros(2,2); [rho, 0; 0, rho]; zeros(2,2)]];
B = [B; zeros(2,2)];
C = [C, zeros(4,2)];

M = expm([A, B; zeros(2,8)]*Tstep);
Ad = M(1:6, 1:6);
Bd = M(1:6, 7:8);
Cd = C;

RES = expm([-A, Sigma*Sigma'; zeros(size(A)), A']*Tstep);
Qd = Ad*RES(1:6, 7:end);
Qd_chol = chol(Qd);

%% Calculate matrices which we need to use in the simulation and controllers
N = 25;
dim_z = 2;
dim_u = 2;
dim_x = 6;
dim_y = 4;
dim_w = 6;

Cz = Cd(1:2, :);
G = eye(dim_x);

H = zeros(dim_z, dim_u, N);
Obs_matrix = zeros(dim_z, dim_x, N+1);
Obs_matrix(:, :, 1) = Cz;

Phi_w = zeros(dim_z, dim_w, N);

for i=1:N
    H(:, :, i) = Obs_matrix(:, :, i)*Bd;
    Phi_w(:, :, i) = Obs_matrix(:, :, i)*G;
    Obs_matrix(:, :, i+1) = Obs_matrix(:, :, i)*Ad;
end
Gamma = zeros(dim_z*N, dim_u*N);

for i=1:N
    for j=1:i
        Gamma(1+(i-1)*dim_z:i*dim_z,1+(j-1)*dim_u:j*dim_u) = H(:, :, i+1-j);
    end
end

Lambda = eye(N*dim_u);
for i=1:N-1
    Lambda(1+i*dim_u:(i+1)*dim_u, 1+(i-1)*dim_u:i*dim_u) = -eye(dim_u);
end

Phi_x = zeros(N*dim_z, dim_x);
for i=1:N
    Phi_x(1+(i-1)*dim_z:i*dim_z, :) = Obs_matrix(:, :, i+1);
end
Phi_w = reshape(Phi_w, dim_w, N*dim_z)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up weights and limits of the problem
% Minimum and maximum flows
u_ub = 500;
u_lb = 0;

u_lb_vec = repmat([u_lb; u_lb]-us, 3*N, 1)';
u_ub_vec = [repmat([u_ub; u_ub]-us, N, 1); inf*ones(2*dim_z*N, 1)];
max_diff_u = 20;
dU_min = repmat([-max_diff_u; -max_diff_u], N, 1);
dU_max = repmat([max_diff_u; max_diff_u], N, 1);

% Weight matrices
Wz = 10*eye(dim_z); 
Wu = 0.1*eye(dim_u);
Wdu = 10*eye(dim_u);

w1 = 25;
w2 = 25;
Ws1 = w1*eye(dim_z);
Ws2 = w2*eye(dim_z);
Wt1 = w1*eye(dim_z);
Wt2 = w2*eye(dim_z);


%% Calculate matrices used in the optimizer
% Here I calculate all parameters which are constant, to save computing
% power

% First, define some necessary conditions for the problem
U_bars = zeros(2,N);
U_bar = reshape(U_bars(:, 1:N), [], 1);

% Target heights
N_sim = (tf-t0)/Ts;
Z_bars = ones(2, N_sim+N+3*N_sim);

%% Determine the setpoints
% Timestep of the first jump
J1 = 150;
J2 = 350;
J3 = 550;
J4 = 750;

% Hij is the wanted height of container j at jump i
H11 = 45;
H12 = 60;
H21 = 28;
H22 = 35;
H31 = 55;
H32 = 77;
H41 = 48;
H42 = 48;

Z_bars(1, 1:J1) = ys(1)*Z_bars(1, 1:J1);
Z_bars(2, 1:J1) = ys(2)*Z_bars(2, 1:J1);
Z_bars(1, J1+1:J2) = H11*Z_bars(1, J1+1:J2);
Z_bars(2, J1+1:J2) = H12*Z_bars(2, J1+1:J2);
Z_bars(1, J2+1:J3) = H21*Z_bars(1, J2+1:J3);
Z_bars(2, J2+1:J3) = H22*Z_bars(2, J2+1:J3);
Z_bars(1, J3+1:J4) = H31*Z_bars(1, J3+1:J4);
Z_bars(2, J3+1:J4) = H32*Z_bars(2, J3+1:J4);

Z_bars(1, J4+1:end) = H41*Z_bars(1, J4+1:end);
Z_bars(2, J4+1:end) = H42*Z_bars(2, J4+1:end);

Z_bars(1, :) = Z_bars(1, :) - ys(1);
Z_bars(2, :) = Z_bars(2, :) - ys(2);


%% Calculate optimization parameters
% phi_z parameters
Wz_bar = kron(eye(N), Wz);
H_z = (Wz_bar*Gamma)'*(Wz_bar*Gamma);
gz_mat = -(Wz_bar*Gamma)'*Wz_bar;
rhoz_mat = Wz_bar'*Wz_bar;

% phi_u parameters
Wu_bar = kron(eye(N), Wu);
H_u = Wu_bar'*Wu_bar;
g_u = -Wu_bar'*Wu_bar*U_bar;
rho_u = 1/2*U_bar'*H_u*U_bar;

% phi_du variables
I0 = [eye(dim_u); zeros(dim_u*(N-1), dim_u)];
Wdu_bar = kron(eye(N), Wdu);
H_du = (Wdu_bar*Lambda)'*(Wdu_bar*Lambda);
Mdu = -(Wdu_bar*Lambda)'*Wdu_bar*I0;
WI = Wdu_bar*I0;


%% Calculate matrices used in the Kalman filter

S = zeros(dim_y, dim_x)';
Rvv = 0.15*eye(4);
P = dare(Ad',Cd',G*Qd*G',Rvv,G*S);
Re = Cd*P*Cd' + Rvv;
Kfx = P*Cd'*inv(Re);
Kfw = S*inv(Re);


%% Define soft constraint variables
A_bar = [
    Lambda, zeros(N*dim_u, 2*N*dim_u);
    Gamma, eye(N*dim_u), zeros(N*dim_u, N*dim_u);
    Gamma, zeros(N*dim_u, N*dim_u), -eye(N*dim_u)
    ];


Ws2_bar = kron(eye(N), Ws2);
Ws1_bar = kron(eye(N), Ws1);
Wt2_bar = kron(eye(N), Wt2);
Wt1_bar = kron(eye(N), Wt1);

H_s = Ws2_bar'*Ws2_bar;
H_t = Wt2_bar'*Wt2_bar;
g_s = Ws1_bar*ones(dim_z*N, 1);
g_t = Wt1_bar*ones(dim_z*N, 1);

%% Caluate noise realizations
W = Qd_chol*randn(6, N_sim);
V = Rvv*randn(4, N_sim);

%% Initiate matrices
U = zeros(dim_u, N_sim);
X = zeros(dim_x, N_sim);
Y = zeros(dim_y, N_sim);
Xhat = zeros(dim_x, N_sim);

xkp1k = zeros(dim_x, 1);

Y(:, 1) = Cd*X(:, 1) + Rvv*randn(4,1);


%% UNCSONTRAINED OPTIMIZATION

for i=1:(tf-t0)/Ts

    Z_bar = Z_bars(:, i:i+N-1);
    Z_bar = reshape(Z_bar, [], 1);

    w = Qd_chol*randn(6, 1);
    X(:, i+1) = Ad*X(:, i) + Bd*U(:, i) + G*W(:, i);
    Y(:, i+1) = Cd*X(:, i+1) + V(:, i);

    [Xhat(:, i), xkp1k, wk] = OneStepKalmanFilterStatic(Ad, Bd, Cd, xkp1k, Y(:, i), U(:, i), Kfx, Kfw);

    Utmp = unconstrainedOptimization(Xhat(:, i), wk, U(:, i), Z_bar, Mdu, H_z, H_u, H_du, g_u, rho_u, WI, gz_mat, rhoz_mat, Phi_x, Phi_w, N);
    U(:, i+1) = Utmp(1:dim_u);

end
% Save the flows
Uu1 = U(1, :)+us(1);
Uu2 = U(2, :)+us(2);
hu1 = Y(1, :)+ys(1);
hu2 = Y(2, :)+ys(2);

%% Plot the controlled flow and heights
fig = figure('Position', [400 100 900 400]);
subplot(1, 2, 1);
plot(t0:Ts:tf, U(1, :)+us(1), '-b');
hold on;
plot([t0, tf], [u_ub, u_ub], '--r');
plot([t0, tf], [u_lb, u_lb], '--r');
hold off;
xlabel('Time [s]');
ylabel('Flow [cm^3/s]');
ylim([-20, 520])
title('F_1');

subplot(1, 2, 2);
plot(t0:Ts:tf, U(2, :)+us(2), '-b');
hold on;
plot([t0, tf], [u_ub, u_ub], '--r');
plot([t0, tf], [u_lb, u_lb], '--r');
hold off;
xlabel('Time [s]');
ylabel('Flow [cm^3/s]');
ylim([-20, 520])
title('F_2');
saveas(fig, '../Exam project/Figures/unconstrained_flow.png')


fig = figure('Position', [400 100 900 400]);
subplot(1, 2, 1)
plot(t0:Ts:tf, Y(1, :)+ys(1), '-b');
hold on;
plot(t0:Ts:tf, Z_bars(1, 1:(tf-t0)/Ts+1)+ys(1), '--r');
hold off;
xlabel('Time [s]');
ylabel('Height [cm]');
title('Tank 1');

subplot(1, 2, 2)
plot(t0:Ts:tf, Y(2,:)+ys(2), '-b');
hold on;
plot(t0:Ts:tf, Z_bars(2, 1:(tf-t0)/Ts+1)+ys(2), '--r');
hold off;
xlabel('Time [s]');
ylabel('Height [cm]');
title('Tank 2');
saveas(fig, '../Exam project/Figures/unconstrained_heights.png')


%% INPUT CONSTRAINED OPTIMIZATION
% Reset matrices and vectors before next simulation
U = zeros(dim_u, N_sim);
X = zeros(dim_x, N_sim);
Y = zeros(dim_y, N_sim);
Xhat = zeros(dim_x, N_sim);
xkp1k = zeros(dim_x, 1);

for i=1:(tf-t0)/Ts

    Z_bar = Z_bars(:, i:i+N-1);
    Z_bar = reshape(Z_bar, [], 1);

    X(:, i+1) = Ad*X(:, i) + Bd*U(:, i) + G*W(:, i);
    Y(:, i+1) = Cd*X(:, i+1) + V(:, i);

    [Xhat(:, i), xkp1k, wk] = OneStepKalmanFilterStatic(Ad, Bd, Cd, xkp1k, Y(:, i), U(:, i), Kfx, Kfw);

    Utmp = constrainedOptimization(Xhat(:, i), wk, U(:, i), u_lb_vec(1:2*N), u_ub_vec(1:2*N), Z_bar, Mdu, H_z, H_u, H_du, g_u, rho_u, WI, gz_mat, rhoz_mat, Phi_x, Phi_w, Lambda, dU_min, dU_max);
    U(:, i+1) = Utmp(1:dim_u);

end

% Save the flows
Ui1 = U(1, :)+us(1);
Ui2 = U(2, :)+us(2);
hi1 = Y(1, :)+ys(1);
hi2 = Y(2, :)+ys(2);

%% Plot the controlled flow and heights
fig = figure('Position', [400 100 900 400]);
subplot(1, 2, 1);
plot(t0:Ts:tf, U(1, :)+us(1), '-b');
hold on;
plot([t0, tf], [u_ub, u_ub], '--r');
plot([t0, tf], [u_lb, u_lb], '--r');
hold off;
xlabel('Time [s]');
ylabel('Flow [cm^3/s]');
ylim([-20, 520])
title('F_1');

subplot(1, 2, 2);
plot(t0:Ts:tf, U(2, :)+us(2), '-b');
hold on;
plot([t0, tf], [u_ub, u_ub], '--r');
plot([t0, tf], [u_lb, u_lb], '--r');
hold off;
xlabel('Time [s]');
ylabel('Flow [cm^3/s]');
ylim([-20, 520])
title('F_2');
saveas(fig, '../Exam project/Figures/input_constrained_flow.png')


fig = figure('Position', [400 100 900 400]);
subplot(1, 2, 1)
plot(t0:Ts:tf, Y(1, :)+ys(1), '-b');
hold on;
plot(t0:Ts:tf, Z_bars(1, 1:(tf-t0)/Ts+1)+ys(1), '--r');
hold off;
xlabel('Time [s]');
ylabel('Height [cm]');
title('Tank 1');

subplot(1, 2, 2)
plot(t0:Ts:tf, Y(2,:)+ys(2), '-b');
hold on;
plot(t0:Ts:tf, Z_bars(2, 1:(tf-t0)/Ts+1)+ys(2), '--r');
hold off;
xlabel('Time [s]');
ylabel('Height [cm]');
title('Tank 2');
saveas(fig, '../Exam project/Figures/input_constrained_heights.png')

%% SOFT CONSTRAINED OPTIMIZATION
% Reset matrices and vectors before next simulation
U = zeros(dim_u, N_sim);
X = zeros(dim_x, N_sim);
Y = zeros(dim_y, N_sim);
Xhat = zeros(dim_x, N_sim);
xkp1k = zeros(dim_x, 1);

% Soft constraint max and min heights
min_h = [20; 30];
max_h = [55; 77];

Rmin = repmat(min_h-ys(1:2), [N, 1]);
Rmax = repmat(max_h-ys(1:2), [N, 1]);

for i=1:N_sim

    Z_bar = Z_bars(:, i:i+N-1);
    Z_bar = reshape(Z_bar, [], 1);


    X(:, i+1) = Ad*X(:, i) + Bd*U(:, i) + G*W(:, i);
    Y(:, i+1) = Cd*X(:, i+1) + V(:, i);

    [Xhat(:, i), xkp1k, wk] = OneStepKalmanFilterStatic(Ad, Bd, Cd, xkp1k, Y(:, i), U(:, i), Kfx, Kfw);

    Utmp = softConstrainedOptimization(Xhat(:, i), wk, U(:, i), u_lb_vec, u_ub_vec, A_bar, Z_bar, Mdu, H_z, H_u, H_du, H_s, H_t, g_u, g_s, g_t, rho_u, I0, WI, dU_min, dU_max, Rmin, Rmax, gz_mat, rhoz_mat, Phi_x, Phi_w, N);
    U(:, i+1) = Utmp(1:dim_u);

end

% Save the flows and heights
Us1 = U(1, :)+us(1);
Us2 = U(2, :)+us(2);
hs1 = Y(1, :)+ys(1);
hs2 = Y(2, :)+ys(2);


%% Plot the controlled flow and heights

fig = figure('Position', [400 100 900 400]);
subplot(1, 2, 1);
plot(t0:Ts:tf, U(1, :)+us(1), '-b');
hold on;
plot([t0, tf], [u_ub, u_ub], '--r');
plot([t0, tf], [u_lb, u_lb], '--r');
hold off;
xlabel('Time [s]');
ylabel('Flow [cm^3/s]');
ylim([-20, 520])
title('F_1');

subplot(1, 2, 2);
plot(t0:Ts:tf, U(2, :)+us(2), '-b');
hold on;
plot([t0, tf], [u_ub, u_ub], '--r');
plot([t0, tf], [u_lb, u_lb], '--r');
hold off;
xlabel('Time [s]');
ylabel('Flow [cm^3/s]');
ylim([-20, 520])
title('F_2');
saveas(fig, '../Exam project/Figures/soft_constrained_flow.png')


fig = figure('Position', [400 100 900 400]);
subplot(1, 2, 1)
plot(t0:Ts:tf, Y(1, :)+ys(1), '-b');
hold on;
plot(t0:Ts:tf, Z_bars(1, 1:(tf-t0)/Ts+1)+ys(1), '--r');
hold off;
xlabel('Time [s]');
ylabel('Height [cm]');
title('Tank 1');

subplot(1, 2, 2)
plot(t0:Ts:tf, Y(2,:)+ys(2), '-b');
hold on;
plot(t0:Ts:tf, Z_bars(2, 1:(tf-t0)/Ts+1)+ys(2), '--r');
hold off;
xlabel('Time [s]');
ylabel('Height [cm]');
title('Tank 2');
saveas(fig, '../Exam project/Figures/soft_constrained_heights.png')


%% Plot everything together
fig = figure('Position', [400 100 900 400]);
subplot(1, 2, 1);
hold on;
plot([t0, tf], [u_ub, u_ub], '--black');
plot([t0, tf], [u_lb, u_lb], '--black');
plot(t0:Ts:tf, Uu1, '-b');
plot(t0:Ts:tf, Ui1, '--g');
plot(t0:Ts:tf, Us1, '-r');
hold off;
xlabel('Time [s]');
ylabel('Flow [cm^3/s]');
ylim([min([0, min(Uu1), min(Ui1), min(Us1)])-20, max([500, max(Uu1), max(Ui1), max(Us1)])+20])
title('F_1');

subplot(1, 2, 2);
hold on;
plot([t0, tf], [u_ub, u_ub], '--black');
plot([t0, tf], [u_lb, u_lb], '--black');
plot(t0:Ts:tf, Uu2, '-b');
plot(t0:Ts:tf, Ui2, '--g');
plot(t0:Ts:tf, Us2, '-r');
hold off;
xlabel('Time [s]');
ylabel('Flow [cm^3/s]');
ylim([min([0, min(Uu2), min(Ui2), min(Us2)])-20, max([500, max(Uu2), max(Ui2), max(Us2)])+20])
title('F_2');
legend("", "", "Unconstrained", "Input constrained", "Soft constrained", "Location", "SouthEast");
saveas(fig, '../Exam project/Figures/all_flows.png')


fig = figure('Position', [400 100 900 400]);
subplot(1, 2, 1)
hold on;
plot(t0:Ts:tf, Z_bars(1, 1:(tf-t0)/Ts+1)+ys(1), '--black');
plot([t0, tf], [max_h(1), max_h(1)], 'Linestyle', '--', 'Color', '#77AC30');
plot([t0, tf], [min_h(1), min_h(1)], 'Linestyle', '--', 'Color', '#77AC30');
plot(t0:Ts:tf, hu1, '-b');
plot(t0:Ts:tf, hi1, '--g');
plot(t0:Ts:tf, hs1, '-r');
hold off;
xlabel('Time [s]');
ylabel('Height [cm]');
title('Tank 1');

subplot(1, 2, 2)
hold on;
plot(t0:Ts:tf, Z_bars(2, 1:(tf-t0)/Ts+1)+ys(2), '--black');
plot([t0, tf], [max_h(2), max_h(2)], 'Linestyle', '--', 'Color', '#77AC30');
plot([t0, tf], [min_h(2), min_h(2)], 'Linestyle', '--', 'Color', '#77AC30');
plot(t0:Ts:tf, hu2, '-b');
plot(t0:Ts:tf, hi2, '--g');
plot(t0:Ts:tf, hs2, '-r');
hold off;
xlabel('Time [s]');
ylabel('Height [cm]');
title('Tank 2');
legend("", "", "", "Unconstrained", "Input constrained", "Soft constrained", "Location", [0.45,0.18,0.15,0.1]);

saveas(fig, '../Exam project/Figures/all_heights.png')


%% PID-Controllers
% P-Controller
% PID params
Kp = 50;
Ki = 0.1;
Kd = 150;

% Reset matrices and vectors before next simulation
U = zeros(dim_u, N_sim);
X = zeros(dim_x, N_sim);
Y = zeros(dim_y, N_sim);
I = 0;


for i=1:N_sim
    w = Qd_chol*randn(6, 1);
    X(:, i+1) = Ad*X(:, i) + Bd*(U(:, i)-us) + G*W(:, i);
    Y(:, i+1) = Cd*X(:, i) + V(:, i);
    if i == 1
        [U(:, i+1), I] = PID_controller(Z_bars(:, i), Y(1:2, i), Y(1:2, i), I, Kp, 0, 0, Tstep, u_ub-us, u_lb-us);
    else
        [U(:, i+1), I] = PID_controller(Z_bars(:, i), Y(1:2, i), Y(1:2, i-1), I, Kp, 0, 0, Tstep, u_ub-us, u_lb-us);
    end
    U(:, i+1) = U(:, i+1) + us;
end
% Save the flows and heights
Up1 = U(1, :);
Up2 = U(2, :);
hp1 = Y(1, :)+ys(1);
hp2 = Y(2, :)+ys(2);


%% PI-Controller
% Reset matrices and vectors before next simulation
U = zeros(dim_u, N_sim);
X = zeros(dim_x, N_sim);
Y = zeros(dim_y, N_sim);
I = 0;

for i=1:N_sim
    w = Qd_chol*randn(6, 1);
    X(:, i+1) = Ad*X(:, i) + Bd*(U(:, i)-us) + G*W(:, i);
    Y(:, i+1) = Cd*X(:, i) + V(:, i);
    if i == 1
        [U(:, i+1), I] = PID_controller(Z_bars(:, i), Y(1:2, i), Y(1:2, i), I, Kp, Ki, 0, Tstep, u_ub-us, u_lb-us);
    else
        [U(:, i+1), I] = PID_controller(Z_bars(:, i), Y(1:2, i), Y(1:2, i-1), I, Kp, Ki, 0, Tstep, u_ub-us, u_lb-us);
    end
    U(:, i+1) = U(:, i+1) + us;
end
% Save the flows and heights
Upi1 = U(1, :);
Upi2 = U(2, :);
hpi1 = Y(1, :)+ys(1);
hpi2 = Y(2, :)+ys(2);


%% PID-Controller
% Reset matrices and vectors before next simulation
U = zeros(dim_u, N_sim);
X = zeros(dim_x, N_sim);
Y = zeros(dim_y, N_sim);
I = 0;


for i=1:N_sim
    w = Qd_chol*randn(6, 1);
    X(:, i+1) = Ad*X(:, i) + Bd*(U(:, i)-us) + G*W(:, i);
    Y(:, i+1) = Cd*X(:, i) + V(:, i);
    if i == 1
        [U(:, i+1), I] = PID_controller(Z_bars(:, i), Y(1:2, i), Y(1:2, i), I, Kp, Ki, Kd, Tstep, u_ub-us, u_lb-us);
    else
        [U(:, i+1), I] = PID_controller(Z_bars(:, i), Y(1:2, i), Y(1:2, i-1), I, Kp, Ki, Kd, Tstep, u_ub-us, u_lb-us);
    end
    U(:, i+1) = U(:, i+1) + us;
end
% Save the flows and heights
Upid1 = U(1, :);
Upid2 = U(2, :);
hpid1 = Y(1, :)+ys(1);
hpid2 = Y(2, :)+ys(2);

%% PID PLOTS
fig = figure('Position', [400 100 900 400]);
subplot(1, 2, 1);
hold on;
plot([t0, tf], [u_ub, u_ub], '--black');
plot([t0, tf], [u_lb, u_lb], '--black');
plot(t0:Ts:tf, Up1, '-b');
plot(t0:Ts:tf, Upi1, '-g');
plot(t0:Ts:tf, Upid1, '-r');
hold off;
xlabel('Time [s]');
ylabel('Flow [cm^3/s]');
ylim([min([0, min(Uu1), min(Ui1), min(Us1)])-20, max([500, max(Uu1), max(Ui1), max(Us1)])+20])
title('F_1');

subplot(1, 2, 2);
hold on;
plot([t0, tf], [u_ub, u_ub], '--black');
plot([t0, tf], [u_lb, u_lb], '--black');
plot(t0:Ts:tf, Up2, '-b');
plot(t0:Ts:tf, Upi2, '-g');
plot(t0:Ts:tf, Upid2, '-r');
hold off;
xlabel('Time [s]');
ylabel('Flow [cm^3/s]');
ylim([min([0, min(Uu2), min(Ui2), min(Us2)])-20, max([500, max(Uu2), max(Ui2), max(Us2)])+20])
title('F_2');
legend("", "", "P", "PI", "PID", "Location", "SouthEast");
saveas(fig, '../Exam project/Figures/PID_all_flows.png')


fig = figure('Position', [400 100 900 400]);
subplot(1, 2, 1)
hold on;
plot(t0:Ts:tf, Z_bars(1, 1:(tf-t0)/Ts+1)+ys(1), '--black');
plot(t0:Ts:tf, hp1, '-b');
plot(t0:Ts:tf, hpi1, '--g');
plot(t0:Ts:tf, hpid1, '-r');
hold off;
xlabel('Time [s]');
ylabel('Height [cm]');
title('Tank 1');

subplot(1, 2, 2)
hold on;
plot(t0:Ts:tf, Z_bars(2, 1:(tf-t0)/Ts+1)+ys(2), '--black');
plot(t0:Ts:tf, hp2, '-b');
plot(t0:Ts:tf, hpi2, '-g');
plot(t0:Ts:tf, hpid2, '-r');
hold off;
xlabel('Time [s]');
ylabel('Height [cm]');
title('Tank 2');
legend("", "P", "PI", "PID", "Location", [0.44,0.18,0.15,0.1]);
saveas(fig, '../Exam project/Figures/PID_all_heights.png')
