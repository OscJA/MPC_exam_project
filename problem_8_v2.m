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
tf = 60*60;

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

%% Markov parameters
N = 15;
dim_z = 2;
dim_u = 2;
dim_x = 6;
dim_y = 4;
dim_w = 6;

Cz = Cd(1:2, :);
G = eye(dim_x);

H = zeros(dim_z, dim_u, N);
% H(:, :, 1) = Ed;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unconstrained optimization
% First define the target values 
R = 40*ones(2, (tf-t0)/Ts);

u_ub = 500;
u_lb = 0;

% I currently set all W to just be identity matrices
Wz = 100*eye(dim_z); 
Wu = 1*eye(dim_u);
Wdu = 5*eye(dim_u);

Phi_x = zeros(N*dim_z, dim_x);
for i=1:N
    Phi_x(1+(i-1)*dim_z:i*dim_z, :) = Obs_matrix(:, :, i+1);
end
% Phi_x = reshape(Obs_matrix(:, :, 2:end), N*dim_z, dim_x);
Phi_w = reshape(Phi_w, dim_w, N*dim_z)';

%% Calculate matrices used in the optimizer
% Here I calculate all parameters which are constant, to save computing
% power

% First, define some necessary conditions for the problem
U_bars = [250*ones(1, (tf-t0)/Ts); 200*ones(1, (tf-t0)/Ts)]-us;
U_bar = reshape(U_bars(:, 1:N), [], 1);

% Target heights
% Z_bars = [35*ones(1, 70), 40*ones(1,160); 49*ones(1,230)]-ys(1:2);
N_sim = (tf-t0)/Ts;
Z_bars = ones(2, N_sim+N);

% Timestep of the first jump
J1 = 150;
J2 = 350;
J3 = 550;
J4 = 750;
% J5 = 800;

% Hij is the wanted height of container j at jump i
H11 = 38;
H12 = 45;
%H21 = 40;
%H22 = 40;
H21 = 35;
H22 = 35;
H31 = 42;
H32 = 42;
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
%Z_bars(1, J4+1:J5) = H41*Z_bars(1, J4+1:J5);
%Z_bars(2, J4+1:J5) = H42*Z_bars(2, J4+1:J5);

Z_bars(1, J4+1:end) = H41*Z_bars(1, J4+1:end);
Z_bars(2, J4+1:end) = H42*Z_bars(2, J4+1:end);


Z_bars(1, :) = Z_bars(1, :) - ys(1);
Z_bars(2, :) = Z_bars(2, :) - ys(2);



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


%% Simulate my optimization

% Initiate matrices
U = zeros(dim_u, N_sim);
X = zeros(dim_x, N_sim);
Y = zeros(dim_y, N_sim);
Xhat = zeros(dim_x, N_sim);

wk = zeros(dim_w, 1);
xkp1k = zeros(dim_x, 1);

Rdd = 5*eye(2);
Y(:, 1) = Cd*X(:, 1) + Rvv*randn(4,1);

u_lb_vec = repmat([u_lb; u_lb]-us, N, 1)';
u_ub_vec = repmat([u_ub; u_ub]-us, N, 1)';
max_diff_u = 20;
dU_min = repmat([-max_diff_u; -max_diff_u], N, 1);
dU_max = repmat([max_diff_u; max_diff_u], N, 1);

for i=1:(tf-t0)/Ts

    Z_bar = Z_bars(:, i:i+N-1);
    Z_bar = reshape(Z_bar, [], 1);

    w = Qd_chol*randn(6, 1);
    X(:, i+1) = Ad*X(:, i) + Bd*U(:, i) + G*w;
    Y(:, i+1) = Cd*X(:, i+1) + Rvv*randn(4,1);

    [Xhat(:, i), xkp1k, wk] = OneStepKalmanFilterStatic(Ad, Bd, Cd, xkp1k, Y(:, i), U(:, i), Kfx, Kfw);

    % Utmp = unconstrainedOptimization(Xhat(:, i), wk, U(:, i), Z_bar, Mdu, H_z, H_u, H_du, g_u, rho_u, WI, gz_mat, rhoz_mat, Phi_x, Phi_w, N);
    Utmp = constrainedOptimization(Xhat(:, i), wk, U(:, i), u_lb_vec, u_ub_vec, Z_bar, Mdu, H_z, H_u, H_du, g_u, rho_u, WI, gz_mat, rhoz_mat, Phi_x, Phi_w, Lambda, dU_min, dU_max);
    U(:, i+1) = Utmp(1:dim_u);


end


fig = figure('Position', [400 100 900 400]);
subplot(1, 2, 1);
plot(t0:Ts:tf, U(1, :)+us(1), '-b');
hold on;
% plot(t0:Ts:(tf_-Ts), U_bars(1, :)+us(1), '--r');
plot([t0, tf], [u_ub, u_ub], '--r');
plot([t0, tf], [u_lb, u_lb], '--r');
hold off;
xlabel('Time [s]');
ylabel('Flow [cm^3/s]');
%legend('Actual', 'Target')
% ylim([min(U(1, :)+us(1))-20, max(U(1, :)+us(1))+20])
ylim([-20, 520])
title('F_1');

subplot(1, 2, 2);
plot(t0:Ts:tf, U(2, :)+us(2), '-b');
hold on;
plot([t0, tf], [u_ub, u_ub], '--r');
plot([t0, tf], [u_lb, u_lb], '--r');
% plot(t0:Ts:(tf_-Ts), U_bars(2, :)+us(2), '--r');
hold off;
xlabel('Time [s]');
ylabel('Flow [cm^3/s]');
% ylim([min(U(2, :)+us(2))-20, max(U(2, :)+us(2))+20])
ylim([-20, 520])
% legend('Actual', 'Target')
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


Yhat = Cd*Xhat;
fig = figure;
fig.Visible = 'off';
subplot(2,1,1);
plot(t0:Tstep:tf, Y(1,:)+ys(1));
hold on;
plot(t0:Tstep:tf-Tstep, Yhat(1,:)+ys(1));
hold off;
xlabel("Time", 'FontSize', 13.5);
ylabel("Height", 'FontSize', 13.5);
title('Tank 1', 'FontSize', 15);

subplot(2,1,2);
plot(t0:Tstep:tf, Y(2,:)+ys(2));
hold on;
plot(t0:Tstep:tf-Tstep, Yhat(2, :)+ys(2));
hold off;
xlabel("Time", 'FontSize', 13.5);
ylabel("Height", 'FontSize', 13.5);
title('Tank 2', 'FontSize', 15);
legend("Simulated", "Predicted")
