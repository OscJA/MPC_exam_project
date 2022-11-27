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
tf_ = 15*60;

%% Discretize & Linearize
hs = ys;
T = sqrt(2*At.*xs)./(ap.*sqrt(rho*g));
% T = (At./ap).*sqrt(hs*4/(2*g));
A = [-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
B = [rho*gam(1) 0;0 rho*gam(2); 0 rho*(1-gam(2)); rho*(1-gam(1)) 0];
C = diag(1./(rho*At));
Cz = C(1:2,:);

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
N = 5;
dim_z = 2;
dim_u = 2;
dim_x = 4;
dim_w = 2;

Cz = C(1:2, :);
G = Ed(1:2, :);

H = zeros(dim_z, dim_u, N);
% H(:, :, 1) = Ed;
Obs_matrix = zeros(dim_z, dim_x, N+1);
Obs_matrix(:, :, 1) = Cz;

Phi_w = zeros(dim_z, dim_w, N);

for i=1:N
    H(:, :, i) = Obs_matrix(:, :, i)*Bd;
    Phi_w(:, :, i) = Obs_matrix(:, :, i)*Ed;
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
R = 40*ones(2, (tf_-t0)/Ts);


% I currently set all W to just be identity matrices
Wz = 100*eye(dim_z); 
Wu = 1*eye(dim_u);
Wdu = 10*eye(dim_u);

Phi_x = reshape(Obs_matrix(:, :, 2:end), N*dim_z, dim_x);
Phi_w = reshape(Phi_w, N*dim_z, dim_w);

%% Calculate matrices used in the optimizer
% Here I calculate all parameters which are constant, to save computing
% power

% First, define some necessary conditions for the problem
% U_bar = repmat([300, 240], 1, N)';
U_bars = [250*ones(1, (tf_-t0)/Ts); 200*ones(1, (tf_-t0)/Ts)];
U_bar = reshape(U_bars(:, 1:N), [], 1);
% Z_bar = repmat([30, 30], 1, N)';
% Z_bars = [repmat([30, 30], 1, 70)'; repmat([40, 30], 1, 230)'];
Z_bars = [35*ones(1, 70), 40*ones(1,160); 49*ones(1,230)];

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

wk = [0; 0];

%% Simulate my optimization

% Initiate matrices
U_realized = zeros(dim_u, (tf_-t0)/Ts);
X = zeros(4, (tf_-t0)/Ts);
X(:, 1) = xs;
Y = zeros(4, (tf_-t0)/Ts);

Rdd = 5*eye(2);
Rvv = 0.1*eye(4);
S = zeros(dim_w, dim_x);

% Find P
P = dlyap(Ad,Ed*Rdd*Ed');

U = [300; 300];
U_realized(:, 1) = U;

for i=1:(tf_-t0)/Ts
    Z_bar = Z_bars(:, i:i+N-1);
    Z_bar = reshape(Z_bar, [], 1);

    xdot = noisyFourTankSystem(0,X(:,i),U(1:dim_u,1),Rdd,d,p);
    X(:, i+1) = X(:, i) + Ts*xdot;
    Y(:, i) = FourTankSystemSensor(X(:, i+1), p)+ Rvv*randn(4,1);

    [xk, wk] = stationaryKalmanFilter(X(:, i), U(1:dim_u), wk, Y(:, i), Ad, Bd, Cd, Ed, P, Rvv, S);
    % U = unconstrainedOptimization(xk, wk, U(1:dim_u), Z_bar, Mdu, H_z, H_u, H_du, g_u, rho_u, WI, gz_mat, rhoz_mat, Phi_x, Phi_w, N);
    U = unconstrainedOptimization(X(:, i), wk, U(1:dim_u), Z_bar, Mdu, H_z, H_u, H_du, g_u, rho_u, WI, gz_mat, rhoz_mat, Phi_x, Phi_w, N);
    U_realized(:, i) = U(1:dim_u);

end

figure;
plot(t0:Ts:(tf_-Ts), U_realized(1, :), '-b');
hold on;
plot(t0:Ts:(tf_-Ts), U_bars(1, :), '--r');
hold off;
legend('Actual', 'Target')
title('Input variable 1 over time');

figure;
plot(t0:Ts:(tf_-Ts), U_realized(2, :), '-b');
hold on;
plot(t0:Ts:(tf_-Ts), U_bars(2, :), '--r');
hold off;
legend('Actual', 'Target')
title('Input variable 2 over time');

figure;
plot(t0:Ts:(tf_-Ts), Y(1, :), '-b');
hold on;
plot(t0:Ts:tf_, Z_bars(1, 1:(tf_-t0)/Ts+1), '--r');
hold off;
title('Tank 1 height');

figure;
plot(t0:Ts:(tf_-Ts), Y(2, :), '-b');
hold on;
plot(t0:Ts:tf_, Z_bars(2, 1:(tf_-t0)/Ts+1), '--r');
hold off;
title('Tank 2 height');

