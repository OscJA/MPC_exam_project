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
ds = [50; 50];
xs0 = [5000; 5000; 5000; 5000]; % [g] Initial guess on xs
xs = fsolve(@FourTankSystemWrap,xs0,[],us,ds,p); % Steady state masses
ys = FourTankSystemSensor(xs,p); % Steady states heights
p = [p; xs; ys; us; ds];

%% Pack the variables for further calculations
At = [A1; A2; A3; A4];
ap = [a1; a2; a3; a4];
gam = [gamma1; gamma2];
Tstep = 4; % Time step (I think)
t0 = 0;
tf = 15*60;

Sigma = [zeros(4,2); eye(2)];

%% Linearization and discretization

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


%% Create Wiener process noise
Nsim = (tf-t0)/Tstep;
dt = Tstep/Nsim; % Nsim euler-maruyama in simulation (and euler steps in CDEKF)
dW = sqrt(dt) * randn(2, Tstep/dt);

%% Simulation

Delta_d = zeros(2, Nsim);
Delta_d(:, 50) = [50, 50];

Rvv = 0.05*eye(4);
Rdd = 0.4*eye(4);

R = Rvv;
G = eye(4);
G(1,1) = 0;
G(2,2) = 0;
N = 5;



X = zeros(6, Nsim);
Y = zeros(4, Nsim);

U = ones(2, Nsim);
U(:, 1:50) = 300*U(:, 1:50);
U(:, 51:end) = 300*U(:, 51:end);

xs = [xs; zeros(2,1)]; % Include the noise variables in xs
X(:, 1) = zeros(6,1);
% X(:, 1) = xs;
% xdot = A*(xs-xs) + B*(U(:,1)-us) + Sigma*dW(:, 1);
% X(:, 2) = X(:, 1) + Tstep*xdot;
% Y(:, 1) = Cd*X(:, 1);
Qd_chol = chol(Qd);

for i=1:Nsim-1
    w = Qd_chol*randn(6, 1);
    X(:, i+1) = Ad*X(:, i) + Bd*(U(:, i)-us) + eye(6)*w;
    % xdot = Ad*(X(:, i)-xs) + Bd*(U(:, i)-us) + Sigma*dW(:, i);
    % X(:, i+1) = X(:, i) + Tstep*xdot;
    X(5:6, i+1) = X(5:6, i+1) + Delta_d(:, i);
    Y(:, i) = Cd*X(:, i) + Rvv*randn(4,1);
end
Y(:, Nsim) = Cd*X(:, Nsim) + Rvv*randn(4,1);
Y = Y + ys;
X = X + xs;

T = t0:Tstep:tf;
[Tkk, Xkk, Xkp1k, Ykk, Ykp1k, P, dkk] = KalmanFilterStatic(Ad, Bd, Cd, T, X(1:4, :)', Y', xs(1:4), ys, us, ds, p);

fig = figure;
plot(t0:Tstep:(tf-Tstep), Y(3,:));
hold on;
plot((t0+Tstep):Tstep:(tf-Tstep), Ykk(:, 3), 'r');
hold off;

fig = figure;
subplot(2,2,3);
plot(t0:Tstep:(tf-Tstep), Y(1,:));
title('Tank 1')

subplot(2,2,4);
plot(t0:Tstep:(tf-Tstep), Y(2,:));
title('Tank 2')

subplot(2,2,1);
plot(t0:Tstep:(tf-Tstep), Y(3,:));
title('Tank 3')

subplot(2,2,2);
plot(t0:Tstep:(tf-Tstep), Y(4,:));
title('Tank 4')

% figure;
% plot(t0:Tstep:(tf-Tstep), Y(1,:), '-b');
% hold on;
% plot(t0+(N+1)*Tstep:Tstep:tf+(N-1)*Tstep, Ypred(1,:), '-r');
% hold off;
% title('Tank 1')

