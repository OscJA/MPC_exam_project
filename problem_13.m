clear;
close all;

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
tf = 20*60;

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

Rvv = 0.15*eye(4);
Rdd = 0.4*eye(4);

G = eye(6);


X = zeros(6, Nsim);
Y = zeros(4, Nsim);
U = zeros(2, Nsim);

xs = [xs; zeros(2,1)]; % Include the noise variables in xs
X(:, 1) = zeros(6,1);
Qd_chol = chol(Qd);

% PID vars
I = 0;

Kp = 100; %*[1, 1];
Ki = 0.5; %*[1, 1];
Kd = 1; %*[1, 1];
u_ub = 600;
u_lb = 0;

% Set the setpoints
Ybar = ones(2, Nsim);
Ybar(:, 1:50) = 0*Ybar(:, 1:50);
Ybar(:, 51:150) = 10*Ybar(:, 51:150);
Ybar(1, 151:end) = 40*Ybar(1, 151:end)-ys(1);
Ybar(2, 151:end) = 40*Ybar(2, 151:end)-ys(2);

U(:, 1) = us;



for i=1:Nsim-1
    w = Qd_chol*randn(6, 1);
    X(:, i+1) = Ad*X(:, i) + Bd*(U(:, i)-us) + G*w;
    Y(:, i) = Cd*X(:, i) + Rvv*randn(4,1);
    if i == 1
        [U(:, i+1), I] = PID_controller(Ybar(:, i), Y(1:2, i), Y(1:2, i), I, Kp, Ki, Kd, Tstep, u_ub-us, u_lb-us);
    else
        [U(:, i+1), I] = PID_controller(Ybar(:, i), Y(1:2, i), Y(1:2, i-1), I, Kp, Ki, Kd, Tstep, u_ub-us, u_lb-us);
    end
    U(:, i+1) = U(:, i+1) + us;
end
Y(:, Nsim) = Cd*X(:, Nsim) + Rvv*randn(4,1);
Y = Y + ys;
X = X + xs;

%% 
fig = figure('Position', [500 250 850 700]);
subplot(1, 2, 1)
plot(t0:Tstep:(tf-Tstep), Y(1,:));
hold on;
plot(t0:Tstep:(tf-Tstep), Ybar(1, :)+ys(1), '--r');
hold off;
xlabel("Time", 'FontSize', 13.5);
ylabel("Height", 'FontSize', 13.5);
title('Tank 1', 'FontSize', 15);

subplot(1, 2, 2)
plot(t0:Tstep:(tf-Tstep), Y(2,:));
hold on;
plot(t0:Tstep:(tf-Tstep), Ybar(2, :)+ys(2), '--r');
hold off;
xlabel("Time", 'FontSize', 13.5);
ylabel("Height", 'FontSize', 13.5);
title('Tank 2', 'FontSize', 15);


%% Ubar
fig = figure('Position', [500 250 850 700]);
subplot(1, 2, 1)
plot(t0:Tstep:(tf-Tstep), U(1,:));
hold on;
plot(t0:Tstep:(tf-Tstep), 0*ones(Nsim,1), '--g');
plot(t0:Tstep:(tf-Tstep), u_ub*ones(Nsim,1), '--r');
hold off;
ylim([-20, u_ub+20])
xlabel("Time", 'FontSize', 13.5);
ylabel("Height", 'FontSize', 13.5);
title('Tank 1', 'FontSize', 15);

subplot(1, 2, 2)
plot(t0:Tstep:(tf-Tstep), U(2,:));
hold on;
plot(t0:Tstep:(tf-Tstep), 0*ones(Nsim,1), '--g');
plot(t0:Tstep:(tf-Tstep), u_ub*ones(Nsim,1), '--r');
hold off;
ylim([-20, u_ub+20])
xlabel("Time", 'FontSize', 13.5);
ylabel("Height", 'FontSize', 13.5);
title('Tank 2', 'FontSize', 15);

