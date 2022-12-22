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

%% Linearization and discretization

hs = ys;
T = sqrt(2*At.*xs)./(ap.*sqrt(rho*g));
% T = (At./ap).*sqrt(hs*4/(2*g));
A = [-1/T(1) 0 1/T(3) 0;0 -1/T(2) 0 1/T(4);0 0 -1/T(3) 0;0 0 0 -1/T(4)];
B = [rho*gam(1) 0;0 rho*gam(2); 0 rho*(1-gam(2)); rho*(1-gam(1)) 0];
C = diag(1./(rho*At));
Cz = C(1:2,:);
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
Cd = C;
Ed = M(1:4, 7:8);

%% Simulation

Rvv = 0.05*eye(4);
Rdd = 0.4*eye(4);

Q = Rdd;
R = Rvv;
% G = [zeros(2,2); eye(2)];
G = eye(4);
G(1,1) = 0;
G(2,2) = 0;
N = 5;
P = dlyap(A,G*Q*G');

t0 = 0;
Ts = 4;
tf_ = 15*60;

S = zeros(4, 4);

X = zeros(4, (tf_-t0)/Ts);
Y = zeros(4, (tf_-t0)/Ts);

U = ones(2, (tf_-t0)/Ts);
U(:, 1:50) = 300*U(:, 1:50);
U(:, 51:end) = 350*U(:, 51:end);

% d = Q*randn(4,1);
% d = d(1:2);
% xdot = Ad*X(:, 1) + Bd*U(:, 1) + Ed*d;
X(:, 1) = xs;
Rdd = 5*eye(2);
xdot = noisyFourTankSystem(0,X(:,1),U(:,1),Rdd,d,p);
xdot = FourTankSystem(0,X(:,1),U(:,1),d,p);
X(:, 2) = X(:, 1) + xdot;
Y(:, 1) = FourTankSystemSensor(X(:, 1), p);
% Y(:, 2) = FourTankSystemSensor(X(:, 2));

[Zs, Rs, xnext, Ps] = DiscreteKalmanFilter(X(:, 1), P, Y(:, 1), Ad, Cd, G, R, Q, S, N);
xnext = xnext(:, 1);


% Movie variables
h = figure('units','normalized','outerposition',[0 0 1 1]);
ax = gca;
ax.NextPlot = 'replaceChildren';
loops = (tf_-t0)/Ts;
Mov(loops) = struct('cdata',[],'colormap',[]);
h.Visible = 'off';

T = t0:Ts:tf_;
Ypred = zeros(4, (tf_-t0)/Ts-1);

for i=2:(tf_-t0)/Ts
    xdot = noisyFourTankSystem(0,X(:,i),U(:,i),Rdd,d,p);
    X(:, i+1) = X(:, i) + Ts*xdot;
    Y(:, i) = FourTankSystemSensor(X(:, i+1), p)+ Rvv*randn(4,1);
    Y_kalman = FourTankSystemSensor(X_kalman(:,:), p);    
    Ypred(:, i-1) = Y_kalman(:, end);
    
    plot(T(1:i),Y(1,1:i),"-b");
    hold on
    Tpred = T(i):Ts:T(i)+Ts*N;
    plot(Tpred, Y_kalman(1,:),"-r");
    ylim([30,50])
    xlim([t0,tf_+Ts*N])
    hold off
    drawnow;
    Mov(i-1) = getframe;
    
end

fig = figure;

subplot(2,2,3);
plot(t0:Ts:(tf_-Ts), Y(1,:));
title('Tank 1')

subplot(2,2,4);
plot(t0:Ts:(tf_-Ts), Y(2,:));
title('Tank 2')

subplot(2,2,1);
plot(t0:Ts:(tf_-Ts), Y(3,:));
title('Tank 3')

subplot(2,2,2);
plot(t0:Ts:(tf_-Ts), Y(4,:));
title('Tank 4')

figure;
plot(t0:Ts:(tf_-Ts), Y(1,:), '-b');
hold on;
plot(t0+(N+1)*Ts:Ts:tf_+(N-1)*Ts, Ypred(1,:), '-r');
hold off;
title('Tank 1')

h.Visible = 'on';
% movie(Mov);
movie(Mov,2,1);

