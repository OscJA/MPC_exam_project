function [Tkk, Xkk, Xkp1k, Ykk, Ykp1k, P, dkk] = KalmanFilterStatic(Ad, Bd, Cd, T, X, Y, xs, ys, us, ds, p)

% --------------------------------------------------------------
% Setup random values
% --------------------------------------------------------------
nx = size(Ad,1); %4;
nx_diff = nx-4;
if nx_diff > 2
    ds = [zeros(nx_diff-2,1); ds];
end
%nu = 2;
G = eye(nx);
% Qz = eye(nu); %diag([1,1]);
S = zeros(nx,4);
%Rvv = eye(nx);
Q = eye(nx);
%Rwv = zeros(nx);
% Wz = eye(nu);
% Wu = eye(nu);
% Wdu = zeros(nu); %eye(nu);
R = eye(4);

% --------------------------------------------------------------
% Setup matices
% --------------------------------------------------------------
N = size(X,1);
Xkk = zeros(N-1,4);
dkk = zeros(N-1,nx_diff);
Xkp1k = zeros(N-1,4);
Ykk = zeros(N-1,4);
Ykp1k = zeros(N-1,4);
Tkk = T(2:end);

% --------------------------------------------------------------
% Static, so we can define this out of the loop
% --------------------------------------------------------------
P = dare(Ad',Cd',G*Q*G',R,G*S);
Re = Cd*P*Cd' + R;
Kfx = P*Cd'*inv(Re);
Kfw = S*inv(Re);

% --------------------------------------------------------------
% Perform 1 step prediction of Kalman filter and save in matrix
% --------------------------------------------------------------
dkkm1 = zeros(nx_diff,1);
xkkm1 = [X(1,:)'-xs;dkkm1];

for i = 1:N-1
    % Innovation
    yk = Y(i+1,:)'-ys;
    ykkm1 = Cd*xkkm1;
    ek = yk-ykkm1;
    xkk = xkkm1 + Kfx * ek;
    ykk = FourTankSystemSensor(xkk(1:4),p);
    
    % Make one-step prediction
    wkk = Kfw * ek;
    xkp1k = Ad*xkk + Bd*(us-us) + wkk;
    ykp1k = FourTankSystemSensor(xkp1k(1:4),p);
    
    % Save into matrices
    dkk(i,:) = xkp1k(5:(5+nx_diff-1))'+ds';
    Xkk(i,:) = xkk(1:4)+xs;
    Xkp1k(i,:) = xkp1k(1:4)'+xs';
    Ykk(i,:) = (ykk(1:4)+ys)';
    Ykp1k(i,:) = (ykp1k+ys)';
    
    % Update for next iteration
    xkkm1 = xkp1k;
end



end