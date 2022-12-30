function [Tkk, Xkk, Xkp1k, Ykk, Ykp1k, P, dkk] = KalmanFilterStatic(Ad, Bd, Cd, T, X, Y, xs, ys, us, ds, Qd, G, R, p)
%% 
% Use the static Kalman filter to predict the states of the system for a
% number of iterations
% Author: Oscar Juul Andersen, s194316
%%

%% Setup matrices
nx = size(Ad,1);
nx_diff = nx-4;


S = zeros(nx,4);
Q = Qd;

N = size(X,1);
Xkk = zeros(N-1,6);
dkk = zeros(N-1,nx_diff);
Xkp1k = zeros(N-1,6);
Ykk = zeros(N-1,4);
Ykp1k = zeros(N-1,4);
Tkk = T(2:end);

%% Calculate static matrices, to save computation time
P = dare(Ad',Cd',G*Q*G',R,G*S);
Re = Cd*P*Cd' + R;
Kfx = P*Cd'*inv(Re);
Kfw = S*inv(Re);

%% Do the Kalman filtering
xkkm1 = X(1,:)'-xs;

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
    Xkk(i,:) = xkk(1:6)+xs;
    Xkp1k(i,:) = xkp1k(1:6)'+xs';
    Ykk(i,:) = (ykk(1:4)+ys)';
    Ykp1k(i,:) = (ykp1k+ys)';
    
    % Update for next iteration
    xkkm1 = xkp1k;
end

end