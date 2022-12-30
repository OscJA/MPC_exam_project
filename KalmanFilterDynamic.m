function [Tkk, Xkk, Xkp1k, Ykk, Ykp1k, Pkk, Pkp1k, dkk] = KalmanFilterDynamic(Ad, Bd, Cd, T, X, Y, xs, ys, us, ds, Q, G, R, p)
%% 
% Use the dynamic Kalman filter to predict the states of the system for a
% number of iterations
% Author: Oscar Juul Andersen, s194316
%%

%% Setup matrices
nx = size(Ad,1);

S = zeros(nx,4);

N = size(X,1);
Xkk = zeros(N-1,6);
Xkp1k = zeros(N-1,6);
dkk = zeros(N-1,nx_diff);
Ykk = zeros(N-1,4);
Ykp1k = zeros(N-1,4);
Tkk = T(2:end);
Pkk = cell(N-1,1);
Pkp1k = cell(N-1,1);
Qkk = Q;


%% Do the Kalman filtering
xkkm1 = X(1,:)'-xs;
pkkm1 = eye(nx);

for i = 1:N-1
    % Filtering
    yk = Y(i+1,:)'-ys;
    ykkm1 = Cd*xkkm1;
    
    ek = yk-ykkm1;
    Rek = Cd*pkkm1*Cd' + R;
    Kfxk = pkkm1*Cd'*inv(Rek);
    Kfwk = S*inv(Rek);
    
    xkk = xkkm1 + Kfxk * ek;
    ykk = FourTankSystemSensor(xkk(1:4),p);
    
    % Make one-step prediction
    wkk = Kfwk * ek;
    pkk = pkkm1 - Kfxk*Rek*Kfxk';
    Qkk = Qkk-Kfwk*Rek*Kfwk';
    
    xkp1k = Ad*xkk + Bd*(us-us) + wkk;
    ykp1k = FourTankSystemSensor(xkp1k(1:4),p);
    pkp1k = Ad*pkk*Ad' + G*Q*G' - Ad*Kfxk*S' - S*Kfxk'*Ad';
    
    % Save into matrices
    dkk(i,:) = xkp1k(5:(5+nx_diff-1))'+ds';
    Xkk(i,:) = xkk+xs;
    Xkp1k(i,:) = xkp1k'+xs';
    Ykk(i,:) = ykk'+ys';
    Ykp1k(i,:) = ykp1k'+ys';
    Pkk{i} = pkk;
    Pkp1k{i} = pkp1k;

    % Update for next iteration
    pkkm1 = pkp1k;
    xkkm1 = xkp1k;

end

end