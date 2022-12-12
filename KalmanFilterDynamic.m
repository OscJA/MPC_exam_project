function [Tkk, Xkk, Xkp1k, Ykk, Ykp1k, Pkk, Pkp1k, dkk] = KalmanFilterDynamic(Ad, Bd, Cd, T, X, Y, xs, ys, us, ds, p)

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
Xkp1k = zeros(N-1,4);
dkk = zeros(N-1,nx_diff);
Ykk = zeros(N-1,4);
Ykp1k = zeros(N-1,4);
Tkk = T(2:end);
Pkk = cell(N-1,1);
Pkp1k = cell(N-1,1);

% --------------------------------------------------------------
% Perform 1 iteration of dynamic Kalman filter and save in matrix
% --------------------------------------------------------------
dkkm1 = zeros(nx_diff,1);
xkkm1 = [X(1,:)'-xs;dkkm1];
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
    Qkk = Q-Kfwk*Rek*Kfwk';
    
    xkp1k = Ad*xkk + Bd*(us-us) + wkk;
    ykp1k = FourTankSystemSensor(xkp1k(1:4),p);
    pkp1k = Ad*pkk*Ad' + Qkk - Ad*Kfxk*S' - S*Kfxk'*Ad';
    
    % Save into matrices
    dkk(i,:) = xkp1k(5:(5+nx_diff-1))'+ds';
    Xkk(i,:) = xkk(1:4)+xs;
    Xkp1k(i,:) = xkp1k(1:4)'+xs';
    Ykk(i,:) = ykk'+ys';
    Ykp1k(i,:) = ykp1k'+ys';
    Pkk{i} = pkk;
    Pkp1k{i} = pkp1k;

    % Update for next iteration
    pkkm1 = pkp1k;
    xkkm1 = xkp1k;

end



end