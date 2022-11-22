function [Z, R, X, P] = DiscreteKalmanFilter(xhat, Pkkm1, yk, C, R)

%% Compute innovation
ek = yk - C*xhat;
Re = C*Pkkm1*C' + R;

%% Compute the filtered state and the filtered process noise

Kfx = Pkkm1*C'*inv(Re);
Kfw = S*inv(Re);
xkk = xhat + Kfw*ek;
Pkk = Pkkm1 - Kfx*Re*Kfx';
what = Kfx*ek;
Qkk = Q - Kfw*Re*Kfw';

end