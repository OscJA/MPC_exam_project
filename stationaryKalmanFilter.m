function [xkk, wkk] = stationaryKalmanFilter(xkm1, Ukm1, wkm1, yk, A, B, C, G, P, R, S)

%% Calculate the noise matrices
Rek = C*P*C'+ R;
Kfx = P*C'/(Rek);
Kfw = S/Rek;

%% Filtering
xkkm1 = A*xkm1 + B*Ukm1 + G*wkm1;
ykkm1 = C*xkkm1;
ek = yk-ykkm1;
xkk = xkkm1 + Kfx*ek;
wkk = Kfw*ek;

end