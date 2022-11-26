function [Z, R, X, P] = DiscreteKalmanFilter(xkkm1, Pkkm1, yk, A, C, G, R, Q, S, N)

%% Compute innovation
ek = yk - C*xkkm1;
Re = C*Pkkm1*C' + R;

%% Compute the filtered state and the filtered process noise

Kfx = Pkkm1*C'*inv(Re);
Kfw = S*inv(Re);
xkk = xkkm1 + Kfw*ek;
Pkk = Pkkm1 - Kfx*Re*Kfx';
what = Kfx*ek;
Qkk = Q - Kfw*Re*Kfw';

%% State predictions
xkp1k = A*xkk + G*what;
Pkp1k = A*Pkkm1*A' + G*Qkk*G' - A*Kfx*S'*G' - G*S*Kfx'*A';

%% Outputs
X = zeros(length(xkkm1), N);
X(:, 1) = xkp1k;
P = zeros(size(Pkkm1,1), size(Pkkm1,2), N);
P(:, :, 1) = Pkp1k;
for i=2:N
    X(:, i) = A*X(:, i-1);
    P(:, :, i) = A*P(:, :, i-1)*A' + G*Q*G';
end

Z = zeros(length(xkkm1), N);
R = zeros(size(Pkkm1,1), size(Pkkm1,2), N);

for i=1:N
    Z(:, i) = C*X(:, i);
    R(:, :, i) = C*P(:, :, i)*C';
end

end