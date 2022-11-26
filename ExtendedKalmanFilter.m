function X = ExtendedKalmanFilter(xkk, A, u, d, p, Ts, N)

% %% One step prediction
% xkkp1 = A*xkk;
% 
% % %% Prediction of state covariance
% % Pkp1k = A*Pkk*A' + Q;
% 
% %% Further predictions!
% X = zeros(length(xkk), N+1);
% X(:, 1:2) = [xkk, xkkp1];
% 
% for i=2:N
%     X(:, i+1) = A*X(:, i);
% end
%% Further predictions (nonlinear)!
xkkp1 = xkk + Ts*FourTankSystem(0,xkk,u,d,p);
X = zeros(length(xkk), N+1);
X(:, 1:2) = [xkk, xkkp1];

for i=2:N
    X(:, i+1) = X(:, i) + Ts*FourTankSystem(0,X(:, i),u,d,p);
end

end
