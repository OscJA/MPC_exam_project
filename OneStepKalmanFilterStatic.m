function [xkk, xkp1k, wkk] = OneStepKalmanFilterStatic(Ad, Bd, Cd, xkkm1, yk, uk, Kfx, Kfw)
%% 
% Perform 1 step prediction of Kalman filter
% Author: Oscar Juul Andersen, s194316
%%

ykkm1 = Cd*xkkm1;
ek = yk-ykkm1;
xkk = xkkm1 + Kfx * ek;

wkk = Kfw * ek;
xkp1k = Ad*xkk + Bd*uk + wkk;


end