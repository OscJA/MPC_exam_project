function [uk, Ikp1] = PID_controller(ybar, yk, ykm1, I, Kp, Ki, Kd, dt, umax, umin)
%% 
% Choose an optimal flow to minimize the input-constrained optimization
% problem
% Author: Oscar Juul Andersen, s194316
%%

ek = ybar - yk;

P = Kp*ek;
D = -Kd*(yk - ykm1)/dt;
uk = P + I + D;
Ikp1 = I + Ki*ek*dt;

% Ensure that the selected flows are within ids bounds
if uk(1) <= umin(1)
    uk(1) = umin(1);
elseif uk(1) >= umax(1)
    uk(1) = umax(1);
end

if uk(2) <= umin(2)
    uk(2) = umin(2);
elseif uk(2) >= umax(2)
    uk(2) = umax(2);
end

end