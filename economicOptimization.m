function U = economicOptimization(xk, wk, ukm1, R, gu, gv, Phi_x, Phi_w, ulb, uub, Lambda, Gamma, I0, dU_min, dU_max, N)
nu = length(ukm1);

% Find b_k
b_k = Phi_x*xk + Phi_w*wk;

% g-vector
g = [gu; gv];


%% Calculate constraints
A = [eye(nu*N), zeros(nu*N);
    -eye(nu*N), zeros(nu*N);
    zeros(nu*N), eye(nu*N);
    Lambda, zeros(nu*N);
    -Lambda, zeros(nu*N);
    Gamma, eye(nu*N)];

b = [ulb;
    -uub;
    zeros(nu*N, 1);
    dU_min+I0*ukm1;
    -dU_max-I0*ukm1;
    R - b_k]';


%% Compute the optimal U
U = linprog(g', -A, -b);

end