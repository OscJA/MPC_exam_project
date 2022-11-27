function U = softConstrainedOptimization(xk, wk, ukm1, ulb, uub, A_bar, Z_bar, Mdu, H_z, H_u, H_du, g_u, rho_u, Lambda, Gamma, I0, WI, dU_min, dU_max, Rmin, Rmax, gz_mat, rhoz_mat, Phi_x, Phi_w, N)

% Find b_k
b_k = Phi_x*xk + Phi_w*wk;
%% phi_z parameters
c_k = Z_bar - b_k;
g_z = gz_mat*c_k;
rho_z = 1/2*c_k'*rhoz_mat*c_k;

%% phi_du parameters
g_du = Mdu*ukm1;
rho_du = 1/2*(WI*ukm1)'*(WI*ukm1);

%% Assemble all the methods
H = H_z + H_u + H_du;
g = g_z + g_u + g_du;
rho = rho_z + rho_u + rho_du;

%% Calculate soft output constraints
bl = [dU_min + I0*ukm1; Rmin - b_k; -inf*ones(length(b_k), 1)];
bu = [dU_max + I0*ukm1; inf*ones(length(b_k), 1); Rmax - b_k];

%% Compute the optimal U
U = qpsolver(H,g,ulb,uub,A_bar,bl,bu,ukm1);

end