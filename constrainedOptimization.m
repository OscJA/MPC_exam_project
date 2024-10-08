function U = constrainedOptimization(xk, wk, ukm1, ulb, uub, Z_bar, Mdu, H_z, H_u, H_du, g_u, rho_u, WI, gz_mat, rhoz_mat, Phi_x, Phi_w, Lambda, dU_min, dU_max)
%% 
% Choose an optimal flow to minimize the input-constrained optimization
% problem
% Author: Oscar Juul Andersen, s194316
%%

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

%% Compute the optimal U
bl = dU_min;
bl(1:length(ukm1)) = bl(1:length(ukm1)) + ukm1;
bu = dU_max;
bu(1:length(ukm1)) = bu(1:length(ukm1)) + ukm1;
A = Lambda;
U = qpsolver(H,g,ulb,uub,A,bl,bu,ukm1);

end