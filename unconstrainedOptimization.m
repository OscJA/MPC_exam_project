function U = unconstrainedOptimization(xk, wk, ukm1, Z_bar, Mdu, H_z, H_u, H_du, g_u, rho_u, WI, gz_mat, rhoz_mat, Phi_x, Phi_w, N)

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
l = zeros(N*length(ukm1), 0);
u = zeros(N*length(ukm1), 0);
bl = zeros(N*length(ukm1), 1);
bu = zeros(N*length(ukm1), 1);
A = zeros(N*length(ukm1), N*length(ukm1));
U = qpsolver(H,g,l,u,A,bl,bu,ukm1);

end