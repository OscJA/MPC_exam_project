function [x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit)
% A interface which coonverts a certain input, into input accepted by the
% quadprog solver

%%
n = length(xinit);

lb = l';
ub = u';

A = [A; -A];
b = [bu; -bl];

Aeq = zeros(0, n);
beq = zeros(0, 1);

% Ensure that H is symetric
H = (H+H')/2;

[x, ~, ~, info] = quadprog(H, g', A, b, Aeq, beq, lb, ub, xinit);


end