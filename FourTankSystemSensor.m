function ys = FourTankSystemSensor(xs, p)
% Convert the mass in the container to the height of water in the container

rho = p(12,1); % Density of water [g/cm3]
A = p(5:8,1); % Tank cross sectional areas [cm2]

ys = xs./(rho*A);

end