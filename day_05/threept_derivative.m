function dCdt = threept_derivative(C_eq, dt) %#ok<FNDEF>

% @parameter C_eq: equally spaced function values
% @parameter dt: uniform grid spacing

% @return dCdt: derivative estimate at every grid point (forward at the
%   first point, central at interior points, backward at the last point)

Neq = length(C_eq);
dCdt = zeros(1, Neq);

dCdt(1) = (-3*C_eq(1) + 4*C_eq(2) - C_eq(3)) / (2*dt);

for i = 2:Neq-1
    dCdt(i) = (C_eq(i+1) - C_eq(i-1)) / (2*dt);
end

dCdt(Neq) = (3*C_eq(Neq) - 4*C_eq(Neq-1) + C_eq(Neq-2)) / (2*dt);

end
