function [s, converged, guess_history] = ...
    shooting_secant(r, s0, s1, epsilon, k_max) %#ok<FNDEF>

% @parameter r: function handle, residual of the far boundary condition
% @parameter s0, s1: two initial guesses for the unknown initial value
% @parameter epsilon: tolerance on |r(s)|
% @parameter k_max: maximum number of iterations

% @return s: converged value of the unknown initial value
% @return converged: shows if convergence was achieved
% @return guess_history: history of guesses

if nargin < 5; k_max = 100; end
if nargin < 4; epsilon = 1e-8; end

guess_history = zeros(k_max + 2, 1);
guess_history(1:2) = [s0; s1];
r0 = r(s0); r1 = r(s1); k = 2;

while abs(r1) >= epsilon && k < k_max + 2
    s2 = s1 - r1 * (s1 - s0) / (r1 - r0);
    s0 = s1; r0 = r1;
    s1 = s2; r1 = r(s1);
    k = k + 1;
    guess_history(k) = s1;
end

s = s1;
converged = abs(r1) < epsilon;
guess_history = guess_history(1:k);

end