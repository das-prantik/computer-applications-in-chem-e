function [k, n, R2, iters] = gauss_newton_powerlaw(Cfit, Rfit, k0, n0, tol, k_max) %#ok<FNDEF>

% @parameter Cfit, Rfit: data to be fit, rate = k * C^n
% @parameter k0, n0: initial guesses for k and n
% @parameter tol: convergence tolerance on the parameter update norm
% @parameter k_max: maximum number of Gauss-Newton iterations

% @return k, n: converged parameter estimates
% @return R2: coefficient of determination for the nonlinear fit
% @return iters: number of iterations performed

k = k0; n = n0;
Nfit = length(Cfit);

for iters = 1:k_max
    JTJ = zeros(2,2);
    JTr = zeros(2,1);

    for i = 1:Nfit
        C = Cfit(i);
        predicted = k * C^n;
        residual = Rfit(i) - predicted;

        J1 = -C^n;
        J2 = -k * C^n * log(C);

        JTJ(1,1) = JTJ(1,1) + J1*J1;
        JTJ(1,2) = JTJ(1,2) + J1*J2;
        JTJ(2,1) = JTJ(2,1) + J2*J1;
        JTJ(2,2) = JTJ(2,2) + J2*J2;

        JTr(1) = JTr(1) + J1*residual;
        JTr(2) = JTr(2) + J2*residual;
    end

    detJ = JTJ(1,1)*JTJ(2,2) - JTJ(1,2)*JTJ(2,1);
    if abs(detJ) < 1e-20
        break;
    end

    delta1 = -(JTJ(2,2)*JTr(1) - JTJ(1,2)*JTr(2)) / detJ;
    delta2 = -(-JTJ(2,1)*JTr(1) + JTJ(1,1)*JTr(2)) / detJ;

    k = k + delta1;
    n = n + delta2;

    if sqrt(delta1^2 + delta2^2) < tol
        break;
    end
end

ratebar = mean(Rfit);
predicted_all = k * Cfit.^n;
SSt = sum((Rfit - ratebar).^2);
SSr = sum((Rfit - predicted_all).^2);
R2 = 1 - SSr/SSt;

end
