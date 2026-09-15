function [dydx, iters] = central_diff_converge(P_fun, x_target, epsilon, k_max) %#ok<FNDEF>
    h = 0.1;
    dydx = (P_fun(x_target + h) - P_fun(x_target - h)) / (2 * h);
    iters = 0;

    while iters < k_max
        iters = iters + 1;
        dydx_prev = dydx;
        h = h / 2;
        dydx = (P_fun(x_target + h) - P_fun(x_target - h)) / (2 * h);

        if abs(dydx - dydx_prev) < epsilon
            break;
        end
    end
end
