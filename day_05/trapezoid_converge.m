function [I, iters, m] = ...
    trapezoid_converge(P_fun, a, b, epsilon, k_max) %#ok<FNDEF>
    m = 6;
    I_prev = 0;
    I = 0;
    iters = 0;

    while iters < k_max
        iters = iters + 1;
        h = (b - a) / m;
        x_nodes = linspace(a, b, m + 1);
        y_nodes = P_fun(x_nodes);

        I = (h / 2) * (y_nodes(1) + 2 * sum(y_nodes(2:end-1))...
            + y_nodes(end));

        if iters > 1 && abs(I - I_prev) < epsilon
            break;
        end

        I_prev = I;
        m = m * 2;
    end
end
