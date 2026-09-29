function [t, Y_march] = runge_kutta_4(F, Y_0, t_0, t_t, h) %#ok<FNDEF>
    n = round((t_t - t_0) / h); t = t_0:h:t_t;
    Y_march = zeros(numel(F), n + 1);
    Y_march(:, 1) = Y_0;

    f = @(tt, yy) cellfun(@(g) g(tt, yy), F(:));

    for i = 1:n
        Y_current = Y_march(:, i); t_current = t(i);

        k1 = f(t_current, Y_current);
        k2 = f(t_current + h/2, Y_current + h*k1/2);
        k3 = f(t_current + h/2, Y_current + h*k2/2);
        k4 = f(t_current + h, Y_current + h*k3);

        Y_march(:, i + 1) = Y_current + h / 6 * (k1 + 2*k2 + 2*k3 + k4);
    end
end