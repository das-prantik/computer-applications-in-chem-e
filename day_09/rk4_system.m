function [x, T] = rk4_system(f, T_0, L, n) %#ok<FNDEF>
    h = L / n;
    x = linspace(0, L, n + 1)';
    T = zeros(n + 1, numel(T_0));
    T(1, :) = T_0(:)';

    for i = 1:n
        T_i = T(i, :)';

        k1 = f(x(i), T_i);
        k2 = f(x(i) + h/2, T_i + h/2 * k1);
        k3 = f(x(i) + h/2, T_i + h/2 * k2);
        k4 = f(x(i) + h, T_i + h * k3);

        T(i + 1, :) = (T_i + h / 6 * (k1 + 2*k2 + 2*k3 + k4))';
    end
end