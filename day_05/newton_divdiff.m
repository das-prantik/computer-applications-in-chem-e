function [P_sym, P_fun, coeffs] = ...
    newton_divdiff(x_data, y_data) %#ok<FNDEF>
    n = length(x_data);
    F = zeros(n, n);
    F(:, 1) = y_data(:);

    for j = 2:n
        for i = 1:n-j+1
            F(i, j) = (F(i+1, j-1) - F(i, j-1)) / (x_data(i+j-1) - x_data(i));
        end
    end

    coeffs = F(1, :);

    syms x
    P_sym = coeffs(1);
    term = sym(1);

    for i = 2:n
        term = term * (x - x_data(i-1));
        P_sym = P_sym + coeffs(i) * term;
    end

    P_sym = simplify(P_sym);
    P_fun = matlabFunction(P_sym);
end
