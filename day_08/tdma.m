function x = tdma(a, b, c, d) %#ok<FNDEF>

% @parameter a: sub-diagonal (a(1) is not used)
% @parameter b: main diagonal
% @parameter c: super-diagonal (c(n) is not used)
% @parameter d: right-hand side
% @return x: solution vector

    n = numel(d);
    c_star = zeros(n, 1); d_star = zeros(n, 1); x = zeros(n, 1);

    % Forward sweep
    c_star(1) = c(1) / b(1);
    d_star(1) = d(1) / b(1);
    for i = 2:n
        denom = b(i) - a(i) * c_star(i-1);
        c_star(i) = c(i) / denom;
        d_star(i) = (d(i) - a(i) * d_star(i-1)) / denom;
    end

    % Back substitution
    x(n) = d_star(n);
    for i = n-1:-1:1
        x(i) = d_star(i) - c_star(i) * x(i+1);
    end
end