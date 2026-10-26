function [z, C] = diffusion_reaction_fd(C_A0, k, D_AB, L, N) %#ok<FNDEF>
    h = L / N; z = linspace(0, L, N + 1)';
    m2h2 = (k / D_AB) * h^2;

    % Unknowns: C_1, ..., C_N (C_0 = C_A0 is known)
    a = ones(N, 1); b = -(2 + m2h2) * ones(N, 1);
    c = ones(N, 1); d = zeros(N, 1);

    d(1) = -C_A0;   % z = 0: Dirichlet boundary moved to the right-hand side
    a(N) = 2;       % z = L: ghost node C_{N+1} = C_{N-1} (zero gradient)

    C = [C_A0; tdma(a, b, c, d)];
end