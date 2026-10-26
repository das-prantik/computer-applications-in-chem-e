function r = shoot_residual(f, T_ho, s, L, n, T_co) %#ok<FNDEF>
    % integrate the IVP with T_c(0) = s and return T_c(L) - T_co
    [~, T] = rk4_system(f, [T_ho; s], L, n);
    r = T(end, 2) - T_co;
end