function C_eq = piecewise_interp(t_raw, C_raw, t_eq) %#ok<FNDEF>
    Nraw = length(t_raw);
    Neq = length(t_eq);
    C_eq = zeros(1, Neq);

    for i = 1:Neq
        tq = t_eq(i);
        for j = 1:Nraw-1
            if tq >= t_raw(j) && tq <= t_raw(j+1)
                C_eq(i) = C_raw(j) + (C_raw(j+1) - C_raw(j)) * ...
                    (tq - t_raw(j)) / (t_raw(j+1) - t_raw(j));
                break;
            end
        end
    end
end
