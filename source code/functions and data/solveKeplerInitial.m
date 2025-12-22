function E0 = solveKeplerInitial(M, e)
    % SOLVEKEPLERINITIAL Newton-Raphson solver for first step
    E0 = M;
    for k = 1:50
        E_new = E0 - (E0 - e*sin(E0) - M) / (1 - e*cos(E0));
        if abs(E_new - E0) < 1e-12
            break;
        end
        E0 = E_new;
    end
end
