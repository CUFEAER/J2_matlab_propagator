function out = E_Calc(E, t, e, ti, T)
% Solves Kepler's Equation for the Eccentric Anomaly (E):
% $f(E) = t + t_i - \frac{T}{2\pi} (E - e \sin E) = 0$
 out = t + ti - T / (2 * pi) * (E - e * sin(E));
end
