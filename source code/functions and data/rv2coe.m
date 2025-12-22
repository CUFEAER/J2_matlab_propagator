function [a, e, i, W0, w0, nu,M] = rv2coe(r, v)
   
    
    mu = 398600; % Earth gravitational parameter (km^3/s^2)
    eps = 1.e-10; % Tolerance for "zero" checks

    % 1. Magnitudes
    rmag = norm(r);
    vmag = norm(v);
    
    % 2. Radial Velocity (vr)
    vr = dot(r, v) / rmag;
    
    % 3. Specific Angular Momentum (h)
    h_vec = cross(r, v);
    hmag = norm(h_vec);
    
    % 4. Inclination (i)
    i = acos(h_vec(3) / hmag);
    
    % 5. Node Line (N)
    N_vec = cross([0 0 1], h_vec);
    Nmag = norm(N_vec);
    
    % 6. RAAN (W0)
    if Nmag ~= 0
        W0 = acos(N_vec(1) / Nmag);
        if N_vec(2) < 0
            W0 = 2*pi - W0; % Quadrant check
        end
    else
        W0 = 0; % Equatorial orbit
        
    end
    
    % 7. Eccentricity Vector (e)
    % E = 1/mu * ((v^2 - mu/r)*R - (r*vr)*V)
    term1 = (vmag^2 - mu/rmag) * r;
    term2 = (rmag * vr) * v;
    e_vec = (1/mu) * (term1 - term2);
    
    % 8. Eccentricity Magnitude
    e = norm(e_vec);
    
    % Argument of Perigee (w0)
    if Nmag ~= 0
        if e > eps
            w0 = acos(dot(N_vec, e_vec) / (Nmag * e));
            if e_vec(3) < 0
                w0 = 2*pi - w0; % Quadrant check
            end
        else
            w0 = 0; % Circular orbit
        end
    else
        w0 = 0; % Equatorial orbit
    end
    
    % 10. True Anomaly (nu)
    if e > eps
        nu = acos(dot(e_vec, r) / (e * rmag));
        if vr < 0
            nu = 2*pi - nu; % Quadrant check
        end
    else
        % Circular orbit case 
        cp = cross(N_vec, r);
        if cp(3) >= 0
            nu = acos(dot(N_vec, r) / (Nmag * rmag));
        else
            nu = 2*pi - acos(dot(N_vec, r) / (Nmag * rmag));
        end
    end
    
    % -- Semi-Major Axis (a) --
    % Derived from Vis-Viva: v^2/2 - mu/r = -mu/2a
    specific_energy = (vmag^2 / 2) - (mu / rmag);
    a = -mu / (2 * specific_energy);
    
   

if e < 1
        % Elliptical Orbit: Calculate Eccentric Anomaly (E) first
        % tan(E/2) = sqrt((1-e)/(1+e)) * tan(nu/2)
        E = 2 * atan(sqrt((1-e)/(1+e)) * tan(nu/2));
        
        % Kepler's Equation: M = E - e*sin(E)
        M = E - e * sin(E);
        
        % Ensure M is in 0..2pi range
        if M < 0, M = M + 2*pi; end
    else
        % Parabolic/Hyperbolic case 
        M = 0; 
end
mu = 398600;
energy = norm(v)^2/2 - mu/norm(r);
if energy >= 0
    error('Unbound orbit! Check input r & v.');
end


end