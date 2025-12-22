function dadt = calc_dadt(r, v, a, Drag_Term)
    % CALC_DADT Calculates rate of change of Semi-Major Axis (da/dt)
    % Inputs:
    %   r          : Current Position Magnitude [km] OR Vector [km]
    %   v          : Current Velocity Magnitude [km/s] OR Vector [km/s]
    %   a          : Current Semi-Major Axis [km]
    %   Drag_Term  : Ballistic Factor (Cd * A / m). Units: [m^2 / kg]
    %
    % Output:
    %   dadt       : Decay rate [km/s] (Negative value)

    % --- Constants ---
    R_e = 6378;          % Earth Radius [km]
    mu  = 398600;        % GM [km^3/s^2]
    
    % Exponential Atmosphere Model (Approx for LEO ~200-800km)
    rho0 = 1.57e-13;     % Reference Density at 500km [kg/m^3]
    h0   = 500;          % Reference Altitude [km]
    H    = 63.4;         % Scale Height [km]

    % --- Handle Inputs ---
    if length(r) == 3, r_mag = norm(r); else, r_mag = r; end
    if length(v) == 3, v_mag = norm(v); else, v_mag = v; end

    % --- 1. Calculate Density (rho) ---
    alt = r_mag - R_e; % Altitude [km]
    
    if alt > 1000
        rho = 0; % Vacuum approximation
    else
        rho = rho0 * exp(-(alt - h0) / H); % [kg/m^3]
    end

    % --- 2. Calculate Drag Acceleration ---
    % a_drag = 0.5 * rho * v^2 * (Cd * A / m)
    
    % Convert v to m/s for physics calc
    v_metric = v_mag * 1000; 
    
    % Acceleration in [m/s^2]
    accel_metric = 0.5 * rho * (v_metric^2) * Drag_Term;
    
    % Convert back to [km/s^2]
    accel_km = accel_metric / 1000;

    % --- 3. Calculate da/dt ---
    % Energy decay: da/dt = -2 * a^2 * v * a_drag / mu
    dadt = -2 * (a^2) * v_mag * accel_km / mu;
end