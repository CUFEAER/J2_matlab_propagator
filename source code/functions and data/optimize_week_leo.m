%% 🚀 STRICT LEO OPTIMIZER (Apogee < 2000 km)
%  - This version finds the best possible 7-day orbit while enforcing
%    a hard constraint that the apogee must remain below 2000 km.

clc; clearvars; close all;

fprintf('=====================================================\n');
fprintf('   🛡️  STRICT LEO 7-DAY OPTIMIZER (Apogee < 2000km) 🛡️\n');
fprintf('=====================================================\n');

R_e = 6378;

% --- 1. SEED WITH THE BEST *POSSIBLE* LEO SHAPE ---
% We create an initial guess that is already at the LEO limit.
perigee_safe = 300;     % Safe for a week-long mission
apogee_limit = 1999.9;  % The hard constraint

rp = R_e + perigee_safe;
ra = R_e + apogee_limit;

start_a = (rp + ra) / 2;
start_e = (ra - rp) / (ra + rp);
start_i = 28.5;  % Ideal latitude match
start_w = 270.0; % Apogee at northernmost point
start_M = 0;

fprintf('Seeding with ideal LEO shape (Perigee: %.0fkm, Apogee: %.0fkm)\n', perigee_safe, apogee_limit);

% --- 2. PHASE 1: RAAN SCAN ---
% Find the best starting rotation for this ideal LEO shape.
fprintf('[Phase 1] Scanning RAAN to find best Earth synchronization...\n');

W_grid = 0:15:360; 
best_score = Inf;
best_W = 0;
warning('off', 'all');

for W = W_grid
    x_test = [start_a, start_e, start_i, W, start_w, start_M];
    % Use a medium time step for the scan
    score = cost_function_strict_leo(x_test, 120); 
    if score < best_score
        best_score = score;
        best_W = W;
    end
end
fprintf('  >> Best starting RAAN found: %.1f deg\n', best_W);

% --- 3. PHASE 2: HIGH-PRECISION REFINEMENT ---
fprintf('[Phase 2] Running high-precision refinement (dt=10s)...\n');
fprintf('  (This may take several minutes)...\n');

x0 = [start_a, start_e, start_i, best_W, start_w, start_M];
options = optimset('Display','iter', 'TolX', 1e-5, 'MaxIter', 250);
fun = @(x) cost_function_strict_leo(x, 10); 
final_x = fminsearch(fun, x0, options);

% --- 4. FINAL VALIDATION ---
fprintf('\n[Final Step] Validating final result...\n');
[final_score, total_minutes, perigee, apogee] = cost_function_strict_leo(final_x, 10);

% Unpack and Normalize
fa = final_x(1); fe = abs(final_x(2)); fi = mod(final_x(3), 360);
fW = mod(final_x(4), 360); fw = mod(final_x(5), 360); fM = mod(final_x(6), 360);

% Calc Nu
kepler_final = @(E) E - fe*sin(E) - deg2rad(fM);
try E_fin = fzero(kepler_final, deg2rad(fM)); catch, E_fin = deg2rad(fM); end
nu_rad = 2 * atan(sqrt((1+fe)/(1-fe)) * tan(E_fin/2));
fnu = mod(rad2deg(nu_rad), 360);
warning('on', 'all');

fprintf('\n\n╔════════════════════════════════════════════════════════╗\n');
fprintf('║          🏆 FINAL STRICT LEO RESULTS 🏆                ║\n');
fprintf('╠════════════════════════════════════════════════════════╣\n');
fprintf('║ Total Dwell Time: %7.2f Minutes (Over 7 Days)        ║\n', total_minutes);
fprintf('║ Average Daily:    %7.2f Minutes / Day                ║\n', total_minutes/7);
fprintf('║ Perigee: %5.0f km   |   Apogee: %5.0f km             ║\n', perigee, apogee);
fprintf('╠════════════════════════════════════════════════════════╣\n');
fprintf('║   ENTER THESE EXACT VALUES IN MODE 2:                  ║\n');
fprintf('║                                                        ║\n');
fprintf('║   a  : %12.4f                                   ║\n', fa);
fprintf('║   e  : %12.8f                                   ║\n', fe);
fprintf('║   i  : %12.4f                                   ║\n', fi);
fprintf('║   W  : %12.4f                                   ║\n', fW);
fprintf('║   w  : %12.4f                                   ║\n', fw);
fprintf('║   nu : %12.4f                                   ║\n', fnu);
fprintf('║                                                        ║\n');
fprintf('╚════════════════════════════════════════════════════════╝\n');


%% ==========================================================
%  🧠 COST FUNCTION (With Hard Apogee Limit)
%  =========================================================
function [cost, minutes, perigee_alt, apogee_alt] = cost_function_strict_leo(x, dt_sim)

    a = x(1); e = abs(x(2)); i = deg2rad(x(3));
    W0 = deg2rad(x(4)); w0 = deg2rad(x(5)); M0 = deg2rad(x(6));
    
    R_e = 6378; mu = 398600; J2 = 1.08263e-3;
    
    perigee_r = a * (1 - e);
    apogee_r  = a * (1 + e);
    perigee_alt = perigee_r - R_e;
    apogee_alt  = apogee_r - R_e;
    minutes = 0;
    
    % --- KEY MODIFICATION: HARD LEO CONSTRAINT ---
    % If the apogee ever exceeds 2000 km, apply a massive penalty.
    if apogee_alt > 2000
        cost = 1e9 + (apogee_alt - 2000) * 10000; % Huge penalty
        return;
    end
    
    % Safety for perigee
    if perigee_alt < 200
        cost = 1e9 + (200 - perigee_alt) * 10000;
        return;
    end
    
    % Physics
    Cd = 2.2; Area = 10e-6; Mass = 1000; Drag_Term = (Cd*Area)/Mass;
    We = 2 * pi * 366.25 / (24 * 3600 * 365.25);
    theta_GMST = getGMST(0); 
    
    target_vec = [R_e*cosd(27)*cosd(31); R_e*cosd(27)*sind(31); R_e*sind(27)];
    
    duration = 168 * 3600; 
    t = 0:dt_sim:duration;
    
    points_in_egypt = 0;
    dist_bonus = 0;
    min_dist_overall = Inf;
    
    a_curr = a; M_curr = M0; E_prev = M0; 
    sin_i = sin(i); cos_i = cos(i);
    
    for k = 1:length(t)
        time = t(k);
        
        n_curr = sqrt(mu/a_curr^3);
        p_curr = a_curr*(1-e^2);
        
        wd = 3/4 * J2 * sqrt(mu) * R_e^2 * (4 - 5 * sin_i^2) / p_curr^3.5;
        Wd = -3 * J2 * sqrt(mu) * R_e^2 * cos_i / (2 * p_curr^3.5);
        
        if k > 1, M_curr = M_curr + n_curr * dt_sim; end
        M_curr = mod(M_curr, 2*pi);
        
        % High precision fzero for final run, fast approx for scanning
        if dt_sim < 20
             eqn = @(E) E - e*sin(E) - M_curr;
             try E = fzero(eqn, E_prev); catch, E = M_curr; end
             E_prev = E;
        else
             E = M_curr;
             for iter=1:3, E = E - (E - e*sin(E) - M_curr)/(1 - e*cos(E)); end
        end

        nu = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
        r_val = p_curr / (1 + e * cos(nu));
        
        % Crash Check
        if r_val - R_e < 80, cost = 1e9; return; end
        
        % Drag
        v_mag = sqrt(mu * (2/r_val - 1/a_curr));
        rho = getDensityLocal(r_val - R_e);
        if rho > 0
             dadt = -2 * a_curr^2 * v_mag * (0.5 * rho * (v_mag*1000)^2 * Drag_Term/1000) / mu;
             a_curr = a_curr + dadt * dt_sim;
        end
        
        u = nu + w0 + wd * time;
        Om = W0 + Wd * time;
        
        X_ECI = r_val * [ cos(Om)*cos(u) - sin(Om)*sin(u)*cos_i;
                          sin(Om)*cos(u) + cos(Om)*sin(u)*cos_i;
                          sin(u)*sin_i ];
                          
        theta_rot = We * time + theta_GMST;
        
        X_ECEF = [ cos(theta_rot) sin(theta_rot) 0;
                  -sin(theta_rot) cos(theta_rot) 0;
                   0              0              1 ] * X_ECI;
                   
        d = norm(X_ECEF - target_vec);
        if d < min_dist_overall, min_dist_overall = d; end
        
        lat = asind(X_ECEF(3)/norm(X_ECEF));
        lon = atan2d(X_ECEF(2), X_ECEF(1));
        
        if (lat >= 22 && lat <= 32) && (lon >= 25 && lon <= 37)
            points_in_egypt = points_in_egypt + 1;
            dist_bonus = dist_bonus + (100 / (d + 1)); 
        end
    end
    
    minutes = (points_in_egypt * dt_sim) / 60;
    
    if points_in_egypt > 0
        cost = - (minutes * 1000) - dist_bonus;
    else
        cost = min_dist_overall;
    end
end

function rho = getDensityLocal(alt)
    if alt > 1000, rho = 0; return; end
    rho0 = 1.57e-13; h0 = 500; H = 63.4; 
    rho = rho0 * exp(-(alt - h0)/H);
end