%% =======================
%  🚀 ORBIT SIM MENU SYSTEM
%  =======================

clc; clearvars; close all;
addpath("functions and data")
%---global constants---
R_e = 6378;    % Earth's equatorial radius [km]
mu = 398600;   % Earth's gravitational parameter [km^3/s^2]
j = 1.08263e-3; % J2 term (related to Earth's oblateness)

% --- Drag Physics Constants (Added) ---
Cd_std = 2.2;       % Standard Drag Coefficient
Area_m2 = 10;       % Cross-section [m^2]
Mass_kg = 1000;     % Mass [kg]


fprintf("╔═══════════════════════════════════════════════════════════════════╗\n");
fprintf("║                SATELLITE ORBIT SIMULATOR (J2 MODEL)               ║\n");
fprintf("║          Developed by CUFE Aerospace Department Students          ║\n");
fprintf("╠═══════════════════════════════════════════════════════════════════╣\n");
fprintf("║ Team Members:                                                     ║\n");
fprintf("║   Ahmed Magdy      Ahmed Sayed                                    ║\n");
fprintf("║   Sohail Sayed     Zeyad Alaa                                     ║\n");
fprintf("╠═══════════════════════════════════════════════════════════════════╣\n");
fprintf("║ Redistribution is free, provided the credited names remain intact ║\n");
fprintf("╚═══════════════════════════════════════════════════════════════════╝\n\n");


%% MODIFICATION START: Added option 6 to the menu
fprintf("Choose input mode:\n");
fprintf("  1) Load elements from TLE file\n");
fprintf("  2) Enter classical orbital elements manually\n");
fprintf("  3) Enter state vectors r & v (km, km/s)\n");
fprintf("  4) Convert state vectors r & v  to COE\n");
fprintf("  5) Run Attitude Dynamics Simulation\n"); % New Option
fprintf("  6) Quit\n");

mode = input("Enter choice (1/2/3/4/5/6): ");

if mode == 6
    disp("Exiting.");
    return;
end

while mode > 6 || mode < 1
    mode = input('please enter one of the choices you have: ');
end
%% MODIFICATION END

%% ================================
% 🚀 MODE 1: Extract from TLE
% ================================
switch mode
    case 1
        fprintf("\n--- TLE MODE SELECTED ---\n");
        fname = input("Enter TLE filename (e.g., 'tleiss.txt'): ", 's');

        % === TLE extraction ===
        [epoch, ndot, i_deg, W0_deg, e, w0_deg, M_deg, n, bstar, epoch_date_formatted] = readTLE(fname);

        fprintf("TLE Loaded Successfully ✔️\n");

        % Convert and compute needed stuff
        M = deg2rad(M_deg);
        E0 = solveKeplerInitial(M, e);
        n_rad = n * 2*pi / 86400;
        ndot_rad_s2 = (ndot * 2) * (2*pi) / (86400^2);
        a = (398600 / n_rad^2)^(1/3);
        
        Drag_Term = (Cd_std * Area_m2) / Mass_kg;


        i = deg2rad(i_deg);
        W0 = deg2rad(W0_deg);
        w0 = deg2rad(w0_deg);



%% ================================
% 🚀 MODE 2: Enter COEs manually
% ================================
    case 2
        fprintf("\n--- MANUAL MODE SELECTED ---\n");

        a = input("Semi-major axis a [km]: ");
        e = input("Eccentricity e: ");
        i_deg = input("Inclination i [deg]: ");
        W0_deg = input("RAAN Ω [deg]: ");
        w0_deg = input("Argument of perigee ω [deg]: ");
        nu0_deg = input("True anomaly ν0 [deg]: ");

        % Conversions
        i = deg2rad(i_deg);
        W0 = deg2rad(W0_deg);
        w0 = deg2rad(w0_deg);
        nu0 = deg2rad(nu0_deg);

        % Compute mean anomaly
        E0 = 2 * atan( tan(nu0/2) * sqrt((1-e)/(1+e)) );
        M = E0 - e*sin(E0);

        % Mean motion
        n_rad = sqrt(mu / a^3);

        ndot_rad_s2 = 0;     % no drag
        
        % --- Drag Setup ---
        % Drag_Term = (Cd * Area_km2) / Mass_kg
        Drag_Term = (Cd_std * (Area_m2)) / Mass_kg;
        
        epoch = 0;


%% ================================
% 🚀 MODE 3: Input r & v manually
% ================================
    case 3
        fprintf("\n--- STATE VECTOR MODE SELECTED ---\n");
        fprintf("Input position vector r (km):\n");
        rx = input("  r_x = ");
        ry = input("  r_y = ");
        rz = input("  r_z = ");

        fprintf("Input velocity vector v (km/s):\n");
        vx = input("  v_x = ");
        vy = input("  v_y = ");
        vz = input("  v_z = ");

        r = [rx ry rz];
        v = [vx vy vz];

        % === Convert to orbital elements ===
        [a, e, i, W0, w0, nu, M] = rv2coe(r, v);

        % Compute eccentric anomaly
        E0 = 2 * atan( tan(nu/2) * sqrt((1-e)/(1+e)) );

        % Mean motion
        mu = 398600;
        n_rad = sqrt(mu / a^3);

        ndot_rad_s2 = 0;
        
        % --- Drag Setup ---
        Drag_Term = (Cd_std * (Area_m2)) / Mass_kg;
        
        epoch = 0;

    case 4
        fprintf("\n--- STATE VECTOR TO COE MODE SELECTED ---\n");
        fprintf("Input position vector r (km):\n");
        rx = input("  r_x = ");
        ry = input("  r_y = ");
        rz = input("  r_z = ");

        fprintf("Input velocity vector v (km/s):\n");
        vx = input("  v_x = ");
        vy = input("  v_y = ");
        vz = input("  v_z = ");

        r = [rx ry rz];
        v = [vx vy vz];

        % === Convert to orbital elements ===
[a, e, i, W0, w0, nu, M] = rv2coe(r, v);

fprintf('\n================ ORBITAL ELEMENTS ================\n');
fprintf(' Semi-major axis (a)        : %12.6f km\n', a);
fprintf(' Eccentricity (e)           : %12.8f\n', e);
fprintf(' Inclination (i)            : %12.6f deg\n', rad2deg(i));
fprintf(' RAAN (Ω)                   : %12.6f deg\n', rad2deg(W0));
fprintf(' Argument of Perigee (ω)    : %12.6f deg\n', rad2deg(w0));
fprintf(' True Anomaly (ν)           : %12.6f deg\n', rad2deg(nu));
fprintf(' Mean Anomaly (M)           : %12.6f deg\n', rad2deg(M));
fprintf('==================================================\n\n');
        return;

%% MODIFICATION START: Add case for Attitude Simulation
    case 5
        fprintf("\n--- ATTITUDE DYNAMICS SIMULATION ---\n");
        duration_att = input("Enter simulation duration (seconds): ");
        AttitudeSimulation(duration_att);
        fprintf("\nAttitude simulation complete. Plots are displayed.\n");
        return; % Exit the main script after showing plots
%% MODIFICATION END

end



%% ================================
% 🚀 Continue your simulation here
% ================================

    fprintf('\n================ ORBITAL ELEMENTS ================\n');
    fprintf(' Semi-major axis (a)        : %12.6f km\n', a);
    fprintf(' Eccentricity (e)           : %12.8f\n', e);
    fprintf(' Inclination (i)            : %12.6f deg\n', rad2deg(i));
    fprintf(' RAAN (Ω)                   : %12.6f deg\n', rad2deg(W0));
    fprintf(' Argument of Perigee (ω)    : %12.6f deg\n', rad2deg(w0));
    fprintf(' Mean Anomaly (M)           : %12.6f deg\n', rad2deg(M));
    fprintf('==================================================\n\n');

fprintf("\nSimulation duration (hours): ");
sim_duration_hours = input('');

fprintf("Time step (seconds): ");
time_step_seconds = input('');

fprintf("\nRunning orbit simulation...\n");
pause(0.5);

% Time vector creation
sim_duration_seconds = sim_duration_hours * 3600;
t = 0:time_step_seconds:sim_duration_seconds;

% Derived Parameters
p = a * (1 - e^2);             % Semi-latus rectum [km]
T = 2 * pi * sqrt(a^3 / mu);   % Orbital period [seconds]

% Earth's sidereal rotation rate [rad/s]
We = 2 * pi * 366.25 / (24 * 3600 * 365.25);

% Calculate GMST Angle for map alignment
theta_GMST = getGMST(epoch);
%% 3.  perturbing forces Calculation (J2 Rates)

% Note: With Drag, 'a' and 'p' change over time. 
% We will recalculate wd and Wd inside the loop now.
% Initial values for reference:
wd = 3/4 * j * sqrt(mu) * R_e^2 * (4 - 5 * sin(i)^2) / p^3.5;
Wd = -3 * j * sqrt(mu) * R_e^2 * cos(i) / (2 * p^3.5);

%% 4. 📐 Solving Kepler's Equation

% 4.2. Initial time offset (ti)
% Corresponds to the initial Mean Anomaly ($M_0 = E_0 - e \sin(E_0)$)
ti = T / (2 * pi) * (E0 - e * sin(E0)); 
%will set to 0 for there is supposed to be none change depending on the situation
%ti=0;

% 4.3. Solve Kepler's Equation for E over time using fzero
E = zeros(size(t));
M_curr = M;
% Arrays for loop storage
nu = zeros(size(t));
r  = zeros(size(t));
xi = zeros(size(t));
yi = zeros(size(t));
zi = zeros(size(t));
E_prev = E0; 

for k = 1:length(t)
    dt = t(k);


    if a <= (R_e + 80) % Stop if altitude < 80km 
        fprintf('\n🚨 SATELLITE RE-ENTERED ATMOSPHERE AT t = %.2f HOURS 🚨\n', t(k)/3600);
        fprintf('Simulation stopped to prevent math errors.\n');
        k_stop = k - 1; % Mark where data ends
        break;
    end

    % 1. Update Mean Motion: Dynamic Drag Update
    % Calculate Mean Motion from CURRENT Semi-Major Axis
    n_curr = sqrt(mu / a^3);
    
    % 2. Update Semi-Major Axis (Kepler 3rd Law: a^3 ~ 1/n^2)
    a_curr = a; % Use the current decaying 'a'
    
    % 3. Update Period and p
    T_curr = 2 * pi / n_curr;
    p_curr = a_curr * (1 - e^2);

    % 4. Update J2 Rates (because p changed)
    wd = 3/4 * j * sqrt(mu) * R_e^2 * (4 - 5 * sin(i)^2) / p_curr^3.5;
    Wd = -3 * j * sqrt(mu) * R_e^2 * cos(i) / (2 * p_curr^3.5);
    % ---------------------------------

     if k == 1
        dt_step = t(1); 
    else
        dt_step = t(k) - t(k-1);
    end
    
    M_curr = M_curr + n_curr * dt_step;
    
    % B. Normalize M to 0..2pi (This FIXES the fzero error)
    M_curr = mod(M_curr, 2*pi);
    
    % C. Solve Kepler's Equation: M = E - e*sin(E)
    % We solve for E such that: E - e*sin(E) - M = 0
    % We normalize the guess (E_prev) as well to keep fzero happy
    guess = mod(E_prev, 2*pi);
    
    % Simple anonymous function 
    kepler_eqn = @(E_sol) E_sol - e * sin(E_sol) - M_curr;
    
    E(k) = fzero(kepler_eqn, guess);
    E_prev = E(k);


    % 4.4. Calculate True Anomaly (nu)
    % $\nu = 2 \arctan\left(\sqrt{\frac{1+e}{1-e}} \tan(E/2)\right)$
    nu(k) = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E(k) / 2));

    % 4.5. Calculate Radius (r)
    % Use p_curr (updated with drag)
    r(k) = p_curr ./ (1 + e .* cos(nu(k)));
    
    % --- ATMOSPHERIC DRAG UPDATE  ---
    % 1. Calculate Velocity Magnitude (Vis-Viva)
    v_mag = sqrt(mu * (2/r(k) - 1/a_curr));
    % 2. Calculate Decay Rate (da/dt)
    dadt = calc_dadt(r(k), v_mag, a_curr, Drag_Term);
    % 3. Update 'a' for the NEXT step
    a = a + dadt * time_step_seconds; 
    % ---------------------------------------
    
    %% 5. 🌐 Position in ECI Frame (3D Orbit) (Moving inside loop for optimization)
    
    % Position in the **Perifocal Frame** (x, y) [km]
    x_k = r(k) .* cos(nu(k));
    y_k = r(k) .* sin(nu(k));
    
    % Time-varying orbital elements due to J2 effect
    W = W0 + Wd * dt; % Updated RAAN
    w = w0 + wd * dt; % Updated Argument of Perigee

    % Transformation matrix (Perifocal to ECI frame)
    % Accounts for W, i, and w.
    % Using new local function To2i
    trn = To2i(W,i,w);

    % Convert Perifocal position (x, y, 0) to ECI (xi, yi, zi)
    X_ECI = trn * [x_k; y_k; 0];
    xi(k) = X_ECI(1);
    yi(k) = X_ECI(2);
    zi(k) = X_ECI(3);

    % --- Compute velocity vector in ECI ---

% Perifocal velocities:
vx_p = -sqrt(mu/p_curr) * sin(E(k));
vy_p =  sqrt(mu/p_curr) * sqrt(1-e^2) * cos(E(k));

% Rotate to ECI
V_ECI = trn * [vx_p; vy_p; 0];

vx(k) = V_ECI(1);
vy(k) = V_ECI(2);
vz(k) = V_ECI(3);

end

% Plotting E and nu for reference 
figure(1);
subplot(2, 1, 1);
plot(t / 3600, rad2deg(E));
title('Eccentric Anomaly (E)'); xlabel('Time [hours]'); ylabel('E [degrees]'); grid on;
subplot(2, 1, 2);
plot(t / 3600, rad2deg(nu));
title('True Anomaly (\nu)'); xlabel('Time [hours]'); ylabel('\nu [degrees]'); grid on;


% Plotting 3D Orbit 
figure(2);
plot3(xi, yi, zi, 'b-');
hold on;
% Simple sphere representation of Earth
[sphere_x, sphere_y, sphere_z] = sphere(20);
surf(R_e * sphere_x, R_e * sphere_y, R_e * sphere_z, 'FaceColor', [0.1 0.1 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
xlabel('ECI X [km]'); ylabel('ECI Y [km]'); zlabel('ECI Z [km]');
title('3D Satellite Orbit (ECI Frame)');
axis equal; grid on;

%% 6. 🗺️ Ground Track Projection (ECEF Frame)

% Convert ECI coordinates to ECEF coordinates to account for Earth's rotation (We)
X_ECEF      = zeros(size(t));
Y_ECEF      = zeros(size(t));
Z_ECEF      = zeros(size(t));
Alt_anim    = zeros(size(t));


for k = 1:length(t)
    % Added theta_GMST to align map correctly
    W_earth_rot = We * t(k) + theta_GMST; % Angle of Earth's rotation
    
    %% MODIFICATION START: Replace manual matrix with call to Ti2e function
    % Rotation matrix from ECI to ECEF: $R_Z(\omega_e t)$
    % R_ECI_to_ECEF = [cos(W_earth_rot), sin(W_earth_rot), 0;
    %                  -sin(W_earth_rot), cos(W_earth_rot), 0;
    %                  0, 0, 1];
    R_ECI_to_ECEF = Ti2e(W_earth_rot);
    %% MODIFICATION END

    X_ECEF_k = R_ECI_to_ECEF * [xi(k); yi(k); zi(k)];

    Alt_anim(k)=norm([xi(k); yi(k); zi(k)])-R_e;

    X_ECEF(k) = X_ECEF_k(1);
    Y_ECEF(k) = X_ECEF_k(2);
    Z_ECEF(k) = X_ECEF_k(3);
end

Ri = sqrt(X_ECEF.^2 + Y_ECEF.^2 + Z_ECEF.^2); % Radius to satellite

% Calculate Latitude (phi) and Longitude (gamma)
phi = asin(Z_ECEF ./ Ri);       % Latitude: $\phi = \arcsin(Z / R)$
gamma = atan2(Y_ECEF, X_ECEF);  % Longitude: $\lambda = \arctan2(Y, X)$

%% ===== LOCAL SOLAR TIME CALCULATION =====

% Earth rotation rate (you already defined We earlier):
% We = 2 * pi * 366.25 / (24 * 3600 * 365.25);

% Local Hour Angle: GMST(t) + longitude
theta_t = theta_GMST + We .* t;      % Total Earth rotation [rad]

% LST in hours:
% Convert theta_t + longitude (gamma) → degrees → divide by 15
LST_hours = mod( rad2deg(theta_t + gamma) / 15 , 24 );

% Optional: display sample values
fprintf("\n===== Local Solar Time Samples =====\n");
fprintf("At t = 0 hours:  LST = %.2f\n", LST_hours(1));
fprintf("At t = %.1f h: LST = %.2f\n", t(end)/3600, LST_hours(end));
fprintf("====================================\n");


lat_deg_all = rad2deg(phi);
lon_deg_all = rad2deg(gamma);

is_in_egypt = (lat_deg_all >= 22 & lat_deg_all <= 32) & ...
              (lon_deg_all >= 25 & lon_deg_all <= 37);

% 3. Calculate Total Time
% Count how many 'dots' are in Egypt and multiply by the time step
points_in_egypt = sum(is_in_egypt);
total_minutes_egypt = (points_in_egypt * time_step_seconds) / 60;

minutes_over_egypt = total_minutes_egypt ;
disp('------------------------------------------------');
disp(['Total Dwell Time over Egypt: ', num2str(minutes_over_egypt), ' minutes']);
disp('------------------------------------------------');

% Plotting Ground Track 
figure(3);
% Load map data (requires 'earth.mat' file)
load('earth.mat');
plot(long, lat, '.' );
hold on;

plot(rad2deg(gamma), rad2deg(phi), '.');
plot(rad2deg(gamma(1)),rad2deg(phi(1)),'go',"MarkerSize",8,"MarkerFaceColor","g")
plot(rad2deg(gamma(end)),rad2deg(phi(end)),'ro',"MarkerSize",8,"MarkerFaceColor","r")
xlabel('Longitude [degrees]');
ylabel('Latitude [degrees]');
title('Satellite Ground Track Projection on Earth');
xlim([-180, 180]); ylim([-90, 90]);
grid on;

%% ===== Animation or not decider =====


fprintf("\nDisplay options:\n");
fprintf("  1) Animate the ground track\n");
fprintf("  2) Show satellite position at a specific date & time\n");
choice_view = input("Select (1 or 2): ");

if choice_view == 1
%% ===== Animation with single-loop tracker (dots only) =====
figure(4); clf; hold on;    % <-- force animation in figure(4)
load('earth.mat');

% Draw Earth's coastlines
plot(long, lat, '.', 'Color', [0.4 0.7 1], 'MarkerSize', 0.5);
hold on;

% Convert full orbit to degrees
lon_full = rad2deg(gamma);
lat_full = rad2deg(phi);

% Unwrap true anomaly to detect loops
nu_unwrapped = unwrap(nu);
loop_idx = floor((nu_unwrapped - nu_unwrapped(1)) / (2*pi)) + 1;
num_loops = max(loop_idx);

% Precompute loop start/end indices
loop_starts = zeros(1, num_loops);
loop_ends   = zeros(1, num_loops);
for kk = 1:num_loops
    idxs = find(loop_idx == kk);
    loop_starts(kk) = idxs(1);
    loop_ends(kk)   = idxs(end);
end

% Satellite markers
plot(lon_full(1),  lat_full(1),  'go','MarkerFaceColor','g','MarkerSize',8);
plot(lon_full(end),lat_full(end),'ro','MarkerFaceColor','r','MarkerSize',8);

% Initialize handles
h_loop = [];
h_night_patch = [];   % night shading patch handle
h_term_line = [];     % terminator line handle
h_sat = plot(lon_full(1), lat_full(1), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

% Animation loop
for j = 1:100:length(t)
    figure(4); % Ensure we're updating the right figure
    
    % ===== DELETE previous NIGHT SHADE + TERMINATOR =====
    % Clear previous shading objects ONLY (keep Earth map)
    if ~isempty(h_night_patch) && isgraphics(h_night_patch)
        delete(h_night_patch);
    end
    if ~isempty(h_term_line) && isgraphics(h_term_line)
        delete(h_term_line);
    end
    
    if mode == 1
        % === Compute Sun position ===
        currentUTC = datetime(epoch_date_formatted) + seconds(t(j));
        [subsolar_lon, subsolar_lat] = subsolarPoint(currentUTC);
        [lonTerm, latTerm] = sunTerminator(subsolar_lon, subsolar_lat);
        
        % === Create CLEAN night shading (like in static mode) ===
        % The key: create ONE polygon for the night side
        h_night_patch = patchCleanNightShading(lonTerm, latTerm);
        
        % === Draw day/night boundary ===
        h_term_line = plot(lonTerm, latTerm, 'w-', 'LineWidth', 1.4);
    end
    
    % Remove previous loop and satellite
    if ~isempty(h_loop) && isgraphics(h_loop)
        delete(h_loop);
    end
    if ~isempty(h_sat) && isgraphics(h_sat)
        delete(h_sat);
    end
    
    curr_loop = loop_idx(j);
    
    % Current loop range
    idx_loop_start = loop_starts(curr_loop);
    idx_loop_end   = loop_ends(curr_loop);
    
    % Plot current loop
    h_loop = plot(lon_full(idx_loop_start:idx_loop_end), ...
                  lat_full(idx_loop_start:idx_loop_end), '.', ...
                  'MarkerSize', 3, 'Color', [0.05 0.6 0.95]);
    
    % Plot satellite
    h_sat = plot(lon_full(j), lat_full(j), 'or', ...
                 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    
    %% ===== ⭐ LOCAL Sidereal TIME CALCULATION for THIS FRAME =====
    theta_t = theta_GMST + We * t(j);     % Earth rotation angle at this time
    LST_hours = mod( rad2deg(theta_t + gamma(j)) / 15 , 24 );
    
    %% ===== Update title =====
    title(sprintf(['Satellite Ground Track - Live\n',...
                   'Lat: %.2f°   Lon: %.2f°   Altitude: %.2fkm   Loop: %d/%d   LST: %.2f hours'], ...
          lat_full(j), lon_full(j), Alt_anim(j), curr_loop, num_loops, LST_hours), ...
          'FontSize', 10);
    
    xlim([-180 180]);
    ylim([-90 90]);
    grid on;
    
    pause(0.05);  % Smooth animation speed
    drawnow;
end
else
    %% ===== STATIC MODE (Show one point in time) =====
if mode == 1
disp("Enter date/time to inspect:");
year  = str2double(input('  Year (YYYY): ', 's'));
month = str2double(input('  Month (1–12): ', 's'));
day   = str2double(input('  Day (1–31): ', 's'));
hour  = str2double(input('  Hour (0–23): ', 's'));
minute= str2double(input('  Minute (0–59): ', 's'));
second= str2double(input('  Second (0–59): ', 's'));


% Ask for time zone offset
fprintf('\nEnter time zone offset from UTC (e.g., +2 for Egypt, -5 for EST):\n');
tz_offset = input('  UTC±HH (e.g., 2, -5, 0): ');

% Create datetime in local time (no timezone)
user_time_local = datetime(year, month, day, hour, minute, second);

% Convert local time to UTC
user_time_utc = user_time_local - hours(tz_offset);
user_time_utc.TimeZone = 'UTC';

% Use UTC time for calculations
user_time = user_time_utc;

fprintf('\nLocal time: %s\n', datestr(user_time_local));
fprintf('UTC time:   %s\n', datestr(user_time_utc));


% Convert to seconds from epoch

t_user = seconds(user_time - epoch_date_formatted);   % epoch from TLE mode

% If manual mode → epoch = 0, user inputs absolute time instead
else
    t_user = input("Enter time since start of simulation [seconds]: ");
end

% Find closest simulation point
[~, idx] = min(abs(t - t_user));
% ===== Print position & velocity at user-chosen time =====
fprintf("\n===== STATE AT USER-SELECTED DATE/TIME =====\n");
fprintf("Time input: %f seconds since epoch\n", t_user);

fprintf("\nECI Position [km]:\n");
fprintf("  x = %.3f\n  y = %.3f\n  z = %.3f\n", xi(idx), yi(idx), zi(idx));

fprintf("\nECI Velocity [km/s]:\n");
fprintf("  vx = %.6f\n  vy = %.6f\n  vz = %.6f\n", vx(idx), vy(idx), vz(idx));
fprintf("=============================================");

Alt=norm([xi(idx), yi(idx), zi(idx)])-R_e;
fprintf("\nAltitude (h) [km] = %.4f\n", Alt);
fprintf("=============================================\n\n");

% Extract that point
lon_pt = rad2deg(gamma(idx));
lat_pt = rad2deg(phi(idx));
if mode == 1
% getting LST
[LAST,LMST]=getLST(user_time,lon_pt);

% ===== LAST & LMST =====
if LMST > 1 && LAST < 1 
    fprintf('======== Printing Local Mean/Apperant Solar Time ========\n')
    fprintf('Local Mean solar time = %3.4f hr\n',LMST)
    fprintf('Local Apperant solar time = %3.4f min\n',LAST*60)
    fprintf("=============================================\n\n");

elseif LMST < 1 && LAST > 1
    fprintf('======== Printing Local Mean/Apperant Solar Time ========\n')
    fprintf('Local Mean solar time = %3.4f min\n',LMST*60)
    fprintf('Local Apperant solar time = %3.4f hr\n',LAST)
    fprintf("=============================================\n\n");

elseif LMST < 1 && LAST < 1
    fprintf('======== Printing Local Mean/Apperant Solar Time ========\n')
    fprintf('Local Mean solar time = %3.4f min\n',LMST*60)
    fprintf('Local Apperant solar time = %3.4f min\n',LAST*60)
    fprintf("=============================================\n\n");

else
    fprintf('======== Printing Local Mean/Apperant Solar Time ========\n')
    fprintf('Local Mean solar time = %3.4f hr\n',LMST)
    fprintf('Local Apperant solar time = %3.4f hr\n',LAST)
    fprintf("=============================================\n\n");
end
end
% ==== Compute Local Sidreal Time ====
theta_t = theta_GMST + We * t(idx);
LST_hours = mod( rad2deg(theta_t + gamma(idx))/15 , 24 );

%% ===== PLOT on FIGURE 4 =====
figure(4); clf; hold on;
load('earth.mat'); 
plot(long, lat, '.', 'Color', [0.4 0.7 1], 'MarkerSize', 0.5);
hold on;

% Plot satellite point
plot(lon_pt, lat_pt, 'ro', 'MarkerFaceColor','r','MarkerSize',9);
if mode == 1
% ========== SUNLIGHT OVERLAY ==========
currentUTC = user_time;
% === Compute Sun position ===
[subsolar_lon, subsolar_lat] = subsolarPoint2(user_time);

% === Terminator curve ===
[lonTerm, latTerm] = sunTerminator(subsolar_lon, subsolar_lat);

% === Shade night ===
shadeNight(lonTerm, latTerm, subsolar_lon,subsolar_lat);

% === Draw day/night boundary ===
plot(lonTerm, latTerm, 'w.', 'LineWidth', 1.4);



% Title including real time + LST
title(sprintf(['Satellite Position @ %s UTC\n',...
               'Lat: %.2f°   Lon: %.2f°\n' ...
               'LST(Local Siderial Time): %.2f hours'], ...
               datestr(user_time), lat_pt, lon_pt, LST_hours));
else
title(sprintf(['Satellite Position after t = %.0f sec\n',...
               'Lat: %.2f°   Lon: %.2f°\n'], ...
               t_user, lat_pt, lon_pt));
end

xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
xlim([-180 180]); ylim([-90 90]);
grid on; 
return;   % Skip animation section

end
