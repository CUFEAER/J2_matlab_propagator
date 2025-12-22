function runAttitudeSimulation(duration)
    % --- Initial Conditions  ---
    wx0 = 0.2;
    wy0 = 0.4;
    wz0 = 0.3;
    
    epsi0  = deg2rad(45);
    theta0 = deg2rad(45);
    phi0   = deg2rad(45);
    
    Ixx = 2;
    Iyy = 1;
    Izz = 0.5;
    
    Mx = 0;
    My = 0;
    Mz = 0;
    
    initial_position = [wx0; wy0; wz0; epsi0; theta0; phi0];
    t_span = [0, duration];
    
    % --- Solve ODE ---
    [t1, y1] = ode45(@(t,y) omega_inertial_angles(t, y, Ixx, Iyy, Izz, Mx, My, Mz), t_span, initial_position);

    % --- Plotting  ---
    figure();
    plot(t1, y1(:, 1:3));
    grid on;
    legend('omega x', 'omega y', 'omega z', 'Location', 'best');
    title('Angular Velocities');

    figure();
    plot(t1, rad2deg(y1(:, 4:6))); % Converted to degrees for readability in plot
    grid on;
    legend('epsi', 'theta', 'phi', 'Location', 'best');
    title('Euler Angles');

    figure();
    plot(t1, y1);
    grid on;
    legend('omega x', 'omega y', 'omega z', 'epsi', 'theta', 'phi', 'Location', 'best');
    title('Combined State');
end

function dydt = omega_inertial_angles(t, y, Ixx, Iyy, Izz, Mx, My, Mz)

    % Unpack state for clarity
    % y(1)=wx, y(2)=wy, y(3)=wz
    % y(4)=epsi, y(5)=theta, y(6)=phi

    % --- DYNAMICS (Fixed missing * symbols) ---
    omega_dot = [ (Mx - y(2)*y(3)*(Izz-Iyy)) / Ixx ; ...
                  (My - y(1)*y(3)*(Ixx-Izz)) / Iyy ; ...
                  (Mz - y(1)*y(2)*(Iyy-Ixx)) / Izz ];

    omega = [y(1); y(2); y(3)];

    % --- KINEMATICS ( ---
    ia = (1/sin(y(5))) * ...
        [ sin(y(6)),            cos(y(6)),             0 ; ...
          cos(y(6))*sin(y(5)), -sin(y(6))*sin(y(5)),   0 ; ...
         -sin(y(6))*cos(y(5)),  -cos(y(6))*cos(y(5)),   sin(y(5)) ];

    iad = ia * omega;

    dydt = [omega_dot; iad];
end