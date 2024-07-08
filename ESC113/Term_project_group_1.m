% Constants and parameters
m = 2500;         % Mass (in kg)
to = 25;          % Initial temperature (in °C)
t_steam = 250;    % Steam temperature (in °C)
w1 = 200;         % Heat transfer rate (in W)
Cp = 2.1;         % Specific heat capacity (in J/kg·°C)
UA = 18.69;       % Overall heat transfer coefficient (in W/°C)
t1=25;            % Initial temperature of container 1
t2=25;            % Initial temperature of container 2
t3=25;            % Initial temperature of container 3
t4=25;            % Initial temperature of container 4

% Define the coefficient matrix A
A = [UA+w1*Cp, 0, 0, 0;
     w1*Cp, -(UA+w1*Cp), 0, 0;
     0, w1*Cp, -(UA+w1*Cp), 0;
     0, 0, w1*Cp, -(UA+w1*Cp)];

% Define the vector V
V = [w1*Cp*to+UA*t_steam;
     -UA*t_steam;
     -UA*t_steam;
     -UA*t_steam];

% Solve the system of equations using Gauss elimination
Tg = gaussEL(A, V); % Steady state temperature by Gauss elimination

% Solve the system of equations using Jacobi method
T1 = [t1; t2; t3; t4]; % Initial guess of temperature (initial temperature of container)

Tj = jacobi(A, V, T1); % Steady state temperature by Jacobi method


% Simulation parameters
t_final = 150;       % Final simulation time (in minutes)
delta_t = 0.1;      % Time step size (in minutes)

% Number of tanks
num_tanks = 4;

% Initialize variables
num_steps = t_final / delta_t + 1;
time = linspace(0, t_final, num_steps);
T = zeros(num_tanks, num_steps);
T(:, 1) = t1;

% Euler's explicit method
for i = 2:num_steps
    
    % Calculate the rate of temperature change for each tank
    dT1_dt = (w1 / m) * (to - T(1,i-1)) + (UA * (t_steam - T(1,i-1)) / (m * Cp));
    dT2_dt = (w1 / m) * (T(1,i-1) - T(2,i-1)) + (UA * (t_steam - T(2,i-1)) / (m * Cp));
    dT3_dt = (w1 / m) * (T(2,i-1) - T(3,i-1)) + (UA * (t_steam - T(3,i-1)) / (m * Cp));
    dT4_dt = (w1 / m) * (T(3,i-1) - T(4,i-1)) + (UA * (t_steam - T(4,i-1)) / (m * Cp));

    % Update the temperature using Euler's method
    T(1,i) = T(1,i-1) + dT1_dt * delta_t;
    T(2,i) = T(2,i-1) + dT2_dt * delta_t;
    T(3,i) = T(3,i-1) + dT3_dt * delta_t;
    T(4,i) = T(4,i-1) + dT4_dt * delta_t;

end

% Plotting
figure;
plot(time, T(1,:), time, T(2,:), time, T(3,:), time, T(4,:))

xlabel('Time (min)');
ylabel('Temperature (°C)');
title('Temperature vs. Time for each Tank');

% Function for solving the system of equations using Jacobi method
function X = jacobi(A, V, x)
    n = size(A, 1);
    x = [25; 25; 25; 25];
    epsl = 1;
    tol = 1e-8;
    i = 0;

    while (epsl > tol)
        i = i + 1;
        Xn = zeros(n, 1);
        for k = 1:n
            s = V(k);
            for j = 1:n
                if j ~= k
                    s = s - A(k, j) * x(j);
                end
            end
            Xn(k, 1) = s / A(k, k);
        end

        m = zeros(n, 1);
        s1 = 0;
        for z = 1:n
            m(z, 1) = Xn(z, 1) - x(z, 1);
            s1 = s1 + m(z, 1)^2;
        end
        epsl = s1^0.5;
        x = Xn;
    end
    
    X = x;
end

% Function for solving the system of equations using Gauss elimination
function x = gaussEL(A, b)
    n = size(A, 1);
    
    % Forward elimination
    for k = 1:n-1
        for i = k+1:n
            m = A(i, k) / A(k, k);
            A(i, k:n) = A(i, k:n) - m * A(k, k:n);
            b(i) = b(i) - m * b(k);
        end
    end

    % Backward substitution
    x = zeros(n, 1);
    x(n) = b(n) / A(n, n);
    for k = n-1:-1:1
        x(k) = (b(k) - A(k, k+1:n) * x(k+1:n)) / A(k, k);
    end
end
