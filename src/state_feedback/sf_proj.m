%% Modeling and Control of CPS - Part II Final Project
clc; clear all; close all;

% Saturation
USE_NU_SAT   = false;       
NU_MAX       = 5.0;       
USE_U_SAT    = false;       
U_MAX        = 8.0;      
 
% Noise
USE_NOISE    = false;
noise_std    = 0.02; 
 
% Output Disturbance 
USE_DIST     = false;
dist_amp     = 0.05;     
dist_freq    = 0.3;        
 
% Trajectory 
traj_type    = 'linear';   % 'linear' | 'circular' | 'sinusoidal'

%% SYSTEM DEFINITIONS & MATRICES
A = [0 1 0 0; 
     0 0 0 0; 
     0 0 0 1; 
     0 0 0 0];
B = [0 0; 
     1 0; 
     0 0; 
     0 1];
C = [1 0 0 0; 
     0 0 1 0];

N = 3; % Number of agents

% State structure: [px; vx; py; vy]
h1 = [0; 0; 0; 0];
h2 = [1; 0; 0; 0];
h3 = [0.5; 0; sqrt(3)/2; 0];
H = [h1, h2, h3];

%% GRAPH TOPOLOGY & PINNING 
% Undirect Graph 

Adj = [0 1 1;
       1 0 1;
       1 1 0];

D = diag(sum(Adj, 2));
L = D - Adj;
Pi_1 = diag([1 0 0]); % Only Agent 1 is pinned to reference

% Coupling Gains
c = 1;    % Inter-agent coupling gain
cl = 3;  % Tracking gain for leader 

% Verify stability of augmented Laplacian
Lp = c * L + cl * Pi_1;
alpha_p = min(real(eig(Lp)));
if alpha_p <= 0.5
    warning('Stability condition alpha_p > 0.5 not met. Increase c or cl.');
end

%% RICCATI CONTROLLER DESIGN 
Q = diag([10, 1, 10, 1]); 
 % Penalty on state error
R = eye(2);               % Penalty on input
P_ricc = are(A, B * inv(R) * B', Q);
K = inv(R) * B' * P_ricc; % Optimal gain matrix

%% OBSERVER DESIGN

Qo = diag([1e-4, 1, 1e-4, 1]); 
Ro = eye(2) * noise_std^2;
Po    = are(A', C' * inv(Ro) * C, Qo);
F = Po * C' * inv(Ro);
co = 1;

%% SIMULATION PARAMETERS
dt = 0.01;      
Tf = 20;        
time = 0:dt:Tf;
Nt = length(time);

% Parametri Moto Lineare (Velocità Costante)
v_lin = [0.4; 0.2];   % Vettore velocità costante [vx; vy]
r0_lin = [0; 0];      % Posizione iniziale del riferimento [x0; y0]
% Moto circolare
omega_c = 0.5; R_c = 4;
% Moto sinusoidale
v_sin = 0.5; A_sin = 3; omega_sin = 0.8;

% State initialization
% X_real: [state_dim, agents, time]
X_real = zeros(4, N, Nt);
X_hat  = zeros(4, N, Nt);
U_phys = zeros(2, N, Nt);
nu_virt = zeros(2, N, Nt);

% Initial Conditions 
X_real(:,:,1) = [ -2,  0, 2;
                   0,  0, 0;
                  -2, -4, 2;
                   0,  0, 0 ];

% Initialize observer states with small error (casual)
X_hat(:,:,1) = X_real(:,:,1) + 0.5 * ones(4, 3); 

%% MAIN SIMULATION LOOP (Forward Euler)
for k = 1:Nt-1
    t = time(k);

 % Reference trajectory
 % computing reference using get_reference
    [xr, u_ff] = get_reference(t, traj_type, v_lin, r0_lin, omega_c, R_c, v_sin, A_sin, omega_sin);
    
    % Measurement with sensor noise
    Y_meas = zeros(2, N);
    for i = 1:N
        noise    = zeros(2,1);
        dist     = zeros(2,1);
        if USE_NOISE
            noise = noise_std * randn(2, 1);
        end
        if USE_DIST 
            % Sinusoidal disturbance: different phase per agent to make it realistic
            dist = dist_amp * sin(dist_freq * t + (i-1)*pi/3) * [1; 1];
        end
            Y_meas(:, i) = C * X_real(:, i, k) + noise + dist;
    end
    
    % distributed protocol & virutal input (nu) 
    % Evaluated using the estimated states X_hat
    nu = zeros(2, N);
    
    for i = 1:N
        % Local error formation for all nodes
        eps_i = zeros(4, 1);
        for j = 1:N
            if Adj(i, j) > 0
                eps_i = eps_i + Adj(i, j) * ((X_hat(:, j, k) - X_hat(:, i, k)) - (H(:, j) - H(:, i)));
            end
        end

        % Cooperative distributed term 

        u_form = c * K * eps_i;
        
        % Tracking term (exclusively on node 1)
        u_track = [0; 0];
        if i == 1
            u_track = cl * K * (xr - X_hat(:, 1, k));
        end
        nu_raw = u_form + u_track + u_ff;
 
        % Saturation on virtual input nu 
        if USE_NU_SAT
            nu_raw = saturate(nu_raw, NU_MAX);
        end
        nu(:, i) = nu_raw;
    end
    nu_virt(:,:,k) = nu;

    % Inverse feedback linearization 
    u1 = fb_lin_agent1(X_hat(:,1,k), nu(:,1));
    u2 = fb_lin_agent2(X_hat(:,2,k), nu(:,2));
    u3 = fb_lin_agent3(X_hat(:,3,k), nu(:,3));
    
    if USE_U_SAT
        u1 = saturate(u1, U_MAX);
        u2 = saturate(u2, U_MAX);
        u3 = saturate(u3, U_MAX);
    end
 
    U_phys(:,1,k) = u1;
    U_phys(:,2,k) = u2;
    U_phys(:,3,k) = u3;
    % Dynamics & Observer update
    for i = 1:N
        % physical dynamics plant
        if i == 1
            x_dot = plant_agent1(X_real(:, i, k), U_phys(:, i, k));
        elseif i == 2
            x_dot = plant_agent2(X_real(:, i, k), U_phys(:, i, k));
        else
            x_dot = plant_agent3(X_real(:, i, k), U_phys(:, i, k));
        end
        
        % Observer dynamics
        x_hat_dot = A * X_hat(:, i, k) + B * nu(:, i) + co*F * (Y_meas(:, i) - C * X_hat(:, i, k));
        
        % integration
        X_real(:, i, k+1) = X_real(:, i, k) + dt * x_dot;
        X_hat(:, i, k+1)  = X_hat(:, i, k) + dt * x_hat_dot;
    end
end

%% COMPUTE ERRORS & METRICS FOR RESULTS TABLE
err_form  = zeros(Nt, 1);
err_track = zeros(Nt, 1);
 
for k = 1:Nt
    xr_curr = get_ref_pos(time(k), traj_type, v_lin, r0_lin, omega_c, R_c, v_sin, A_sin);
    err_track(k) = norm(X_real([1,3], 1, k) - xr_curr);
    e12 = norm((X_real([1,3],2,k) - X_real([1,3],1,k)) - (H([1,3],2) - H([1,3],1)));
    e13 = norm((X_real([1,3],3,k) - X_real([1,3],1,k)) - (H([1,3],3) - H([1,3],1)));
    err_form(k) = e12 + e13;
end
 
% Convergence time (threshold: error < 0.1 m)
THRESH = 0.1;
conv_idx_form  = find(err_form  < THRESH, 1, 'first');
conv_idx_track = find(err_track < THRESH, 1, 'first');
t_conv_form    = NaN; if ~isempty(conv_idx_form),  t_conv_form  = time(conv_idx_form);  end
t_conv_track   = NaN; if ~isempty(conv_idx_track), t_conv_track = time(conv_idx_track); end
 
% Steady-state window: last 10 s
ss_idx = time >= (Tf - 10);
ss_err_form  = mean(err_form(ss_idx));
ss_err_track = mean(err_track(ss_idx));
peak_err_form  = max(err_form);
peak_err_track = max(err_track);
 
% Max physical inputs per agent
max_u = zeros(2, N);
for i = 1:N
    max_u(:,i) = max(abs(squeeze(U_phys(:,i,1:end-1))), [], 2);
end
 
fprintf('\n========== RESULTS TABLE ==========\n');
fprintf('Trajectory: %s\n', traj_type);
fprintf('Saturation nu: %d (max=%.1f)  |  u: %d (max=%.1f)\n', USE_NU_SAT, NU_MAX, USE_U_SAT, U_MAX);
fprintf('Noise: %d (std=%.3f)  |  Disturbance: %d (amp=%.3f)\n', USE_NOISE, noise_std, USE_DIST, dist_amp);
fprintf('---\n');
fprintf('Peak tracking error:    %.4f m\n', peak_err_track);
fprintf('Peak formation error:   %.4f m\n', peak_err_form);
fprintf('Conv. time (tracking):  %.2f s\n', t_conv_track);
fprintf('Conv. time (formation): %.2f s\n', t_conv_form);
fprintf('SS tracking error:      %.4f m\n', ss_err_track);
fprintf('SS formation error:     %.4f m\n', ss_err_form);
fprintf('Max |u| Agent1: [%.2f, %.2f]\n', max_u(1,1), max_u(2,1));
fprintf('Max |u| Agent2: [%.2f, %.2f]\n', max_u(1,2), max_u(2,2));
fprintf('Max |u| Agent3: [%.2f, %.2f]\n', max_u(1,3), max_u(2,3));
fprintf('====================================\n\n');

%% 8. PLOTS
 
% Agent Trajectories
figure('Name', 'Agent Trajectories (XY Plane)');
hold on; grid on;
colors = {'b', 'r', 'g'};
for i = 1:N
    plot(squeeze(X_real(1,i,:)), squeeze(X_real(3,i,:)), colors{i}, 'LineWidth', 1.5, 'DisplayName', sprintf('Agent %d', i));
end
plot_ref_traj(time, traj_type, v_lin, r0_lin, omega_c, R_c, v_sin, A_sin, omega_sin);
xlabel('X [m]'); ylabel('Y [m]');
title(['Nonlinear Multi-Agent Formation — ', upper(traj_type)]);
legend('Location', 'best'); axis equal;
 
% Formation & Tracking Errors
figure('Name', 'Formation Errors Over Time');
subplot(2,1,1);
plot(time, err_track, 'b', 'LineWidth', 1.5); hold on; grid on;
yline(THRESH, 'k--', 'Threshold');
xlabel('Time [s]'); ylabel('Error [m]');
title('Leader Tracking Error'); legend('err\_track', 'Threshold');
 
subplot(2,1,2);
plot(time, err_form, 'r', 'LineWidth', 1.5); hold on; grid on;
yline(THRESH, 'k--', 'Threshold');
xlabel('Time [s]'); ylabel('Error [m]');
title('Formation Structure Error'); legend('err\_form', 'Threshold');
 
% Physical Inputs u per agent
figure('Name', 'Physical Inputs U_phys');
agent_labels = {'Agent 1', 'Agent 2', 'Agent 3'};
ax_labels    = {'u_x', 'u_y'};
for i = 1:N
    subplot(N, 2, 2*(i-1)+1);
    plot(time(1:end-1), squeeze(U_phys(1,i,1:end-1)), colors{i}, 'LineWidth', 1.2);
    if USE_U_SAT
        yline( U_MAX, 'k--'); yline(-U_MAX, 'k--');
    end
    grid on;
    xlabel('Time [s]'); ylabel('u_x [N]');
    title([agent_labels{i}, ' — u_x']);
 
    subplot(N, 2, 2*(i-1)+2);
    plot(time(1:end-1), squeeze(U_phys(2,i,1:end-1)), colors{i}, 'LineWidth', 1.2);
    if USE_U_SAT
        yline( U_MAX, 'k--'); yline(-U_MAX, 'k--');
    end
    grid on;
    xlabel('Time [s]'); ylabel('u_y [N]');
    title([agent_labels{i}, ' — u_y']);
end
sgtitle('Physical Inputs (dashed = saturation bound)');
 
% Virtual Inputs nu
figure('Name', 'Virtual Inputs nu');
for i = 1:N
    subplot(N, 2, 2*(i-1)+1);
    plot(time(1:end-1), squeeze(nu_virt(1,i,1:end-1)), colors{i}, 'LineWidth', 1.2);
    if USE_NU_SAT
        yline( NU_MAX, 'k--'); yline(-NU_MAX, 'k--');
    end
    grid on;
    xlabel('Time [s]'); ylabel('\nu_x [m/s^2]');
    title([agent_labels{i}, ' — \nu_x']);
 
    subplot(N, 2, 2*(i-1)+2);
    plot(time(1:end-1), squeeze(nu_virt(2,i,1:end-1)), colors{i}, 'LineWidth', 1.2);
    if USE_NU_SAT
        yline( NU_MAX, 'k--'); yline(-NU_MAX, 'k--');
    end
    grid on;
    xlabel('Time [s]'); ylabel('\nu_y [m/s^2]');
    title([agent_labels{i}, ' — \nu_y']);
end
sgtitle('Virtual Inputs (dashed = saturation bound)');
 
% Observer estimation error
figure('Name', 'Observer Estimation Error');
obs_err = zeros(N, Nt);
for i = 1:N
    for k = 1:Nt
        obs_err(i,k) = norm(X_real(:,i,k) - X_hat(:,i,k));
    end
end
hold on; grid on;
for i = 1:N
    plot(time, obs_err(i,:), colors{i}, 'LineWidth', 1.5, 'DisplayName', agent_labels{i});
end
xlabel('Time [s]'); ylabel('||x - x_{hat}|| [m]', 'Interpreter', 'tex');
title('Luenberger Observer Estimation Error');
legend('Location','best');
 
% Animation
figure('Name', 'Multi-Agent Formation Animation');
hold on; grid on;
plot_ref_traj(time, traj_type, v_lin, r0_lin, omega_c, R_c, v_sin, A_sin, omega_sin);
for i = 1:N
    plot(squeeze(X_real(1,i,:)), squeeze(X_real(3,i,:)), [colors{i}, ':'], 'LineWidth', 0.5);
end
h_robots = gobjects(N,1);
for i = 1:N
    h_robots(i) = plot(X_real(1,i,1), X_real(3,i,1), 'o', 'Color', colors{i}, ...
                       'MarkerFaceColor', colors{i}, 'MarkerSize', 8);
end
link_x = [X_real(1,1,1), X_real(1,2,1), X_real(1,3,1), X_real(1,1,1)];
link_y = [X_real(3,1,1), X_real(3,2,1), X_real(3,3,1), X_real(3,1,1)];
h_shape = plot(link_x, link_y, 'k-', 'LineWidth', 1.5);
legend([h_robots; h_shape], {'Agent 1 (Leader)', 'Agent 2', 'Agent 3', 'Formation'}, 'Location','best');
margin = 3;
xlim([min(X_real(1,:,:),[],'all')-margin, max(X_real(1,:,:),[],'all')+margin]);
ylim([min(X_real(3,:,:),[],'all')-margin, max(X_real(3,:,:),[],'all')+margin]);
title(['Animation — ', upper(traj_type)]); xlabel('X [m]'); ylabel('Y [m]'); axis equal;
 
for k = 1:10:Nt
    for i = 1:N
        set(h_robots(i), 'XData', X_real(1,i,k), 'YData', X_real(3,i,k));
    end
    set(h_shape, 'XData', [X_real(1,1,k), X_real(1,2,k), X_real(1,3,k), X_real(1,1,k)], ...
                 'YData', [X_real(3,1,k), X_real(3,2,k), X_real(3,3,k), X_real(3,1,k)]);
    drawnow; pause(0.005);
end

%% LOCAL FUNCTIONS: FEEDBACK LINEARIZATION 

function v = saturate(v, vmax)
    v = max(min(v, vmax), -vmax);
end
 
function [xr, u_ff] = get_reference(t, traj_type, v_lin, r0_lin, omega_c, R_c, v_sin, A_sin, omega_sin)
    if strcmp(traj_type, 'linear')
        xr    = [r0_lin(1)+v_lin(1)*t; v_lin(1); r0_lin(2)+v_lin(2)*t; v_lin(2)];
        u_ff  = [0; 0];
    elseif strcmp(traj_type, 'circular')
        xr    = [R_c*cos(omega_c*t); -R_c*omega_c*sin(omega_c*t); ...
                 R_c*sin(omega_c*t);  R_c*omega_c*cos(omega_c*t)];
        u_ff  = [-R_c*omega_c^2*cos(omega_c*t); -R_c*omega_c^2*sin(omega_c*t)];
    else
        xr    = [v_sin*t; v_sin; A_sin*sin(omega_sin*t); A_sin*omega_sin*cos(omega_sin*t)];
        u_ff  = [0; -A_sin*omega_sin^2*sin(omega_sin*t)];
    end
end
 
function pos = get_ref_pos(t, traj_type, v_lin, r0_lin, omega_c, R_c, v_sin, A_sin)
    if strcmp(traj_type, 'linear')
        pos = [r0_lin(1)+v_lin(1)*t; r0_lin(2)+v_lin(2)*t];
    elseif strcmp(traj_type, 'circular')
        pos = [R_c*cos(omega_c*t); R_c*sin(omega_c*t)];
    else
        pos = [v_sin*t; A_sin*sin(omega_sin*t)];
    end
end
 
function plot_ref_traj(time, traj_type, v_lin, r0_lin, omega_c, R_c, v_sin, A_sin, omega_sin)
    if strcmp(traj_type, 'linear')
        plot(r0_lin(1)+v_lin(1)*time, r0_lin(2)+v_lin(2)*time, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
    elseif strcmp(traj_type, 'circular')
        plot(R_c*cos(omega_c*time), R_c*sin(omega_c*time), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
    else
        plot(v_sin*time, A_sin*sin(omega_sin*time), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
    end
end

function u_phys = fb_lin_agent1(x_hat, nu)
    % Agent 1
    v_x = x_hat(2);
    v_y = x_hat(4);
    phi = [ -0.35*v_x - 0.08*v_x*abs(v_x);
            -0.45*v_y - 0.10*v_y*abs(v_y) ];
    B_mat = diag([1/1.2, 1/1.2]);
    %input-output feedback linearizing 
    u_phys = B_mat \ (nu - phi);
end

function u_phys = fb_lin_agent2(x_hat, nu)
    % Agent 2
    p_x = x_hat(1);
    v_x = x_hat(2);
    p_y = x_hat(3);
    v_y = x_hat(4);
    phi = [ -0.40*v_x;
            -0.55*v_y ];
    B_mat = diag([ 1/(1 + 0.25*sin(p_x)^2), 1/(1 + 0.30*cos(p_y)^2) ]);
    u_phys = B_mat \ (nu - phi);
end

function u_phys = fb_lin_agent3(x_hat, nu)
    % Agent 3
    p_x = x_hat(1);
    v_x = x_hat(2);
    p_y = x_hat(3);
    v_y = x_hat(4);
    phi = [ -0.30*v_x + 0.15*sin(p_y)*v_y;
            -0.35*v_y + 0.15*cos(p_x)*v_x ];
    B_mat = [ 1 + 0.20*cos(p_y)^2, 0.08*sin(p_x);
              0.08*sin(p_y),       1 + 0.20*sin(p_x)^2 ];
    u_phys = B_mat \ (nu - phi);
end

%% LOCAL FUNCTIONS: NONLINEAR PLANT DYNAMICS

function x_dot = plant_agent1(x, u)
    v_x = x(2); v_y = x(4);
    phi = [ -0.35*v_x - 0.08*v_x*abs(v_x);
            -0.45*v_y - 0.10*v_y*abs(v_y) ];
    B_mat = diag([1/1.2, 1/1.2]);
    accel = phi + B_mat * u;
    x_dot = [ v_x; accel(1); v_y; accel(2) ];
end

function x_dot = plant_agent2(x, u)
    p_x = x(1); v_x = x(2); p_y = x(3); v_y = x(4);
    phi = [ -0.40*v_x;
            -0.55*v_y ];
    B_mat = diag([ 1/(1 + 0.25*sin(p_x)^2), 1/(1 + 0.30*cos(p_y)^2) ]);
    accel = phi + B_mat * u;
    x_dot = [ v_x; accel(1); v_y; accel(2) ];
end

function x_dot = plant_agent3(x, u)
    p_x = x(1); v_x = x(2); p_y = x(3); v_y = x(4);
    phi = [ -0.30*v_x + 0.15*sin(p_y)*v_y;
            -0.35*v_y + 0.15*cos(p_x)*v_x ];
    B_mat = [ 1 + 0.20*cos(p_y)^2, 0.08*sin(p_x);
              0.08*sin(p_y),       1 + 0.20*sin(p_x)^2 ];
    accel = phi + B_mat * u;
    x_dot = [ v_x; accel(1); v_y; accel(2) ];
end