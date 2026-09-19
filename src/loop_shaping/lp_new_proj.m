%% Part II - Formation Control of Nonlinear Planar Mobile Robots (Loop Shaping)
clear; clc; close all;

%% Setting parameters
% Loop-Shaping parameters
p.Kp = 1;         
p.Kd = 1;          
p.gamma = 2.0;       

% Controller
p.USE_LEAD = false;      
p.tau_lead = 0.02;  

p.SHOW_FREQUENCY_ANALYSIS = true;

% Saturation
p.USE_NU_SAT = true;     
p.NU_MAX = 5.0;       
p.USE_U_SAT = true;    
p.U_MAX = 8.0;       
 
%noise and disturbances
p.USE_NOISE = true;       
p.noise_std = 0.02;       

p.USE_DIST = true;
p.dist_amp = 0.05;      
p.dist_freq = 5;       
 
% trajectory selector
p.traj_type = 'linear';   % linear | circular | sinusoidal

% Trajectory parameters
p.v_lin = [0.4; 0.2];   
p.r0_lin = [0; 0];
p.omega_c = 0.5;        
p.R_c = 4;
p.v_sin = 0.5;          
p.A_sin = 3; 
p.omega_sin = 0.8;

%% formation setup
N = 3;
% Grafo indiretto connesso: 1<->2, 2<->3, 3<->1
A_adj = [0 1 1; 1 0 1; 1 1 0]; 
D_deg = diag(sum(A_adj, 2));
L = D_deg - A_adj;

% laplaciana ancorata (Pinning sull'agente 1)
Pi1 = diag([1 zeros(1,N-1)]); 
p.Lp = L + p.gamma * Pi1;

% analisi modale
[Vp, Lambdap] = eig(p.Lp); 
lambda_vals = diag(Lambdap);
% vettori di formazione (triangolo in movimento)
p.H = [0, 1, 0.5; 
       0, 0, sqrt(3)/2]; 

%% Bode & Nyquist
if p.SHOW_FREQUENCY_ANALYSIS
    s_tf = tf('s');
    G = 1 / (s_tf^2);

    if p.USE_LEAD
        K_ctrl = p.Kp + (p.Kd * s_tf) / (1 + p.tau_lead * s_tf);
        ctrl_title = sprintf('Lead Compensator  (Kp=%.1f, Kd=%.1f, \\tau=%.2f)', p.Kp, p.Kd, p.tau_lead);
    else
        K_ctrl = p.Kp + p.Kd * s_tf;
        ctrl_title = sprintf('PD Puro  (Kp=%.1f, Kd=%.1f)', p.Kp, p.Kd);
    end

    colors = {[0.839 0.153 0.157], [0.122 0.467 0.706], [0.173 0.627 0.173]};
    lw         = 1.8;
    fs_title   = 13;
    fs_label   = 11;

    % Pre-calcolo
    L_loops     = cell(1, length(lambda_vals));
    S_sens_all  = cell(1, length(lambda_vals));
    W_sens_all  = cell(1, length(lambda_vals));
    leg_entries = cell(1, length(lambda_vals));
    for i = 1:length(lambda_vals)
        mu_i           = lambda_vals(i);
        L_loops{i}     = G * K_ctrl * mu_i;
        S_sens_all{i}  = 1 / (1 + L_loops{i});
        W_sens_all{i}  = L_loops{i} / (1 + L_loops{i});
        leg_entries{i} = sprintf('\\mu = %.2f', mu_i);
    end

    hFig = figure('Name', 'Frequency Domain Analysis', ...
                  'Position', [100 80 1000 620], 'Color', 'w');
    tg = uitabgroup(hFig);

    % Bode
    tab1 = uitab(tg, 'Title', 'Bode  L(s)');
    ax1a = axes('Parent', tab1, 'Position', [0.10 0.54 0.82 0.36]);
    hold(ax1a, 'on'); grid(ax1a, 'on'); box(ax1a, 'on');
    ax1b = axes('Parent', tab1, 'Position', [0.10 0.10 0.82 0.36]);
    hold(ax1b, 'on'); grid(ax1b, 'on'); box(ax1b, 'on');

    for i = 1:length(lambda_vals)
        [mag, ph, wout] = bode(L_loops{i});
        mag_dB = 20*log10(squeeze(mag));
        ph_deg = squeeze(ph);
        semilogx(ax1a, wout, mag_dB, 'Color', colors{i}, 'LineWidth', lw, ...
                 'DisplayName', leg_entries{i});
        semilogx(ax1b, wout, ph_deg, 'Color', colors{i}, 'LineWidth', lw, ...
                 'DisplayName', leg_entries{i});
    end
    yline(ax1a,   0, 'k--', '0 dB',   'LineWidth', 0.8, 'LabelHorizontalAlignment', 'left');
    yline(ax1b, -180, 'k--', '-180°', 'LineWidth', 0.8, 'LabelHorizontalAlignment', 'left');
    ylabel(ax1a, 'Magnitudine (dB)', 'FontSize', fs_label);
    ylabel(ax1b, 'Fase (deg)',       'FontSize', fs_label);
    xlabel(ax1b, 'Frequenza (rad/s)','FontSize', fs_label);
    title(ax1a, sprintf('Bode  L(s) = G(s)\\cdotK(s)\\cdot\\mu_i  —  %s', ctrl_title), ...
          'FontSize', fs_title, 'FontWeight', 'bold');
    legend(ax1a, 'Location', 'southwest', 'FontSize', 10);

    % Nyquist
    tab2 = uitab(tg, 'Title', 'Nyquist  L(s)');
    ax2  = axes('Parent', tab2, 'Position', [0.10 0.10 0.82 0.80]);
    hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');

    for i = 1:length(lambda_vals)
        [re, im, ~] = nyquist(L_loops{i});
        re = squeeze(re); im = squeeze(im);
        plot(ax2,  re,  im, 'Color', colors{i}, 'LineWidth', lw,      'DisplayName', leg_entries{i});
        plot(ax2,  re, -im, '--', 'Color', colors{i}, 'LineWidth', lw*0.7, 'HandleVisibility', 'off');
    end
    plot(ax2, -1, 0, 'kx', 'MarkerSize', 14, 'LineWidth', 2.5, 'DisplayName', 'Punto critico (-1,0)');
    xline(ax2, 0, 'Color', [0.75 0.75 0.75]);
    yline(ax2, 0, 'Color', [0.75 0.75 0.75]);
    xlim(ax2, [-3 1.5]); ylim(ax2, [-2.5 2.5]);
    xlabel(ax2, 'Re\{L(j\omega)\}', 'FontSize', fs_label);
    ylabel(ax2, 'Im\{L(j\omega)\}', 'FontSize', fs_label);
    title(ax2,  sprintf('Nyquist  L(s)  —  %s', ctrl_title), ...
          'FontSize', fs_title, 'FontWeight', 'bold');
    legend(ax2, 'Location', 'best', 'FontSize', 10);

    % Sensitivity 
    tab3 = uitab(tg, 'Title', 'Sensitività  S(s)');
    ax3  = axes('Parent', tab3, 'Position', [0.10 0.12 0.82 0.76]);
    hold(ax3, 'on'); grid(ax3, 'on'); box(ax3, 'on');

    for i = 1:length(lambda_vals)
        [mag_s, ~, w_s] = bode(S_sens_all{i});
        mag_s_dB = 20*log10(squeeze(mag_s));
        semilogx(ax3, w_s, mag_s_dB, 'Color', colors{i}, 'LineWidth', lw, ...
                 'DisplayName', leg_entries{i});
    end
    xlabel(ax3, 'Frequenza (rad/s)',    'FontSize', fs_label);
    ylabel(ax3, '|S(j\omega)|  (dB)',   'FontSize', fs_label);
    title(ax3,  sprintf('Sensitività  S(s)  —  %s', ctrl_title), ...
          'FontSize', fs_title, 'FontWeight', 'bold');
    legend(ax3, 'Location', 'southeast', 'FontSize', 10);

    % Sensitività Complementare
    tab4 = uitab(tg, 'Title', 'Sens. Complementare  W(s)');
    ax4  = axes('Parent', tab4, 'Position', [0.10 0.12 0.82 0.76]);
    hold(ax4, 'on'); grid(ax4, 'on'); box(ax4, 'on');

    for i = 1:length(lambda_vals)
        [mag_w, ~, w_w] = bode(W_sens_all{i});
        mag_w_dB = 20*log10(squeeze(mag_w));
        semilogx(ax4, w_w, mag_w_dB, 'Color', colors{i}, 'LineWidth', lw, ...
                 'DisplayName', leg_entries{i});
    end
    yline(ax4, -3, 'k--', '-3 dB  (f_{bw})', 'LineWidth', 1.2, ...
          'LabelHorizontalAlignment', 'left', 'FontSize', 10);
    xlabel(ax4, 'Frequenza (rad/s)',         'FontSize', fs_label);
    ylabel(ax4, '|W(j\omega)|  (dB)',        'FontSize', fs_label);
    title(ax4,  sprintf('Sensitività Complementare  W(s)  —  %s', ctrl_title), ...
          'FontSize', fs_title, 'FontWeight', 'bold');
    legend(ax4, 'Location', 'southwest', 'FontSize', 10);
end
%% SIMULATION
% Condizioni iniziali [px, py, vx, vy] per i 3 agenti (12 stati fisici)
X0_agents = [ -2; -1;  0;  0; 
               2; -2;  0;  0; 
               0;  3;  0;  0];   

% Condizioni iniziali per il filtro Lead Compensator (6 stati aggiuntivi z)
X0_filter = zeros(6, 1);

% Vettore di stato iniziale completo (18 elementi)
X0 = [X0_agents; X0_filter];    

tspan = [0, 100];

options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);
[t, X_out] = ode45(@(t,X) formation_dynamics(t, X, p), tspan, X0, options);

%% POST-PROCESSING ED ESTRAZIONE DATI
p1_x = X_out(:,1); p1_y = X_out(:,2);
p2_x = X_out(:,5); p2_y = X_out(:,6);
p3_x = X_out(:,9); p3_y = X_out(:,10);

% Ricalcolo riferimento per plot
ref_x = zeros(length(t),1); ref_y = zeros(length(t),1);
for k = 1:length(t)
    [r, ~] = get_reference(t(k), p);
    ref_x(k) = r(1); ref_y(k) = r(2);
end

% Errori di formazione
err_21_x = (p2_x - p1_x) - (p.H(1,2) - p.H(1,1));
err_21_y = (p2_y - p1_y) - (p.H(2,2) - p.H(2,1));
err_31_x = (p3_x - p1_x) - (p.H(1,3) - p.H(1,1));
err_31_y = (p3_y - p1_y) - (p.H(2,3) - p.H(2,1));

%% PLOT STATICI E ANIMAZIONE
% ricalcolo u per verifica della saturazione
Nt = length(t);
U1 = zeros(2, Nt); U2 = zeros(2, Nt); U3 = zeros(2, Nt);

for k = 1:Nt
    tk = t(k);
    Xk = X_out(k, :)';
    
    % Chiamiamo la dinamica per ottenere le accelerazioni effettive (dv)
    dXk = formation_dynamics(tk, Xk, p);
    dv1 = dXk(3:4);
    dv2 = dXk(7:8);
    dv3 = dXk(11:12);
    
    % Estraiamo stati fisici necessari per la IOFL inversa
    v1 = Xk(3:4);
    p2 = Xk(5:6); v2 = Xk(7:8);
    p3 = Xk(9:10); v3 = Xk(11:12);
    
    % Inversione della dinamica per ricalcolare esattamente l'input u applicato
    phi1 = [-0.35*v1(1) - 0.08*v1(1)*abs(v1(1)); -0.45*v1(2) - 0.10*v1(2)*abs(v1(2))];
    B1 = diag([1/1.2, 1/1.2]);
    U1(:, k) = B1 \ (dv1 - phi1);
    
    phi2 = [-0.40*v2(1); -0.55*v2(2)];
    B2 = diag([1 / (1 + 0.25*sin(p2(1))^2), 1 / (1 + 0.30*cos(p2(2))^2)]);
    U2(:, k) = B2 \ (dv2 - phi2);
    
    phi3 = [-0.30*v3(1) + 0.15*sin(p3(2))*v3(2); -0.35*v3(2) + 0.15*cos(p3(1))*v3(1)];
    B3 = [1 + 0.20*cos(p3(2))^2, 0.08*sin(p3(1)); 0.08*sin(p3(2)), 1 + 0.20*sin(p3(1))^2];
    U3(:, k) = B3 \ (dv3 - phi3);
end

hFigSim = figure('Name', 'Risultati Simulazione Temporale', 'Position', [150 100 1100 700], 'Color', 'w');
tgSim = uitabgroup(hFigSim);

% Sforzo di controllo
tabU = uitab(tgSim, 'Title', 'Sforzo di Controllo (U)');

axU1x = subplot(3,2,1, 'Parent', tabU);
plot(axU1x, t, U1(1,:), 'r', 'LineWidth', 1.2); grid(axU1x, 'on'); hold(axU1x, 'on');
if p.USE_U_SAT, yline(axU1x, p.U_MAX, 'k--', 'Sat Max'); yline(axU1x, -p.U_MAX, 'k--'); end
title(axU1x, 'Agente 1 - Input u_x'); ylabel(axU1x, 'Forza [N]');

axU1y = subplot(3,2,2, 'Parent', tabU);
plot(axU1y, t, U1(2,:), 'r', 'LineWidth', 1.2); grid(axU1y, 'on'); hold(axU1y, 'on');
if p.USE_U_SAT, yline(axU1y, p.U_MAX, 'k--', 'Sat Max'); yline(axU1y, -p.U_MAX, 'k--'); end
title(axU1y, 'Agente 1 - Input u_y');

axU2x = subplot(3,2,3, 'Parent', tabU);
plot(axU2x, t, U2(1,:), 'b', 'LineWidth', 1.2); grid(axU2x, 'on'); hold(axU2x, 'on');
if p.USE_U_SAT, yline(axU2x, p.U_MAX, 'k--'); yline(axU2x, -p.U_MAX, 'k--'); end
title(axU2x, 'Agente 2 - Input u_x'); ylabel(axU2x, 'Forza [N]');

axU2y = subplot(3,2,4, 'Parent', tabU);
plot(axU2y, t, U2(2,:), 'b', 'LineWidth', 1.2); grid(axU2y, 'on'); hold(axU2y, 'on');
if p.USE_U_SAT, yline(axU2y, p.U_MAX, 'k--'); yline(axU2y, -p.U_MAX, 'k--'); end
title(axU2y, 'Agente 2 - Input u_y');

axU3x = subplot(3,2,5, 'Parent', tabU);
plot(axU3x, t, U3(1,:), 'g', 'LineWidth', 1.2); grid(axU3x, 'on'); hold(axU3x, 'on');
if p.USE_U_SAT, yline(axU3x, p.U_MAX, 'k--'); yline(axU3x, -p.U_MAX, 'k--'); end
title(axU3x, 'Agente 3 - Input u_x'); xlabel(axU3x, 'Tempo [s]'); ylabel(axU3x, 'Forza [N]');

axU3y = subplot(3,2,6, 'Parent', tabU);
plot(axU3y, t, U3(2,:), 'g', 'LineWidth', 1.2); grid(axU3y, 'on'); hold(axU3y, 'on');
if p.USE_U_SAT, yline(axU3y, p.U_MAX, 'k--'); yline(axU3y, -p.U_MAX, 'k--'); end
title(axU3y, 'Agente 3 - Input u_y'); xlabel(axU3y, 'Tempo [s]');

%Traiettorie X-Y
tabXY = uitab(tgSim, 'Title', 'Traiettorie X-Y');
axXY = axes('Parent', tabXY);
plot(axXY, p1_x, p1_y, 'r', 'LineWidth', 1.5); hold(axXY, 'on'); grid(axXY, 'on');
plot(axXY, p2_x, p2_y, 'b', 'LineWidth', 1.5);
plot(axXY, p3_x, p3_y, 'g', 'LineWidth', 1.5);
plot(axXY, ref_x, ref_y, 'k--', 'LineWidth', 1);

step = floor(length(t)/5);
for k = 1:step:length(t)
    plot(axXY, [p1_x(k) p2_x(k) p3_x(k) p1_x(k)], [p1_y(k) p2_y(k) p3_y(k) p1_y(k)], 'k-', 'LineWidth', 0.5);
end
title(axXY, ['Traiettorie Non Lineari — ', upper(p.traj_type)]);
legend(axXY, 'Agente 1', 'Agente 2', 'Agente 3', 'Riferimento', 'Location', 'best');
xlabel(axXY, 'Posizione X [m]'); ylabel(axXY, 'Posizione Y [m]'); axis(axXY, 'equal');

% Errori di formazione
tabErr = uitab(tgSim, 'Title', 'Errori di Formazione');
axEx = subplot(2,1,1, 'Parent', tabErr);
plot(axEx, t, err_21_x, 'b', t, err_31_x, 'g', 'LineWidth', 1.5); grid(axEx, 'on');
title(axEx, 'Errore di Formazione Lungo X'); legend(axEx, 'Agent 2 vs 1', 'Agent 3 vs 1');

axEy = subplot(2,1,2, 'Parent', tabErr);
plot(axEy, t, err_21_y, 'b', t, err_31_y, 'g', 'LineWidth', 1.5); grid(axEy, 'on');
title(axEy, 'Errore di Formazione Lungo Y'); xlabel(axEy, 'Tempo [s]');

% Animazione
tabAnim = uitab(tgSim, 'Title', 'Animazione Live');
axAnim = axes('Parent', tabAnim);
hold(axAnim, 'on'); grid(axAnim, 'on'); axis(axAnim, 'equal');
plot(axAnim, ref_x, ref_y, 'k--', 'LineWidth', 1);

margine = 2;
xlim(axAnim, [min([p1_x; p2_x; p3_x]) - margine, max([p1_x; p2_x; p3_x]) + margine]);
ylim(axAnim, [min([p1_y; p2_y; p3_y]) - margine, max([p1_y; p2_y; p3_y]) + margine]);

h_scia1 = plot(axAnim, p1_x(1), p1_y(1), 'r-', 'LineWidth', 0.8);
h_scia2 = plot(axAnim, p2_x(1), p2_y(1), 'b-', 'LineWidth', 0.8);
h_scia3 = plot(axAnim, p3_x(1), p3_y(1), 'g-', 'LineWidth', 0.8);

h_triangolo = plot(axAnim, [p1_x(1) p2_x(1) p3_x(1) p1_x(1)], [p1_y(1) p2_y(1) p3_y(1) p1_y(1)], 'k-', 'LineWidth', 1.5);

h_agente1 = plot(axAnim, p1_x(1), p1_y(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
h_agente2 = plot(axAnim, p2_x(1), p2_y(1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
h_agente3 = plot(axAnim, p3_x(1), p3_y(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

tgSim.SelectedTab = tabAnim;

step_animazione = max(1, floor(length(t) / 300));
for k = 1:step_animazione:length(t)
    if ~isvalid(hFigSim), break; end 
    
    set(h_agente1, 'XData', p1_x(k), 'YData', p1_y(k));
    set(h_agente2, 'XData', p2_x(k), 'YData', p2_y(k));
    set(h_agente3, 'XData', p3_x(k), 'YData', p3_y(k));
    set(h_triangolo, 'XData', [p1_x(k) p2_x(k) p3_x(k) p1_x(k)], 'YData', [p1_y(k) p2_y(k) p3_y(k) p1_y(k)]);
    set(h_scia1, 'XData', p1_x(1:k), 'YData', p1_y(1:k));
    set(h_scia2, 'XData', p2_x(1:k), 'YData', p2_y(1:k));
    set(h_scia3, 'XData', p3_x(1:k), 'YData', p3_y(1:k));
    drawnow;
    pause(0.01); % Decommenta se vuoi rallentare l'animazione
end

%% Function

function [r, dr] = get_reference(t, p)
    if strcmp(p.traj_type, 'linear')
        r  = p.r0_lin + p.v_lin * t;
        dr = p.v_lin;
    elseif strcmp(p.traj_type, 'circular')
        r  = [p.R_c*cos(p.omega_c*t); p.R_c*sin(p.omega_c*t)];
        dr = [-p.R_c*p.omega_c*sin(p.omega_c*t); p.R_c*p.omega_c*cos(p.omega_c*t)];
    else % sinusoidal
        r  = [p.v_sin*t; p.A_sin*sin(p.omega_sin*t)];
        dr = [p.v_sin; p.A_sin*p.omega_sin*cos(p.omega_sin*t)];
    end
end

function v = saturate(v, vmax)
    for idx = 1:length(v)
        v(idx) = max(min(v(idx), vmax), -vmax);
    end
end

function dX = formation_dynamics(t, X, p)
    
    %Variabili di stato fisiche 
    p1 = X(1:2); v1 = X(3:4);
    p2 = X(5:6); v2 = X(7:8);
    p3 = X(9:10); v3 = X(11:12);
    
    % Estrazione stati del filtro (Lead Compensator)
    z_x = [X(13); X(15); X(17)];
    z_y = [X(14); X(16); X(18)];
    
    %  Sensori 
    y1_meas = p1; y2_meas = p2; y3_meas = p3;
    
    if p.USE_NOISE
        noise = p.noise_std * [sin(100*t + 1); cos(110*t + 2)];
        y1_meas = y1_meas + noise;
        y2_meas = y2_meas + noise;
        y3_meas = y3_meas + noise;
    end
    
    if p.USE_DIST
        disturbo = p.dist_amp * [sin(p.dist_freq*t); cos(p.dist_freq*t)];
        y2_meas = y2_meas + disturbo;
    end
    
    % generazione riferimento r(t) 
    [r_t, dr_t] = get_reference(t, p);
    
    % Calcolo errori locali
    eta1 = y1_meas - r_t - (p.H(:,1) - p.H(:,1));
    eta2 = y2_meas - r_t - (p.H(:,2) - p.H(:,1)); 
    eta3 = y3_meas - r_t - (p.H(:,3) - p.H(:,1));
    
    % Derivate degli errori 
    deta1 = v1 - dr_t; deta2 = v2 - dr_t; deta3 = v3 - dr_t;
    
    eta_x = [eta1(1); eta2(1); eta3(1)];
    eta_y = [eta1(2); eta2(2); eta3(2)];
    deta_x = [deta1(1); deta2(1); deta3(1)];
    deta_y = [deta1(2); deta2(2); deta3(2)];
    
    % Errore distribuito 
    epsilon_p_x = -p.Lp * eta_x;
    epsilon_p_y = -p.Lp * eta_y;
    d_epsilon_p_x = -p.Lp * deta_x;
    d_epsilon_p_y = -p.Lp * deta_y;
    
    dz_x = zeros(3,1);
    dz_y = zeros(3,1);

    if p.USE_LEAD
        dz_x = (p.Kd * d_epsilon_p_x - z_x) / p.tau_lead;
        dz_y = (p.Kd * d_epsilon_p_y - z_y) / p.tau_lead;
        
        nu_x = p.Kp * epsilon_p_x + z_x;
        nu_y = p.Kp * epsilon_p_y + z_y;
    else
        nu_x = p.Kp * epsilon_p_x + p.Kd * d_epsilon_p_x;
        nu_y = p.Kp * epsilon_p_y + p.Kd * d_epsilon_p_y;
    end
    
    nu1 = [nu_x(1); nu_y(1)];
    nu2 = [nu_x(2); nu_y(2)];
    nu3 = [nu_x(3); nu_y(3)];
    
    if p.USE_NU_SAT
        nu1 = saturate(nu1, p.NU_MAX);
        nu2 = saturate(nu2, p.NU_MAX);
        nu3 = saturate(nu3, p.NU_MAX);
    end
    
    % FL
    % agente 1
    phi1 = [-0.35*v1(1) - 0.08*v1(1)*abs(v1(1)); ...
            -0.45*v1(2) - 0.10*v1(2)*abs(v1(2))];
    B1 = diag([1/1.2, 1/1.2]);
    u1 = B1 \ (nu1 - phi1);
    
    % angente 2
    phi2 = [-0.40*v2(1); -0.55*v2(2)];
    B2 = diag([1 / (1 + 0.25*sin(p2(1))^2), 1 / (1 + 0.30*cos(p2(2))^2)]);
    u2 = B2 \ (nu2 - phi2);
    
    % agetne 3
    phi3 = [-0.30*v3(1) + 0.15*sin(p3(2))*v3(2); ...
            -0.35*v3(2) + 0.15*cos(p3(1))*v3(1)];
    B3 = [1 + 0.20*cos(p3(2))^2,  0.08*sin(p3(1)); ...
          0.08*sin(p3(2)),        1 + 0.20*sin(p3(1))^2];
    u3 = B3 \ (nu3 - phi3);
    
    % Saturazione input fisico
    if p.USE_U_SAT
        u1 = saturate(u1, p.U_MAX);
        u2 = saturate(u2, p.U_MAX);
        u3 = saturate(u3, p.U_MAX);
    end
    
    dX = zeros(18,1); % Vettore derivate da 18 elementi
    
    % Dinamica fisica 
    dX(1:2)   = v1; dX(3:4)   = phi1 + B1*u1;
    dX(5:6)   = v2; dX(7:8)   = phi2 + B2*u2;
    dX(9:10)  = v3; dX(11:12) = phi3 + B3*u3;
    
    % Dinamica del lead Compensator 
    dX(13) = dz_x(1); dX(14) = dz_y(1);
    dX(15) = dz_x(2); dX(16) = dz_y(2);
    dX(17) = dz_x(3); dX(18) = dz_y(3);
end