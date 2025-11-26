% UNIVERSIDADE FEDERAL DA BAHIA (UFBA)
% DEPARTAMENTO DE ENGENHARIA ELÉTRICA E DA COMPUTAÇÃO (DEEC)
% PROJETO DE CONTROLE DE SISTEMAS
% GABRIEL CORREIA DOS SANTOS

% DESAFIO IV - CONTROLE ADAPTATIVO AUTO-AJUSTÁVEL

%% CONFIGURAÇÃO INICIAL
clear; clc; close all;

%% =============================================================================
% SIMULAÇÃO DE VALIDAÇÃO (QUESTÃO 1, 2 e 3)
% =============================================================================

format shortG

disp('============================================================');
disp('RESOLUÇÃO DAS QUESTÕES 1, 2 E 3');
disp('============================================================');

%% --- QUESTÃO 1: DEFINIÇÃO DO PERÍODO DE AMOSTRAGEM ---
disp(' ');
disp('--- Questão 1: Critério de Amostragem ---');

% Planta Contínua: G(s) = (-0.2s + 1) / (s^2 + 2s + 1)
Gs = tf([-0.2, 1], [1, 2, 1]);

% Step Info (SettlingTime padrão do MATLAB é 2%, compatível com Python)
SI = stepinfo(Gs, 'SettlingTimeThreshold', 0.02);

% Cálculo do Undershoot (MATLAB retorna positivo, Python tbm)
% PeakTime e outros já vêm na struct
ts_ma = SI.SettlingTime;

% Exibição formatada igual ao Python
disp('SI = struct with fields:');
disp(' ');
fprintf('         RiseTime: %.4f\n', SI.RiseTime);
fprintf('     SettlingTime: %.4f\n', SI.SettlingTime);
fprintf('      SettlingMin: %.4f\n', SI.SettlingMin);
fprintf('      SettlingMax: %.4f\n', SI.SettlingMax);
fprintf('        Overshoot: %.4f\n', SI.Overshoot);
fprintf('       Undershoot: %.4f\n', SI.Undershoot);
fprintf('             Peak: %.4f\n', SI.Peak);
fprintf('         PeakTime: %.4f\n', SI.PeakTime);

ts_mf_desejado = ts_ma / 2;
Ta = ts_mf_desejado / 10;

disp(' ');
fprintf('Tempo de Acomodação (Malha Aberta): %.4f s\n', ts_ma);
fprintf('Tempo de Acomodação Desejado (MF):  %.4f s\n', ts_mf_desejado);
fprintf('Período de Amostragem Escolhido (Ta): %.4f s\n', Ta);

%% --- QUESTÃO 2: DISCRETIZAÇÃO DA PLANTA ---
disp(' ');
disp('--- Questão 2: Modelo Discreto G(z) ---');

Gz = c2d(Gs, Ta, 'tustin');
[num_z, den_z] = tfdata(Gz, 'v');

% Normalização (garantir den(1) = 1)
num_z = num_z / den_z(1);
den_z = den_z / den_z(1);

A_coeffs = den_z;       % [1, a1, a0]
B_coeffs = num_z(2:end); % [b0, b1] (Ignora z^0 do Tustin)

disp('Pd =');
Gz
fprintf('Sample time: %.5f seconds\n', Ta);

disp(' ');
disp('Model Properties');
FZ = tf([1 0], [1], Ta, 'Variable', 'z');
disp('FZ =');
disp(FZ);

%% --- QUESTÃO 3: PROJETO DO CONTROLADOR ---
disp('--- Questão 3: Controlador RST por Alocação de Polos ---');

% 1. Especificações
xi = 0.7;
wn = 4 / ts_mf_desejado;

% 2. Modelo de Referência
Hs = tf([wn^2], [1, 2*xi*wn, wn^2]);
disp(' ');
disp('H =');
disp(Hs);

Hz = c2d(Hs, Ta, 'tustin');
disp('Hd =');
disp(Hz);

[num_hz, den_hz] = tfdata(Hz, 'v');
Am_coeffs = den_hz;
B_ref_coeffs = num_hz;

% Prints dos vetores de coeficientes
disp('nmf ='); disp(B_ref_coeffs);
disp('dmf ='); disp(Am_coeffs);
disp('nma ='); disp(B_coeffs);
disp('dma ='); disp(A_coeffs);

% 3. Observador e Polinômio Desejado
p_obs = 0.8;
Ao_coeffs = conv([1, -p_obs], [1, -p_obs]);
pol_des = conv(Am_coeffs, Ao_coeffs);

disp('pol_des ='); disp(pol_des);

% --- MATRIZ DE SYLVESTER (Construção para exibição das colunas) ---
% A_aug = A * (z-1)
A_aug_disp = conv(A_coeffs, [1, -1]);
% B aumentado = [0, b0, b1]
B_aug_disp = [0, B_coeffs];

% Montagem manual das colunas para printar igual ao Python
col1 = [A_aug_disp, 0];
col2 = [0, A_aug_disp];
col3 = [B_aug_disp, 0, 0];
col4 = [0, B_aug_disp, 0];
col5 = [0, 0, B_aug_disp];

disp('col1 ='); disp(col1);
disp('col2 ='); disp(col2);
disp('col3 ='); disp(col3);
disp('col4 ='); disp(col4);
disp('col5 ='); disp(col5);

% Resolução Oficial
[R_poly, S_poly, T_poly] = projetar_controlador_rst(A_coeffs, B_coeffs, Am_coeffs, Ao_coeffs);

% Solução (Vetor x da equação M*x = P)
% x = [r0', r1', s0, s1, s2]
% Note: R_poly já é convolucionado, então precisamos recuperar r' para mostrar no 'Sol'
% Mas podemos montar o vetor Sol com os dados que temos:
% R' = R / (z-1). Como R = conv(r', [1, -1]), fazemos deconv
[r_prime, ~] = deconv(R_poly, [1, -1]);
Sol = [r_prime(:); S_poly(:)];

disp('Sol ='); disp(Sol);

% Exibição dos Polinômios Finais e TFs
disp(' ');
disp('S =');
disp(tf(S_poly, [1], Ta, 'Variable', 'z'));

disp('R =');
disp(tf(R_poly, [1], Ta, 'Variable', 'z'));

disp('To =');
disp(tf(Ao_coeffs, [1], Ta, 'Variable', 'z'));

disp('T =');
disp(tf(T_poly, [1], Ta, 'Variable', 'z'));

% Vetores C e F
disp('C =');
fprintf('Num: '); disp(S_poly);
fprintf('Den: '); disp(R_poly);

disp(' ');
disp('F =');
fprintf('Num: '); disp(T_poly);
fprintf('Den: '); disp(S_poly);


%% =============================================================================
% SIMULAÇÃO DE VALIDAÇÃO (QUESTÃO 3)
% =============================================================================
disp('Simulando resposta em malha fechada...');

t_final = 50;
steps = floor(t_final / Ta);
time = linspace(0, t_final, steps);

y_sim = zeros(1, steps);
u_sim = zeros(1, steps);
r_sim = zeros(1, steps);
qy_sim = zeros(1, steps);
qu_sim = zeros(1, steps);

u_buf = zeros(1, 10);
y_buf = zeros(1, 10);

[Ap, Bp, Cp, Dp] = tf2ss(cell2mat(Gs.num), cell2mat(Gs.den));
x_p = [0; 0];

for k = 1:steps
    t = (k-1) * Ta;
    r_val = 1.0;
    r_sim(k) = r_val;

    % Perturbações
    dy = 0; if t >= 15, dy = 0.2; end
    du = 0; if t >= 35, du = 0.2; end
    qy_sim(k) = dy; qu_sim(k) = du;

    % Controle
    tr_term = sum(T_poly * r_val);
    sy_term = 0;
    for i = 1:length(S_poly)
        if i <= length(y_buf), sy_term = sy_term + S_poly(i) * y_buf(i); end
    end
    ru_past = 0;
    for i = 2:length(R_poly)
        if (i-1) <= length(u_buf), ru_past = ru_past + R_poly(i) * u_buf(i-1); end
    end
    
    u_val = (tr_term - sy_term - ru_past) / R_poly(1);
    u_sim(k) = u_val;
    u_buf = [u_val, u_buf(1:end-1)];

    % Planta
    u_app = u_val + du;
    n_sub = 20; dt_sub = Ta/n_sub;
    for s = 1:n_sub
        dx = Ap * x_p + Bp * u_app;
        x_p = x_p + dx * dt_sub;
    end
    y_val = (Cp * x_p) + Dp * u_app + dy;
    y_sim(k) = y_val;
    y_buf = [y_val, y_buf(1:end-1)];
end

% PLOTAGEM FIGURA 1
figure('Name', 'Questão 3', 'Color', 'w');
subplot(2,1,1);
plot(time, y_sim, 'b', 'LineWidth', 2); hold on;
plot(time, r_sim, 'k--');
grid on; title('Questões 1-3: Controle Fixo por Alocação de Polos');
ylabel('Amplitude'); legend('Saída y(t)', 'Referência');

subplot(2,1,2);
plot(time, u_sim, 'r', 'LineWidth', 1.5);
grid on; ylabel('Sinal de Controle'); xlabel('Tempo (s)'); legend('Controle u(t)');

% PLOTAGEM FIGURA 2
figure('Name', 'Perturbações', 'Color', 'w');
subplot(2,1,1);
plot(time, qy_sim, 'b', 'LineWidth', 1.5);
grid on; title('Perfil da Perturbação de Saída (t=15s)'); ylabel('Perturbação Dy'); ylim([-0.05, 0.3]);

subplot(2,1,2);
plot(time, qu_sim, 'b', 'LineWidth', 1.5);
grid on; title('Perfil da Perturbação de Entrada (t=35s)'); ylabel('Perturbação Du'); xlabel('Tempo (s)'); ylim([-0.05, 0.3]);


%% =============================================================================
% 3. EXECUÇÃO DA QUESTÃO 4
% =============================================================================
disp(' ');
disp('============================================================');
disp('QUESTÃO 4: Adaptação sem Perturbação');
disp('============================================================');

[t1, y1, u1, th1, true_th] = run_adaptive_simulation(4, 1.05, 'positional', Gs, Ta, Am_coeffs, Ao_coeffs);
[t2, y2, u2, th2, ~]       = run_adaptive_simulation(4, 0.95, 'positional', Gs, Ta, Am_coeffs, Ao_coeffs);

figure('Name', 'Questão 4', 'Color', 'w');
subplot(2, 1, 1);
plot(t1, y1, 'b', 'LineWidth', 1.5); hold on;
plot(t2, y2, 'r--', 'LineWidth', 1.5);
plot(t1, ones(size(t1)), 'k:');
title('Questão 4: Resposta Adaptativa (Posicional)');
legend('Inic 1.05x', 'Inic 0.95x', 'Referência'); grid on; ylabel('Saída');

subplot(2, 1, 2);
plot(t1, th1(:, 1), 'b', 'LineWidth', 1.5); hold on;
plot(t2, th2(:, 1), 'r--', 'LineWidth', 1.5);
yline(true_th(1), 'k:', 'Real');
title('Convergência do Parâmetro a1');
legend('Est. (1.05x)', 'Est. (0.95x)', 'Real'); grid on; ylabel('Valor de a1');

%disp('OBSERVAÇÕES Q4:');
%disp('- O sistema converge para a referência independentemente do erro inicial.');

%% =============================================================================
% 4. EXECUÇÃO DA QUESTÃO 5
% =============================================================================
disp(' ');
disp('============================================================');
disp('QUESTÃO 5: Rejeição de Perturbação e Regressores');
disp('============================================================');

[t_pos, y_pos, u_pos, th_pos, ~] = run_adaptive_simulation(5, 1.05, 'positional', Gs, Ta, Am_coeffs, Ao_coeffs);
[t_inc, y_inc, u_inc, th_inc, ~] = run_adaptive_simulation(5, 1.05, 'incremental', Gs, Ta, Am_coeffs, Ao_coeffs);

figure('Name', 'Questão 5', 'Color', 'w');
subplot(2, 1, 1);
plot(t_pos, y_pos, 'b', 'LineWidth', 1.5); hold on;
plot(t_inc, y_inc, 'g--', 'LineWidth', 1.5);
plot(t_pos, ones(size(t_pos)), 'k:');
title('Questão 5: Rejeição de Perturbação (t=35s)');
legend('Posicional', 'Incremental', 'Referência'); grid on; ylabel('Saída');

subplot(2, 1, 2);
plot(t_pos, th_pos(:, 3), 'b', 'LineWidth', 1.5); hold on;
plot(t_inc, th_inc(:, 3), 'g--', 'LineWidth', 1.5);
xline(35, 'k:', 'Perturbação');
title('Estimativa de b0 durante Perturbação');
legend('Posicional (Viés)', 'Incremental (Robusto)', 'Perturbação'); grid on; ylabel('Valor de b0');

%disp('OBSERVAÇÕES Q5:');
%disp('- Ambos rejeitam perturbação devido à ação integral.');
%disp('- Posicional sofre viés (drift). Incremental é robusto.');


%% =============================================================================
% FUNÇÕES LOCAIS
% =============================================================================

function [R, S, T] = projetar_controlador_rst(A_poly, B_poly, Am_poly, Ao_poly)
    A_aug = conv(A_poly, [1, -1]);
    P_cl = conv(Am_poly, Ao_poly);
    
    bb = [0, B_poly(:)']; 
    aa = A_aug(:)';
    
    M = zeros(5, 5);
    M(1:4, 1) = aa; M(2:5, 2) = aa;
    M(1:3, 3) = bb; M(2:4, 4) = bb; M(3:5, 5) = bb;
    
    rhs = zeros(5, 1);
    rhs(1:length(P_cl)) = P_cl(:);
    
    if rank(M) < 5, R=[]; S=[]; T=[]; return; end
    
    x = M \ rhs;
    r_prime = x(1:2)'; s = x(3:5)';
    
    R = conv(r_prime, [1, -1]);
    S = s;
    
    gain_B = sum(B_poly); gain_Am = sum(Am_poly);
    if abs(gain_B) < 1e-6, gain_B = 1e-6; end
    T = Ao_poly * (gain_Am / gain_B);
end

function [theta, P] = rls_update(theta, P, target_val, phi, lambda_factor)
    phi = phi(:);
    y_hat = phi' * theta;
    epsilon = target_val - y_hat;
    denom = lambda_factor + phi' * P * phi;
    K = (P * phi) / denom;
    theta = theta + K * epsilon;
    P = (P - K * phi' * P) / lambda_factor;
end

function [time, y_hist, u_hist, theta_hist, theta_true] = run_adaptive_simulation(q_num, init_scale, reg_type, Gs, Ta, Am, Ao)
    fprintf('Iniciando: Q%d, Escala=%.2f, Regressor=%s...\n', q_num, init_scale, reg_type);
    
    Gz = c2d(Gs, Ta, 'tustin');
    [num, den] = tfdata(Gz, 'v');
    num = num / den(1); den = den / den(1);
    
    a_true = den(2:end); b_true = num(2:end);
    theta_true = [a_true(:); b_true(:)];
    
    theta_est = theta_true * init_scale;
    P = eye(4) * 1000;
    
    t_final = 50; steps = floor(t_final / Ta);
    time = linspace(0, t_final, steps);
    
    y_hist = zeros(1, steps); u_hist = zeros(1, steps); theta_hist = zeros(steps, 4);
    [Ap, Bp, Cp, Dp] = tf2ss(cell2mat(Gs.num), cell2mat(Gs.den));
    x_plant = [0; 0];
    u_buf = zeros(1, 10); y_buf = zeros(1, 10);
    R_poly = [1]; S_poly = [0]; T_poly = [1];
    
    for k = 1:steps
        t = (k-1) * Ta;
        di = 0; if q_num == 5 && t >= 35, di = 0.2; end
        
        u_applied = u_buf(1) + di;
        for s=1:10, dx = Ap*x_plant + Bp*u_applied; x_plant = x_plant + dx*(Ta/10); end
        y_val = (Cp * x_plant) + Dp * u_applied;
        y_hist(k) = y_val; y_buf = [y_val, y_buf(1:end-1)];
        
        if k > 4
            if strcmp(reg_type, 'positional')
                phi = [-y_buf(2); -y_buf(3); u_buf(2); u_buf(3)]; target = y_val;
            else
                phi = [-(y_buf(2)-y_buf(3)); -(y_buf(3)-y_buf(4)); (u_buf(2)-u_buf(3)); (u_buf(3)-u_buf(4))];
                target = y_buf(1) - y_buf(2);
            end
            [theta_est, P] = rls_update(theta_est, P, target, phi, 1.0);
            [r_new, s_new, t_new] = projetar_controlador_rst([1; theta_est(1:2)], theta_est(3:4), Am, Ao);
            if ~isempty(r_new), R_poly = r_new; S_poly = s_new; T_poly = t_new; end
        end
        theta_hist(k, :) = theta_est';
        
        tr = sum(T_poly); sy = 0;
        for i = 1:length(S_poly), if i<=length(y_buf), sy=sy+S_poly(i)*y_buf(i); end, end
        ru = 0; for i = 2:length(R_poly), if (i-1)<=length(u_buf), ru=ru+R_poly(i)*u_buf(i-1); end, end
        u_val = (tr - sy - ru) / R_poly(1);
        if u_val > 10, u_val=10; elseif u_val < -10, u_val=-10; end
        u_hist(k) = u_val; u_buf = [u_val, u_buf(1:end-1)];
    end
end

% Mil milhão de anos depois...