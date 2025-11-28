% UNIVERSIDADE FEDERAL DA BAHIA (UFBA)
% DEPARTAMENTO DE ENGENHARIA ELÉTRICA E DA COMPUTAÇÃO (DEEC)
% PROJETO DE CONTROLE DE SISTEMAS
% GABRIEL CORREIA DOS SANTOS

% DESAFIO III - SISTEMAS DE CONTROLE AMOSTRADOS

%% CONFIGURAÇÃO INICIAL
clear; clc; close all;
clear all; close all; clc;

% =========================================================================
% Sintonia:
% - Alocação de Polos: Tau = 0.8s (Equilibrado)
% - IMC: Tau = 0.5s (Rápido)
% =========================================================================

% =========================================================================
% 1. SISTEMA E AMOSTRAGEM
% =========================================================================
s = tf('s');

% Planta Original
P = (0.2 * (10 - s)) / ((s + 1)^2);

% Cálculo do Tempo de Amostragem (Ta)
info = stepinfo(P);
ts_ma = info.SettlingTime;
Ta = (ts_ma / 2) / 10;

fprintf('Ta calculado: %.4fs\n', Ta);

% =========================================================================
% 2. PREPARAÇÃO PARA O PROJETO
% =========================================================================
z = tf('z', Ta);

% Definições auxiliares para a fórmula do IMC
% s via Tustin
s_tustin = (2/Ta) * (z - 1) / (z + 1);

% Parte Não-Invertível (Zero Instável)
B_menosz_tf = 10 - s_tustin;
B_gain = dcgain(B_menosz_tf);

% Inversa da parte invertível (Pn_tilde^-1)
Pn_til_s = 0.2 / ((s + 1)^2);
Pn_til_z = c2d(Pn_til_s, Ta, 'tustin');
Pn_til_inv = 1 / Pn_til_z;

% Planta Completa Discretizada (para simulação)
P_z = c2d(P, Ta, 'tustin');

% =========================================================================
% 3. PROJETO DOS CONTROLADORES
% =========================================================================

% --- A. Alocação de Polos (Tau = 0.8s) ---
tau_pp = 0.8;
[R_pp, S_pp, T_pp, C_pp_tf] = calcular_rst_final(tau_pp, s, Ta, B_gain, Pn_til_inv, B_menosz_tf);

fprintf('\n--- [Q3] Alocação de Polos (Tau=%.1fs) ---\n', tau_pp);
fprintf('R: %s\n', mat2str(R_pp, 4));
fprintf('S: %s\n', mat2str(S_pp, 4));
fprintf('T: %s\n', mat2str(T_pp, 4));
fprintf('Função de Transferência C(z):\n');
display(C_pp_tf);

% --- B. IMC (Tau = 0.5s) ---
tau_imc = 0.5;
[R_imc, S_imc, T_imc, C_imc_tf] = calcular_rst_final(tau_imc, s, Ta, B_gain, Pn_til_inv, B_menosz_tf);

fprintf('\n--- [Q4] IMC (Tau=%.1fs) ---\n', tau_imc);
fprintf('R: %s\n', mat2str(R_imc, 4));
fprintf('S: %s\n', mat2str(S_imc, 4));
fprintf('T: %s\n', mat2str(T_imc, 4));
fprintf('Função de Transferência C(z):\n');
display(C_imc_tf);

% =========================================================================
% 4. SIMULAÇÃO
% =========================================================================
tsim = 50;
steps = floor(tsim / Ta);
t_vec = (0:steps-1) * Ta;

% Sinais
ref = ones(1, steps);
ref(1:floor(1/Ta)) = 0; 

du = zeros(1, steps);
du(floor(35/Ta):end) = 0.2; 

dy = zeros(1, steps);
dy(floor(15/Ta):end) = 0.2; 

% Simulação
[y_pp, u_pp] = simular_loop(R_pp, S_pp, T_pp, P_z, steps, ref, du, dy);
[y_imc, u_imc] = simular_loop(R_imc, S_imc, T_imc, P_z, steps, ref, du, dy);

% =========================================================================
% 5. PLOTAGEM
% =========================================================================

% --- FIGURA 1: COMPARATIVA (Q5) ---
figure(1);
subplot(2,1,1);
plot(t_vec, ref, 'k--', 'LineWidth', 1); hold on;
plot(t_vec, y_pp, 'r-', 'LineWidth', 1.5);
plot(t_vec, y_imc, 'g-.', 'LineWidth', 1.5);
xline(15, ':', 'Color', [0.5 0.5 0.5]);
xline(35, ':', 'Color', [0.5 0.5 0.5]);
title('Comparação Final [Q5]: Alocação vs IMC');
ylabel('Saída y(t)');
legend('Referência', ['Alocação (\tau=' num2str(tau_pp) 's)'], ['IMC (\tau=' num2str(tau_imc) 's)'], 'Location', 'best');
ylim([0 1.5]); % só coloquei isso aqui para os gráficos ficarem todos iguais 
grid on;

subplot(2,1,2);
plot(t_vec, u_pp, 'r-', 'LineWidth', 1.5); hold on;
plot(t_vec, u_imc, 'g-.', 'LineWidth', 1.5);
ylabel('Controle u(t)');
xlabel('Tempo (s)');
legend('Controle PP', 'Controle IMC');
grid on;

% --- FIGURA 2: DETALHE ALOCAÇÃO (Q2/Q3) ---
figure(2);
subplot(2,1,1);
plot(t_vec, ref, 'k--', 'LineWidth', 1); hold on;
plot(t_vec, y_pp, 'r-', 'LineWidth', 1.5);
xline(15, ':', 'Color', [0.5 0.5 0.5]);
text(15.5, 0.5, 'Pert. Saída', 'Color', [0.4 0.4 0.4]);
xline(35, ':', 'Color', [0.5 0.5 0.5]);
text(35.5, 0.5, 'Pert. Entrada', 'Color', [0.4 0.4 0.4]);
title(['Detalhe [Q2/Q3]: Alocação de Polos (\tau=' num2str(tau_pp) 's)']);
ylabel('Saída y(t)');
legend('Referência', 'Saída Controlada');
ylim([0 1.5]); 
grid on;

subplot(2,1,2);
plot(t_vec, u_pp, 'r-', 'LineWidth', 1.5);
ylabel('Controle u(t)');
xlabel('Tempo (s)');
legend('Sinal de Controle u(t)');
grid on;

% --- FIGURA 3: DETALHE IMC (Q4) ---
figure(3);
subplot(2,1,1);
plot(t_vec, ref, 'k--', 'LineWidth', 1); hold on;
plot(t_vec, y_imc, 'g-.', 'LineWidth', 1.5);
xline(15, ':', 'Color', [0.5 0.5 0.5]);
text(15.5, 0.5, 'Pert. Saída', 'Color', [0.4 0.4 0.4]);
xline(35, ':', 'Color', [0.5 0.5 0.5]);
text(35.5, 0.5, 'Pert. Entrada', 'Color', [0.4 0.4 0.4]);
title(['Detalhe [Q4]: IMC (\tau=' num2str(tau_imc) 's)']);
ylabel('Saída y(t)');
legend('Referência', 'Saída Controlada');
ylim([0 1.5]); % só coloquei isso aqui para os gráficos ficarem todos iguais 
grid on;

subplot(2,1,2);
plot(t_vec, u_imc, 'g-.', 'LineWidth', 1.5);
ylabel('Controle u(t)');
xlabel('Tempo (s)');
legend('Sinal de Controle u(t)');
grid on;

% --- FIGURA 4: LUGAR DAS RAÍZES ---
figure(4);
subplot(1,2,1);
rlocus(P_z * C_pp_tf);
zgrid; title(['Lugar das Raízes: Alocação (\tau=' num2str(tau_pp) 's)']); axis equal;

subplot(1,2,2);
rlocus(P_z * C_imc_tf);
zgrid; title(['Lugar das Raízes: IMC (\tau=' num2str(tau_imc) 's)']); axis equal;

fprintf('\nSimulação concluída.\n');

% =========================================================================
% FUNÇÕES AUXILIARES
% =========================================================================

function [R, S, T, C_tf] = calcular_rst_final(tau_desejado, s, Ta, B_gain, Pn_til_inv, B_menosz_tf)
    % Filtro F(z)
    F_s = 1 / ((tau_desejado * s + 1)^2);
    F_z = c2d(F_s, Ta, 'tustin');
    
    Frz = F_z / B_gain;
    
    % Cálculo do Controlador
    % TOLERÂNCIA DE 0.1 -> Isso força o cancelamento dos polos parasitas em z=-1
    C_tf = minreal((Frz * Pn_til_inv) / (1 - Frz * B_menosz_tf), 0.1);
    
    % Extração dos Polinômios
    [num, den] = tfdata(C_tf, 'v');
    
    % Normalização
    norm = den(1);
    R = den / norm;
    S = num / norm;
    T = S;
end

function [y, u] = simular_loop(R, S, T, P_z, steps, ref, du, dy)
    [Bn, An] = tfdata(P_z, 'v');
    norm_p = An(1);
    Bn = Bn / norm_p;
    An = An / norm_p;
    
    An_rest = An(2:end); 
    Rm = R(2:end);
    
    y = zeros(1, steps);
    u = zeros(1, steps);
    
    u_hist = zeros(length(Rm), 1);
    y_hist = zeros(length(S), 1);
    r_hist = zeros(length(T), 1);
    
    u_plant_hist = zeros(length(Bn), 1);
    y_plant_hist = zeros(length(An_rest), 1);
    
    for k = 1:steps
        r_hist = [ref(k); r_hist(1:end-1)];
        
        u_val = (T * r_hist) - (S * y_hist) - (Rm * u_hist);
        u(k) = u_val;
        
        u_real = u_val + du(k);
        u_plant_hist = [u_real; u_plant_hist(1:end-1)];
        
        y_clean = (Bn * u_plant_hist) - (An_rest * y_plant_hist);
        y(k) = y_clean + dy(k);
        
        u_hist = [u(k); u_hist(1:end-1)];
        y_hist = [y(k); y_hist(1:end-1)];
        y_plant_hist = [y_clean; y_plant_hist(1:end-1)];
    end
end

% Mil milhão de anos depois...