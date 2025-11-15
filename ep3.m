%% Script para PTC3314 - Exercício de Simulação 3
% Baseado no NUSP 12544308

clear;
clc;
close all;

%% ------------------- Constantes e Definições -------------------
nusp = 12544308;
f0 = 1.308e9; % Frequência central f0 = 1,308 GHz

% Definição da faixa de frequência da simulação
f_start = f0 - 500e6; % f0 - 500 MHz
f_stop = f0 + 500e6;  % f0 + 500 MHz
passo = 0.5e6;       % Passo de 0.5 MHz
f = f_start:passo:f_stop; % Vetor de frequências

% Impedâncias (usando valores relativos a eta0)
eta_1 = 1;     % Ar
eta_f = 1/4;   % Dielétrico final (sqrt(16) = 4)

% Limite para cálculo da Largura de Banda (BW)
BW_limit = 0.001; %

% Pre-alocação de vetores para os resultados
rho1_sq_q1 = zeros(size(f));
rho1_sq_q2 = zeros(size(f));

%% ------------------- Loop de Simulação -------------------
% A função de impedância de entrada é:
% Z_in = Z0 * (Z_L + 1j * Z0 * tan(beta*d)) / (Z0 + 1j * Z_L * tan(beta*d))
% Onde beta*d = (2*pi*f/v) * (v / (4*f0)) = (pi * f) / (2 * f0)

for i = 1:length(f)
    freq_atual = f(i);
    
    % Argumento da tangente (comprimento elétrico)
    % beta*d = (pi * f) / (2 * f0)
    beta_d = (pi * freq_atual) / (2 * f0);
    
    % --- Questão 1: Uma camada ---
    % eta2 = sqrt(eta1 * etaf) = sqrt(1 * 1/4) = 1/2
    eta_2_q1 = 1/2; 
    
    % Impedância de entrada vista da interface 1
    eta_in_q1 = eta_2_q1 * (eta_f + 1j * eta_2_q1 * tan(beta_d)) / ...
                          (eta_2_q1 + 1j * eta_f * tan(beta_d));
    
    % Coeficiente de reflexão na interface 1
    rho1_q1 = (eta_in_q1 - eta_1) / (eta_in_q1 + eta_1);
    rho1_sq_q1(i) = abs(rho1_q1)^2;
    
    
    % --- Questão 2: Duas camadas ---
    eta_2_q2 = 1/sqrt(2); %
    eta_3_q2 = 1/sqrt(8); %
    
    % Passo 1: Impedância de entrada vista da interface 2 (entrada da camada 3)
    eta_in_2 = eta_3_q2 * (eta_f + 1j * eta_3_q2 * tan(beta_d)) / ...
                           (eta_3_q2 + 1j * eta_f * tan(beta_d));
    
    % Passo 2: Impedância de entrada vista da interface 1 (entrada da camada 2)
    % A carga para a camada 2 é eta_in_2
    eta_in_1 = eta_2_q2 * (eta_in_2 + 1j * eta_2_q2 * tan(beta_d)) / ...
                           (eta_2_q2 + 1j * eta_in_2 * tan(beta_d));
                           
    % Coeficiente de reflexão na interface 1
    rho1_q2 = (eta_in_1 - eta_1) / (eta_in_1 + eta_1);
    rho1_sq_q2(i) = abs(rho1_q2)^2;
end

%% ------------------- Resultados e Gráficos -------------------

% --- Questão 1 ---
fprintf('--------- Questão 1 (f0 = %.3f GHz) ---------\n', f0/1e9);
% Valores específicos
f0_idx = find(f == f0);
fprintf('  |rho1(f0 - 500 MHz)|^2 = %.5f\n', rho1_sq_q1(1));
fprintf('  |rho1(f0)|^2           = %.5f (ideal: 0)\n', rho1_sq_q1(f0_idx));
fprintf('  |rho1(f0 + 500 MHz)|^2 = %.5f\n', rho1_sq_q1(end));

% Cálculo da Largura de Banda (BW)
indices_bw_q1 = find(rho1_sq_q1 <= BW_limit);
f_low_q1 = f(indices_bw_q1(1)) / 1e6;   % em MHz
f_high_q1 = f(indices_bw_q1(end)) / 1e6; % em MHz
BW_q1 = f_high_q1 - f_low_q1;
fprintf('  Largura de Banda (<= %.3f): %.1f MHz (de %.1f a %.1f MHz)\n', ...
         BW_limit, BW_q1, f_low_q1, f_high_q1);

% Gráfico 1
figure;
plot(f / 1e6, rho1_sq_q1, 'b', 'LineWidth', 1.5);
hold on;
% Linhas de referência
plot([f_low_q1, f_high_q1], [BW_limit, BW_limit], 'r--', 'LineWidth', 1);
title(sprintf('Questão 1: Casador de 1 Camada (f_0 = %.3f GHz)', f0/1e9));
xlabel('Frequência (MHz)');
ylabel('|\rho_1|^2');
grid on;
ylim([0 0.18]); % Ajuste do eixo Y para melhor visualização
legend('|\rho_1|^2', sprintf('Limite BW (%.3f)', BW_limit), 'Location', 'north');

% --- Questão 2 ---
fprintf('\n--------- Questão 2 (f0 = %.3f GHz) ---------\n', f0/1e9);
% Valores específicos
fprintf('  |rho1(f0 - 500 MHz)|^2 = %.5f\n', rho1_sq_q2(1));
fprintf('  |rho1(f0)|^2           = %.5f (ideal: 0)\n', rho1_sq_q2(f0_idx));
fprintf('  |rho1(f0 + 500 MHz)|^2 = %.5f\n', rho1_sq_q2(end));

% Cálculo da Largura de Banda (BW)
indices_bw_q2 = find(rho1_sq_q2 <= BW_limit);
f_low_q2 = f(indices_bw_q2(1)) / 1e6;   % em MHz
f_high_q2 = f(indices_bw_q2(end)) / 1e6; % em MHz
BW_q2 = f_high_q2 - f_low_q2;
fprintf('  Largura de Banda (<= %.3f): %.1f MHz (de %.1f a %.1f MHz)\n', ...
         BW_limit, BW_q2, f_low_q2, f_high_q2);

% Gráfico 2
figure;
plot(f / 1e6, rho1_sq_q2, 'b', 'LineWidth', 1.5);
hold on;
% Linhas de referência
plot([f_low_q2, f_high_q2], [BW_limit, BW_limit], 'r--', 'LineWidth', 1);
title(sprintf('Questão 2: Casador de 2 Camadas (f_0 = %.3f GHz)', f0/1e9));
xlabel('Frequência (MHz)');
ylabel('|\rho_1|^2');
grid on;
ylim([0 0.06]); % Ajuste do eixo Y (note que o máximo é bem menor)
legend('|\rho_1|^2', sprintf('Limite BW (%.3f)', BW_limit), 'Location', 'north');

fprintf('\nAnálise: BW Q1 = %.1f MHz, BW Q2 = %.1f MHz.\n', BW_q1, BW_q2);
fprintf('O casador de 2 camadas (binomial) tem uma largura de banda %.2f vezes maior.\n', BW_q2 / BW_q1);