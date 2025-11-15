%% Script para PTC3314 - Exercício de Simulação 4
% Baseado no NUSP 12544308
% MODIFICADO: Eixo Y do Gráfico (a) ajustado para [0 100]

clear;
clc;
close all;

%% ------------------- 1. Constantes e Parâmetros -------------------
% Parâmetros do NUSP
nusp = 12544308;
mnp = mod(nusp, 1000); % 308
epsilon_r1 = 3 + mnp / 1000.0; % 3.308 (Vidro)
epsilon_r2 = 1.0; % Ar

% Constantes Físicas
f = 1.0e9; % 1.0 GHz
c = 3.0e8; % Velocidade da luz (m/s)
eta_0 = 376.73; % Impedância do vácuo (Ohms)
H_y1_plus = 1.0e-3; % 1.0 mA_ef/m

% Parâmetros dos Meios
eta_1 = eta_0 / sqrt(epsilon_r1); % Impedância do vidro
eta_2 = eta_0 / sqrt(epsilon_r2); % Impedância do ar
omega = 2 * pi * f;
k1 = omega / c * sqrt(epsilon_r1); % Número de onda no vidro
k2 = omega / c * sqrt(epsilon_r2); % Número de onda no ar

% Vetor de Ângulos
passo = 0.1;
theta_1_deg = (0:passo:90)';
theta_1_rad = deg2rad(theta_1_deg);

fprintf('Parâmetros Iniciais (NUSP: %d)\n', nusp);
fprintf('Epsilon_r1 (vidro) = %.3f\n', epsilon_r1);
fprintf('eta_1 = %.2f Ohms\n', eta_1);
fprintf('eta_2 = %.2f Ohms\n', eta_2);

%% ------------------- a) Ângulo de Transmissão -------------------

% Cálculo do Ângulo Crítico
sin_theta_c = sqrt(epsilon_r2 / epsilon_r1);
theta_c_rad = asin(sin_theta_c);
theta_c_deg = rad2deg(theta_c_rad);

fprintf('\n--------- Questão (a) ---------\n');
fprintf('Ângulo Crítico (theta_c): %.2f graus\n', theta_c_deg);

% Cálculo de theta_2 (Lei de Snell)
sin_theta_2 = sqrt(epsilon_r1 / epsilon_r2) * sin(theta_1_rad);

% Para theta_1 > theta_c, a parte real de theta_2 é 90 graus
theta_2_rad = zeros(size(theta_1_rad));
idx_pre_c = theta_1_rad <= theta_c_rad;
idx_post_c = theta_1_rad > theta_c_rad;

theta_2_rad(idx_pre_c) = asin(sin_theta_2(idx_pre_c));
theta_2_rad(idx_post_c) = pi / 2; % Parte real
theta_2_deg = rad2deg(theta_2_rad);

% Gráfico (a)
figure;
plot(theta_1_deg, theta_2_deg, 'b', 'LineWidth', 1.5);
hold on;
xline(theta_c_deg, 'r--', 'LineWidth', 1.5);
legend('Ângulo Transmitido $\theta_2$', sprintf('Ângulo Crítico $\\theta_c = %.2f^\\circ$', theta_c_deg), 'Interpreter', 'latex');
title('Gráfico (a): $\theta_2$ vs $\theta_1$', 'Interpreter', 'latex');
xlabel('Ângulo de Incidência $\theta_1$ (graus)', 'Interpreter', 'latex');
ylabel('Ângulo de Transmissão $\theta_2$ (graus)', 'Interpreter', 'latex');
ylim([0 100]); % <--- MODIFICAÇÃO SOLICITADA
grid on;

%% ------------------- b) Impedância de Carga Z_L -------------------

% Z_L = eta_2 * cos(theta_2) (Polarização TM / E no plano)
Z_L = zeros(size(theta_1_rad), 'like', 1i); % Vetor complexo

% Antes de theta_c (Z_L é real)
Z_L(idx_pre_c) = eta_2 * cos(theta_2_rad(idx_pre_c));

% Depois de theta_c (Z_L é puramente imaginário)
% cos(theta_2) = -j * sqrt( (e1/e2)sin^2(t1) - 1 )
argumento_sqrt = (epsilon_r1 / epsilon_r2) * sin(theta_1_rad(idx_post_c)).^2 - 1;
Z_L(idx_post_c) = -1i * eta_2 * sqrt(argumento_sqrt);

% Valores Específicos
Z_L_0 = Z_L(1);
[~, idx_c] = min(abs(theta_1_deg - theta_c_deg));
Z_L_theta_c = Z_L(idx_c);
Z_L_90 = Z_L(end);

fprintf('\n--------- Questão (b) ---------\n');
fprintf('Z_L(0 deg)   = %.2f + j(%.2f) Ohms\n', real(Z_L_0), imag(Z_L_0));
fprintf('Z_L(theta_c) = %.2f + j(%.2f) Ohms\n', real(Z_L_theta_c), imag(Z_L_theta_c));
fprintf('Z_L(90 deg)  = %.2f + j(%.2f) Ohms\n', real(Z_L_90), imag(Z_L_90));

% Gráfico (b)
figure;
plot(theta_1_deg, real(Z_L), 'b', 'LineWidth', 1.5);
hold on;
plot(theta_1_deg, imag(Z_L), 'g', 'LineWidth', 1.5);
xline(theta_c_deg, 'r--', 'LineWidth', 1.5);
title('Gráfico (b): Impedância de Carga $Z_L$ vs $\theta_1$', 'Interpreter', 'latex');
xlabel('Ângulo de Incidência $\theta_1$ (graus)', 'Interpreter', 'latex');
ylabel('Impedância $Z_L$ ($\Omega$)', 'Interpreter', 'latex');
legend('Parte Real $Re(Z_L)$', 'Parte Imaginária $Im(Z_L)$', sprintf('Ângulo Crítico $\\theta_c = %.2f^\\circ$', theta_c_deg), 'Interpreter', 'latex');
grid on;

%% ------------------- c) Relações de Potência -------------------

% Z_z1 = eta_1 * cos(theta_1) (Polarização TM)
Z_z1 = eta_1 * cos(theta_1_rad);

% Coeficiente de reflexão (rho_0)
rho_0 = (Z_L - Z_z1) ./ (Z_L + Z_z1);

% Relações de Potência
R_power_ratio = abs(rho_0).^2; % Potência Refletida |Nz1-|/|Nz1+|
T_power_ratio = 1 - R_power_ratio; % Potência Transmitida |Nz2+|/|Nz1+|

% Ângulo de Brewster (theta_p)
tan_theta_p = sqrt(epsilon_r2 / epsilon_r1);
theta_p_rad = atan(tan_theta_p);
theta_p_deg = rad2deg(theta_p_rad);

% Encontrar índices para valores específicos
[~, idx_0] = min(abs(theta_1_deg - 0.0));
[~, idx_p] = min(abs(theta_1_deg - theta_p_deg));
[~, idx_40] = min(abs(theta_1_deg - 40.0));

fprintf('\n--------- Questão (c) ---------\n');
fprintf('Ângulo de Brewster (theta_p): %.2f graus\n', theta_p_deg);
fprintf('Valores em theta_1 = 0 deg:\n');
fprintf('  Transmitida (T): %.3f\n', T_power_ratio(idx_0));
fprintf('  Refletida (R):   %.3f\n', R_power_ratio(idx_0));
fprintf('Valores em theta_1 = theta_p (%.2f deg):\n', theta_p_deg);
fprintf('  Transmitida (T): %.3f\n', T_power_ratio(idx_p));
fprintf('  Refletida (R):   %.3f (ideal: 0)\n', R_power_ratio(idx_p));
fprintf('Valores em theta_1 = 40 deg:\n');
fprintf('  Transmitida (T): %.3f\n', T_power_ratio(idx_40));
fprintf('  Refletida (R):   %.3f\n', R_power_ratio(idx_40));

% Gráfico (c)
figure;
plot(theta_1_deg, R_power_ratio, 'r', 'LineWidth', 1.5);
hold on;
plot(theta_1_deg, T_power_ratio, 'g', 'LineWidth', 1.5);
xline(theta_c_deg, 'k:', 'LineWidth', 1.5);
xline(theta_p_deg, 'b--', 'LineWidth', 1.5);
title('Gráfico (c): Relações de Potência vs $\theta_1$', 'Interpreter', 'latex');
xlabel('Ângulo de Incidência $\theta_1$ (graus)', 'Interpreter', 'latex');
ylabel('Relação de Potência');
legend('Refletida $\frac{|N_{z1}^{-}|}{|N_{z1}^{+}|}$', 'Transmitida $\frac{|N_{z2}^{+}|}{|N_{z1}^{+}|}$', ...
       sprintf('Ângulo Crítico $\\theta_c = %.2f^\\circ$', theta_c_deg), ...
       sprintf('Ângulo Brewster $\\theta_p = %.2f^\\circ$', theta_p_deg), ...
       'Interpreter', 'latex', 'Location', 'east');
grid on;
ylim([-0.05, 1.05]);

%% ------------------- d) Campo |H_y(z)| para 60 graus -------------------

theta_1_d_deg = 60.0; % 60 > theta_c, então é TIR
theta_1_d_rad = deg2rad(theta_1_d_deg);

% Parâmetros em 60 graus
Z_z1_d = eta_1 * cos(theta_1_d_rad);
argumento_sqrt_d = (epsilon_r1 / epsilon_r2) * sin(theta_1_d_rad)^2 - 1;
Z_L_d = -1i * eta_2 * sqrt(argumento_sqrt_d);
rho_0_d = (Z_L_d - Z_z1_d) / (Z_L_d + Z_z1_d);

% Vetor de posição z (em metros)
z_neg = linspace(-0.20, 0, 400); % -20 cm a 0
z_pos = linspace(0, 0.10, 200); % 0 a +10 cm

% z < 0 (Vidro): Onda Estacionária
beta_z1_d = k1 * cos(theta_1_d_rad);
Hy_neg = H_y1_plus * abs(exp(-1i * beta_z1_d * z_neg) - rho_0_d * exp(1i * beta_z1_d * z_neg));

% z > 0 (Ar): Onda Evanescente
alpha_z_d = k2 * sqrt(argumento_sqrt_d);
Hy_at_0 = H_y1_plus * abs(1 - rho_0_d);
Hy_pos = Hy_at_0 * exp(-alpha_z_d * z_pos);

% Valores Específicos
Hy_at_5cm = Hy_at_0 * exp(-alpha_z_d * 0.05);

fprintf('\n--------- Questão (d) ---------\n');
fprintf('Campo em theta_1 = 60 graus (TIR):\n');
fprintf('  |H_y(z=0)|   = %.3f mA_ef/m\n', Hy_at_0 * 1000);
fprintf('  |H_y(z=5cm)| = %.3f mA_ef/m\n', Hy_at_5cm * 1000);

% Gráfico (d)
figure;
plot(z_neg * 100, Hy_neg * 1000, 'b', 'LineWidth', 1.5); % z em cm, H em mA
hold on;
plot(z_pos * 100, Hy_pos * 1000, 'g', 'LineWidth', 1.5); % z em cm, H em mA
xline(0, 'k--', 'LineWidth', 1);
scatter(5.0, Hy_at_5cm * 1000, 'r', 'filled');
text(5.0, Hy_at_5cm * 1000 + 0.01, sprintf(' (5.0 cm, %.3f mA/m)', Hy_at_5cm * 1000), 'HorizontalAlignment', 'left');
title('Gráfico (d): Magnitude do Campo $|\dot{H}_y(z)|$ em $\theta_1 = 60^\circ$', 'Interpreter', 'latex');
xlabel('Posição z (cm)');
ylabel('$|\dot{H}_y|$ ($mA_{ef}/m$)', 'Interpreter', 'latex');
legend('Vidro (z<0) - Onda Estacionária', 'Ar (z>0) - Onda Evanescente', 'Interface z=0');
grid on;

fprintf('\nSimulação concluída.\n');