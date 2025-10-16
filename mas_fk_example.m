%% Ejemplo FK: 
% -
% Modelo mecánico M-R-A. Suponemos un sistema "real" con 
% dinámicas no modeladas y variabilidad en parámetros de modelo. 
% Simplificamos el sistema un caso ideal.
% Estará instrumentado con acelerómetro y sensor de profundidad z = (velocidad, distancia)'  
clc; clear; close all; 
addpath("sims")

%% Parámetros Físicos 
m = 1; % [kg ]
c = 1; % [N.s/m]
k = 1; % [N/m]

%% Variables de estado
% Definimos las variables de estado 

A = [0 1; -k/m -c/m]; 
B = [0 , 1/m]'; % Entrada en Fuerza. Conversión aceleración
H = [1, 0; -k/m, -c/m]; % Misma C en variables de estado. Matriz de Observación
D_meas = [0; 1/m]; % Causalidad del sistema 

G = ss(A, B, H, D_meas); % Sistema Continuo
G.InputName = "Force"; 
G.InputUnit = "N"; 

G.OutputName = {'position'; 'acceleration'}; 

[y,t] = step(G); % Real Esperado en implementación digital. 

pos_mod = y(:, 1); 
vel_mod = y(:, 2); 

%% Simulación en Simscape 
simout = sim("mass_fk_example_sim");  % output format: meas_acc, meas_pos, real_acc, real_pos
simt = simout.simout.Time ; 

% Measured
meas_acc = simout.simout.Data(:, 1); 
meas_pos = simout.simout.Data(:, 2); 

real_acc = simout.simout.Data(:, 3); 
real_pos = simout.simout.Data(:, 4); 

u        = simout.simout.Data(:, 5); 

Ts = simt(2) - simt(1); % Tiempo en discretizacion 

%% Discretizacion de modelo 
Ad = [1 Ts; -k/m*Ts (1 - c/m*Ts)]; 
Bd = [0 Ts/m]'; 


%% Comparacion entre real + medidas 
subplot(2, 1, 1); 
plot(simt, meas_acc, "-", simt, real_acc, "r-"); 
grid on; 
title("Comparativa de aceleraciones"); 
xlabel("time [s]"); ylabel("Acc [m/s^2]"); 
legend("meas", "real");

subplot(2, 1, 2); 
plot(simt, meas_pos, simt, real_pos); 
grid on; 
title("Comparativa de posiciones"); 
xlabel("time [s]"); ylabel("Pos [m]"); 
legend("meas", "real");


%% Viz 
figure; 
plot(simt, real_pos, t, pos_mod)
grid on; 
title("Posisiton for MRA system"); 
xlabel("Time [s]"); 
ylabel("pos [m]"); 

legend("model", "real sin ruido")

%% Filtro de Kalman 
N = length(simt); 

% inicialización de vectores - ejemplo offline 
x_hat_priori = zeros(2, N); 
x_hat_posteriori = zeros(2, N); 
x_hat_posteriori(:, 1) = [0; 2]; 

P_priori = zeros(2,2, N); 
P_posterior = zeros(2,2, N); 

x_hat_posteriori(:,1) = [meas_pos(2); 0];
P_posterior(:,:,1) = diag([0.1, 0.1]); 

K = zeros(2, 2, N); 

% Matrices de Kalman 
Q = diag([1e-4, 1e-3]);   % proceso: permite cambios moderados en pos/vel
R = diag([0.7, 0.3]);   % mediciones: confío más en pos, menos en acc



%% Implementacion digital de filtro de Kalman 
for k = 1:(N-1)
    % Predicción 
    x_hat_priori(:, k + 1) = Ad * x_hat_posteriori(:, k) + Bd * u(k+1);
    P_priori(:,:, k + 1) = Ad * P_posterior(:, :, k) * Ad' + Q; % Covarianza - Predecida 

    % Observacion correctiva 
    K(:, :, k + 1) = P_priori(:,:, k + 1) * H' / (H * P_priori(:,:, k + 1) * H' + R);
    
    % Medicion
    z = [meas_pos(k); meas_acc(k)]; 
    
    x_hat_posteriori(:, k + 1) = x_hat_priori(:, k + 1) + K(:, : , k + 1) * (z - (H * x_hat_priori(:, k + 1) + D_meas * u(k))); % Estimación corregida
    P_posterior(:, :, k + 1) = (eye(2) - K(:, : ,k + 1) * H) * P_priori(:,:, k + 1); % Covarianza - Corregida 
end 

%%  Vizualización de resultados 
figure;

pos_priori = x_hat_priori(1, :)'; 
pos_posteriori= x_hat_posteriori(1, :)'; 

% subplot(2,1,1);
plot(real_pos, 'k'); hold on;
plot(pos_priori, '--b');
plot(pos_posteriori, 'r');

title("Comparativa de Posición");
hold on; 
plot(meas_pos)
hold off; 
legend("Real", "A priori", "A posteriori", "Measured");
grid on;


%% Dispersión ? 
theta = linspace(0, 2*pi, 50);   % Para generar elipse

T = [];
Y = [];
Z = [];

for k = 1:N
    P = P_priori(:,:,k);      % <-- Cambia a P_posteriori(:,:,k) si quieres el filtrado
    [V, D] = eig(P);
    L = V * sqrt(D);
    ellipse = L * [cos(theta); sin(theta)];

    T = [T; simt(k) * ones(1, length(theta))];
    Y = [Y; ellipse(1,:)];
    Z = [Z; ellipse(2,:)];
end

figure;
surf(T, Y, Z, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
colormap turbo;
title('Evolución de la covarianza (túnel)');
xlabel('Tiempo [s]');
ylabel('Δ Posición');
zlabel('Δ Velocidad');
grid on;
view(3);
