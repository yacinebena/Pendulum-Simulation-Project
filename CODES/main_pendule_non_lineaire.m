% Définir les constantes
g = 9.81;  % Accélération due à la gravité (m/s^2)
l = 1.0;   % Longueur du pendule (m)

% Définir la fonction f qui représente le système d'EDO non linéaire
function dydt = systeme_pendule_non_lineaire(t, y)
    g = 9.81;  % Accélération due à la gravité (m/s^2)
    l = 1.0;   % Longueur du pendule (m)
    theta = y(1);
    omega = y(2);
    dtheta_dt = omega;
    domega_dt = -(g/l) * sin(theta);
    dydt = [dtheta_dt; domega_dt];
endfunction

% Conditions initiales
theta0 = pi / 4;  % Angle initial (45 degrés)
omega0 = 0.0;     % Vitesse angulaire initiale (rad/s)
y0 = [theta0; omega0];

% Domaine de la solution
intervalle_temps = [0 100];  % Résoudre de t = 0 à t = 10 secondes
pas = 0.01;                 % Taille du pas
tol = 1e-6;                 % Tolérance pour les méthodes implicites

% Résoudre avec différentes méthodes
[temps_euler, solution_euler] = euler_explicite(@systeme_pendule_non_lineaire, intervalle_temps, y0, pas);
[temps_euler_imp, solution_euler_imp] = euler_implicite(@systeme_pendule_non_lineaire, intervalle_temps, y0, pas, tol);
[temps_rk4, solution_rk4] = RK4(@systeme_pendule_non_lineaire, intervalle_temps, y0, pas);
[temps_cn, solution_cn] = crank_nicolson(@systeme_pendule_non_lineaire, intervalle_temps, y0, pas, tol);

% Tracer les résultats
figure;
subplot(2, 1, 1);
plot(temps_euler, solution_euler(1, :), 'b', 'DisplayName', 'Euler explicite');
hold on;
plot(temps_euler_imp, solution_euler_imp(1, :), 'r', 'DisplayName', 'Euler implicite');
plot(temps_rk4, solution_rk4(1, :), 'g', 'DisplayName', 'RK4');
plot(temps_cn, solution_cn(1, :), 'm', 'DisplayName', 'Crank-Nicolson');
xlabel('Temps (s)');
ylabel('Angle θ(t) (rad)');
title('Évolution de l’angle du pendule');
legend('show');
grid on;

subplot(2, 1, 2);
plot(temps_euler, solution_euler(2, :), 'b', 'DisplayName', 'Euler explicite');
hold on;
plot(temps_euler_imp, solution_euler_imp(2, :), 'r', 'DisplayName', 'Euler implicite');
plot(temps_rk4, solution_rk4(2, :), 'g', 'DisplayName', 'RK4');
plot(temps_cn, solution_cn(2, :), 'm', 'DisplayName', 'Crank-Nicolson');
xlabel('Temps (s)');
ylabel('Vitesse angulaire ω(t) (rad/s)');
title('Évolution de la vitesse angulaire du pendule');
legend('show');
grid on;



% Affichage des figures
shg;

