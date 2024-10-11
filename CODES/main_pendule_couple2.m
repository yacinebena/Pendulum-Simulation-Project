% Définition des paramètres
g = 9.81; % Accélération due à la gravité
k = 1.0;  % Constante du système
l = 1.0;  % Constante de frottement
m = 1.0;  % Constante de couplage
dt = 0.01; % Pas de temps
T = 5;   % Temps total de simulation
tol = 1e-6; % Tolérance pour les méthodes implicites

% Conditions initiales
theta1_0 = pi/4; % Angle initial de theta1
theta2_0 = pi/6; % Angle initial de theta2
omega1_0 = 0.0;  % Vitesse angulaire initiale de theta1
omega2_0 = 0.0;  % Vitesse angulaire initiale de theta2

% Conditions initiales sous forme de vecteur
y0 = [theta1_0; omega1_0; theta2_0; omega2_0];

% Fonction définissant le système d'équations différentielles
f = @(t, y) [
    y(2); % d(theta1)/dt = omega1
    -(g * k * y(1) + l * y(2) + m * l^2 * (y(1) - y(3))); % d(omega1)/dt = alpha1
    y(4); % d(theta2)/dt = omega2
    -(g * k * y(3) + l * y(4) + m * l^2 * (y(3) - y(1)))  % d(omega2)/dt = alpha2
];

% Intervalle de temps
intervalle_temps = [0 T];

% Résolution avec Euler explicite
[temps_euler_exp, sol_euler_exp] = euler_explicite(f, intervalle_temps, y0, dt);

% Résolution avec RK4
[temps_rk4, sol_rk4] = RK4(f, intervalle_temps, y0, dt);

% Résolution avec Euler implicite
[temps_euler_imp, sol_euler_imp] = euler_implicite(f, intervalle_temps, y0, dt, tol);

% Résolution avec Crank-Nicolson
[temps_crank, sol_crank] = crank_nicolson(f, intervalle_temps, y0, dt, tol);

% Affichage des résultats
figure;
hold on;
plot(temps_euler_exp, sol_euler_exp(1, :), 'r', 'DisplayName', 'Euler Explicite - \theta_1');
plot(temps_euler_exp, sol_euler_exp(3, :), 'b', 'DisplayName', 'Euler Explicite - \theta_2');
plot(temps_rk4, sol_rk4(1, :), 'r--', 'DisplayName', 'RK4 - \theta_1');
plot(temps_rk4, sol_rk4(3, :), 'b--', 'DisplayName', 'RK4 - \theta_2');
plot(temps_euler_imp, sol_euler_imp(1, :), 'g', 'DisplayName', 'Euler Implicite - \theta_1');
plot(temps_euler_imp, sol_euler_imp(3, :), 'c', 'DisplayName', 'Euler Implicite - \theta_2');
plot(temps_crank, sol_crank(1, :), 'g--', 'DisplayName', 'Crank-Nicolson - \theta_1');
plot(temps_crank, sol_crank(3, :), 'c--', 'DisplayName', 'Crank-Nicolson - \theta_2');
xlabel('Temps (s)');
ylabel('Angle (rad)');
legend show;
title('Comparaison des méthodes numériques pour la résolution des pendules couplés');
grid on;
hold off;

