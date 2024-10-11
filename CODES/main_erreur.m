% Définir les constantes
g = 9.81;  % Accélération due à la gravité (m/s^2)
l = 1.0;   % Longueur du pendule (m)

% Définir la fonction f qui représente le système d'EDO
function dydt = systeme_pendule(t, y)
    g = 9.81;  % Accélération due à la gravité (m/s^2)
    l = 1.0;   % Longueur du pendule (m)
    theta = y(1);
    X = y(2);
    dtheta_dt = X;
    dX_dt = -(g/l) * sin(theta);
    dydt = [dtheta_dt; dX_dt];
endfunction

% Solution exacte pour le pendule simple (petites oscillations)
function theta_exact = solution_exact(t, theta0, g, l)
    omega = sqrt(g / l);  % Fréquence angulaire
    theta_exact = theta0 * cos(omega * t);
endfunction

% Conditions initiales
theta0 = 1e-13;  % Angle initial (45 degrés)
X0 = 0.0;         % Vitesse angulaire initiale (rad/s)
y0 = [theta0; X0];

% Domaine de la solution
intervalle_temps = [0 5];  % Résoudre de t = 0 à t = 10 secondes
pas = 0.001;                 % Taille du pas
tol = 1e-14;                 % Tolérance pour les méthodes implicites

% Résoudre avec différentes méthodes
[temps_euler, solution_euler] = euler_explicite(@systeme_pendule, intervalle_temps, y0, pas);
[temps_euler_imp, solution_euler_imp] = euler_implicite(@systeme_pendule, intervalle_temps, y0, pas, tol);
[temps_rk4, solution_rk4] = RK4(@systeme_pendule, intervalle_temps, y0, pas);
[temps_cn, solution_cn] = crank_nicolson(@systeme_pendule, intervalle_temps, y0, pas, tol);

% Calculer la solution exacte
theta_exact = solution_exact(temps_rk4, theta0, g, l);

% Extraire les résultats
theta_euler = solution_euler(1, :);
X_euler = solution_euler(2, :);

theta_euler_imp = solution_euler_imp(1, :);
X_euler_imp = solution_euler_imp(2, :);

theta_rk4 = solution_rk4(1, :);
X_rk4 = solution_rk4(2, :);

theta_cn = solution_cn(1, :);
X_cn = solution_cn(2, :);

% Calculer les erreurs
erreur_euler = norm(theta_euler - solution_exact(temps_euler, theta0, g, l)) / norm(solution_exact(temps_euler, theta0, g, l));
erreur_euler_imp = norm(theta_euler_imp - solution_exact(temps_euler_imp, theta0, g, l)) / norm(solution_exact(temps_euler_imp, theta0, g, l));
erreur_rk4 = norm(theta_rk4 - solution_exact(temps_rk4, theta0, g, l)) / norm(solution_exact(temps_rk4, theta0, g, l));
erreur_cn = norm(theta_cn - solution_exact(temps_cn, theta0, g, l)) / norm(solution_exact(temps_cn, theta0, g, l));

% Afficher les erreurs
printf('Erreur relative pour Euler explicite : %e\n', erreur_euler);
printf('Erreur relative pour Euler implicite : %e\n', erreur_euler_imp);
printf('Erreur relative pour RK4 : %e\n', erreur_rk4);
printf('Erreur relative pour Crank-Nicolson : %e\n', erreur_cn);

% Tracer les résultats
figure;

% Tracé de θ(t) pour chaque méthode
subplot(2, 1, 1);
plot(temps_euler, theta_euler, 'b', 'DisplayName', 'Euler explicite');
hold on;
plot(temps_euler_imp, theta_euler_imp, 'r', 'DisplayName', 'Euler implicite');
plot(temps_rk4, theta_rk4, 'g', 'DisplayName', 'RK4');
plot(temps_cn, theta_cn, 'm', 'DisplayName', 'Crank-Nicolson');
plot(temps_rk4, solution_exact(temps_rk4, theta0, g, l), 'k--', 'DisplayName', 'Solution exacte');
xlabel('Temps (s)');
ylabel('Angle θ(t) (rad)');
title('Évolution de θ(t) au cours du temps');
legend('show');
grid on;

% Tracé de X(t) pour chaque méthode
subplot(2, 1, 2);
plot(temps_euler, X_euler, 'b', 'DisplayName', 'Euler explicite');
hold on;
plot(temps_euler_imp, X_euler_imp, 'r', 'DisplayName', 'Euler implicite');
plot(temps_rk4, X_rk4, 'g', 'DisplayName', 'RK4');
plot(temps_cn, X_cn, 'm', 'DisplayName', 'Crank-Nicolson');
xlabel('Temps (s)');
ylabel('Vitesse angulaire X(t) (rad/s)');
title('Évolution de X(t) au cours du temps');
legend('show');
grid on;

% Affichage des figures
shg;

