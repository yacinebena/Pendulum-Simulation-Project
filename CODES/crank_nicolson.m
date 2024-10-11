function [temps, solution] = crank_nicolson(f, intervalle_temps, y0, pas, tol)
    % Initialiser les variables
    temps = intervalle_temps(1):pas:intervalle_temps(2); % Générer les valeurs de temps
    n_steps = length(temps);        % Nombre de pas de temps
    solution = zeros(length(y0), n_steps); % Initialiser la solution
    solution(:,1) = y0; % Condition initiale

    % Boucle sur chaque pas de temps
    for i = 1:(n_steps-1)
        t_n = temps(i);
        y_n = solution(:,i);

        % Estimation initiale pour y_{n+1} en utilisant Euler explicite
        y_next = y_n + pas * f(t_n, y_n);

        % Résoudre l'équation implicite par itération (point fixe)
        diff = inf;
        while diff > tol
            y_next_old = y_next;
            y_next = y_n + (pas / 2) * (f(t_n, y_n) + f(t_n + pas, y_next_old));
            diff = norm(y_next - y_next_old);
        endwhile

        % Mettre à jour la solution
        solution(:,i+1) = y_next;
    end
end

