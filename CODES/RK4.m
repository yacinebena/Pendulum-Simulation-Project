function [temps, solution] = RK4(f, intervalle_temps, y0, pas)
    % Initialiser les variables
    temps = intervalle_temps(1):pas:intervalle_temps(2); % Générer les valeurs de temps
    n_steps = length(temps);        % Nombre de pas de temps
    solution = zeros(length(y0), n_steps); % Initialiser la solution
    solution(:,1) = y0; % Condition initiale

    % Boucle sur chaque pas de temps
    for i = 1:(n_steps-1)
        t = temps(i);
        y = solution(:,i);

        % Calculer k1, k2, k3, k4
        k1 = f(t, y);
        k2 = f(t + pas/2, y + pas/2 * k1);
        k3 = f(t + pas/2, y + pas/2 * k2);
        k4 = f(t + pas, y + pas * k3);

        % Calculer la prochaine valeur de y
        solution(:,i+1) = y + pas/6 * (k1 + 2*k2 + 2*k3 + k4);
    end
end

