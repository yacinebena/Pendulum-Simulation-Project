function [temps, solution] = euler_explicite(f, intervalle_temps, y0, pas)
    % Initialiser les variables
    temps = intervalle_temps(1):pas:intervalle_temps(2); % Générer les valeurs de temps
    n_steps = length(temps);        % Nombre de pas de temps
    solution = zeros(length(y0), n_steps); % Initialiser la solution
    solution(:,1) = y0; % Condition initiale

    % Boucle sur chaque pas de temps
    for i = 1:(n_steps-1)
        t = temps(i);
        y = solution(:,i);

        % Calculer la prochaine valeur de y en utilisant la méthode d'Euler explicite
        solution(:,i+1) = y + pas * f(t, y);
    end
end

