% Plot line of best fit for harvesting and wave speed
b = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]; % Harvesting parameter (years^-1)
c = [2.013, 1.906, 1.793, 1.674, 1.547, 1.410, 1.258, 1.087, 0.884, 0.621]; % Wave speed (km/years)

% Plot the data
figure;
plot(b, c, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'Computed Wave Speeds'); s
hold on;

% Fit a line to the data
coeffs = polyfit(b, c, 1); % Linear fit (1st degree polynomial)
c_fit = polyval(coeffs, b); % Evaluate the fit at the data points

% Plot the line of best fit
plot(b, c_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Line of Best Fit'); % Red line for best fit

% Add labels, title, and legend
xlabel('Harvesting');
ylabel('Wave speed');
title('Wave Speed vs. Harvesting');
legend('Location', 'northeast');
grid on;
hold off;
