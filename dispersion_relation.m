%% Plot the dispersion relation

a = linspace(0.01, 5, 500); % Avoid a=0 to prevent division by zero

% Compute 'c'
c = a + 1 ./ a;

% Plot the function
figure;
plot(a, c, 'LineWidth', 1.5);
xlabel('a', FontSize=12);
ylabel('c', FontSize=12);
title('Plot of c = a + 1/a', FontSize=14);
grid on;
xlim([0 5]); ylim([0,10]); % Optional, focus on range of 'a'
