%% Phase plane

%% system we want to plot: 
%% u' = w
%% w' = -cw - u * (1-u)

f = @(t,Y) [Y(2); -3*Y(2) - Y(1)*(1-Y(1))];

y1 = linspace(-0.5, 1.25, 20);
y2 = linspace(-0.75, 0.5, 20);

[x,y] = meshgrid(y1, y2);

u = zeros(size(x));
v = zeros(size(x));

t = 0;
for i = 1:numel(x)
    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

quiver(x,y,u,v,'r');
xlabel('u');
ylabel('w');
axis tight equal;

% Define desired initial conditions
initial_conditions = [
    0.99, 0;       % (u0, w0) = (1, 0)
    0.5, -0.5;     % (u0, w0) = (0.5, -2)
    0.25, -0.5;
    0.7, 0.1;
];

% Prepare figure
hold on
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);

% Loop over desired initial conditions
for ic = 1:size(initial_conditions, 1)
    u0 = initial_conditions(ic, 1);
    w0 = initial_conditions(ic, 2);
    [ts, ys] = ode15s(f, [0, 20], [u0; w0], options);
    
    % Plot the trajectory
    plot(ys(:, 1), ys(:, 2), 'LineWidth', 1.5) % (u, w) phase plot
    plot(ys(1, 1), ys(1, 2), 'bo') % Initial point
    plot(ys(end, 1), ys(end, 2), 'ks') % Final point
end

% Finalize plot
xlabel('u', FontSize=12);
ylabel('w', FontSize=12);
axis tight equal;
title('Phase Plane for Fisher-KPP with c=3', FontSize=16);
hold off;
