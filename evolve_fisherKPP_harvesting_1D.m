%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical simulation of the Fisher-KPP equation in 1D
% Closed domain using Crank-Nicholson and Adams-Bashforth.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear;

%% Set up parameters

% System parameters
par.D = 1;      % Diffusion coefficient
par.r = 1;      % Intrinsic growth rate
par.b = 0.1; %

% Numerical parameters
nx = 1000;           % Number of grid points
Lx = 385;           % Length of domain
h = Lx / (nx-1);        % Spatial step size
x = linspace(0, Lx, nx); % Domain

tf = 100;           % Final time
dt = 0.01;           % Time step
iter = ceil(tf / dt); % Number of time steps
plot_img = 100;      % Plot every plot_img steps

C = nx/2;

% Diffusion operator (2nd derivative matrix, Neumann BCs)
L2 = ComputeLinearOperator_1D_Fisher(nx, h);
L2 = par.D * L2; % multiply by diffusion coefficient
% L2(1:C, :) = par.D * L2(1:C, :);
% L2(C+1:end, :) = (par.D/12) * L2(C+1:end, :);


% Implicit matrix for Crank-Nicholson scheme
DV = speye(nx) - (dt / 2) * L2;

% Initial condition
% u = exp(-0.05 * ((x - 1) / 10).^2)'; % Gaussian bump
u = exp(-x)';

%% Evolve in time

for k =1 
    % Crank-Nicholson for single time step to get history
    fu = fisher_KPP_harvesting_nonlin(u, par);
    u = DV \ (u + (dt/2)*L2*u + dt*fu);
end

fu_minus = fu;

for k = 2:iter
    % nonlinear reaction term
    fu = fisher_KPP_nonlin(u, par);

    % Crank-Nicholson + Adams-Bashforth for subsequent steps
    u = DV \ (u + (dt/2)*L2*u  + 3*dt*fu./2 - dt*fu_minus./2);

    fu_minus = fu;

    % Plot results at intervals
    if mod(k, plot_img) == 0
        figure(1);
        plot(x, u, 'b', 'LineWidth', 2);
        xlabel('x');
        ylabel('u(x,t)');
        title(sprintf('Time = %.2f', k * dt));
        set(gca, 'LineWidth', 2, 'FontSize', 16);
        box on; xlim([x(1), x(end)]); ylim([0, 1.1]);
        drawnow; % pause(0.2);
    end
end




