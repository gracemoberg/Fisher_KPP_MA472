%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical simulation of the Fisher-KPP equation in 1D
% Closed domain using Crank-Nicholson and Adams-Bashforth.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear;

%% Set up parameters

% System parameters
par.D = 1;      % Diffusion coefficient
par.r = 1;      % Intrinsic growth rate
par.b = 1; %

% nx_steps = 200:100:1500;
% 
% for g = nx_steps

    % Numerical parameters
    nx = 1000;           % Number of grid points
    Lx = 500;           % Length of domain
    h = Lx / (nx-1);        % Spatial step size
    x = linspace(0, Lx, nx); % Domain
    
    tf = 800;           % Final time
    dt = 0.01;           % Time step
    iter = ceil(tf / dt); % Number of time steps
    plot_img = 100;      % Plot every plot_img steps
    
    % wave speed initializations
    threshold = 0.00001;  % Threshold for front detection
    point1 = floor(nx / 4);  % Spatial point 1 (quarter of the domain)
    point2 = floor(3* nx / 4);  % Spatial point 2 (three-quarters of domain)
    
    % Initialize variables to store crossing times
    time_point1 = [];
    time_point2 = [];

    C1 = floor(nx/3);
    C2 = floor(2*nx/3);
    
    % Diffusion operator (2nd derivative matrix, Neumann BCs)
    L2 = ComputeLinearOperator_1D_Fisher(nx, h);
    L2 = par.D * L2; % multiply by diffusion coefficient
    % L2(1:C, :) = par.D * L2(1:C, :);
    % L2(C+1:end, :) = (par.D/2) * L2(C+1:end, :);
    % L2(1:C1, : ) = par.D * L2(1:C1, :);
    % L2(C1+1:C2, :) = 2*par.D/3 *L2(C1+1:C2, :);
    % L2(C2+1:end, :) = par.D/3 * L2(C2+1:end, :);


    % Implicit matrix for Crank-Nicholson scheme
    DV = speye(nx) - (dt / 2) * L2;
    
    % Initial condition
    % u = exp(-0.5 * ((x - 10) / 10).^2)'; % Gaussian bump
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
        fu = fisher_KPP_harvesting_nonlin(u, par);
    
        % Crank-Nicholson + Adams-Bashforth for subsequent steps
        u = DV \ (u + (dt/2)*L2*u  + 3*dt*fu./2 - dt*fu_minus./2);
    
        fu_minus = fu;
    
        % Check if the wavefront has reached x1 or x2
       if isempty(time_point1) && u(point1) >= threshold
            time_point1 = (k - 1) * dt;  % Record rise time for point 1
        end
    
        if isempty(time_point2) && u(point2) >= threshold
            time_point2 = (k - 1) * dt;  % Record rise time for point 2
        end
    
        % Break early if both threshold crossings are recorded
        if ~isempty(time_point1) && ~isempty(time_point2)
            break;
        end
    
        % Plot results at intervals
        % if mod(k, plot_img) == 0
        %     figure(1);
        %     plot(x, u, 'b', 'LineWidth', 2);
        %     xlabel('x');
        %     ylabel('u(x,t)');
        %     % title(sprintf('Time = %.2f', k * dt));
        %     set(gca, 'LineWidth', 2, 'FontSize', 16);
        %     box on; xlim([x(1), x(end)]); ylim([0, 1.1]);
        %     drawnow; % pause(0.5);
        % end
    end
    
    if ~isempty(time_point1) && ~isempty(time_point2)
        distance = abs(x(point2) - x(point1));  % Spatial distance
        time_difference = abs(time_point2 - time_point1);  % Time difference for wavefront
        wave_speed = distance / time_difference;  % Wave speed calculation
        fprintf('Wave speed for nx = %d: %.3f units/time\n', wave_speed);
    else
        disp('Wave did not cross threshold at both points.');
    end

% end






