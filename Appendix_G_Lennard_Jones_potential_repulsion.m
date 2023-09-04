% Title: Random Motion of Multiple Particles with Lennard-Jones Repulsion
% Author: Gopi Reddy Rajendran
% Description:
% This MATLAB code simulates the motion of multiple particles in a 1D space with Lennard-Jones repulsion forces.
% The particles are initially distributed within a 1D box of size 'L'.
% Each particle undergoes a random walk with Lennard-Jones repulsion interactions.
% The simulation is carried out over a specified time ('T') with a given time step ('dt').
% The particles' positions are updated at each time step, and the simulation results are visualized as both a video and individual PNG frames.

% Clear the workspace and command window
clc; % Clear command window
clear; % Clear workspace
close all; % Close all figure windows

% Simulation parameters
N = 2; % Number of particles
d = 3; % Distance of each particle (not used in this simulation)
dt = 0.01; % Time step for simulation
dr = 0.1; % Initial displacement of particles from grid points
V = 0; % Initial velocity of particles (not used in this simulation)
drcut = 5; % Cutoff distance for Lennard-Jones potential
T = 4.2; % Total simulation time
nSteps = T / dt; % Number of time steps in simulation
plotDelay = 0 / 1; % Time delay between plot updates (not used in this simulation)

% Determine the size of the simulation box
L = 10;

% Calculate the cutoff distance squared for efficiency
drcut2 = drcut * drcut;

% Generate initial positions and velocities for particles
% within a small range of the grid points
for i = 1:N
    xcurr = [3 4]; % Initial positions (customize as needed)
    vx0(i) = 2 * (rand - 0.5) * V; % Initial velocities (not used in this simulation)
    xprev(i) = xcurr(i) - vx0(i) * dt; % Previous positions
end

% Create a folder for saving simulation results
if ~exist('LJ_Repulsion', 'dir')
    mkdir('LJ_Repulsion')
end

% Create VideoWriter object for saving animation
outputVideo = VideoWriter('LJ_Repulsion.avi');
outputVideo.FrameRate = 30; % Frames per second
open(outputVideo);

% Create an initial scatter plot for the first frame
fig = figure('visible', 'off');
scatter(xcurr, 0, 500, 'b.'); % Initial positions as blue dots
title({'Lennard-Jones Repulsion Simulation'; 'step = 0'});
axis([0 L -L L]);
drawnow;
frame = getframe(gcf);
writeVideo(outputVideo, frame);
saveas(fig, fullfile('LJ_Repulsion', sprintf('Initial.png')))
close(gcf);

% Main simulation loop
for n = 1:nSteps
    % Calculate distances between particles and their separations in x direction
    dx = repmat(xcurr', [1, N]) - repmat(xcurr, [N, 1]);
    
    % Apply periodic boundary conditions to distances
    dx = dx - L * round(dx / L);
    
    % Calculate the squared distances between particles
    dr2 = dx.^2;
    [row, col, DR2] = find(triu(dr2, 1));
    linearIndex = (col - 1) * N + row;
    dx = dx(linearIndex);
    
    % Calculate Lennard-Jones forces
    invDR2 = zeros(size(DR2));
    cfIndex = (DR2 < drcut2);
    invDR2(cfIndex) = 1 ./ DR2(cfIndex);
    f = 48 * invDR2.^4 .* (invDR2.^3 - 0.5);
    fx = f .* dx;
    fx = full(sparse(row, col, fx, N, N));
    fx = sum(-fx' + fx, 2);
    
    % Update particle positions and velocities
    for i = 1:N
        xnew(i) = 2 * xcurr(i) - xprev(i) + fx(i) * dt^2;
        
        % Apply periodic boundary conditions
        if xnew(i) < 0
            xnew(i) = xnew(i) + L;
        end
        if xnew(i) > L
            xnew(i) = xnew(i) - L;
        end
        
        vx(i) = (xnew(i) - xprev(i)) / (2 * dt);
    end
    
    % Update previous positions
    for i = 1:N
        xprev(i) = xcurr(i);
        xcurr(i) = xnew(i);
    end
    
    % Plot the current configuration and save it as a PNG file
    fig = figure('visible', 'off');
    scatter(xcurr, 0, 500, 'b.'); % Current positions as blue dots
    title({'Lennard-Jones Repulsion Simulation'; ['step = ', num2str(n)]});
    axis([0 L -L L]);
    drawnow;
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);
    
    % Save the plot periodically (adjust as needed)
    if mod(n, 10) == 0
        saveas(fig, fullfile('LJ_Repulsion', sprintf('step%d.png', n)))
    end
    
    close(gcf);
end

% Finish video
close(outputVideo);