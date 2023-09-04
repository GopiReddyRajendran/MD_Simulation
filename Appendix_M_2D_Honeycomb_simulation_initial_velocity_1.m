% Author: Gopi Reddy Rajendran
% Description: This MATLAB code simulates a 2D molecular dynamics (MD)
% simulation with initial velocity as 1 of particles with a honeycomb pattern, where particles interact via a Lennard-Jones potential.

% Clear the workspace and command window
clc; % Clear the command window
clear; % Clear the workspace
close all; % Close all figure windows

% Simulation parameters
N = 320; % Number of particles
dt = 0.01; % Time step for the simulation
V = 1; % Initial velocity of particles
drcut = 4; % Cutoff distance for Lennard-Jones potential
T = 50; % Total simulation time
nSteps = T / dt; % Number of time steps in the simulation

% Calculate the cutoff distance squared for efficiency
drcut2 = drcut * drcut;

% Initialize arrays to store particle positions and velocities
xprev = zeros(1, N);
yprev = zeros(1, N);
xcurr = zeros(1, N);
ycurr = zeros(1, N);
xnew = zeros(1, N);
ynew = zeros(1, N);
vx0 = zeros(1, N);
vy0 = zeros(1, N);
vx = zeros(1, N);
vy = zeros(1, N);

% Generate initial positions for particles in a honeycomb pattern
% and random initial velocities within a small range
gx = [0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30];
gy = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19];

% Initialize initial positions and velocities for particles
for i = 1:N
    xcurr(i) = gx(i);
    ycurr(i) = gy(i);
    vx0(i) = 2 * (rand - 0.5) * V;
    vy0(i) = 2 * (rand - 0.5) * V;
    xprev(i) = xcurr(i) - vx0(i) * dt;
    yprev(i) = ycurr(i) - vy0(i) * dt;
end

% Create the "Honeycomb0" folder if it doesn't exist
if ~exist('Honeycomb1', 'dir')
    mkdir('Honeycomb1')
end

% Create VideoWriter object for creating a simulation video
outputVideo = VideoWriter('Honeycomb1.avi');
outputVideo.FrameRate = 30; % Frames per second
open(outputVideo);

% Create the initial scatter plot of particles and save it as an image
fig = figure('visible', 'off');
scatter(xcurr, ycurr, 50, 'b.');
title({'2D Molecular Dynamics Simulation (Honeycomb Pattern)'; 'step = 0'});
axis([0 32 0 20]);
drawnow;
frame = getframe(gcf);
writeVideo(outputVideo, frame);
saveas(fig, fullfile('Honeycomb1', sprintf('Initial.png')))
close(gcf);

% Main simulation loop
for n = 1:nSteps
    % Particle interaction calculations
    
    % Calculate distances between particles and their separations in x and y directions
    dx = repmat(xcurr', [1, N]) - repmat(xcurr, [N, 1]);
    dy = repmat(ycurr', [1, N]) - repmat(ycurr, [N, 1]);
    
    % Apply periodic boundary conditions to distances
    dx = dx - 32 * round(dx / 32);
    dy = dy - 20 * round(dy / 20);
    
    % Calculate the squared distances between particles
    dr2 = dx.^2 + dy.^2;
    
    % Find pairs of particles within the cutoff distance
    [row, col, DR2] = find(triu(dr2, 1));
    linearIndex = (col - 1) * N + row;
    dx = dx(linearIndex);
    dy = dy(linearIndex);
    
    invDR2 = zeros(size(DR2));
    cfIndex = (DR2 < drcut2);
    invDR2(cfIndex) = 1.0 ./ DR2(cfIndex);
    
    % Calculate Lennard-Jones forces between pairs
    f = 48 * invDR2.^4 .* (invDR2.^3 - 0.5);
    fx = f .* dx;
    fy = f .* dy;
    
    % Distribute the forces back to particles
    fx = full(sparse(row, col, fx, N, N));
    fy = full(sparse(row, col, fy, N, N));
    fx = sum(-fx' + fx, 2);
    fy = sum(-fy' + fy, 2);
    
    % Update particle positions and velocities
    for i = 1:N
        xnew(i) = 2 * xcurr(i) - xprev(i) + fx(i) * dt^2;
        ynew(i) = 2 * ycurr(i) - yprev(i) + fy(i) * dt^2;
        
        % Apply periodic boundary conditions
        if xnew(i) < 0
            xnew(i) = xnew(i) + 32;
        end
        if xnew(i) > 32
            xnew(i) = xnew(i) - 32;
        end
        if ynew(i) < 0
            ynew(i) = ynew(i) + 20;
        end
        if ynew(i) > 20
            ynew(i) = ynew(i) - 20;
        end
        
        vx(i) = (xnew(i) - xprev(i)) / (2 * dt);
        vy(i) = (ynew(i) - yprev(i)) / (2 * dt);
    end
    
    % Update previous positions
    xprev = xcurr;
    yprev = ycurr;
    
    % Update current positions
    xcurr = xnew;
    ycurr = ynew;
    
    % Plot the current configuration and save it as a PNG file
    fig = figure('visible', 'off');
    scatter(xcurr, ycurr, 50, 'b.');
    title({'2D Molecular Dynamics Simulation (Honeycomb Pattern)'; ['step = ', num2str(n)]});
    axis([0 32 0 20]);
    drawnow;
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);
    
    % Save an image every 10 time steps
    if mod(n, 10) == 0
        saveas(fig, fullfile('Honeycomb1', sprintf('step%d.png', n)))
    end
    
    close(gcf);
end

% Finish video
close(outputVideo);
