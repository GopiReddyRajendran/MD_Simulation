% Author: Gopi Reddy Rajendran
% Description: 
%   This MATLAB code simulates a 1D molecular dynamics system with Lennard-Jones
%   interactions. It initializes particles, calculates forces, and records the
%   simulation as both a video and individual PNG frames
%
%   The code generates an initial configuration of particles within a specified range,
%   applies periodic boundary conditions, calculates inter-particle forces using the
%   Lennard-Jones potential, and updates particle positions and velocities. The simulation
%   progresses through multiple time steps, and at each step, it plots the current
%   configuration, saving snapshots of the system as PNG images.
%
%   The simulation allows you to explore the behavior of particles in a 1D system with
%   Lennard-Jones forces and observe how they evolve over time.

% Clear the workspace and command window
clc; % Clear command window
clear; % Clear workspace
close all; % Close all figure windows

% Simulation parameters
N = 5; % Number of particles
d = 3; % Distance of each particle (not used in this simulation)
dt = 0.01; % Time step for simulation
dr = 0.1; % Initial displacement of particles from grid points
V = 1; % Initial velocity of particles
drcut = 3; % Cutoff distance for Lennard-Jones potential
T = 5; % Total simulation time
nSteps = T/dt; % Number of time steps in the simulation
plotDelay = 0.1; % Time delay between plot updates

% Determine the size of the simulation box
L = 10;

% Calculate the cutoff distance squared for efficiency
drcut2 = drcut*drcut;

% Generate initial positions and velocities for particles
% within a small range of the grid points
gx = [0.5 2.5 4.5 6.5 8.5];
for i = 1:N
    xcurr(i) = gx(i) + 2*(rand-0.5)*dr; % Initial positions with random variations
    vx0(i) = 2*(rand-0.5)*V; % Initial velocities with random variations
    xprev(i) = xcurr(i) - vx0(i)*dt; % Previous positions
end

% Create the "1D_Result" folder if it doesn't exist
if ~exist('1D_Result', 'dir')
    mkdir('1D_Result')
end

% Create VideoWriter object for saving animation
outputVideo = VideoWriter('1D_Result.avi');
outputVideo.FrameRate = 30; % Frames per second
open(outputVideo);

% Create an initial scatter plot for the first frame
fig = figure('visible', 'off');
scatter(xcurr,0,500,'b.'); % Initial positions as blue dots
title({'1-Dimensional Molecular Dynamics Simulation'; 'step = 0'});
axis([0 L -L L]);
drawnow;
frame = getframe(gcf);
writeVideo(outputVideo, frame);
saveas(fig, fullfile('1D_Result', sprintf('Initial.png')))
close(gcf);

% Main simulation loop
for n = 1:nSteps
    % Calculate distances between particles and their separations in x direction
    dx = repmat(xcurr',[1,N])-repmat(xcurr,[N,1]);
    
    % Apply periodic boundary conditions to distances
    dx = dx - L*round(dx/L);
    
    % Calculate the squared distances between particles
    dr2 = dx.^2;
    [row,col,DR2] = find(triu(dr2,1));
    linearIndex = (col-1)*N+row;
    dx = dx(linearIndex);
    invDR2 = zeros(size(DR2));
    cfIndex = (DR2<drcut2);
    invDR2(cfIndex) = 1./DR2(cfIndex);
    f = 48*invDR2.^4.*(invDR2.^3 - 0.5);
    fx = f.*dx;
    fx = full(sparse(row,col,fx,N,N));
    fx = sum(-fx'+fx,2);
    
    % Update particle positions and velocities
    for i = 1:N
        xnew(i) = 2*xcurr(i) - xprev(i) + fx(i)*dt^2;
        
        % Apply periodic boundary conditions
        if xnew(i) < 0
            xnew(i) = xnew(i) + L;
        end
        if xnew(i) > L
            xnew(i) = xnew(i) - L;
        end
        vx(i) = (xnew(i) - xprev(i))/(2*dt);
    end
    for i = 1:N
        xprev(i) = xcurr(i);
        xcurr(i) = xnew(i);
    end
    
    % Plot the current configuration and save it as a PNG file
    fig = figure('visible', 'off');
    scatter(xcurr,0,500,'b.'); % Current positions as blue dots
    title({'1-Dimensional Molecular Dynamics Simulation';['step = ',num2str(n)]});
    axis([0 L -L L]);
    drawnow;
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);
    
    % Save the plot periodically (adjust as needed)
    if mod(n,10) == 0
        saveas(fig, fullfile('1D_Result', sprintf('step%d.png', n)))
    end
    close(gcf);    
end

% Finish video
close(outputVideo);
