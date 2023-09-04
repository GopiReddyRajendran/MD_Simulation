% Author: Gopi Reddy Rajendran
% Description: This code simulates the 2D molecular dynamics of particles in
% a square simulation box. The particles interact through a Lennard-Jones
% potential, and their positions and velocities are updated over time using
% numerical integration. The code generates a video of the particle motion,
% as well as histograms of the speed and velocity distributions.

% Clear the workspace and command window
clc; % Clear command window
clear; % Clear workspace
close all; % Close all figure windows

% Simulation parameters
N = 100; % Number of particles
d = 3; % Distance of each particle
dt = 0.01; % Time step for simulation
dr = 0.1; % Initial displacement of particles from grid points
V = 1; % Initial velocity of particles
drcut = 4; % Cutoff distance for Lennard-Jones potential
T = 50; % Total simulation time
nSteps = T / dt; % Number of time steps in simulation

% Round N to the nearest square number to create a square grid of particles
N = round(sqrt(N))^2;

% Determine the size of the simulation box
L = round(d * sqrt(N));

% Calculate the cutoff distance squared for efficiency
drcut2 = drcut * drcut;

% Generate initial positions and velocities for particles
% within a small range of the grid points
[gx, gy] = meshgrid(d / 2:d:L - (d / 2), d / 2:d:L - (d / 2));
gx = gx(1:N);
gy = gy(1:N);
for i = 1:N
    xcurr(i) = gx(i) + 2 * (rand - 0.5) * dr;
    ycurr(i) = gy(i) + 2 * (rand - 0.5) * dr;
    vx0(i) = 2 * (rand - 0.5) * V;
    vy0(i) = 2 * (rand - 0.5) * V;
    xprev(i) = xcurr(i) - vx0(i) * dt;
    yprev(i) = ycurr(i) - vy0(i) * dt;
end

% Create the "plots" folder if it doesn't exist
if ~exist('2D_Initial_Velocity_1', 'dir')
    mkdir('2D_Initial_Velocity_1')
end
% Create VideoWriter object
outputVideo = VideoWriter('2D_Initial_Velocity_1.avi');
outputVideo.FrameRate = 30; % Frames per second
open(outputVideo);

fig = figure('visible', 'off');
scatter(xcurr, ycurr, 500, 'b.');
title({'2D Molecular Dynamics Simulation (Initial Velocity = 1)'; 'step = 0'});
axis([0 L 0 L]);
drawnow;
frame = getframe(gcf);
writeVideo(outputVideo, frame);
saveas(fig, fullfile('2D_Initial_Velocity_1', sprintf('Initial.png')))
close(gcf);

% Main simulation loop
for n = 1:nSteps
    % Calculate distances between particles and their separations in x and y directions
    dx = repmat(xcurr', [1, N]) - repmat(xcurr, [N, 1]);
    dy = repmat(ycurr', [1, N]) - repmat(ycurr, [N, 1]);
    % Apply periodic boundary conditions to distances
    dx = dx - L * round(dx / L);
    dy = dy - L * round(dy / L);
    % Calculate the squared distances between particles
    dr2 = dx.^2 + dy.^2;
    [row, col, DR2] = find(triu(dr2, 1));
    linearIndex = (col - 1) * N + row;
    dx = dx(linearIndex);
    dy = dy(linearIndex);
    invDR2 = zeros(size(DR2));
    cfIndex = (DR2 < drcut2);
    invDR2(cfIndex) = 1. / DR2(cfIndex);
    f = 48 * invDR2.^4 .* (invDR2.^3 - 0.5);
    fx = f .* dx;
    fy = f .* dy;
    fx = full(sparse(row, col, fx, N, N));
    fy = full(sparse(row, col, fy, N, N));
    fx = sum(-fx' + fx, 2);
    fy = sum(-fy' + fy, 2);
    for i = 1:N
        xnew(i) = 2 * xcurr(i) - xprev(i) + fx(i) * dt^2;
        ynew(i) = 2 * ycurr(i) - yprev(i) + fy(i) * dt^2;
        if xnew(i) < 0
            xnew(i) = xnew(i) + L;
        end
        if xnew(i) > L
            xnew(i) = xnew(i) - L;
        end
        if ynew(i) < 0
            ynew(i) = ynew(i) + L;
        end
        if ynew(i) > L
            ynew(i) = ynew(i) - L;
        end
        vx(i) = (xnew(i) - xprev(i)) / (2 * dt);
        vy(i) = (ynew(i) - yprev(i)) / (2 * dt);
    end
    for i = 1:N
        xprev(i) = xcurr(i);
        yprev(i) = ycurr(i);
        xcurr(i) = xnew(i);
        ycurr(i) = ynew(i);
    end

    % Plot the current configuration and save it as a PNG file
    fig = figure('visible', 'off');
    scatter(xcurr, ycurr, 500, 'b.');
    title({'2D Molecular Dynamics Simulation (Initial Velocity = 1)'; ['step = ', num2str(n)]});
    axis([0 L 0 L]);
    drawnow;
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);
    if mod(n, 10) == 0
        saveas(fig, fullfile('2D_Initial_Velocity_1', sprintf('step%d.png', n)))
    end
    close(gcf);
end
% Finish video
close(outputVideo);

% Speed distribution
speed = sqrt(vx.^2 + vy.^2);
figure;
histfit(speed);
title('Speed Distribution');
xlabel('Speed');
ylabel('Probability');
saveas(gcf, fullfile('2D_Initial_Velocity_1', sprintf('Speed Distribution.png')))

% Velocity distribution
figure;
histfit(vx);
title('Velocity Distribution');
xlabel('Vx');
ylabel('Probability');
saveas(gcf, fullfile('2D_Initial_Velocity_1', sprintf('Velocity Distribution(x-axis).png')))

% Velocity distribution
figure;
histfit(vy);
title('Velocity Distribution');
xlabel('Vy');
ylabel('Probability');
saveas(gcf, fullfile('2D_Initial_Velocity_1', sprintf('Velocity Distribution(y-axis).png')))