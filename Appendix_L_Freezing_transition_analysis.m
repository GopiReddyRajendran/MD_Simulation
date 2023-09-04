% Author: Gopi Reddy Rajendran
% Description:
% This code simulates a 2D molecular dynamics system with
% Lennard-Jones interactions and records various properties
% such as particle positions, velocities, energies, and
% temperature over time.

% Clear the workspace and command window
clc; % clear command window
clear; % clear workspace
close all; % close all figure windows

% Simulation parameters
N = 100; % number of particles
d = 1; % distance of each particle
dt = 0.01; % time step for simulation
dr = 0.1; % initial displacement of particles from grid points
V = 0; % initial velocity of particles
drcut = 3; % cutoff distance for Lennard-Jones potential
T = 50; % total simulation time
nSteps = T/dt; % number of time steps in simulation
kB = 1.380649e-23;  % Boltzmann constant (J/K)

% Round N to the nearest square number to create a square grid of particles
N = round(sqrt(N))^2;

% Determine the size of the simulation box
L = round(d*sqrt(N));

% Calculate the cutoff distance squared for efficiency
drcut2 = drcut*drcut;

% Generate initial positions and velocities for particles
% within a small range of the grid points
[gx,gy] = meshgrid(d/2:d:L-(d/2),d/2:d:L-(d/2));
gx = gx(1:N);
gy = gy(1:N);
for i = 1:N
    xcurr(i) = gx(i) + 2*(rand-0.5)*dr;
    ycurr(i) = gy(i) + 2*(rand-0.5)*dr;
    vx0(i) = 2*(rand-0.5)*V;
    vy0(i) = 2*(rand-0.5)*V;
    xprev(i) = xcurr(i) - vx0(i)*dt;
    yprev(i) = ycurr(i) - vy0(i)*dt;
    XatomsOrig(i) = xcurr(i);
    YatomsOrig(i) = ycurr(i);
end

% Create the "plots" folder if it doesn't exist
if ~exist('2D_Initial_Velocity_0', 'dir')
    mkdir('2D_Initial_Velocity_0')
end

% Create VideoWriter object
outputVideo = VideoWriter('2D_Initial_Velocity_0.avi');
outputVideo.FrameRate = 30; % Frames per second
open(outputVideo);

% Initial plot
fig = figure('visible', 'off');
scatter(xcurr, ycurr, 500, 'b.');
title({'2D Molecular Dynamics Simulation (Initial Velocity = 0)'; 'step = 0'});
axis([0 L 0 L]);
drawnow;
frame = getframe(gcf);
writeVideo(outputVideo, frame);
saveas(fig, fullfile('2D_Initial_Velocity_0', sprintf('Initial.png')))
close(gcf);

% Main simulation loop
for n = 1:nSteps
    % Calculate distances between particles and their separations in x and y directions
    dx = repmat(xcurr',[1,N])-repmat(xcurr,[N,1]);
    dy = repmat(ycurr',[1,N])-repmat(ycurr,[N,1]);
    
    % Apply periodic boundary conditions to distances
    dx = dx - L*round(dx/L);
    dy = dy - L*round(dy/L);
    
    % Calculate the squared distances between particles
    dr2 = dx.^2 + dy.^2;
    [row,col,DR2] = find(triu(dr2,1));
    linearIndex = (col-1)*N+row;
    dx = dx(linearIndex);
    dy = dy(linearIndex);
    invDR2 = zeros(size(DR2));
    cfIndex = (DR2<drcut2);
    invDR2(cfIndex) = 1./DR2(cfIndex);
    f = 48*invDR2.^4.*(invDR2.^3 - 0.5);
    fx = f.*dx;
    fy = f.*dy;
    fx = full(sparse(row,col,fx,N,N));
    fy = full(sparse(row,col,fy,N,N));
    fx = sum(-fx'+fx,2);
    fy = sum(-fy'+fy,2);
    
    % Update particle positions and velocities
    for i = 1:N
        xnew(i) = 2*xcurr(i) - xprev(i) + fx(i)*dt^2;
        ynew(i) = 2*ycurr(i) - yprev(i) + fy(i)*dt^2;
        
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
        
        MSD(i) = (((xnew(i)-XatomsOrig(i))^2)+((ynew(i)-YatomsOrig(i))^2));
        
        vx(i) = (xnew(i) - xprev(i))/(2*dt);
        vy(i) = (ynew(i) - yprev(i))/(2*dt);
        KE(i) = 0.5 * ((vx(i)^2) + (vy(i)^2));
        PE(i) = sqrt(((sum(fx(i)))^2)+((sum(fy(i)))^2));
    end
    
    AvgMSD(n) = (sum(MSD)/N);
    
    % Average Energy Calculations
    KineticEnergy(n) = (sum(KE)/N);
    PotentialEnergy(n) = (sum(PE)/N);
    TotalEnergy(n) = KineticEnergy(n) + PotentialEnergy(n);
    Temp(n) = (2 * KineticEnergy(n)) / (3 * N * kB);
    
    % Update particle positions and velocities for the next time step
    for i = 1:N
        xprev(i) = xcurr(i);
        yprev(i) = ycurr(i);
        xcurr(i) = xnew(i);
        ycurr(i) = ynew(i);
    end
    
    % Plot the current configuration and save it as a PNG file
    fig = figure('visible', 'off');
    scatter(xcurr, ycurr, 500, 'b.');
    title({'2D Molecular Dynamics Simulation (Initial Velocity = 0)'; ['step = ', num2str(n)]});
    axis([0 L 0 L]);
    drawnow;
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);
    if mod(n, 10) == 0
        saveas(fig, fullfile('2D_Initial_Velocity_0', sprintf('step%d.png', n)))
    end
    close(gcf);
end

% Finish video
close(outputVideo);

% Plotting the temperature versus time
time = (1:nSteps);
figure;
plot(time, KineticEnergy);
xlabel('Time');
ylabel('Kinetic Energy');
title('Kinetic Energy versus Time');
saveas(gcf, fullfile('2D_Initial_Velocity_0', sprintf('Kinetic Energy.png')))

% Plotting the potential energy versus time
figure;
plot(time, PotentialEnergy);
xlabel('Time');
ylabel('Potential Energy');
title('Potential Energy versus Time');
saveas(gcf, fullfile('2D_Initial_Velocity_0', sprintf('Potential Energy.png')))

% Plotting the total energy versus time
figure;
plot(time, TotalEnergy);
xlabel('Time');
ylabel('Total Energy');
title('Total Energy versus Time');
saveas(gcf, fullfile('2D_Initial_Velocity_0', sprintf('Total Energy.png')))

% Plotting the displacement versus time
figure;
plot(time, AvgMSD);
xlabel('Time');
ylabel('Displacement');
title('Displacement versus Time');
saveas(gcf, fullfile('2D_Initial_Velocity_0', sprintf('Displacement.png')))

% Plotting the temperature versus time
figure;
plot(time, Temp);
xlabel('Time');
ylabel('Temperature');
title('Temperature versus Time');
saveas(gcf, fullfile('2D_Initial_Velocity_0', sprintf('Temperature.png')))
