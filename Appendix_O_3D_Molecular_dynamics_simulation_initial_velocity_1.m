% Author: Gopi Reddy Rajendran
% Description: This MATLAB code simulates a 3D molecular dynamics (MD)
% simulation with initial velocity as 1 of particles interacting with a Lennard-Jones potential.

% Clear the workspace and command window
clc; % clear command window
clear; % clear workspace
close all; % close all figure windows

% Simulation parameters
N = 100; % number of particles
d = 1; % distance of each particle
dt = 0.01; % time step for simulation
dr = 0.1; % initial displacement of particles from grid points
V = 1; % initial velocity of particles
drcut = 4; % cutoff distance for Lennard-Jones potential
T = 50; % total simulation time
nSteps = T/dt; % number of time steps in simulation
plotDelay = 0.05; % time delay between plot updatesclc;

% MD Initialising 
% Round N to the nearest square number to create a square grid of particles
N = round((N^(1/3)))^3;

% Determine the size of the simulation box
L = round(d*(N^(1/3)));

% Calculate the cutoff distance squared for efficiency
drcut2 = drcut*drcut;

xprev = zeros(1,N);     yprev = zeros(1,N);     zprev = zeros(1,N);
xcurr = zeros(1,N);     ycurr = zeros(1,N);     zcurr = zeros(1,N);
xnew = zeros(1,N);     ynew = zeros(1,N);     znew = zeros(1,N);


vx0 = zeros(1,N);     vy0 = zeros(1,N);     vz0 = zeros(1,N);
vx = zeros(1,N);     vy = zeros(1,N);     vz = zeros(1,N);

% Generate initial positions and velocities for particles
% within a small range of the grid points
[gx,gy,gz] = meshgrid(d/2:d:L-(d/2),d/2:d:L-(d/2),d/2:d:L-(d/2));
gx = gx(1:N);
gy = gy(1:N);
gz = gz(1:N);
for i = 1:N
    xcurr(i) = gx(i) + 2*(rand-0.5)*dr;
    ycurr(i) = gy(i) + 2*(rand-0.5)*dr;
    zcurr(i) = gz(i) + 2*(rand-0.5)*dr;
    vx0(i) = 2*(rand-0.5)*V;
    vy0(i) = 2*(rand-0.5)*V;
    vz0(i) = 2*(rand-0.5)*V;
    xprev(i) = xcurr(i) - vx0(i)*dt;
    yprev(i) = ycurr(i) - vy0(i)*dt;
    zprev(i) = zcurr(i) - vz0(i)*dt;
end

% Create the "plots" folder if it doesn't exist
if ~exist('3D_1', 'dir')
    mkdir('3D_1')
end
% Create VideoWriter object
outputVideo = VideoWriter('3D1.avi');
outputVideo.FrameRate = 10; % Frames per second
open(outputVideo);

fig = figure('visible', 'off');
scatter3(xcurr,ycurr,zcurr,500,'b.');
title({'3D Molecular Dynamics Simulation';'step = 0'});
axis([0 L 0 L 0 L]);
drawnow;
frame = getframe(gcf);
writeVideo(outputVideo, frame);
saveas(fig, fullfile('3D_1', sprintf('Initial.png')))
close(gcf);

% MD Simulation
for n = 1:nSteps
    
    dx = repmat(xcurr',[1,N])-repmat(xcurr,[N,1]);
    dy = repmat(ycurr',[1,N])-repmat(ycurr,[N,1]);
    dz = repmat(zcurr',[1,N])-repmat(zcurr,[N,1]);

    dx = dx - L*round(dx/L);
    dy = dy - L*round(dy/L);
    dz = dz - L*round(dz/L);

    dr2 = dx.^2 + dy.^2 + dz.^2;

    [row,col,DR2] = find(triu(dr2,1));
    linearIndex = (col-1)*N+row;
    dx = dx(linearIndex);
    dy = dy(linearIndex);
    dz = dz(linearIndex);

    invDR2 = zeros(size(DR2));
    cfIndex = (DR2<drcut2);
    invDR2(cfIndex) = 1./DR2(cfIndex);
    f = 48*invDR2.^4.*(invDR2.^3 - 0.5);

    fx = f.*dx;
    fy = f.*dy;
    fz = f.*dz;

    fx = full(sparse(row,col,fx,N,N));
    fy = full(sparse(row,col,fy,N,N));
    fz = full(sparse(row,col,fz,N,N));

    fx = sum(-fx'+fx,2);
    fy = sum(-fy'+fy,2);
    fz = sum(-fz'+fz,2);
    for i = 1:N
        xnew(i) = 2*xcurr(i) - xprev(i) + fx(i)*dt^2;
        ynew(i) = 2*ycurr(i) - yprev(i) + fy(i)*dt^2;
        znew(i) = 2*zcurr(i) - zprev(i) + fz(i)*dt^2;
        
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
        if znew(i) < 0
            znew(i) = znew(i) + L;
        end
        if znew(i) > L
            znew(i) = znew(i) - L;
        end
        
        vx(i) = (xnew(i) - xprev(i))/(2*dt);
        vy(i) = (ynew(i) - yprev(i))/(2*dt);
        vz(i) = (znew(i) - zprev(i))/(2*dt);
    end
    
    for i = 1:N
           xprev(i) = xcurr(i);
           yprev(i) = ycurr(i);
           zprev(i) = zcurr(i);

           xcurr(i) = xnew(i);
           ycurr(i) = ynew(i);
           zcurr(i) = znew(i);
            
    end
    
    % Plot the current configuration and save it as a PNG file
    fig = figure('visible', 'off');
    scatter3(xcurr,ycurr,zcurr,500,'b.');
    title({'3D Molecular Dynamics Simulation';['step = ',num2str(n)]});
    axis([0 L 0 L 0 L]);
    drawnow;
    frame = getframe(gcf);
    writeVideo(outputVideo, frame);
    if mod(n,10) == 0
        saveas(fig, fullfile('3D_1', sprintf('step%d.png', n)))
    end
close(gcf);    
end
% Finish video
close(outputVideo);