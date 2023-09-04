% Title: Random Motion of Multiple Blue Discs in 1D
% Author: Gopi Reddy Rajendran
% Description:
% This MATLAB code simulates the random motion of multiple blue discs in 1D using a simple random walk model.
% The initial positions of the discs are set using linspace and are plotted as blue dots with MarkerSize of 50.
% Each disc then moves to a new position that is randomly chosen within a range of 'delta_r' units from its current position. The new positions are plotted as red dots with MarkerSize of 50, and the plot is updated using drawnow to show the motion of the discs in real-time.
% Code:
% Clear the workspace and command window
clear;
clc;
% Set the number of discs
N = 10;
% Set the range of positions for the discs
position = linspace(-20, 20, N);
% Set the range of displacement for the discs
delta_r = 0.8;
% Plot the initial positions of the discs as blue dots
plot(position, 0, 'b.', 'MarkerSize', 50);
hold on;
% Loop through each disc and update its position to a new random position
for i = 1:N
new_position(i) = position(i) + 2*(rand-0.5)*delta_r;
% Plot the new position of the disc as a red dot
plot(new_position(i), 0, 'r.', 'MarkerSize', 50);
% Set plot limits
axis([-25 25 -0.5 0.5]);
title('Random Motion of multiple Blue Disc'); % Title for the plot
% Update the plot to show the motion of the disc in real-time
drawnow;
end
% Release the plot hold
hold off;
% End of code