% Title: Random Motion of a Blue Disc in 1D
% Author: Gopi Reddy Rajendran
% Description:
% This MATLAB code simulates the random motion of a blue disc in 1D using a simple random walk model.
% The initial position of the disc is set to 'position' and is plotted as a blue dot with MarkerSize of 50.
% The disc then moves to a new position that is randomly chosen within a range of 'delta_r' units from its current position. The new position is plotted as a red dot with MarkerSize of 50, and the plot is updated using drawnow to show the motion of the disc in real-time.
% Code:
% Clear the workspace and command window
clear;
clc;
% Set the initial position of the disc
position = 5;
% Set the range of displacement for the disc
delta_r = 0.8;
% Plot the initial position of the disc as a blue dot
plot(position, 1, 'b.', 'MarkerSize', 50);
hold on;
% Update the position of the disc to a new random position
new_position = position + 2*(rand-0.5)*delta_r;
% Plot the new position of the disc as a red dot
plot(new_position, 1, 'r.', 'MarkerSize', 50);
% Set plot limits
axis ([4 6 0 2])
title('Random Motion of a Blue Disc'); % Title for the plot
% Update the plot to show the motion of the disc in real-time
drawnow;
% Release the plot hold
hold off;
% End of code
