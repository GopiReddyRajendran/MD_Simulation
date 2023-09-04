% Title: Ten equally spaced blue disc in a row
% Author: Gopi Reddy Rajendran
% Description:
% This MATLAB code generates a scatter plot of data points with x and y coordinates using the scatter function.
% The scatter function is used to create a scatter plot, where x and y are the coordinates of the data points, and 1000 specifies the size of the markers, and 'b.' specifies the color and marker type of the markers (in this case, blue dots).
% The x values are generated using the linspace function to create a linearly spaced vector from 1 to 10 with 10 points.
% The y values are set to 0 for all data points.
% The axis function is used to set the limits of the plot to x-axis from 0 to 11 and y-axis from -0.5 to 0.5.
% Code:
% Clear the workspace and command window
clear;
clc;
% Initialize variables
x = linspace(1, 10, 10);
y = 0;
% Generate scatter plot
scatter(x, y, 1000, 'b.');
% Set plot limits
axis([0 11 -0.5 0.5]);
% Add title and labels to the plot
title('Ten equally spaced blue disc in a row');
% End of code
