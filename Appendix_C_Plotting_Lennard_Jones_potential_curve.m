% Title: Lennard-Jones Potential Curve
% Author: Gopi Reddy Rajendran
% Description:
% This MATLAB code generates a plot of the Lennard-Jones potential energy curve.
% The Lennard-Jones potential is given by the formula: U(r) = 4 * epsilon * ((sigma/r)^12 - (sigma/r)^6)
% where epsilon is the depth of the potential well, sigma is the distance at which the potential is zero, and r is the distance between the two particles.
% In this example, we plot the potential energy (U) as a function of the reduced distance (r/sigma) and the reduced potential energy (U/epsilon), where sigma and epsilon are constants.
% Code:
% Clear the workspace and command window
clear;
clc;
% Parameters for the Lennard-Jones potential
epsilon = 1.0; % depth of the potential well
sigma = 1.0; % distance at which the potential is zero
% Generate a range of distances (r) from 0.5 to 5.0 using linspace
r = linspace(0.5,5);
% Calculate the Lennard-Jones potential energy for each distance in the range
U = 4 * epsilon * ((sigma ./ r).^12 - (sigma ./ r).^6);
% Plot the Lennard-Jones potential energy versus the reduced distance (r/sigma) and the reduced potential energy (U/epsilon)
plot(r/sigma, U/epsilon, 'b-', 'LineWidth', 2);
xlabel('V(r)/ε'); % Label for the x-axis
ylabel('r/σ'); % Label for the y-axis
axis([0 4 -2 2]);
title('Lennard-Jones Potential Energy Curve'); % Title for the plot
grid on; % Add a grid to the plot
% End of code
