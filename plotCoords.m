function plotCoords(idx)

% Load file generated from c++ code
M = csvread('illumination_coords.txt');

x = M(:,1);
y = M(:,2);

if nargin== 1
    x = x(idx);
    y = y(idx);
end

th = linspace(0, 2*pi, 100);

plot(x, y, '.-', cos(th), sin(th), 'r');
axis image

set(gcf, 'color', 'w')