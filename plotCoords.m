% Load file generated from c++ code
M = csvread('illumination_coords.txt');

x = M(:,1);
y = M(:,2);

th = linspace(0, 2*pi, 100);

plot(x, y, '.-', cos(th), sin(th), 'r');
axis image

set(gcf, 'color', 'w')