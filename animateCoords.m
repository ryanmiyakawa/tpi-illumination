function animateCoords(idx)

% Load file generated from c++ code
M = csvread('illumination_coords.txt');

x = M(:,1);
y = M(:,2);


th = linspace(0, 2*pi, 100);
for k = 1:length(x)*2
    plot(x(mod(k, length(x))), y(mod(k, length(x))), '.-', cos(th), sin(th), 'r');
    pause(0.05);
end
axis image

set(gcf, 'color', 'w')