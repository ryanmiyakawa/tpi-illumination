lobes = 5;
alpha = 1.7;
a = 4;
b = 3;

N = 1000;
t = linspace(0, 2*pi * (lobes + 1), N+1);
t = t(1:end-1);


x = cos(t) + 1/alpha * cos(t/(lobes + 1));
y = sin(t) + 1/alpha * sin(t/(lobes + 1));

plot(x,y,'.');
axis image