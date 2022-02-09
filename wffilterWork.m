
load('scanpath.mat');

n=length(scanpath);

x=scanpath(:,1)'; y=scanpath(:,2)';

figure(1); plot(x,y,'.'); axis image;

x=x+10; y=y+10;

s=0.4; %scaling factor

looptime=1; % [ s ]

f_cutoff=1000; % [ Hz ]

cutoff=s*f_cutoff*looptime;

t=1:n; t=t-n/2;

filter=exp(-(t/cutoff).^2);

figure(2); plot(filter);
fx=fftshift(fft(x));
fy=fftshift(fft(y));
fx=fx.*filter;
fy=fy.*filter;
u=abs(ifft(fx))-10;
v=abs(ifft(fy))-10;
figure(3); plot(u,v,'.'); axis image;