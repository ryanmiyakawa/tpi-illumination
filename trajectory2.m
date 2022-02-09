%trajectory

sigmamin=0.2; sigmamax=0.8; pen=0.1; nsamples=10000;   %from function call
sigmamin=sigmamin+0.5*pen;sigmamax=sigmamax-0.5*pen;
deltasigma=sigmamax-sigmamin;
n=round(10*deltasigma)+2;               %number of rings
r=sigmamin:deltasigma/(n-1):sigmamax;
dphi=pen./r;
arclength=(2*pi-dphi).*r; step=sum(arclength)/nsamples; nsteps=round(arclength/step);

c=1;
phi=0:(2*pi-dphi(c))/(nsteps(1)-1):2*pi-dphi(c);phi=phi+1.5*pi+0.5*dphi(c);
u=r(c).*cos(phi);
v=r(c).*sin(phi);

for c=2:n
    
   phi=0:(2*pi-dphi(c))/(nsteps(c)-1):2*pi-dphi(c);phi=phi+1.5*pi+0.5*dphi(c);
   x=r(c).*cos(phi);
   y=r(c).*sin(phi);
   u=[u x];
   v=[v y];
   
end

%unit circle

phi=0:0.1:2*pi;
x=cos(phi);
y=sin(phi);

figure(1);plot(u,v,x,y);axis image;

scanpath= [u' v'];

%save('scanpath.mat','scanpath');