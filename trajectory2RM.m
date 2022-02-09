%trajectory

sigmamin    = 0.2; 
sigmamax    = 0.8; 
pen         = 0.2; % ?? maybe the effective width of a sample point, used to compute ring radii

nsamples    = 10000;   %from function call

sigmamin    = sigmamin + 0.5*pen;
sigmamax    = sigmamax - 0.5*pen;

deltasigma  = sigmamax - sigmamin;

numRings    = round(10*deltasigma) + 2;  % number of rings

ringRadii   = sigmamin:deltasigma/(numRings-1):sigmamax;
dphi        = pen./ringRadii;

arclength   = (2*pi-dphi).*ringRadii; % part of ring before skipping to next ring
step        = sum(arclength)/nsamples; 
nsteps      = round(arclength/step);

c = 1;
phi         =  0 : (2*pi-dphi(c))/(nsteps(1)-1) : 2*pi-dphi(c);

phi=phi+1.5*pi+0.5*dphi(c);
u=ringRadii(c).*cos(phi);
v=ringRadii(c).*sin(phi);

for c=2:numRings
    
   phi=0:(2*pi-dphi(c))/(nsteps(c)-1):2*pi-dphi(c);phi=phi+1.5*pi+0.5*dphi(c);
   x=ringRadii(c).*cos(phi);
   y=ringRadii(c).*sin(phi);
   u=[u x];
   v=[v y];
   
end

%unit circle

phi=0:0.1:2*pi;
x=cos(phi);
y=sin(phi);

figure(1);plot(u,v, '.', x,y);axis image;

scanpath= [u' v'];

%save('scanpath.mat','scanpath');