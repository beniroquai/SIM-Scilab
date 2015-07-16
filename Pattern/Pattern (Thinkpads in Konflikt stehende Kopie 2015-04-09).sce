//Generation of Pattern

//1. linSI.m
//
// Clear previous variables and graphics windows
//
clear ;
//plotlibmode

//
// Set parameters
//
fc = 40;
n = 480;
i= 1;

nphases = 3; // number of phases
nangles = 3; // number of angles

theta=linspace(0,%pi,nangles+1);
theta=theta(1:nangles);
phi=linspace(0,2*%pi,nphases+1);
phi=phi(1:nphases);

I = zeros(n,n,nphases*nphases);

[X1,Y1] = meshgrid(1:n,1:n);

ki = 2*%pi*(fc)/n;


for itheta=1:nphases
for iphi=1:nphases
    rottheta= [cos(theta(itheta)) -sin(theta(itheta));
sin(theta(itheta)) cos(theta(itheta)) ];
    
    k = rottheta*[ki; 0];

kx = k(1); ky = k(2);
I(:,:,(itheta-1)*3+iphi) = (1-cos(kx*X1+ky*Y1+phi(iphi)))/2;

pattern = abs(I(:,:,(itheta-1)*3+iphi));
WriteImage(pattern, 'phase'+string(itheta)+'_'+string(iphi)+'.jpg');
end
end


