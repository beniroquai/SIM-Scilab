//Image Reconstruction 

path = './acq';
i=0;
maxpic = 9

ki = 2*%pi*(fc)/n;





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
n = 400;
count = 9

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
    
    
    //Patterngeneration
    k = rottheta*[ki; 0];

    kx = k(1); ky = k(2);
    I(:,:,(itheta-1)*3+iphi) = (1-cos(kx*X1+ky*Y1+phi(iphi)))/2;

    pattern = abs(I(:,:,(itheta-1)*3+iphi));




// Reconstruction



    
    //reading set of taken pictures
    data(:;:;i) = ReadImage(path+string(x)+'.jpg');
    
    //bringing pictures in reciprocal space and shift them to the 
    //moves the zero frequency component to the center of the spectrum
    DIbar(:;:;i) = fftshift(fft2(ifftshift(data(:;:;i)));
    //Illumination Patterned Object in Reciprocal Space

i=i+1;
end

nx = kx*n/(2*pi); ny = ky*n/(2*pi);
ind = [(0:(nphases-1)/2) -(nphases-1)/2:1:-1];
xx=repmat(x,n,1);
yy=repmat(y,1,n);



fo  r j = 1:nphases
temp_separated = zeros(n,n,nphases);
Ir(:,:)=exp(-1i*(kx*ind(j)*xx+ky*ind(j)*yy));
for k = 1:nphases
temp_separated(:,:,k) = inv_phase_matrix(j,k).*DIbars(:,:,k);
sp(:,:,(itheta-1)*3+j) = sp(:,:,(itheta-1)*3+j)+temp_separated(:,:,k);
end
//figure; colormap('gray'); imagesc(abs(sp(:,:,(itheta-1)*3+j))); axis('square');title('1')
hv(:,:) = OTF(n,n,-nx*ind(j),-ny*ind(j),fc);

//shift the OTF by taking inverse Fourier transform and exponential factors

replc(:,:) = fft2(ifftshift(ifft2(sp(:,:,(itheta-1)*3+j)).*Ir(:,:)));
scalefactor = abs(cos((j-1)*pi/3));

figure; colormap('gray'); imagesc(abs(replc(:,:))); axis('square');title('2')
rc(:,:,(itheta-1)*3+j) = replc(:,:).*(scalefactor*conj(hv(:,:)));
hs = hs + abs(scalefactor*hv(:,:)).^2;
end
end




// Deconvolution and reconstruction with a Wiener

for t = 1:nphases*nphases
dr = dr+rc(:,:,t)./ ( hs + .005*length(itheta)*(.0000001)^2);
end
figure; colormap('gray'); imagesc(abs(dr(:,:))); axis('square'); title('Reconstruction of The Object in Reciprocal Space')

//Triangular function

[k_x, k_y]=meshgrid(-n/2+1:n/2, -n/2+1:n/2);
k_r = sqrt(k_x.^2+k_y.^2);
k_max = .9*.9*fc*((nphases-1)/2+1);
bhs = cos(pi*k_r/(2*k_max));
indi = find( k_r > k_max );
bhs(indi) = 0;
figure; colormap('gray'); imagesc(abs(dr.*bhs)); axis('square'); title('Apodization of The Object in Reciprocal Space')
drr = dr.*bhs;
fimage = ifft2(ifftshift(drr));

figure; colormap('gray'); imagesc(abs(fimage)); axis('square'); title('Reconstruction of The Object in Real Space')



