//1. Reconstruction.sci
//
// Clear previous variables and graphics windows
//

clear;
stacksize(100000000);

// Set Camera Parameters PIKE_F_145B_C_fiber.pdf

cpixel_width = 1388;     //  Number of Pixel width
cpixel_height = 1038;     //  Number of Pixel height 
pixelsize =  6.45;       //  Pixel_Size [mu m]
cs_width = cpixel_width * pixelsize;   // width of Camera-Sensor [mm]
cs_height = cpixel_height * pixelsize;   // hight of Camera-Sensor [mm]


// Set Projector Parameters
ps_width = 12;   // width of Camera-Sensor [mm]
ps_height = 7;    // hight of Camera-Sensor [mm]
ppixel_width = 480;     //  Number of Pixel width
ppixel_height = 320;     //  Number of Pixel height 
ppixel_size = 7,6;       //  Pixel_Size [mu m]

// Set Parameters for Grating
start_angle = 11;   // Initial phase (!=0 for better result)
nphases = 3; // number of phases
nangles = 3; // number of angles
pixelperiode = 5;   // How many Pixels for a periode of a grating


// Set Lens Parameters
na = .25;   // numerical aperture (1.4 for oil immersion lens)
mag = 10;     // magnification of Lens (@ TL=200mm)
lambda = 1000; // wavelength in nanometers
rho = 0,61*lambda/na;

fc = 10; 


// Set Image Parameters
datasets = 1;
impath = 'C:\Users\Bene\Dropbox\Dokumente\Optical Design\Superresolution\Aufnahmen\Aufnahme_Matlab\';
extension = '.png';

//
// Pattern (theta) and phase angles (phi)
//
theta=linspace(0,%pi,nangles+1);
theta=theta(1:nangles);
phi=linspace(0,2*%pi,nphases+1);
phi=phi(1:nphases);

// LOAD IMAGES 
//I(:,:,(itheta-1)*3+iphi) = (1-cos(kx*X1+ky*Y1+phi(iphi)))/2;


file = impath + string(1) + '_' + string(1) + extension;
pixelnumber =  size(imread(file));
n = pixelnumber(1);
m = pixelnumber(2);


// Parameter for Reconstruction of the superresolved image
rc = zeros(n,n,nphases^2);
sp = zeros(n,n,nphases*nphases);
I = zeros(n,n,nphases*nphases);
rimage = zeros(n,n,3);
hv = zeros(n,n);
replc = zeros(n,n);
dr = zeros(n,n);

x=1:n;
y=(1:m)';

raw = zeros(n,m,datasets*nphases*nangles);
spectrum = zeros(n,m,datasets*nphases*nangles);


ki = 2*%pi*(fc)/n; //???????????????



//
// Generate arrays X1 and Y1 with n rows and n columns. For array X1 (Y1),
// rows (columns) are identical and columns (rows) have values ranging from 1 to n.
//
// Example for n=5:
//
// X1=[1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5]
// Y1=[1 1 1 1 1; 2 2 2 2 2; 3 3 3 3 3; 4 4 4 4 4; 5 5 5 5 5]
//
[X1,Y1] = meshgrid(1:n,m:n); //?????????




//
// Information components separation
//
//Use this to store the values of DIbar, which is the B matrix ( AX = B ) 
DIbars = zeros(n,m,nphases);
//The A matrix of the equations( AX = B ), it depends on the shift phases
//phi
phase_matrix = [1 1 1;1 exp(%i*2*%pi/3) exp(-%i*2*%pi/3);1 exp(%i*4*%pi/3) exp(-%i*4*%pi/3);];
//Inversed A matrix
inv_phase_matrix = inv(phase_matrix);



for j = 1:nangles
    k= [cos(j) -sin(j); sin(j) cos(j) ]*[ki; 0]; //?????????????????????
    kx = k(1);
    ky = k(2);
    
    nx = kx*n/(2*%pi); ny = ky*m/(2*%pi);
    ind = [(0:(nphases-1)/2) -(nphases-1)/2:1:-1];
    xx=repmat(x,n,1);
    yy=repmat(y,1,n);
    
    for i = 1:nphases
        file = impath + string(i) + '_' + string(j) + extension;
        I=imread(file)
        Igray=RGB2Gray(I)
        //Igray=Igray./max(max(Igray));
        raw(:,:,i)=  Igray;
        fourier = fftshift(fftw(ifftshift(raw(:,:,1))))
        
        //Show/Save FFT of Image
        norm_spectrum = max(max(abs(fftw(raw(:,:,1)))));
        imshow(abs(fourier)./norm_spectrum * 255);
        imwrite(abs((fourier)./norm_spectrum * 255),impath + string(i) + '_' + string(j) + '_fft'  + extension)
        //-------
    end
    
    nx = kx*n/(2*%pi); ny = ky*m/(2*%pi);
    ind = [(0:(nphases-1)/2) -(nphases-1)/2:1:-1];
    xx=repmat(x,n,1);
    yy=repmat(y,1,n);


// This part is solving the equations to separate the components

for jj= 1:nphases
temp_separated = zeros(n,m,nphases);
Ir(:,:)=exp(-%i*(kx*ind(jj)*xx+ky*ind(jj)*yy));
for kk = 1:nphases
temp_separated(:,:,kk) = inv_phase_matrix(jj,kk).*DIbars(:,:,kk);
sp(:,:,(itheta-1)*3+jj) = sp(:,:,(itheta-1)*3+jj)+temp_separated(:,:,kk);
end

//%shift the OTF by taking inverse Fourier transform and exponential factors

replc(:,:) = fft2(ifftshift(ifft2(sp(:,:,(itheta-1)*3+j)).*Ir(:,:)));
scalefactor = abs(cos((j-1)*pi/3));


rc(:,:,(itheta-1)*3+j) = replc(:,:).*(scalefactor*conj(hv(:,:)));
hs = hs + abs(scalefactor*hv(:,:)).^2;
end
end

end


