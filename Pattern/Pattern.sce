//Generation of Pattern

//1. linSI.m
//
// Clear previous variables and graphics windows
//
clear ;
//
// Set parameters
//

impath = 'C:\Users\Bene\Dropbox\Dokumente\Optical Design\Superresolution\Scilab\Pattern\'

// Microscope Objective

M_Microscope = 10;       // Magnification @TL 200mm
NA = 2;
lambda = 555e-9;
f_mo = 200e-3/M_Microscope;

// Magnification
f_TL = 400e-3  ;         // Focal Length Tube Lens
f_PL = 400e-3 ;         // Focal Length Projection Lens

mag_p = f_mo/f_PL;      // Demagnification of DLP Display into Probe Plane
mag_i = f_TL/f_mo;      // Magnification of Probe into Camera


// Periode
period = 10;             // How many pixel form a periode




rho_p = 0.61*lambda/NA;
fc = 1/rho_p;


// Set Camera Parameters PIKE_F_145B_C_fiber.pdf

pixel_camera_width = 1388;     //  Number of Pixel width
pixel_camera_height = 1038;     //  Number of Pixel height 
pixelsize_camera =  6.45e-6;       //  Pixel_Size [mu m]
cs_width = pixelsize_camera * pixel_camera_width;   // width of Camera-Sensor [mm]
cs_height = pixelsize_camera * pixel_camera_height;   // hight of Camera-Sensor [mm]


// Set Projector Parameters
ps_width = 1e-3;   // width of Camera-Sensor [mm]
ps_height = 72e-3;    // hight of Camera-Sensor [mm]
pixel_projector_width = 320;// 480;     //  Number of Pixel width
pixel_projector_height = 320;     //  Number of Pix<<l height 
pixelsize_projector = 7.6e-6;       //  Pixel_Size [mu m]

// find best Period to fullfill resolution enhancement

period_best = ceil(((rho_p/mag_p)/pixelsize_projector));
disp(period_best, 'Best Period for DLP to fit: ')
fc_dlp  = floor(1/(period_best*pixelsize_projector))

// Set number of pattern, angle, rotation

nphases = 3; // number of phases
nangles = 3; // number of angles


// ------------------ START -----------

n = pixel_projector_width;
m = pixel_projector_height;

theta=linspace(0,%pi,nangles+1);
theta=theta(1:nangles);
phi=linspace(0,2*%pi,nphases+1);
phi=phi(1:nphases);

I = zeros(m,n,nphases*nangles);



[X1,Y1] = meshgrid(1:n,1:m);

ki = 2*%pi*(fc_dlp)/1000;


for itheta=1:nphases
for iphi=1:nphases
    rottheta= [cos(theta(itheta)) -sin(theta(itheta));
sin(theta(itheta)) cos(theta(itheta)) ];
    
k = rottheta*[ki; 0];

kx = k(1); ky = k(2);
I(:,:,(itheta-1)*3+iphi) = (1/2+squarewave((kx*X1+ky*Y1+phi(iphi)))/2);
imshow(I(:,:,(itheta-1)*3+iphi))

pattern = I(:,:,(itheta-1)*nangles+iphi);
WriteImage(pattern, impath+'period_'+string(period_best)+'_phase'+string(itheta)+'_'+string(iphi)+'.jpg');


end
end

