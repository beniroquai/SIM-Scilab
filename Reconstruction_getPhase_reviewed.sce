//1. Reconstruction.sci
//
// Clear previous variables and graphics windows
//
xdel(winsid()); clear; clc 
ieee(2)  // no more warning division by zero
stacksize(200000000);


//---- REGISTER FITTING FUNCTIONS

function angle = atan2(y, x)

if (x>0) then 
    angle = atan(y/x);    
elseif (y>=0 & x<0) then 
        angle = atan(y/x) + %pi;
elseif (y< 0 & x<0) then 
        angle = atan(y/x) - %pi;
elseif (y> 0 & x==0) then 
        angle = %pi/2;
elseif (y< 0 & x==0) then 
        angle = -%pi/2;
else 
    angle = 0;
end

endfunction

function h = OTF( m, n, off_x, off_y, fc )
    h=zeros(m,n);
    for k=1:m
        for l=1:n
            q = sqrt((k-m/2-1-off_y)^2+(l-n/2-1-off_x)^2);
            if q>fc
                h(k,l)=0;
            else
                b = acos( q/fc );
                h(k,l)=(2*b-sin(2*b))/%pi;
                %h(k,l)=1;
            end
        end
    end
endfunction

function y = fun2fit(x, c)
    y= c(4)+c(1)*sin(2*%pi*x*c(2) - c(3))//*2*%pi*c(2))
    //y= c(4)+c(1)*sin(x*c(2) + c(3))
    // y =  c(1) * sin(c(2)*2*%pi *x+c(3)) + c(4);
endfunction 

function e = myerror(c, x, y, frequency)
    e = fun2fit(x, c, frequency) - y;
endfunction

function y = fun2fit(x, c, frequency)
    //y= c(4)+c(1)*sin(2*%pi*x*frequency + c(3))
    y= c(4)+c(1)*sin(2*%pi*x*c(2) - c(3))//*2*%pi*c(2))
    // y =  c(1) * sin(c(2)*2*%pi *x+c(3)) + c(4);
endfunction 


// Set Camera Parameters PIKE_F_145B_C_fiber.pdf

cpixel_width = 1388;     //  Number of Pixel width
cpixel_height = 1038;     //  Number of Pixel height 
pixelsize =  6.45;       //  Pixel_Size [mu m]
cs_width = cpixel_width * pixelsize;   // width of Camera-Sensor [mm]
cs_height = cpixel_height * pixelsize;   // hight of Camera-Sensor [mm]

//fc = 20; 
fc_est = [0 0 0] //analyitisch bestimmte Ortsfrequenz aus Bildern


// Set Projector Parameters
ps_width = 12;   // width of Camera-Sensor [mm]
ps_height = 7;    // hight of Camera-Sensor [mm]
ppixel_width = 480;     //  Number of Pixel width
ppixel_height = 320;     //  Number of Pix<<l height 
ppixel_size = 7.6e-6;       //  Pixel_Size [mu m]

// Set Lens Parameters
na = .25;   // numerical aperture (1.4 for oil immersion lens)
f_mo = 20;  // Brennweite Mikroskopobjektiv in mm
f_tl = 400; // Brennweite Tubuslinse in m
mag = f_tl/f_mo;     // magnification of Lens (@ TL=400mm)/fMO=20mm
lambda = 550e-9; // wavelength in nanometers
rhop = 0.61*lambda./na; //Zerstreukreis von Objektiv, beugungsbegrenzt
rhop_Sensor = rhop*mag;
rhop_SensorInPixel = round(rhop_Sensor/ppixel_size)

fc_sensor = 1/rhop_Sensor*1e-3     // Grenzfrequenz Objektstrukturen auf Sensor 
fc_objekt = 1/rhop*1e-3;            // Grenzfrequenz Objektstrukturen in Objektebene

fc = 2*na/lambda;

//NEED TO BE REPLACED BY FILTERING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// Set Parameters for Grating
init_theta = 0;    // Initial rotational Angle in Degrees
init_theta = init_theta/180*%pi;

init_phi = 5;       // Initial phase of Illumination pattern/sinosiodal pattern
init_phi = init_phi/180*%pi;

nphi = 3; // number of phases
ntheta = 3; // number of angles
pixelperiode = 5;   // How many Pixels for a periode of a grating

// Generate Vorzeichenmatrix
ind = [(0:(nphi-1)/2) -(nphi-1)/2:1:-1];

//Radius for Zero-Filtering ->Estimate Frequency of Grating
R = 15 ; // radius
// Conversion into rad

theta=linspace(init_theta, init_theta+%pi, ntheta+1);
theta=theta(1:ntheta);
theta_real = theta;
phi=linspace(init_phi, init_phi+2*%pi,nphi+1);
phi=phi(1:nphi);
phi_spectrum = phi;
phi_fit = zeros((ntheta-1)*3+nphi)




// Set Image Parameters
datasets = 1;
//impath = 'C:\Users\Bene\Dropbox\Dokumente\Optical Design\Superresolution\Aufnahmen\Aufnahme_Matlab\';
impath = 'C:\Users\Bene\Dropbox\Dokumente\Optical Design\Superresolution\Aufnahmen\Aufnahme_Nr_4\ROI\';
impath = 'C:\Users\Bene\Dropbox\Dokumente\Optical Design\Superresolution\Aufnahmen\Aufnahme Nr8\ROI\';
impath = 'C:\Users\Bene\Dropbox\Dokumente\Optical Design\Superresolution\Aufnahmen\Aufnahme Nr 10_Sin8\'
impath = 'C:\Users\Bene\Pictures\Picasa\Exports\Fluro_11Neu\'
extension = '.jpg';

file = impath + string(1) + '_' + string(1) + extension;
pixelnumber =  size(imread(file));
n = pixelnumber(1);
m = pixelnumber(2);


// Parameter for Reconstruction of the superresolved image
rc = zeros(n,n,nphi^2);
sp = zeros(n,n,nphi*nphi);
I = zeros(n,n,nphi*nphi);
rimage = zeros(n,n,3);
hv = zeros(n,n);
replc = zeros(n,n);
dr = zeros(n,n);


spectrum = zeros(n,m,datasets*nphi*ntheta);



//PARAMETERS for getting Phase

//Region of Interest for determining Rotation Angle in Fitting the crosssection to a sine
windowsize = 64    ;

//
// Information components separation
//
//Use this to store the values of DIbar, which is the B matrix ( AX = B ) 
DIbars = zeros(n,m,nphi);



// Generate Filter for Zero Frequency Filtering
n2 = n/2;
[x,y] = meshgrid(-n2:n2-1) ;
M = sqrt(x.^2 + y.^2) < R ;
M = bool2s(~M);


IM_temp = zeros(1, ntheta)
RE_temp = zeros(1, ntheta)

for j = 1:ntheta

    f_temp = 0;

    for i = 1:nphi
        file = impath + string(j) + '_' + string(i) + extension;
        Igray=RGB2Gray(imread(file));
       
        //Igray=Igray./max(max(Igray));
        //raw=  Igray;
        //fourier = fftshift(fftw(ifftshift(Igray)));
        fourier = fftshift(fftw(double(Igray)));
        spectrum(:,:,(j-1)*ntheta+i) = fourier;

        X = M.*fourier;     //Filter Zero-Order
        X_abs=(abs(X)./max(abs(X)));
        
        //Show/Save FFT of Image
        norm_spectrum = max(max(abs(fft2(Igray))));
        imshow(abs(fourier)./norm_spectrum * 2.^8);
        imwrite(X_abs*2, impath + 'zero_filtered_' + string(j) + '_' + string(i) + '_fft'  + extension)
        //-------
        
        
        
       

X = X(1: n/2-5, :)
        [b,a]=max(abs(X));  //Get Coordinates of -1. Diff-Order

        //Get "Real-Part" as the Real-part of the Rotation angle
        RE = (n/2-a(2));
        //Get "Imaginary-Part" as the Imag-part of the Rotation angle
        IM = -(n/2-a(1));
        
        RE_temp(j) = RE;
        IM_temp(j) = IM;

        //---------------------------------------GET PHASE ------------------------------------------------------------------------------------------------------

        // Gitterfrequenz aus Fourierspektrum bestimmen

        f_temp = f_temp+  sqrt(RE^2+IM^2);



        rot_angle = atan2(IM, RE);   //equals atan2 already, Calculate Angle of the Rotatation Vector
        //Estimate Phase of Shifted Spectra Shroff et al

        phi_fit((j-1)*3+i) = atan2(imag(fourier(a(2),a(1))), real(fourier(a(2),a(1))))


        if phi_fit((j-1)*3+i) <= 0 then 
            phi_fit((j-1)*3+i)= phi_fit((j-1)*3+i) + %pi;
        end
        disp(   phi_fit((j-1)*3+i)/%pi*180, "Die Phasenverschiebung:" )
        disp( f_temp, "Ortsfrequenz:")
    end       

    fc_est(j) = f_temp/nphi


    //Take Care of the Symmetrie around the Complex plane

    //while(rot_angle<0) rot_angle = rot_angle +2*%pi;
    //end

//while(rot_angle > %pi)
//rot_angle = rot_angle - 2*%pi;
//end
//while(rot_angle < %pi)
//    rot_angle = rot_angle + %pi;
//end

    theta_real(j)=rot_angle;

    disp(rot_angle/%pi*180,' ist ',j,'Der Rot_Angle für Datenset ')
    //disp(IM,' Imaginär:  ',RE,'Realteil: ')


end

//WARUM```?????
//Durch spiegelverkehrte Aufnahme!
//theta_real = theta

theta_temp = theta(2)
theta(2)=theta(3)
theta(3)=theta_temp

theta_temp = theta_real(2)
theta_real(2)=theta_real(3)
theta_real(3)=theta_temp

//theta = theta_real;



//------------??????????????????????????? WARUM ist das nötig?

//---------------------------------------GET PHASE ------------------------------------------------------------------------------------------------------

j=1 //only for one rotation angle = 0°
frequency_est = 0;      // initliaze frequency_estimate

//Get Crosssection X/Y
crosssection = zeros(nphi, n);
x_pos = linspace(1, n, n);
y_pos = x_pos;

    x_min = 1;
    y_min = n/2;

for i=1:ntheta
    file = impath + string(j) + '_' + string(i) + extension;        // Read Image with Rot = 0°, Phase 1,2,3
    I=imread(file);
    Igray=RGB2Gray(I);      
    [n,m] = size(Igray);


    x_max = x_min + windowsize;

    for k = 1:(n)
        //Get Y-Pixel Coordinate for Pixel-Crossection at given Rot-Angle of Sine-Pattern
        y_pos(k) = round(tan(theta(j))*x_pos(k))+y_min;
        //disp(y_pos(k))
        if (y_pos(k) < 1) break;
        elseif (y_pos(k) > n) break;
        end

        //Save Grayvalue at Pixel in crossection-vector
        crosssection(i, k) = Igray(y_pos(k), x_pos(k));
    end 

    // Estimate Frequency of Crosssection Signal

    max_freq = abs(fft(crosssection(i,:)))
    max_freq = max_freq(R:50)
    [u, v] = max(max_freq)
    frequency_est = frequency_est+v+R;
end

//Assume, Frequency of Pattern is always the same!
frequency_est = 1/(frequency_est/(%pi*nphi));       // 2*Pixelperiode equals 1/f



for i=1:nphi
    value_sub = crosssection(i, x_min:x_max);

    // Measured data in vectors x and y
    x= linspace(0, windowsize, windowsize+1)
    x=x';
    y = value_sub'; 

    //y =  c(1).* cos(2*%pi*c(2)*x + c(4)) + c(3);
    // This is our guessed coefficients. First attempt.
    amplitude = (max(value_sub)-min(value_sub))/2
    frequency = frequency_est//.09//2*%pi*0.1   //max(fourier)
    frequency = 1/(fc_est(i)/%pi);
    //frequency = 1/9;
    phase = phi(i)*2*%pi*frequency // phi(i);
    offset = min(value_sub)+amplitude

    c0 = [amplitude frequency phase offset]';
    y0 = fun2fit(x, c0); 





    // Launch the optimization function. Provide data and coefficients
    [f, copt] = leastsq(list(myerror, x, y, frequency), c0) 

    yopt0 = fun2fit(x, copt);

    // Plot the original points and the starting point in our problem
    set("current_figure", i)
    plot(x, y0, 'ro-')
    plot(x, value_sub, 'gx-')
    plot(x, yopt0, 'bo-')
    legend(['Starting point'; 'Original data'; 'After 1st. least square fit']);
    disp(copt, 'Optimierungsparameter')

    while copt(3)<0 
        copt(3) = 2*%pi+copt(3);
    end
    get_phi = copt(3)//(2*%pi*copt(2));
    get_phi = modulo(get_phi, (2*%pi));

    phi_spectrum(1,i) = get_phi
end


//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



//---------------------------------------Tausche Simulierte Daten gegen erfasste Parameter------------------------------------------------------------------------------------------------------


//fc_est = [20 20 20]
//phi = phi_fit
//phi = phi_spectrum;


//The A matrix of the equations( AX = B ), it depends on the shift phases
//phi



A = [1 exp(%i*phi(1)) exp(-%i*phi(1)); 1 exp(%i*phi(2)) exp(-%i*phi(2));1 exp(%i*phi(3)) exp(-%i*phi(3))];
A_inv = inv(A);

A = [1 1 1;1 exp(%i*2*%pi/3) exp(-%i*2*%pi/3);1 exp(%i*4*%pi/3) exp(-%i*4*%pi/3);];
//Inversed A matrix
A_inv = inv(A);

//Solving the Equation-System


D = zeros(n,n,nphi*ntheta);


//solve EQS for one Theta!


for k =1:ntheta     //For Rotationangles

    for i=1:nphi


        for j=1:nphi
            //Summe(Inverse Phasenmatrikx * Spektrum) einer Zeile
            D(:,:,((k-1)*ntheta)+i) = A_inv(i,j)*spectrum(:,:,(((k-1)*ntheta)+j))+ D(:,:,((k-1)*ntheta)+i);

            filename_new = impath+"spectrum_"+string(j)+"_"+string(i)+"_"+string(k)+".png";
            imwrite(abs(spectrum(:,:,(((k-1)*ntheta)+i)))./max(abs(spectrum(:,:,(((k-1)*ntheta)+i)))), filename_new)
        end
        filename_new = impath+"D_"+string(k)+"_"+string(i)+".png";
        imwrite(abs(D(:,:,((k-1)*ntheta)+i))./max(abs(D(:,:,((k-1)*ntheta)+i))), filename_new)
    end


end




//------------ Verschieben der Spektren an richitge Stelle - Entsprechend Gitterfrequenz/Phase





hs = zeros(n, n); // Temporäre matrix für OTF pro Phase 


for j=1:ntheta
    //?????????????????????????????????????????????????????????????????????????????????????????

    //ki = 2*%pi*frequency/n; // Verschieben der Spektren in Pixeleinheiten, anhand der Maxima aus dem Spektrum
    ki = mean(fc_est(j))*2*%pi/n;
    ki = mean(80/(2*%pi))/n;

    //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Kann man sich hier sparen, sind 1:1 die RE/IM aus Max(Spektrum!!!!!!!)


    //?????????????????????????????????????????????????????????????????????????????????????????



    k = [cos(theta(j)) -sin(theta(j)); sin(theta(j)) cos(theta(j)) ]*[ki; 0]; // Rotate Phaseshift according to acquired rotation angle

I_ramp = zeros(n,m); //Reserve RAM for Phase-Ramp
x=1:n;//Phasenrampe in X/Richtung 
y=x';//Phasenrampe in Y/Richtung
x_ramp=repmat(x,n,1); //Phasenrampe in X-Ebene
y_ramp=repmat(y,1,n); //Phasenrampe in Y/Ebene
ind = [(0:(nphi-1)/2) -(nphi-1)/2:1:-1];
xx=repmat(x,n,1); //Grauverlauf??
yy=repmat(y,1,n);

Ir = zeros(n,n);
    kx = k(1); 
    ky = k(2);
    nx = kx*n/(2*%pi); 
    ny = ky*n/(2*%pi);

    for i=1:nphi
        //Erzeugen der OTF
        //hv(:,:) = OTF(n, n,-nx*ind(i),-ny*ind(i), fc);
        Ir(:,:)=exp(-%i*(kx*ind(j)*xx+ky*ind(j)*yy));
        //Verschiebe die Spektren an die richtige Stelle
        I_ramp(:,:)=exp(-%i*(kx*ind(i)*x_ramp+ky*ind(i)*y_ramp));               // Phasenrampe um Phasenshift zu erzeugen! Entsprechend der Phasenverschiebung der Sinusmuster

        filename_new = impath+"PhasenRampe_"+string(j)+"_"+string(i)+".png"
        imwrite(atan(imag(I_ramp(:,:)),real(I_ramp(:,:)))./max(atan(imag(I_ramp(:,:)),real(I_ramp(:,:)))), filename_new)


        filename_new = impath+"before_shifted_"+string(j)+"_"+string(i)+".png"
        imwrite(10*abs(D(:,:,((j-1)*ntheta)+i))./max(abs(D(:,:,((j-1)*ntheta)+i))), filename_new)

        D(:,:,((j-1)*ntheta)+i) = (fftw(ifftshift(ifft(D(:,:,((j-1)*ntheta)+i)).*I_ramp)));
        imshow(10*abs(D(:,:,((j-1)*ntheta)+i))./max(abs(D(:,:,((j-1)*ntheta)+i))));
        filename_new = impath+"shifted_"+string(j)+"_"+string(i)+".png"
        imwrite(10*abs(D(:,:,((j-1)*ntheta)+i))./max(abs(D(:,:,((j-1)*ntheta)+i))), filename_new)
        
        //D(:,:,((j-1)*ntheta)+i) =D(:,:,((j-1)*ntheta)+i).*conj(hv(:,:)); //*Scalefactor??????
        
        hs = hs + abs(hv(:,:)).^2;
    end    
end

dr = zeros(n,n);

for t = 1:nphi*ntheta
dr = dr+D(:,:,t)./ ( hs + .005*length(nphi)*(.0000001)^2);
end
fimage = ifftshift(ifft(dr));


imshow(abs(fimage)./max(abs(fimage)))
w=getdate()
//mprintf("hour:%d,minute:%d,second:%d",w(7),w(8),w(9))
filename_new = impath+"renconstructed_"+string(w(7))+"_"+string(w(8))+"_"+string(w(9))+".png"
imwrite(abs(fimage)./max(abs(fimage)), filename_new)



//Experimental
// Deconvolution and reconstruction with a Wiener

for t = 1:nphi*ntheta
dr = dr+rc(:,:,t)./ ( hs + .005*length(ntheta)*(.0000001)^2);
end
figure; colormap('gray'); imagesc(abs(dr(:,:))); axis('square'); title('Reconstruction of The Object in Reciprocal Space')
%
%Triangular function
%
[k_x, k_y]=meshgrid(-n/2+1:n/2, -n/2+1:n/2);
k_r = sqrt(k_x.^2+k_y.^2);
k_max = .9*.9*fc*((nphases-1)/2+1);
bhs = cos(pi*k_r/(2*k_max));
indi = find( k_r > k_max );
bhs(indi) = 0;
figure; colormap('gray'); imagesc(abs(dr.*bhs)); axis('square'); title('Apodization of The Object in Reciprocal Space')
drr = dr.*bhs;
fimage = ifft2(ifftshift(drr));

figure; imshow(abs(fimage)); 

