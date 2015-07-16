//1. Reconstruction.sci
//
// Clear previous variables and graphics windows
//
xdel(winsid()); clear; clc 

stacksize(100000000);


//---- REGISTER FITTING FUNCTIONS

function y = fun2fit(x, c)
    y= c(4)+c(1)*sin(2*%pi*x*c(2) + c(3))
    // y =  c(1) * sin(c(2)*2*%pi *x+c(3)) + c(4);
endfunction 

function e = myerror(c, x, y)
    e = fun2fit(x, c) - y;
endfunction 



// Set Camera Parameters PIKE_F_145B_C_fiber.pdf

cpixel_width = 1388;     //  Number of Pixel width
cpixel_height = 1038;     //  Number of Pixel height 
pixelsize =  6.45;       //  Pixel_Size [mu m]
cs_width = cpixel_width * pixelsize;   // width of Camera-Sensor [mm]
cs_height = cpixel_height * pixelsize;   // hight of Camera-Sensor [mm]

fc = 20; 
fc_est = [0 0 0] //analyitisch bestimmte Ortsfrequenz aus Bildern


// Set Projector Parameters
ps_width = 12;   // width of Camera-Sensor [mm]
ps_height = 7;    // hight of Camera-Sensor [mm]
ppixel_width = 480;     //  Number of Pixel width
ppixel_height = 320;     //  Number of Pix<<l height 
ppixel_size = 7,6;       //  Pixel_Size [mu m]

// Set Lens Parameters
na = .25;   // numerical aperture (1.4 for oil immersion lens)
mag = 10;     // magnification of Lens (@ TL=200mm)
lambda = 1000; // wavelength in nanometers
rho = 0,61*lambda/na;



//NEED TO BE REPLACED BY FILTERING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// Set Parameters for Grating
init_theta = 0;    // Initial rotational Angle in Degrees
init_theta = init_theta/180*%pi;

init_phi = 0;       // Initial phase of Illumination pattern/sinosiodal pattern
init_phi = init_phi/180*%pi;

nphi = 3; // number of phases
ntheta = 3; // number of angles
pixelperiode = 5;   // How many Pixels for a periode of a grating


// Conversion into rad

theta=linspace(init_theta, init_theta+%pi, ntheta+1);

theta=theta(1:ntheta);
theta_real = theta;
phi=linspace(init_phi, init_phi+2*%pi,nphi+1);
phi=phi(1:nphi);
phi_real = phi;
phi_est = zeros((ntheta-1)*3+nphi)




// Set Image Parameters
datasets = 1;
//impath = 'C:\Users\Bene\Dropbox\Dokumente\Optical Design\Superresolution\Aufnahmen\Aufnahme_Matlab\';
impath = 'C:\Users\Bene\Dropbox\Dokumente\Optical Design\Superresolution\Aufnahmen\Aufnahme_Nr_4\ROI\';
impath = 'C:\Users\Bene\Pictures\Picasa\Exports\Fluro_11Neu\'
extension = '.jpg';


// LOAD IMAGES 
//I(:,:,(itheta-1)*3+iphi) = (1-cos(kx*X1+ky*Y1+phi(iphi)))/2;
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



raw = zeros(n,m,datasets*nphi*ntheta);
spectrum = zeros(n,m,datasets*nphi*ntheta);


ki = 2*%pi*(fc)/n; //???????????????


//PARAMETERS for getting Phase

//Region of Interest for determining Rotation Angle
windowsize = 15;

//
// Information components separation
//
//Use this to store the values of DIbar, which is the B matrix ( AX = B ) 
DIbars = zeros(n,m,nphi);


R = 30 ; // radius
//n2 = floor(n/2) ;
n2 = n/2;
[x,y] = meshgrid(-n2:n2-1) ;
M = sqrt(x.^2 + y.^2) < R ;
M = bool2s(~M);

for j = 1:ntheta
    
    fc_det = 0;
    
    for i = 1:nphi
        file = impath + string(j) + '_' + string(i) + extension;
        I=imread(file);
        Igray=RGB2Gray((I));
        //Igray=Igray./max(max(Igray));
        raw(:,:,i)=  Igray;
        fourier = fftshift(fftw(ifftshift(raw(:,:,i))));
        spectrum(:,:,(j-1)*ntheta+i) = fourier;

        X = M.*fourier;     //Filter Zero-Order
        X_abs=(abs(X)./max(abs(X)));




        //Show/Save FFT of Image
        norm_spectrum = max(max(abs(fftw(raw(:,:,1)))));
        imshow(abs(fourier)./norm_spectrum * 2.^8);
        imwrite(abs((fourier)./norm_spectrum * 2.^8),impath + 'fft_'  + string(j) + '_' + string(i)  + extension)
        imwrite(abs((atan(imag(fourier),real(fourier)))./max(atan(imag(fourier),real(fourier)))),impath + 'angle_' + string(j) + '_' + string(i) + '_fft'  + extension)
        imwrite(X_abs*255, impath + 'zero_filtered_' + string(j) + '_' + string(i) + '_fft'  + extension)
        //-------

        [b,a]=max(abs(X));  //Get Coordinates of -1. Diff-Order
        [sx,sy]=size(abs(X));    //Get Size of Spectrum

        //Get "Real-Part" as the Real-part of the Rotation angle
        RE = (sx/2-a(2));
        //Get "Imaginary-Part" as the Imag-part of the Rotation angle
        IM = -(sy/2-a(1));
        
        // Gitterfrequenz aus Fourierspektrum bestimmen
        
        fc_det = fc_det+  sqrt(RE^2+IM^2);
        
        

        rot_angle = atan(IM, RE);   //equals atan2 already, Calculate Angle of the Rotatation Vector
        //Estimate Phase of Shifted Spectra

        phi_est((j-1)*3+i) = atan(imag(fourier(a(2),a(1))), real(fourier(a(2),a(1))))
        
        if phi_est <= 0 then 
            phi_est = phi_est + %pi;
        end
        disp(   phi_est((j-1)*3+i)/%pi*180, "Die Phasenverschiebung:" )
        disp( fc_det, "Ortsfrequenz:")
    end       
   
fc_est(j) = fc_det/nphi


    //Take Care of the Symmetrie around the Complex plane
    //if(RE<0 & IM >=0) then rot_angle = rot_angle  + %pi;
    //elseif(RE<0 & IM <0) then rot_angle = rot_angle  - (%pi);
    //else rot_angle = rot_angle;
    //end 


    //rot_angle=abs(rot_angle-%pi)
    //if rot_angle > %pi then rot_angle = rot_angle - %pi
    //end

    //if(rot_angle>= %pi) then rot_angle = rot_angle - %pi;
    //end
    //if(rot_angle<= -%pi) then rot_angle = rot_angle + %pi;
    //end
    theta_real(j)=rot_angle;

    disp(rot_angle/%pi*180,' ist ',j,'Der Rot_Angle für Datenset ')
    //disp(IM,' Imaginär:  ',RE,'Realteil: ')


end



//---------------------------------------GET PHASE 

j=1

for i=1:ntheta


    file = impath + string(j) + '_' + string(i) + extension;
    I=imread(file);
    Igray=RGB2Gray(I);      
    [n,m] = size(Igray);
    x_min = n/2;
    x_max = n/2+windowsize;


        //Get Crosssection X/Y
        value_positive = linspace(0, n/2, n/2+1);
        x_positive = linspace(0,n/2,n/2+1);
        y_positive = linspace(0,m/2,m/2+1);
        
        value_negative = linspace(0, n/2, n/2+1);
        x_negative = linspace(0,n/2,n/2+1);
        y_negative = linspace(0,m/2,m/2+1);
        
        value = linspace(0,m,m+1)
        
        for i = 1:(max(x_positive))
            y_positive(i) = round(tan(theta_real(j))*x_positive(i));
                if (y_positive(i) >= m/2) break;
            end
            value_positive(i) = Igray(y_positive(i)+m/2, x_positive(i)+n/2);
            value_negative(i) = Igray(-(y_positive(i))+m/2, -x_positive(i)+n/2);
        end 
    //inverse negative row
    value_negative=value_negative(:,$:-1:1);
    //stitch positive/negative values together
    value=cat(2, value_negative, value_positive);


    x_sub=linspace(x_min, x_max, windowsize+1);
    value_sub = value(x_min:x_max);

    // Measured data in vectors x and y
    x = x_sub';
    x= linspace(0, windowsize, windowsize+1)
    x=x';
    y = value_sub'; 


    //y =  c(1).* cos(2*%pi*c(2)*x + c(4)) + c(3);
    // This is our guessed coefficients. First attempt.
    amplitude = (max(value_sub)-min(value_sub))/2
    frequency = .09//2*%pi*0.1   //max(fourier)
    phase = 0
    offset = min(value_sub)+amplitude

    c0 = [amplitude frequency phase offset]';
    y0 = fun2fit(x, c0); 

    // Plot the original points and the starting point in our problem
    set("current_figure", 0)
    plot(x, y0, 'ro') 
    plot(x, value_sub, 'gx')


    // Launch the optimization function. Provide data and coefficients
    [f, copt] = leastsq(list(myerror, x, y), c0) 
    yopt0 = fun2fit(x, copt);
    plot(x, yopt0, 'bo-')
    legend(['Original data'; 'Starting point'; 'After 1st. least square fit']);
    disp(copt, 'Optimierungsparameter')

    while copt(3)<0 
        copt(3) = %pi+copt(3);
    end
    get_phi = copt(3);
    get_phi = modulo(get_phi, (%pi));
    get_phi_deg = (get_phi/%pi)*180
    disp(get_phi, get_phi_deg,  'Phi;');
    phi_real(1,i) = get_phi
end



//---------------------------------------Tausche Simulierte Daten gegen erfasste Parameter



phi = phi_est
theta = theta_real;
//phi = phi_real;


//The A matrix of the equations( AX = B ), it depends on the shift phases
//phi



A = [1 exp(%i*phi(1)) exp(-%i*phi(1)); 1 exp(%i*phi(2)) exp(-%i*phi(2));1 exp(%i*phi(3)) exp(-%i*phi(3))];
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

//GET Phaseshift
//j=1
value_positive = linspace(0, n/2, n/2+1);
x_positive = linspace(0,n/2,n/2+1);
y_positive = linspace(0,m/2,m/2+1);

//D=H*G, OTF times patterned illuminated object
I1_bar = fft(ifft(D(:,:,2))*exp(-%i*2*%pi*(phi(2))));


//------------ Solve

//
//D_0 = A_inv(1,1)*l(1)+ A_inv(1,2)*l(2)+ A_inv(1,3)*l(3);
//D_p = A_inv(1,1)*l(1)+ A_inv(1,2)*l(2)+ A_inv(1,3)*l(3);
//D_n = A_inv(1,1)*l(1)+ A_inv(1,2)*l(2)+ A_inv(1,3)*l(3);
//D_0 = D_0+D_p+D_n;  
//imshow(abs(D(:,:,1))./max(abs(D(:,:,1))))
//imshow(abs(D(:,:,2))./max(abs(D(:,:,2))))
//imshow(abs(D(:,:,3))./max(abs(D(:,:,3))))
//D_0=D(:,:,1)+D(:,:,1)+D(:,:,1);
//D_x = fftshift(ifft(ifftshift(D_0)));
//imshow(abs(D_x)./max(abs(D_x)));



//------------ Verschieben der Spektren an richitge Stelle - Entsprechend Gitterfrequenz/Phase

for j=1:ntheta
//?????????????????????????????????????????????????????????????????????????????????????????

ki = 2*%pi*(fc_est(i)/2)/n; //?????????????

//?????????????????????????????????????????????????????????????????????????????????????????

ind = [(0:(nphi-1)/2) -(nphi-1)/2:1:-1]; //create Vorzeichen for phaseshift

I_ramp = zeros(n,m); //Reserve RAM for Phase-Ramp

x=1:n;//Phasenrampe in X/Richtung 
y=x';//Phasenrampe in Y/Richtung
x_ramp=repmat(x,n,1); //Phasenrampe in X-Ebene
y_ramp=repmat(y,1,n); //Phasenrampe in Y/Ebene



    k = [cos(theta(j)) -sin(theta(j)); sin(theta(j)) cos(theta(j)) ]*[ki; 0]; // Rotate Phaseshift according to acquired rotation angle

    kx = k(1); 
    ky = k(2);
    nx = kx*n/(2*%pi); 
    ny = ky*n/(2*%pi);


    for i=1:nphi

        I_ramp(:,:)=exp(-%i*(kx*ind(i)*x_ramp+ky*ind(i)*y_ramp)); // Phasenrampe um Phasenshift zu erzeugen! Entsprechend der Phasenverschiebung der Sinusmuster
        //imshow((atan(imag(I_ramp),real(I_ramp)))./max(atan(imag(I_ramp),real(I_ramp))))
        filename_new = impath+"before_shifted_"+string(j)+"_"+string(i)+".png"
        imwrite(abs(D(:,:,((j-1)*ntheta)+i))./max(abs(D(:,:,((j-1)*ntheta)+i))), filename_new)
        
        D(:,:,((j-1)*ntheta)+i) = fft(ifftshift(ifft(D(:,:,((j-1)*ntheta)+i)).*I_ramp));
        imshow(abs(D(:,:,((j-1)*ntheta)+i))./max(abs(D(:,:,((j-1)*ntheta)+i))));
        filename_new = impath+"shifted_"+string(j)+"_"+string(i)+".png"
        imwrite(abs(D(:,:,((j-1)*ntheta)+i))./max(abs(D(:,:,((j-1)*ntheta)+i))), filename_new)
    end    
end


//D=H*G, OTF times patterned illuminated object


I1_bar = fft(ifft(D(:,:,nphi*1))*exp(-%i*2*%pi*(phi(1))));
I2_bar = fft(ifft(D(:,:,nphi*2))*exp(-%i*2*%pi*(phi(2))));
I3_bar = fft(ifft(D(:,:,nphi*3))*exp(-%i*2*%pi*(phi(3))));
I1_barbar = fft(ifft(D(:,:,nphi*1))*exp(%i*2*%pi*(2*phi(1))));
I2_barbar = fft(ifft(D(:,:,nphi*2))*exp(%i*2*%pi*(2*phi(2))));
I3_barbar = fft(ifft(D(:,:,nphi*3))*exp(%i*2*%pi*(2*phi(3))));


// Combine---------------
Ic2 = I1_bar + I1_barbar;
Ic4 = I2_bar + I2_barbar;
Ic6 = I3_bar + I3_barbar;

Ic1 = D(:,:,ntheta*1);
Ic2 = D(:,:,ntheta*2);
Ic3 = D(:,:,ntheta*3);


I_rec1 = Ic1+Ic2;
I_rec2 = Ic2+Ic4;
I_rec3 = Ic3+Ic6;

I_rec1 = ifft(fftshift(I_rec1));
I_rec2 = ifft(fftshift(I_rec2));
I_rec3 = ifft(fftshift(I_rec3));


//figure; imshow(abs(I_rec1)./max(abs(I_rec1)));  title('Reconstruction of The Object in Reciprocal Space')
//figure; imshow(abs(I_rec2)./max(abs(I_rec1)));  title('Reconstruction of The Object in Reciprocal Space')
//figure;  imshow(abs(I_rec3)./max(abs(I_rec1))); title('Reconstruction of The Object in Reciprocal Space')
//
I_rec = I_rec1 + I_rec2 + I_rec3
imshow(abs(I_rec)./max(abs(I_rec)))
w=getdate()
//mprintf("hour:%d,minute:%d,second:%d",w(7),w(8),w(9))
filename_new = impath+"renconstructed_"+string(w(7))+"_"+string(w(8))+"_"+string(w(9))+".png"
imwrite(abs(I_rec)./max(abs(I_rec)), filename_new)


for j=1:nphi
    itheta = j
    for k = 1:nphi
        temp_separated(:,:,k) = A_inv(j,k).*spectrum(:,:,k);
        sp(:,:,(itheta-1)*3+j) = sp(:,:,(itheta-1)*3+j)+temp_separated(:,:,k);
    end
    norm_sp = max(abs(sp(:,:,(itheta-1)*3+j)));
    //imshow(abs(sp(:,:,(itheta-1)*3+j))./norm_sp*2.^8);

end

