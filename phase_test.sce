xdel(winsid()); clear; clc 
//funcprot(0);

function y = fun2fit(x, c)
    y= c(4)+c(1)*sin(2*%pi*x*c(2) + c(3))
 // y =  c(1) * sin(c(2)*2*%pi *x+c(3)) + c(4);
endfunction 

function e = myerror(c, x, y)
  e = fun2fit(x, c) - y;
endfunction 

//Region of Intersion
windowsize = 15;



impath = 'C:\Users\Bene\Dropbox\Dokumente\Optical Design\Superresolution\Aufnahmen\Aufnahme_Nr_2\ROI\';
extension = '.png';
theta_real = [0.    1.0471976    2.0943951 ]

j=1

for i=1:3


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
phi = copt(3);
phi = modulo(phi, (%pi));
phi_deg = (phi/%pi)*180
disp(phi, phi_deg,  'Phi;')
end
