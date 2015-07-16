// Delete figures, memory and screen
xdel(winsid()); clear; clc 
function y = fun2fit(x, c)
  y =  c(1) * sin(c(2) *x+c(3)) + c(4);
endfunction 

// This is how we define the error to be minimized.
// Function leastsq iterates through this error until it's
// not posible to reduce it any more.
function e = myerror(c, x, y)
  e = fun2fit(x, c) - y;
endfunction 

// Declare the objective function, or function to be minimized
//exec('C:\Users\Bene\Dropbox\Dokumente\Optical Design\Superresolution\Scilab\nonlinear_fit_1.sci') 

// Define our initial data or measured points
//x = [175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215]';
x=linspace(0,11, 12)
x=x'
y = [77,78,78,79,88,96,99,100,99,92,84,78]';  

// This is our guessed coefficients. First attempt.
c0 = [27 1 3 90]';
y0 = fun2fit(x, c0); 

// Plot the original points and the starting point in our problem
set("current_figure", 0)
plot(x, y, 'g*', x, y0, 'ro') 

for i=1:100
// Launch the optimization function. Provide data and coefficients
[f, copt] = leastsq(list(myerror, x, y), c0) 
c0=copt
end

// Plot results after optimization
yopt0 = fun2fit(x, copt);
plot(x, yopt0, 'bo-')
legend(['Original data'; 'Starting point'; 'After 1st. least square fit']);



