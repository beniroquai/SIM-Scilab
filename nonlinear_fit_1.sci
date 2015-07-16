// This is the function to be approximated.
// The error function calls this function iteratively
function y = fun2fit(x, c)
  y =  c(1)*x .* cos(c(2)*x) + c(3);
endfunction 

// This is how we define the error to be minimized.
// Function leastsq iterates through this error until it's
// not posible to reduce it any more.
function e = myerror(c, x, y)
  e = fun2fit(x, c) - y;
endfunction 
