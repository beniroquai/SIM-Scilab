function y=yth(t, x)
   y  = 27*sin(x(1)*t+x(2))+90
endfunction



// we have the m measures (ti, yi):
m = 41;
tm = [175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215]';
ym = [77,78,78,79,88,96,99,100,99,92,84,78,78,79,82,90,95,100,98,90,84,79,79,83,88,98,100,103,102,98,90,83,80,81,88,95,99,103,101,95,91]';
tm = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5]';
ym = [0.79, 0.59, 0.47, 0.36, 0.29, 0.23, 0.17, 0.15, 0.12, 0.08]';
// measure weights (here all equal to 1...)
wm = ones(m,1);

// and we want to find the parameters x such that the model fits the given
// data in the least square sense:
//
//  minimize  f(x) = sum_i  wm(i)^2 ( yth(tm(i),x) - ym(i) )^2

// initial parameters guess
x0 = [1.5 ; 0.8];

// in the first examples, we define the function fun and dfun
// in scilab language
function e=myfun(x, tm, ym, wm)
   e = wm.*( yth(tm, x) - ym )
endfunction

function g=mydfun(x, tm, ym, wm)
   v = wm.*sin(-x(2)*tm)
   g = [v , -x(1)*tm.*v]
endfunction

// now we could call leastsq:

// 1- the simplest call
[f,xopt, gopt] = leastsq(list(myfun,tm,ym,wm),x0)

// 2- we provide the Jacobian
[f,xopt, gopt] = leastsq(list(myfun,tm,ym,wm),mydfun,x0)

// a small graphic (before showing other calling features)
tt = linspace(0,1.1*max(tm),100)';
yy = yth(tt, xopt);
scf();
plot(tm, ym, "kx")
plot(tt, yy, "b-")
legend(["measure points", "fitted curve"]);
xtitle("a simple fit with leastsq")

// 3- how to get some information (we use imp=1)
[f,xopt, gopt] = leastsq(1,list(myfun,tm,ym,wm),mydfun,x0)

// 4- using the conjugate gradient (instead of quasi Newton)
[f,xopt, gopt] = leastsq(1,list(myfun,tm,ym,wm),mydfun,x0,"gc")

// 5- how to provide bound constraints (not useful here !)
xinf = [-%inf,-%inf];
xsup = [%inf, %inf];
// without Jacobian:
[f,xopt, gopt] = leastsq(list(myfun,tm,ym,wm),"b",xinf,xsup,x0)
// with Jacobian :
[f,xopt, gopt] = leastsq(list(myfun,tm,ym,wm),mydfun,"b",xinf,xsup,x0)

// 6- playing with some stopping parameters of the algorithm
//    (allows only 40 function calls, 8 iterations and set epsg=0.01, epsf=0.1)
[f,xopt, gopt] = leastsq(1,list(myfun,tm,ym,wm),mydfun,x0,"ar",40,8,0.01,0.1)
