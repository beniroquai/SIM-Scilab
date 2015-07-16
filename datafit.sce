//generate the data
function y=FF(x, p),y=p(1)*(x-p(2))+p(3)*x.*x,endfunction
X=[];Y=[];
pg=[34;12;14] //parameter used to generate data
for x=0:.1:3, Y=[Y,FF(x,pg)+100*(rand()-.5)];X=[X,x];end
Z=[Y;X];

//The criterion function
function e=G(p, z),
  y=z(1),x=z(2);
  e=y-FF(x,p),
endfunction

//Solve the problem
p0=[3;5;10]    
[p,err]=datafit(G,Z,p0);

scf(0);clf()
plot2d(X,FF(X,pg),5) //the curve without noise
plot2d(X,Y,-1)  // the noisy data
plot2d(X,FF(X,p),12) //the solution

//the gradient of the criterion function
function s=DG(p, z),
  a=p(1),b=p(2),c=p(3),y=z(1),x=z(2),
  s=-[x-b,-a,x*x]
endfunction

[p,err]=datafit(G,DG,Z,p0);
scf(1);clf()
plot2d(X,FF(X,pg),5) //the curve without noise
plot2d(X,Y,-1)  // the noisy data
plot2d(X,FF(X,p),12) //the solution

// Add some bounds on the estimate of the parameters
// We want positive estimation (the result will not change)
[p,err]=datafit(G,DG,Z,'b',[0;0;0],[%inf;%inf;%inf],p0,algo='gc');
scf(1);clf()
plot2d(X,FF(X,pg),5) //the curve without noise
plot2d(X,Y,-1)  // the noisy data
plot2d(X,FF(X,p),12) //the solution
