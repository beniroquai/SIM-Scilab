function [x,y]=rotation(x0,y0,angle)
    angle=%pi*angle/180//convertion degré<->radian
    M=[cos(angle) -sin(angle); sin(angle) cos(angle)]//rotation matrix
    vec=M*[x0;y0]//rotation
    x=vec(1);y=vec(2)
endfunction

function im=imrotation(im0,angle)
    im0=double(im0)
    // angle reduction to [090[
    angle=pmodulo(angle,360)
    while angle>90
        angle=angle-90
        im0=im(:,$:-1:1)'//rotation 90 degree
    end
    //rotation 0<angle<90
    if angle<>0 then 
        [nx,ny]=size(im0)
        [x,y]=rotation(0,ny,angle),tx=fix(-x)+1//x  shift
        [x,y]=rotation(nx,0,angle),txinf=fix(nx-x+1)//xinf shift
        [x,y]=rotation(nx,ny,angle),ty=fix(y-ny)+1//y shift
        im=zeros(nx+txinf-tx,ny+ty)
        for x=1:nx+tx
            for y=1:ny+ty
                [x0,y0]=rotation(x-tx,y,-angle)//after rotation of 'angle' pixel (x0,y0) ->(x,y)
                x1=fix(x0),x2=x1+1,y1=fix(y0),y2=y1+1// (x0,y0) is in the square [x1,x2]*[y1,y2]
                px=x2-x0,py=y2-y0//px% of x1, (1-px)% of x2,py% of y1, (1-py)% of y2 
                if (x1>0)&(x2<nx+1)&(y1>0)&(y2<ny+1) then 
                    im(x,y)=px*py*im0(x1,y1)+px*(1-py)*im0(x1,y2)+(1-px)*py*im0(x2,y1)+(1-px)*(1-py)*im0(x2,y2)
                end
            end
        end
        
    end
    im=uint8(im)
endfunction

function im=imrotate(im0,angle)
    num_type=type(im0)
    im=im0;//in case of  errors
    if type(angle)==1 then 
        if (num_type==8)|(num_type==1) then //grayscale image 
            im=imrotation(im0,angle)
        elseif num_type==17 then //RGB image
            R=im0(:,:,1),R=imrotation(R,angle)//red
            G=im0(:,:,2),G=imrotation(G,angle)//green
            B=im0(:,:,3),B=imrotation(B,angle)//blue
            [p,n]=size(R),im=hypermat([p,n,3])//initialize result
            im(:,:,1)=R, im(:,:,2)=G, im(:,:,3)=B//final image
        else printf('wrong type (%d) for first argument ''im0''\n',numtype)
        end
    else  printf('wrong type (%d) for second argument ''angle''\n',type(angle))
    end
    im=uint8(im)
endfunction


// replace image.png by your test image and uncomment code below :
//[im0]=imread('image.png');
//im0=uint8(im0);// convert image to 8bits unsigned 
//im0=imresize(im0,0.5);// smaller test image to avoid optimize cpu time!
//im=imrotation(im0,30);
//imshow(im);   

 
