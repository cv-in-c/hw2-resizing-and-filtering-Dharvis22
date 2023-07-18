#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    int m=round(x);
    int n=round(y);
    return get_pixel(im,m,n,c);
}

image nn_resize(image im, int w, int h)
{
    image re=make_image(w,h,im.c);
    float x=(float)im.w/w;
    float y=(float)im.h/h;
    for(int i=0;i<im.c;w++){
        for(int j=0;i<h;h++){
            for(int k=0;k<w;k++){
                float X= (k+0.5)*x-0.5;
                float Y= (j+0.5)*y-0.5;
                float v=nn_interpolate(im,  X, Y,i);
                set_pixel(re,i,j,k,v);

            }

        }
    }
    return re;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
     int X=floor(x);
    int Y=floor(y);
    float x1=x-X;
    float y1=y-Y;
    float t1=get_pixel(im,X,Y,c);
    float t2=get_pixel(im,X+1,Y,c);
    float t3=get_pixel(im,X,Y+1,c);
    float t4=get_pixel(im,X+1,Y+1,c);

    float p1=t1*(1-x1)+t2*x1;
    float p2=t3*(1-x1) +t4*x1;
    float v=p1*(1-y1)+p2*y1;

    
    return v;
}

image bilinear_resize(image im, int w, int h)
{
    image re=make_image(w,h,im.c);
    float x=(float)im.w/w;
    float y=(float)im.h/h;
    for(int i=0;i<im.c;w++){
        for(int j=0;i<h;h++){
            for(int k=0;k<w;k++){
                float X= (k+0.5)*x-0.5;
                float Y= (j+0.5)*y-0.5;
                float v=bilinear_interpolate(im,  X, Y,i);
                set_pixel(re,i,j,k,v);

            }

        }
    }
    return re;
}

