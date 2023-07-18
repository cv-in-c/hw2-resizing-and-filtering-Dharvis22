#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    float sum=0;
    for(int i=0;i<im.c;i++){
        for(int j=0;j<im.h;j++){
            for(int k=0;k<im.w;k++){
                sum+=get_pixel(im,k,j,i);

            }
        }
    }
    for(int i=0;i<im.c;i++){
        for(int j=0;j<im.h;j++){
            for(int k=0;k<im.w;k++){
                set_pixel(im,k,j,i,get_pixel(im,k,j,i)/sum);
            }
        }
    }
}

image make_box_filter(int w)
{
     image f=make_image(w,w,1);
    float v=1.0/(w*w);
    for(int i=0;i<f.h;i++){
        for(int j=0;j<f.w;j++){
            set_pixel(f,j,i,0,v);
        }
    }

    return f;
}

image convolve_image(image im, image filter, int preserve)
{
    assert(im.c==filter.c||filter.c==1);
   int preserve_channels=(preserve&&im.c==filter.c);
   image convolved =make_image(im.w,im.h,preserve_channels ?im.c:1);
   int half_width=filter.w/2;
   int half_height=filter.h/2;
   for(int i=0;i<convolved.c;i++){
    for(int j=0;j<convolved.h;j++){
        for(int k=0;k<convolved.w;k++){
            float s=0;
            for(int fy=0;fy<filter.h;++fy){
                for(int fx=0;fx<filter.w;fx++){
                    int X=k+fx-half_width;
                    int Y=j+fy- half_height;
                    float filter_v=get_pixel(filter,fx,fy,0);
                    float image_v=get_pixel(im,X,Y,i);
                    s+=filter_v*image_v;


                }
            }
            if(preserve_channels){
                set_pixel(convolved,k,j,i,s);

            }
            else{
                set_pixel(convolved,k,j,0,s);
            }
        }
    }
   }

    return convolved;
}
   


image make_highpass_filter()
{
    image filter=make_image(3,3,1);
    set_pixel(filter ,0,0,0,0);
    set_pixel(filter,1,0,0,-1);
    set_pixel(filter,2,0,0,0);
    set_pixel(filter,0,1,0,-1);
    set_pixel(filter,1,1,0,4);
    set_pixel(filter,2,1,0,-1);
    set_pixel(filter,0,2,0,0);
    set_pixel(filter,1,2,0,-1);
    set_pixel(filter,2,2,0,0);
    
    return filter;
    
}

image make_sharpen_filter()
{
    image filter=make_image(3,3,1);
    set_pixel(filter ,0,0,0,0);
    set_pixel(filter,1,0,0,-1);
    set_pixel(filter,2,0,0,0);
    set_pixel(filter,0,1,0,-1);
    set_pixel(filter,1,1,0,5);
    set_pixel(filter,2,1,0,-1);
    set_pixel(filter,0,2,0,0);
    set_pixel(filter,1,2,0,-1);
    set_pixel(filter,2,2,0,0);
    
    return filter;
}
    

image make_emboss_filter()
{
   
     image filter=make_image(3,3,1);
    set_pixel(filter ,0,0,0,-2);
    set_pixel(filter,1,0,0,-1);
    set_pixel(filter,2,0,0,0);
    set_pixel(filter,0,1,0,-1);
    set_pixel(filter,1,1,0,1);
    set_pixel(filter,2,1,0,1);
    set_pixel(filter,0,2,0,0);
    set_pixel(filter,1,2,0,1);
    set_pixel(filter,2,2,0,2);
    
    return filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: Filters like highpass,sharpen etc. should be preserved so that the output image maintains the color information. 

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Post-processing techniques like unsharp masking or edge enhancement can be applied to improve the overall appearance.


image make_gaussian_filter(float sigma)
{
    int size=(int)(6*sigma);
    if(size%2==0)
        size+= 1;

    image filter=make_image(size,size,1);

    int center=size/2;
    float sum=0;

    for (int i=0;i < filter.h; ++i) {
        for (int j=0; j<filter.w; ++j) {
            float x=j-center;
            float y =i-center;
            float v =expf(-(x * x + y * y)/ 2 * sigma * sigma));
            set_pixel(filter, j, i, 0, v);
            sum+= v;
              }
    }

    for(int i = 0; i<filter.h; ++i) {
        for(int j = 0; j<filter.w; ++j) {
            float v = get_pixel(filter,j,i,0)/sum;
            set_pixel(filter,j,i,0,v);
        }
    }

    return filter;
}

image add_image(image a, image b)
{
    assert(a.w==b.w&&a.h==b.h&&a.c==b.c);

    image result =make_image(a.w,a.h,a.c);

    for (int i = 0; i < result.c; i++) {
        for (int j = 0; j < result.h; j++) {
            for (int k = 0; k < result.w;k++) {
                float va = get_pixel(a,k,j,i);
                float vb = get_pixel(b,k,j,i);
                set_pixel(result,k,j,i,va+vb);
            }
        }
    }

    return result;
}

image sub_image(image a, image b)
{
    assert(a.w==b.w&&a.h==b.h&&a.c==b.c);

    image result =make_image(a.w,a.h,a.c);

    for (int i=0; i<result.c;i++) {
        for (int j = 0; j < result.h;j++) {
            for (int k= 0; k< result.w;k++) {
                float va = get_pixel(a,i,j,k);
                float vb=get_pixel(b,i,j,k);
                set_pixel(result,k,j,i,va-vb);
            }
        }
    }
    return result;
}

image make_gx_filter()
{
     image filter=make_image(3,3,1);

    set_pixel(filter,0,0,0,-1);
    set_pixel(filter,1,0,0,0);
    set_pixel(filter,2,0,0,1);

    set_pixel(filter,0,1,0,-2);
    set_pixel(filter,1,1,0,0);
    set_pixel(filter,2,1,0,2);

    set_pixel(filter,0,2,0,-1);
    set_pixel(filter,1,2,0,0);
    set_pixel(filter,2,2,0,1);

    return filter;
}

image make_gy_filter()
{
   mage filter=make_image(3,3,1);

    set_pixel(filter,0,0,0,-1);
    set_pixel(filter,1,0,0,-2);
    set_pixel(filter,2,0,0,-1);

    set_pixel(filter,0,1,0,0);
    set_pixel(filter,1,1,0,0);
    set_pixel(filter,2,1,0,0);

    set_pixel(filter,0,2,0,1);
    set_pixel(filter,1,2,0,2);
    set_pixel(filter,2,2,0,1);

    return filter;

void feature_normalize(image im)
{
     float min_v=INFINITY;
    float max_v=-INFINITY;

    for (int c=0; c<im.c; c++) {
        for (int h=0; h<im.h; h++) {
            for (int w=0; w<im.w; w++) {
                float v=get_pixel(im,w,h,c);
                if(v<min_v)
                    min_v=v;
                if(v>max_v)
                    max_v=v;
            }
        }
    }
float range=max_v-min_v;
    if(range==0)
        range=1;

    for(int c=0; c<im.c; ++c){
        for (int h=0; h<im.h; ++h){
            for (int w=0; w<im.w; ++w){
                float v =(get_pixel(im, w, h, c) - min_v) / range;
                set_pixel(im,w,h,c,v);
            }
        }
    }
}

}

image *sobel_image(image im)
{
    
   image *results =calloc(2, sizeof(image));

    image gx_filter=make_gx_filter();
    image gy_filter=make_gy_filter();

    image gx=convolve_image(im, gx_filter, 0);
    image gy=convolve_image(im, gy_filter, 0);

    results[0] = make_image(im.w, im.h, 1);
    results[1] = make_image(im.w, im.h, 1);

    for (int i=0; i< im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {
            float gx_v = get_pixel(gx, j, i, 0);
            float gy_v = get_pixel(gy, j, i, 0);
            float mag = sqrtf(gx_v * gx_v + gy_v * gy_v);
            float theta = atan2f(gy_v, gx_v);
            set_pixel(results[0], j, i, 0, mag);
            set_pixel(results[1], j, i, 0, theta);
        }
    }

    free_image(gx_filter);
    free_image(gy_filter);
    free_image(gx);
    free_image(gy);

    return results;
}


image colorize_sobel(image im)
{
   
      image *sobel=sobel_image(im);

    image mag=sobel[0];
    image theta=sobel[1];

    feature_normalize(mag);

    image colorized=make_image(im.w, im.h, 3);

    for (int i = 0; i < im.h; ++i) {
        for (int j = 0; j < im.w; ++j) {
            float hue=(get_pixel(theta, j, i, 0)+TWOPI)/(2*TWOPI);
            float saturation=get_pixel(mag,j,i,0);
            float value = 1;
            set_pixel(colorized,j,i,0,hue);
            set_pixel(colorized,j,i,1,saturation);
            set_pixel(colorized,j,i,2,value);
        }
    }

    hsv_to_rgb(colorized);

    free_image(mag);
    free_image(theta);
    free(sobel);
    return colorized;
}
