function a=a_matrix_imgx(pix)
    a =[(0.01+(20*((255-pix)/255))) ,   (0.05+(1.3*((255-pix)/255)));...
        (0.05 +(1.3*((255-pix)/255))), (  0.01+(2*((255-pix)/255)))];


    
 kappapom=max(eig(a))/min(eig(a));
 eigenos=eig(a);

end