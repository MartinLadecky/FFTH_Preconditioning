function a=a_matrix_img(pix)
if pix < 130
    a =[0.04 ,   0.008;...
        0.008 ,   0.04];

else 
    a =[2,  1.3;...
        1.3, 2];
    
    
 kappapom=max(eig(a))/min(eig(a));
 eigenos=eig(a);

end