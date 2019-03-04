function a=a_matrix_img(pix)
if pix < 130
    a =[0.01 ,   0.005;...
        0.005 ,   0.01];

else 
    a =[2,  1.3;...
        1.3, 2];
    
    
 kappapom=max(eig(a))/min(eig(a));
 eigenos=eig(a);

end