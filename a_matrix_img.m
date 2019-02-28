function a=a_matrix_img(pix)
if pix < 130
    a =[0.01 ,   0.000;...
        0.000 ,   0.01];

else 
    a =[2,  0.9;...
        0.9, 2];
    
    
 kappapom=max(eig(a))/min(eig(a));
 eigenos=eig(a)

end