function[a,kappapom,eigenos]=a_matrix_img_aniso(pix,par)

if pix < 130
    a =[ 1,  0 ;...
         0, par ];
else 
    a =[10,0;...
         0,  10*par];
end


  kappapom=max(eig(a))/min(eig(a));
 eigenos=eig(a)
end
