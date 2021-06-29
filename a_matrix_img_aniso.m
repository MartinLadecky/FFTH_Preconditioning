function[a,kappapom,eigenos]=a_matrix_img_aniso(pix,par)

if pix < 130
    a =[ 1, 0.5 ;...
         0.5, 5 ];
else 
    a =[ 4,  2;...
         2,  10*par];
end


  kappapom=max(eig(a))/min(eig(a));
 eigenos=eig(a);
end
