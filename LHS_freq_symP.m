function [Mc] = LHS_freq_symP(A,c,G,M_m)
% transformace Fx a derivace G*F
% grad(c)

GFx=G.*M_m.*c; 
%iF[grad(c)]
FGFx=fftshift(ifft2(ifftshift(GFx)));
%A(x)(iF[grad(c)])
AFGFx=cat(3,A(:,:,1,1).*FGFx(:,:,1)+A(:,:,1,2).*FGFx(:,:,2),...
            A(:,:,2,1).*FGFx(:,:,1)+A(:,:,2,2).*FGFx(:,:,2));
%div(F[A(x)(iF[grad(c)])])
GAFGFx=G.*fftshift(fft2(ifftshift(AFGFx)));
Mc=M_m.*(GAFGFx(:,:,1)+GAFGFx(:,:,2)); 

end