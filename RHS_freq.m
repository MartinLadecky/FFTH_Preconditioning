function [RHS_0] = RHS_freq(A,E,G_nfft)
AE=cat(3,A(:,:,1,1).*E(1)+A(:,:,1,2).*E(2),...
         A(:,:,2,1).*E(1)+A(:,:,2,2).*E(2));
%  ifftshift(AE(:,:,1))
%  ifftshift(AE(:,:,2))
%  ifftshift(AE)
%  return
FAE=fftshift(fft2(ifftshift(AE)));
% return
  
GFAE=G_nfft.*FAE; 

RHS_0=GFAE(:,:,1)+GFAE(:,:,2);

end