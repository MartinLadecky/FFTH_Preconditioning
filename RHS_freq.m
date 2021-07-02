function [RHS_0] = RHS_freq(A,E,G)

AE=cat(3,A(:,:,1,1).*E(1)+A(:,:,1,2).*E(2),...
         A(:,:,2,1).*E(1)+A(:,:,2,2).*E(2));

FAE=fftshift(fft2(ifftshift(AE)));

GFAE=conj(G).*FAE; 

RHS_0=GFAE(:,:,1)+GFAE(:,:,2);

end