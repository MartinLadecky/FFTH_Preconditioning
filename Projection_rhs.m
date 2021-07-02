function [GPAFGFx] = Projection_rhs(A,E,G,M_f)
%% aplication of Material properties
AFGFx=cat(3,A(:,:,1,1).*E(1)+A(:,:,1,2).*E(2),...
            A(:,:,2,1).*E(1)+A(:,:,2,2).*E(2));

%% Projection operator


GFAFGFx_p=conj(G).*fftshift(fft2(ifftshift(AFGFx)));
GFAFGFx=GFAFGFx_p(:,:,1)+GFAFGFx_p(:,:,2);

GC_refG_inv=M_f.^-1;%GC_refG.^-1;
GC_refG_inv((end+1)/2,(end+1)/2)=0;


GPAFGFx=G.*GC_refG_inv.*GFAFGFx;
GPAFGFx((end+1)/2,(end+1)/2,:)=0;

end