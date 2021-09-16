function [GPAFGFx] = Projection_plain_Cref(GFx,G,M_f,C_ref)

C_ref=chol(C_ref);

FGFx=fftshift(ifft2(ifftshift(GFx)));
%% Projection operator with Cref on 
CFGFx=cat(3,C_ref(1,1).*FGFx(:,:,1)+C_ref(1,2).*FGFx(:,:,2),...
            C_ref(2,1).*FGFx(:,:,1)+C_ref(2,2).*FGFx(:,:,2));

        
%% Projection operator

GFAFGFx_p=conj(G).*fftshift(fft2(ifftshift(CFGFx)));
%GFAFGFx_p=G.*GFx;
GFAFGFx=GFAFGFx_p(:,:,1)+GFAFGFx_p(:,:,2);

GC_refG_inv=M_f.^-1;
GC_refG_inv((end+1)/2,(end+1)/2)=0;


GPAFGFx=G.*GC_refG_inv.*GFAFGFx;
GPAFGFx((end+1)/2,(end+1)/2,:)=0;

FGPAFGFx=fftshift(ifft2(ifftshift(GPAFGFx)));

GPAFGFx=cat(3,C_ref(1,1)'.*FGPAFGFx(:,:,1)+C_ref(1,2)'.*FGPAFGFx(:,:,2),...
            C_ref(2,1)'.*FGPAFGFx(:,:,1)+C_ref(2,2)'.*FGPAFGFx(:,:,2));

end
