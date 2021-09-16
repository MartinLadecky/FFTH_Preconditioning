function [GPAFGFx] = Projection_Cref(A,GFx,G,M_f,C_ref)

FGFx=fftshift(ifft2(ifftshift(GFx)));
%% aplication of Material properties
AFGFx=cat(3,A(:,:,1,1).*FGFx(:,:,1)+A(:,:,1,2).*FGFx(:,:,2),...
            A(:,:,2,1).*FGFx(:,:,1)+A(:,:,2,2).*FGFx(:,:,2));


CAFGFx=cat(3,C_ref(1,1).*AFGFx(:,:,1)+C_ref(1,2).*AFGFx(:,:,2),...
    C_ref(2,1).*AFGFx(:,:,1)+C_ref(2,2).*AFGFx(:,:,2));
%% Projection operator

GFAFGFx_p=conj(G).*fftshift(fft2(ifftshift(CAFGFx)));
GFAFGFx=GFAFGFx_p(:,:,1)+GFAFGFx_p(:,:,2);

GC_refG_inv=M_f.^-1;
GC_refG_inv((end+1)/2,(end+1)/2)=0;

%% 
GPAFGFx=G.*GC_refG_inv.*GFAFGFx;

GPAFGFx((end+1)/2,(end+1)/2,:)=0;

end