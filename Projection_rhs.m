function [FBFGPAFGFx] = Projection_rhs(A,E,G,C_ref,M_f)
% transformace Fx a derivace G*F
% grad(c)
L=chol(C_ref);
%% aplication of Material properties
AFGFx=cat(3,A(:,:,1,1).*E(1)+A(:,:,1,2).*E(2),...
            A(:,:,2,1).*E(1)+A(:,:,2,2).*E(2));
%div(F[A(x)(iF[grad(c)])])


%% Projection operator
%   G*((G'*C*G)^-1)*G';  %
% BtAFGFx_p=cat(3,L(1,1).*AFGFx(:,:,1)+L(2,1).*AFGFx(:,:,2),...
%             L(1,2).*AFGFx(:,:,1)+L(2,2).*AFGFx(:,:,2));

GFAFGFx_p=G.*fftshift(fft2(ifftshift(AFGFx)));
GFAFGFx=GFAFGFx_p(:,:,1)+GFAFGFx_p(:,:,2);

% Inverse operator ((G'*C*G)^-1)
% GC_refG=-(C_ref(1,1).*(G(:,:,1).^2)+C_ref(2,2).*(G(:,:,2).^2)...
%                    +2*C_ref(1,2).*(G(:,:,1).*(G(:,:,2))));
% GC_refG((end+1)/2,(end+1)/2)=1;
GC_refG_inv=M_f.^-1;%GC_refG.^-1;
GC_refG_inv((end+1)/2,(end+1)/2)=0;
%GC_refG_inv;
%% 
GPAFGFx=G.*GC_refG_inv.*GFAFGFx;

GPAFGFx((end+1)/2,(end+1)/2,:)=0;

FGPAFGFx=fftshift(ifft2(ifftshift(GPAFGFx)));

BFGPAFGFx=cat(3,C_ref(1,1).*FGPAFGFx(:,:,1)+C_ref(1,2).*FGPAFGFx(:,:,2),...
            C_ref(2,1).*FGPAFGFx(:,:,1)+C_ref(2,2).*FGPAFGFx(:,:,2));
FBFGPAFGFx=fftshift(fft2(ifftshift(BFGPAFGFx)));


end