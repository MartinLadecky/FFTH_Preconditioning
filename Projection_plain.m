function [GPAFGFx] = Projection_plain(GFx,G,M_f)
%% Projection operator

GFAFGFx_p=conj(G).*GFx;
%GFAFGFx_p=G.*GFx;
GFAFGFx=GFAFGFx_p(:,:,1)+GFAFGFx_p(:,:,2);

GC_refG_inv=M_f.^-1;
GC_refG_inv((end+1)/2,(end+1)/2)=0;


GPAFGFx=G.*GC_refG_inv.*GFAFGFx;
GPAFGFx((end+1)/2,(end+1)/2,:)=0;

end