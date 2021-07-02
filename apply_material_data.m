function [FAFGFx] = apply_material_data(A,GFx)

if  size(GFx,1) == 2

AFGFx=cat(3,A(:,:,1,1).*GFx(1)+A(:,:,1,2).*GFx(2),...
            A(:,:,2,1).*GFx(1)+A(:,:,2,2).*GFx(2));
FAFGFx=fftshift(fft2(ifftshift(AFGFx)));

else
    
%% Projection operator
FGFx=fftshift(ifft2(ifftshift(GFx)));
AFGFx=cat(3,A(:,:,1,1).*FGFx(:,:,1)+A(:,:,1,2).*FGFx(:,:,2),...
            A(:,:,2,1).*FGFx(:,:,1)+A(:,:,2,2).*FGFx(:,:,2));
FAFGFx=fftshift(fft2(ifftshift(AFGFx)));
end

