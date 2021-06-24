function [val] = scalar_product_grad_energy(a,b,A)

b=fftshift(ifft2(ifftshift(b)));
%% aplication of Material properties
Ab=cat(3,A(:,:,1,1).*b(:,:,1)+A(:,:,1,2).*b(:,:,2),...
            A(:,:,2,1).*b(:,:,1)+A(:,:,2,2).*b(:,:,2));
%div(F[A(x)(iF[grad(c)])])


%% Projection operator
%   G*((G'*C*G)^-1)*G';  %

Ab=fftshift(fft2(ifftshift(Ab)));



val_1=sum(sum(a(:,:,1).*(Ab(:,:,1).')'));
val_2=sum(sum((a(:,:,2).')'.*Ab(:,:,2)));
val=val_1+val_2;
end