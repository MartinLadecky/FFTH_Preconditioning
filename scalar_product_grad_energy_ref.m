function [val] = scalar_product_grad_energy_ref(a,b,A)

b=fftshift(ifft2(ifftshift(b))); % iFFT

Ab=cat(3,A(1,1).*b(:,:,1)+A(1,2).*b(:,:,2),...
            A(2,1).*b(:,:,1)+A(2,2).*b(:,:,2));
        
Ab=fftshift(fft2(ifftshift(Ab))); % FFT

val_1=(a(:,:,1).')'.*Ab(:,:,1);
val_2=(a(:,:,2).')'.*Ab(:,:,2);
val=sum(sum(val_1+val_2));
end