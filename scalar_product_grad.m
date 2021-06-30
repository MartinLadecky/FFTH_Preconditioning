function [val] = scalar_product_grad(a,b)
a=conj(a);
val_1=a(:,:,1).*b(:,:,1);
val_2=a(:,:,2).*b(:,:,2);
val=sum(sum(val_1+val_2));
end