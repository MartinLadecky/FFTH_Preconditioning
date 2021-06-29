function [val] = scalar_product_grad(a,b)

val_1=(b(:,:,1).')'.*a(:,:,1);
val_2=(b(:,:,2).')'.*a(:,:,2);
val=sum(sum(val_1+val_2));
end