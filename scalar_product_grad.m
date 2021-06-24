function [val] = scalar_product_grad(a,b)

val_1=sum(sum(a(:,:,1).*(b(:,:,1).')'));
val_2=sum(sum((a(:,:,2).')'.*b(:,:,2)));
val=val_1+val_2;
end