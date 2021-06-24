function [val] = scalar_product_energy(a,b,A)
val=sum(sum((a.')'.* A.*b));
end