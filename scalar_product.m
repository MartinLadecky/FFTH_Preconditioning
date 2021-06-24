function [val] = scalar_product(a,b)
val=sum(sum((a.')'.*b));
end