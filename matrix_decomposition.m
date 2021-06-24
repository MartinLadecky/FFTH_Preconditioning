a= 5
b=3
c=2
A = [a 0 0 c 0 0;
     0 a 0 0 c 0;
     0 0 a 0 0 c;
     c 0 0 b 0 0;
     0 c 0 0 b 0;
     0 0 c 0 0 b]
 eig(A)
 
L = chol(A)
L'*L

B= [a c;
     c b]
 Lb = chol(B)


