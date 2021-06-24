function [A_0E] = Hom_parameter_grad(C_grad,A,G,E)
%%  A_0E=int(A(x)(E+grad(u)))dx/(volume of omega)
%% input
% A   [N_bf2,N_bf1,2,2] -matrix of material parameters in every point of grid
% G_n [N_bf2,N_bf1,2]   -matrix of coeficients of 1st derivative
% C   [N_bf2,N_bf1]     -solution in frequencies 
% E   [1,2]             -vector [0;1] or [1;0]
%% Output
% A_0*E [2,1]           -Part of solution A_0 
%% 
Gnu=C_grad; % grad(c)
gradu=fftshift(ifft2(ifftshift(Gnu))); % iF[grad(c)]

Egradu=cat(3,gradu(:,:,1)+E(1),gradu(:,:,2)+E(2)); % E+iF[grad(c)]

a_1E=cat(3,A(:,:,1,1).*Egradu(:,:,1)+A(:,:,1,2).*Egradu(:,:,2),...
           A(:,:,2,1).*Egradu(:,:,1)+A(:,:,2,2).*Egradu(:,:,2)); %A(x)(E+iF[grad(c)])
%        Egradu
%        return

A_0E= ([sum(sum(a_1E(:,:,1))),sum(sum(a_1E(:,:,2)))]./numel(C_grad(:,:,1)))'; % num inegration
end