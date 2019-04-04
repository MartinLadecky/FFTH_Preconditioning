function [M] = M_mean(N_bf1,N_bf2,C)
%% Input
% N_bf1 []- Number of bf in x dir,
% N_bf2 []- Number of bf in y dir,
%% Output
% G_n -Matrix of coeficients of 1st derivative 
%       saved in "meshgrid size" 
% G_n [:,:,1]- der in x direction
% G_n [:,:,2]- der in y direction
%% 
[m(:,:,1),m(:,:,2)] =meshgrid(-(N_bf1-1)/2:(N_bf1-1)/2,...
                              -(N_bf2-1)/2:(N_bf2-1)/2);
M=zeros(N_bf2,N_bf1);
for i=1:N_bf2
    for j=1:N_bf1
       M(i,j)=1./(pi*sqrt([m(i,j,1),m(i,j,2)]*C*[m(i,j,1),m(i,j,2)]'));
    end
end
M((N_bf2+1)/2,(N_bf1+1)/2)=1; % Zero division ==>0 
end