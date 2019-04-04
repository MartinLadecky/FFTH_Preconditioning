function [G] = G_clasic(N_bf1,N_bf2)
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
G=zeros(N_bf2,N_bf1,2);
for i=1:N_bf2
    for j=1:N_bf1
       G(i,j,:)=pi*1i*m(i,j,:);
    end
end
G((N_bf2+1)/2,(N_bf1+1)/2,:)=0; % Zero division ==>0 
end