
%function [A_0,st,t,Nbf]=Hom_solver(N)
%% HOMOGENIZATION PROBLEM
%  sol PDE -div(A(x)grad(u))=div(A(x)E))
%  then homogenized A_0E=int(A(x)(E+grad(u)))dx/(volume of omega)
%% Input
% N   [1]    -number of points in sample
%% Output
% A_0 [2,2] -Homogenized material parameters
% st  [1]   -Number of iteration steps
% t   [1]   -Time
%% BEGIN
clc;clear;

tic
E_0=[1 0;0 1];
A_0=zeros(2,2);
loop1=1;
counter=1;
figure
%[1,5,10,14,17,19,20]
for N_1=[3]% number of points in x_1-1   3^loop

    N_2=N_1; % number of points in x_2
    N=N_2*N_1
    %% Mesh parameters
    h_1=2*pi/(N_1); % step in x_2
    h_2=2*pi/(N_2); % step in x_2

    %% Coordinates
    x=zeros(N_2,N_1,2);
    [x(:,:,1),x(:,:,2)]=meshgrid(-pi+h_1/2:h_1:pi-h_1/2,-pi+h_2/2:h_2:pi-h_2/2);
    %% Material coeficient matrix and data analysis
     Min_eig = 100000;
     Max_eig= 0;
     A=zeros(N_2,N_1,2,2);
     for i=1:N_2
         for j=1:N_1    
             A(i,j,:,:)=a_matrix(x(i,j,:));

             pom = zeros(2);
             pom(1,1) = A(i,j,1,1); pom(1,2) = A(i,j,1,2);
             pom(2,1) = A(i,j,2,1); pom(2,2) = A(i,j,2,2);        
             p = eig(pom);
             min_eig = min(p);
             max_eig = max(p);
             if (min_eig<Min_eig) 
                 Min_eig=min_eig; 
             end
             if (max_eig>Max_eig)
                 Max_eig=max_eig;
             end        
         end       
     end
     Max_kappa = Max_eig/Min_eig;
     [Min_eig,Max_eig,Max_kappa]
     if Min_eig<0

         return
     end

    d=[mean(mean(A(:,:,1,1))) mean(mean(A(:,:,1,2))) ;...
       mean(mean(A(:,:,2,1))) mean(mean(A(:,:,2,2)))];

    %% Derivatives matrices
    G  =G_clasic(N_1,N_2);
    G_n=G_normed(N_1,N_2,[1 0;0 1]);
    G_m=G_normed(N_1,N_2,d);

    %% Preconditioning with constant frequencie
    M_f_const= -(d(1,1).*(G(:,:,1).^2)+d(2,2).*(G(:,:,2).^2)...
                       +2*d(1,2).*(G(:,:,1).*(G(:,:,2))));                   
    M_f_const((end+1)/2,(end+1)/2)=-1;

    [M_n] = M_mean(N_1,N_2,[1 0;0 1]);
    [M_m] = M_mean(N_1,N_2,d);


    %% Conditions:
    steps =300;
    toler = 1e-10;
    c_000=zeros(N_2,N_1);

    %% SOLVERs with M_1 preconditionig
    c_0 = c_000;
    K=zeros(N,N);
    PK=zeros(N,N);
    counter=1;
    for j = 1:N_1 
        for i = 1:N_2
            v_1=zeros(N_2,N_1);
            v_1(i,j)=1;
            row=LHS(A,v_1,G);
            Prow=PLHS(A,v_1,G,M_m);
            K(counter,:)=reshape(row,[N,1]);
            PK(counter,:)=reshape(Prow,[N,1]);     
            counter=counter+1;      
        end
    end
    %K=-K
    e_1=-eig(K);
    %e_1=e_1/max(e_1);

    e_2=-eig(PK);
    %e_2=e_2/max(e_2);

    x_0=ones(1,N)';

    %figure(1)
    %plot(e_1,x_0+2+0.1*loop1, 'b.','LineWidth',1)
    %hold on;
    %plot(e_2,x_0+1+0.1*loop1, 'r.','LineWidth',1)
    
    figure(5)
    hold on;
    e1_size=size(e_1);
    scatter((1:1: e1_size(1)),sort(e_1)+1, 'bo','LineWidth',1)
    e2_size=size(e_2);
    scatter((1:1: e2_size(1)),sort(e_2)+1, 'rx','LineWidth',1)
    set(gca,'YScale','log')
    loop1=loop1+1;
end
%figure(1)
%legend('K','PK')
%ylim([1, 4])
%xlim([0, 1])
figure(5)
legend('K','PK')
