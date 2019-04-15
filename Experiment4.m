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
for loop=(5:2:10)%[1,5,10,11]%,20,21,22
N_1=2*(loop^2)+1%2*(loop^2)+1% number of points in x_1-1

N_2=N_1; % number of points in x_2

%% Mesh parameters
h_1=2*pi/(N_1); % step in x_2
h_2=2*pi/(N_2); % step in x_2

%% Coordinates
x=zeros(N_2,N_1,2);
[x(:,:,1),x(:,:,2)]=meshgrid(-pi+h_1/2:h_1:pi-h_1/2,-pi+h_2/2:h_2:pi-h_2/2);

%% Material coeficient matrix

 Pixels = imread('structure_3.png');

 pixa=round(linspace(1,size(Pixels,2),N_1));
 piya=round(linspace(1,size(Pixels,1),N_2));
 
  
  Min_eig = 100000;
 Max_eig= 0;
 A=zeros(N_2,N_1,2,2);  
  for i=1:N_2
     for j=1:N_1    
          A(i,j,:,:)=a_matrix_img(Pixels(piya(i),pixa(j)));
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

%% Material ananlysis
d=[mean(mean(A(:,:,1,1))) mean(mean(A(:,:,1,2))) ;...
   mean(mean(A(:,:,2,1))) mean(mean(A(:,:,2,2)))];

%% Derivatives
G=G_clasic(N_1,N_2);
G_n=G_matrix(N_1,N_2);
G_m=G_mean(N_1,N_2,d);

%% Preconditioning
[M_m] = M_mean(N_1,N_2,d);

M_fGn_const= -(d(1,1).*(G_n(:,:,1).^2)+d(2,2).*(G_n(:,:,2).^2)...
                   +2*d(1,2).*(G_n(:,:,1).*(G_n(:,:,2))));  
               
M_fGn_const((end+1)/2,(end+1)/2)=1;

M_fG_const= -(d(1,1).*(G(:,:,1).^2)+d(2,2).*(G(:,:,2).^2)...
                   +2*d(1,2).*(G(:,:,1).*(G(:,:,2))));  
               
M_fG_const((end+1)/2,(end+1)/2)=-1;

%% Conditions:
steps = 300;
toler = 1e-6;
c_000=rand(N_2,N_1);
%% SOLVER without preconditionig and G
c_0 = c_000;
for k=1:1
    E=E_0(:,k); 
    tic;
    [C,st,norm_evol1]=CG_solver(A,G_m,c_0,E,steps,toler,M_m);
    T1(k,counter) = toc;
    S1(k,counter) = st;
    A_3(:,k)=Hom_parameter(C,A,G,E); % Compute homogenized parameter
end
%% SOLVER with symetric preconditionig M and G
c_0 = c_000;
for k=1:1
    E=E_0(:,k); 
    tic;
    [C,st,norm_evol2]=CG_solver_symPrec(A,G,c_0,E,steps,toler,M_m);

    T2(k,counter) = toc;
    S2(k,counter) = st;
    A_3(:,k)=Hom_parameter(C,A,G,E); % Compute homogenized parameter
end

%% SOLVER with constant preconditionig
c_0=c_000;
for k=1:1
    E=E_0(:,k);
    tic;
    [C,st,norm_evol3]=CGP_solver_left(A,G,c_0,E,steps,toler,M_fG_const);% with preconditioning
    T3(k,counter)=toc;
    S3(k,counter) = st;
    A_3(:,k)=Hom_parameter(C,A,G,E);% Compute homogenized parameter

end


% error 
  NoP(counter)=N_1*N_2;
  counter=counter+1;
end
%% 
%save('experiment_data/exp2/NoP.mat','NoP');
% save('experiment_data/exp2/S1.mat','S1');
% save('experiment_data/exp2/S2.mat','S2');
% save('experiment_data/exp2/S3.mat','S3');
% save('experiment_data/exp2/T1.mat','T1');
% save('experiment_data/exp2/T2.mat','T2');
% save('experiment_data/exp2/T3.mat','T3'); 

 figure 
 hold on
 plot((1:S1(1,end)) ,norm_evol1)
 plot((1:S2(1,end)),norm_evol2)
 plot((1:S3(1,end)),norm_evol3)
set(gca, 'XScale', 'log', 'YScale', 'log');
legend('in G','Symetric','Left')
%% Plot
% Plot steps
figure 
 hold on
 plot(NoP,S1(1,:),'r')
 plot(NoP,S2(1,:),'--k')
 plot(NoP,S3(1,:))
legend('in G','Symetric','Left')
 
% Plot times
 figure 
 hold on
 plot(NoP,T1(1,:),'r')
 plot(NoP,T2(1,:),'b')
 plot(NoP,T3(1,:),'black')
legend('in G','Symetric','Left')

 %% Plot A function