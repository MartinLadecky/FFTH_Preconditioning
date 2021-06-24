%% Example 3: Effect of phase anisotropy
%  experiment 4.3 in the paper
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

counter=1;
for loop=1:2:20
  par=loop 
N_1=301; %  number of points in x_1
N_2=N_1; % number of points in x_2

%% Mesh parameters
h_1=2*pi/(N_1); % step in x_2
h_2=2*pi/(N_2); % step in x_2

%% Coordinates
x=zeros(N_2,N_1,2);
[x(:,:,1),x(:,:,2)]=meshgrid(-pi+h_1/2:h_1:pi-h_1/2,-pi+h_2/2:h_2:pi-h_2/2);
%% Derivatives
G_n=G_matrix(N_1,N_2);
%% Material coeficient matrix

 Pixels = imread('structure_3.png');

 pixa=round(linspace(1,size(Pixels,2),N_1));
 piya=round(linspace(1,size(Pixels,1),N_2));
 
  
  Min_eig = 100000;
 Max_eig= 0;
 A=zeros(N_2,N_1,2,2);  
  for i=1:N_2
     for j=1:N_1    
         [A(i,j,:,:),kappapom(i,j),eigenos(i,j,:)]=a_matrix_img_aniso(Pixels(piya(i),pixa(j)),par);
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
 Pomer1=max(max(kappapom));
 Pomer2=min(min(kappapom));
 Pomer12(counter)=Pomer1./Pomer2
 
 
d=[mean(mean(A(:,:,1,1))) mean(mean(A(:,:,1,2))) ;...
   mean(mean(A(:,:,2,1))) mean(mean(A(:,:,2,2)))];

%% Gradient matrices
G  =G_clasic(N_1,N_2);
G_n=G_normed(N_1,N_2,[1 0;0 1]);
G_m=G_normed(N_1,N_2,d);

%% Preconditioning with constant frequency: identity 
[M_n] = M_mean_half_inv(N_1,N_2,[1 0;0 1]);

%% Preconditioning with constant frequency: mean value 
[M_m] = M_mean_half_inv(N_1,N_2,d);

%% Preconditioning with significant frequencie
[U1,U2] = precon_decomposition(A,G,N_1,N_2);

%% Conditions:
steps = 1000;
toler = 1e-6;
c_000=rand(N_2,N_1);

%% SOLVER without preconditionig
c_0 = c_000;
for k=1:1
    E=E_0(:,k); 
    tic;
    %[C,st]=CG_solver(A,G_n,c_0,E,steps,toler); % without preconditioning
    [C1,st,norm_evol1]=solver_CG_Gn(A,G_n,c_0,E,steps,toler,M_n);
    T1(k,counter) = toc;
    S1(k,counter) = st;
    %A_01(:,k)=Hom_parameter(C,A,G_n,E); % Compute homogenized parameter
end

%% SOLVER with constant preconditionig
c_0=c_000;
for k=1:1
    E=E_0(:,k);
    tic;
   % [C,st]=CGP_solver_constant(A,G_n,c_0,E,steps,toler,GG0);% with preconditioning
    [C2,st,norm_evol2]=solver_CG_Gn(A,G_m,c_0,E,steps,toler,M_m);
    T2(k,counter) = toc;
    S2(k,counter) = st;
    %A_02(:,k)=Hom_parameter(C,A,G_n,E);% Compute homogenized parameter
end
%% SOLVER with efficient enhanced preconditionig
c_0=c_000;
for k=1:1
    tic;
    E=E_0(:,k);
    %[C,st]=CGP_solver_1f_v3(A,G_n,c_0,E,steps,toler,GG0,GG1,U1,U2); % with better preconditioning
    [C3,st,norm_evol3]=solver_PCG_1f_left(A,G,c_0,E,steps,toler,U1,U2); % with better preconditioning
    T3(k,counter)=toc;
    S3(k,counter) = st;
    %A_03(:,k)=Hom_parameter(C,A,G_n,E);% Compute homogenized parameter
end
%% error 
   NoP(1,counter)=par;
   Kappa(counter,:)=[Min_eig,Max_eig,Max_kappa]
  counter=counter+1;      
    
end

%% Plot 
% Plot steps
 figure 
 hold on
 plot(NoP,S1(1,:),'r')
 plot(NoP,S2(1,:),'b')
 plot(NoP,S3(1,:),'black')
 legend('G_n','G_m','M_2')
 
% Plot times
 figure 
 hold on
 plot(NoP,T1(1,:),'r')
 plot(NoP,T2(1,:),'b')
 plot(NoP,T3(1,:),'black')
 legend('G_n','G_m','M_2')
%% Save plot data 
 if ~exist('experiment_data/exp3', 'dir')
       mkdir('experiment_data/exp3')
end
save('experiment_data/exp3/NoP.mat','NoP');
save('experiment_data/exp3/S1.mat','S1');
save('experiment_data/exp3/S2.mat','S2');
save('experiment_data/exp3/S3.mat','S3');
save('experiment_data/exp3/T1.mat','T1');
save('experiment_data/exp3/T2.mat','T2');
save('experiment_data/exp3/T3.mat','T3'); 
 
 
 