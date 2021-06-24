%% Example 1: Significant frequency in material
%  experiment 4.1 in the paper
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

for loop=[1,2,3,4,5,10,11,12]%[1,5,10,14,17,19,20]
N_1=2*(loop^2)+1% number of points in x_1-1   3^loop
N_2=N_1; % number of points in x_2

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
steps =300;
toler = 1e-6;
c_000=rand(N_2,N_1);

%% SOLVERs with M_0 preconditionig
% Preconditioner incorporated in G matricies
c_0 = c_000;
for k=1:1 
    E=E_0(:,k);
    tic;
    [C1,st,norm_evol1]=CG_solver(A,G,c_0,E,steps,toler,M_n);
    T1(k,counter) = toc;
    S1(k,counter) = st;
    A_1(:,k)=Hom_parameter(C1,A,G,E); % Compute homogenized parameter
end
%% SOLVERs with M_1 preconditionig
c_0 = c_000;
for k=1:1 
    E=E_0(:,k);
    tic;
    [C2,st,norm_evol2]=solver_CG_Gn(A,G_m,c_0,E,steps,toler,M_m);
    T2(k,counter) = toc;
    S2(k,counter) = st;
    A_2(:,k)=Hom_parameter(C2,A,G,E); % Compute homogenized parameter
end
%% SOLVER with M_2 preconditionig
c_0=c_000;
for k=1:1
    E=E_0(:,k);
    tic;
    [C3,st,norm_evol3]=solver_PCG_1f_left(A,G,c_0,E,steps,toler,U1,U2); % with better preconditioning
    T3(k,counter)=toc;
    S3(k,counter) = st;
    A_3(:,k)=Hom_parameter(C3,A,G,E);% Compute homogenized parameter
end
%% Plot
 figure 
 hold on
 plot((1:S1(1,end)) ,norm_evol1)
 plot((1:S2(1,end)),norm_evol2)
% loglog((1:S22(1,end)),norm_evol22)
 loglog((1:S3(1,end)),norm_evol3)
set(gca, 'XScale', 'linear', 'YScale', 'log');
legend('G','G_m','M_2')
title('Norm of residua')
%% error 
  NoP(counter)=N_1*N_2;
  counter=counter+1;
end


NoS1=(1:S1(1,end));
NoS2=(1:S2(1,end));
NoS3=(1:S3(1,end));

%% Plot results
 figure 
 hold on
 loglog((1:S1(1,end)) ,norm_evol1)
 loglog((1:S2(1,end)),norm_evol2)
 loglog((1:S3(1,end)),norm_evol3)
set(gca, 'XScale', 'log', 'YScale', 'log');
legend('$\widehat{K}^{\textrm{ref}}_{0}$','$\widehat{K}^{\textrm{ref}}_{1}$','$\widehat{K}^{\textrm{ref}}_{2}$','Interpreter','latex')

% Plot steps
figure 
 hold on
 plot(NoP,S1(1,:),'r')
 plot(NoP,S2(1,:),'b')
 plot(NoP,S3(1,:),'black')
 % plot(NoP,S4(1,:),'black')
legend('$\widehat{K}^{\textrm{ref}}_{0}$','$\widehat{K}^{\textrm{ref}}_{1}$','$\widehat{K}^{\textrm{ref}}_{2}$','Interpreter','latex')
axis([ 0 1e6 0 50])
%yticks([0 5 10 15 20 25 30])
%xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
ylabel('number of iterations')
xlabel('grid size - $ | | $ ','Interpreter','latex')
saveas(gcf,'exp1_steps.eps')

% Plot times
 figure 
 hold on
 plot(NoP,T1(1,:),'r')
 plot(NoP,T2(1,:),'b')
 plot(NoP,T3(1,:),'black')
legend('$\widehat{K}_{0}$','Interpreter','latex')
 
 %% Save plot data 
if ~exist('experiment_data/exp1', 'dir')
       mkdir('experiment_data/exp1')
end

save('experiment_data/exp1/NoP.mat','NoP');
save('experiment_data/exp1/S1.mat','S1');
save('experiment_data/exp1/S2.mat','S2');
save('experiment_data/exp1/S3.mat','S3');
save('experiment_data/exp1/T1.mat','T1');
save('experiment_data/exp1/T2.mat','T2');
save('experiment_data/exp1/T3.mat','T3');

if ~exist('experiment_data/expPAMM1', 'dir')
       mkdir('experiment_data/expPAMM1')
end
save('experiment_data/expPAMM1/NoS1.mat','NoS1'); 
save('experiment_data/expPAMM1/NoS2.mat','NoS2');
save('experiment_data/expPAMM1/NoS3.mat','NoS3');

save('experiment_data/expPAMM1/S1.mat','S1');
 save('experiment_data/expPAMM1/S2.mat','S2');
 save('experiment_data/expPAMM1/S3.mat','S3');
save('experiment_data/expPAMM1/T1.mat','norm_evol1');
save('experiment_data/expPAMM1/T2.mat','norm_evol2');
save('experiment_data/expPAMM1/T3.mat','norm_evol3'); 
 