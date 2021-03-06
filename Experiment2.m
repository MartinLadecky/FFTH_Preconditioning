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
for loop=[1,5,10,14,17,19,20]
N_1=2*(loop^2)+1% number of points in x_1-1
N_2=N_1; % number of points in x_2

%% Mesh parameters
h_1=2*pi/(N_1); % step in x_2
h_2=2*pi/(N_2); % step in x_2

%% Coordinates
x=zeros(N_2,N_1,2);
[x(:,:,1),x(:,:,2)]=meshgrid(-pi+h_1/2:h_1:pi-h_1/2,-pi+h_2/2:h_2:pi-h_2/2);

%% Material coeficient matrix and data analysis

 Pixels = imread('structure_4.png');
  %structure = matfile('structure11.mat');
  %Pixels = structure.newimage;

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
 
d=[mean(mean(A(:,:,1,1))) mean(mean(A(:,:,1,2))) ;...
   mean(mean(A(:,:,2,1))) mean(mean(A(:,:,2,2)))];

%% Deri%% Derivatives matrices
G  =G_clasic(N_1,N_2);
G_n=G_normed(N_1,N_2,[1 0;0 1]);
G_m=G_normed(N_1,N_2,d);
%% Preconditioning with constant frequencie
M_f_const= -(d(1,1).*(G(:,:,1).^2)+d(2,2).*(G(:,:,2).^2)...
                   +2*d(1,2).*(G(:,:,1).*(G(:,:,2))));                   
M_f_const((end+1)/2,(end+1)/2)=-1;

[M_n] = M_mean(N_1,N_2,[1 0;0 1]);
[M_m] = M_mean(N_1,N_2,d);


%% Preconditioning with significant frequencie
A_freq11=fft2(ifftshift(A(:,:,1,1)));
A_freq22=fft2(ifftshift(A(:,:,2,2)));
A_freq12=fft2(ifftshift(A(:,:,1,2)));
    % 0 frequency
d11 = A_freq11(1,2); % nad diagonalu v matici FAF !
d22 = A_freq22(1,2);
d12 = A_freq12(1,2);

       
dd11 = A_freq11(1,1);
dd22 = A_freq22(1,1);
dd12 = A_freq12(1,1);

G_1 = G(:,:,1);
G_2 = G(:,:,2);
G_1p = [G_1(:,2:end),G_1(:,1)];% posunuta matica derivaci
G_2p = [G_2(:,2:end),G_2(:,1)];

GG0 = (dd11*G_1.*G_1 + dd22*G_2.*G_2 + dd12*(G_2.*G_1 + G_1.*G_2))/N_1/N_2;
GG0((end+1)/2,(end+1)/2)=-1; %Predpodmienovacia 0 f
GG0 = -GG0;

GG1 = -conj((d11*conj(G_1).*G_1p + d22*conj(G_2).*G_2p + ...
      d12*conj(G_1).*G_2p + d12*conj(G_2).*G_1p))/N_1/N_2;
GG1 = -GG1;%Predpodmienovacia 1 f

P = zeros(N_1,N_1,N_2); % sada predpodminovacich matic
for r = 1:N_2    
    m1 = GG0(r,:); % radek Predpod matice 0f
    m2 = GG1(r,:); % radek Predpod matice 1f  
    p = spdiags([m2(1:N_1-1)'],[-1],N_1,N_1); % zoradi m2 riadok do poddiagonaly, riedka matica
    p = p + p' + diag(m1);     
    P(:,:,r) = p;% zaradi to do glob matice
end  

U1 = zeros(N_2,N_1);
U2 = U1;
for r = 1:N_2    
    m1 = GG0(r,:);
    m2 = GG1(r,:);    
    p = spdiags([m2(1:N_1-1)'],[-1],N_1,N_1);
%     p(1,N_1) = m2(N_1);  % ZDE OPRAVA NEBO NE
    p = p + p' + diag(m1); 
    pom = chol(p);  % p = pom'*pom % chol rozklad
    U1(r,:) = diag(pom); % diagonala chol
    U2(r,1:N_1-1) = diag(pom,1); % mimodiagonala chol
end

%% Conditions:
steps = 500;
toler = 1e-6;
c_000=rand(N_2,N_1);
%% SOLVERs with M_0 preconditionig
% Preconditioner incorporated in G matricies
c_0 = c_000;
for k=1:1 
    E=E_0(:,k);
    tic;
    [C1,st,norm_evol1]=CG_solver(A,G_n,c_0,E,steps,toler,M_n);
    T1(k,counter) = toc;
    S1(k,counter) = st;
    A_1(:,k)=Hom_parameter(C1,A,G,E) % Compute homogenized parameter
end
%% SOLVERs with M_1 preconditionig
c_0 = c_000;
for k=1:1 
    E=E_0(:,k);
    tic;
    [C2,st,norm_evol2]=CG_solver(A,G_m,c_0,E,steps,toler,M_m);
    %[C22,st,norm_evol2]=CGP_solver_left(A,G,c_0,E,steps,toler,M_f_const);
    T2(k,counter) = toc;
    S2(k,counter) = st;
    A_2(:,k)=Hom_parameter(C2,A,G,E) % Compute homogenized parameter
end
%% SOLVER with M_2 preconditionig
c_0=c_000;
for k=1:1
    E=E_0(:,k);
    tic;
    [C3,st,norm_evol3]=CGP_solver_1f_left(A,G,c_0,E,steps,toler,U1,U2); % with better preconditioning
    %[C33,st]=CGP_solver_1f_v3(A,G_n,c_0,E,steps,toler,GG0,GG1,U1,U2);
    T3(k,counter)=toc;
    S3(k,counter) = st;
    A_3(:,k)=Hom_parameter(C3,A,G,E)% Compute homogenized parameter

end
%% error 
  NoP(counter)=N_1*N_2;
  counter=counter+1;
end
%% 
save('experiment_data/exp2/NoP.mat','NoP');
save('experiment_data/exp2/S1.mat','S1');
save('experiment_data/exp2/S2.mat','S2');
save('experiment_data/exp2/S3.mat','S3');
save('experiment_data/exp2/T1.mat','T1');
save('experiment_data/exp2/T2.mat','T2');
save('experiment_data/exp2/T3.mat','T3'); 

NoS1=(1:S1(1,end));
NoS2=(1:S2(1,end));
NoS3=(1:S3(1,end));

save('experiment_data/expPAMM/NoP.mat','NoP'); 

save('experiment_data/expPAMM/NoS1.mat','NoS1'); 
save('experiment_data/expPAMM/NoS2.mat','NoS2');
save('experiment_data/expPAMM/NoS3.mat','NoS3');


save('experiment_data/expPAMM/S1.mat','S1');
 save('experiment_data/expPAMM/S2.mat','S2');
 save('experiment_data/expPAMM/S3.mat','S3');
save('experiment_data/expPAMM/T1.mat','norm_evol1');
save('experiment_data/expPAMM/T2.mat','norm_evol2');
save('experiment_data/expPAMM/T3.mat','norm_evol3'); 

%% Plot
 figure 
 hold on
 plot((1:S1(1,end)) ,norm_evol1)
 plot((1:S2(1,end)),norm_evol2)
 plot((1:S3(1,end)),norm_evol3)
set(gca, 'XScale', 'log', 'YScale', 'log');
legend('G_n','G_m','M_2')

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
 
 
