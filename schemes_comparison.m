
%% BEGIN
clc;clear;




counter=1;
steps = 200;
for loop=[7]
N_1=2*(loop^2)+1

N_2=N_1; % number of points in x_2

%% Mesh parameters
h_1=2*pi/(N_1); % step in x_2
h_2=2*pi/(N_2); % step in x_2

%% Coordinates
x=zeros(N_2,N_1,2);
[x(:,:,1),x(:,:,2)]=meshgrid(-pi+h_1/2:h_1:pi-h_1/2,-pi+h_2/2:h_2:pi-h_2/2);

%% Material coeficient matrix

Pixels = imread('structure_4.png');

pixa=round(linspace(1,size(Pixels,2),N_1));
piya=round(linspace(1,size(Pixels,1),N_2));
par=1000
  

C=zeros(N_2,N_1,2,2);
    for i=1:N_2
        for j=1:N_1    
            C(i,j,:,:)=a_matrix_img_aniso(Pixels(piya(i),pixa(j)),par);
          % C(i,j,:,:)=a_matrix(x(i,j,:));
          % C(i,j,:,:)=a_matrix_img_aniso(Pixels(piya(i),pixa(j)),par);          
        end       
    end

 

%% Material ananlysis
d=[mean(mean(C(:,:,1,1))) mean(mean(C(:,:,1,2)));...
   mean(mean(C(:,:,2,1))) mean(mean(C(:,:,2,2)))];
% d=[1 0;
%    0 1];

C_ref=zeros(N_2,N_1,2,2); 
C_ref_inv=zeros(N_2,N_1,2,2); 
    for i=1:N_2
        for j=1:N_1    
            C_ref(i,j,:,:)=d;
            C_ref_inv(i,j,:,:)=d^-1;
        end
    end

%% Gradient operator
G = G_clasic(N_1,N_2);

%% Preconditioner

M=(d(1,1).*(G(:,:,1).^2)+d(2,2).*(G(:,:,2).^2)...
                   +2*d(1,2).*(G(:,:,1).*(G(:,:,2))));  
               
M((end+1)/2,(end+1)/2)=1;

%% Conditions:
tau=0.25;
toler = 1e-6;

c_0=zeros(N_2,N_1);
E_0=[1 0;0 1];
A_0=zeros(2,2);

E=E_0(:,1);
%% Gradient-Based solvers

disp('Gradient-Based Conjugage Gradient solver')

   [C_GB_CG,st_GB_CG, norm_evol_GB_CG_grad, norm_evol_GB_CG_energy, estim_GB_CG,  sol_norm_GB_CG]...
       =solver_GB_CG(C,G,c_0,E,M,d,steps,toler,tau);

    S_GB_CG(1,counter) = st_GB_CG+1;
    
    A_(:,1)=Hom_parameter_grad(C_GB_CG,C,G,E) % Compute homogenized parameter
    A_GB_CG(counter)=A_(1,1);


disp('Modifield Gradient-Based Conjugage Gradient solver')


    [C_GB_CG_mod,st_GB_CG_mod,norm_evol_GB_CG_mod_grad,norm_evol_GB_CG_mod_energy, estim_GB_CG_mod,  sol_norm_GB_CG_mod]...
        =solver_GB_CG_modif(C,G,c_0,E,M,d,steps,toler,tau);

    S_GB_CG_mod(1,counter) =st_GB_CG_mod+1;

    A_(:,1)=Hom_parameter_grad(C_GB_CG_mod,C,G,E) % Compute homogenized parameter
    A_GB_CG_mod(counter)=A_(1,1);


disp('Gradient-Based Preconditioned Conjugage Gradient solver')


    [C_GB_PCG,st_GB_PCG, norm_evol_GB_PCG_rr, norm_evol_GB_PCG_energy, estim_GB_PCG,  sol_norm_GB_PCG]...
     =solver_GB_PCG(C,G,c_0,E,M,steps,toler,tau);


    S_GB_PCG(1,counter) = st_GB_PCG+1;

    A_(:,1)=Hom_parameter_grad(C_GB_PCG,C,G,E) % Compute homogenized parameter
    A_GB_PCG(counter)=A_(1,1);


%% Displacement-Based Preconditioned Conjugage Gradient solver
disp('Displacement-Based Preconditioned Conjugage Gradient solver')

toler = 1e-10;

    [C_DB_PCG,st_DB_PCG,norm_evol_DB_PCG_rr,norm_evol_DB_PCG_energy, norm_evol_DB_PCG_grad, estim_DB_PCG,  sol_norm_DB_PCG]...
        =solver_DB_PCG(C,G,c_0,E,M,steps,toler,tau);% with preconditioning

    S_DB_PCG(1,counter) = st_DB_PCG+1;
    
    A_(:,1)=Hom_parameter(C_DB_PCG,C,G,E)% Compute homogenized parameter
    A_DB_PCG(counter)=A_(1,1);  


NoP(counter)=N_1*N_2;
counter=counter+1;
end

% % Plot estimates
 rel_estim_DB_PCG=(estim_DB_PCG./estim_DB_PCG(1)).^(1/2);
 rel_estim_GB_CG=estim_GB_CG./estim_GB_CG(1);


%% Plot residuals
figure 
hold on

    plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_rr,'--xb')
    plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_energy,'--ob')

    plot((1:S_GB_CG(1,end)),norm_evol_GB_CG_grad,'-.xr')
    plot((1:S_GB_CG(1,end)),norm_evol_GB_CG_energy,'-.or')

    plot((1:S_GB_CG_mod(1,end)),norm_evol_GB_CG_mod_grad,'-.xg')
    plot((1:S_GB_CG_mod(1,end)),norm_evol_GB_CG_mod_energy,'-.og')

    plot((1:S_GB_PCG(1,end)),norm_evol_GB_PCG_rr,'-.xk')
    plot((1:S_GB_PCG(1,end)),norm_evol_GB_PCG_energy,'-.ok')

set(gca, 'XScale', 'linear', 'YScale', 'log');
legend( 'DB PCG || r ||',       'DB PCG || r ||_M', ...
        'GB CG || r^* ||',      'GB CG || r^* ||_M',...
        'GB CG mod || r^{m} ||','GB CG mod || r ||_M',...
        'GB PCG || r^{@} ||',   'GB PCG || r ||_M') %'DB PCG || Dr ||',
title('Residuals ') 
%% Plot residuals with estimates 
figure 
hold on

    plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_rr,'--xb')
    plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_energy,'--ob')

    plot((1:S_GB_CG(1,end)),norm_evol_GB_CG_grad,'-.xr')
    plot((1:S_GB_CG(1,end)),norm_evol_GB_CG_energy,'-.or')

    plot((1:S_GB_CG_mod(1,end)),norm_evol_GB_CG_mod_grad,'-.xg')
    plot((1:S_GB_CG_mod(1,end)),norm_evol_GB_CG_mod_energy,'-.og')

    plot((1:S_GB_PCG(1,end)),norm_evol_GB_PCG_rr,'-.xk')
    plot((1:S_GB_PCG(1,end)),norm_evol_GB_PCG_energy,'-.ok')

    plot((1:numel(estim_DB_PCG)),abs(rel_estim_DB_PCG),'.')
    plot((1:numel(estim_DB_PCG)),abs(rel_estim_DB_PCG./(1-tau)),'.')
 
set(gca, 'XScale', 'linear', 'YScale', 'log');
legend( 'DB PCG || r ||',       'DB PCG || r ||_M', ...
        'GB CG || r^* ||',      'GB CG || r^* ||_M',...
        'GB CG mod || r^{m} ||','GB CG mod || r ||_M',...
        'GB PCG || r^{@} ||',   'GB PCG || r ||_M',...
        'DB estim || e_k ||_K',' DB estim || e_k ||_K') %'DB PCG || Dr ||',
title('Residuals with  erorr K-norm estimates') 


%% Plot solution norm
figure 
hold on

    plot((1:S_DB_PCG(1,end)),abs(sol_norm_DB_PCG(end)-sol_norm_DB_PCG),'-.^')
    plot((1:S_GB_CG(1,end)),abs(sol_norm_DB_PCG(end)-sol_norm_GB_CG),'-.*')

    plot((1:S_GB_CG_mod(1,end)),abs(sol_norm_DB_PCG(end)-sol_norm_GB_CG_mod),'-.o')
  
    plot((1:S_GB_PCG(1,end)),abs(sol_norm_DB_PCG(end)-sol_norm_GB_PCG),'-.^')
 
set(gca, 'XScale', 'linear', 'YScale', 'log');
legend('DB PCG || Du ||_C','GB CG || Du ||_C','GB mCG || Du ||_C','GB PCG || Du ||_C')
title('solution norm: relative')

%% Plot solution norm
figure 
hold on

    plot((1:S_DB_PCG(1,end)),sol_norm_DB_PCG,'-.^')
    plot((1:S_GB_CG(1,end)),sol_norm_GB_CG,'-.*')

    plot((1:S_GB_CG_mod(1,end)),sol_norm_GB_CG_mod,'-.o')

    plot((1:S_GB_PCG(1,end)),sol_norm_GB_PCG,'-.^')
 
set(gca, 'XScale', 'linear', 'YScale', 'log');
legend('DB PCG ','GB CG ','GB CG modif','GB PCG ')
title('solution norm: absolut')

%% Plot hom_mat prop
% A_refer=real(A3);%9.403399951113226;
% figure 
%  hold on
% 
%  plot(NoP,abs(A_refer-real(A3))','^')
%   plot(NoP,abs(A_refer-real(Ag))','-.^') 
%   
%  plot(NoP,abs(A_refer-real(Ap))','-.*')
%  plot(NoP,abs(A_refer-real(Apc))','-.>') 
% 
%  set(gca, 'XScale', 'linear', 'YScale', 'log');
% 
% legend('DB','DB modif','GB orig','GB modif')
%  title(' Hom parameter error')
% figure 
%  hold on
%  plot(NoP,imag(A1)','x')
%  plot(NoP,imag(A2)','o')
%  plot(NoP,imag(A3)','^')
%  plot(NoP,imag(Ap)','-.*')
%  plot(NoP,imag(Ag)','-.^') 
% legend('in G','Symetric','Left','GB','grad norm')


close all

 %% Plot A function