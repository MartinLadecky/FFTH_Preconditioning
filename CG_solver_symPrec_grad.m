function [c_1,st,norm_evol,norm_sol] = CG_solver_symPrec_grad(A,G,c_0,E,steps,toler,M_m)
%% input
% A   [N_bf2,N_bf1,2,2]- matrix of material parameters in every point of grid
% G_i [N_bf2,N_bf1,2]   -matrix of coeficients of 1st derivative with incorporated preconditioner M.^(-1/2) from both sides
% M   [N_bf2,N_bf1,2]   -preconditioning matrix M.^(-1/2)
% c_0 [N_bf2,N_bf1]     -initial solution
% E   [2,1]             -vector [0;1] or [1;0]
% toler [1]             -relative tolerance
% steps [1]             -max number of steps
%% Output
% c_0 [N_bf,1]          -vector of solution
% st  [1]               -number of stepts
%% 
% toler           % relative tolerance
% steps           % -max number of steps
norm_evol=0;

u_0=c_0./M_m;
M_0 = LHS_freq_symP(A,u_0,G,M_m);

b_0 = M_m.*RHS_freq(A,E,G);

r_0 = b_0-M_0; % 

z_0=r_0./M_m;
grad_z_0=G.*(z_0./(M_m.^-2));
nr0 =sqrt(scalar_product_grad(grad_z_0,grad_z_0));

%nr0 = sqrt(sum(sum(abs(z_0.*z_0))));
%nr0 = sqrt(sum(sum(abs(r_0.*(M_m.^-2).*r_0))));
p_0 = r_0;
for st = 1:steps
    M_1 = LHS_freq_symP(A,p_0,G,M_m);
    
    r_0r_0= sum(sum((r_0.')'.*r_0));
    p_0Ap_0=sum(sum((p_0.')'.*M_1));
    alfa_0 =r_0r_0/p_0Ap_0;
    
    u_1 = u_0 + alfa_0*p_0;
    x_1=u_1.*M_m;
    grad_x_1=G.*(x_1);
    norm_sol(st)=sqrt(scalar_product_grad(grad_x_1,grad_x_1));
    r_1 = r_0 - alfa_0*M_1;
    
    z_1=r_1./(M_m);
    grad_z_1=G.*(z_1./(M_m.^-2));
    nr1 =sqrt(scalar_product_grad(grad_z_1,grad_z_1));
    norm_evol(st)=nr1/nr0;
   %norm_evol(st)=sqrt(sum(sum(abs(z_1.*z_1))))/nr0;%norm(r_1.*(M_m.^-2).*r_1,'fro')/nr0;
    
    if (norm_evol(st)<toler)
         %c_1 = u_1; 
         break; 
    end
    
    beta_0 =sum(sum((r_1.')'.*r_1))/sum(sum((r_0.')'.*r_0));
    p_1 = r_1 + beta_0*p_0;
    %p_1.*(M_m);
    %%
    p_0 = p_1;
    r_0 = r_1;
    u_0 = u_1;
%    norm_evol(st)=norm(c_1,'fro');
end
c_1 =u_1.*M_m;
end