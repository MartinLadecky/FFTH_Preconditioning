function [c_1,st] = CGP_solver_constant(A,G_n,c_0,E,steps,toler,M_f)
%% input
% A   [N_bf2,N_bf1,2,2] -matrix of material parameters in every point of grid
% G_n [N_bf2,N_bf1,2]   -matrix of coeficients of 1st derivative
% c_0 [N_bf2,N_bf1]     -initial solution
% E   [2,1]             -vector [0;1] or [1;0]
% toler [1]             -relative tolerance
% steps [1]             -max number of steps
%% Output
% c_0 [N_bf2,N_bf1]     -vector of solution
% st  [1]               -number of stepts
%% 
% toler           % relative tolerance
% steps           % -max number of steps
c_1 = c_0;
M_0 = LHS_freq(A,c_1,G_n); % System matrix * initial solution(x_0=c_0)
b_0 = RHS_freq(A,E,G_n);     % Right hand side vector

r_0 = b_0-M_0; % x_0=0
r_0 = r_0./M_f; % solve lin system rM_0=M_f^(-1)*r_0: rM_0 is idagonal matrix
nr0 = norm(r_0,'fro');
if (norm(r_0,'fro')/nr0<toler)
    st = 0;
   return;
end
p_0 = r_0;
for st = 1:steps
    M_1 = LHS_freq(A,p_0,G_n) ;
    M_1 = M_1./M_f;
    alfa_0 = real(sum(sum((r_0.')'.*r_0))/sum(sum((p_0.')'.*M_1)));
    x_1 = c_1 + alfa_0.*p_0;
    r_1 = r_0-alfa_0*M_1;
    if (norm(r_1,'fro')/nr0<toler)
         c_1 = x_1; 
        break; 
    end    
    beta_0 = real(sum(sum((r_1.')'.*r_1))/sum(sum((r_0.')'.*r_0)));
    p_1 = r_1 + beta_0*p_0;    
    %% 
    p_0 = p_1;
    r_0 = r_1;
    c_1 = x_1;     
end
%x_0=fftshift(ifft2(ifftshift(c_0)));
end