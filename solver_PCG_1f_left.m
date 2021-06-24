function [c_1,st,norm_evol] = solver_PCG_1f_left(A,G,c_0,E,steps,toler,U1,U2)
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
% steps           % max number of steps

[n2,n1] = size(c_0);

M_0 = LHS_freq(A,c_0,G); % System matrix * initial solution(x_0=c_0)
b_0 = RHS_freq(A,E,G);% Right hand side vector


r_0 = b_0-M_0; % x_0=0

nr0 = norm(r_0,'fro');

%%  NEW method
r_pom=r_0;
for s = 1:n1-1
        pom = conj(U2(:,s))./U1(:,s);    %delenie poddiag diag celnom, pre vsetky bloky naraz    
        r_pom(:,s+1) = r_pom(:,s+1) - r_pom(:,s).*pom; 
end
    r_pom= r_pom./U1;% vydelim diagonalov
for s = n1:-1:2
        pom = U2(:,s-1)./U1(:,s);
        %pom;
        r_pom(:,s-1) = r_pom(:,s-1) - r_pom(:,s).*pom; 
end
    z_0 = r_pom./U1;
  
%% 
p_0 = z_0;

for st = 1:steps
    Ap_0 = LHS_freq(A,p_0,G);     
    
    alfa_0 = (sum(sum((z_0.')'.*r_0))/sum(sum((p_0.')'.*Ap_0)));
    
    c_1 = c_0 + alfa_0.*p_0;
    r_1 = r_0-alfa_0*Ap_0;
    
    norm_evol(st)=norm(r_1,'fro')/nr0;
    if (norm_evol(st)<toler)
         break; 
    end 
    % Solve auxiliary system   
    r_pom=r_1;
    for s = 1:n1-1
        pom = conj(U2(:,s))./U1(:,s);
        r_pom(:,s+1) = r_pom(:,s+1) - r_pom(:,s).*pom; 
    end
    r_pom= r_pom./U1;
    for s = n1:-1:2
        pom = U2(:,s-1)./U1(:,s);
        r_pom(:,s-1) = r_pom(:,s-1) - r_pom(:,s).*pom; 
    end
    z_1 = r_pom./U1;  

    beta_1 = (sum(sum((z_1.')'.*r_1))/sum(sum((z_0.')'.*r_0)));
    p_1 = z_1 + beta_1*p_0;    
    %% 
    p_0 = p_1;
    r_0 = r_1;
    c_0 = c_1;
    z_0 = z_1;
end
%x_0=fftshift(ifft2(ifftshift(c_0)));
end