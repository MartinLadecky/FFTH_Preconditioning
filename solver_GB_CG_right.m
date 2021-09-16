function [c_1,st,norm_evol_rr, norm_evol_rz,  norm_sol,e_norm_error]...
    = solver_GB_CG_right(A,G,c_0,E,M,C_ref,steps,toler)
    %% input
    
    grad_c_0=G.*c_0; % Gradient [N_bf2,N_bf1,2]
    c_0=grad_c_0;

    norm_sol(1)=sqrt(scalar_product_grad_energy(c_0,c_0,A));

    M_0 = Projection_Cref(A,grad_c_0,G,M,C_ref); 
    b_0 = Projection_rhs_Cref(A,E,G,M,C_ref) ;
    r_0 = b_0-M_0;
   
    
    nr0=sqrt(scalar_product_grad(r_0,r_0));
    norm_evol_rr(1)=nr0/nr0;
    
    nz0r0 =sqrt(scalar_product_grad_energy(r_0,r_0,C_ref));
    norm_evol_rz(1)=nz0r0/nz0r0;
    
    
    p_0 = r_0;
    
    for st = 1:steps
        Ap_0 = Projection_Cref(A,p_0,G,M,C_ref);

        
        r_0r_0=scalar_product_grad(r_0,r_0);
        p_0Ap_0=scalar_product_grad(p_0,Ap_0);
        alfa_0 = r_0r_0/p_0Ap_0;
        
        c_1 = c_0 + alfa_0.*p_0; 
        norm_sol(st+1)=sqrt(scalar_product_grad_energy(c_1,c_1,A));

        r_1 = r_0 - alfa_0*Ap_0 ;
        
        nr1=sqrt(scalar_product_grad(r_1,r_1));
        norm_evol_rr(st+1)=nr1/nr0;   
        

        nz1r1=sqrt(scalar_product_grad_energy(r_1,r_1,C_ref));
        norm_evol_rz(st+1)=nz1r1/nz0r0;
        
            if ( norm_evol_rr(st+1)<toler) 
                break; 
            end  
        
    
        r_1r_1=scalar_product_grad(r_1,r_1);
        beta_1 =  r_1r_1/r_0r_0;
        p_1 = r_1 + beta_1*p_0;
        %% reiniticialisation
        p_0 = p_1;
        r_0 = r_1;
        
        c_0 = c_1;
        
    end
end

