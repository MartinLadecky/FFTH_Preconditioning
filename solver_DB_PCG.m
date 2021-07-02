function [c_1,st,norm_evol_rr,norm_evol_energy,norm_evol_grad, norm_sol]...
    = solver_DB_PCG(A,G,c_0,E,M,steps,toler)
    %% input
    norm_sol(1)=sqrt(scalar_product_grad_energy(G.*c_0,G.*c_0,A));
    M_0 = LHS_freq(A,c_0,G); 
    b_0 = RHS_freq(A,E,G);  

    r_0 = b_0-M_0; 
    
    nr0 =sqrt(scalar_product(r_0,r_0));
    norm_evol_rr(1)=nr0/nr0;
    
    nDr0Dr0=sqrt(scalar_product_grad(G.*r_0,G.*r_0));
    norm_evol_grad(1)=nDr0Dr0/nDr0Dr0;

    z_0 = r_0./M; 
    
    nz0r0 = sqrt(scalar_product(r_0,z_0));
    norm_evol_energy(1)=nz0r0/nz0r0;
    
    p_0 = z_0;

    for st = 1:steps
        Ap_0 = LHS_freq(A,p_0,G);
     
        z_0r_0= scalar_product(z_0,r_0 );
        p_0Ap_0=scalar_product(p_0,Ap_0);
        alfa_0 = z_0r_0/p_0Ap_0;
        
        c_1 = c_0 + alfa_0.*p_0;
        
        norm_sol(st+1)=sqrt(scalar_product_grad_energy(G.*c_1,G.*c_1,A));
        r_1 = r_0-alfa_0*Ap_0;
        
        nr1 =sqrt(scalar_product(r_1,r_1 ));
        norm_evol_rr(st+1)=nr1/nr0;

        
        z_1=r_1./M;
        
        z_1r_1=scalar_product(z_1,r_1 );
        norm_evol_energy(st+1)=sqrt(z_1r_1)/nz0r0;
        
        
        nDr1Dr1=sqrt(scalar_product_grad(G.*r_1,G.*r_1));
        norm_evol_grad(st+1)=nDr1Dr1/nDr0Dr0;

            if ( norm_evol_rr(st+1)<toler)
                break; 
            end  

        beta_1 =z_1r_1 /z_0r_0;
        p_1 = z_1 + beta_1*p_0;
        %% reiniticialisation
       
        p_0 = p_1;
        r_0 = r_1;
        z_0 = z_1; 
        c_0 = c_1;

    end
end






