function [c_1,st,norm_evol_rr,norm_evol_energy,norm_evol_grad, estim,norm_sol] = solver_DB_PCG(A,G,c_0,E,steps,toler,M_f,tau)
    %% input
    norm_sol(1)=sqrt(scalar_product_grad_energy(G.*c_0,G.*c_0,A));
    M_0 = LHS_freq(A,c_0,G); % System matrix * initial solution(x_0=c_0)
    b_0 = RHS_freq(A,E,G);  % Right hand side vector

    r_0 = b_0-M_0; % x_0=0
    
    nr0 =norm(r_0,'fro');
    norm_evol_rr(1)=nr0/nr0;
    
    nDr0Dr0=sqrt(scalar_product_grad(G.*r_0,G.*r_0));
    norm_evol_grad(1)=nDr0Dr0/nDr0Dr0;

    z_0 = r_0./M_f; % solve lin system rM_0=M_f^(-1)*r_0: rM_0 is idagonal matrix


    nz0r0 = sqrt(scalar_product(r_0,z_0));
    norm_evol_energy(1)=nz0r0/nz0r0;
    

    p_0 = z_0;

    k=1;
    d=0;
    estim=0;
    for st = 1:steps
        Ap_0 = LHS_freq(A,p_0,G);
     
        z_0r_0= sum(sum((z_0.')'.*r_0));
        p_0Ap_0=sum(sum((p_0.')'.*Ap_0));
        alfa_0 = z_0r_0/p_0Ap_0;
        

        c_1 = c_0 + alfa_0.*p_0;
        
        norm_sol(st+1)=sqrt(scalar_product_grad_energy(G.*c_1,G.*c_1,A));
        r_1 = r_0-alfa_0*Ap_0;
        

         z_1=r_1./M_f;
         z_1r_1=real(sum(sum((z_1.')'.*r_1)));
       	%nr1=sqrt(sum(sum(abs(r_1.*z_1))));
        nr1 =norm(r_1,'fro');
        norm_evol_rr(st+1)=nr1/nr0;
        nDr1Dr1=sqrt(scalar_product_grad(G.*r_1,G.*r_1));
        norm_evol_grad(st+1)=nDr1Dr1/nDr0Dr0;

        norm_evol_energy(st+1)=sqrt(z_1r_1)/nz0r0;
            if ( norm_evol_rr(st+1)<toler)
                %c_1 = c_0; 
                break; 
            end  
        
        
        G.*z_1;
        
        
        beta_1 =z_1r_1 /z_0r_0;
        p_1 = z_1 + beta_1*p_0;
        
        Delta(st)=real(alfa_0*z_1r_1);
        curve(st)=0;
        curve=curve+Delta(st);
        if st >1
            S=findS(curve,Delta,k);
            num = S*Delta(st);
            den = sum(Delta(k:st-1));
            
            while (d>= 0) && (num/den<= tau)
                delay(k)=d;
                estim(k)=den;
                k=k+1;
                d=d-1;
                den=sum(Delta(k:st-1));
            end
            d=d+1;
        end
        
        
        %% 
        p_0 = p_1;
        r_0 = r_1;
        z_0 = z_1; 
        c_0 = c_1;
    end
end


function [S]=findS(curve,Delta,k)
    ind=find((curve(k)./curve) <= 1e-4,1,'last');
    if isempty(ind)
       ind = 1 ;
    end
    S = max(curve(ind:end-1)./Delta(ind:end-1));
end








