function RepForce = F_rep(O_i, O_object, rho_0,eta)

%Force Calculation
RepForce={1:3, rand(3,1)}; %memory preallocation

for i= 1:3
   
    rho_q= norm((O_i{i})-O_object);
    
    if rho_q <= rho_0 %radius of influence check
        grad_rho_q= ((O_i{i}) - O_object)/rho_q ;
        RepForce{i}= eta(i) * ( (1/rho_q) - (1/rho_0) ) *  (1/(rho_q^2)) *grad_rho_q;
        
    elseif rho_q > rho_0 
        
        RepForce{i}=0;
        
    end   
end
end