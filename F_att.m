function AttForce = F_att(O_i, O_f, d, zeta)
 

%Force Calculation
AttForce={1:3, rand(3,1)}; %memory preallocation

for i=1:3
    
    if norm(O_i{i}-O_f{i}) <= d
        AttForce{i}= -zeta(i)*(O_i{i}-O_f{i}); %Parabolic
        
    elseif norm(O_i{i}-O_f{i}) > d    
       AttForce{i}= -d*zeta(i)*((O_i{i}-O_f{i}))/(norm(O_i{i}-O_f{i})); %Conic
    end
    
end
end
