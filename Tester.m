%Test
eta=[1 1 1];

O_w= [1.5; -1; 1.25]; %Wall Obstacle 
rho_w=0.5;

O_m= [2; 0; 2.25]; %Mount Obstacle 
rho_m=0.25;

O_c= [2; 0; 1.25]; %Cylinder Obstacle 
rho_c=0.75;

q_1=[0 0 0];
O_1=O_positions(q_1);
Fep_w=F_rep(O_1, O_w, rho_w,eta);
Fep_m=F_rep(O_1, O_m, rho_m,eta);
Fep_c=F_rep(O_1, O_c, rho_c,eta)

for i=1:3   
Total{i}=Fep_w{i}+Fep_m{i}+Fep_c{i}
end
