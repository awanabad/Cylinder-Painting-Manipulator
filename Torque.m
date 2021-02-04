function Tao=Torque(q_start,q_end)

d_at=0.3; %switching distance (Conic to Parabolic)
zeta=[0.5 1.5 2.1]; %Attractive Force Tuning Const
eta=[1.5 0.7 0.1]; %Repulsive Force Tuning Const

%Obstacles (wall, mount and gas cylinder)
O_w= [1.5; -1; 1.25]; %Wall Obstacle 
rho_w=0.8;
O_m= [2; 0; 2.25]; %Mount Obstacle 
rho_m=0.25;
O_c= [2; 0; 1.25]; %Cylinder Obstacle 
rho_c=0.55;

%Starting and Ending Origins 
O_start=O_positions(q_start);
O_end=O_positions(q_end);
%Determining Jacobian 
J_v=J_matrix(q_start); %(Pre-transposed in Function)

%Calculating Attractive Force
Fat=F_att(O_start, O_end, d_at, zeta);

%Calculating and Summing Repulsive Force
Frep_w=F_rep(O_start, O_w, rho_w,eta);
Frep_m=F_rep(O_start, O_m, rho_m,eta);
Frep_c=F_rep(O_start, O_c, rho_c,eta);

Frep={1:3, rand(3,1)}; %memory preallocation

for i=1:3   
    Frep{i}=Frep_w{i}+Frep_m{i}+Frep_c{i};
end

%Torque Calculation
Tao_temp={1:3, rand(3,1)}; %memory preallocation

for j=1:3
    if Fat{j}==0 %Since cell element shows single 0 for 0 matrix
        Fat{j} = [0; 0; 0];
   end
    if Frep{j}==0
        Frep{j} = [0; 0; 0];
    end

    Tao_rep=(J_v{j}*Frep{j});
    Tao_att=(J_v{j}*Fat{j});
    Tao_temp{j}=Tao_rep+Tao_att;
end

Tao=Tao_temp{1}+Tao_temp{2}+Tao_temp{3};

end
 