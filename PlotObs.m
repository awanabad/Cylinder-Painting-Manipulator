function F= PlotObs()
 
%Obstacles (wall, mount and gas cylinder)
O_w= [1.5; -1; 1.25]; %Wall Obstacle 
rho_w=0.8;
O_m= [2; 0; 2.25]; %Mount Obstacle 
rho_m=0.25;
O_c= [2; 0; 1.25]; %Cylinder Obstacle 
rho_c=0.55;

Rho=[rho_w, rho_c, rho_m];
Origin={O_w, O_c, O_m};


[A,B,C] = sphere;         %plot obstacles as spheres
    
for i=1:numel(Origin)
    
surf(Rho(i)*A+Origin{i}(1), Rho(i)*B+Origin{i}(2), Rho(i)*C+Origin{i}(3));   
hold on 
axis equal    
end