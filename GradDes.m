function GD = GradDes(q_start, q_end)
%Gradient Descent Constants
eps=0.01; %region of convergence
alpha=0.01; %step size

q_next=q_start; 
count=0; %while loop counter

while norm(q_next - q_end) > eps   
  T_q=(Torque(q_next,q_end))';
  q_next=q_next+ alpha*(T_q/norm(T_q)); 
    
  scatter3( q_next(1),  q_next(2),  q_next(3) );
  hold on
  axis equal
  count=count+1;
  
  Percentage= 100*(abs(q_next(1)) +  abs(q_next(2)) +  abs(q_next(3)))/abs((q_end(1)) +  abs(q_end(2)) +  abs(q_end(3))); 
  disp('Step Number');
  disp(count);
  disp('Percentage Done (%)(value sometimes questionable lol)');
  disp(round(Percentage, 2,'decimal'));
  
end

 disp('Step Number');
 disp(count+1);
 disp('Percentage Goal Reached (%)');
 disp(100);
 disp('Gradient Descent Plot Complete');
 
 PlotObs(); %Plotting Obstacles
 
end