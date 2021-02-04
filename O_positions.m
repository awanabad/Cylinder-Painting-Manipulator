function O = O_positions(q)
%T Matrices
T=T_matrix(q);

%Origins
O={1:3, rand(1,3)}; %memory preallocation
for j=1:numel(T)
    O_temp=[T{j}(1,4), T{j}(2,4) T{j}(3,4)];
    O{j}=O_temp';
end

end
