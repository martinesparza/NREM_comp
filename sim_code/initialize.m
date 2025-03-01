function [ C ] = initialize( text,duration_vec,amplitude_vec,max_sim,dt)
% Funcion para initializar un celula de memoria
rows = length(duration_vec);
z_axis = length(amplitude_vec);
C = cell(rows,2);

for i = duration_vec
    iterations = i/dt;
    margin = 10000;
    loc = find(duration_vec == i);
    C{loc,1} = strcat(text,{' '},num2str(i));
    C{loc,2} = zeros(max_sim, iterations + (2*margin) + 1,z_axis);

end
end

