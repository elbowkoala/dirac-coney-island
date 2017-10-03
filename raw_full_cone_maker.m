result_comp = zeros(size(cones(:,:,1)));

for i = 1:num_scans
        result_comp = result_comp + result(:,:,i);
end
