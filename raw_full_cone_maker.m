raw_full_cone = zeros(size(cones(:,:,1)));

for i = 1:num_scans
    raw_full_cone = raw_full_cone + cones(:,:,i);
end
