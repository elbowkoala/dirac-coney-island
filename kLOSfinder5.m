
function [ cone_ksym, k_offset ] = kLOSfinder5(cone, bin_E, bin_k)
    
    cone_b = Binning_2d(cone, bin_E, bin_k);
    cone_f = imgaussfilt(cone_b,1);
    
    flip_dots = [];
    NN = floor(2/3*size(cone_f,1));
    for nn = 1:(size(cone_f,1) - NN)
        cone_s = cone_f(nn:nn+NN,:);  
        cone_s = cone_s ./ sum(cone_s(:));
        cone_n = cone_s;
        %cone_n = (cone_s - min(cone_s(:)))/(max(cone_s(:))-min(cone_s(:)));
        flip_dots(nn) = sum(dot( cone_n, flipud(cone_n)));
    end
    LOS = find(flip_dots==max(flip_dots));
    cone_ksym = cone_b(LOS:LOS+NN,:);
    k_offset = LOS-1;
    
end
