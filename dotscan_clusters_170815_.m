tic; 

dot_dirac_Es = zeros(1,num_scans);
dot_dirac_ks = zeros(1,num_scans);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_E = 5;
bin_k = 2;

rfc_sigma = 5;

window_E_coor = 19; %%%%%%%%%%%%%%% hand-pick point from rfc_bg%%%%%%%%%%%
window_k_coor = 15.5; %%%%%%%%%%%%% hand-pick point from rfc_bg %%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Making the 0th-Order Template %%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_full_cone = zeros(size(cone_range_K,2),size(cone_range_E,2));

for i = 1:1:num_scans;
    raw_full_cone = raw_full_cone + cones(:,:,i);
end

[rfc_window, window_k_offset] = kLOSfinder5(raw_full_cone, bin_E, bin_k);
rfc_b = rfc_window(11:end-10,55:end-61);
rfc_bg = mat2gray(rfc_b);
%rfc_bg = mat2gray(rfc_bg(11:end-10,31:77));
window_k_coor = size(rfc_bg,1)/2;




disp('Template step completed '), toc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Round 1 Scan with Large Window Including Bulk Bands %%%%%%%%%%%%%

scan_maxes = cat(2,[]);

for i = 1:1:num_scans;  

    scone = cones(:,:,i);
    scone = imgaussfilt(scone, rfc_sigma);
    scone_b = Binning_2d(scone, bin_E, bin_k);
    scone_bg = mat2gray(scone_b);

    scan = normxcorr2( rfc_bg, scone_bg );
    scan = scan(:,51:150);
    [scan_max, scan_max_i] = max(scan(:));
    scan_maxes(i) = scan_max;

    [y_peak, x_peak] = find(scan==max(scan(:)));
    y_offset = y_peak - size(rfc_bg,1);
    x_offset = x_peak - size(rfc_bg,2) + 50;

    dot_dirac_Es(i) = (x_offset + window_E_coor -1)*bin_E + bin_E/2; %%%%%%!!!!%%
    dot_dirac_ks(i) = (y_offset + window_k_coor -1)*bin_k + bin_k/2;    
end
disp('Scanning Round 1 completed'), toc


E_map = reshape(dot_dirac_Es, 31, 31);
k_map = reshape(dot_dirac_ks, 31, 31);

EE_sd = mean([std(E_map(:,2:end)-E_map(:,1:end-1)),std(E_map(2:end,:)-E_map(1:end-1,:))])
kk_sd = mean([std(k_map(:,2:end)-k_map(:,1:end-1)),std(k_map(2:end,:)-k_map(1:end-1,:))])


figure, subplot(2,2,[1 2]), imagesc(rfc_bg), axis xy, title('window1');
subplot(2,2,3), imagesc(E_map), axis xy, title('Emap1');
subplot(2,2,4), imagesc(k_map), axis xy, title('kmap1');

figure, subplot(1,2,1), imagesc(E_map), axis xy, title('Emap');
subplot(1,2,2), imagesc(k_map), axis xy, title('kmap');
%}
