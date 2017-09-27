tic; 

%rfc_big_Es = zeros(1,num_scans);
%rfc_big_ks = zeros(1,num_scans);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_E = 5;
bin_k = 2;

rfc_sigma = 5;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Making the 0th-Order Template %%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_full_cone = zeros(size(cone_range_K,2),size(cone_range_E,2));

for i = 1:1:num_scans;
    raw_full_cone = raw_full_cone + cones(:,:,i);
end

rfc_b = Binning_2d(raw_full_cone, bin_E, bin_k);
rfc_big_window_Kpix = (50:130);
rfc_big_window_Epix = (40:130);
rfc_big_window = rfc_b(rfc_big_window_Kpix, rfc_big_window_Epix); 
rfc_big_bwg = mat2gray(rfc_big_window);
%rfc_bg = mat2gray(rfc_bg(11:end-10,31:77));
window_big_k_coor = size(rfc_big_bwg,1)/2;
window_big_E_coor = 35; %%%%%%%%%%%%%%% hand-pick point from rfc_bg%%%%%%%%%%%

disp('Template step completed '), toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Round 1 Scan with Large Window Including Bulk Bands %%%%%%%%%%%%%

scan_maxes = cat(2,[]);
figure,
NNN=0;
for i = 12% [1,round(961*rand),round(961*rand),round(961*rand)]%1:1:num_scans;  
   % NNN=NNN+1;
    if rem(i,100) == 0
        disp(['Now on scan ',num2str(i)])
    end
    
    scone = cones(:,:,i);
    scone = imgaussfilt(scone, rfc_sigma);
    scone_b = Binning_2d(scone, bin_E, bin_k);
    scone_bg = mat2gray(scone_b);

    scann = normxcorr2( rfc_big_bwg, scone_bg );
    scan = scann(:,71:end);
    [scan_max, scan_max_i] = max(scan(:));
    scan_maxes(i) = scan_max;
    
    [y_peak, x_peak] = find(scan==max(scan(:)));
    y_offset = y_peak - size(rfc_big_bwg,1);
    x_offset = x_peak - size(rfc_big_bwg,2) + 70;
    
    %plot([1,size(rfc_bwg,2),size(rfc_bwg,2),1,1]+x_offset,[1,1,size(rfc_bwg,1),size(rfc_bwg,1),1]+y_offset,'r'), hold off;
    
    rfc_big_Es(i) = (x_offset + window_big_E_coor -1)*bin_E + bin_E/2; %%%%%%!!!!%%
    rfc_big_ks(i) = (y_offset + window_big_k_coor  -1)*bin_k + bin_k/2;    
    
    scan_Eaxis = x_peak - 30 : x_peak + 30;
    scan_Kaxis = y_peak - 50: y_peak + 50;
    
    subplot(1,2,1), %imagesc(scan_Kaxis,scan_Eaxis,rot90(scan,-1)), axis xy, title(num2str(i)), hold on;
    imagesc(rot90(scan,-1)), axis xy, hold on;
    plot(size(scan,1)-y_peak,x_peak,'b*'), hold off;
    subplot(1,2,2), imagesc(rot90(scone_bg,-1)), axis xy, hold on;
    plot(flipud([1,size(rfc_bwg,1),size(rfc_bwg,1),1,1]+size(scone_bg,1)-size(rfc_bwg,1)-y_offset),[1,1,size(rfc_bwg,2),size(rfc_bwg,2),1]+x_offset,'r'), hold off;
    
end
disp('Scanning Round 1 completed'), toc


E_map = reshape(rfc_big_Es, 31, 31);
k_map = reshape(rfc_big_ks, 31, 31);

EE_sd = mean([std(E_map(:,2:end)-E_map(:,1:end-1)),std(E_map(2:end,:)-E_map(1:end-1,:))])
kk_sd = mean([std(k_map(:,2:end)-k_map(:,1:end-1)),std(k_map(2:end,:)-k_map(1:end-1,:))])

%{
figure, subplot(2,2,[1 2]), imagesc(rfc_bwg), axis xy, title('window1');
subplot(2,2,3), imagesc(E_map), axis xy, title('Emap1');
subplot(2,2,4), imagesc(k_map), axis xy, title('kmap1');

figure, subplot(1,2,1), imagesc(E_map), axis xy, title('Emap');
subplot(1,2,2), imagesc(k_map), axis xy, title('kmap');
%}
