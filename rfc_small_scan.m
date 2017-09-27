tic; 

rfc_small_Es = zeros(1,num_scans);
rfc_small_ks = zeros(1,num_scans);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_E = 5;
bin_k = 2;

rfc_sigma = 3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Making the 0th-Order Template %%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_full_cone = zeros(size(cone_range_K,2),size(cone_range_E,2));

for i = 1:1:num_scans;
    raw_full_cone = raw_full_cone + cones(:,:,i);
end

rfc_b = Binning_2d(raw_full_cone, bin_E, bin_k);
rfc_small_window_Kpix = (70:110);
rfc_small_window_Epix = (55:95);
rfc_small_window = rfc_b(rfc_small_window_Kpix, rfc_small_window_Epix); 
rfc_small_bwg = mat2gray(rfc_small_window);

window_small_k_coor = size(rfc_small_bwg,1)/2;
window_small_E_coor = 19; %%%%%%%%%%%%%%% hand-pick point from rfc_bg%%%%%%%%%%%

disp('Template step completed '), toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Round 1 Scan with Large Window Including Bulk Bands %%%%%%%%%%%%%

scan_maxes = cat(2,[]);
%figure,
NNN=0;
for i =  1:num_scans%[round(961*rand),round(961*rand),round(961*rand), round(961*rand)]%1:1:num_scans;  
    NNN=NNN+1;
    if rem(i,100) == 0
        disp(['Now on scan ',num2str(i)])
    end
    
    scone = cones(:,:,i);
    Epix = (rfc_big_Es(i) - bin_E/2)/bin_E;
    kpix = (rfc_big_ks(i) - bin_k/2)/bin_k;
    Epixhr = size(rfc_small_bwg,2)/2 + 15; 
    kpixhr = size(rfc_small_bwg,1)/2 + 10;
    
    scone = imgaussfilt(scone, rfc_sigma);
    scone_b = Binning_2d(scone, bin_E, bin_k);
    scone_b = scone_b( round(kpix-kpixhr):round(kpix+kpixhr), round(Epix-Epixhr):round(Epix+Epixhr));                    
    scone_bg = mat2gray(scone_b);

    scann = normxcorr2( rfc_small_bwg, scone_bg );
    scan=scann(:,size(rfc_small_bwg,2)+1:end-size(rfc_small_bwg,2));
    [scan_max, scan_max_i] = max(scan(:));
    scan_maxes(i) = scan_max;
    
    [y_peak, x_peak] = find(scan==max(scan(:)));
    y_offset = y_peak - size(rfc_small_bwg,1);
    x_offset = x_peak - size(rfc_small_bwg,2) + size(rfc_small_bwg,2);
        
    rfc_small_Es(i) = (x_offset + window_small_E_coor + round(Epix-Epixhr) -1)*bin_E + bin_E/2; %%%%%%!!!!%%
    rfc_small_ks(i) = (y_offset + window_small_k_coor + round(kpix-kpixhr) -1)*bin_k + bin_k/2;    
    %{
    subplot(2,4,NNN), %imagesc(scan_Kaxis,scan_Eaxis,rot90(scan,-1)), axis xy, title(num2str(i)), hold on;
    imagesc(rot90(scan,-1)), axis xy, hold on;
    plot(size(scan,1)-y_peak,x_peak,'b*'), hold off;
    subplot(2,4,NNN+4), imagesc(rot90(scone_bg,-1)), axis xy, hold on;
    %plot([1,size(rfc_small_bwg,2),size(rfc_small_bwg,2),1,1]+x_offset,[1,1,size(rfc_small_bwg,1),size(rfc_small_bwg,1),1]+y_offset,'r'), hold off;
    plot([1,1,size(rfc_small_bwg,1),size(rfc_small_bwg,1),1]+y_offset,[1,size(rfc_small_bwg,2),size(rfc_small_bwg,2),1,1]+x_offset,'r'), hold off;
    %plot(flipud([1,size(rfc_small_bwg,1),size(rfc_small_bwg,1),1,1]+size(scone_bg,1)-size(rfc_bwg,1)+y_offset),[1,1,size(rfc_bwg,2),size(rfc_bwg,2),1]+x_offset,'r'), hold off;
    %}
end
disp('Scanning Round 1 completed'), toc


E_map = reshape(rfc_small_Es, 31, 31);
k_map = reshape(rfc_small_ks, 31, 31);

EE_sd = mean([std(E_map(:,2:end)-E_map(:,1:end-1)),std(E_map(2:end,:)-E_map(1:end-1,:))])
kk_sd = mean([std(k_map(:,2:end)-k_map(:,1:end-1)),std(k_map(2:end,:)-k_map(1:end-1,:))])

%{
figure, subplot(2,2,[1 2]), imagesc(rfc_bwg), axis xy, title('window1');
subplot(2,2,3), imagesc(E_map), axis xy, title('Emap1');
subplot(2,2,4), imagesc(k_map), axis xy, title('kmap1');

figure, subplot(1,2,1), imagesc(E_map), axis xy, title('Emap');
subplot(1,2,2), imagesc(k_map), axis xy, title('kmap');
%}
