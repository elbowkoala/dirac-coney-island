tic; 

rfc_FL_Es = zeros(1,num_scans);
rfc_FL_ks = zeros(1,num_scans);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_E = 3;
bin_k = 2;

rfc_FL_sigma = 10;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Making the 0th-Order Template %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
raw_full_cone = zeros(size(cone_range_K,2),size(cone_range_E,2));

for i = 1:1:num_scans;
    raw_full_cone = raw_full_cone + cones(:,:,i);
end
%}
rfc_b = Binning_2d(raw_full_cone, bin_E, bin_k);
rfc_FL_window_Kpix = (30:150);
rfc_FL_window_Epix = (190:230);
rfc_FL_window = rfc_b(rfc_FL_window_Kpix, rfc_FL_window_Epix); 
rfc_FL_bwg = mat2gray(rfc_FL_window);

window_FL_k_coor = size(rfc_FL_bwg,1)/2;
window_FL_E_coor = 26; %%%%%%%%%%%%%%% hand-pick point from rfc_bg%%%%%%%%%%%

disp('Template step completed '), toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Round 1 Scan with Large Window Including Bulk Bands %%%%%%%%%%%%%

scan_maxes_FL = zeros(1,num_scans);
%scan_spread = zeros(1,num_scans);
%figure,
NNN=0;
for i =  1:num_scans%round(961*rand)%[round(961*rand),round(961*rand),round(961*rand),round(961*rand),round(961*rand), round(961*rand)]%1:1:num_scans;  
    NNN=NNN+1;
    if rem(i,100) == 0
        disp(['Now on scan ',num2str(i)])
    end
    
    scone = cones(:,:,i);
    
    scone = imgaussfilt(scone, rfc_FL_sigma);
    scone_b = Binning_2d(scone, bin_E, bin_k);
    scone_b = scone_b(:,round(100*5/bin_E)+1:end);
    scone_bg = mat2gray(scone_b);

    scann = normxcorr2( rfc_FL_bwg, scone_bg );
    scan=scann(size(rfc_FL_bwg,1)+1:end-size(rfc_FL_bwg,1),size(rfc_FL_bwg,2)+1:end-size(rfc_FL_bwg,2));
    [scan_max, scan_max_i] = max(scan(:));
    [smax_y,smax_x] = ind2sub(size(scan),scan_max_i);
    scan_maxes_FL(i) = scan_max;
    
    [y_peak, x_peak] = find(scan==max(scan(:)));
    y_offset = y_peak - size(rfc_FL_bwg,1) + size(rfc_FL_bwg,1);
    x_offset = x_peak - size(rfc_FL_bwg,2) + size(rfc_FL_bwg,2);
        
    rfc_FL_Es(i) = (x_offset + round(100*5/bin_E) + window_FL_E_coor  -1)*bin_E + bin_E/2; %%%%%%!!!!%%
    rfc_FL_ks(i) = (y_offset + window_FL_k_coor  -1)*bin_k + bin_k/2;    
   
    scan_gray = mat2gray(scan);
    %{
    subplot(1,4,1), %imagesc(scan_Kaxis,scan_Eaxis,rot90(scan,-1)), axis xy, title(num2str(i)), hold on;
    imagesc(rot90(scan,-1)), axis xy, hold on;
    plot(size(scan,1)-y_peak,x_peak,'b*'), hold off;
        title({['Sc max=',num2str(round(scan_maxes_small(i),3))];['spread=',num2str(round(scan_spread(i)))]});
        
    subplot(1,4,2), imagesc(rot90(rfc_FL_bwg,-1)), axis xy;
    title(['DPI=',num2str(round(DPI_small(i),2))]);

    
    subplot(1,4,3), imagesc(rot90(scone_bg,-1)), axis xy, hold on;
    %plot([1,size(rfc_small_bwg,2),size(rfc_small_bwg,2),1,1]+x_offset,[1,1,size(rfc_small_bwg,1),size(rfc_small_bwg,1),1]+y_offset,'r'), hold off;
    plot([0,0,-size(rfc_FL_bwg,1),-size(rfc_FL_bwg,1),1]+size(scone_bg,1)-y_offset,[1,size(rfc_FL_bwg,2),size(rfc_FL_bwg,2),1,1]+x_offset,'r'), hold on;
    %plot([size(scone_bg,1)-y_offset-window_FL_k_coor],[x_offset+window_FL_E_coor],'r*'), hold on;
    plot([0,-size(rfc_FL_bwg,1)]+size(scone_bg,1)-y_offset,[x_offset+window_FL_E_coor,x_offset+window_FL_E_coor],'w'), hold off;
    %plot(flipud([1,size(rfc_small_bwg,1),size(rfc_small_bwg,1),1,1]+size(scone_bg,1)-size(rfc_bwg,1)+y_offset),[1,1,size(rfc_bwg,2),size(rfc_bwg,2),1]+x_offset,'r'), hold off;
    title(['E=',num2str(rfc_FL_Es(i)),' k=',num2str(rfc_FL_ks(i))]);
    
    subplot(1,4,4), imagesc(rot90(imgaussfilt(cones(:,:,i),rfc_FL_sigma),-1)), axis xy, hold on;
    plot([1,300],[rfc_FL_Es(i),rfc_FL_Es(i)],'w'), hold off;
    title(['i=',num2str(i)]);
    %}
end
disp('Scanning Round 1 completed'), toc

%{
E_FL_map = reshape(rfc_FL_Es, 31, 31);
k_FL_map = reshape(rfc_FL_ks, 31, 31);

EE_small_sd = mean([std(E_FL_map(:,2:end)-E_FL_map(:,1:end-1)),std(E_FL_map(2:end,:)-E_FL_map(1:end-1,:))])
kk_small_sd = mean([std(k_FL_map(:,2:end)-k_FL_map(:,1:end-1)),std(k_FL_map(2:end,:)-k_FL_map(1:end-1,:))])


figure, subplot(2,2,[1 2]), imagesc(rfc_bwg), axis xy, title('window1');
subplot(2,2,3), imagesc(E_map), axis xy, title('Emap1');
subplot(2,2,4), imagesc(k_map), axis xy, title('kmap1');

figure, subplot(1,2,1), imagesc(E_map), axis xy, title('Emap');
subplot(1,2,2), imagesc(k_map), axis xy, title('kmap');
%}
