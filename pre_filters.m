
pre_filter_vecs = cat(1,[]); 


pre_filter_vecs(1,:) = [172, 188, scan_spread];
pre_filter_range(1,:) = [204,276];
pre_filter_vecs(2,:) = rfc_small_ks;
pre_filter_range(2,:) = [172,188];
pre_filter_vecs(3,:) = DPI_big;
pre_filter_range(2,:)
%pre_filter_vecs(4,:) = DPIs;

pre_filter_ranges = [204,276; 172,188; 215, 434];%[0,.0022;  172, 188;  44,63;  1200, 2050];

allfiltered = ones(1,961);
figure;
for NN = 1:size(pre_filter_ranges,1)
    filtered = pre_filter_vecs(NN,:);
    filtered(filtered < pre_filter_ranges(NN,1)) = 0;
    filtered(filtered > pre_filter_ranges(NN,2)) = 0;
    filtered(filtered~=0) = 1;
    allfiltered = allfiltered .* filtered;
    
    subplot(size(pre_filter_ranges,1),1,NN)
    hist(pre_filter_vecs(NN,:),25), hold on;
    plot([pre_filter_ranges(NN,1),pre_filter_ranges(NN,1)],[0,max(hist(pre_filter_vecs(NN,:),25))],'r'), hold on;
    plot([pre_filter_ranges(NN,2),pre_filter_ranges(NN,2)],[0,max(hist(pre_filter_vecs(NN,:),25))],'r'), hold off;   
end



figure, imagesc(reshape(allfiltered,31,31)), axis xy
%{
pre_filter2_vec = ABEK_ks;
pre_filter2_range = [172,188];

pre_filter3_vec = combadges;
pre_filter3_range = [.25, .40];

pre_filter4_vec = ABEK_MCSs;
pre_filter4_range = [44, 63];

filtered1 = pre_filter1_vec;
filtered1(filtered1< pre_filter1_range(1)) = NaN;
filtered1(filtered1 > pre_filter1_range(2)) = NaN;

filtered2 = pre_filter2_vec;
filtered2(filtered2< pre_filter2_range(1)) = NaN;
filtered2(filtered2 > pre_filter2_range(2)) = NaN;

filtered3 = pre_filter3_vec;
filtered3(filtered3< pre_filter3_range(1)) = NaN;
filtered3(filtered3 > pre_filter3_range(2)) = NaN;

allfiltered = filtered1 .* filtered2 .* filtered3;

figure, 
subplot(1,4,1), imagesc(reshape(filtered1,31,31)), axis xy;
subplot(1,4,2), imagesc(reshape(filtered2,31,31)), axis xy;
subplot(1,4,3), imagesc(reshape(filtered3,31,31)), axis xy;
subplot(1,4,4), imagesc(reshape(allfiltered,31,31)), axis xy

allfiltered_bw = allfiltered;
allfiltered_bw(isnan(allfiltered_bw)) = 0;
allfiltered_bw(allfiltered_bw~=0) = 1;
figure, imagesc(reshape(allfiltered_bw,31,31)), axis xy

figure, 
subplot(3,1,1)
hist(pre_filter1_vec,25), hold on;
plot([pre_filter1_range(1),pre_filter1_range(1)],[0,max(hist(pre_filter1_vec,25))],'r'), hold on;
plot([pre_filter1_range(2),pre_filter1_range(2)],[0,max(hist(pre_filter1_vec,25))],'r'), hold off;

subplot(3,1,2)
hist(pre_filter2_vec,25), hold on;
plot([pre_filter2_range(1),pre_filter2_range(1)],[0,max(hist(pre_filter2_vec,25))],'r'), hold on;
plot([pre_filter2_range(2),pre_filter2_range(2)],[0,max(hist(pre_filter2_vec,25))],'r'), hold off;

subplot(3,1,3)
hist(pre_filter3_vec,25), hold on;
plot([pre_filter3_range(1),pre_filter3_range(1)],[0,max(hist(pre_filter3_vec,25))],'r'), hold on;
plot([pre_filter3_range(2),pre_filter3_range(2)],[0,max(hist(pre_filter3_vec,25))],'r'), hold off;
%}

