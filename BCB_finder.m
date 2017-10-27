rfc_ncorr_scan;
ss = imgaussfilt(scone,5);
figure, 
rsss = [5, 10, 15];
rssn = 1;
for rss = rsss    
    sss = ss(rfc_ks_b4(i)-rss:rfc_ks_b4(i)+rss,:);
    ssss = sum(sss);
    smin = find(ssss(:,400:550)==min(ssss(:,400:550)));
    subplot(length(rsss)+1,1,rssn), plot(ssss), hold on, plot([400+smin-1],[ssss(400+smin-1)],'r*')
    rssn = rssn+1;
end
    
subplot(length(rsss)+1,1,rssn), plot(sum(ss));
suptitle(['i=',num2str(i)])
figure, imagesc(imgaussfilt(cones(:,:,i),6)), axis xy
title(['k=',num2str(rfc_ks_b4(i))])