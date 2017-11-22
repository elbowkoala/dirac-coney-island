
wannasee = 0;
scone_sigma = 5;

tic
kLOS = zeros(1,961);
for i = 1:961% [868,908,472,470,325,865]% [round(961*rand),round(961*rand),round(961*rand),round(961*rand),round(961*rand), round(961*rand)]%1:1:num_scans;  
    if rem(i,100) == 0
        disp(['Now on scan ',num2str(i), ';  time: ',num2str(toc)])
    end
    if DPI_big(i) < 200
        disp('Scan is low DPI; skipping')
        continue
    end
    scone = result1i(:,:,i); 
    sconee = scone(:,1:round(rfc_FL_Es(i)));
    sconee = imgaussfilt(sconee, scone_sigma);
    sconeee = sconee;
    sconeee(sconee<(mean(sconee(:))-0*std(sconee(:))))=0;
    sconeee(sconee>(mean(sconee(:))+2*std(sconee(:)))) = mean(scone(:))+2*std(sconee(:));
    %figure, imagesc(sconeee), axis xy
    
    klosrange = (100:250);
    normcorr = zeros(1,length(klosrange));
    krow_n = 1;
    for krow =  klosrange
        hwidth = 40;%min(krow-1, (size(sconee,1)-krow-1));
        f_img1 = mat2gray(sconee(krow-hwidth:krow+hwidth, :));
        f_img2 = flipud(f_img1);
           
        f_mean = mean(f_img1(:));
        f_denom = sum(sum((f_img1-f_mean).^2));
        
        normcorr(krow_n) = sum(sum((f_img1-f_mean).*(f_img2-f_mean)))/f_denom;
        krow_n = krow_n+1;
    end
    %figure, imagesc(normcorr), axis xy
    
    kLOS(i) = klosrange(1) + find(normcorr==max(normcorr)) -1;
    
    if wannasee == 1
        disp(['i=',num2str(i),' kLOS=',num2str(kLOS(i))])
        
        figure, imagesc(imgaussfilt(scone,scone_sigma)), axis xy, hold on;
        plot([0,size(sconee,2)],[kLOS(i),kLOS(i)],'w'), hold off
        title(['i=',num2str(i),' kLOS=',num2str(kLOS(i))])
    end
    
end

