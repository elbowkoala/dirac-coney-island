
wannasee = 1;
scone_sigma = 8;

tic

kLOS_sat = zeros(1,961);
for i = [round(961*rand),round(961*rand),round(961*rand),round(961*rand),round(961*rand), round(961*rand)]%1:1:num_scans;  
    if rem(i,100) == 0
        disp(['Now on scan ',num2str(i), ';  time: ',num2str(toc)])
    end
    cone = cones(:,:,i); 
    scone = cone(:,round(2*pre_dos_Es(i)-BCB_Es(i)):round(BCB_Es(i)));
    scone = imgaussfilt(scone, scone_sigma);
    scone(scone>(mean(scone(:))+2*std(scone(:)))) = mean(scone(:))+2*std(scone(:));

    klosrange = (100:250);
    normcorr = zeros(1,length(klosrange));
    krow_n = 1;
    for krow =  klosrange
        hwidth = min(krow-1, (size(scone,1)-krow-1));
        f_img1 = mat2gray(scone(krow-hwidth:krow+hwidth, :));
        f_img2 = flipud(f_img1);
           
        f_mean = mean(f_img1(:));
        f_denom = sum(sum((f_img1-f_mean).^2));
        
        normcorr(krow_n) = sum(sum((f_img1-f_mean).*(f_img2-f_mean)))/f_denom;
        krow_n = krow_n+1;
    end
    
    kLOS_sat(i) = klosrange(1) + find(normcorr==max(normcorr)) -1;
    
    if wannasee == 1
        figure, imagesc(imgaussfilt(cone,scone_sigma)), axis xy, hold on;
        plot([0,size(cone,2)],[kLOS_sat(i),kLOS_sat(i)],'w'), hold off
    end
    
end

