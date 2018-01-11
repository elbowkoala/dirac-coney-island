
wannasee = 0;
scone_sigma = 8;

tic

kLOS_sat = zeros(1,961);
for i = 1:961%[round(961*rand),round(961*rand),round(961*rand),round(961*rand),round(961*rand), round(961*rand)]%1:1:num_scans;  
    if rem(i,100) == 0
        disp(['Now on scan ',num2str(i), ';  time: ',num2str(toc)])
    end
    cone = result1i(:,:,i); 
    scone = medfilt2(imgaussfilt(cone, scone_sigma),[10,10]);
    scone_mod = scone(:,200:600);%round(2*pre_dos_Es(i)-BCB_Es(i)):round(BCB_Es(i)));
    scone_mod(scone_mod>(mean(scone_mod(:))+2*std(scone_mod(:)))) = mean(scone_mod(:))+2*std(scone_mod(:));

    klosrange = (100:250);
    normcorr = zeros(1,length(klosrange));
    krow_n = 1;
    for krow =  klosrange
        hwidth = min(krow-1, (size(scone_mod,1)-krow-1));
        f_img1 = mat2gray(scone_mod(krow-hwidth:krow+hwidth, :));
        f_img2 = flipud(f_img1);
           
        f_mean = mean(f_img1(:));
        f_denom = sum(sum((f_img1-f_mean).^2));
        
        normcorr(krow_n) = sum(sum((f_img1-f_mean).*(f_img2-f_mean)))/f_denom;
        krow_n = krow_n+1;
    end
    
    kLOS_sat(i) = klosrange(1) + find(normcorr==max(normcorr)) -1;
    
    if wannasee == 1
        figure, imagesc([1:size(scone,1)],fliplr([1:size(scone,2)]),rot90(scone)), axis xy, hold on;
        plot([kLOS_sat(i),kLOS_sat(i)],[0,size(scone,2)],'w')
        %plot([0,size(scone,2)],[kLOS_sat(i),kLOS_sat(i)],'w'), hold off
        title(num2str(i))
    end
    
end

