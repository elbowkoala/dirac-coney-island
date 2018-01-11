
above_box_eV_ranges = [+.6,+.25];%[-.075,-.085; -.075,-.1; -.075,-.125; -.075,-.15; -.075,-.175];


for beVi = 1:size(above_box_eV_ranges,1)
    above_box_eV_range = above_box_eV_ranges(beVi,:)
    above_box_invA_range = [-.35,-.25]; %Actually use whole k-range of cones data
    above_be2 = above_box_eV_range(2)/pix2eV;
    above_be1 = above_box_eV_range(1)/pix2eV;
    above_bk1 = above_box_invA_range(1)/pix2invA;
    above_bk2 = above_box_invA_range(2)/pix2invA;

    above_Iperpix = zeros(1,961);
     for i = find(kLOS>0)
            above_res1i = result1i(:,:,i);
            if kLOS(i) + above_bk1 < 1
                above_res1i_box = above_res1i(1:round(kLOS(i)+above_bk2),round(AVE_FL-above_be1):round(AVE_FL-above_be2));
            elseif kLOS(i) + above_bk2 > size(above_res1i,1)
                above_res1i_box = above_res1i(round(kLOS(i)+above_bk1):end,round(AVE_FL-above_be1):round(AVE_FL-above_be2));
            else
                above_res1i_box = above_res1i(round(kLOS(i)+above_bk1):round(kLOS(i)+above_bk2), round(AVE_FL-above_be1):round(AVE_FL-above_be2));%res1i(round(kLOS(i)+bk1):round(kLOS(i)+bk2), round(AVE_FL-be1):round(AVE_FL-be2));
            end
            above_number_of_pixels = size(above_res1i_box,1)*size(above_res1i_box,2);
            above_Iperpix(i) =  sum(sum(above_res1i_box))/above_number_of_pixels;
    end
    ave_above_Iperpix = mean(above_Iperpix)%(above_Iperpix>0))
    % panel_mean = zeros(1,size(region_list,1)); 
    figure,
    imagesc(imgaussfilt(result1i(:,:,random_scan),5)), axis xy, hold on;
    if size(above_res1i_box,1) == size(result1i(:,:,1),1)
        plot([round(AVE_FL-above_be1),round(AVE_FL-above_be1),round(AVE_FL-above_be2),round(AVE_FL-above_be2),round(AVE_FL-above_be1)],...
        [1,size(res1i,1),size(res1i,1),1,1],'r', 'DisplayName','ave above FL'), hold on;
    else
        plot([round(AVE_FL-above_be1),round(AVE_FL-above_be1),round(AVE_FL-above_be2),round(AVE_FL-above_be2),round(AVE_FL-above_be1)],...
        [round(kLOS(i)+above_bk1),round(kLOS(i)+above_bk2),round(kLOS(i)+above_bk2),round(kLOS(i)+above_bk1),round(kLOS(i)+above_bk1)],'r','DisplayName','ave above FL'), hold on;
    end
    title({['average Iperpix all scans = ',num2str(ave_above_Iperpix)]; ['E range [',num2str(above_box_eV_range(1)),', ',num2str(above_box_eV_range(2)),']eV']});
end