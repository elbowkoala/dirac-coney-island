above_box_eV_range = [-.175,-.075];%
above_box_invA_range = [-.2,.2]; %Actually use whole k-range of cones data

above_boxtop_eV = above_box_eV_range(1);%.24;% 0.2;
above_boxbot_eV = above_box_eV_range(2);%.30;%0.30;

above_be2 = round(above_boxtop_eV/pix2eV);
above_be1 = round(above_boxbot_eV/pix2eV);

above_bk1 = above_box_invA_range(1)/pix2invA;
above_bk2 = above_box_invA_range(2)/pix2invA;

above_Iperpix = zeros(1,961);
for i = find(kLOS>0)
    res1i = result1i(:,:,i);
%     if ((kLOS(i) + bk1) < 0) || ((kLOS(i)+bk2) > 300)
%             continue
%         end
    res1i_box = res1i(:, round(AVE_FL-above_be1):round(AVE_FL-above_be2));%res1i(round(kLOS(i)+bk1):round(kLOS(i)+bk2), round(AVE_FL-be1):round(AVE_FL-be2));
    number_of_pixels = size(res1i_box,1)*size(res1i_box,2);
    above_Iperpix(i) =  sum(sum(res1i_box))/number_of_pixels;
end

figure, 
random_scan = round(961*rand);
imagesc(imgaussfilt(result1i(:,:,random_scan),5)), axis xy, hold on;
if size(res1i_box,1) == size(res1i,1)
    plot([round(AVE_FL-above_be1),round(AVE_FL-above_be1),round(AVE_FL-above_be2),round(AVE_FL-above_be2),round(AVE_FL-above_be1)],...
    [1,size(res1i,1),size(res1i,1),1,1],'r'), hold off;
else
    plot([round(AVE_FL-above_be1),round(AVE_FL-above_be1),round(AVE_FL-above_be2),round(AVE_FL-above_be2),round(AVE_FL-above_be1)],...
    [round(kLOS(i)+above_bk1),round(kLOS(i)+above_bk2),round(kLOS(i)+above_bk2),round(kLOS(i)+above_bk1),round(kLOS(i)+above_bk1)],'r'), hold off;
end
title(['scan ',num2str(random_scan)])
above_aveIperpix = mean(above_Iperpix(above_Iperpix>0))



side_box_eV_range = [0.4,0.45];%
side_box_invA_range = [.18,.23];%

side_boxtop_eV = side_box_eV_range(1);%.24;% 0.2;
side_boxbot_eV = side_box_eV_range(2);%.30;%0.30;

side_be2 = round(side_boxtop_eV/pix2eV);
side_be1 = round(side_boxbot_eV/pix2eV);

side_bk1 = side_box_invA_range(1)/pix2invA;
side_bk2 = side_box_invA_range(2)/pix2invA;

side_Iperpix = zeros(1,961);
for i = find(kLOS>0)
    res1i = result1i(:,:,i);
    if ((kLOS(i) + side_bk1) < 0) || ((kLOS(i)+side_bk2) > 300)
            continue
        end
    res1i_box = res1i(round(kLOS(i)+side_bk1):round(kLOS(i)+side_bk2), round(AVE_FL-side_be1):round(AVE_FL-side_be2));
    number_of_pixels = size(res1i_box,1)*size(res1i_box,2);
    side_Iperpix(i) =  sum(sum(res1i_box))/number_of_pixels;
end

figure, 
ax = subplot(1,2,1);
random_scan = round(961*rand);
imagesc(imgaussfilt(result1i(:,:,random_scan),5)), axis xy, hold on;
plot([round(AVE_FL-side_be1),round(AVE_FL-side_be1),round(AVE_FL-side_be2),round(AVE_FL-side_be2),round(AVE_FL-side_be1)],...
    [round(kLOS(i)+side_bk1),round(kLOS(i)+side_bk2),round(kLOS(i)+side_bk2),round(kLOS(i)+side_bk1),round(kLOS(i)+side_bk1)],'r'), hold off;
title(['scan ',num2str(random_scan)])
ax = subplot(1,2,2);
hh = histogram(side_Iperpix(side_Iperpix>0),12);


side_aveIperpix = mean(side_Iperpix(side_Iperpix>0))