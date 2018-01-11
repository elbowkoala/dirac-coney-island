kcutoff = 100;
figure
ii = 1;
for i = 1:5
    I = round(961*rand);
    while dos_Es(I) == 0
       I = round(961*rand);
    end
    ksympixel = kLOS(I);
    cone = result1i(:,:,I);
    cone_mod = imgaussfilt(cone,1);
    kweights = abs([1:size(cone,1)]-ksympixel);
    kweights(kweights>kcutoff) = 0;
    DOS = sum(repmat(kweights, size(cone,2), 1)'.*cone_mod);
    
    int_width = 40/2;
    int_dos_x = int_width+1:size(cone,2)-int_width;
    int_dos = zeros(size(int_dos_x));
    int_dos_i = 1;
    for iii = int_dos_x
        int_dos(int_dos_i+int_width) = trapz(DOS(iii-int_width:iii+int_width));
        int_dos_i = int_dos_i + 1;
    end
    
    if rem(ii,5) == 0
        plot_place = (floor(ii/5)-1)*3*5 + 5;
    else
        plot_place = (floor(ii/5))*3*5 + rem(ii,5);
    end
    subplot(3,5,plot_place)
    imagesc(imgaussfilt(cone,5)), axis xy, hold on;
    plot([dos_Es(I),dos_Es(I)],[300,-1200],'r')
    title(num2str(I),'Fontsize',5)

    subplot(3,5,plot_place+5)
    plot(DOS), hold on;
    plot([dos_Es(I),dos_Es(I)],[0,5000],'r')
    xlim([1,800])
    ylim([0,.012*trapz(DOS)])
    
    subplot(3,5,plot_place+10)
    plot(int_dos), hold on
    plot([dos_Es(I),dos_Es(I)],[0,400000],'r')
    xlim([1,800])
    ylim([0,.006*trapz(int_dos)])
    
    ii = ii + 1;
end

% 
% figure,
% subplot(3,1,1), imagesc(imgaussfilt(cone,5)), axis xy
% subplot(3,1,2), plot(DOS), xlim([1,800])
% subplot(3,1,3), plot(region1_areas), xlim([1-int_width,800-int_width])
%    
dosfit_width = 50/2;
dosfit_m = zeros(1,size(int_dos,2));
for II = 100:700
    dosfit_x = II-dosfit_width:II+dosfit_width;
    dosfit_y = int_dos(dosfit_x);
    linefit = polyfit(dosfit_x, dosfit_y, 1);
    dosfit_m(II) = linefit(1);
end

figure
subplot(2,1,1)
imagesc(imgaussfilt(cone,5)), axis xy
subplot(2,1,2)
yyaxis left
plot(int_dos)
yyaxis right
plot(dosfit_m), hold on;
plot([1,800],[0,0],'r')



