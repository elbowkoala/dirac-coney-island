

%Scan1=importdata('scan1.txt');
Scan1 = imgaussfilt(cones(:,:,1),2.5);
Scan1 = cones(:,:,1);
Scan1_b = Binning_2d(raw_full_cone,5,2);


%Cent_FL=334;
%Cent_K=188;

Cent_FL = round(mean(fermi_));
Cent_FL_b = round(Cent_FL / bin_E + (bin_E/2));
Cent_K = round(mean(ABEK_ks));
Cent_K_b = round(Cent_K / bin_k + (bin_k/2));

% Cent_K=188; scan6
% Cent_K=188; scan5
% Cent_K=188; scan4
% Cent_K=188; scan3
% Cent_K=186; scan2
% Cent_K=184; scan1


%background=ones(1,498);
%{
background = ones(1,size(Scan1,2));
for j=1:size(Scan1,2);
background(j)=mean(Scan1(Cent_K-80:Cent_K-40,j));
end
for i=1:size(Scan1,1);
Scan1(i,:)=Scan1(i,:)-background(:)';
end
%}

K_width = 100;
K_width_b = round(K_width / bin_k);
K_range =Cent_K-K_width:Cent_K+K_width;
K_range_b = Cent_K_b - K_width_b : Cent_K_b + K_width_b;
E_range=Cent_FL-550 :Cent_FL+50;
E_range_b = Cent_FL_b - 550/bin_E : Cent_FL_b + 50/bin_E;

E_resolve=1.599./(496*2); %fixed resolution 8/16/17 (0.17595-0.15177)./10; %in eV
E_resolve_b = E_resolve * bin_E;
X_axis=0.512*sin((K_range-Cent_K)*0.04631./180*3.1415*14/30)*sqrt(110-4); %mtm for each pixel, fixed res 8/16/17 (added *14/30, changed 117 to 110)
X_axis_b = 0.512*sin((K_range_b-Cent_K_b)*0.04631./180*3.1415*14/30)*sqrt(110-4) * bin_k;
Y_axis=(E_range-Cent_FL)*E_resolve;
Y_axis_b = (E_range_b - Cent_FL_b)*E_resolve_b;



plot_x=ones(length(E_range),1)*X_axis;
plot_y=ones(K_width*2+1,1)*Y_axis;
plot_x_b = ones(length(E_range_b),1)*X_axis_b;
plot_y_b = ones(K_width_b*2+1,1)*Y_axis_b;


figure
plot_yy = plot_y';
imagesc(X_axis,flipud(Y_axis'), rot90(norman(Scan1(K_range,E_range),0,3),1)), axis xy, hold on;
shading flat,
colormap gray,
title('Scan 1')
xlabel('k (A^{-1})')
ylabel('E-E_{B} (eV)')
set(gca, 'TickDir', 'out')


figure
plot_yy_b = plot_y_b';
imagesc(plot_x_b(1,:),flipud(plot_yy_b(:,1)),rot90(norman(Scan1_b(K_range_b,E_range_b),0,5),1)), axis xy
colormap jet
set(gca, 'TickDir', 'out')
%hold on;
%plot([plot_x_b(1),plot_x_b(end)],[0,0],'r'), hold off
shading flat,
title('Scan 1')
%}
%axis([-0.2,0.2,-0.35,0.1])
%caxis([40,3000]) 
% caxis([40,3000]) scan6
% caxis([40,3000]) scan5
% caxis([40,3000]) scan4
% caxis([40,3000]) scan3
% caxis([0,6700])  scan2
% caxis([15,700])  scan1

