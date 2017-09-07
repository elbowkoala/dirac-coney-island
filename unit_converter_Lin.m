i = 833;

%Scan1=importdata('scan1.txt');
Scan1 = cones(:,:,i);


%Cent_FL=334;
%Cent_K=188;

Cent_FL = round(fermi_(i));
Cent_K = round(draw_ks_906(i));

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

K_width=100;
K_range=Cent_K-K_width:Cent_K+K_width;
E_range=Cent_FL-550:Cent_FL+50;

E_resolve=1.599./(496*2); %fixed resolution 8/16/17 (0.17595-0.15177)./10; %in eV
X_axis=0.512*sin((K_range-Cent_K)*0.04631./180*3.1415*14/30)*sqrt(110-4); %mtm for each pixel, fixed res 8/16/17 (added *14/30, changed 117 to 110)
Y_axis=(E_range-Cent_FL)*E_resolve;

plot_x=ones(length(E_range),1)*X_axis;
plot_y=ones(K_width*2+1,1)*Y_axis;



figure
pcolor(plot_x, plot_y', Scan1(K_range,E_range)'), hold on;
shading flat,
title('Scan 1')
xlabel('k (A^{-1})')
ylabel('E-E_{B} (eV)')
%axis([-0.2,0.2,-0.35,0.1])
%caxis([40,3000]) 
% caxis([40,3000]) scan6
% caxis([40,3000]) scan5
% caxis([40,3000]) scan4
% caxis([40,3000]) scan3
% caxis([0,6700])  scan2
% caxis([15,700])  scan1

