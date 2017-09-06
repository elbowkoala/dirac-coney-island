%function [DPI_map] = DPI_map_maker(dirac_Es,dirac_ks,Epixwidth,kpixwidth)
dirac_Es = draw_Es_905_tails;
dirac_ks = draw_ks_905_tails;
Epixwidth = 100;
kpixwidth = 100;

DPIs = zeros(1,num_scans);

for i = 1:num_scans
    cone = cones(:,:,i);
    
    %rfc_DPE = round(dot_dirac_Es(i));
    %rfc_DPk = round(dot_dirac_ks(i));
    %rfc_DPI(i) = sum(sum(cone(rfc_DPk-50:rfc_DPk+50, rfc_DPE-50:rfc_DPE+50)));
    
    DPE = round(dirac_Es(i));
    DPk = round(dirac_ks(i));
    DPIs(i) = sum(sum(cone(DPk-round(kpixwidth/2):DPk+round(kpixwidth/2),...
        DPE-round(Epixwidth/2):DPE+round(Epixwidth/2))));
end

DPI_map = reshape(DPIs,31,31);


%figure, imagesc(line_DPI_map_829), axis xy, title('line scan 170824: DP Intensity')


%{
figure, subplot(2,2,1), imagesc(E_map), axis xy, title('rfc E map');
subplot(2,2,2), imagesc(rfc_DPI_map), axis xy, title('rfc DPI map');
subplot(2,2,3), imagesc(reshape(ABEK_Es,31,31)), axis xy, title('ABEK E map');
subplot(2,2,4), imagesc(ABEK_DPI_map), axis xy, title('ABEK DPI map');
%}
    