load('new_fermi_finds.mat')
dirac_Es = draw_Es_905_tails;
dirac_ks = draw_ks_905_tails;

aaaa = reshape(region_list{iii},1,961);
aaaa = find(aaaa==1);
[E_interp,K_interp] = meshgrid(-149.5:.5:300, -80:.5:80);
regional_arpes = zeros(size(E_interp));

DE_eV_ave = 0;
for i = aaaa
    cone = cones(:,:,i);
    [E_coor,K_coor] = meshgrid(1:size(cone,2),1:size(cone,1));
    E_coor = E_coor - dirac_Es(i);
    K_coor = K_coor - dirac_ks(i);
    interp_arpes = interp2(E_coor,K_coor,cone,E_interp,K_interp);
    regional_arpes = regional_arpes + interp_arpes;
    
    DPE_pix = dirac_Es(i);
    FL_pix = fermi_(i);
    DPK_pix = dirac_ks(i);
    DE_eV = (FL_pix - DPE_pix)*E_resolve;
    DE_eV_ave = DE_eV_ave + DE_eV;
end
DE_eV_ave = DE_eV_ave / length(aaaa);

E_bin=10;
K_bin=3;

E_binned = E_interp(1,1:E_bin:end-1);

if K_bin == 1
    K_binned = K_interp(:,1);
else    
    K_binned = K_interp(1:K_bin:end,1);
end
Eaxis=E_binned;

regional_arpes_binned = Binning_2d(regional_arpes,E_bin,K_bin);
[region_arpes_symmetrized,Kaxis]=Symmetrized_spectra(region_arpes_binned,K_binned);

figure, imagesc(Eaxis, Kaxis, regional_arpes), axis xy
xticks([0])
xticklabel({num2str(DE_eV_ave)})
