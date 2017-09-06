load('new_fermi_finds.mat')
dirac_Es = draw_Es_905_tails;
dirac_ks = draw_ks_905_tails;
iii = 5;
aaaa = reshape(region_list{iii},1,961);
aaaa = find(aaaa==1);
[E_interp,K_interp] = meshgrid(-149.5:.5:350, -80:.5:80);
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
    DE_eV = (-FL_pix + DPE_pix)*E_resolve;
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
X_axis=0.512*sin((K_range-Cent_K)*0.04631./180*3.1415*14/30)*sqrt(110-4); %mtm for each pixel, fixed res 8/16/17 (added *14/30, changed 117 to 110)
K_resolve_approx = 0.512*0.04631/180*3.1415*14/30*sqrt(110-4);

regional_arpes_binned = Binning_2d(regional_arpes,E_bin,K_bin);
[region_arpes_symmetrized,Kaxis]=Symmetrized_spectra(region_arpes_binned,K_binned);

figure 
%subplot(1,2,2)
imagesc(Kaxis, flip(Eaxis), rot90(regional_arpes)), axis xy
yticks([0,-(DE_eV_ave/E_resolve)])
yticklabels({num2str(DE_eV_ave),'0'})
ylabel_str = ' E - E$_{F}$, (eV)';% \int_{0}^{2} x^2\sin(x) dx $$'
ylabel(ylabel_str,'Interpreter','latex')
xticks([-50,0,50])
xticklabels({num2str(-50*K_resolve_approx),'0',num2str(50*K_resolve_approx)})
xlabel('Momentum ({\AA}$^{-1}$)','Interpreter','latex')
