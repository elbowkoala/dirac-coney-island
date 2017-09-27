pix2eV = (1.599/(2*496));
pix2invA = 0.512*sqrt(110-4)*0.04631*(14/30)*(3.1415/180);
load('new_fermi_finds.mat')
%%%Run DPI_map_maker.m with correct E,k arrays%%%


E_map_pix = E_map;%reshape(dot_dirac_Es,31,31);%reshape(ABEK_Es,31,31);
FL_map_pix = reshape(fermi_,31,31);
K_map_pix = k_map;%reshape(ABEK_ks,31,31);
V_map_pix = (reshape(ABEK_As,31,31));
B_map_pix = reshape(ABEK_Bs,31,31);
MC_map = norman(reshape(TBDs,31,31),0,5);

%DPI_map_maker;
%DPI_map_pix = DPI_map;  

E_map_eV = - (FL_map_pix - E_map_pix) .*  pix2eV;
K_map_invA = (K_map_pix-mean(K_map_pix(:))) .* pix2invA;
V_map_eVA = V_map_pix * (bin_E/bin_k) * .8107;

figure;
subplot(2,3,1)
imagesc([1:31].*15,[1:31].*15,E_map_eV), axis xy
title('Energy Map')
xlabel('X Position (um)')
ylabel('Y Position (um)')
colormap jet
c = colorbar();
c.Label.String = 'DP Energy (eV)';
c.Label.FontSize = 8;

subplot(2,3,2)
imagesc([1:31].*15,[1:31].*15,K_map_invA), axis xy
title('Momentum Map')
xlabel('X Position (um)')
ylabel('Y Position (um)')
colormap jet
c = colorbar();
c.Label.String = 'DP Momentum (A^{-1})';
c.Label.FontSize = 8;

subplot(2,3,4)
imagesc(V_map_eVA), axis xy
title('Dirac Velocity Map')
xlabel('X Position')
ylabel('Y Position')
c = colorbar();
c.Label.String = 'Dirac Velocity (eV-A)';
c.Label.FontSize = 8;

subplot(2,3,5)
imagesc(B_map_pix), axis xy
title('Square Term (pix units)')
xlabel('X Position')
ylabel('Y Position')
c = colorbar();
c.Label.String = 'Square term (arb pix units)';
c.Label.FontSize = 8;

%{
subplot(2,3,3)
imagesc(DPI_map), axis xy
title('DP Intensity Map')
xlabel('X Position')
ylabel('Y Position')
c = colorbar();
c.Label.String = 'DP Region Intensity (arb. units)';
c.Label.FontSize = 8;
%}
subplot(2,3,6)
imagesc(MC_map), axis xy
title('Max Corr Map')
xlabel('X Position')
ylabel('Y Position')
c = colorbar();
c.Label.String = 'Max Corr (arb. units)';
c.Label.FontSize = 8;

