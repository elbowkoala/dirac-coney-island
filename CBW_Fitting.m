

i = 34;
cone_sigma = 2.5;
CBW_Ksize = 50;

cone = cones(:,:,i);
coneb = Binning_2d(cone,bin_E,bin_k);
%[fass,fass_k_off] = kLOSfinder5(cone,bin_E,bin_k);
conebf = imgaussfilt(coneb,cone_sigma);

FL_bpix = round((fermi_(i) - bin_E/2)/bin_E);
DE_bpix = round((ABEK_Es(i)-bin_E/2)/bin_E);
DK_bpix = round((ABEK_ks(i)-bin_k/2)/bin_k);
CBW = conebf(DK_bpix-CBW_Ksize:DK_bpix+CBW_Ksize,DE_bpix:FL_bpix);
CBWB = zeros(size(CBW));

A_it_it = ABEK_As(i);
B_it_it = ABEK_Bs(i);
combadge_yp = A_it_it*(draw_x0) + B_it_it*((draw_x0).^2) + E_0;
combadge_yn = -A_it_it*(draw_x0) + B_it_it*((draw_x0).^2) + E_0;

combadge_itp = horzcat(combadge_yp,draw_x0+K_0);
combadge_itn = horzcat(combadge_yn,draw_x0+K_0);

combadge_itp_curt = combadge_itp(max(1,round(-A_it_it/(2*B_it_it))+K_0):end,:);
combadge_itn_curt = combadge_itn(1:min(length(draw_x0),round(K_0+A_it_it/(2*B_it_it))),:);

combadge_itp_it = reshape(combadge_itp_curt',1,[]);
combadge_itn_it = reshape(combadge_itn_curt',1,[]);

combadge_ITP = insertShape(draw_box, 'Line', combadge_itp_it);
combadge_ITN = insertShape(draw_box, 'Line', combadge_itn_it);
combadge_ITP = combadge_ITP(:,:,1);
combadge_ITN = combadge_ITN(:,:,1);

combadge_IT = combadge_ITP + combadge_ITN; 
combadge_IT = imgaussfilt(combadge_IT, krillin_sigma); %to not penalize intensity from a good surface band that's slightly broad
combadge = combadge_IT;
combadge(combadge~=0)=1;
combadge = abs(1-combadge);

starfleet = bwconncomp(combadge);
combadge(starfleet.PixelIdxList{1})=0;
combadge(starfleet.PixelIdxList{2})=0;
combadge(starfleet.PixelIdxList{3})=0;



parab_x = (1:size(CBW,1))';
parab_A = 0.05;
parab_E0 = 40;
parab_K0 = CBW_Ksize;
parab = parab_A*(parab_x-parab_K0-1).^2 + parab_E0;
parab_pts = horzcat(parab,parab_x);
parab_pts = reshape(parab_pts',1,[]);

CB = insertShape(CBWB,'Line',parab_pts);
figure, imagesc(CB), axis xy

