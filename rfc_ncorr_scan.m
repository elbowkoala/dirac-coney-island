tic; 
%load 'cones.mat'
%load 'ssd_big_scan_171003.mat'

wannasee = 0;
scan_is =  i;% [round(961*rand), round(961*rand),round(961*rand),round(961*rand)] ;
%{
rfc_Es_b4 = zeros(1,num_scans);
rfc_ks_b4 = zeros(1,num_scans);
rfc_wedges_b4 = zeros(1,num_scans); 
rfc_corrs_b4 = zeros(1,num_scans);
rfc_corrspreads_b4 = zeros(1,num_scans);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_E = 5;
bin_k = 2;

scone_sigma = 5;
ncorr_TH = 0.95;

rfc_patch_E_b_range = [51:95];  
rfc_patch_k_b_range = [51:130];

wedge_Emin = 34;
wedge_Emax = 45;
wedge_kmin = 33;
wedge_kmax = 47;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Making the 0th-Order Template %%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_full_cone = zeros(size(cone_range_K,2),size(cone_range_E,2));

for i = 1:1:num_scans
    raw_full_cone = raw_full_cone + cones(:,:,i);
end

rfc_b = Binning_2d(raw_full_cone, bin_E, bin_k);
rfc_patch = rfc_b(rfc_patch_k_b_range, rfc_patch_E_b_range);
%rfc_patch = rfc_b(40:140,65:210);
%rfc_patch = rfc_b(:,195:235);
%[rfc_patch_sym, rfc_patch_symkaxis] = Symmetrized_spectra(rfc_patch,(1:size(rfc_patch,1))');
rfc_patch_gray = mat2gray(rfc_patch);
template = rfc_patch_gray;%IT_processor(rfc_patch_gray,5,5);

template_DPK_coor = 41;
template_DPE_coor = 23;

wedge = insertShape(zeros(size(template)),'Line',[wedge_Emax,wedge_kmin, wedge_Emin,round(mean([wedge_kmin,wedge_kmax])), wedge_Emax, wedge_kmax]);
wedge = wedge(:,:,1);
wedge(wedge~=0) = 1;
for w_i = 30:size(template,2)
    wedge_col = wedge(:,w_i);
    pntA = find(wedge_col==1,1,'first');
    pntB = find(wedge_col==1,1,'last');
    wedge_col(pntA:pntB) = 1;
    wedge(:,w_i) = wedge_col;
end

mean_ssd_big_k = mean(ssd_big_ks);

tic
for i = scan_is%[round(961*rand),round(961*rand),round(961*rand),round(961*rand),round(961*rand), round(961*rand)]%1:1:num_scans;  
    if rem(i,100) == 0
        disp(['Now on scan ',num2str(i), ';  time: ',num2str(toc)])
    end
    b_E = round((ssd_big_Es(i)-bin_E/2)/bin_E);
    b_k = round((mean_ssd_big_k-bin_k/2)/bin_k);
    scone_b_E_range = [round( b_E - size(template,2)/2 - 20) : round(b_E + size(template,2)/2 + 20)];
    scone_b_k_range = [round( b_k - size(template,1)/2 - 10) : round(b_k + size(template,1)/2 + 10)];
    
    
    
    scone = cones(:,:,i);
    
    scone = imgaussfilt(scone, scone_sigma);
    scone_b = Binning_2d(scone, bin_E, bin_k);
    scone_b = scone_b(scone_b_k_range,scone_b_E_range);
    
    scone_bg = mat2gray(scone_b);    
    image = scone_bg;

    N = size(template,1);
    M = size(template,2);
    
    normcorr = zeros([(size(image,1)-N+1), (size(image,2)-M+1)]);
    f_mean = zeros(size(normcorr));
    f_denom = zeros(size(normcorr));
    t_mean = nanmean(template(:));
    t_denom = sum(nansum((template-t_mean).^2));


    for x = 1:(size(image,2) - M +1)
        for y = 1:(size(image,1) - N +1)
            f_mean(y,x) = mean( mean( image(y:y+N-1, x:x+M-1)));
%              f_denom = 0;
            %f_denom = sum(nansum((f_image(y:y+N-1,x:x+M-1)-f_mean).^2));
            %NC = nansum( dot( (IT_template(1:N,1:M) - t_mean),(f_image(y:y+N-1,x:x+M-1)-f_mean)));
            %normcorr(y,x) = sum(NC) / sqrt(t_denom * f_denom);

            for m = 1:size(template,2)
                %normcol = (IT_template(1:N,m) - t_mean) .* (f_image(y:y+N-1,x+m-1) - f_mean);
                %normcorr(y,x) = normcorr(y,x) + nansum(normcol);
                %f_denomcol = (f_image(y:y+N-1,x+m-1)-f_mean).^2;
                %f_denom = f_denom + nansum(f_denomcol);

                for n = 1:size(template,1)
                    if isnan(template(n,m))==1
                        continue
                    end
                    normcorr(y,x) = normcorr(y,x) + (template(n,m)-t_mean)*(image(y+n-1,x+m-1)-f_mean(y,x));
                    %t_denom = t_denom + (IT_template(n,m)-t_mean).^2;
                    f_denom(y,x) = f_denom(y,x) + (image(y+n-1,x+m-1)-f_mean(y,x)).^2;
                    %normcorr(y,x) = normcorr(y,x) + (fass_image(y+n-1, x+m-1) - IT_template(n,m)).^2;
                end

            end
            %if (x == 5) && (y==5)
            %    disp(['AB_IT_i=',num2str(AB_IT_i),'t_denom=',num2str(t_denom)])
            %end
            normcorr(y,x) = normcorr(y,x) / sqrt(t_denom * f_denom(y,x));

        end
    end

    [corrmax_row,corrmax_col] = find(normcorr==max(normcorr(:)));
    rfc_corrs_b4(i) = max(normcorr(:));
    
    DPE_b = corrmax_col + template_DPE_coor + scone_b_E_range(1)-1 - 1;
    DPK_b = corrmax_row + template_DPK_coor +scone_b_k_range(1)-1 - 1;

    rfc_Es_b4(i) = DPE_b * bin_E + bin_E/2;
    rfc_ks_b4(i) = DPK_b * bin_k + bin_k/2;
    
    overlay = zeros(size(image));
    wedge_overlay = zeros(size(image));
    overlay((1:size(template,1))+corrmax_row-1,(1:size(template,2))+corrmax_col-1) = template;
    wedge_overlay((1:size(template,1))+corrmax_row-1,(1:size(template,2))+corrmax_col-1) = wedge;
    
    rfc_wedges_b4(i) = sum(dot(wedge_overlay,image))/(nnz(wedge_overlay));
    [wedge_row,wedge_col] = find(wedge~=0);
    
    
    normcorr_g = mat2gray(normcorr);
    %rfc_corrspread = length(find(normcorr_g > normcorrspread_TH));
    rfc_corrspreads_b4(i) = length(find(normcorr_g > ncorr_TH));
    
    if wannasee == 1
        
        figure

        subplot(2,2,1)
        imagesc((1:size(template,1)), fliplr(1:size(template,2)),rot90(template,1)), axis xy, hold on;
        plot([wedge_row],[wedge_col],'r.'), hold off
        title(['maxcorr = ',num2str(rfc_corrs_b4(i))])


        subplot(2,2,2)
        imagesc((1:size(image,1)),fliplr(1:size(image,2)),rot90(image,1)), axis xy
        title(['DPI = ',num2str(round(DPI_big(i)))])

        subplot(2,2,3)
        imagesc((1:size(normcorr,1)),fliplr(1:size(normcorr,2)),rot90(normcorr,1)), axis xy, hold on
        plot(corrmax_row,corrmax_col,'r+'), hold off
        title(['spread = ',num2str(round(rfc_corrspreads_b4(i),3))])


        subplot(2,2,4)
        imagesc((1:size(image,1)),fliplr(1:size(image,2)),rot90(overlay+image,1)), axis xy
        title(['ws = ',num2str(round(rfc_wedges_b4(i),3))])

        suptitle(['i=',num2str(i),'  E=',num2str(rfc_Es_b4(i)),'  k=',num2str(rfc_ks_b4(i))])
        pause(.01)
    end
end
%disp('Scanning Round 1 completed'), toc
toc

%E_small_map = reshape(rfc_Es_b4, 31, 31);
%k_small_map = reshape(rfc_ks_b4, 31, 31);

%EE_small_sd = mean([std(E_small_map(:,2:end)-E_small_map(:,1:end-1)),std(E_small_map(2:end,:)-E_small_map(1:end-1,:))]);
%kk_small_sd = mean([std(k_small_map(:,2:end)-k_small_map(:,1:end-1)),std(k_small_map(2:end,:)-k_small_map(1:end-1,:))]);

