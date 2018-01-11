load_data = 1;
if load_data == 1
    load 'cones.mat';
    load 'ssd_big_scan_171003.mat';
    load 'rfc_big_scan_170927.mat';
    load rfc_FL_scan_170927.mat;
    load rfc_ncorr_scan_171019w.mat;
    load kLOS_171122.mat;
    load pre_dos_Es;
    load BCB_E_scan.mat;
    bad_DPI = find(DPI_big < 200);
end

global x_x y_y FL_param
global x__x y__y %bcb_find dos_start
wannasee = 1;

dos_k_cutoff = 90;
ssKpm = 5;
bcbmarj = 25;
a1_range = [0:10:50];
bcb_sigma = 8;
dos_sigma = 3;


%final_bcb_params = zeros(4,961);
final_bcb_Rsquareds_180109 = zeros(1,961);
bcb_finds_180109 = zeros(1,961);
%final_dos_params = zeros(4,961);
final_dos_Rsquareds_180109 = zeros(1,961);
dos_Es_180109 = zeros(1,961);
dos_Ks_180109 = zeros(1,961);
scan_errors_180109 = zeros(1,961);
%CBS_180109 = zeros(1,961);

tic;
for i =  1:961%25%round(961*rand)
	if ismember(i, bad_DPI) == 1
        disp(['Scan ',num2str(i),' is bad DPI']);
        continue
    end 
    
    try 
        %rfc_ncorr_scan
        
        %cone = result1i(:,:,i);
        cone = zeros(300,800);
        cone = cone + result1i(:,:,i);
        %{
        [s_row, s_col] = ind2sub([31,31],i);
        [s_sites] = Neighbor_sites(s_row,s_col,31,31);
        
        for ii = s_sites'
            [ss_row, ss_col] = ind2sub([31,31],ii);
            site_dist = sum(abs([s_row,s_col] - [ss_row,ss_col]));
            if site_dist == 2
                continue
            end
            cone = cone + result1i(:,:,ii);
            disp('added'),disp(num2str(ii))
        end
        %}
        ss = imgaussfilt(cone,bcb_sigma);   
        %ss(ss>(mean(ss(:))+2.5*std(ss(:)))) = mean(ss(:))+2.5*std(ss(:));
        ssK = kLOS(i);%rfc_ks_after(i);
        sss = 1000*sum( ss(ssK-ssKpm : ssK+ssKpm,:));
        smin = 399+find(sss(:,400:550)==min(sss(:,400:550)));
        
        center_EDC_mid_min = find(sss==min(sss(400:500)));
        [pks,locs] = findpeaks(sss,'MinPeakProminence',.1*max(sss));
        %CBS_180109(i) = locs(find(locs>center_EDC_mid_min,1,'first'));
        %figure, findpeaks(sss,[1:length(sss)],'MinPeakProminence',.1*max(sss),'Annotate','extents')
        
        
        y_y = sss(smin-bcbmarj:end);
        x_x = [1:length(y_y)];
        FL_param = rfc_FL_Es(i) - (smin-bcbmarj);
        
        SS_tot = sum((y_y-mean(y_y)).^2);
        %a1 = bcbmarj+10;%find(y_y(bcbmarj:bcbmarj+50) == max(y_y(bcbmarj:bcbmarj+50)));
        a2 = sss(smin);
        a3 = max(y_y)-a2;
        a4 = 0;

        a1_trials = zeros(4,length(a1_range));
        a1_Rs = zeros(1,length(a1_range));
        a1_n = 1;
        for a1 = bcbmarj + a1_range
            a0 = [a1, a2, a3, a4];
            afinal = round(fminsearch(@bcbdevsum,a0));
            afinal(1) = round((afinal(1)));

            R_squared = 1 - bcbdevsum(afinal)/SS_tot;
            a1_trials(:,a1_n) = afinal';
            a1_Rs(a1_n) = R_squared;
            a1_n = a1_n+1;
        end
        best_a1_n = find(a1_Rs==max(a1_Rs));
        afinal = a1_trials(:,best_a1_n)';
        bcb_find = smin - bcbmarj + afinal(1);

        y_fit = zeros(size(x_x));
        y_fit(1:afinal(1)) = afinal(2)^2;
        y_fit(afinal(1)+1:round(FL_param)) = afinal(2)^2+afinal(3)^2;
        y_fit(round(FL_param)+1:end) = afinal(4)^2;
        
        bcb_finds_180109(i) = bcb_find;
        %final_bcb_params(:,i) = afinal';
        final_bcb_Rsquareds_180109(i) = R_squared;
        
        dos_start = 199+find(sss(200:250)==max(sss(200:250)),1,'last');%round(pre_dos_Es(i)-150);
        dos_end = round(bcb_finds_180109(i));%min([500,round(bcb_finds(i))]);
        
        %%Now use found bcb to do DOS fit%%
        ssKs = ssK + [-10:1:10];
        kLOS_Rsquareds = zeros(size(ssKs));
        kLOS_intercepts = zeros(size(ssKs));
        kLOS_y__fit = zeros(length(ssKs),(dos_end-dos_start+1));
        kLOS_params = zeros(length(ssKs),4);
        for ssK_i = 1%:length(ssKs)
            
            %kwts = abs([1:size(ss,1)]-ssKs(ssK_i));
            kwts = abs([1:size(ss,1)]-ssK);
            kwts(kwts>dos_k_cutoff) = 0;
            DOS = sum(repmat(kwts, size(ss,2), 1)' .*imgaussfilt(cone,dos_sigma));

            %dos_start = find(DOS(200:250) == max(DOS(200:250))) + 199 + 20;
            y__y = DOS(dos_start : dos_end);% min(500,bcb_find));
            x__x = (1:length(y__y));

            b1 = round(length(x__x)/2);
            b2 = y__y(1)/b1;
            b3 = y__y(1);
            %b4 = b2;
            b4 = -b3;
            b0 = [b1 b2 b3 b4];% b5];
            bfinal =  (fminsearch(@dosdevsum,b0));

            bfinal(1) = round(abs(bfinal(1)));
            bfinal(2) = abs(bfinal(2));
            bfinal(4) = -2*bfinal(2)*bfinal(1)+bfinal(3);

            y__fit = zeros(size(x__x));
            y__fit(1:bfinal(1)) = -bfinal(2)*x__x(1:bfinal(1)) + bfinal(3);
            y__fit(bfinal(1)+1:end) = bfinal(2)*x__x(bfinal(1)+1:end) + bfinal(4);

            intercept = (bfinal(4) - bfinal(3))/(-bfinal(2) - bfinal(2));
            SS__tot = sum((y__y-mean(y__y)).^2);
            R__squared = 1 - (dosdevsum(bfinal))/SS__tot;
            
            kLOS_y__fit(ssK_i,:) = y__fit;
            kLOS_params(ssK_i,:) = bfinal;
            kLOS_Rsquareds(ssK_i) = R__squared; 
            kLOS_intercepts(ssK_i) = intercept;
            
            %{
            figure, 
            subplot(1,2,1)
            plot(x__x, y__y,'k'), hold on;
            plot(x__x, y__fit, 'r'), hold off;
            title(['E=',num2str(round(dos_E,2)),'  Rsq=',num2str(round(R__squared,4))])
            
            subplot(1,2,2)
            imagesc(ss), axis xy, hold on;
            plot([0,size(ss,2)],[kLOS(i),kLOS(i)],'w'), hold on;
            plot([0,size(ss,2)],[ssKs(ssK_i),ssKs(ssK_i)],'r'), hold off
            title(['kLOS=',num2str(ssKs(ssK_i))])
            %}
        end
        %max_kLOS_i = find(kLOS_Rsquareds==max(kLOS_Rsquareds));
        max_kLOS_i = 1;
        %dos_K = ssKs(max_kLOS_i);
        dos_K = ssK;
        %dos_E = kLOS_intercepts(max_kLOS_i) + dos_start - 1; 
        dos_E = kLOS_intercepts(1)+dos_start-1;
        
        dos_Es_180109(i) = dos_E;
        dos_Ks_180109(i) = dos_K; 
        final_dos_params(:,i) = kLOS_params(max_kLOS_i,:)';
        final_dos_Rsquareds_180109(i) = kLOS_Rsquareds(max_kLOS_i);
        
        
        
        

        if wannasee == 1
            figure,
            subplot(2,2,1), plot(kLOS_Rsquareds)
            %plot(x_x+smin-bcbmarj,y_y,'k'), hold on;
            plot([1:length(sss)],sss,'k'), hold on;
            plot(x_x+smin-bcbmarj,y_fit,'r'), hold off;
            xlim([0,size(ss,2)])
            title(['BCB=',num2str(bcb_finds_180109(i))])

    %         subplot(2,2,2),
    %         plot(x__x, y__y,'k'), hold on;
    %         plot(x__x, [kLOS_y__fit(round(length(ssKs)/2),:)], 'r'), hold off;
    %         title(['preE=',num2str(round(pre_dos_Es(i),2)),'  Rsq=',num2str(kLOS_Rsquareds(round(length(ssKs)/2)))])
    %             
            subplot(2,2,3), imagesc(ss), axis xy, hold on;
            %plot([0,size(ss,2)],[kLOS(i),kLOS(i)],'w'), hold on;
            plot([bcb_finds_180109(i),bcb_finds_180109(i)],[1,size(ss,1)],'w'), hold on;
            plot(dos_Es_180109(i), dos_Ks_180109(i), 'r*'), hold on;
            %plot([CBS_180109,CBS_180109],[1,300],'c'), hold on;
            plot([rfc_FL_Es(i),rfc_FL_Es(i)],[1,300],'r'), hold off;
            title(['DPI=',num2str(round(DPI_big(i)))])

            %plot([0,size(ss,2)],[dos_K,dos_K],'r'), hold on;

            subplot(2,2,2),
            plot(x__x, y__y,'k'), hold on;
            plot(x__x, [kLOS_y__fit(max_kLOS_i,:)], 'r'), hold off;
            title(['E=',num2str(round(dos_E,2)),'  Rsq=',num2str(round(final_dos_Rsquareds_180109(i),3))])

            subplot(2,2,4)
            imagesc(imgaussfilt(cones(:,:,i),5)), axis xy, hold on;
            plot(pre_dos_Es(i),kLOS(i),'w*'), hold on;
            plot(dos_Es_180109(i), dos_Ks_180109(i), 'r*'), hold on;
            plot([bcb_finds_180109(i),bcb_finds_180109(i)],[1,size(ss,1)],'w'), hold off;

            suptitle(['i=',num2str(i)])

            
%             figure,
%             subplot(2,2,1)
%             plot((1:length(sss)), sss, 'k'), hold on;%,(YY+i-400),'k'), hold on;
%             plot(smin,sss(smin),'r*'), hold on;
%             plot((x_x+smin-bcbmarj), y_fit, 'r'), hold off
%             xlim([1,size(ss,2)])
%             title(['Rsq=',num2str(R_squared)])
% 
%             subplot(2,2,2), plot(DOS)
%             title(['BCB=',num2str(bcb_find)])
% 
%             subplot(2,2,3)
%             %imagesc(imgaussfilt(cones(:,:,i),5)), axis xy, hold on;
%             %imagesc(imgaussfilt(cones(:,:,i),bcb_sigma)), axis xy, hold on;
%             %plot([BCB_Es(i),BCB_Es(i)],[1,size(ss,1)],'r'), hold on;
%             %plot([smin-bcbmarj+afinal(1),smin-bcbmarj+afinal(1)],[1,size(ss,1)],'w'), hold on;
%             %plot(dos_E,kLOS(i),'w*'), hold on;
%             plot([500,500],[1,size(ss,1)],'w'), hold off;
% 
%             subplot(2,2,4)
%             plot(x__x, y__y,'k'), hold on;
%             plot(x__x, y__fit, 'r'), hold off;
%             title(['E=',num2str(round(dos_E,2)),'  Rsq=',num2str(round(R__squared,2))])
% 
%             suptitle(['i=',num2str(i),'   DPI=',num2str(DPI_big(i))])
%             %figure, imagesc(ss), axis xy, title(['i=',num2str(i),' used for bcb fit'])
        end
        
    catch ME
        disp('Error occured with scan '); disp(i)
        scan_errors_180109(i) = 1;
    end
end
toc