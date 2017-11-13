load_data = 0;
if load_data == 1
    load 'cones.mat';
    load 'ssd_big_scan_171003.mat';
    load 'rfc_big_scan_170927.mat';
    load rfc_FL_scan_170927.mat;
    load rfc_ncorr_scan_171019w.mat;
    load kLOS;
    bad_DPI = find(DPI_big < 200);
end

global x_x y_y FL_param
global x__x y__y %bcb_find dos_start
wannasee = 0;

dos_k_cutoff = 90;
ssKpm = 5;
bcbmarj = 25;
a1_range = [0:10:50];
bcb_sigma = 7;
dos_sigma = 5;


final_bcb_params = zeros(4,961);
final_bcb_Rsquareds = zeros(1,961);
final_dos_params = zeros(4,961);
final_dos_Rsquareds = zeros(1,961);
bcb_finds = zeros(1,961);
dos_Es = zeros(1,961);
scan_errors = zeros(1,961);


tic;
for i = 1:961%round(961*rand)
	if ismember(i, bad_DPI) == 1
        disp(['Scan ',num2str(i),' is bad DPI']);
        continue
    end 
    
    try 
        %rfc_ncorr_scan
        cone = cones(:,:,i);
        ss = imgaussfilt(cone,bcb_sigma);   
        ssK = kLOS(i);%rfc_ks_after(i);
        sss = 1000*sum( ss(ssK-ssKpm : ssK+ssKpm,:));
        smin = find(sss(:,400:550)==min(sss(:,400:550)));
        
        y_y = sss(smin+399-bcbmarj:end);
        x_x = [1:length(y_y)];
        FL_param = rfc_FL_Es(i) - (smin+399-bcbmarj);
        
        SS_tot = sum((y_y-mean(y_y)).^2);
        %a1 = bcbmarj+10;%find(y_y(bcbmarj:bcbmarj+50) == max(y_y(bcbmarj:bcbmarj+50)));
        a2 = sss(smin+399);
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
        bcb_find = smin + 399 - bcbmarj + afinal(1);

        y_fit = zeros(size(x_x));
        y_fit(1:afinal(1)) = afinal(2)^2;
        y_fit(afinal(1)+1:round(FL_param)) = afinal(2)^2+afinal(3)^2;
        y_fit(round(FL_param)+1:end) = afinal(4)^2;
        
        bcb_finds(i) = bcb_find;
        final_bcb_params(:,i) = afinal';
        final_bcb_Rsquareds(i) = R_squared;
        
        
        %%Now use found bcb to do DOS fit%%
        kwts = abs([1:size(ss,1)]-ssK);
        kwts(kwts>dos_k_cutoff) = 0;
        DOS = sum(repmat(kwts, size(ss,2), 1)' .*imgaussfilt(cone,dos_sigma));

        dos_start = find(DOS(200:250) == max(DOS(200:250))) + 199 + 20;
        y__y = DOS(dos_start : min([500,bcb_find]));% min(500,bcb_find));
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
        dos_E = intercept + dos_start - 1; 

        dos_Es(i) = dos_E;
        final_dos_params(:,i) = bfinal';
        final_dos_Rsquareds(i) = R__squared;
        

        if wannasee == 1
            figure,
            subplot(2,2,1)
            plot((1:length(sss)), sss, 'k'), hold on;%,(YY+i-400),'k'), hold on;
            plot(smin+399,sss(smin+399),'r*'), hold on;
            plot((x_x+smin+399-bcbmarj), y_fit, 'r'), hold off
            xlim([1,size(ss,2)])
            title(['Rsq=',num2str(R_squared)])

            subplot(2,2,2), plot(DOS)
            title(['BCB=',num2str(bcb_find)])

            subplot(2,2,3)
            %imagesc(imgaussfilt(cones(:,:,i),5)), axis xy, hold on;
            imagesc(imgaussfilt(cones(:,:,i),bcb_sigma)), axis xy, hold on;
            plot([smin+399-bcbmarj+afinal(1),smin+399-bcbmarj+afinal(1)],[1,size(ss,1)],'w'), hold on;
            plot(dos_E,kLOS(i),'w*'), hold off

            subplot(2,2,4)
            plot(x__x, y__y,'k'), hold on;
            plot(x__x, y__fit, 'r'), hold off;
            title(['E=',num2str(round(dos_E,2)),'  Rsq=',num2str(round(R__squared,2))])

            suptitle(['i=',num2str(i),'   DPI=',num2str(DPI_big(i))])
            %figure, imagesc(ss), axis xy, title(['i=',num2str(i),' used for bcb fit'])
        end
        
    catch ME
        disp('Error occured with scan '); disp(i)
        scan_errors(i) = 1;
    end
end
toc