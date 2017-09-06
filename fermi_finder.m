%max_L = zeros(1,num_scans);
max_R = zeros(1,num_scans);
%min_C = zeros(1,num_scans);
fermi_= zeros(1,num_scans);
%bbbb = zeros(1,num_scans);
%cond_band_beg = zeros(1,num_scans);

for i = 1:num_scans;

    duhh = sum(imgaussfilt(cones(:,:,i),5),1);
    duhhh = smooth(duhh,20,'rlowess');
    duhhh = duhhh./max(duhhh);

    if rem(i,100)==0;
        disp(i)
    end
    %max_L(i) = find(duhh==max(duhh(50:150)));
    %max_R(i) = find(duhhh==max(duhhh(600:700)),1,'last');
    fermi_(i) = find(duhhh(600:700)>=0.5,1,'last') + 599;
    %fermi(i) = find(duhh(max_R(i):end)<=(duhh(max_R(i))/2),1) + max_R(i) -1;
    %min_C(i) = find(duhh==min(duhh(max_L(i):max_R(i))));
    %cond_band_beg(i) = find(duhh(min_C(i):max_R(i))>duhh(max_L(i)),1,'first') + min_C(i) - 1;
end

figure, 
plot(fermi_)
%subplot(3,1,1), plot(duhmaxAs,'r');
%subplot(3,1,2), plot(duhmaxBs,'b');
%subplot(3,1,3), plot(duhminCs,'p');