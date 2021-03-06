function [processed_IT] = IT_processor(input_IT, E_B1, E_B2)
%%%E_B1, E_B2 are pixel values of where the linear weights start/stop %%%

%processed_IT_ = imgaussfilt(input_IT,draw_sigma);
processed_IT_ = (input_IT);

for i = 1 : E_B1
    processed_IT_(:,i) = processed_IT_(:,i) * (1 - ((E_B1 - i) / (E_B1)));
end

E_end = size(input_IT,2);
for i = E_end-E_B2 : E_end
    processed_IT_(:,i) = processed_IT_(:,i) * (E_end-i) / (E_B2);
end
%{
for i = E_0 : E_B2
    %%%add negative weight to region inside cone between DP and beginning
    %%%of conduction band
    %processed_IT_col = processed_IT_(:,i);
    %processed_IT_col(find(processed_IT_col==0))=-mean(processed_IT(:));
    %processed_IT_(:,i) = processed_IT_col;
end
%}

processed_IT = processed_IT_;% - mean(processed_IT_(:));
end

