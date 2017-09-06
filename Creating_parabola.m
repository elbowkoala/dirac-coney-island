function [ output_matrix ] = Creating_parabola( A0,B0,E0,K0,ex_width,in_matrix )
%CREATING_PARABOLA Summary of this function goes here
%   Detailed explanation goes here
    output_matrix=zeros(size(in_matrix));
    m_center=round(size(in_matrix,1)/2);
    for i= 1:size(in_matrix,1)
        k_pos=(i-m_center-K0);
        first_intervals = [round(A0*(k_pos-1)+B0*(k_pos-1)^2+E0),round(A0*(k_pos+1)+B0*(k_pos+1)^2+E0)];
        second_intervals = [round(-A0*(k_pos-1)+B0*(k_pos-1)^2+E0),round(-A0*(k_pos+1)+B0*(k_pos+1)^2+E0)];
        
        good_cols=unique([min(first_intervals)-ex_width:max(first_intervals)+ex_width,min(second_intervals)-ex_width:max(second_intervals)+ex_width]);
        
        good_cols=good_cols(good_cols > 0 & good_cols <= size(in_matrix,2));
        output_matrix(i,good_cols)=1;
        
        
    end
    

end

