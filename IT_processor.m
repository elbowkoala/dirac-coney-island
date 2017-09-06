function [processed_IT] = IT_processor(input_IT, draw_sigma)
processed_IT_ = imgaussfilt(input_IT,draw_sigma);
processed_IT_ = mat2gray(processed_IT_);
processed_IT = processed_IT_ ./ sum(processed_IT_(:));
end
