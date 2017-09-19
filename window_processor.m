function processed_window = window_processor(input_window)
processed_window = input_window ./ sum(input_window(:));
%processed_window = processed_window_ - mean(processed_window_);

end