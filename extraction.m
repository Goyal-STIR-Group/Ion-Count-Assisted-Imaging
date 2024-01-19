function [voltages_test, scan_test]= extraction(scan,voltages)
% Function for extracting the start of image frames
% Authors: Akshay Agarwal, Xinglin He
% Version: 1.0, January 19, 2024
end_index = numel(scan);
scan_ds = downsample(scan, 1e4); % Downsampling the scan signal to make searching for the start of image frame faster   
for i = 1:numel(scan_ds)
    if(scan_ds(i:i+80)<-4.69) % Searching for a long, flat region in the scan signal, which corresponds to the end of the previous frame
        break;
    end
end
index_frame_search_start = (i-1)*1e4 + 1 + 8e5; % Upsampling the index to the original scan voltage
frame_search_array = scan(index_frame_search_start:index_frame_search_start + 3e5); % 3e5 chosen to guarantee that this chunk contains the start of the frame
frame_start_index = find(frame_search_array>-4.69,1) + index_frame_search_start-2; 
voltages_test = voltages(frame_start_index:end_index);
scan_test = scan(frame_start_index:end_index);
end