%% Code for forming ICAM and conventional images using SE detector and scan data
% Authors: Akshay Agarwal, Xinglin He
% Version: 1.0, January 19, 2024
%% Pixel allocation
clear;clc;
tic;
% Constants
line_period = 104602.277; % Period of one horizontal scan, computed from knowledge of pixel dwell time, number of pixels, and other manufacturer-provided timing info
pixel_index_start = 2000; % Time after the start of line that the first pixel is scanned. obtained from analysis of the scan waveform
pixel_dur_indices = 200; % Duration of one pixel in time indices. Pixel dwell time is 2 us, and sampling time is 10 ns.
no_lines = 512;
no_frames = 4; % number of complete frames in one file
filename_base = "image2_data";
no_file = 4; % number of datasets for one image

% Precompute constant values of line start indices
line_start_offsets = pixel_index_start + round((0:no_lines-1) * line_period);

% Preallocate arrays
pk_sum_image = zeros(512, 512); % This array will store the total SE pulse peak voltage at each pixel acoss all frames
no_peaks_image = zeros(512, 512); % This array will store the number of SE pulses at each pixel across all frames
frame_pk_voltages = zeros(16,512,512); % This array will store the SE pulse peak voltage at each pixel acoss in each frame separately for analysis
frame_pk_cts = zeros(15,512,512); % This array will store the number of SE pulses at each pixel acoss in each frame separately for analysis

% Loop to process datasets
for n = 1:no_file
    filename = strcat(filename_base, num2str(n, "%d"), ".mat");
    pulses = importdata(filename);
    [voltages_test,scan_test] = extraction(pulses(2, :), pulses(1, :)); % The helper function 'extraction.m' extracts the first complete image frame from the current file
    % Loop over frames
    for k = 1:no_frames
        k
        % Extract data from the start of each frame until the end
        if k > 1
            [voltages_test,scan_test] = extraction(scan_test, voltages_test); 
        end

        % Preallocate arrays for this frame
        pk_sum = zeros(512, 512); % Temporary storage for sum of SE pulse peak voltage in current frame
        no_peaks = zeros(512, 512); % Temporary storage for number of SE pulses in current frame
        
        % Loop over lines
        for j = 1:no_lines
            for i = 1:512 % looping over each pixel
                pixel_voltages = voltages_test(line_start_offsets(j)-pixel_dur_indices+(i-1)*pixel_dur_indices:line_start_offsets(j)-pixel_dur_indices+i*pixel_dur_indices); % Extracting the voltages for that pixel
                [~, pks] = peakfinder5(pixel_voltages,0.07,[],1,true,false,100e-9,0.08,50e-9); %sfunction similar to findpeaks that returns the SE pulse peak values
                if ~isempty(pks)
                    pk_sum(i, j) = sum(pks);
                    no_peaks(i, j) = numel(pks);
                end
            end
        end
        % Accumulate results for this frame
        %pk_sum_image = pk_sum_image + pk_sum;
        %no_peaks_image = no_peaks_image + no_peaks;
        frame_pk_voltages(k + 4*(n-1),:,:) = pk_sum;
        frame_pk_cts(k + 4*(n-1),:,:) = no_peaks;
    end
end
toc;
%% Forming images with all frames
c1 = 0.163; % mean voltage response from one SE
beam_current = 0.105e-12; % incident beam current
tau = 2e-6; % Pixel dwell time
lambda_rate = beam_current/1.6e-19*1e-6; % Rate of incident ions per us
lambda = beam_current/1.6e-19*tau*no_frames*no_file; % Number of incident ions per pixel
icam_image = zeros(512,512); % Array to store ICAM image
image_pk_sum = squeeze(sum(frame_pk_voltages,1)); % Summing over the peak voltages for each frame to get the total peak voltages
image_no_peaks = squeeze(sum(frame_pk_cts,1)); % Summing over number of pulses for each frame to get the total number of peaks
scaled_pk_sum_image = image_pk_sum/c1;
for i = 1:512
    for j = 1:512
        if image_no_peaks(i,j)>0
            fun_icam = @(eta) eta - scaled_pk_sum_image(i,j)/(image_no_peaks(i,j)/exp(-lambda_rate*(1-exp(-eta))*0.13) + lambda*exp(-eta)); % implementing the ICAM estimator
            icam_image(i,j) = fzero(fun_icam,0.5); % Computing ICAM estimate for current pixel by finding the root
        end
    end
end
conv_image = scaled_pk_sum_image/lambda; % Forming conventional image by scaling the sum of peak voltages
%% Plotting images
figure;
imagesc(icam_image');
colormap gray
pbaspect([1 1 1]);
caxis([0 8]);
figure;
imagesc(conv_image');
colormap gray
pbaspect([1 1 1]);
caxis([0 8]);
c = colorbar();
c.FontSize = 16;
%% Forming images with different numbers of frames
% Using different numbers of frames to make images for the analyses in
% Figures 2, 3, and 4.
icam_frame_images = zeros(16,512,512); 
conv_frame_images = zeros(16,512,512);
file_order = randperm(16); % Creating a random order of frames
%file_order = 16:-1:1;
frame_pk_voltages_rand = frame_pk_voltages(file_order,:,:); % Pulse peak voltages in the random frame order
frame_pk_cts_rand = frame_pk_cts(file_order,:,:); % Number of pulses in random frame order
for ii = 1:16 % Looping over frames (in random order) to create images
    ii
    scaled_pk_sum_image = squeeze(sum(frame_pk_voltages_rand(1:ii,:,:),1))/c1;
    image_no_peaks = squeeze(sum(frame_pk_cts_rand(1:ii,:,:),1));
    for i = 1:512
        for j = 1:512
            if image_no_peaks(i,j)>0
                fun_icam = @(eta) eta - scaled_pk_sum_image(i,j)/(image_no_peaks(i,j)/exp(-lambda_rate*(1-exp(-eta))*0.13) + lambda_rate*ii*tau/1e-6*exp(-eta));
                icam_frame_images(ii,i,j) = fzero(fun_icam,0.5);
            end
        end
    end
    conv_frame_images(ii,:,:) = scaled_pk_sum_image/(lambda_rate*ii*tau/1e-6);
end
%% Variance for images with different numbers of frames
for i = 1:16
    conv_image = conv_frame_images(i,:,:);
    conv_var(i) = var(conv_image(:));
    ml_image = icam_frame_images(i,:,:);
    mli_var(i) = var(ml_image(:));
end
%% Plotting variances
loglog((1:16)*0.11*2*6.25,conv_var.^0.5,'x','Color','red','MarkerSize',6,'LineWidth',2);
hold on;
loglog((1:16)*0.11*2*6.25,mli_var.^0.5,'o','Color','blue','MarkerSize',6,'LineWidth',2);
xlim([1 22]);
ylim([1 3]);
ax = gca;
xlabel('dose (ions/pixel)','FontSize',18);
ylabel('std. dev.','FontSize',18);
set(ax,'XTick',0:5:20,'FontSize',18);
set(ax,'YTick',0:0.5:10,'FontSize',18);
lgd = legend('conv','ICA');
lgd.Box = 'off';
grid on;
%% Plotting images with different numbers of frames
icam_image_show = squeeze(icam_frame_images(6,:,:));
conv_image_show = squeeze(conv_frame_images(6,:,:));
figure;
imagesc(conv_image_show(:,100:388)');
colormap gray
pbaspect([16 9 1]);
caxis([0 8]);
figure;
imagesc(icam_image_show(:,100:388)');
colormap gray
pbaspect([16 9 1]);
caxis([0 8]);
%% Creating and plotting histograms for images with different numbers of frames
no_bins = 256;
for ii = 1:16
    mli_image_analysis = squeeze(icam_frame_images(ii,:,:));
    conv_image_analysis = squeeze(conv_frame_images(ii,:,:));
    %conv_image_analysis = conv_image_analysis/mean(mean(conv_image_analysis))*mean(mean(mli_image_analysis));
    conv_hist(ii,:) = histcounts(conv_image_analysis,no_bins,"BinLimits",[0,10]); % Setting limits of histogram to SE yields of 0 and 10
    mli_hist(ii,:) = histcounts(mli_image_analysis,no_bins,"BinLimits",[0,10]);
    mli_std(ii) = var(mli_image_analysis(:)).^0.5;
    conv_std(ii) = var(conv_image_analysis(:)).^0.5;
end
hist_bins = linspace(0,10,256);
plot(hist_bins,conv_hist(16,:),'LineWidth',2,'Color','red');
hold on;
plot(hist_bins,mli_hist(16,:),'LineWidth',2,'Color','blue');
hold on;
plot(hist_bins,mli_hist(7,:),'--','LineWidth',2,'Color','blue');
xlim([0 8]);
ylim([0 7000]);
ax = gca;
xlabel('pixel SE yield','FontSize',18);
ylabel('# counts','FontSize',18);
set(ax,'XTick',0:1:10,'FontSize',18);
set(ax,'YTick',0:1000:8000,'FontSize',18);
lgd = legend('conv, $\lambda=22$', 'ICA, $\lambda=22$','ICA, $\lambda=9.63$','Interpreter','latex');
ldg.Location = 'NorthWest';
lgd.Box = 'off';
lgd.FontSize = 14;
%% Creating and plotting histograms and gold and silicon regions for analysis in Figure 3
no_bins = 100;
conv_hist = zeros(16,no_bins);
mli_hist = zeros(16,no_bins);
for ii = 1:16
    mli_image_analysis = squeeze(icam_frame_images(ii,50:200,50:200)); % Change first index to 300:450 for silicon region
    conv_image_analysis = squeeze(conv_frame_images(ii,50:200,50:200));
    [conv_hist(ii,:),edges,~] = histcounts(conv_image_analysis,no_bins,"BinLimits",[0,10]);
    mli_hist(ii,:) = histcounts(mli_image_analysis,no_bins,"BinLimits",[0,10]);
    mli_std(ii) = var(mli_image_analysis(:)).^0.5;
    conv_std(ii) = var(conv_image_analysis(:)).^0.5;
end
hist_bins = linspace(0,10,no_bins);
plot(edges(1:end-1) + 0.05,conv_hist(16,:),'LineStyle','-','LineWidth',2,'Color','red');
hold on;
plot(edges(1:end-1) + 0.05,mli_hist(16,:),'LineStyle','--','LineWidth',2,'Color','blue');
plot(edges(1:end-1) + 0.05,mli_hist(9,:),'LineStyle',':','LineWidth',2,'Color','blue');
xlim([0 6]);
ylim([0 2000]);
ax = gca;
xlabel('pixel SE yield','FontSize',18);
ylabel('# counts','FontSize',18);
set(ax,'XTick',0:1:10,'FontSize',18);
set(ax,'YTick',0:500:2000,'FontSize',18);
lgd = legend('conv, $\lambda=21$', 'ICA, $\lambda=21$','ICA, $\lambda=11.8$','Interpreter','latex');
ldg.Location = 'NorthWest';
lgd.Box = 'off';
lgd.FontSize = 14;
grid on;

%% SNR estimation using Thong's method for full image

im_sz=512; % Image size
sz_fit=3; % Number of points that will be fitted to extract noise-free autocorr at zero offset

test_im=squeeze(icam_image);
test_im_mean=mean(mean(test_im));
test_im=test_im-test_im_mean; % Subtracting mean, to make SNR insensitive to constant intensity shifts
noise_corr=xcorr2(test_im,test_im)/im_sz/im_sz; % autocorr of image

noise_corr = noise_corr/noise_corr(511,512); % normalizing to maximum of autocorr
[~,corr_sz]=size(noise_corr);
% next, we will fit the autocorr near 0 offset (number of fit points
% determined by sz_fit) to estimate the noise-free autocorr at zero offset.
fit_pts=noise_corr(im_sz-sz_fit:im_sz-1,im_sz); 
x_fit=im_sz-sz_fit:im_sz-1;
fit_interp=polyfit(x_fit',fit_pts,1);
x_interp=im_sz-sz_fit:im_sz;
nf_fit=fit_interp(1)*x_interp+fit_interp(2);
phi_nf=nf_fit(end); % Extracting noise free autocorr from extrapolated fit
phi=noise_corr(im_sz,im_sz);
phi_noise=phi-phi_nf;
SNR=phi_nf/phi_noise; % Calculating SNR
%% SNR and FRC estimation from two 8 frame images
% Here, we will divide up the 16 frames into random permutations of 2 sets
% of 8 frames. We will make one conventional and one ICAM image from each set, and compute the SNR
% for each image. We will also compute the fourier ring correlation
% resolution from the two half-dataset images. 
n = 100; % number of times computation will be repeated. Each repetition has a different permutation of frames used to make the half-dataset images
no_frames = 16;
snr_conv = zeros(1,n);
snr_icam = zeros(1,n);
icam_res = zeros(1,n); % To store FRC resolution for each repetition for ICAM
conv_res = zeros(1,n); % To store FRC resolution for each repetition for conv
sz_im = 512;
ft_array = linspace(-0.5,0.5,sz_im);
[ft_X, ft_Y] = meshgrid(ft_array,ft_array); % Fourier space grid
filter_array = ft_X.^2 + ft_Y.^2;
ring_radius = 0.01; % Radius of fourier space ring for computing FRC
no_rings = 0.5/ring_radius;
frc_conv = zeros(n,no_rings);
frc_mli = zeros(n,no_rings);
f_grid = (0.5:1:no_rings)*ring_radius;
for ii = 1:n
    ii
    file_order = randperm(16);
    mli_half_images = zeros(2,512,512);
    conv_half_images = zeros(2,512,512);
    file_order = randperm(16);
    frame_pk_voltages_half = frame_pk_voltages(file_order(1:no_frames),:,:);
    frame_pk_cts_half = frame_pk_cts(file_order(1:no_frames),:,:);
    % image 1
    scaled_pk_sum_image = squeeze(sum(frame_pk_voltages_half(1:round(no_frames/2),:,:),1))/c1;
    image_no_peaks = squeeze(sum(frame_pk_cts_half(1:round(no_frames/2),:,:),1));
    for i = 1:512
        for j = 1:512
            if image_no_peaks(i,j)>0
                fun_icam = @(eta) eta - scaled_pk_sum_image(i,j)/(image_no_peaks(i,j)/exp(-lambda_rate*(1-exp(-eta))*0.13) + lambda_rate*no_frames/2*tau/1e-6*exp(-eta));
                mli_half_images(1,i,j) = fzero(fun_icam,0.5);
            end
        end
    end
    conv_half_images(1,:,:) = scaled_pk_sum_image/(lambda_rate*no_frames/2*tau/1e-6);
    % image 2
    scaled_pk_sum_image = squeeze(sum(frame_pk_voltages_half(round(no_frames/2)+1:no_frames,:,:),1))/c1;
    image_no_peaks = squeeze(sum(frame_pk_cts_half(round(no_frames/2)+1:no_frames,:,:),1));
    for i = 1:512
        for j = 1:512
            if image_no_peaks(i,j)>0
                fun_icam = @(eta) eta - scaled_pk_sum_image(i,j)/(image_no_peaks(i,j)/exp(-lambda_rate*(1-exp(-eta))*0.13) + lambda_rate*no_frames/2*tau/1e-6*exp(-eta));
                mli_half_images(2,i,j) = fzero(fun_icam,0.5);
            end
        end
    end
    conv_half_images(2,:,:) = scaled_pk_sum_image/(lambda_rate*no_frames/2*tau/1e-6);

    conv_half_1_0mean = conv_half_images(1,:,:) - mean(mean(conv_half_images(1,:,:)));
    conv_half_2_0mean = conv_half_images(2,:,:) - mean(mean(conv_half_images(2,:,:)));
    mli_half_1_0mean = mli_half_images(1,:,:) - mean(mean(mli_half_images(1,:,:)));
    mli_half_2_0mean = mli_half_images(2,:,:) - mean(mean(mli_half_images(2,:,:)));
    ncc_conv = sum(sum(conv_half_1_0mean.*conv_half_2_0mean))/(sqrt(sum(sum(conv_half_1_0mean.^2)))*sqrt(sum(sum(conv_half_2_0mean.^2))));
    ncc_mli = sum(sum(mli_half_1_0mean.*mli_half_2_0mean))/(sqrt(sum(sum(mli_half_1_0mean.^2)))*sqrt(sum(sum(mli_half_2_0mean.^2))));
    snr_conv(ii) = ncc_conv/(1-ncc_conv);
    snr_icam(ii) = ncc_mli/(1-ncc_mli);
    noise_img1_conv = squeeze(conv_half_images(1,:,:));
    noise_img2_conv = squeeze(conv_half_images(2,:,:));
    ft_noise_img1_conv = fftshift(fftshift(fft2(noise_img1_conv),2),1);
    ft_noise_img2_conv = fftshift(fftshift(fft2(noise_img2_conv),2),1);
    noise_img1_mli = squeeze(mli_half_images(1,:,:));
    noise_img2_mli = squeeze(mli_half_images(2,:,:));
    ft_noise_img1_mli = fftshift(fftshift(fft2(noise_img1_mli),2),1);
    ft_noise_img2_mli = fftshift(fftshift(fft2(noise_img2_mli),2),1);
    % Computing FRC value for each ring
    for jj = 1:no_rings
        inner_rad = (jj-1)*ring_radius;
        outer_rad = jj*ring_radius;
        ft_filter = filter_array >= inner_rad.^2 & filter_array <= outer_rad.^2;
        ft_filter_shifted = fftshift(fftshift(ft_filter,2),1);
        filter_real = real(ifft2(ft_filter_shifted));
        ft_noise_img1_filtered = ft_filter.*ft_noise_img1_conv;
        ft_noise_img2_filtered = ft_filter.*ft_noise_img2_conv;
        frc_conv(ii,jj) = sum(ft_noise_img1_filtered(:).*conj(ft_noise_img2_filtered(:)))/sqrt(sum(ft_noise_img1_filtered(:).*conj(ft_noise_img1_filtered(:))).*sum(ft_noise_img2_filtered(:).*conj(ft_noise_img2_filtered(:))));
        frc_conv(ii,jj) = real(frc_conv(ii,jj));
        ft_noise_img1_filtered = ft_filter.*ft_noise_img1_mli;
        ft_noise_img2_filtered = ft_filter.*ft_noise_img2_mli;
        frc_mli(ii,jj) = sum(ft_noise_img1_filtered(:).*conj(ft_noise_img2_filtered(:)))/sqrt(sum(ft_noise_img1_filtered(:).*conj(ft_noise_img1_filtered(:))).*sum(ft_noise_img2_filtered(:).*conj(ft_noise_img2_filtered(:))));
        frc_mli(ii,jj) = real(frc_mli(ii,jj));
    end
end
% Finding resolution from FRC. Resolution is defined as the spatial freq.
% at which FRC = 1/3
for ii = 1:n
    mli_pt = InterX([2*f_grid;frc_mli(ii,:)],[2*f_grid;1/3*ones(1,no_rings)]);
    icam_res(ii) = 1/mli_pt(1)*20;
    conv_pt = InterX([2*f_grid;frc_conv(ii,:)],[2*f_grid;1/3*ones(1,no_rings)]);
    conv_res(ii) = 1/conv_pt(1)*20;
end
%% Calculating mean SNR and FRC and plotting FRC
mean_snr_conv = mean(snr_conv);
mean_snr_mli = mean(snr_icam);
mean_frc_conv = mean(frc_conv,1);
mean_frc_mli = mean(frc_mli,1);
mli_pt = InterX([2*f_grid;mean_frc_mli],[2*f_grid;1/3*ones(1,no_rings)]);
mli_res_mean = 1/mli_pt(1)*20;
conv_pt = InterX([2*f_grid;mean_frc_conv],[2*f_grid;1/3*ones(1,no_rings)]);
conv_res_mean = 1/conv_pt(1)*20;
errorbar(f_grid*2,mean_frc_conv,var(frc_conv,1).^0.5, 'Color','red','LineWidth',2);
hold on;
errorbar(f_grid*2,mean_frc_mli,var(frc_mli,1).^0.5,'Color','blue','LineWidth',2);
hold on;
plot(f_grid*2,1/3*ones(1,no_rings),'Color','black','LineWidth',2,'LineStyle','--');
xlim([0 1]);
ylim([0 1]);
ax = gca;
ax.LineWidth = 2;
%xlabel('spatial frequency (nm$^{-1}$)','FontSize',24,'Interpreter','latex');
xlabel('spatial frequency (inverse pixels)','FontSize',32,'Interpreter','latex');
ylabel('FRC','FontSize',24,'interpreter','latex');
set(ax,'Xtick',0:0.2:1,'FontSize',70);
%set(ax, 'TickLabelInterpreter', 'latex', 'XTickLabel', {'0','1/200','1/100', '3/200','1/50','1/40' ,'3/100', '7/200','1/25','9/200', '1/20'},'FontSize',18);
set(ax,'YTick',0:0.2:1,'FontSize',70);
lgd = legend('conv', 'ICA','Interpreter','latex');
ldg.Location = 'NorthWest';
lgd.Box = 'off';
lgd.FontSize = 70;
grid on;