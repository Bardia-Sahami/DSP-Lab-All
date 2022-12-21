%% 4-1-a:
data1 = imread('lena.png');
figure(1)
imshow(data1);
title('Lena');

%% 4-1-b:
data2 = im2double(data1);

%% 4-1-c:
figure(2)
imhist(data2);

%% 4-1-d:
data3 = histeq(data2);
figure(3)
imshowpair(data1, data3,'montage');
title('Lena                                     Lena - Chromatic');

%% 4-1-e:
figure(4)
subplot(1,2,1);
imhist(data2);
subplot(1,2,2);
imhist(data3);

%% 4-2-a:
data4 = imread('Image02.jpg');
figure(5)
imshow(data4);

%% 4-2-b:
data5 = imnoise(data4,"gaussian",0, 0.2);
figure(6)
imshow(data5);

%% 4-2-c:
filter1 = ones(3,3)/9;
data6 = imfilter(data5,filter1);
figure(7)
imshowpair(data5, data6,'montage');

%% 4-2-d:
filter2 = ones(5,5)/25;
data7 = imfilter(data5,filter2);
figure(8)
imshowpair(data5, data7,'montage');

%% 4-2-e:
p = 10/100;
data8 = imnoise(data4,"salt & pepper",p);
figure(9)
imshowpair(data4, data8,'montage');

%% 4-2-f:
filter3 = ones(3,3)/9;
data9 = imfilter(data8,filter3);
figure(10)
imshowpair(data8, data9,'montage');

%% 4-2-g:
load('filter.mat');
filter4 = ftrans2(Num);
figure(11)
subplot(2,1,1)
freqz(Num);
subplot(2,1,2)
freqz2(filter4);

%% 4-2-h:
data10 = imfilter(data5,filter4);
data11 = imfilter(data8,filter4);

figure(12)
subplot(2,1,1)
imshowpair(data5, data10,'montage');
subplot(2,1,2)
imshowpair(data8, data11,'montage');

%% 4-2-i:
% At the bottom of the code

%% 4-2-j:
NoisyMedIm = MedFilt(data8,3,3);
figure(13)
imshowpair(data8, NoisyMedIm,'montage');

%% 4-2-k:
% On Report

%% 4-3-a and 4-3-b:
I = imread('Image03.jpg'); % image
I1 = squeeze(I(:,:,2)); % convert 3dimentional array to 2dimentional array
figure(14)
imagesc(I1)
%%%%%%%% matlab wavelet transform %%%%%%%%%%%%%
% [LoD,HiD] = wfilters('haar','d');
% [cA,cH,cV,cD] = dwt2(I1,LoD,HiD, 'modeh','symh'); % our picture must be have 2 dimentional array

%%%%%%%%% wavelet transform by db1 %%%%%%%%%
[cA,cH,cV,cD] = dwt2(I1, 'db2'); % our picture must be have 2 dimentional array

figure(15);
subplot(2,2,1)
imagesc(cA)
colormap gray
title('Approximation')
subplot(2,2,2)
imagesc(cH)
% imshow(cH)
colormap gray
title('Horizontal')
subplot(2,2,3)
imagesc(cV)
% imshow(cV)
colormap gray
title('Vertical')
subplot(2,2,4)
imagesc(cD)
% imshow(cD)
colormap gray
title('Diagonal')



xrecD = idwt2([], [], [], cD, 'db1');
xrecHD = idwt2([], cH, [], cD, 'db1');
figure(16);
subplot(1,2,1);
imagesc(xrecD)
title('Diagonal')

subplot(1,2,2);
imagesc(xrecHD)
title('Horizontal-Diagonal')

xrecH = idwt2([], cH, [], [], 'db1');
figure(17);
subplot(1,2,1);
imshow(xrecH)

I2 = xrecH(1:size(I1,1), 1:size(I1, 2));
figure(17);
subplot(1,2,2);
imagesc(I2);
colormap(gray);



%%%%%%% combined cH,cV,cD images %%%%%%%%
combined_im = cat(3,cH,cV,cD);
figure(18);
imagesc(combined_im)
imshow(combined_im);

[cAsym,cHsym,cVsym,cDsym] = dwt2(I1,'db2','mode','sym');


%% 4-4-a:
data13  = imread('Image04.png');
filter5 = fspecial('motion',15,20*pi/180);
data14  = imfilter(data13,filter5);
figure(19)
imshow(data14);

%% 4-4-b:
for ii = 1:6
    estimated_nsr = [0, 0.001 0.01, 0.1, 1, 5];
    wnr           = deconvwnr(data14, filter5, estimated_nsr(ii));
    f = figure(20);
    f.WindowState = 'maximized';
    subplot(3,2,ii);
    imshow(wnr);
    title("Restoration of Blurred, Noisy Image Using NSR = " ...
        + num2str(estimated_nsr(ii)));
end

%% 4-4-c:
data15 = imnoise(data14,"gaussian",0, 0.01);
figure(21)
imshowpair(data14, data15,'montage');

%% 4-4-d:
for ii = 1:6
    estimated_nsr = [0, 0.001 0.01, 0.1, 1, 5];
    wnr           = deconvwnr(data15, filter5, estimated_nsr(ii));
    f = figure(22);
    f.WindowState = 'maximized';
    subplot(3,2,ii);
    imshow(wnr);
    title("Restoration of Blurred, Noisy Image Using NSR = " ...
        + num2str(estimated_nsr(ii)));
end

%% 4-5-a:
data16 = im2double(imread('glass.tif'));
figure(23)
imshow(data16);
title('glass');

%% 4-5-b:
data17 = fftshift(fft2(data16));

figure(24);
mesh(10*log10(abs(data17)));
grid on;
title('DFT Mag');

figure(25);
mesh(angle(data17));
grid on;
title('DFT phase');

%% 4-5-c:
% At the bottom of the file.

%% 4-5-d:
Output_image = FFT_LP_2D(data16,0.1*pi);

f = figure(26);
f.WindowState = 'maximized';
subplot(1,2,1);
imshow(data16); 
title('Original Glass');

subplot(1,2,2);
imshow(Output_image); 
title('Low pass Filtered Glass');

%% 4-5-e:
data18 = imresize(data16,1/4,'nearest');

f = figure(27);
f.WindowState = 'maximized';
subplot(1,3,1);
imshow(data16);
title('Original Glass');

subplot(1,3,2)
imshow(data18);
title('Down Sampled Glass with coeff. 1/4');

data19 = FFT_LP_2D(data18,0.65*pi);
subplot(1,3,3)
imshow(data19);
title('Law Passed Down Sampled Glass with 0.65pi');






%% 4-2-i:
function denoiseIm = MedFilt(img,k1,k2)
[m,n,l] = size(img) ;
k11 = floor(k1/2) ;
k22 = floor(k2/2);
denoiseIm = img;
for i = 1 : m
    if(mod(k1,2) == 0 ) || (mod(k2,2) == 0)
        disp('even filter size input in my Median Filter design')
        break;
    end
    for j = 1 : n
        for p = 1 : l
            temp = img(max(i - k11,1):min(i+k11,m) , max(j - k22,1):min(j + k22,n) , p);
            temp = reshape(temp,[],1) ;
            denoiseIm(i,j,p) = median(temp);
        end
    end
end
end

%% 4-5-c:
function Output_image = FFT_LP_2D(input_image, cutoff_frequency)
if cutoff_frequency < 0 || cutoff_frequency > pi
    disp("The frequency you entered is in false range");
else
    [ix,iy,iz] = size(input_image);
    hr = (ix-1)/2;
    hc = (iy-1)/2;
    [x, y] = meshgrid(-hc:hc, -hr:hr);

    mg = sqrt((x/hc).^2 + (y/hr).^2);
    lp = double(mg <= cutoff_frequency);

    IM = fftshift(fft2(double(input_image)));
    IP = zeros(size(IM));
    for z = 1:iz
        IP(:,:,z) = IM(:,:,z) .* lp;
    end
    Output_image = abs(ifft2(ifftshift(IP),'symmetric'));
end
end