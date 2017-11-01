%%=======================Extraction===========================
% reason:
% the code is for data extraction from dsp and decode, then extract features from audio data 
% It is used to applying my signal processing knowledge to practical programming
txtname = '427';
fname = [txtname,'.txt'];
fid = fopen( fname );
Data = fscanf(fid,'%x');
fclose(fid);

% ----------------------------------
% 
fs = 16e3;
len = length(Data);
n = 128; % 128 bands
L = floor(length(Data)/(n+1)/2); 
s_tmp = zeros(1,L,'double');
s1 = zeros( n, L, 'double');

% 4.20 for mapping 0~2 to 0~1
mx = 0;
mi = 1;

counts = 0;
ct = 0;
i = 0;
for k = 1:2:len
    
    if k+1>len 
        break
    end
    %¼ì²â¿ªÍ·
    if ((Data(k)==170) && (Data(k+1)==170)) || (ct == 128)
        k = k+2;          
        flag = 0;
        counts = [counts ct];
        ct = 0;
        i = i+1;
        continue
    end
    
    temp = (bitshift(Data(k),8) + Data(k+1));
    % 4.17 for mapping from 0~2 to 0~255
    % s_tmp(i) = u16int2int16(temp)/16384*255/2;
    % 4.18 for mapping 0~2 to 0~1 to 0~255
    ct = ct+1;
    s_tmp(ct) = u16int2int16(temp)/262144; % shift 4, all positive for func
    s1(129-ct,i) = s_tmp(ct);
    if s_tmp(ct) > mx
        mx = s_tmp(ct);
    end
    if s_tmp(ct) ~=0
        if s_tmp(ct) < mi
            mi = s_tmp(ct);
        end
    end
        
   
    k
    disp('processing...')
end

%s1(s1==0) = mi/10;
offset = abs(10*log10(mi/1000))*1.1;
%--------------------Figuaration-------------------
% 4.18 for mapping 0~2 to 0~1 to dB
s2 = s1;
s2(s2==0) = mi/1000;
s2 = 10*log10(s2);
s2 = s2 + offset;
figure
imagesc(s2)
colorbar
title('Original WOLA')

% 
counts = counts(3:length(counts));
figure
plot(counts)

%=====================Spectrum Subtraction==================
% s1 WOLA signal E^2
% plot(s1)
% noise around before 150 (0.6 s)
% Spectrum subtraction

E_noise = sum(s1(:,1:50),2) / 150;
mag_s  = s1 - E_noise;
mag_s( mag_s<0 ) = 0;
%---------------------Figuration------------------
noise = E_noise;
noise(noise==0) = 1e-8;
noise = 10*log10(noise);
noise = noise+offset;
figure
plot(noise)
title('estimated average Noise')
view(90,90)

s3 = mag_s;
s3(s3==0) = mi/1000;
s3 = 10*log10(s3);
s3 = s3 + offset;
figure
imagesc(s3)
colorbar
title('Processed WOLA')

s1 = flipud(s1);
mag_s = flipud(mag_s);
%=====================Spectral flux================
flux_o = spectralFlux(s1, n, L);
%flux_o = flux_o / max(flux_o);
figure
plot(flux_o)
title('Spectral Flux (original)')

flux_p = spectralFlux(mag_s, n, L);
%flux_p = flux_p / max(flux_p);
figure
plot(flux_p)
title('Spectral Flux (processed)')

%=====================KL Distance===================
%-----------------KL for original WOLA--------------
for i = 2:L
    [KLD_s(i), KLD_v(:,i)] = KLdistance( sqrt(s1(:,i)), sqrt(s1(:,(i-1))) );
end
figure
plot(KLD_s)
title('KL-Distance (scalar)')
figure
plot(KLD_v)
title('KL-Distance (vector)')

%=====================Centroid======================
%----------------original WOLA centroid-------------
Freq = zeros(n,1,'double'); % remove DC
hop = fs/n;
for i = 1:n
    Freq(i) = i*hop - hop/2;
end

[C, C_normal] = spectralCentroid(sqrt(s1), Freq, fs);

figure
plot(C)
title('Spectral Centroid')
ylabel('Hz')
figure
plot(C_normal)
title('normalized Spectral Centroid')

%=====================Energy=======================
Eng_o = zeros(1,L);
for j = 1:L
    Eng_o(j) = sum(s1(:,j));
end
figure
plot(Eng_o)
title('Energy (original)')
Eng_o = 10*log10(Eng_o);
Eng_o = Eng_o + offset;
figure
plot(Eng_o)
title('Energy in dB (original)')

Eng_p = zeros(1,L);
for j = 1:L
    Eng_p(j) = sum(mag_s(:,j));
end
figure
plot(Eng_p)
title('Energy (processed)')
Eng_p = 10*log10(Eng_p);
Eng_p = Eng_p + offset;
figure
plot(Eng_p)
title('Energy in dB (processed)')

%%====================RECOVER========================

x = WOLA_recover(s1, fs, n, L);
y = WOLA_recover(mag_s, fs, n, L);

%%===========4.25 UPDATE=============================
%------------normalized spectral flux----------------
flux_o_new = spectralFlux_new(s1, n, L);
flux_o_new = flux_o_new / max(flux_o_new);
figure
plot(flux_o_new)
title('Spectral Flux new (original)')

flux_p_new = spectralFlux_new(mag_s, n, L);
flux_p_new = flux_p_new / max(flux_p_new);
figure
plot(flux_p_new)
title('Spectral Flux new (processed)')

%-----------reduced centroid for original-------------
for i = 4:4:n
    Freq_new(i/4) = Freq(i);
    s1_new(i/4,:) = s1(i,:);
    mag_s_new(i/4,:) = mag_s(i,:);
end
Freq_new  = Freq_new';
[C_new, C_normal_new] = spectralCentroid(sqrt(s1_new), Freq_new, fs);
[C_new, C_normal_new] = spectralCentroid(sqrt(mag_s_new), Freq_new, fs);

figure
plot(C_new)
title('Spectral Centroid new')
ylabel('Hz')
figure
plot(C_normal_new)
title('normalized Spectral Centroid new')

