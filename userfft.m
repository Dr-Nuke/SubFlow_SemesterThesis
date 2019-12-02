close all
% % dt = 1/100; % sampling rate
% % et = 20; % end of the interval
% % t = 0:dt:et; % sampling range
% % yy1=2;
% % yy2=3*(t-5);
% % yy3=sin(2*pi*t).*randn(1,size(yy2,2));
% % yy4=(-t.^2)/10;
% % y=yy1 +yy2+ yy3+yy4;
% % 
% % 
% % % y = 2*sin(2*pi*t); % sample the signal
% % noise = randn(1,size(y,2)); % random noise
% % ey = y + noise; % samples with noise
% % Y = fft(y); % Fourier transform of original signal
% % eY = fft(ey); % Fourier transform of noisy signal
% % fY = fix(eY/500)*500; % set numbers < 150 to zero
% % ifY = ifft(fY); % inverse Fourier transform of fixed data
% % cy = real(ifY); % remove imaginary parts
% % 
% % plot(t,y)
% % 
% % figure
% % subplot(2,1,1); % first of two plots
% % plot(t,ey); grid on % plot with grid
% % axis([0 et -10 10]); % adjust scaling
% % xlabel('Time (s)'); % time expressed in seconds
% % ylabel('Amplitude'); % amplitude as function of time
% % 
% % subplot(2,1,2); % second of two plots
% % plot(t,cy); grid on % plot with grid
% % axis([0 et -10 10]); % adjust scaling
% % xlabel('Time (s)'); % time expressed in seconds
% % ylabel('Amplitude'); % amplitude as function of time
% % 
% % figure
% % n = size(y,2)/2; % 2nd half are complex conjugates
% % plot(Y/n,'r+') % Fourier transform of original
% % hold on % put more on same plot
% % plot(eY/n,'bx') % Fourier transform of noisy signal

%x = rand(1,10); % suppose 10 samples of a random signal
x=sig;
y = fft(x); % Fourier transform of the 
iy = ifft(y); % inverse Fourier transform
x2 = real(iy); % chop off tiny imaginary parts
norm(x-x2); % compare original with inverse of transformed

%dt = 1/100; % sampling rate
dt=1/2500;
%et = 4; % end of the interval
et=1;
t = dt:dt:et; % sampling range
%y = 3*sin(4*2*pi*t) + 5*sin(2*2*pi*t); % sample the signal
%y=signal;
y=sig;

subplot(2,1,1); % first of two plots
plot(t,y); grid on % plot with grid
axis([0 et -8 8]); % adjust scaling
xlabel('Time (s)'); % time expressed in seconds
ylabel('Amplitude'); % amplitude as function of time

Y = fft(y); % compute Fourier transform
n = size(y,2)/2; % 2nd half are complex conjugates
amp_spec = abs(Y)/n; % absolute value and normalize

subplot(2,1,2); % second of two plots
freq = (0:79)/(2*n*dt); % abscissa viewing window
plot(freq,amp_spec(1:80)); grid on % plot amplitude spectrum
xlabel('Frequency (Hz)'); % 1 Herz = number of cycles/second
ylabel('Amplitude'); % amplitude as function of frequency

noise = randn(1,size(y,2)); % random noise
ey = y + noise; % samples with noise
eY = fft(ey); % Fourier transform of noisy signal
n = size(ey,2)/2; % use size for scaling
amp_spec = abs(eY)/n; % compute amplitude spectrum

figure % plots in new window
subplot(2,1,1); % first of two plots
plot(t,ey); grid on % plot noisy signal with grid
axis([0 et -8 8]); % scale axes for viewing
xlabel('Time (s)'); % time expressed in seconds
ylabel('Amplitude'); % amplitude as function of time
subplot(2,1,2); % second of two plots
freq = (0:79)/(2*n*dt); % abscissa viewing window
plot(freq,amp_spec(1:80)); grid on % plot amplitude spectrum
xlabel('Frequency (Hz)'); % 1 Herz = number of cycles per second
ylabel('Amplitude'); % amplitude as function of frequency

figure % new window for plot
plot(Y/n,'r+') % Fourier transform of original
hold on % put more on same plot
plot(eY/n,'bx') % Fourier transform of noisy signal

fY = fix(eY/100)*100; % set numbers < 100 to zero
ifY = ifft(fY); % inverse Fourier transform of fixed data
cy = real(ifY); % remove imaginary parts

figure % new window for plot
plot(t,cy); grid on % plot corrected signal
axis([0 et -8 8]); % adjust scale for viewing
xlabel('Time (s)'); % time expressed in seconds
ylabel('Amplitude'); % amplitude as function of time