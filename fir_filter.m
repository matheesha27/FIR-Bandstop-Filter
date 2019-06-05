%FIR Bandstop Filter
 
close all;
clear all;
 
indexNo = 150009;
%Get A, B, C values from the index number
A = fix(mod(indexNo,1000)/100);
B = fix(mod(indexNo,100)/10);
C = mod(indexNo,10);
 
%Define filter specifications
Ap = 0.05 + 0.01*A;     %Passband ripple
Aa = 40 + B;            %Minimum stopband attenuation
Omega_p1 = C*100 + 300; %Lower passband frequency
Omega_p2 = C*100 + 850; %Upper passband frequency
Omega_a1 = C*100 + 400; %Lower stopband frequency
Omega_a2 = C*100 + 700; %Upper stopband frequency
Omega_s = 2*(C*100 + 1200); %Sampling frequency
T = 2*pi/Omega_s;
 
Bt = min((Omega_a1 - Omega_p1),(Omega_p2 - Omega_a2));  %Transition width
Omega_c1 = Omega_p1 + Bt*0.5;   %Lower cut-off frequency
Omega_c2 = Omega_p2 - Bt*0.5;   %Upper cut-off frequency
 
%Choose delta value
delta_p = (10^(0.05*Ap)-1)/(10^(0.05*Ap)+1);
delta_a = 10^(-0.05*Aa);
delta = min(delta_p, delta_a);
 
%Actual Stopband attenuation
Aaa = -20*log10(delta);
 
%Choose parameter alpha
if Aaa<=21
    alpha = 0;
elseif (21<=Aaa) && (Aaa<=50)
    alpha = 0.5842*(Aaa-21)^0.4 + 0.07886*(Aaa-21);
elseif Aaa>50
    alpha = 0.1102*(Aaa-8.7);
end
 
%Choose parameter D
if Aaa<=21
    D = 0.9222;
else
    D = (Aaa-7.95)/14.36;
end
 
%Choose the lowest odd value of N
if mod(ceil((Omega_s*D)/Bt+1),2) == 0;
    N = ceil((Omega_s*D)/Bt+1) + 1;
else
    N = ceil((Omega_s*D)/Bt+1);
end
 
%Plot the Window function (wk) from Kaiser window function
nr = -(N-1)/2 : 1 : (N-1)/2;     %define the range where wk is non-zero.
beta = alpha*(1 - ((2*nr)/(N-1)).^2).^0.5;
I_beta = 1; I_alpha = 1;
for k = 1 : 1 : 100
    I_beta = I_beta + ((1/factorial(k))*(beta/2).^k).^2;
    I_alpha = I_alpha + ((1/factorial(k))*(alpha/2).^k).^2;
end
wk = I_beta./I_alpha;
figure;
stem(nr,wk);
title('Windowing function');
xlabel('n'); ylabel('w[n]');
grid on;
 
%Compute h[n]
n1 = -(N-1)/2 : 1 : -1; %Range for negative values
n2 = 1 : 1 : (N-1)/2;   %Range for positive values
h1 = (((1/pi)./n1).*(sin(Omega_c1*n1*T) - sin(Omega_c2*n1*T)));
h2 = (((1/pi)./n2).*(sin(Omega_c1*n2*T) - sin(Omega_c2*n2*T)));
h0 = 1 + (2*(Omega_c1-Omega_c2))/Omega_s;
n = [n1, 0, n2];
hn = [h1, h0, h2];      %h[n] array
figure;
stem(n,hn);
title('Impulse Response of Ideal Bandstop Filter');
xlabel('n'); ylabel('h[n]');
grid on;
 
%Compute Digital filter
filt = hn.*wk;
 
figure;
stem(n, filt);
title('Impulse Response of Non-causal Filter');
xlabel('n'); ylabel('H');
grid on;
 
figure;
n_new = 0:1:(N-1);
stem(n_new, filt);
title('Impulse Response of Causal Filter');
xlabel('n'); ylabel('H');
grid on;
 
%Magnitude Response
fvtool(filt)
 
%compute the Omega values and plot the Excitation in time domain
w1 = Omega_p1/2;
w2 = (Omega_a1+Omega_a2)/2;
w3 = (Omega_p2+Omega_s/2)/2;
 
ns = 0:1:300;   %No. of samples
xnT = sin(w1*ns*T)+sin(w2*ns*T)+sin(w3*ns*T); %Excitation function
figure;
plot(ns, xnT);
title('Excitation');
xlabel('time(s)');
ylabel('x(nT)');
 
y = conv2(xnT,filt);
figure;
plot([1:length(y)]*T*(length(xnT))/(length(y)),y);
title('Filtered Signal');
xlabel('time(s)');
 
N_FFT = 2^nextpow2(numel(ns));   %Next power of 2 from length of y
xnT_FFT = fft (xnT, N_FFT)/numel(ns);
f = (Omega_s)/2* linspace(0, 1, N_FFT/2+1);
figure;
plot(f, 2*abs(xnT_FFT(1:N_FFT/2+1)));
title('Discrete Fourier Transform Excitation');
xlabel('Frequency(rad/s)');
 
N_FFT = 2^nextpow2(numel(ns));   %Next power of 2 from length of y
Y_FFT = fft(y, N_FFT)/numel(ns);
f = (Omega_s)/2*linspace(0, 1, N_FFT/2+1);
figure;
plot(f, 2*abs(Y_FFT(1:N_FFT/2+1)));
title('Discrete Fourier Transform Filtered from Kaiser Window Bandstop Filter');
xlabel('Frequency(rad/s)');
 
%Ideal Filter
xnT=sin(w1*ns*T)+sin(w3*ns*T);
Y=conv2(xnT,filt);
N_FFT = 2^nextpow2(numel(ns));   %Next power of 2 from length of y
Y_FFT = fft(Y, N_FFT)/numel(ns);
f = (Omega_s)/2*linspace(0, 1, N_FFT/2+1);
figure;
plot(f, 2*abs(Y_FFT(1:N_FFT/2+1)));
title('DFT Filtered from Ideal Bandstop Filter');
xlabel('Frequency(rad/s)');
 
%%%%%%%%%%%%%end%%%%%%%%%%%%%

