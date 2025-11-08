# Design-of-FIR-Filters-using-hanning-window

# DESIGN OF FIR DIGITAL FILTER 

# AIM:           
  To generate design of low pass FIR digital filter using SCILAB 

# APPARATUS REQUIRED:
  PC Installed with SCILAB 

# PROGRAM 
# High Pass FIR Filter
```
clc;
clear;
close;

// High Pass FIR Filter using Hanning Window
N = 31;              // Filter length
fc = 0.3;            // Normalized cutoff frequency
n = 0:N-1;
alpha = (N-1)/2;

// Hanning window
w = 0.5 - 0.5*cos(2*%pi*n/(N-1));

// Ideal High Pass Filter Impulse Response
hd = zeros(1, N);
for i = 1:N
    if (i-1) == alpha then
        hd(i) = 1 - 2*fc;
    else
        hd(i) = -sin(2*%pi*fc*(i-1-alpha)) / (%pi*(i-1-alpha));
    end
end

// Multiply by window
h = hd .* w;

// Frequency Response
[H, f] = frmag(h, 512);

// Plot
figure;
subplot(2,1,1);
plot(f, 20*log10(abs(H)));
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
title('HIGH PASS FIR FILTER (Hanning Window)');

subplot(2,1,2);
plot(f, atan(imag(H), real(H)));
xlabel('Normalized Frequency');
ylabel('Phase (radians)');
title('Phase Response');
```
# Low Pass FIR Filter
```
clc;
clear;
close;

// Low Pass FIR Filter using Hanning Window
N = 31;              // Filter length
fc = 0.3;            // Normalized cutoff frequency (fc = Fcutoff / (Fs/2))
n = 0:N-1;
alpha = (N-1)/2;

// Hanning window
w = 0.5 - 0.5*cos(2*%pi*n/(N-1));

// Ideal Low Pass Filter Impulse Response
hd = zeros(1, N);
for i = 1:N
    if (i-1) == alpha then
        hd(i) = 2*fc;
    else
        hd(i) = sin(2*%pi*fc*(i-1-alpha)) / (%pi*(i-1-alpha));
    end
end

// Multiply by window
h = hd .* w;

// Frequency Response
[H, f] = frmag(h, 512);

// Plot
figure;
subplot(2,1,1);
plot(f, 20*log10(abs(H)));
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
title('LOW PASS FIR FILTER (Hanning Window)');

subplot(2,1,2);
plot(f, atan(imag(H), real(H)));
xlabel('Normalized Frequency');
ylabel('Phase (radians)');
title('Phase Response');
```
# OUTPUT
# High Pass FIR Filter
<img width="759" height="599" alt="image" src="https://github.com/user-attachments/assets/232a1984-28ba-4448-bd35-0de2be8bb120" />

# Low Pass FIR Filter
<img width="762" height="601" alt="image" src="https://github.com/user-attachments/assets/b5189b2f-5549-4304-9080-17977807f6f0" />

# RESULT
  Design of low pass FIR digital filter using SCILAB is generated.
