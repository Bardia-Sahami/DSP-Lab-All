%part 1:
t = linspace(0, 2, 201); %or t = 0:0.01:2
A = 5; f = 1; 
y = A*sin(2*pi*f*t);
figure(1);
subplot(1,1,1);
plot(t, y)
title('Sin Function');
xlabel('Time (Second)');
ylabel('Sin(2?ft)')

%part 2:
t1 = 0:0.01:2;
A = 5; f = 1; 
y1 = A*sin(2*pi*f*t1);
t2 = -0.5:0.005:0.5;
r = 0.5 + (0.5+0.5).*rand(1,201); %random
y2 = y1 + r; %sin with noise
figure(1);
subplot(2,1,1);
plot(t1,y1)
title('Normal Sin');xlabel('Time (Second)');ylabel('Sin(2?ft)');
subplot(2,1,2);
plot(t2,y2)
hold on;
title('Noisy Sin');xlabel('Time (Second)');ylabel('Sin(2?ft)');

%part 3:
t1 = 0:0.01:2;
A = 5; f = 1; 
y1 = A*sin(2*pi*f*t1);
t2 = -0.5:0.005:0.5;
r = 0.5 + (0.5+0.5).*rand(1,201); %random
y2 = y1 + r; %sin with noise
N1 = 0; N2 = 20;
len = N2 - N1 + 1; %window length
window = ones(1,len)/len;
y3 = conv(y1, window);
t3 = 0:0.01:(0.01*(length(y3)-1));
figure(2);
plot(t3, y3)

%part 4:
t1 = 0:0.01:2;
A = 5; f = 1; 
y1 = A*sin(2*pi*f*t1);
t2 = -0.5:0.005:0.5;
r = 0.5 + (0.5+0.5).*rand(1,201); %random
y2 = y1 + r; %sin with noise
N1 = 0; N2 = 20;
len = N2 - N1 + 1; %window length
b = ones(1,len)/len;
a = 1;
y4 = filter(b, a, y2);
t4 = 0:0.01:(0.01*(length(y4)-1));
figure(3);
plot(t4, y4)

%part 5:
w0 = pi/4; n = 101;
[y5,t5] = singen(w0,n);
figure(4);
stem(t5,y5, '.')

%part 6:
t6 = 0:0.01:4;
y6 = cos(2*pi*t6) + cos(8*pi*t6) + cos(12*pi*t6);
figure(5);
plot(t6,y6)
t7 = 0:0.2:4; %fs = 5 so Ts = 0.2
y7 = cos(2*pi*t7) + cos(8*pi*t7) + cos(12*pi*t7);
hold on;
stem(t7,y7)
b = zeros (1, 21); b(1) = 1;
a = zeros (1, 21); a(1) = 1;
y8 = filter(b, a, y7);
t8 = 0:0.01:(0.01*(length(y8)-1));
figure(6);
plot(t8,y8)


%part 8:
t00 = -5:0.01:5;
y00 = (sinc(5*t00)).^2;
figure(7);
subplot(3,2,[1, 2]);
plot(t00,y00)
title('Main');

t01 = -5:0.25:5;
y01 = (sinc(5*t01)).^2;
subplot(3,2,3);
plot(t01,y01)
title('4 Hz');

t02 = -5:0.2:5;
y02 = (sinc(5*t02)).^2;
subplot(3,2,4);
plot(t02,y02)
title('5 Hz');

t03 = -5:0.1:5;
y03 = (sinc(5*t03)).^2;
subplot(3,2,5);
plot(t03,y03)
title('10 Hz');

t04 = -5:0.05:5;
y04 = (sinc(5*t04)).^2;
subplot(3,2,6);
plot(t04,y04)
title('20 Hz');

%part 9:
t9 = -1:0.01:1;
f = [pi/16 5*pi/16 9*pi/16 13*pi/16];
y9 = zeros (1, length(t9));
for fs=f
    x = cos(2*fs*t9);
    y9 = y9 + x;
end
figure(8);
L = length(y9);
y9_sp = abs(fftshift(fft(y9))/L);
subplot(2,1,1);
plot(t9,y9)
title('Sin Signal'); xlabel('Time');ylabel('Spectrum');
subplot(2,1,2);
plot(normalFreq(L),y9_sp)
title('Spectrum'); xlabel('Freq');ylabel('Spectrum');
