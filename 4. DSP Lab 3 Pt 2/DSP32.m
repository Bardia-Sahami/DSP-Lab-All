%% 3-2-a:
N = 300;
f1 = 4;
f2 = 8;
f3 = 12;
fs2 = 400;
t1 = 0:(1/fs2):(2-1/fs2);
t2 = 2:(1/fs2):(4-1/fs2);
t3 = 4:(1/fs2):(6-1/fs2);
x = [cos(2*pi*f1*t1),cos(2*pi*f2*t2),cos(2*pi*f3*t3)];

% H1:
B1 = [0.969531 -1.923772 0.969531];
A1 = [1 -1.923772 0.939063];

% H2:
B2 = [0.996088 -1.976468 0.996088];
A2 = [1 -1.976468 0.992177];

n = 0:N;
H1 = impz(B1,A1,n);
H2 = impz(B2,A2,n);

figure(1);
subplot(2,2,1);
plot(n,abs(H1), 'LineWidth', 1.5);
grid on;
xlabel('sample');
ylabel('|H1(z)|');
title('Impulse Response, magnitude');

subplot(2,2,2);
plot(n, abs(H2), 'LineWidth', 1.5);
grid on;
xlabel('sample');
ylabel('|H2(z)|');
title('Impulse Response, magnitude');

subplot(2,2,3);
plot(n,(180/pi)*angle(H1), 'LineWidth', 1.5);
grid on;
xlabel('freq');
ylabel('\angle H1(z)');
title('Impulse Response, angle');

subplot(2,2,4);
plot(n, (180/pi)*angle(H2), 'LineWidth', 1.5);
grid on;
xlabel('freq');
ylabel('\angle H2(z)');
title('Impulse Response, angle');

% H1:
figure(2)
freqz(B1,A1,1024,'whole');
% H2:
figure(3)
freqz(B2,A2,1024,'whole');

% Canonical Presentation
csys1  = canon(tf (B1 , A1),'companion');
csys2  = canon(tf (B1 , A1),'companion');

% y = zeros(1,N);
% for k = 1:3
%   delta_f = 4;
%   beta = tan(pi*delta_f/fs);
%   a1 = 2*x(k)/(1+beta);
%   a2 = (1-beta)/(1+beta);
%   b1 = (-2*x(k))/(1+beta) ;
%   b2 = 1/(1+beta) ;
%   for i = (1 + N/3*(j-1)) : (N/3*j)
%       if (i == 1 || i == (1 + N/3) || i == (1 + N/3*2))
%           y(i) = (1/(1+beta))*x(i) ;
%           w1 = y(1) ;
%           m1 = xx(1);
%       elseif (i == 2 || i == (2 + N/3) || i == (2 + N/3*2))
%           y(i) = -a1*w1 + (1/(1+beta))*x(i) + b1*m1 ;
%           w2 = w1;
%           w1 = y(2) ;
%           m2 = m1;
%           m1 = xx(2);
%       else
%           y(i) = -a1*w1 -a2*w2 + (1/(1+beta))*x(i) + b1*m1 + b2*m2 ;
%           w2 = w1;
%           w1 = y(i);
%           m2 = m1;
%           m1 = xx(i);
%       end
%   end  
% end

%% 3-2-b:
% On Report

% S1 = stepinfo(tf (B1 , A1),'SettlingTimeThreshold',0.01);
% st1 = S1.SettlingTime;
% st1

% figure(4)
% step(tf (B1 , A1));
% hold on;
% grid on;
% xlabel('n'); 
% ylabel('h1[n]'); 
% 
% figure(40)
% step(tf (B2 , A2));
% hold on;
% grid on;
% xlabel('n'); 
% ylabel('h2[n]'); 

figure(4);
subplot(2,1,1);
plot(n, H1, 'LineWidth', 1.5);
grid on;
xlabel('sample (n)');
ylabel('H1(z)');
title('Impulse Response');

subplot(2,1,2);
plot(n, H2, 'LineWidth', 1.5);
grid on;
xlabel('sample (n)');
ylabel('H2(z)');
title('Impulse Response');



%% 3-2-c:
y1 = filter(B1,A1,x);

figure(5);
subplot(2,1,1);
plot(linspace(0,6,length(x)), x,'r', 'LineWidth', 1.5);
grid on;
xlabel('t (sec)');
ylabel('x(t)');
title('Input Signal');

subplot(2,1,2);
plot(linspace(0,6,length(y1)), y1, 'c', 'LineWidth', 1.5);
grid on;
xlabel('t (sec)');
ylabel('y(t)');
title('Notch Filter Output, \Deltaf = 4');


y2 = filter(B2,A2,x);

figure(50);
subplot(2,1,1);
plot(linspace(0,6,length(x)), x,'r', 'LineWidth', 1.5);
grid on;
xlabel('t (sec)');
ylabel('x(t)');
title('Input Signal');

subplot(2,1,2);
plot(linspace(0,6,length(y2)), y2, 'c', 'LineWidth', 1.5);
grid on;
xlabel('t (sec)');
ylabel('y(t)');
title('Notch Filter Output, \Deltaf = 5');

%% 3-2-d:
figure(6);
[max_H1, indx] = max(H1);
plot(n,abs(H1), 'm', 'LineWidth', 1.5);
hold on;
plot(n(indx), max_H1, 'ro', 'LineWidth', 2);
grid on;
xlabel('freq');
ylabel('|H1(z)|');
title('Impulse Response, magnitude');
legend('|H1(z)|', 'max(|H1(z)|)');

figure(60);
[max_H2, indx] = max(H2);
plot(n,abs(H2), 'm', 'LineWidth', 1.5);
hold on;
plot(n(indx),max_H2, 'ro','LineWidth', 2);
grid on;
xlabel('freq');
ylabel('|H2(z)|');
title('Impulse Response, magnitude');
legend('|H2(z)|', 'max(|H2(z)|)');


%% 3-2-e:
% already presented in 3-2-c and 3-2-d for H2
y2 = filter(B2,A2,x);

figure(7);
subplot(2,1,1);
plot(linspace(0,6,length(x)), x, 'r', 'LineWidth', 1.5);
grid on;
xlabel('t (sec)');
ylabel('x(n)');
title('Input Signal');
subplot(2,1,2);
plot(linspace(0,6,length(y2)), y2, 'm', 'LineWidth', 1.5);
grid on;
xlabel('t (sec)');
ylabel('y(t)');
title('Notch Filter Output, \Deltaf = 0.5');

figure(8);
[max_H2, indx] = max(H2);
plot(n,abs(H2), 'm', 'LineWidth', 1.5);
hold on;
plot(n(indx),max_H2, 'ro','LineWidth', 2);
grid on;
xlabel('freq');
ylabel('|H2(z)|');
title('Impulse Response, magnitude');
legend('|H2(z)|', 'max(|H2(z)|)');

%% 3-2-f:
[h1, w1] = freqz(B1,A1,2048);
[h2, w2] = freqz(B2,A2,2048);

figure(9)
subplot(2,1,1);
plot(w1/pi*fs2/2,abs(h1),'c', 'LineWidth', 1.5);
xlim([0, 20]);
grid on;
xlabel('freq'); 
ylabel('|H_1(f)|'); 
title('Frequency Response of H_1(f)');
subplot(2,1,2);
plot(w2/pi*fs2/2,abs(h2), 'c', 'LineWidth', 1.5);
xlim([0, 20]);
grid on;
xlabel('freq'); 
ylabel('|H_2(f)|'); 
title('Frequency Response of H_2(f)');

%% 3-2-g:
tmp1 = find(abs(abs(h1) - 0.7*max(abs(h1))) < 0.015);
tmp2 = find(abs(abs(h2) - 0.7*max(abs(h2))) < 0.1);

figure(10)
plot(w1/pi*fs2/2, abs(h1), 'c', 'LineWidth', 1.5);
hold on;
plot(w1(tmp1(1))/pi*fs2/2, abs(h1(tmp1(1))), 'ro', 'LineWidth', 1.5);
hold on;
plot(w1(tmp1(end))/pi*fs2/2,abs(h1(tmp1(end))), 'bo', 'LineWidth', 1.5);
hold on;
plot(w2/pi*fs2/2, abs(h2), 'm--', 'LineWidth', 1.5);
hold on;
plot(w2(tmp2(1))/pi*fs2/2, abs(h2(tmp2(1))), 'go', 'LineWidth', 1.5);
hold on;
plot(w2(tmp2(end))/pi*fs2/2,abs(h2(tmp2(end))), 'ko', 'LineWidth', 1.5);
xlim([0,20]);
grid on;
xlabel('freq'); 
ylabel('|H(f)|'); 
title('Notch Filter Response');
legend('|H_1(f)|','f_L_1','f_H_1','|H_2(f)|','f_L_2','f_H_2');

%% 3-2-h:
% H3:
B3 = [0.030469 0 -0.030469];
A3 = [1 -1.923772 0.939063];
% H4:
B4 = [0.003912 0 -0.003912];
A4 = [1 -1.976468 0.992177];

n = 0:N;
H3 = impz(B3,A3,n);
H4 = impz(B4,A4,n);

figure(11);
subplot(2,2,1);
plot(n,abs(H3), 'm', 'LineWidth', 1.5);
grid on;
xlabel('sample');
ylabel('|H3(z)|');
title('Impulse Response, magnitude');

subplot(2,2,3);
plot(n,(180/pi)*angle(H3), 'm', 'LineWidth', 1.5);
grid on;
xlabel('freq');
ylabel('\angle H3(z)');
title('Impulse Response, angle');

subplot(2,2,2);
plot(n, abs(H4), 'm', 'LineWidth', 1.5);
grid on;
xlabel('sample');
ylabel('|H4(z)|');
title('Impulse Response, magnitude');

subplot(2,2,4);
plot(n, (180/pi)*angle(H4), 'm', 'LineWidth', 1.5);
grid on;
xlabel('freq');
ylabel('\angle H4(z)');
title('Impulse Response, angle');

% H3:
figure(12)
freqz(B3,A3,1024,'whole');
% H2:
figure(13)
freqz(B4,A4,1024,'whole');

% Canonical Presentation
Gp3 = tf (B3 , A3);
Gp4 = tf (B4 , A4);
GpssObs3  = canon(Gp3,'companion');
GpssObs4  = canon(Gp4,'companion');

% y = zeros(1,N);
% for k = 1:3
%   delta_f = 4;
%   beta = tan(pi*delta_f/fs);
%   a1 = 2*x(k)/(1+beta);
%   a2 = (1-beta)/(1+beta);
%   b1 = (-2*x(k))/(1+beta) ;
%   b2 = 1/(1+beta) ;
%   for i = (1 + N/3*(j-1)) : (N/3*j)
%       if (i == 1 || i == (1 + N/3) || i == (1 + N/3*2))
%           y(i) = (1/(1+beta))*x(i) ;
%           w1 = y(1) ;
%           m1 = xx(1);
%       elseif (i == 2 || i == (2 + N/3) || i == (2 + N/3*2))
%           y(i) = -a1*w1 + (1/(1+beta))*x(i) + b1*m1 ;
%           w2 = w1;
%           w1 = y(2) ;
%           m2 = m1;
%           m1 = xx(2);
%       else
%           y(i) = -a1*w1 -a2*w2 + (1/(1+beta))*x(i) + b1*m1 + b2*m2 ;
%           w2 = w1;
%           w1 = y(i);
%           m2 = m1;
%           m1 = xx(i);
%       end
%   end  
% end

% 3-2-b:
figure(14);
subplot(2,1,1);
plot(n, H3, 'LineWidth', 1.5);
grid on;
xlabel('sample (n)');
ylabel('H1(z)');
title('Impulse Response');

subplot(2,1,2);
plot(n, H4, 'LineWidth', 1.5);
grid on;
xlabel('sample (n)');
ylabel('H2(z)');
title('Impulse Response');

% 3-2-c:
y3 = filter(B3,A3,x);

figure(15);
subplot(2,1,1);
plot(linspace(0,6,length(x)), x, 'r', 'LineWidth', 1.5);
grid on;
xlabel('t (sec)');
ylabel('x(t)');
title('Input Signal');

subplot(2,1,2);
plot(linspace(0,6,length(y3)), y3, 'c', 'LineWidth', 1.5);
grid on;
xlabel('t (sec)');
ylabel('y(t)');
title('Peak Filter Output, \Deltaf = 4');

% 3-2-d:
figure(16);
[max_H3, indx] = max(H3);
plot(n,abs(H3), 'm', 'LineWidth', 1.5);
hold on;
plot(n(indx), max_H3, 'ro', 'LineWidth', 2);
grid on;
xlabel('freq');
ylabel('|H3(z)|');
title('Impulse Response, magnitude');
legend('|H3(z)|', 'max(|H3(z)|)');

% 3-2-e:
y4 = filter(B4,A4,x);

figure(17);
subplot(2,1,1);
plot(linspace(0,6,length(x)), x, 'r', 'LineWidth', 1.5);
grid on;
xlabel('t (sec)');
ylabel('x(t)');
title('Input Signal');
subplot(2,1,2);
plot(linspace(0,6,length(y4)), y4, 'm', 'LineWidth', 1.5);
grid on;
xlabel('t (sec)');
ylabel('y(t)');
title('Peak Filter Output, \Deltaf = 0.5');

figure(18);
[max_H4, indx] = max(H4);
plot(n,abs(H4), 'm', 'LineWidth', 1.5);
hold on;
plot(n(indx),max_H4, 'ro', 'LineWidth', 1.5);
grid on;
xlabel('freq');
ylabel('|H4(z)|');
title('Impulse Response, magnitude');
legend('|H4(z)|', 'max(|H4(z)|)');

% 3-2-f:
[h3, w3] = freqz(B3,A3,2048);
[h4, w4] = freqz(B4,A4,2048);

figure(19)
subplot(2,1,1);
plot(w3/pi*fs2/2,abs(h3), 'c', 'LineWidth', 1.5);
xlim([0, 20]);
grid on;
xlabel('freq'); 
ylabel('|H_3(f)|'); 
title('Frequency Response of H_3(f)');
subplot(2,1,2);
plot(w4/pi*fs2/2,abs(h4), 'c', 'LineWidth', 1.5);
xlim([0, 20]);
grid on;
xlabel('freq'); 
ylabel('|H_4(f)|'); 
title('Frequency Response of H_4(f)');

% 3-2-g:
tmp1 = find(abs(abs(h3) - 0.7*max(abs(h3))) < 0.01);
tmp2 = find(abs(abs(h4) - 0.7*max(abs(h4))) < 0.1);

figure(20)
plot(w3/pi*fs2/2, abs(h3), 'c', 'LineWidth', 1.5);
hold on;
plot(w3(tmp1(1))/pi*fs2/2, abs(h3(tmp1(1))), 'ro', 'LineWidth', 1.5);
hold on;
plot(w3(tmp1(end))/pi*fs2/2,abs(h3(tmp1(end))), 'bo', 'LineWidth', 1.5);
hold on;
plot(w4/pi*fs2/2, abs(h4), 'm--', 'LineWidth', 1.5);
hold on;
plot(w4(tmp2(1))/pi*fs2/2, abs(h4(tmp2(1))), 'go', 'LineWidth', 1.5);
hold on;
plot(w4(tmp2(end))/pi*fs2/2,abs(h4(tmp2(end))), 'ko', 'LineWidth', 1.5);
xlim([0,20]);
grid on;
xlabel('freq'); 
ylabel('|H(f)|'); 
title('Peak Filter Response');
legend('|H_3(f)|','f_L_3','f_H_3','|H_4(f)|','f_L_4','f_H_4');