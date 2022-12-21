%%3.3.ALL (B/C/D, change parameters!)
%compressor
w0 = .15*pi;
n = 1:600;
N = numel(n); % number of samples
ro = .2;
c0 = .5;
% ro = .1;
% c0 = 1.5;
% ro = 1/2;
% c0 = 1.3;
% ro = 1/4;
% ro = 1/2;
% c0 = 1.3;
lambda = .9;
c = zeros(1,N);
f = zeros(1,N);

x = zeros(1, N);
n1 = 1:200;
n2 = 200+1:400;
n3 = 400+1:600;
A1 = 2; A2 = 4; A3 = .5;
x(n1) = A1*cos(w0*n1);
x(n2) = A2*cos(w0*n2);
x(n3) = A3*cos(w0*n3);

w=c0;
for n=1:N
   c(n) = lambda*w+(1-lambda)*abs(x(n)); 
   w = c(n);
end

%%%% Compressor Gain %%%%%
for n=1:N
   if c(n) >= c0
       f(n) = (c(n)/c0)^(ro-1);
   else
       f(n) = 1;
   end
end

f = smooth(f, 7); % moving avrage function
f = transpose(f); % for scaler dot
y = f .* x; % output


figure(1);
subplot(2, 1, 1);
plot(x); 
xlabel('Time'); 
ylabel('Amplitude'); 
title('Input');

figure(1);
subplot(2, 1, 2);
plot(y); 
xlabel('Time'); 
ylabel('Amplitude'); 
title('Output');

figure(2);
subplot(2, 1, 1);
plot(c); 
xlabel('Time'); 
ylabel('Amplitude'); 
title('Control Signal');

figure(2);
subplot(2, 1, 2);
plot(f); 
xlabel('Time'); 
ylabel('Amplitude'); 
title('Compressor Gain');


%Expander
w0 = .15*pi;
n = 1:600;
N = numel(n); % number of samples
ro = 2;
% ro = 10;
c0 = .5;
lambda = .9;
c = zeros(1,N);
f = zeros(1,N);

x = zeros(1, N);
n1 = 1:200;
n2 = 200+1:400;
n3 = 400+1:600;
A1 = 2; A2 = 4; A3 = .5;
x(n1) = A1*cos(w0*n1);
x(n2) = A2*cos(w0*n2);
x(n3) = A3*cos(w0*n3);

w=c0;
for n=1:N
   c(n) = lambda*w+(1-lambda)*abs(x(n)); 
   w = c(n);
end

%%%% Expander Gain %%%%%
for n=1:N
   if c(n) < c0
       f(n) = (c(n)/c0)^(ro-1);
   else
       f(n) = 1;
   end
end

f = smooth(f, 7); % moving avrage function
f = transpose(f); % for scaler dot
y = f .* x; % output


figure(3);
subplot(2, 1, 1);
plot(x); 
xlabel('Time'); 
ylabel('Amplitude'); 
title('Input');

figure(3);
subplot(2, 1, 2);
plot(y); 
xlabel('Time'); 
ylabel('Amplitude'); 
title('Output');

figure(4);
subplot(2, 1, 1);
plot(c); 
xlabel('Time'); 
ylabel('Amplitude'); 
title('Control Signal');

figure(4);
subplot(2, 1, 2);
plot(f); 
xlabel('Time'); 
ylabel('Amplitude'); 
title('Expander Gain');