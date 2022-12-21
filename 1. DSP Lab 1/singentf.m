function y =singentf(w0,n) 
    delta = zeros(1,n); delta(1) = 1;
    a = [1, -2*cos(w0), 1];
    b = [0, sin(w0)];
    y = filter(b,a,delta);
end