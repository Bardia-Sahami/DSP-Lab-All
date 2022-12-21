function a = normalFreq(Ls)
    fnormal=-(Ls-1)/2:(Ls-1)/2;
%     a=fftshift(fnormal);
%     a( a >= .5 ) = a( a >= .5 )-1;
    a=fnormal;
end