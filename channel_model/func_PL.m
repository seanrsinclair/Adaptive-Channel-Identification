function [ PL ] = func_PL( d, d0, d1, L )
%func_PL: calculate Path Loss in dB
[m,k] = size(d);
for i = m:-1:1
    for j = k:-1:1
        if (d(i,j) > d1)
            PL(i,j) = -L-35*log10(d(i,j));
        elseif (d(i,j) > d0)
            PL(i,j) = -L-15*log10(d1)-20*log10(d(i,j));
        else
            PL(i,j) = -L-15*log10(d1)-20*log10(d0);
        end
    end
end
end

