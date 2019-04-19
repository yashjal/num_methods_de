function e = expTS(x,term)
e = ones(size(x,1),1);
k = 1;
for n=1:lim
    e = e + (x.^n)./k;
    k = k*(k+1);
end

