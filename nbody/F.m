function a = F(i,r,n,G,m,r0)
% computes F_i, but i here is actually i*n
    ri = r(i:i+n-1);
    N = length(r);
    a = G*m(1)*(r0-ri)./(norm(r0-ri))^3;
    mcnt = 2;
    for j=1:n:N
        if i~=j
            rj = r(j:j+n-1);
            a = a + G*m(mcnt)*(rj-ri)./(norm(rj-ri))^3;    
        end
        mcnt = mcnt + 1;
    end
end

