% abs stability regions of SBDF2 and ETDRK2 for u'=aiu-b(1+i)u, b>0

f = @(a,b) ((4/3)-(4/3)*b*(1+1i)-sqrt(((4/3)*b*(1+1i)-(4/3))^2-4*(1-(2/3)*a*1i)*((1/3)-(2/3)*b*(1+1i))))/(2*(1-(2/3)*a*1i));
g = @(a,b) exp(a*1i)+((exp(a*1i)-1)*(-b*(1+1i)))/(a*1i)-(1/a^2)*(exp(a*1i)-1-a*1i)*((-b*(1+1i))*(exp(a*1i)+((exp(a*1i)-1)*(-b*(1+1i)))/(a*1i))+b*(1+1i));

for x=-20:0.2:20
    if x == 0
        continue
    end
    for y=0.01:0.2:10
        if abs(g(x,y)) <= 1
            plot(x,y,'b.'); hold on
        else
            plot(x,y,'r.'); hold on
        end
    end
end

