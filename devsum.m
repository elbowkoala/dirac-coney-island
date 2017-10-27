
a3 = ((max(XX) - min(XX))/10)^2;
a2 = (max(XX)+min(XX))/2;
a1 = max(YY)*a3;
a0 = [a1, a2, a3];
afinal = fminsearch(@devsum,a0);

function DD = devsum(a)
global XX YY
DD = sum((YY - a(1))./((XX-a(2)).^2 + a(3)).^2);
end