




function DD = trydevsum(a)
global XX YY
a(1) = abs(a(1));
a(4) = abs(a(4));
a(3) = abs(a(3));
a(6) = abs(a(6));

DD = sum((YY - (a(1)./((XX-a(2)).^2+a(3))) - (a(4)./((XX-a(5)).^2+a(6)))).^2);
%DD = sum((YY - a(1)./((XX-a(2)).^2 + a(3))).^2);% - b(1)./((XX-b(2)).^2+b(3))).^2);
end