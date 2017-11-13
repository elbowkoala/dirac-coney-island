

function DD = bcbdevsum(a)

global x_x y_y FL_param
FL_p = round(FL_param);
a(1) = round((a(1)));
%a(2) = abs(a(2));
%a(3) = abs(a(3));
%a(4) = abs(a(4));

y_fit = zeros(size(x_x));
y_fit(1:a(1)) = a(2)^2;
y_fit(a(1)+1:(FL_p)) = a(2)^2+a(3)^2;
y_fit((FL_p)+1:end) = a(4)^2;

DD = sum((y_y - y_fit).^2);

%sum((y_y.*f1 - a(3).*f1).^2) + ...
%         sum((y_y.*f2 - a(4).*f2).^2) + ...
%         sum((y_y.*f3 - a(5).*f3).^2);

%f3 = f1;
%f1(a(1):a(2)) = 1;
%f2(a(2)+1:a(3)) = 1;
%f3(a(3)+1:end) = 1;

%DD = sum(((y_y.*f1)-a(4).*f1).^2) + sum(((y_y.*f2)-a(5).*f2).^2) + sum(((YY.*f3)-a(6).*f3).^2);

%DD = sum((YY - a(1)./((XX-a(2)).^2+a(3)) - a(4)./((XX-a(5)).^2+a(6))).^2);
%DD = sum((YY - a(1)./((XX-a(2)).^2 + a(3))).^2);% - b(1)./((XX-b(2)).^2+b(3))).^2);
end