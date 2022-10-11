function p =  perimeter_of_an_ellipse(semi_major_axis, semi_minor_axis)
% A method for computing the perimeter of an ellipse.
% Essentially, the problem of calculating the distance from the equator to the pole is to calculate the quadrants and the perimeter of an ellipse.
% There is no unique solution to this, so here are several formulas that can be used. 
% Two solutions with infinite series are the most accurate.
% Vectorized.
% Version: 1.00 13 Oct 22
% Input:   a - ellipse major semi-axis [meter]
%          b - ellipse minor semi-axis [meter]
% Output:  p - ellipse perimeter [meter]
% Reference: Kawase, Kazushige. "A general formula for calculating meridian
% arc length and its application to coordinate conversion in the Gauss-Krüger projection." Bulletin of the Geospatial Information Authority of Japan 59 (2011): 1-13. https://yadda.icm.edu.pl/baztech/element/bwmeta1.element.baztech-3b526150-d9d6-4f81-a65a-b1bce37472ab/c/Weintrit.pdf
%
% Keywords: geometric geodesy, perimeter ellipse, map projection,
% coordinate conversions
%
% Copyright (c) 2022, Matej Varga
% All rights reserved.
% Email: mvarga1989@gmail.com

a= semi_major_axis;
b= semi_minor_axis;
e= sqrt((a.^2-b.^2)./a.^2);
h=(a-b).^2./(a+b).^2;

% approximation formulas
p.app1= 2*pi*sqrt((a.^2+b.^2)./2);
p.app2 =pi.*sqrt(2*(a.^2+b.^2)-(a-b).^2/2);
p.app3= pi.*(3/2*(a+b)-sqrt(a.*b));
p.app4= pi.*[3.*(a+b)-sqrt((3.*a+b).*(a+3.*b))];
p.app5 = pi*(a+b)*(1+h/8).^2;
p.app6 = pi*(a+b)*[1+3*h./(10+sqrt(4-3.*h))];
p.app7 = pi*(a+b)*(64-3*h.^2)./(64-16.*h);
p.app8 = pi*(a+b)*(256-48*h-21.*h.^2)./(256-112.*h+3.*h.^2);
p.app9 = pi*(a+b)*(3-sqrt(1-h))./2;

% exact formulas
i=[1:30];
supp.inf1_ser = 1-sum((factorial(2.*i)).^2.* e.^(2.*i)./(2.^i.*factorial(i)).^4./(2.*i-1));
p.inf1  = 2.*a*pi.*supp.inf1_ser;

%
n=[0:50];
for ind=1:numel(n);
supp.inf2_ser(ind) = [0.5.^n(ind)./factorial(n(ind))].^2.*h.^n(ind);
end

p.inf2  = pi.*(a+b).*sum(supp.inf2_ser);

end
