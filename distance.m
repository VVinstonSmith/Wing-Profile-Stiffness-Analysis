function dist = distance(p1, p2, p0)

x0=p0(1);
y0=p0(2);
x1=p1(1);
y1=p1(2);
x2=p2(1);
y2=p2(2);

a = y2-y1;
b = -(x2-x1);
c = (x2-x1)*y1 - (y2-y1)*x1;

dist = (a*x0+b*y0+c)/sqrt(a^2+b^2);

end

