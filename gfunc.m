function y = gfunc(x)

y = ones(size(x))*40;
i = find(x);
y(i) = sin(40*pi*x(i))./(pi*x(i));