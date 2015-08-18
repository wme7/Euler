
double = 0;
square = 0;
for i = 1:1E6
x = rand;
tic; s=x*x; double=double+toc;
tic; S=x^2; square=double+toc;
end

disp([s,S]);
disp([double/60,square/60]);