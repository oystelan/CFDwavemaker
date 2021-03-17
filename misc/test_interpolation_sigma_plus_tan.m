
clear;
depth = 100;

eta = 7;
z = -depth:1:eta+1;
val = exp(0.1*z);

subplot(1,2,1)
plot(val,z)

nl = 15;

s = -1:1/(nl-1):0;



a = 7*pi/18;
b = 1.0;


for i=1:nl,   
    s_mark(i)= (-(tan(((-s(i))*a))^b)/(tan(a)^b));
end

s2z = eta + s_mark .* (depth + eta);

subplot(1,2,1)
hold on;
val2 = exp(0.1.*s2z);
plot(val2,s2z,'x')
subplot(1,2,2)
plot(val2,s,'x-');

z2 = -65:10:5;
% Transform z2 coordinates to tan/sigma coordinates
s2 = (z2-eta)./(eta+depth);
for i=1:length(s2),   
    s2_mark(i) = - atan((-s2(i)*(tan(a)^b))^(1/b))./a;
end

val3 = interp1(s,val2,s2_mark);
subplot(1,2,1)
plot(val3,z2,'o')
subplot(1,2,2)
hold on;
plot(val3,s2_mark,'o')
