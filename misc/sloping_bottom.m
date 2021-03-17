clear;
dw = 0.05;
hs = 10;
tp = 15;
w = [0.45:dw:0.7];
S = jonsnor(w./(2*pi),hs,tp)./(2*pi);
a = sqrt(2*dw*S);
d = 50;
g = 9.81;
k = w.^2./g;
alpha = 0.1/3;
k = k(1);

i = 1;
for z = -d:10:0,
    p(i) = real((exp(-k.*(z+2*d))+((1-1j*alpha)/(1+1j*alpha)).*exp(k.*z))./...
        (1+exp(-2.*k.*d)));
    
    p_dz(i) = real((-k.*exp(-k.*(z+2*d))+((1-1j*alpha)/(1+1j*alpha)).*k.*exp(k.*z))./...
        (1+exp(-2.*k.*d)));
    zz(i) = z;
    i = i + 1;
    
end

subplot(1,2,1)
plot(p,zz)
subplot(1,2,2)
plot(p_dz,zz)