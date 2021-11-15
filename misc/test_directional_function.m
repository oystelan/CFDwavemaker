clear;


w = [0.1 0.1];
theta = [0 30*pi/180];
theta = [0 0]
h = 30;

%calculate wavenumber;
for i = 1:length(w),
    omega = w(i);
    W = @(k1)omega-sqrt(9.81*k1.*tanh(k1.*h));
    k(i) = abs(fzero(W,0.0));
end

% 
kx = k.*cos(theta);
ky = k.*sin(theta);


i = 1
m = 2
gamma_nm = cos(theta(i)-theta(m));
k_nm_plus = sqrt(k(i) * k(i) + k(m) * k(m) + (2. * k(i) * k(m) * gamma_nm));
k_nm_minus = sqrt(k(i) * k(i) + k(m) * k(m) - (2. * k(i) * k(m) * gamma_nm));

k_nm_plus2 = 2*k(1)
return

w_nm_plus = k_nm_plus *tanh(k_nm_plus*h)
w_nm_minus = k_nm_minus *tanh(k_nm_minus*h)

% knm alternative

kx2_min = abs(kx(i)-kx(m));
ky2_min = abs(ky(i)-ky(m));

kx2_plu = abs(kx(i)+kx(m));
ky2_plu = abs(ky(i)+ky(m));

k2_min = sqrt(kx2_min^2 + ky2_min^2);
k2_plu = sqrt(kx2_plu^2 + ky2_plu^2);


w_nm_minus2 = k2_min * tanh(k2_min*h)
w_nm_plus2 = k2_plu * tanh(k2_plu*h)

