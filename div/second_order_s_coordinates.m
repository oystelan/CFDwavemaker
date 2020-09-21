clear;
dw = 0.01;
hs = 10;
tp = 15;
w = [0.45:dw:0.7];
S = jonsnor(w./(2*pi),hs,tp)./(2*pi);

a = sqrt(2*dw*S);
ampl_target = 13;
a = ampl_target*a./sum(a);



% sum(a)
% return
h = 300.;
g = 9.81;
k = w.^2./g;

t = 2.;
x = 0.;

signk1 = [1 -1 1 -1];
signk2 = [1 -1 -1 1];

for i = 1:length(k),
    eta1(i) = a(i)*exp(1i*(k(i)*x-w(i)*t));
    for j = 1:length(k),
        for s = 1:4,
            w1 = signk1(s)*w(i);
            w2 = signk2(s)*w(j);
            k1 = signk1(s)*k(i);
            k2 = signk2(s)*k(j);
            
            if (k1 + k2) == 0,
                eta2_temp(s) = 0.;
            else
                %                                 w12 = w(i) + w(j);
                w12 = w1 + w2;
                sigma12 = sqrt(g*abs(k1+k2)*tanh(abs(k1+k2)*h));
                %                 sigma12 = sqrt(g*abs(k1)+abs(k2)*tanh(abs(k1)+abs(k2)*h));
                
%                 chs12 = cosh(abs(k1+k2)*h*(1+s))
                h12(s) = (1i*w12*g*g)/(sigma12^2-w12^2)*...
                    ((k1*k2)/(w1*w2)-(w1*w2 + w1*w1 + w2*w2)/(2*(g*g))+...
                    (k1*k1*w2+k2*k2*w1)/(2*w1*w2*w12));
                
                %                 h12(s) = (1i*w12*g*g)/(sigma12^2-w12^2);
                
                c12(s) = ((1i*w12)/g)*h12(s)+(w1*w1+w2*w2+w1*w2)/(2*g)-...
                    g*((k1*k2)/(2*w1*w2));
                
                e1(s) = exp(1i*(k1*x-w1*t));
                e2(s) = exp(1i*(k2*x-w2*t));
                
                a1(s) = a(i)/2;
                a2(s) = a(j)/2;
                
                eta2_temp(s) = a1(s)*a2(s)*c12(s)*e1(s)*e2(s);
                
                
            end
            
        end
        %         h12
        %         return
        %             return
        eta2_sum(i,j) = real(sum(eta2_temp(1:2)));
        eta2_diff(i,j) = real(sum(eta2_temp(3:4)));
        %         return
    end
end
% figure('Position',[100 100 900 400])
% subplot(1,2,1)
% surf(k,k,eta2_sum)
% view(0,90),xlabel('k1'),ylabel('k2')
% colorbar;
% subplot(1,2,2)
% surf(k,k,eta2_diff)
% view(0,90),xlabel('k1'),ylabel('k2')
% colorbar;
eta2_sum = sum(sum(eta2_sum));
eta2_diff = sum(sum(eta2_diff));
eta2 = sum(sum(eta2_sum+eta2_diff));

eta = sum(real(eta1))+eta2;

z = -30;

ss = ((h+z)/(h+eta))-1;
for i = 1:length(k),
    shs1 = sinh(abs(k(i))*h*(1+ss))/cosh(abs(k(i))*h);
    chs1 = cosh(abs(k(i))*h*(1+ss))/cosh(abs(k(i))*h);
    
    u1(i) = a(i)*((g*k(i))/w(i))*chs1*exp(1i*(k(i)*x-w(i)*t));
    ww1(i) = a(i)*(-1i*(g*k(i))/w(i))*shs1*exp(1i*(k(i)*x-w(i)*t));
    for j = 1:length(k),
        for s = 1:4,
            w1 = signk1(s)*w(i);
            w2 = signk2(s)*w(j);
            k1 = signk1(s)*k(i);
            k2 = signk2(s)*k(j);
            
            shs2 = sinh(abs(k2)*h*(1+ss))/cosh(abs(k2)*h);
            chs2 = cosh(abs(k2)*h*(1+ss))/cosh(abs(k2)*h);
            
            if (k1 + k2) == 0,
                u2_temp(s) = 0.;
                w2_temp(s) = 0.;
            else
                %                                 w12 = w(i) + w(j);
                w12 = w1 + w2;
                sigma12 = sqrt(g*abs(k1+k2)*tanh(abs(k1+k2)*h));
                %                 sigma12 = sqrt(g*abs(k1)+abs(k2)*tanh(abs(k1)+abs(k2)*h));
                
                
                chs12 = cosh(abs(k1+k2)*h*(1+ss))/cosh(abs(k1+k2)*h);
                shs12 = sinh(abs(k1+k2)*h*(1+ss))/cosh(abs(k1+k2)*h);
                
                h12(s) = (1i*w12*g*g)/(sigma12^2-w12^2)*...
                    ((k1*k2)/(w1*w2)-(w1*w2 + w1*w1 + w2*w2)/(2*(g*g))+...
                    (k1*k1*w2+k2*k2*w1)/(2*w1*w2*w12));
                
                %                 h12(s) = (1i*w12*g*g)/(sigma12^2-w12^2);
                
                u12(s) = 1i*(k1+k2)*h12(s)*chs12+ (g*(h+z))/(2*(h+eta)) * ...
                    ((k1*k1*shs1)/w1 + (k2*k2*shs2)/w2);
                
                ww12(s) = abs(k1+k2)*h12(s)*shs12 - 1i*(g*(h+z))/(2*(h+eta))* ...
                    ((k1*abs(k1)*chs1)/w1 + (k2*abs(k2)*chs2)/w2);
                
                e1(s) = exp(1i*(k1*x-w1*t));
                e2(s) = exp(1i*(k2*x-w2*t));
                
                a1(s) = a(i)/2;
                a2(s) = a(j)/2;
                
                u2_temp(s) = a1(s)*a2(s)*u12(s)*e1(s)*e2(s);
                w2_temp(s) = a1(s)*a2(s)*ww12(s)*e1(s)*e2(s);
                
                
            end
            
        end
        %         h12
        %         return
        %             return
        u2_sum(i,j) = real(sum(u2_temp(1:2)));
        u2_diff(i,j) = real(sum(u2_temp(3:4)));
        w2_sum(i,j) = real(sum(w2_temp(1:2)));
        w2_diff(i,j) = real(sum(w2_temp(3:4)));
        %         return
    end
end
figure('Position',[100 100 900 700])
subplot(2,2,1)
surf(k,k,u2_sum)
view(0,90),xlabel('k1'),ylabel('k2')
title('u2_{sum}')
colorbar;
subplot(2,2,2)
surf(k,k,u2_diff)
view(0,90),xlabel('k1'),ylabel('k2')
title('u2_{diff}')
colorbar;
subplot(2,2,3)
surf(k,k,w2_sum)
view(0,90),xlabel('k1'),ylabel('k2')
title('w2_{sum}')
colorbar;
subplot(2,2,4)
surf(k,k,w2_diff)
view(0,90),xlabel('k1'),ylabel('k2')
title('w2_{diff}')
colorbar;
u2_sum = sum(sum(u2_sum));
u2_diff = sum(sum(u2_diff));
u2 = sum(sum(u2_sum+u2_diff));

w2_sum = sum(sum(w2_sum));
w2_diff = sum(sum(w2_diff));
w2 = sum(sum(w2_sum+w2_diff));


% return

% /* Second order wave elevation */
% traditional implementation


eta2b_sum = zeros(length(k));
eta2b_diff = zeros(length(k));
u2b_sum = zeros(length(k));
u2b_diff = zeros(length(k));
w2b_sum = zeros(length(k));
w2b_diff = zeros(length(k));
for i = 1:length(k)
    phi_i = k(i)*x - w(i) * t;
    Rn = k(i) * tanh(k(i) * h);
    
    eta2b_sum(i,i) = 0.5*a(i)^2*k(i)*cos(2*phi_i);
    
    for m = (i+1):length(k),
        gamma_nm = 1;
        k_nm_plus = sqrt(k(i) * k(i) + k(m) * k(m) + (2. * k(i) * k(m) * gamma_nm));
        k_nm_minus = sqrt(k(i) * k(i) + k(m) * k(m) - (2. * k(i) * k(m) * gamma_nm));
        
        Rm = k(m) * tanh(k(m) * h);
        
        
        D_nm_plus = (sqrt(Rn) + sqrt(Rm)) * (sqrt(Rm) * (k(i) * k(i) - Rn * Rn) + sqrt(Rn) * (k(m) * k(m) - Rm * Rm)) +...
            2. * (sqrt(Rn) + sqrt(Rm))^2. * (k(i) * k(m) * gamma_nm - Rn * Rm) / ((sqrt(Rn) + sqrt(Rm))^2. - k_nm_plus * tanh(k_nm_plus * h));
        D_nm_minus = (sqrt(Rn) - sqrt(Rm)) * (sqrt(Rm) * (k(i) * k(i) - Rn * Rn) - sqrt(Rn) * (k(m) * k(m) - Rm * Rm)) +...
            2. * (sqrt(Rn) - sqrt(Rm))^2. * (k(i) * k(m) * gamma_nm + Rn * Rm) / ((sqrt(Rn) - sqrt(Rm))^2. - k_nm_minus * tanh(k_nm_minus * h));
        
        %          Catch NaN when two equal frequency components interact
        if (w(i) == w(m))
            D_nm_minus = 0.;
        end
        
        alpha_nm_minus = (((w(i) / w(m)) + (w(m) / w(i))) + (g * g / (w(i) * w(m))) *...
            ((D_nm_minus - k(i) * k(m) * (gamma_nm + tanh(k(i) * h) * tanh(k(m) * h))) / (w(i) * w(m))));
        alpha_nm_plus = (((w(i) / w(m)) + (w(m) / w(i))) + (g * g / (w(i) * w(m))) *...
            ((D_nm_plus - k(i) * k(m) * (gamma_nm - tanh(k(i) * h) * tanh(k(m) * h))) / (w(i) * w(m))));
        
       
        beta_nm_minus = D_nm_minus / (2 * k(i) * k(m) * (w(i) - w(m)));
		beta_nm_plus = D_nm_plus / (2 * k(i) * k(m) * (w(i) + w(m)));
%         beta_nm_minus = 2*(w(i)-w(m))/((w(i)-w(m))^2-g*k_nm_minus);
%         beta_nm_plus = 0*(w(i)+w(m))/((w(i)+w(m))^2-g*k_nm_plus);
%         return
        phi_m = k(m) * x - w(m) * t;
        
        eta2b_sum(i,m) = ((a(i) * a(m) * w(i) * w(m)) / (2. * g)) * (alpha_nm_plus * cos(phi_i + phi_m));
        eta2b_diff(i,m) = ((a(i) * a(m) * w(i) * w(m)) / (2. * g)) * (alpha_nm_minus * cos(phi_i - phi_m));
        
        u2b_sum(i,m) = (a(i) * a(m) * w(i) * w(m))*(k(i)+k(m))*beta_nm_plus*...
            cos(phi_i+phi_m)*cosh(k_nm_plus*(z+h))/cosh(k_nm_plus*h);
        
        u2b_diff(i,m) = (a(i) * a(m) * w(i) * w(m))*(k(i)-k(m))*beta_nm_minus*...
            cos(phi_i-phi_m)*cosh(k_nm_minus*(z+h))/cosh(k_nm_minus*h);
        
        w2b_sum(i,m) = (a(i) * a(m) * w(i) * w(m))*(k_nm_plus)*beta_nm_plus*...
            sin(phi_i+phi_m)*sinh(k_nm_plus*(z+h))/cosh(k_nm_plus*h);
        
        w2b_diff(i,m) = (a(i) * a(m) * w(i) * w(m))*(k_nm_minus)*beta_nm_minus*...
            sin(phi_i-phi_m)*sinh(k_nm_minus*(z+h))/cosh(k_nm_minus*h);
    end
end

% figure('Position',[100 100 900 400])
% subplot(1,2,1)
% surf(k,k,eta2b_sum)
% view(0,90),xlabel('k1'),ylabel('k2')
% colorbar;
% subplot(1,2,2)
% surf(k,k,eta2b_diff)
% view(0,90),xlabel('k1'),ylabel('k2')
% colorbar;
eta2b_sum = sum(sum(eta2b_sum));
eta2b_diff = sum(sum(eta2b_diff));
eta2b = eta2b_sum+eta2b_diff;

% figure('Position',[100 100 900 700])
% subplot(2,2,1)
% surf(k,k,u2b_sum)
% view(0,90),xlabel('k1'),ylabel('k2')
% title('u2b_{sum}')
% colorbar;
% subplot(2,2,2)
% surf(k,k,u2b_diff)
% view(0,90),xlabel('k1'),ylabel('k2')
% title('u2b_{diff}')
% colorbar;
% subplot(2,2,3)
% surf(k,k,w2b_sum)
% view(0,90),xlabel('k1'),ylabel('k2')
% title('w2b_{sum}')
% colorbar;
% subplot(2,2,4)
% surf(k,k,w2b_diff)
% view(0,90),xlabel('k1'),ylabel('k2')
% title('w2b_{diff}')
% colorbar;
u2b_sum = sum(sum(u2b_sum));
u2b_diff = sum(sum(u2b_diff));
u2b = u2b_sum+u2b_diff;

w2b_sum = sum(sum(w2b_sum));
w2b_diff = sum(sum(w2b_diff));
w2b = w2b_sum+w2b_diff;



% Sharma & Dean approach 1981
% /* Second order wave elevation */
% traditional implementation



eta2c_sum = zeros(length(k));
eta2c_diff = zeros(length(k));
for i = 1:length(k)
    phi_i = k(i)*x - w(i) * t;
    Rn = k(i) * tanh(k(i) * h);
    
    for m = 1:length(k),
        gamma_nm = 1;
        %         k_nm_plus = sqrt(k(i) * k(i) + k(m) * k(m) + (2. * k(i) * k(m) * gamma_nm));
        %         k_nm_minus = sqrt(k(i) * k(i) + k(m) * k(m) - (2. * k(i) * k(m) * gamma_nm));
        
        k_nm_plus = abs(k(i)+k(m));
        k_nm_minus = abs(k(i)-k(m));
        
        Rm = k(m) * tanh(k(m) * h);
        
        D_nm_plus = (2. * (sqrt(Rn) + sqrt(Rm))^2. * (k(i) * k(m) - Rn * Rm))/...
            ((sqrt(Rn) + sqrt(Rm))^2. - k_nm_plus * tanh(k_nm_plus * h))+...
            ((sqrt(Rn) + sqrt(Rm)) * (sqrt(Rn)*(k(m)^2 - Rm^2)+sqrt(Rm)*(k(i)^2-Rn^2)))/...
            ((sqrt(Rn) + sqrt(Rm))^2. - k_nm_plus * tanh(k_nm_plus * h));
        
        
        D_nm_minus = (2. * (sqrt(Rn) - sqrt(Rm))^2. * (k(i) * k(m) + Rn * Rm))/...
            ((sqrt(Rn) - sqrt(Rm))^2. - k_nm_minus * tanh(k_nm_minus * h))+...
            ((sqrt(Rn) - sqrt(Rm)) * (-sqrt(Rn)*(k(m)^2 - Rm^2)+sqrt(Rm)*(k(i)^2-Rn^2)))/...
            ((sqrt(Rn) - sqrt(Rm))^2. - k_nm_minus * tanh(k_nm_minus * h));
        
        %          Catch NaN when two equal frequency components interact
        if (w(i) == w(m))
            D_nm_minus = 0.;
        end
        
        %         alpha_nm_minus = (((w(i) / w(m)) + (w(m) / w(i))) + (g * g / (w(i) * w(m))) *...
        %         ((D_nm_minus - k(i) * k(m) * (gamma_nm + tanh(k(i) * h) * tanh(k(m) * h))) / (w(i) * w(m))));
        
        alpha_nm_minus = (D_nm_minus-(k(i)*k(m)+Rn*Rm))/sqrt(Rn*Rm) + (Rn + Rm);
        
        %     alpha_nm_plus = (((w(i) / w(m)) + (w(m) / w(i))) + (g * g / (w(i) * w(m))) *...
        %         ((D_nm_plus - k(i) * k(m) * (gamma_nm - tanh(k(i) * h) * tanh(k(m) * h))) / (w(i) * w(m))));
        
        alpha_nm_plus = (D_nm_plus-(k(i)*k(m)-Rn*Rm))/sqrt(Rn*Rm) + (Rn + Rm);
        
        phi_m = k(m) * x - w(m) * t;
        
        eta2c_sum(i,m) = 0.25* a(i) * a(m)  * (alpha_nm_plus * cos(phi_i + phi_m));
        eta2c_diff(i,m) = 0.25* a(i) * a(m) * (alpha_nm_minus * cos(phi_i - phi_m));
        
        
        u2c_sum(i,m) = 0.25*(a(i)*g/w(i))*(a(m)*g/w(m))*(k(i)+k(m))*((D_nm_plus/(w(i)+w(m)))*cos(phi_i+phi_m)*...
            cosh(k_nm_plus*(z+h))/cosh(k_nm_plus*h));
        if (w(i) == w(m))
            u2c_diff(i,m) = 0.;
        else
            u2c_diff(i,m) = 0.25*(a(i)*g/w(i))*(a(m)*g/w(m))*(k(i)-k(m))*((D_nm_minus/(w(i)-w(m)))*cos(phi_i-phi_m)*...
                cosh(k_nm_minus*(z+h))/cosh(k_nm_minus*h));
        end
        
        w2c_sum(i,m) = 0.25*(a(i)*g/w(i))*(a(m)*g/w(m))*(k_nm_plus)*((D_nm_plus/(w(i)+w(m)))*sin(phi_i+phi_m)*...
            sinh(k_nm_plus*(z+h))/cosh(k_nm_plus*h));
        if (w(i) == w(m))
            w2c_diff(i,m) = 0.;
        else
            w2c_diff(i,m) = 0.25*(a(i)*g/w(i))*(a(m)*g/w(m))*(k_nm_minus)*((D_nm_minus/(w(i)-w(m)))*sin(phi_i-phi_m)*...
                sinh(k_nm_minus*(z+h))/cosh(k_nm_minus*h));
        end
    end
end

% figure('Position',[100 100 900 400])
% subplot(1,2,1)
% surf(k,k,eta2c_sum)
% view(0,90),xlabel('k1'),ylabel('k2')
% colorbar;
% subplot(1,2,2)
% surf(k,k,eta2c_diff)
% view(0,90),xlabel('k1'),ylabel('k2')
% colorbar;
eta2c_sum = sum(sum(eta2c_sum));
eta2c_diff = sum(sum(eta2c_diff));
eta2c = eta2c_sum+eta2c_diff;

figure('Position',[100 100 900 700])
subplot(2,2,1)
surf(k,k,u2c_sum)
view(0,90),xlabel('k1'),ylabel('k2')
title('u2c_{sum}')
colorbar;
subplot(2,2,2)
surf(k,k,u2c_diff)
view(0,90),xlabel('k1'),ylabel('k2')
title('u2c_{diff}')
colorbar;
subplot(2,2,3)
surf(k,k,w2c_sum)
view(0,90),xlabel('k1'),ylabel('k2')
title('w2c_{sum}')
colorbar;
subplot(2,2,4)
surf(k,k,w2c_diff)
view(0,90),xlabel('k1'),ylabel('k2')
title('w2c_{diff}')
colorbar;
u2c_sum = sum(sum(u2c_sum))
u2c_diff = sum(sum(u2c_diff))
u2c = u2c_sum+u2c_diff
w2c_sum = sum(sum(w2c_sum))
w2c_diff = sum(sum(w2c_diff))
w2c = w2c_sum+w2c_diff

% close all;
% [eta2_sum, eta2b_sum, eta2c_sum]
% [eta2_diff, eta2b_diff, eta2c_diff]

[u2_sum u2b_sum u2c_sum]
[u2_diff u2b_diff u2c_diff]

[w2_sum w2b_sum w2c_sum]
[w2_diff w2b_diff w2c_diff]
