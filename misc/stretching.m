clear;

% sfac = 1.0;
%
% xc = 0.5;
%
% nx = 10;
%
% xcoord = sfac.^[1:nx]


% % A program to test 1D grid refinement.
L=5.0;
xc=1.0;
A=2
s=[0:0.1:1];n=size(s);
% % Alternative 1
% for i=1:n(2),
%     x(i)=L*s(i)+A*(xc-L*s(i))*s(i)*(1-s(i));
% end

% % Alternative 2
% for i=1:n(2),
%  x(i)= cosd(s(i)*90);
% end

% Alternative 2
aa = 0;
for i=1:n(2),
    sfac = 3.0;
    x(i)= 1+(n(2)-(i-1))*sfac;
    aa = aa + x(i);
end
x = x./aa
x = cumsum(x)
% return
% double sfac = 3.0;
% double dd[nl];
% double ddsum = 0.;
% for (int ii=0; ii<nl; ii++){
%     //fprintf(stdout,"%d",ii);
%     dd[ii] = (1.+(nl-ii)*sfac);
%     ddsum+=dd[ii];
%     }
%     
%     return dd[layerno]/ddsum;

% Alternative 2 NICE
a_end = 0.01;
a_potens = 1;
for i=1:n(2),   
    x(i)= 1-(tand((1-s(i))*a_end)^a_potens)/(tand(a_end)^a_potens);
    % and the inverse
    s(i) = atan(((1- x(i))*(tan(a_end*pi/180)^a_potens))^(1/a_potens))/(a_end*pi/180);
end




plot(s,x,'LineWidth',2);hold on
for i=1:n(2),
    plot([s(i),s(i)],[0.0, x(i)]);
    plot([0.0,s(i)],[x(i), x(i)]);
end
xlabel('\xi','Fontsize',18)
ylabel('x','Fontsize',18)
set(gca,'Box','on');
set(gca,'Fontsize',18, 'LineWidth',2)
% plot([0,1],[0,L],'r')
text(0.05,L-0.1,['A=',num2str(A),...
    ' x=',num2str(xc)],'Fontsize',18)
hold off; % print -depsc exampleplot