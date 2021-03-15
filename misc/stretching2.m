clear;

% sfac = 1.0;
%
% xc = 0.5;
%
% nx = 10;
%
% xcoord = sfac.^[1:nx]
nl = 15;

s = -1:1/(nl-1):0;

a = 7*pi/18;
b = 1.0;


for i=1:nl,   
    x(i)= (-(tan(((-s(i))*a))^b)/(tan(a)^b));
end




plot(s,x,'LineWidth',2);hold on
for i=1:nl,
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