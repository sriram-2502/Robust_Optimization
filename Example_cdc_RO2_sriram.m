clc;
clear;
close all;

%%
tspan=[0 500];
z0=0*randn(8,1);
l0=1e-6*ones(4,1);
r0=1e-6*ones(8,1);
y0=[z0;l0;r0];
global t0
t0=10;

% change it to false to have moving targets (sriram)
stationary = true;

[t,y]=ode15s(@(t,y) odeloc(t,y,stationary),tspan,y0);
om=pi/30;
om=pi/50;


if (stationary)
    an = 0; % stationary
else
    an=[0*t(t<=t0);om*(t(t>=t0)-t0)]; % moving targets
end

b1=[1;1]*5;
b2=[-1;1]*5;
b3=[-1;-1]*5;
b5=[1;-1]*5;
b4=[0; -1*sqrt(2)]*5;

% moving targets
C1=[[cos(an),sin(an)]*b1, [-sin(an),cos(an)]*b1];
C2=[[cos(an),sin(an)]*b2, [-sin(an),cos(an)]*b2];
C3=[[cos(an),sin(an)]*b3, [-sin(an),cos(an)]*b3];
C4=[[cos(an),sin(an)]*b4, [-sin(an),cos(an)]*b4];
C5=[[cos(an),sin(an)]*b5,[-sin(an),cos(an)]*b5];

n=size(y,1);
y(n,:)

%% original plots
figure(1);
plot(t,y(:,1:8));
xlabel('Time sec');
ylabel('Agents time response')

figure(3)
plot(y(1:end,1),y(1:end,2),y(1:end,3),y(1:end,4),y(1:end,5),y(1:end,6),y(1:end,7),y(1:end,8))

grid
figure(2);
hold on
plot(y(20,1),y(20,2),'o','MarkerSize',10);
plot(y(985,1),y(985,2),'o','MarkerSize',10,'MarkerFaceColor','b');
plot(y(985,3),y(985,4),'o','MarkerSize',10,'MarkerFaceColor','r');
plot(y(985,5),y(985,6),'o','MarkerSize',10,'MarkerFaceColor','y');
plot(y(985,7),y(985,8),'o','MarkerSize',10,'MarkerFaceColor','g');
plot(y(end,1),y(end,2),'ko','MarkerSize',10);
plot(y(end,3),y(end,4),'ko','MarkerSize',10);
plot(y(end,5),y(end,6),'ko','MarkerSize',10);
plot(y(end,7),y(end,8),'ko','MarkerSize',10);
hold off
grid


%% gif
figure(5);
clear Frame
t_len = 100;
for j=1:t_len, k=max(find(t<=j));
if(stationary)
% plot stationary target
    plot([y(k,1);C1(1,1)],[y(k,2);C1(1,2)],[y(k,1);C2(1,1)],[y(k,2);C2(1,2)],...
        [y(k,3);C3(1,1)],[y(k,4); C3(1,2)],[y(k,5);C4(1,1)],[y(k,6); C4(1,2)],...
        [y(k,7);C1(1,1)],[y(k,8); C1(1,2)],[y(k,7);C5(1,1)],[y(k,8); C5(1,2)],...
        [y(k,1);y(k,3)],[y(k,2);y(k,4)],[y(k,3);y(k,5)],[y(k,4);y(k,6)], ...
        [y(k,7);y(k,1)],[y(k,8);y(k,2)],[y(k,3);y(k,7)],[y(k,4);y(k,8)], ...
        [y(k,7);y(k,5)],[y(k,8);y(k,6)],...
        y(k,1),y(k,2),'o',y(k,3),y(k,4),'o',y(k,5),y(k,6),'o',y(k,7),y(k,8),'o',...
        C1(1,1),C1(1,2),'sk',C2(1,1),C2(1,2),'sk',C3(1,1),C3(1,2),'sk',C4(1,1),C4(1,2),'sk',C5(1,1),C5(1,2),'sk');
else
    % plot moving targets gif
    plot([y(k,1);C1(k,1)],[y(k,2);C1(k,2)],[y(k,1);C2(k,1)],[y(k,2);C2(k,2)],[y(k,3);C3(k,1)],[y(k,4); C3(k,2)],[y(k,5);C4(k,1)],[y(k,6); C4(k,2)],[y(k,7);C1(k,1)],[y(k,8); C1(k,2)],[y(k,7);C5(k,1)],[y(k,8); C5(k,2)],[y(k,1);y(k,3)],[y(k,2);y(k,4)],[y(k,3);y(k,5)],[y(k,4);y(k,6)],[y(k,7);y(k,1)],[y(k,8);y(k,2)],[y(k,3);y(k,7)],[y(k,4);y(k,8)],[y(k,7);y(k,5)],[y(k,8);y(k,6)],y(k,1),y(k,2),'o',y(k,3),y(k,4),'o',y(k,5),y(k,6),'o',y(k,7),y(k,8),'o',C1(k,1),C1(k,2),'sk',C2(k,1),C2(k,2),'sk',C3(k,1),C3(k,2),'sk',C4(k,1),C4(k,2),'sk',C5(k,1),C5(k,2),'sk');
end

axis([-8,8,-8,8]);grid on

pause(.1);
text(6,-6,string(j))
Frame(k) = getframe;
end
grid off

%% plot initial configuration (sriram)
figure(6)
colors = colororder;
blue = colors(1,:);
red = colors(2,:);
yellow = colors(3,:);
green = colors(5,:);
gray = [.7 .7 .7]; % Obstacle color -> Grey

s1 = subplot(2,2,1);
box(s1,'on')
grid(s1,'on');
hold on
start = 1;
last = 985;

% plot constraint line
% plot shaded area for uncertainty
% fill vertices of a triangle in each side of the constraint line
u1 = fill([1.25,8,8],[1.25,-6.3,-5.5],gray,'FaceAlpha',0.5,'EdgeColor','none');
fill([1.25,-8,-8],[1.25,10.5,9.7],gray,'FaceAlpha',0.5,'EdgeColor','none')
fill([1.25,-8,-8],[1.25,10.5,11.3],gray,'FaceAlpha',0.5,'EdgeColor','none')
fill([1.25,8,8],[1.25,-4.7,-5.5],gray,'FaceAlpha',0.5,'EdgeColor','none')

xx = -8:0.01:8;
rho = -1; %uncertainty
yy = rho*xx + 5/2;
c = plot(xx,yy,'--k','LineWidth',2);
rho = -1+0.1; %uncertainty
yy = rho*xx + 5/2;
%plot(xx,yy)
rho = -1-0.1; %uncertainty
yy = rho*xx + 5/2;
%plot(xx,yy)

% plot line for traj
plot(y(start:last,1),y(start:last,2),'color',blue,'LineWidth',2);
plot(y(start:last,3),y(start:last,4),'color',red,'LineWidth',2);
plot(y(start:last,5),y(start:last,6),'color',yellow,'LineWidth',2);
plot(y(start:last,7),y(start:last,8),'color',green,'LineWidth',2);

% % plot trail for agents
% idx = 0;
% tspan = start:10:last;
% for i = tspan
%     alpha1 = 0.99^(length(tspan)-idx);
%     scatter(y(i,1),y(i,2),'o','MarkerEdgecolor','none','MarkerFaceColor', blue,'MarkerFaceAlpha',alpha1);
%     scatter(y(i,3),y(i,4),'o','MarkerEdgecolor','none','MarkerFaceColor', red,'MarkerFaceAlpha',alpha1);
%     scatter(y(i,5),y(i,6),'o','MarkerEdgecolor','none','MarkerFaceColor', yellow,'MarkerFaceAlpha',alpha1);
%     scatter(y(i,7),y(i,8),'o','MarkerEdgecolor','none','MarkerFaceColor', green,'MarkerFaceAlpha',alpha1);
%     idx = idx+1;
% end

% plot line between agents
k = last;  
plot([y(k,1);y(k,3)],[y(k,2);y(k,4)],'color','k'); hold on;
plot([y(k,3);y(k,5)],[y(k,4);y(k,6)],'color','k'); hold on;
plot([y(k,7);y(k,1)],[y(k,8);y(k,2)],'color','k'); hold on;
plot([y(k,3);y(k,7)],[y(k,4);y(k,8)],'color','k'); hold on;
plot([y(k,7);y(k,5)],[y(k,8);y(k,6)],'color','k'); hold on;

% plot line between agents and targets
if(stationary)
    plot([y(k,1);C1(1,1)],[y(k,2);C1(1,2)],'-','color',red); hold on;
    plot([y(k,7);C1(1,1)],[y(k,8); C1(1,2)],'-','color',green); hold on;
    plot([y(k,1);C2(1,1)],[y(k,2);C2(1,2)],'-','color',red); hold on;
    plot([y(k,3);C3(1,1)],[y(k,4);C3(1,2)],'-','color',blue); hold on;
    plot([y(k,5);C4(1,1)],[y(k,6); C4(1,2)],'-','color',yellow); hold on;
    plot([y(k,7);C5(1,1)],[y(k,8); C5(1,2)],'-','color',green); hold on;

    % plot targets
    alpha2 = 1; sz = 100;
    t1=scatter(C1(1,1),C1(1,2),sz,'s','MarkerEdgecolor',red,'MarkerFaceColor',green, 'MarkerFaceAlpha', alpha2, 'LineWidth',2); hold on
    t2=scatter(C2(1,1),C2(1,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',red, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    t3=scatter(C3(1,1),C3(1,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',blue, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    t4=scatter(C4(1,1),C4(1,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',yellow, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    t5=scatter(C5(1,1),C5(1,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',green, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on


else
    plot([y(k,1);C1(k,1)],[y(k,2);C1(k,2)],'-','color',red); hold on;
    plot([y(k,7);C1(k,1)],[y(k,8); C1(k,2)],'-','color',green); hold on;
    plot([y(k,1);C2(k,1)],[y(k,2);C2(k,2)],'-','color',red); hold on;
    plot([y(k,3);C3(k,1)],[y(k,4);C3(k,2)],'-','color',blue); hold on;
    plot([y(k,5);C4(k,1)],[y(k,6); C4(k,2)],'-','color',yellow); hold on;
    plot([y(k,7);C5(k,1)],[y(k,8); C5(k,2)],'-','color',green); hold on;

    %plot targets
    alpha2 = 1; sz = 100;
    scatter(C1(k,1),C1(k,2),sz,'s','MarkerEdgecolor',red,'MarkerFaceColor',green, 'MarkerFaceAlpha', alpha2, 'LineWidth',2); hold on
    scatter(C2(k,1),C2(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',red, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C3(k,1),C3(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',blue, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C4(k,1),C4(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',yellow, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C5(k,1),C5(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',green, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
end

% plot start position
alpha = 0.5; sz = 50;
scatter(y(start,1),y(start,2),sz,'MarkerFaceColor',blue,'MarkerEdgeColor','k','MarkerFaceAlpha', alpha);
scatter(y(start,3),y(start,4),sz,'MarkerFaceColor',red,'MarkerEdgeColor','k','MarkerFaceAlpha', alpha);
scatter(y(start,5),y(start,6),sz,'MarkerFaceColor',yellow,'MarkerEdgeColor','k','MarkerFaceAlpha', alpha);
scatter(y(start,7),y(start,8),sz,'MarkerFaceColor',green,'MarkerEdgeColor','k','MarkerFaceAlpha', alpha);


% plot final position
a0=scatter(0,0,sz,'MarkerFaceColor','w','MarkerEdgeColor','k');
a1=scatter(y(last,1),y(last,2),sz,'MarkerFaceColor',blue,'MarkerEdgeColor','k');
a2=scatter(y(last,3),y(last,4),sz,'MarkerFaceColor',red,'MarkerEdgeColor','k');
a3=scatter(y(last,5),y(last,6),sz,'MarkerFaceColor',yellow,'MarkerEdgeColor','k');
a4=scatter(y(last,7),y(last,8),sz,'MarkerFaceColor',green,'MarkerEdgeColor','k');

% legend([a0,a1,a2,a3,a4,t1,t2,t3,t4,t5,c,u1], ...
%     'start','agent 1','agent 2','agent 3','agent 4',...
%     'target 1','target 2','target 3','target 4','target 5',...
%     'constraint','uncertainty')

legend([a0,c,u1], 'start','constraint','uncertainty')

axis square
xlim([-8,8])
ylim([-8,8])
axes = gca;
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;

%% Change in uncertainty value (sriram)
figure(6)
s2 = subplot(2,2,2);
box(s2,'on')
grid(s2,'on');
hold on

% plot constraint line
% plot shaded area for uncertainty
% fill vertices of a triangle in each side of the constraint line
fill([1.25,-8,-8],[1.25,10.5,1.25],gray,'FaceAlpha',0.5,'EdgeColor','none')
fill([1.25,1.25,10.5],[1.25,-8,-8],gray,'FaceAlpha',0.5,'EdgeColor','none')
fill([1.25,-8,1.25],[1.25,10.5,10],gray,'FaceAlpha',0.5,'EdgeColor','none')
fill([1.25,10.5,10.5],[1.25,-8,1.25],gray,'FaceAlpha',0.5,'EdgeColor','none')

xx1 = -12:0.01:12;
rho = -1; %uncertainty
yy1 = rho*xx1 + 5/2; % original constraint
plot(xx1,yy1,'--k','LineWidth',2)

% plot changes in uncertainy as a shaded area
rho = 0; %uncertainty
xx2 = -12:0.01:12;
yy2 = 1.25*ones(1,length(xx2)); % vertical line passing through 1.25
%plot(xx2,yy2)

rho = 0; %uncertainty
yy = -12:0.01:12;
xx = 1.25*ones(length(yy)); % horizontal line passing through 1.25
%plot(xx,yy)

start = 985;
last = length(y);


% plot line for traj
plot(y(start:last,1),y(start:last,2),'color',blue,'LineWidth',2);
plot(y(start:last,3),y(start:last,4),'color',red,'LineWidth',2);
plot(y(start:last,5),y(start:last,6),'color',yellow,'LineWidth',2);
plot(y(start:last,7),y(start:last,8),'color',green,'LineWidth',2);

% plot trail for agents
% idx = 0;
% tspan = start:10:last;
% for i = tspan
%     alpha1 = 0.99^(length(tspan)-idx);
%     scatter(y(i,1),y(i,2),'o','MarkerEdgecolor','none','MarkerFaceColor', blue,'MarkerFaceAlpha',alpha1);
%     scatter(y(i,3),y(i,4),'o','MarkerEdgecolor','none','MarkerFaceColor', red,'MarkerFaceAlpha',alpha1);
%     scatter(y(i,5),y(i,6),'o','MarkerEdgecolor','none','MarkerFaceColor', yellow,'MarkerFaceAlpha',alpha1);
%     scatter(y(i,7),y(i,8),'o','MarkerEdgecolor','none','MarkerFaceColor', green,'MarkerFaceAlpha',alpha1);
%     idx = idx+1;
% end

% plot line between agents
k = last;  
plot([y(k,1);y(k,3)],[y(k,2);y(k,4)],'color','k'); hold on;
plot([y(k,3);y(k,5)],[y(k,4);y(k,6)],'color','k'); hold on;
plot([y(k,7);y(k,1)],[y(k,8);y(k,2)],'color','k'); hold on;
plot([y(k,3);y(k,7)],[y(k,4);y(k,8)],'color','k'); hold on;
plot([y(k,7);y(k,5)],[y(k,8);y(k,6)],'color','k'); hold on;

% plot line between agents and targets
if(stationary)
    plot([y(k,1);C1(1,1)],[y(k,2);C1(1,2)],'-','color',red); hold on;
    plot([y(k,7);C1(1,1)],[y(k,8); C1(1,2)],'-','color',green); hold on;
    plot([y(k,1);C2(1,1)],[y(k,2);C2(1,2)],'-','color',red); hold on;
    plot([y(k,3);C3(1,1)],[y(k,4);C3(1,2)],'-','color',blue); hold on;
    plot([y(k,5);C4(1,1)],[y(k,6); C4(1,2)],'-','color',yellow); hold on;
    plot([y(k,7);C5(1,1)],[y(k,8); C5(1,2)],'-','color',green); hold on;

    % plot targets
    alpha2 = 1; sz = 100;
    scatter(C1(1,1),C1(1,2),sz,'s','MarkerEdgecolor',red,'MarkerFaceColor',green, 'MarkerFaceAlpha', alpha2, 'LineWidth',2); hold on
    scatter(C2(1,1),C2(1,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',red, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C3(1,1),C3(1,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',blue, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C4(1,1),C4(1,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',yellow, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C5(1,1),C5(1,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',green, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on


else
    plot([y(k,1);C1(k,1)],[y(k,2);C1(k,2)],'-','color',red); hold on;
    plot([y(k,7);C1(k,1)],[y(k,8); C1(k,2)],'-','color',green); hold on;
    plot([y(k,1);C2(k,1)],[y(k,2);C2(k,2)],'-','color',red); hold on;
    plot([y(k,3);C3(k,1)],[y(k,4);C3(k,2)],'-','color',blue); hold on;
    plot([y(k,5);C4(k,1)],[y(k,6); C4(k,2)],'-','color',yellow); hold on;
    plot([y(k,7);C5(k,1)],[y(k,8); C5(k,2)],'-','color',green); hold on;

    % plot targets
    alpha2 = 1; sz = 100;
    scatter(C1(k,1),C1(k,2),sz,'s','MarkerEdgecolor',red,'MarkerFaceColor',green, 'MarkerFaceAlpha', alpha2, 'LineWidth',2); hold on
    scatter(C2(k,1),C2(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',red, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C3(k,1),C3(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',blue, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C4(k,1),C4(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',yellow, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C5(k,1),C5(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',green, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
end

% plot start position
alpha = 0.5; sz = 50;
scatter(y(start,1),y(start,2),sz,'MarkerFaceColor',blue,'MarkerEdgeColor','k','MarkerFaceAlpha', alpha);
scatter(y(start,3),y(start,4),sz,'MarkerFaceColor',red,'MarkerEdgeColor','k','MarkerFaceAlpha', alpha);
scatter(y(start,5),y(start,6),sz,'MarkerFaceColor',yellow,'MarkerEdgeColor','k','MarkerFaceAlpha', alpha);
scatter(y(start,7),y(start,8),sz,'MarkerFaceColor',green,'MarkerEdgeColor','k','MarkerFaceAlpha', alpha);


% plot final position
scatter(y(last,1),y(last,2),sz,'MarkerFaceColor',blue,'MarkerEdgeColor','k');
scatter(y(last,3),y(last,4),sz,'MarkerFaceColor',red,'MarkerEdgeColor','k');
scatter(y(last,5),y(last,6),sz,'MarkerFaceColor',yellow,'MarkerEdgeColor','k');
scatter(y(last,7),y(last,8),sz,'MarkerFaceColor',green,'MarkerEdgeColor','k');


axis square
xlim([-8,8])
ylim([-8,8])

axes = gca;
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;


%% zoom plot for one agent under change in uncertainty (sriram)
figure(6)
s3 = subplot(2,2,3);
box(s3,'on')
grid(s3,'on');
hold on

% plot constraint line
% plot shaded area for uncertainty
% fill vertices of a triangle in each side of the constraint line
fill([1.25,-8,-8],[1.25,10.5,1.25],gray,'FaceAlpha',0.5,'EdgeColor','none')
fill([1.25,1.25,10.5],[1.25,-8,-8],gray,'FaceAlpha',0.5,'EdgeColor','none')
fill([1.25,-8,1.25],[1.25,10.5,10],gray,'FaceAlpha',0.5,'EdgeColor','none')
fill([1.25,10.5,10.5],[1.25,-8,1.25],gray,'FaceAlpha',0.5,'EdgeColor','none')

xx1 = -12:0.01:12;
rho = -1; %uncertainty
yy1 = rho*xx1 + 5/2; % original constraint
plot(xx1,yy1,'--k','LineWidth',2)

% plot changes in uncertainy as a shaded area
rho = 0; %uncertainty
xx2 = -12:0.01:12;
yy2 = 1.25*ones(1,length(xx2)); % vertical line passing through 1.25
%plot(xx2,yy2)

rho = 0; %uncertainty
yy = -12:0.01:12;
xx = 1.25*ones(length(yy)); % horizontal line passing through 1.25
%plot(xx,yy)


start = 985;
last = length(y);

% plot line for traj
plot(y(start:last,1),y(start:last,2),'color',blue,'LineWidth',2);
plot(y(start:last,3),y(start:last,4),'color',red,'LineWidth',2);
plot(y(start:last,5),y(start:last,6),'color',yellow,'LineWidth',2);
plot(y(start:last,7),y(start:last,8),'color',green,'LineWidth',2);

% plot trail for agents
% idx = 0;
% tspan = start:10:last;
% for i = tspan
%     alpha1 = 0.99^(length(tspan)-idx);
%     scatter(y(i,1),y(i,2),'o','MarkerEdgecolor','none','MarkerFaceColor', blue,'MarkerFaceAlpha',alpha1);
%     scatter(y(i,3),y(i,4),'o','MarkerEdgecolor','none','MarkerFaceColor', red,'MarkerFaceAlpha',alpha1);
%     scatter(y(i,5),y(i,6),'o','MarkerEdgecolor','none','MarkerFaceColor', yellow,'MarkerFaceAlpha',alpha1);
%     scatter(y(i,7),y(i,8),'o','MarkerEdgecolor','none','MarkerFaceColor', green,'MarkerFaceAlpha',alpha1);
%     idx = idx+1;
% end

% plot line between agents
k = last;  
plot([y(k,1);y(k,3)],[y(k,2);y(k,4)],'color','k'); hold on;
plot([y(k,3);y(k,5)],[y(k,4);y(k,6)],'color','k'); hold on;
plot([y(k,7);y(k,1)],[y(k,8);y(k,2)],'color','k'); hold on;
plot([y(k,3);y(k,7)],[y(k,4);y(k,8)],'color','k'); hold on;
plot([y(k,7);y(k,5)],[y(k,8);y(k,6)],'color','k'); hold on;

% plot line between agents and targets
if(stationary)
    plot([y(k,1);C1(1,1)],[y(k,2);C1(1,2)],'-','color',red); hold on;
    plot([y(k,7);C1(1,1)],[y(k,8); C1(1,2)],'-','color',green); hold on;
    plot([y(k,1);C2(1,1)],[y(k,2);C2(1,2)],'-','color',red); hold on;
    plot([y(k,3);C3(1,1)],[y(k,4);C3(1,2)],'-','color',blue); hold on;
    plot([y(k,5);C4(1,1)],[y(k,6); C4(1,2)],'-','color',yellow); hold on;
    plot([y(k,7);C5(1,1)],[y(k,8); C5(1,2)],'-','color',green); hold on;

    % plot targets
    alpha2 = 1; sz = 100;
    scatter(C1(1,1),C1(1,2),sz,'s','MarkerEdgecolor',red,'MarkerFaceColor',green, 'MarkerFaceAlpha', alpha2, 'LineWidth',2); hold on
    scatter(C2(1,1),C2(1,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',red, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C3(1,1),C3(1,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',blue, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C4(1,1),C4(1,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',yellow, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C5(1,1),C5(1,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',green, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on


else
    plot([y(k,1);C1(k,1)],[y(k,2);C1(k,2)],'-','color',red); hold on;
    plot([y(k,7);C1(k,1)],[y(k,8); C1(k,2)],'-','color',green); hold on;
    plot([y(k,1);C2(k,1)],[y(k,2);C2(k,2)],'-','color',red); hold on;
    plot([y(k,3);C3(k,1)],[y(k,4);C3(k,2)],'-','color',blue); hold on;
    plot([y(k,5);C4(k,1)],[y(k,6); C4(k,2)],'-','color',yellow); hold on;
    plot([y(k,7);C5(k,1)],[y(k,8); C5(k,2)],'-','color',green); hold on;

    % plot targets
    alpha2 = 1; sz = 100;
    scatter(C1(k,1),C1(k,2),sz,'s','MarkerEdgecolor',red,'MarkerFaceColor',green, 'MarkerFaceAlpha', alpha2, 'LineWidth',2); hold on
    scatter(C2(k,1),C2(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',red, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C3(k,1),C3(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',blue, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C4(k,1),C4(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',yellow, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
    scatter(C5(k,1),C5(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',green, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
end

% plot start position
alpha = 0.5; sz = 50;
scatter(y(start,1),y(start,2),sz,'MarkerFaceColor',blue,'MarkerEdgeColor','k','MarkerFaceAlpha', alpha);
scatter(y(start,3),y(start,4),sz,'MarkerFaceColor',red,'MarkerEdgeColor','k','MarkerFaceAlpha', alpha);
scatter(y(start,5),y(start,6),sz,'MarkerFaceColor',yellow,'MarkerEdgeColor','k','MarkerFaceAlpha', alpha);
scatter(y(start,7),y(start,8),sz,'MarkerFaceColor',green,'MarkerEdgeColor','k','MarkerFaceAlpha', alpha);

% plot final position
scatter(y(last,1),y(last,2),sz,'MarkerFaceColor',blue,'MarkerEdgeColor','k');
scatter(y(last,3),y(last,4),sz,'MarkerFaceColor',red,'MarkerEdgeColor','k');
scatter(y(last,5),y(last,6),sz,'MarkerFaceColor',yellow,'MarkerEdgeColor','k');
scatter(y(last,7),y(last,8),sz,'MarkerFaceColor',green,'MarkerEdgeColor','k');


axis square
xlim([-0.285735880983456 2.33039197224147])
ylim([-0.814016768681192 1.80211108454373])

axes = gca;
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;

%% moving targets

%%
tspan=[0 500];
z0=0*randn(8,1);
l0=1e-6*ones(4,1);
r0=1e-6*ones(8,1);
y0=[z0;l0;r0];
global t0
t0=10;
stationary = false;
[t,y]=ode15s(@(t,y) odeloc(t,y,stationary),tspan,y0);
om=pi/30;
om=pi/50;



an=[0*t(t<=t0);om*(t(t>=t0)-t0)];
b1=[1;1]*5;
b2=[-1;1]*5;
b3=[-1;-1]*5;
b5=[1;-1]*5;
b4=[0; -1*sqrt(2)]*5;

C1=[[cos(an),sin(an)]*b1, [-sin(an),cos(an)]*b1];
C2=[[cos(an),sin(an)]*b2, [-sin(an),cos(an)]*b2];
C3=[[cos(an),sin(an)]*b3, [-sin(an),cos(an)]*b3];
C4=[[cos(an),sin(an)]*b4, [-sin(an),cos(an)]*b4];
C5=[[cos(an),sin(an)]*b5,[-sin(an),cos(an)]*b5];


n=size(y,1);
y(n,:)

%% yellow trail
s4 = subplot(2,2,4);
box(s4,'on')
grid(s4,'on');
grid on
tlen = 100;
tstart = 1;
timestamps1 = tstart:5:tlen;
idx = 1;

% plot constraint line
% plot shaded area for uncertainty
% fill vertices of a triangle in each side of the constraint line
fill([1.25,8,8],[1.25,-6.3,-5.5],gray,'FaceAlpha',0.5,'EdgeColor','none')
fill([1.25,-8,-8],[1.25,10.5,9.7],gray,'FaceAlpha',0.5,'EdgeColor','none')
fill([1.25,-8,-8],[1.25,10.5,11.3],gray,'FaceAlpha',0.5,'EdgeColor','none')
fill([1.25,8,8],[1.25,-4.7,-5.5],gray,'FaceAlpha',0.5,'EdgeColor','none')

xx = -8:0.01:8;
rho = -1; %uncertainty
yy = rho*xx + 5/2;
plot(xx,yy,'--k','LineWidth',2)
rho = -1+0.1; %uncertainty
yy = rho*xx + 5/2;
%plot(xx,yy)
rho = -1-0.1; %uncertainty
yy = rho*xx + 5/2;
%plot(xx,yy)


% trail for control points (less than agent trail)
timestamps2 = tstart:1:tlen;

% plot trail for agents
k_control = [];
for j=timestamps1, k=max(find(t<=j));
    % start fading after each control point
    if any(j == timestamps2)
        alpha1 = 1;
        idx = 1;
        k_control = [k_control, k];
    else
        alpha1 = 0.99^(20/idx); % last marker with alpha 1
    end
    %scatter(y(k,5),y(k,6),'o','MarkerEdgecolor','none','MarkerFaceColor', yellow,'MarkerFaceAlpha',alpha1); hold on
    % update index for trail
    idx = idx + 1;
end

start=1;
% plot traj line for agents
for j=timestamps1, k=max(find(t<=j));
    % start fading after each control point
    plot(y(start:k,5),y(start:k,6),'Color',yellow,'LineWidth',2); hold on
end

% overlap agents at control points
sz=50;
for k = k_control
    scatter(y(k,5),y(k,6),sz,'o','MarkerEdgecolor','k','MarkerFaceColor', yellow,'MarkerFaceAlpha',alpha1); hold on
end

% plot trail for control points and lines
% idx = 1; alpha2 = 1; sz = 100;
% for j=timestamps2, k=max(find(t<=j)); 
%     plot([y(k,5);C4(k,1)],[y(k,6); C4(k,2)],'color',yellow); hold on;    
%     scatter(C4(k,1),C4(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',yellow, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
%     
%     % update index for trail
%     idx = idx + 1;
% end

% idx = 1; alpha1 = 1; sz = 100;
% timestamps2 = tlen:1:tlen;
% for j=timestamps2, k=max(find(t<=j)); 
%     % agent end pos
%     scatter(y(k,5),y(k,6),'o','MarkerEdgecolor','k','MarkerFaceColor', yellow,'MarkerFaceAlpha',alpha1); hold on 
% 
%     % line between agents and control points
%     %plot([y(k,5);C4(k,1)],[y(k,6); C4(k,2)],'color',yellow); hold on; 
%     
%     % control points
%     %scatter(C4(k,1),C4(k,2),sz,'s','MarkerEdgecolor','k','MarkerFaceColor',yellow, 'MarkerFaceAlpha',alpha2, 'LineWidth',1); hold on
% 
%     % update index for trail
%     idx = idx + 1;
% end

axis square
xlim([-4,4])
ylim([-4,4])

axes = gca;
set(axes,'FontSize',15);
xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
box on
axes.LineWidth=2;


%% functions
function dw=odeloc(t,w,stationary)
global t0
z=w(1:8);
l=w(9:12);
u1=w(13);
v1=w(14);
u2=w(15);
v2=w(16);
u3=w(17);
v3=w(18);
u4=w(19);
v4=w(20);
om=pi/50;
cc=.001;


%stationary = false;
if (stationary)
    an = 0; % stationary
else
    an=[0*t(t<=t0);om*(t(t>=t0)-t0)]; % moving targets
end

if t<=300, q0=0.1;
else q0=1;
end

b1=[1;1]*5;
b2=[-1;1]*5;
b3=[-1;-1]*5;
b5=[1;-1]*5;
b4=[0; -1*sqrt(2)]*5;

% move targets in a circle
c1=[[cos(an),sin(an)]*b1, [-sin(an),cos(an)]*b1];
c2=[[cos(an),sin(an)]*b2, [-sin(an),cos(an)]*b2];
c3=[[cos(an),sin(an)]*b3, [-sin(an),cos(an)]*b3];
c4=[[cos(an),sin(an)]*b4, [-sin(an),cos(an)]*b4];
c5=[[cos(an),sin(an)]*b5, [-sin(an),cos(an)]*b5];


c=[c1';c2';c3';c4';c5'];
I=eye(2);O=zeros(2);
Z1=[-4*I  I   O   I;
    I  -4*I  I   I;
    O    I -3*I  I;
    I    I   I  -5*I];
Z2=[I I O O O;
    O O I O O;
    O O O I O;
    I O O O I];

Z2*c;
dl1=30*(z(1)+z(2)+u1*[1 -1]*[z(1);z(2)]-2.5-v1*(u1^2-q0));
%dl1=10*(z(1)+z(2)-2.5);
dl(1,:)=((l(1)>0)|(dl1>0))*dl1;

dl2=30*(z(3)+z(4)+u2*[1 -1]*[z(3);z(4)]-2.5-v2*(u2^2-q0));
%dl2=10*(z(3)+z(4)-5);
dl(2,:)=((l(2)>0)|(dl2>0))*dl2;

dl3=30*(z(5)+z(6)+u3*[1 -1]*[z(5);z(6)]-2.5-v3*(u3^2-q0));
%dl3=1*(z(5)+z(6)-2.5);
dl(3,:)=((l(3)>0)|(dl3>0))*dl3;

dl4=30*(z(7)+z(8)+u4*[1 -1]*[z(7);z(8)]-2.5-v4*(u4^2-q0));
%dl4=10*(z(7)+z(8)-2.5);
dl(4,:)=((l(4)>0)|(dl4>0))*dl4;

du4=10*([1 -1]*[z(7);z(8)]-2*v4*u4);
dv4=70*((v4>0)|((l(4)+cc)*(u4^2-q0)>0))*(l(4)+cc)*(u4^2-q0);
du3=10*([1 -1]*[z(5);z(6)]-2*v3*u3);
dv3=70*((v3>0)|((l(3)+cc)*(u3^2-q0)>0))*(l(3)+cc)*(u3^2-q0);
du2=10*([1 -1]*[z(3);z(4)]-2*v2*u2);
dv2=70*((v2>0)|((l(2)+cc)*(u2^2-q0)>0))*(l(2)+cc)*(u2^2-q0);
du1=10*([1 -1]*[z(1);z(2)]-2*v1*u1);
dv1=70*((v1>0)|((l(1)+cc)*(u1^2-q0)>0))*(l(1)+cc)*(u1^2-q0);
dz=1*(Z1*z+Z2*c)-10*kron((l+cc),[1;1])-30*[1;-1;0;0;0;0;0;0]*u1*(l(1)+cc)-30*[0;0;1;-1;0;0;0;0]*u2*(l(2)+cc)-30*[0;0;0;0;1;-1;0;0]*u3*(l(3)+cc)-30*[0;0;0;0;0;0;1;-1]*u4*(l(4)+cc);


dw=[dz;dl;du1;dv1;du2;dv2;du3;dv3;du4;dv4];
end