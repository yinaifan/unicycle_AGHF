%% Main for dynamcal unicycle motion generation using AGHF

%--------------------------------------------------------------------------
% see paper "AffineGeometricHeatFlowandMotionPlanningforDynamicSystems" ShenyuLiu,YinaiFan,Mohamed-AliBelabbas
% 11th IFAC Symposium on Nonlinear Control Systems NOLCOS 2019
%
% code auther: Yinai Fan
%--------------------------------------------------------------------------

clear;
clc;

%---------------------------- AGHF parameters ----------------------------

% smax for steady state
smax = 0.002;
% # of grids in t
tgrids = 100;
% # of grids in s
sgrids = 100;
% # of grids of integration for final trajectory extraction
intgrids = 10000;
% initial and final value for states
X0 = [0;0;0;0;0];
Xf = [0;-1;0;0;0];
N = length(X0);
M =2;
% motion duration
T = 1;
% penalty value (the lambda in paper)
k = 50000;

% solve for trajectory, see "unicycle_drift_AGHF" function file for
% implementation details
[p, dp, xfull, xdrift, bufull, budrift, sol, s, t, tint, cost] = unicycle_drift_AGHF(smax, tgrids, intgrids, sgrids, k, X0, Xf, T, N, M);

%% plotting

figure('Name','actuated curve length vs s')
plot(s,cost);
title('actuated curve length vs s')
xlabel('s')
ylabel('curve length')
grid on;

u1 = sol(:,:,1);    
u2 = sol(:,:,2);    
u3 = sol(:,:,3);    
u4 = sol(:,:,4);    
u5 = sol(:,:,5);   

figure('Name','3D configuration space curve deformation')
X=u1(1,:);
Y=u2(1,:);
Z=u3(1,:);
h1=plot3(X,Y,Z,'k','LineWidth',2);
hold on
plot3(X,Y,Z,'r','LineWidth',2);
legend('x-y-\theta trajectory','initial guess')
axis([min(u1(end,:))-0.1,max(u1(end,:))+0.1,min(u2(end,:))-0.1,max(u2(end,:))+0.1,min(u3(end,:))-0.1,max(u3(end,:))+0.1]);
grid ON;
title('3D configuration space curve deformation');
xlabel('x');
ylabel('y');
zlabel('\theta');
%pause;

for i=1:sgrids
    X=u1(i,:);h1.XDataSource='X';
    Y=u2(i,:);h1.YDataSource='Y';
    Z=u3(i,:);h1.ZDataSource='Z';
    refreshdata(h1,'caller');
    drawnow;
    pause(0.05);
end

figure('Name','xy-plane trajectory deformation')

X=u1(1,:);
Y=u2(1,:);
h1=plot(X,Y,'k','LineWidth',2);
hold on
plot(X,Y,'r','LineWidth',2);
legend('x-y trajectory','initial guess')
axis([min(u1(end,:))-0.1,max(u1(end,:))+0.1,min(u2(end,:))-0.1,max(u2(end,:))+0.1]);
grid ON;
title('xy-plane trajectory and unicycle deformation');
xlabel('x');
ylabel('y');
%pause;

for i=1:sgrids
    
    X=u1(i,:);h1.XDataSource='X';
    Y=u2(i,:);h1.YDataSource='Y';
    
    refreshdata(h1,'caller');
    drawnow;
    pause(0.05);
    
end

%% animation

figure('Name','unicycle animation');
plot(p(1,:),p(2,:),'k','LineWidth',2);
hold on;

plot(interp1(tint,xdrift(1,:),t),interp1(tint,xdrift(2,:),t),':c','LineWidth',2);

legend('AGHF solution', 'integrated path','AutoUpdate','off')
%plot actual trajectory
Xs=interp1(tint,xdrift(1,:),t);
Ys=interp1(tint,xdrift(2,:),t);
Zs=interp1(tint,xdrift(3,:),t);
Z1s=interp1(tint,xdrift(4,:),t);
Z2s=interp1(tint,xdrift(5,:),t);

w = 0.2;
h = 0.1;
Rot_s = [cos(Zs(1)) -sin(Zs(1));
       sin(Zs(1)) cos(Zs(1))];
p1s = [Xs(1),Ys(1)]' + Rot_s*[w/2;h/2];
p2s = [Xs(1),Ys(1)]' + Rot_s*[-w/2;h/2];
p3s = [Xs(1),Ys(1)]' + Rot_s*[-w/2;-h/2];
p4s = [Xs(1),Ys(1)]' + Rot_s*[w/2;-h/2];
h_box_patch=patch([p1s(1) p2s(1) p3s(1) p4s(1)], [p1s(2) p2s(2) p3s(2) p4s(2)],[0.4,0.4,1]);
h_v=quiver(Xs(1),Ys(1),Z1s(1)*cos(Zs(1)),Z1s(1)*sin(Zs(1)),'r');
qfactor = 0.05;
qhead = 10;
qwidth = 1.5;
set(h_v,'AutoScale','off','AutoScaleFactor',qfactor,'MaxHeadSize',qhead,'linewidth',qwidth)

axis equal
axis([-1.1,1.1,-1.6,0.6]);
grid ON;
xlabel('x');
ylabel('y');
%pause;

ind_plot = 0:10:tgrids;
ind_plot(1) = 1;
ind_count = 1;
for i=1:tgrids
    delete(h_box_patch);
    delete(h_v);
    Rot_s = [cos(Zs(i)) -sin(Zs(i));
        sin(Zs(i)) cos(Zs(i))];
    p1s = [Xs(i),Ys(i)]' + Rot_s*[w/2;h/2];
    p2s = [Xs(i),Ys(i)]' + Rot_s*[-w/2;h/2];
    p3s = [Xs(i),Ys(i)]' + Rot_s*[-w/2;-h/2];
    p4s = [Xs(i),Ys(i)]' + Rot_s*[w/2;-h/2];
    pxs = [p1s(1) p2s(1) p3s(1) p4s(1) p1s(1)];
    pys = [p1s(2) p2s(2) p3s(2) p4s(2) p1s(2)];
    
    h_box_patch=patch(pxs(1:4),pys(1:4),[0.4,0.4,1]);
    h_v=quiver(Xs(i),Ys(i),Z1s(i)*cos(Zs(i)),Z1s(i)*sin(Zs(i)),'r');
    set(h_v,'AutoScaleFactor',qfactor,'MaxHeadSize',qhead,'linewidth',qwidth)
    
    if i == ind_plot(ind_count)
        patch(pxs(1:4),pys(1:4),[0.4,0.4,1],'EdgeAlpha', 0.4+(i/tgrids)*0.6, 'FaceAlpha', 0.4+(i/tgrids)*0.6);
        U = Z1s(i)*cos(Zs(i));
        V = Z1s(i)*sin(Zs(i));
        if (abs(U)>0.0001)&&(abs(V)>0.0001)
            quiver(Xs(i),Ys(i),U,V,'r','AutoScaleFactor',qfactor,'MaxHeadSize',qhead,'linewidth',qwidth);
        end
        ind_count = ind_count + 1;
    end

    drawnow;
    pause(0.005)
end

delete(h_box_patch);
delete(h_v);


