function [p, dp, xfull, xdrift, bufull, budrift, sol, t, x, xnew, cost] = unicycle_drift_AGHF(tmax, xpoints, intgrids, tpoints, k, X0, Xf, T)
% unicycle_drift_AGHF Solve the AGHF for dynamical unicycle
% code auther: Yinai Fan
%
% system: dot_x = f(x) + F(x)*u
% Lagragian of actuated curve length: (x_dot-f(x))'*G(x_dot-f(x))
%
% The AGHF is equivalent as the Euler Lagrange equation:
% dx/ds = AGHF  <===> dx/ds = - (dL/dx-d(dL/d(dot_x))/dt)
% In this implementation, the EL equation is used for solving PDE
%
% Inputs:
% tmax -- smax
% xpoints -- number of grids in time t
% intgrids -- number of grids for integration of solution
% tpoints -- number of grids in s
% k -- penalty (lambda is used in paper)
% X0, Xf -- initial and final value of states
% T -- motion duration
%
% Outputs:
% p, dp -- steady state PDE solution for x, dot_x of the AGHF
% xfull -- state trajectory integrated by using full control (virtual control included)
% xdrift -- solution of the motion, state trajectory integrated by using actual control
% bufull -- full control (virtual control included)
% budrift -- actual control
% sol -- raw solution of pdepe
% t -- s grids vector for solving AGHF,  x -- time grids vector for solving AGHF
% xnew -- time vector for integration of solutions
% cost -- actuated curve length values along s grids



%----------------------- pdepe setup ------------------------------
m = 0;
x = linspace(0,T,xpoints);                  % discretization of the curve in time
t = [0 logspace(-4,log10(tmax),tpoints-1)]; % discretization of s interval in log scale

% global variable for the Euler Lagrange function and Metric G
global get_EL get_G

% generate functions for EL, metric G, drift term f, control matrix F
[get_EL, get_G, get_f, get_F] = generate_fh_of_model(); 

opts = odeset('AbsTol',1e-4);
tic;
display('solving PDE...');
% Solve AGHF
sol = pdepe(m,@(x,t,u,DuDx) mypdexpde(x,t,u,DuDx,k),@(x) mypdexic(x, X0, Xf, T),@(xl,ul,xr,ur,t) mypdexbc(xl,ul,xr,ur,t, X0, Xf),x,t,opts);
% The solution "sol" is of form sol(t,x,i)
toc;


%------------------ Extract Controls and state trajectory------------------
tic;
display('Extracting Trajectory...');

% # of integration grids
xpoints_new = intgrids;
% initialize controls and states
budrift = zeros(2,xpoints_new); % actual control
bufull = zeros(5,xpoints_new); % full control (virtual control included)
xnew = linspace(0,T,xpoints_new); % time grids for integration

% use spline to interpolate state from steady state pde solution
p1=spline(x,sol(end,:,1),xnew);
p2=spline(x,sol(end,:,2),xnew);
p3=spline(x,sol(end,:,3),xnew);
p4=spline(x,sol(end,:,4),xnew);
p5=spline(x,sol(end,:,5),xnew);
p = [p1;p2;p3;p4;p5];
% initial value of drift term f(x(0))
drift0 = get_f(p(1,1), p(2,1), p(3,1), p(4,1), p(5,1));
% manualy defferentiate to get solution of x_dot
dp = [drift0 diff(p.').'/(T/(xpoints_new-1))];

% control extraction
for i = 1 : xpoints_new
    % f(x(t))
    drift = get_f(p(1,i),p(2,i),p(3,i),p(4,i),p(5,i));
    % get [Fc F]
    Ffull = get_F(p(1,i),p(2,i),p(3,i),p(4,i),p(5,i));
    % get full control
    bufull(:,i) = Ffull^(-1) * ( dp(:,i) - drift );
    % get actual control
    budrift(:,i)= bufull(4:5,i);
end


%Integrate to get state trajectory

% initialize state trajctroy with full control
xfull   =   zeros(5,xpoints_new);
%Initial state
xfull(1,1)  =   p1(1); 
xfull(2,1)  =   p2(1);
xfull(3,1)  =   p3(1);
xfull(4,1)  =   p4(1);
xfull(5,1)  =   p5(1);
% initialize state trajctroy with actual control
xdrift   =   zeros(5,xpoints_new);
%Initial state
xdrift(1,1)  =   p1(1);
xdrift(2,1)  =   p2(1);
xdrift(3,1)  =   p3(1);
xdrift(4,1)  =   p4(1);
xdrift(5,1)  =   p5(1);

% Euler integration to get trajectory (can be replaced by ode45 to get better result)
for i =1: xpoints_new-1
    cBfull = get_F(xfull(1,i),xfull(2,i),xfull(3,i),xfull(4,i),xfull(5,i));
    cBdrift = [0 0 0 1 0;
               0 0 0 0 1]';
    drift = get_f(xdrift(1,i),xdrift(2,i),xdrift(3,i),xdrift(4,i),xdrift(5,i));
    drift_full = get_f(xfull(1,i),xfull(2,i),xfull(3,i),xfull(4,i),xfull(5,i));
    xdrift(:,i+1) = xdrift(:,i) + ( drift + cBdrift*budrift(:,i) )*(xnew(i+1)-xnew(i));
    xfull(:,i+1) = xfull(:,i) + ( drift_full + cBfull*bufull(:,i) )*(xnew(i+1)-xnew(i));
end

% initialize cost--actuated curve length
cost = zeros(1,tpoints);
% calculate the atuated curve length for each grid of s
for j = 1:tpoints
    [x_temp, dx_temp]=pdeval(m,x,sol(j,:,1),x);
    [y_temp, dy_temp]=pdeval(m,x,sol(j,:,2),x);
    [th_temp, dth_temp]=pdeval(m,x,sol(j,:,3),x);
    [z1_temp, dz1_temp]=pdeval(m,x,sol(j,:,4),x);
    [z2_temp, dz2_temp]=pdeval(m,x,sol(j,:,5),x);
    
    %for each s grid, integrate  (x_dot-f(x))'*G(x_dot-f(x)) to get curve
    %length
    for i = 1:xpoints
        X = [x_temp(i) y_temp(i) th_temp(i) z1_temp(i) z2_temp(i)]';
        dX = [dx_temp(i) dy_temp(i) dth_temp(i) dz1_temp(i) dz2_temp(i)]';
        f = get_f(X(1),X(2),X(3),X(4),X(5));
        G = get_G(X(1),X(2),X(3),X(4),X(5),k);
        cost(j) = cost(j) + (dX-f)'*G*(dX-f)*(T/(xpoints-1));
    end
end
toc;
display('Done!!!!!');

end

% The followings are the PDE, initial condition and boundary condition for pdepe:


function [c,f,s] = mypdexpde(x,t,u,DuDx,k)    % Define PDE; right-hand-side of AGHF

global get_EL

% evaluate EL
EL = get_EL(u(1),u(2),u(3),u(4),u(5),DuDx(1),DuDx(2),DuDx(3),DuDx(4),DuDx(5),k);
% dL/dx
pLx = EL(:,1);
% dL/d(dot_x)
pLxd = EL(:,2);

f = pLxd;
s = -pLx;

c = [1;1;1;1;1];

% for details of f, s, c, please see the documentation for pdepe

end
% --------------------------------------------------------------------------
 
function u0 = mypdexic(x, X0, Xf, T)    % Initial condition for AGHF

% straight line connecting X0 and Xf
u0=X0+(Xf-X0)*(x/T);
% add some sinusoidal deviation to x1 initial guess
u0(1)=X0(1)+(Xf(1)-X0(1))*(x/T) + 0.01*sin(2*pi*x/T);

end
% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = mypdexbc(xl,ul,xr,ur,t, X0, Xf)    % Boundary condition for AGHF

pl = ul-X0;
ql = [0;0;0;0;0];
pr = ur-Xf;
qr = [0;0;0;0;0];

end

% generate functions for calculating EL equation:
% Output:
% get_EL -- function handle for calculationg EL equation, it returns dL/dx
% and % dL/d(dot_x), it takes x and x_dot as input
% get_G -- function handle for calculating the metric G(x)
% get_f -- function handle for calculating the drift term f(x)
% get_F -- function handle for calculating the full control matrix [Fc(x) F(x)]
function [get_EL, get_G, get_f, get_F] = generate_fh_of_model() 
tic;
display('generating functions for pdepe...');

syms x y th z1 z2 k dx dy dth dz1 dz2 real
X = [x y th z1 z2]';
Xd = [dx dy dth dz1 dz2]';

% drift term
f = [cos(th)*z1 sin(th)*z1 z2 0 0]';
% [Fc F]
F= eye(5);

kb = .01; % barrier function gain
pb = 1; % order of barrier function
% penalty matrix (D matrix in the paper)
K = diag([k k k 1 1]);

% b is the penalty term for state constraint
%b=0; % no state contraints
b = kb/((2*pi/4)^2 - z2^2)^pb; %contraint on angular velocity, |z2|<pi/2
%b = kb/(2^2 - z1^2)^pb; %contraint on translational velocity, |z1|<2
%b = kb/(z1+2)^pb + kb/(2-z1)^pb; %contraint on translational velocity |z1|<2

% actualted curve length
H = (F.')^(-1)*K*F^(-1);
% actualted curve length with state constraints penalty
G =(b+1)*H;
L =  simplify( (Xd - f).' * G * (Xd - f) );

% taking derivatives symbolically to get the EL euqation terms
pLx = [diff(L,X(1)) diff(L,X(2)) diff(L,X(3)) diff(L,X(4)) diff(L,X(5))]';
pLxd = [diff(L,Xd(1)) diff(L,Xd(2)) diff(L,Xd(3)) diff(L,Xd(4)) diff(L,Xd(5))]';

% generate functions
get_EL = matlabFunction([pLx, pLxd],'vars', [X',Xd',k]);
get_G = matlabFunction(G,'vars', [X',k]);
get_f = matlabFunction(f,'vars',X');
get_F = matlabFunction(F,'vars',X');

toc;
end