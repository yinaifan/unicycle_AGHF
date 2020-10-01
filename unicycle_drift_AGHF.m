function [p, dp, xfull, xdrift, bufull, budrift, sol, t, x, xnew, cost] = unicycle_drift_AGHF(tmax, xpoints, intgrids, tpoints, k, X0, Xf, T, N, M)
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
% please modify this function if system is changed
[get_EL, get_G, get_f, get_F] = generate_fh_of_model(N,M); 

opts = odeset('AbsTol',1e-4);
tic;
display('solving PDE...');
% Solve AGHF
sol = pdepe(m,@(x,t,u,DuDx) mypdexpde(x,t,u,DuDx,k,N),@(x) mypdexic(x, X0, Xf, T),@(xl,ul,xr,ur,t) mypdexbc(xl,ul,xr,ur,t, X0, Xf, N),x,t,opts);
% The solution "sol" is of form sol(t,x,i)
toc;

%------------------ Extract Controls and state trajectory------------------
tic;
display('Extracting Trajectory...');

% # of integration grids
xpoints_new = intgrids;
% initialize controls and states
budrift = zeros(M,xpoints_new); % actual control
bufull = zeros(N,xpoints_new); % full control (control on unadmissible directions included!!!)
xnew = linspace(0,T,xpoints_new); % time grids for integration

% use spline to interpolate state from steady state pde solution, p
p = zeros(N,xpoints_new);
for i = 1:N
    p(i,:)=spline(x,sol(end,:,i),xnew);
end
% initial value of drift term f(x(0))
drift0 = get_f(p(:,1));
% manualy defferentiate to get solution of x_dot, represented by dp
dp = [drift0 diff(p.').'/(T/(xpoints_new-1))];

% control extraction
for i = 1 : xpoints_new
    % f(x(t)), note: the F_d(x) in paper
    drift = get_f(p(:,1));
    % get [Fc F], note: the F_bar matrix in paper
    Ffull = get_F(p(:,1));
    % get full control
    bufull(:,i) = Ffull^(-1) * ( dp(:,i) - drift );
    % get actual control
    budrift(:,i)= bufull(N-M+1:end,i);
end


%Integrate to get state trajectory

% initialize state trajctroy with full control
xfull = zeros(N,xpoints_new);
%Initial state
xfull(:,1) = p(:,1); 
% initialize state trajctroy with actual control
xdrift = zeros(N,xpoints_new);
%Initial state
xdrift(:,1)= p(:,1);

% Euler integration to get trajectory (can be replaced by ode45 to get better result)
for i =1: xpoints_new-1
    %[Fc F]
    cBfull = get_F(xfull(:,i));
    % just F
    cBdrift = cBfull(:,(end-M+1):end);
    drift = get_f(xdrift(:,i));
    drift_full = get_f(xfull(:,i));
    xdrift(:,i+1) = xdrift(:,i) + ( drift + cBdrift*budrift(:,i) )*(xnew(i+1)-xnew(i));
    xfull(:,i+1) = xfull(:,i) + ( drift_full + cBfull*bufull(:,i) )*(xnew(i+1)-xnew(i));
end

% initialize cost -- actuated curve length
cost = zeros(1,tpoints);
% calculate the atuated curve length for each grid of s
X_temp = zeros(N,xpoints);
dX_temp = zeros(N,xpoints);
for j = 1:tpoints
    for kk = 1:N
    [X_temp(kk,:), dX_temp(kk,:)]=pdeval(m,x,sol(j,:,kk),x);
    end
    %for each s grid, integrate  (x_dot-f(x))'*G(x_dot-f(x)) to get curve
    %length
    for i = 1:xpoints
        X = X_temp(:,i);
        dX = dX_temp(:,i);
        f = get_f(X);
        G = get_G(X,k);
        cost(j) = cost(j) + (dX-f)'*G*(dX-f)*(T/(xpoints-1));
    end
end
toc;
display('Done!!!!!');

end

% The followings are the PDE, initial condition and boundary condition for pdepe:


function [c,f,s] = mypdexpde(x,t,u,DuDx,k,N)    % Define PDE; right-hand-side of AGHF

global get_EL

% evaluate EL
EL = get_EL(u,DuDx,k);
% dL/dx
pLx = EL(:,1);
% dL/d(dot_x)
pLxd = EL(:,2);

f = pLxd;
s = -pLx;

c = ones(N,1);

% for details of f, s, c, please see the documentation for pdepe

end
% --------------------------------------------------------------------------
 
function u0 = mypdexic(x, X0, Xf, T)    % Initial condition for AGHF

% straight line connecting X0 and Xf
u0=X0+(Xf-X0)*(x/T);
% add some sinusoidal deviation to x1 initial guess (can be discarded)
u0(1)=X0(1)+(Xf(1)-X0(1))*(x/T) + 0.01*sin(2*pi*x/T);

end
% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = mypdexbc(xl,ul,xr,ur,t, X0, Xf, N)    % Boundary condition for AGHF

pl = ul-X0;
ql = zeros(N,1);
pr = ur-Xf;
qr = zeros(N,1);

end

% generate functions for calculating EL equation:
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!! Need to modify this function if system is changed !!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Inputs:
% N,M -- # of states and inputs
% Output:
% get_EL -- function handle for calculationg EL equation, it returns dL/dx
% and % dL/d(dot_x), it takes x and x_dot as input
% get_G -- function handle for calculating the metric G(x)
% get_f -- function handle for calculating the drift term f(x)
% get_F -- function handle for calculating the full control matrix [Fc(x) F(x)]
function [get_EL, get_G, get_f, get_F] = generate_fh_of_model(N,M) 
tic;
display('generating functions for pdepe...');

% ------------------------- system setup ---------------------------------

% make the state and state deirvative simbolic
syms x y th z1 z2 k dx dy dth dz1 dz2 real
X = [x y th z1 z2]';
Xd = [dx dy dth dz1 dz2]';

% drift term f(x) note: Fd(x) is used in the paper
f = [cos(th)*z1 sin(th)*z1 z2 0 0]';
% [Fc F], note: the F_bar matrix in paper
F= sym(eye(N));
% penalty matrix (D matrix in the paper)
D = diag([k*ones(1,N-M) ones(1,M)]);
% the metric, before adding the state constraint barrier functions
H = (F.')^(-1)*D*F^(-1);


% -------------------- state constraint barrier function ------------------
% B is a vector of barrier function for 1 state constraint, each scaler state
% constraint need 1 barrier function
B = [];
% barrier function parameters
kb = .01; % barrier function gain
pb = 1; % order of barrier function

% b is one penalty term for state constraint
%b=0; % no state contraints
b = kb/((2*pi/4)^2 - z2^2)^pb; %contraint on angular velocity, |z2|<pi/2
B = [B b]; % attach the barrier function to B
%b = kb/(2^2 - z1^2)^pb; %contraint on translational velocity, |z1|<2
%b = kb/(z1+2)^pb + kb/(2-z1)^pb; %contraint on translational velocity |z1|<2


% -------------------- Metric and curve Length ----------------------------

% the metric with state constraint barier functions
G =(sum(B)+1)*H;
% actuated curve length
L =  simplify( (Xd - f).' * G * (Xd - f) );

% -------------------- Function Generations -------------------------------

% taking derivatives symbolically to get the EL euqation terms
% pLx -- dL/dx
% pLxd -- dL/d(x_dot)
pLx = sym('pLx', [N 1],'real');
pLxd = sym('pLxd', [N 1],'real');
for i=1:N
    pLx(i) =  diff(L,X(i));
    pLxd(i) = diff(L,Xd(i));
end

% generate functions
get_EL = matlabFunction([pLx, pLxd],'vars', {X,Xd,k});
get_G = matlabFunction(G,'vars', {X,k});
get_f = matlabFunction(f,'vars',{X});
get_F = matlabFunction(F,'vars',{X});

toc;
end