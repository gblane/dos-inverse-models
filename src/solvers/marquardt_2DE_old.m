function [X, info, perf, costf, mod, vid] = marquardt_2DE_old(fun_marq2DE_six_param, fpar, x0, opts, data, flags, en)
% marquardt_2DE_old Inverse two-layer model using Marquardt's method (old version).
%
% [X, info, perf, costf, mod, vid] = marquardt_2DE_old(fun_marq2DE_six_param, fpar, x0, opts, data, flags, en)
%
% Written by Giles Blaney, Ph.D. (Originally by Bertan H.)
%
% Inputs:
%   fun_marq2DE_six_param - Function handle for the forward model [-]
%   fpar                  - Fitting parameters indices [-]
%   x0                    - Initial guess for fitting parameters [mixed]
%   opts                  - Optimization options [-]
%   data                  - Experimental data [mixed]
%   flags                 - Vector of booleans controlling parameters to fit [-]
%   en                    - Zeroth order Bessel function roots [-]
%
% Outputs:
%   X     - Recovered parameters [mixed]
%   info  - Performance information [-]
%   perf  - Performance metrics [-]
%   costf - Cost function value [unitless]
%   mod   - Forward model output [mixed]
%   vid   - Video/visualization information [-]

x0(end+1)=1;

phi_data=data(:, 2);
AC_data=data(:, 3);
phi_err=data(:, 4);
AC_err=data(:, 5);

phi_err=sqrt(phi_err.^2+phi_err(1)^2);
phi_data=phi_data-phi_data(1);

AC_err=sqrt((AC_err/AC_data(1)).^2+(AC_err(1)*AC_data/AC_data(1)^2).^2);
AC_data=AC_data/AC_data(1);

data(:, 2)=phi_data;
data(:, 3)=AC_data;
data(:, 4)=phi_err;
data(:, 5)=AC_err;

costf = zeros(opts(4)+1,1); % Keeping track of the cost function
% vid = zeros(opts(4)+1,1);
%%%%%%
% [x,n,f,J] = checkfun(fun_marq2DE,fpar,x0,opts,data);
[x,n,f,J] = checkfun(fun_marq2DE_six_param,fpar,x0,opts,data,flags,en); %written by Angelo edited by Giles Blaney
costf(1)=(f'*f); %This is the cost function (chisquare) from iteration #1

%  Initial values
A = J'*J;           % Creating a square Jacobian matrix that is symmetric and positive semi-definite
F = (f'*f)/2;       % Madsen Eq (3.1b). Definition of Least squares.
g = J'*f;           % Gradient: Comes from the definition of least squares - dF/dx - Madsen Eq (3.3) and (3.4a)
ng = norm(g,inf);               
mu = opts(1) * max(diag(A));    % initial mu-value is set to max of the elements in the diagonal of A   
kmax = opts(4);                 % iteration limit
Trace = nargout > 2;            
if  Trace                       % if there is more than 2 output arguments (Trace will always be 1)
    X = x*ones(1,kmax+1); X = X';
    perf = [F; ng; mu]*ones(1,kmax+1); perf = perf';
end
k = 1;   
nu = 2;   
nh = 0;   
stop = 0;

while   ~stop
    if  ng <= opts(2)                  
        stop = 1;                        % Stopping criteria with gradient
    else
        h = (A + mu*eye(n))\(-g);        % This is the Tikhonov-type pseudo inverse for LM step (Miller notes p.151)
        nh = norm(h);   
        nx = opts(3) + max(svd(x));      % opts(3) plus maximum singular value of x. max(svd(x)) is the same as norm(x)       
        if nh <= opts(3)*nx  
            stop = 2;                    % Stopping criteria with step-size 
        elseif  nh >= nx/eps   
            stop = 4;                           % Almost singular ?
        end                                     
    end
    if  ~stop
        xnew = x + h;  % This is the Newton's step: xk+1 = xk - f(x)/f'(x) 
        h = xnew - x;   
        dL = (h'*(mu*h - g))/2;          % dL :: gain predicted by the model :: L(0)- L(h) :: given in Eq (3.7b)
%         [fnew,Jnew,mod,fram] = feval(fun_marq2DE, xnew,fpar,data);
        [fnew,Jnew,mod,fram] = feval(fun_marq2DE_six_param, xnew,fpar,data,flags,en);% written by Angelo edited by Giles Blaney
%         [fn,Jn,mod,fram] = fun_marq2DE(xnew,fpar,data);
        Fnew = (fnew'*fnew)/2;   
        dF = F - Fnew;                            % dF :: Gain of the actual function :: F(x)- F(x+h)
        if  (dL > 0) & (dF > 0)                 % Update x and modify mu
            x = xnew;   
            F = Fnew;  
            J = Jnew;  
            f = fnew;
            A = J'*J;   
            g = J'*f;   
            ng = norm(g,inf);
            mu = mu * max(1/3, 1 - (2*dF/dL - 1)^3);   % updating the value of mu and nu based on Madsen 2.21 (p.14)
            nu = 2;
        else                                           % Same  x, increase  mu
            mu = mu*nu;  
            nu = 2*nu;
        end
        k = k + 1;
        costf(k) = F;           % Assigning values to the cost 'vector'
%         vid(k) = fram;
        vid=[];
        if  Trace  
            X(k,:) = x';   
            perf(k,:) = [F ng mu]'; 
        end
        if  k > kmax  
            stop = 3;           % Stopping criteria with max number of iterations 
        end
    end
end
%  Set return values
if  Trace
    X = X(1:k,:);   
    perf = perf(1:k,:);
else  
    X = x;  
end
info = [F  ng  nh  mu/max(diag(A))  k-1  stop];

X=X(:, 1:end-1);

% function  [x,n,f,J] = checkfun(fun_marq2DE,fpar,x0,opts,data)
function  [x,n,f,J] = checkfun(fun_marq2DE_six_param,fpar,x0,opts,data,flags,en) %written by Angelo
% CHECKFUN checks the initialization and sizing (# of rows and columns) of
% each parameter in FUN_MARQ2DE
sx = size(x0);   n = max(sx);
if  (min(sx) > 1)
    error('x0  should be a vector'), end
x = x0(:);   
% [f J] = feval(fun_marq2DE,x,fpar,data);
[f J] = feval(fun_marq2DE_six_param,x,fpar,data,flags,en); %written by Angelo
sf = size(f);   sJ = size(J);
if  sf(2) ~= 1, error('f  must be a column vector'), end
if  sJ(1) ~= sf(1), error('row numbers in  f  and  J  do not match'), end
if  sJ(2) ~= n, error('number of columns in  J  does not match  x'), end
%  Thresholds
if  length(opts) < 4, error('opts  must have 4 elements'), end
if  length(find(opts(1:4) <= 0)), error('The elements in  opts  must be strictly positive'), end
