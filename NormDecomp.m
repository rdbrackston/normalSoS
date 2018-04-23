function [ U ] = NormDecomp( f, vars, varargin )
%NormDecomp Function to perform the normal decomposition of a vector field
%   This function uses the Sum of Squares (SOS) optimisation method to
%   generate a Lyapunov function U, such that the vector field f is
%   decomposed into two orthogonal components: f = -grad(U) + g.

% Use varargin to set the number of iterations
if nargin == 3
    iters = varargin{1};
    p = 2;
elseif nargin == 4
    iters = varargin{1};
    p = varargin{2};
else
    iters = 1;
    p = 2;
end

syms alph epsl
n = length(vars);

% =============================================
%  Initialize the sum of squares program and form the variables
prog = sosprogram(vars);

% Evaluate the correct extended basis for V
basis = minimalbasis(f, vars);
% basis = extendedbasis(f, vars);
disp('Chosen basis as:')
disp(basis.')


[prog,V] = sospolyvar(prog,basis,'wscoeff');

% =============================================
%  Set up the lower bounding polynomial, bnd
o = zeros(n,1);
if n>1
    bnd = vars(1);    % Initialise the type of bnd
    % First add the maximum even order polynomial for each x[ii]
    for ii=1:n
        for ib=1:length(basis)
            tmp = double(feval(symengine,'degree',basis(ib),vars(ii)));
            if mod(tmp,2)==0 && tmp>0
                o(ii) = max([tmp ,o(ii)]);
            end
        end
        bnd = bnd + vars(ii)^o(ii);
    end
    
    % Now add any fully even terms of mixed x[ii], e.g. x[1]^2*x[2]^2
    for ib=1:length(basis)
        iter = 0;
        for ii=1:n
            tmp = double(feval(symengine,'degree',basis(ib),vars(ii)));
            if (mod(tmp,2)==0 && tmp>0)
                iter = iter + 1;
            end
        end
        if iter>1;    bnd = bnd + basis(ib);    end
    end
    bnd = bnd - vars(1);
    [~,trms] = coeffs(bnd);
    bnd = monomials(trms,1);
    b = length(bnd);
else
    b = 1;
    for ib=1:length(basis)
        o(1) = max([feval(symengine,'degree',basis(ib),vars(1)),o(1)]);
    end
    bnd = vars(1)^o(1);
    bnd = [bnd];
end
disp('Chosen bound as:')
disp(bnd.')

% Evaluate grad(V)
for iv=1:length(vars)
    gV(iv) = diff(V,vars(iv));
end


% =============================================
% Constraint 1 : positive definiteness of V
% V(x) - E*bnd >= 0
[prog,E] = sospolymatrixvar(prog,monomials(vars(1),0),[1,b]);
prog = sosineq(prog,V-E*bnd);
for ib=1:b
    prog = sosineq(prog,E(b));
end

% Constraint 2: matrix inequality imposes
% grad(U).gU <= 0. which also imposes grad(U).f <=0
expr = [-gV*f,           gV;
        gV.', eye(length(vars))];
prog = sosmatrixineq(prog,expr);

% =============================================
% Set objective to maximise sum over E then call solver
prog = sossetobj(prog, -sum(E));
prog = sossolve(prog);

U = sosgetsol(prog,V);


%% Now iteratively improve it

for ii=1:iters

    % =============================================
    %  Initialize the second sum of squares program and form the variables
    prog2 = sosprogram(vars, [epsl,alph]);
    [prog2, V] = sospolyvar(prog2, basis,'wscoeff');

    % Evaluate grad(V) and grad(U)
    for iv=1:length(vars)
        gV(iv) = diff(V,vars(iv));
        gU(iv) = diff(U,vars(iv));
    end

    % Constraint 1: V(x) - epsilon(x1^2 + x2^2 + ...) >= 0. For positive definiteness of V
    prog2  = sosineq(prog2, V-epsl*sum(vars.^p));

    % Constraint 2: matrix inequality imposes
    % grad(U).gU <= 0. which also imposes grad(U).f <=0
    expr = [-gV*f,           gV;
            gV.', eye(length(vars))];
    prog2 = sosmatrixineq(prog2,expr);

    % Constraint 3: Wynn inequality
    expr2 = gV*(f+2*gU.') - alph*gU*f - (1+alph)*sum(gU.^2);
    prog2 = sosineq(prog2,expr2);

    % Constraint 4: Alpha, gamma greater than zero
    prog2 = sosineq(prog2,alph);
    prog2 = sosineq(prog2,epsl);

    % Set objective and solve
    prog2 = sossetobj(prog2, alph);
    prog2 = sossolve(prog2);

    fprintf('Alpha = %e',sosgetsol(prog2,alph));
    U = sosgetsol(prog2, V);

end

end


function [basis] = minimalbasis(f, vars)

    % Evaluate the correct minimal basis for V
    basis = 1;
    % For each function f_i, find the polynomial terms and multiply by x_i
    for iv=1:length(vars)
        [coeff,monom] = coeffs(f(iv),vars,'all');
        tmp = coeff.*monom.*vars(iv);
        basis = basis + sum(reshape(tmp,[],1));
    end
    % Now reduce the basis to a column vector with all coefficients equal to 1
    [coeff,monom] = coeffs(basis,vars,'all');
    coeff = double(logical(coeff));
    basis = coeff.*monom;
    basis = reshape(basis,[],1);
    basis = basis(basis~=0);

end


function [basis] = extendedbasis(f, vars)

    % Start from the minimal basis
    basis = minimalbasis(f,vars);
    n = length(vars);

    % Find maximum total degree
    d = 0;
    for ib=1:length(basis)
        d = max([double(feval(symengine, 'degree', basis(ib))), d]);
    end
    
    % Now find maximum individual degree for each xᵢ, o[ii]
    o = zeros(n,1);
    for ii=1:n
        for ib=1:length(basis)
            o(ii) = max([double(feval(symengine,'degree',basis(ib),vars(ii)))...
                            ,o(ii)]);
        end
    end

    % Add all mixed terms up to total degree d, and with maximum individual degree o[ii]
    basis = monomials(vars,0:d);
    for ii=1:length(basis)
        for ix=1:n
            if feval(symengine, 'degree', basis(ii), vars(ix)) > o(ix)
                basis(ii) = 0;
            end
        end
    end
    basis = basis(basis~=0);
% 
end

