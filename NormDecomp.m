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

syms epsl alph

% =============================================
%  Initialize the sum of squares program and form the variables
prog = sosprogram(vars,epsl);

% Evaluate the correct extended basis for V
basis = extendedbasis(f, vars);

[prog,V] = sospolyvar(prog,basis,'wscoeff');

% Evaluate grad(V)
for iv=1:length(vars)
    gV(iv) = diff(V,vars(iv));
end

% =============================================
% Constraint 1 : positive definiteness of V
% V(x) - epsilon(x1^2 + x2^2 + ... xn^2) >= 0
prog = sosineq(prog,V-epsl*sum(vars.^p));

% Constraint 2: matrix inequality imposes
% grad(U).gU <= 0. which also imposes grad(U).f <=0
expr = [-gV*f,           gV;
        gV.', eye(length(vars))];
prog = sosmatrixineq(prog,expr);

% =============================================
% Set objective to maximise epsilon then call solver
prog = sossetobj(prog, -epsl);
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
    
%     basis = 1.0;
% 
%     % Loop over the elements of f
%     for (i,fi) in enumerate(f)
%         fTmp = [];
%         for elem in fi
%             fTmp += coefficient(elem)*elem;
%         end
%         basis = basis + fTmp*x[i];
%     end
%     basis = monomials(basis);

end


function [basis] = extendedbasis(f, vars)

    % Start from the minimal basis
    basis = minimalbasis(f,vars);
    disp(basis)
    n = length(vars);

    % Find maximum total degree
    d = feval(symengine, 'degree', basis);
    disp(d)
    
    % Now find maximum individual degree for each xᵢ, o[ii]
    o = zeros(n,1);
    for ii=1:n
        o(ii) = feval(symengine, 'degree', basis, vars(ii));
    end
% 
%     % Now find maximum individual degree for each xᵢ, o[ii]
%     o = zeros(Int,1,n)
%     for ii=1:n
%         o[ii] = maximum([degree(Vi,x[ii]) for Vi in basis]);
%     end
% 
    % Add all mixed terms up to total degree d, and with maximum individual degree o[ii]
    basis = monomials(vars,0:d);
    for ii=1:length(basis)
        for ix=1:length(vars)
            if feval(symengine, 'degree', basis(ii), vars(ix)) > o(ix)
                basis(ii) = 0;
            end
        end
    end
    basis = basis(basis~=0);
    disp('Found minimal basis as:')
    disp(basis.')
% 
end

