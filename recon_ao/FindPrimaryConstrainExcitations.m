function primaryOut = FindPrimaryConstrainExcitations(primaryIn,M_PrimaryToExcitations,primaryHeadroom,maximizeVec,constraintEqA,constraintEqb,lambda,options)
% Find display primaries that satisfy constraints on cone excitations
%
% Syntax:
%
% Description:
%
%   Bounds on primaries are [0,1].
%
% Inputs:
%
% Outputs:
%
% Optional key/value pairs
%
% See also:

% History
%
% 9/16/22  dhb  Wrote initial version

% Examples
%{
    primaryIn = [0.5 0.5 0.5]';
    M_PrimaryToExcitations = ...
       [0.4155    1.0703    0.1670
        0.1265    0.9730    0.2068
        0.0039    0.0155    0.1559];
    maximizeVec = [1 0 0];
    constraintEqA = [1 -1 0 ; 0 0 1];
    constraintEqb = [0 0.1]';
    lambda = 0.01;
    primaryHeadroom = 0.05;
    primaryOut = ...
        FindPrimaryConstrainExcitations(primaryIn,M_PrimaryToExcitations,primaryHeadroom,maximizeVec,constraintEqA,constraintEqb,lambda);
    excitationsOut = M_PrimaryToExcitations*primaryOut;
    fprintf('Primaries input: %0.4f, %0.4f, %0.4f\n',primaryIn(1),primaryIn(2),primaryIn(3));
    fprintf('Primaries out:   %0.4f, %0.4f, %0.4f\n',primaryOut(1),primaryOut(2),primaryOut(3));
    fprintf('Excitations out: %0.4f, %0.4f, %0.4f\n',excitationsOut(1),excitationsOut(2),excitationsOut(3));
    fprintf('L - M (expect 0): %0.4f\n',excitationsOut(1)-excitationsOut(2));
%}

foptions = optimset('fmincon');
foptions = optimset(foptions,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
if (any(primaryIn < primaryHeadroom) || any(primaryIn > 1-primaryHeadroom))
    error('Passed primaries not in gamut');
end

x0 = primaryIn;
vlb = zeros(size(x0)) + primaryHeadroom;
vub = ones(size(x0)) - primaryHeadroom;

primaryOut0 = fmincon(@(x)FitFunction(x,M_PrimaryToExcitations,maximizeVec,constraintEqA,constraintEqb,lambda),x0,[],[],[],[],vlb,vub,[],foptions);
primaryOut = fmincon(@(x)FitFunction(x,M_PrimaryToExcitations,maximizeVec,constraintEqA,constraintEqb,lambda),primaryOut0,[],[],[],[],vlb,vub,[],foptions);

end

function f = FitFunction(x,M_PrimaryToExcitations,maximizeVec,constraintEqMat,constraintEqb,lambda)

weightedExcitations = maximizeVec*M_PrimaryToExcitations*x;

constraintVals = constraintEqMat*M_PrimaryToExcitations*x-constraintEqb;
constraintVal = sqrt(mean(constraintVals(:).^2));

f = -1000*(lambda*weightedExcitations - (1-lambda)*constraintVal);

end
