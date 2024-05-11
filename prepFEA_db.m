function [passiveElms,edofMat,iK,jK,U,F,freeDofs,Kp] = prepFEA_db(nelX,nelY,sc)
nodeNumbers = reshape(1:(1+nelX)*(1+nelY),1+nelY,1+nelX);                  % Get edofmat, stiffness matrix triplets, and initialise U
edofMat = repmat(reshape(2*nodeNumbers(1:end-1,1:end-1)+1,[],1),1,8)+...
    repmat([0,1,2*nelY+[2,3,0,1],-2,-1],nelX*nelY,1);
iK = reshape(kron(edofMat,ones(8,1))',[],1); 
jK = reshape(kron(edofMat,ones(1,8))',[],1);
numDof = 2*(nelY+1)*(nelX+1); allDofs = 1:numDof; U = zeros(numDof,1);
% ---------------------------------------------- Define loads and supports 
[px,py] = deal(ceil(nelX/(40*sc))*sc,ceil(nelY/(10*sc))*sc);               % Passive block relative sizes
elementX = ceil(nelX/2-px/2)+(1:px);                                       % Elements in X dir subjected to traction
elementY = ceil(nelY/2-py/2)+(1:py);                                       % Elements in Y dir subjected to traction
[elementYGrid, elementXGrid] = ndgrid(elementY,elementX);                  % Elements grid subjected to traction
forceElms = (elementXGrid(:)-1)*nelY+elementYGrid(:);                      % Elements subjected to traction
sF = repmat(-1*[0 1 0 1 0 1 0 1],[numel(forceElms),1]);                    % Integrate traction over element surface
iF = edofMat(forceElms,:);                                                 % Nodes subjected to traction
F = sparse(iF(:),ones(size(iF(:))),sF(:),numDof,1); F = F/sum(abs(F(:)));  % Build RHS force vector with unit magnitude
Kp = sparse(numDof,numDof);                                                % Empty multipoint constraint penalty stiffness matrix
fixedDofs = [1:2*(nelY+1),2*(nelY+1)*(nelX+1)-2*(nelY+1) + (1:2*(nelY+1))];% Fixed dof in each end of cantilever
freeDofs = setdiff(allDofs,fixedDofs);                                     % Get free dofs from fixed dofs
% ------------------- Set passive elements at load and rollor support point
passiveElms = forceElms;                                                   % Compute passive elements from grid ids
end