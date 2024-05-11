function [passiveElms,edofMat,iK,jK,U,F,freeDofs,Kp]=prepFEA_mbb(nelX,nelY,sc) 
nodeNumbers = reshape(1:(1+nelX)*(1+nelY),1+nelY,1+nelX);                  % Get edofmat, stiffness matrix triplets, and initialise U
edofMat = repmat(reshape(2*nodeNumbers(1:end-1,1:end-1)+1,[],1),1,8)+...
    repmat([0,1,2*nelY+[2,3,0,1],-2,-1],nelX*nelY,1);
iK = reshape(kron(edofMat,ones(8,1))',[],1); 
jK = reshape(kron(edofMat,ones(1,8))',[],1);
numDof= 2*(nelY+1)*(nelX+1); allDofs = 1:numDof; U = zeros(numDof,1);
% ---------------------------------------------- Define loads and supports 
[px,py] = deal(ceil(nelX/(30*sc))*sc,ceil(nelY/(30*sc))*sc);               % Passive block sizes (coarse scaled)
forceElms = 1:nelY:px*nelY;                                                % Elements subjected to traction
iF = edofMat(forceElms,:);                                                 % Nodes subjected to traction
sF = repmat(-1*[0,0,0,0,0,1,0,1],[numel(forceElms),1]);                    % Integrate traction over element surface
F = sparse(iF(:),ones(size(iF(:))),sF(:),numDof,1); F = F/sum(abs(F(:)));  % Build RHS force vector with unit magnitude
[iB,jB] = meshgrid(numDof -flip(0:2*(nelY+1):2*(nelY+1)*px)-2*(nelY+1)*px);% Distributed Y-roller support (multipoint constraint)
Kp = sparse(iB(:),jB(:),1e4,numDof,numDof);                                %    with penalty method (penalty value = 1e4)  
fixedDofs = 1:2:2*(nelY+1);                                                % Symmetry plane
freeDofs = setdiff(allDofs,fixedDofs);                                     % Get free dofs from fixed dofs
% ------------------- Set passive elements at load and rollor support point
[yInd1, xInd1] = ndgrid(1:py,1:px);                                        % Load area (element grid ids)
[yInd2, xInd2] = ndgrid(nelY+1-(1:py),sort(nelX+1-(xInd1)) - px);          % Support area (element grid ids)
passiveElms = ([xInd1(:);xInd2(:)]-1)*nelY+[yInd1(:); yInd2(:)];           % Compute passive elements from grid ids
end