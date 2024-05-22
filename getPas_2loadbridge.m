function [passiveElms]=getPas_2loadbridge(nelX,nelY,sc)
sizeBlockLoad=1/20; nelY=nelY*sc; nelX=nelX*sc;
yInd = nelY+1-(1:sizeBlockLoad*nelY);
xbc_ind1 = 1:(sizeBlockLoad*nelX);           % left corner
xbc_ind2 = (nelX+1-(xbc_ind1));      % right corner
xl_ind1 = xbc_ind1 + (1/4-sizeBlockLoad/2)*nelX; % load point 1
xl_ind2 = xbc_ind1 + (3/4-sizeBlockLoad/2)*nelX; % load point 2

xInd = [xbc_ind1(:);xbc_ind2(:);xl_ind1(:);xl_ind2(:)];
passiveElms = reshape((xInd(:)-1)*nelY+yInd(:)',[],1);
%rhoPassive = ones(size(indPassive,1),1);
end