function [rhoPhys,TO] = deHomTop808(nelX,nelY,volFrac,rMin,wMin,wMax,...
    dMin,deHomFrq,eval,TO)
%% An 808 line dehomogenisation code for multi-scale TopOpt, May 03 2024
fprintf('\nInitialising...'); T = tic;
% ----------------------------------------------------- MATERIAL PROPERTIES
[E,EMin,nu] = deal(1,1e-9,1/3);                                            % Young's modulus of stiff and isotropic background material and Poisson's ratio
% ----------------------------------------- PREPARE FINITE ELEMENT ANALYSIS
[KE0,BMatrixMid] = preIntegrateKE();                                       % Get pre-intergrated stiffness components and stress-displacement matrix
[pasE,edofMat,iK,jK,U,F,freeDofs,Kp] = prepFEA(nelX,nelY,1);               % Get passive elements, edofmat, dof and element info, and penalty matrix
ne = nelX*nelY;                                                            % Number of elements
actE = setdiff(1:ne,pasE);                                                 % Active elements
dv = reshape(actE(:) + ne*(0:3),[],1);                                     % Active design variables
% ------------------------------------ PREPARE FILTER AND ROBUST PROJECTION
rsc = 3.0;                                                                 % Indicator field radius scale wrt. rMin
[dy,dx] = ndgrid(-ceil(rsc*rMin)+1:ceil(rsc*rMin)-1); d= sqrt(dx.^2+dy.^2);% Convolution kernel size
[hS,hW] = deal(max(0,rsc*rMin-d),max(0,rMin-d));                           % Convolution kernel for width and indicator field
HsS = conv2(ones(nelY,nelX),hS,'same');                                    % Matrix of filter weights for indicator
HsW = conv2(ones(nelY,nelX),hW,'same');                                    %    and width field
[dEta,beta,betaFinal] = deal(0.05,0.1,32);                                 % Projection range and steepness; initial and final
eta = [dEta 0 -dEta] + 0.5;                                                % Projection threshold for robust design
dilatedVolfrac = volFrac;                                                  % Initialise dilated volume fraction
% ----------------------- INITIALISE ORIENTATIONS FROM PRINCIPAL DIRECTIONS
C = getIsotropicElasticity(E,EMin,nu,ones(ne, 1));                         % Get elastic properties from solid
U = doFEA(F,U,C,iK,jK,freeDofs,KE0,Kp);                                    % Call FE-analysis to get displacements
stress = (reshape(C(1,1,[1,2,3,2,4,5,3,5,6]),3,3)*BMatrixMid*U(edofMat'))';% Element stress from solid region
a = atan2(2*stress(:,3),(stress(:,1)-stress(:,2)))/2;                      % Principal stress directions
a(not(a > 0)) = a(not(a > 0)) + pi;                                        % Correct angle invariants to [0,2pi]
cScale = volFrac/(F'*U);                                                   % Solid compliance scale
% ------------------------------------------------------- INITIALISE DESIGN
s = ones(ne,1);                                                            % Initialise indicator field
w = max(min((1 - sqrt(1 - volFrac))./s.*ones(ne,2),wMax),wMin);            % Initialise thickness to volume fraction and clamp to bounds
[w(pasE,:),s(pasE,:)] = deal(1.0);                                         % Insert solid material in passive elements
[wTilde,wBar,sTilde,sPhys] = deal(w,w,s,repmat(s,[1,1,3]));                % Initialise modified fields
% ---------------------------------------- DEHOMOGENISATION PRECOMPUTATIONS
deHomGrid = prepPhasor(nelX,nelY,dMin,wMin);                               % Get struct of precomputed sizes, parameters and index sets
align = prepPhasorKernels(deHomGrid);                                      % Get struct of element-wise neighborhood candidates and distances
alignItr = 20;                                                             % Set alignment iterations during optimization
rhoPhys = zeros(deHomGrid.fsize(1),deHomGrid.fsize(2));                    % Initialise the dehomogenised field
pasRhoPhys=prepFEA(deHomGrid.fsize(2),deHomGrid.fsize(1),deHomGrid.fscale);% Get passive elements for dehomogenised field
% --------------------------------------------------------- INITIALISE PLOT
figure; tiledlayout('flow','TileSpacing','tight','Padding','tight');       % Tiled figure
tmpPlt=ones(nelY,nelX); colormap(flip(gray)); nexttile;                    % Temporary plot data
h(1) = imagesc(tmpPlt); hold on;                                           % Rank-2 density plot
[h(2),h(3)]=deal(quiver(tmpPlt,tmpPlt,0),quiver(tmpPlt,tmpPlt,0)); hold off% Quiver plots of tangents
for i=2:3;h(i).LineWidth=2;h(i).Alignment='center';h(i).ShowArrowHead=0;end% Set linewidth, align to centre, and remove head
if deHomFrq>0; nexttile; h(4)=imagesc(rhoPhys); end                        % Dehomgenisation density
for i=1:1+~~deHomFrq;nexttile(i);clim([0 1]);axis off;daspect([1,1,1]); end
% --------------------------------------------- INITIALISE OC AND ITERATION
[x,x0,xOld0,xOld1] = deal([w(:); a(:); s(:)]);                             % Initialise design vector
xmin0 = reshape(repmat([wMin wMin -4*pi 0.0],ne,1),[],1);                  % Min. box constraints
xmax0 = reshape(repmat([wMax wMax  4*pi 1.0],ne,1),[],1);                  % Max. box constraints
xmove = reshape(repmat([0.2  0.2   0.05 0.2],ne,1),[],1);                  % Design variables move limit
sWeight = 0.05;                                                            % Indicator penalty objective weight
stopCrit = 1e-3; loopMax = 300; updateFreq = 50; updateLoop = 0;           % Stop criteria wrt. change, max loops, and continuations loop frequency
if exist('TO','var'); loopMax=0; wPhys=TO.w; N=TO.N; f=TO.f; J=TO.J; end   % If optimisation data is provided, skip
loop = 0; change = 1; Ttot = toc(T); fprintf(' time: %.3f\n',Ttot);        % Misc.
%% ====================================== MULTI-SCALE TOPOLOGY OPTIMISATION
while change > stopCrit && loop < loopMax
    [T,loop,updateLoop] = deal(tic,loop + 1,updateLoop + 1);               % Start optimisation loop
    % -------------------------------- CONTINUATION ON PROJECTION SHARPNESS
    if (beta<betaFinal && (not(updateLoop<updateFreq) || change<stopCrit))
        beta = min(betaFinal,max((beta*2),1)); restartAs = 1;              % Update beta and restart asymptotes
        updateFreq = max(updateFreq - 2,30); updateLoop = 0;               % Update frequency and restart update iterator
        fprintf(' |- beta update: %.2f\n',beta);                           % Print update
    else; restartAs = 0;
    end
    % ------------------------------------------------------------ UPDATE W
    sBar = conv2(reshape(s,[nelY,nelX]),hS,'same')./HsS;                   % Filter indicator field
    wBar(:) = convn(reshape(w,[nelY,nelX,2]),hW,'same')./HsW;              % Filter thickness field
    [wTilde(actE,:),sTilde(actE)] = deal(wBar(actE,:),sBar(actE));         % Modify active elements only
    sPhys(:) = (tanh(beta*eta)+tanh(beta*(sTilde(:)-eta)))./...            % Robust projected indicator field
        (tanh(beta*eta)+tanh(beta*(1-eta)));
    wPhys = wTilde.*sPhys;                                                 % Compute projected thickness field
    % ------------------------------------------------------- GET MU FROM W
    rho = 1 - prod(1 - wPhys,2);                                           % Physical densities
    [muPhys,dmudw] = getMuFromW(rho(:,:,1),wPhys(:,:,1));                  % Compute mu from w (eroded design) and obtain sensitivities
    % --------------------------------------------------------- FE-ANALYSIS
    [C,dCdmu1,dCdmu2,dCda] = getConstitutiveRank2(a,muPhys,E,nu,EMin);     % Get elastic properties and sensitivities wrt. mu1, mu2, and a
    U = doFEA(F,U,C,iK,jK,freeDofs,KE0,Kp);                                % Call FE-analysis to get state field; displacements
    % --------------------- COMPLIANCE OBJECTIVE FUNCTION AND SENSETIVITIES
    J = F'*U;                                                              % Global compliance
    Ue2 = reshape(pagemtimes(reshape(U(edofMat'),8,1,[]),...               % Precompute U(e)'.*U(e) for sum(K(e).*Ue2) e.i. U(e)'*K(e)*U(e)
          reshape(U(edofMat'),1,8,[])),64,[]);
    dKdMu1 = sum(pagemtimes(KE0,dCdmu1),3);                                % Get state sensitivities wrt. mu1, dKdmu1 = sum(KE0i*dCidmu1,i=1..6)
    dKdMu2 = sum(pagemtimes(KE0,dCdmu2),3);                                % Get state sensitivities wrt. mu2, dKdmu2 = sum(KE0i*dCidmu2,i=1..6)
    dKda   = sum(pagemtimes(KE0,dCda),3);                                  % Get state sensitivities wrt. a,   dKda   = sum(KE0i*dCida,i=1..6)y
    dJdmu = -[sum(dKdMu1.*Ue2,1)', sum(dKdMu2.*Ue2,1)'];                   % Compliance sensitivities wrt. mu1 and mu2
    dJda  =  -sum(dKda.*Ue2,1)';                                           % Compliance sensitivities wrt. a
    % ----------- INDICATOR FIELD VOLUME PENALTY FUNCTION AND SENSETIVITIES
    S = mean(sPhys(:,3));                                                  % Weighted indicator field volume fraction (dilated design)
    dSds = ones(ne,1)./ne;                                                 % Sensitivities wrt. s
    % --------------- VOLUME FRACTION CONSTRAINT FUNCTION AND SENSETIVITIES
    f = squeeze(mean(rho));                                                % Domain volume fraction
    dfdw = (1-fliplr(wPhys(:,:,3)))./ne;                                   % Dilated volume fraction sensitivities wrt. w
    if mod(loop,20)==0; dilatedVolfrac=(f(3)/f(2))*volFrac; end            % Update dilated volume fraction wrt. nominal volume fraction
    % -------------------------------- CHAIN AND PRODUCT RULE SENSITIVITIES
    [dmudw(pasE,:,:),dfdw(pasE,:,:),dSds(pasE,:)] = deal(0.0);             % Correct for passive elements
    dJdw = squeeze(sum(dJdmu .* dmudw,2));                                 % mu to w sensitivities
    [dJds,dfds] = deal(sum(dJdw.*wTilde,2), sum(dfdw.*wTilde,2));          % Indicator field sensitivities wrt. indicator field
    [dJdw,dfdw] = deal(dJdw.*sPhys(:,:,1), dfdw.*sPhys(:,:,3));            %    and widths
    dHds = beta.*(1-tanh(beta*(sTilde(:)-eta([1,3,3]))).^2)./(tanh(beta*...% Robust projection of sensitivities wrt. indicator field
        eta([1,3,3]))+tanh(beta*(1-eta([1,3,3])))) .* [dJds,dfds,dSds];
    dFdw = convn(reshape([dJdw dfdw],[nelY,nelX,4])./HsW,hW,'same');       % Filter sensitivities wrt. widths
    dFds = convn(reshape(dHds,[nelY,nelX,3])./HsS,hS,'same');              % Filter sensitivities wrt. indicator field
    [dJdw,dfdw] = deal(dFdw(:,:,[1,2]),dFdw(:,:,[3,4]));                   % Extract filtered sensitivities
    [dJds,dfds,dSds] = deal(dFds(:,:,1),dFds(:,:,2),dFds(:,:,3));
    % ---------------------- OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES
    if loop == 1; xOld0 = x; xOld1 = xOld0; as = []; end                   % Initialise asymtotes and historic design variables
    [obj,g] = deal(J*cScale + S*sWeight, f(3)/dilatedVolfrac-1.0);         % Set and scale objective and dilated volume fraction constraint
    dObjdx =[dJdw(:)*cScale;dJda(:)*cScale;dJds(:)*cScale+dSds(:)*sWeight];% Set and scale sensitivities
    dgdx = [dfdw(:); zeros(ne,1); dfds(:)]./dilatedVolfrac;                % Set constraint sensitivities
    x=[w(:); a(:); s(:)];                                                  % Update design vector
    [xmin,xmax] = deal(max(xmin0,x-xmove),min(xmax0,x+xmove));             %    and bounds
    [x0(dv),as] = ocUpdate(loop,x(dv),xmax(dv),xmin(dv),dObjdx(dv),...     % Call OC. (Only active design variables (dv) are considered)
        g,dgdx(dv),[0.7 1.2],xOld0(dv),xOld1(dv),as,beta,restartAs);
    xOld1 = xOld0; xOld0 = x; x = x0;                                      % Update historic design variables
    [w(:), a(:), s(:)] = deal(x(1:2*ne), x(1+2*ne:3*ne), x(1+3*ne:4*ne));  % Distribute updated design variables
    change = max(abs(x(:)-xOld0(:)));                                      % Compute design change
    % ----------------------------------------- ON-THE-FLY DEHOMOGENISATION
    N = cat(3,[cos(a+pi/2),sin(a+pi/2)],[cos(a),sin(a)]);                  % Layer direction vectors
    if mod(loop,deHomFrq) == 0; tp = tic;
        rhoPhys = phasorDehomogenise(deHomGrid,wMin,wPhys(:,:,2),N,...     % Call phasor-based dehomogenisation
            alignItr,align);
        rhoPhys(pasRhoPhys) = 1.0;                                         % Correct passive solid elements
    else; tp = tic;
    end; phT = toc(tp);
    % ---------------------------------------------- PRINT RESULTS AND PLOT
    fprintf([' Itr: %3i Obj: %.4f J: %.3f S: %.3f Vol: %.3f (ph: %.3f)',...% Print stats
        ' ch: %.3f Time: %.3f (ph: %.3f)\n'],loop,obj,J,S,...
        f(2),mean(rhoPhys(:)),change,toc(T),phT);
    Tp = [1,-1,-1,1].*reshape(N,[],4).*wPhys(:,[2,2,1,1],2).^0.75; % Scale tangents for plotting
    h(1).CData(:) = rho(:,:,2);                                            % Update Rank-2 density plot
    for i=1:2; h(i+1).UData(:)= Tp(:,i*2-1); h(i+1).VData(:) =Tp(:,i*2); end % Update quiver plots of Rank-2 tangents
    if mod(loop,deHomFrq) == 0 && deHomFrq>0; h(4).CData = rhoPhys; end    % Update dehomgenisation density
    drawnow; Ttot = Ttot + toc(T);                                         % Accumulate time
end
% --------------- POST EVALUATE MULTI-SCALE RESULT FROM INTERMEDIATE DESIGN
if loopMax > 0
    wPhys = wTilde.*sPhys(:,:,2);                                          % Correct intermediate thickness field if not provided
    rho = 1 - prod(1 - wPhys,2);                                           % Get relative denisty and
    muPhys = getMuFromW(rho,wPhys);                                        %    compute mu from w from intermediate design
    C = getConstitutiveRank2(a,muPhys,E,nu,EMin);                          % Get elastic properties and sensitivities wrt. mu1, mu2, and a
    U = doFEA(F,U,C,iK,jK,freeDofs,KE0,Kp);                                % Call FE-analysis to get displacements
    [J,f] = deal(F'*U, mean(rho));                                         % Global compliance and volume fraction
    fprintf(['Multi-scale structure, intermediate design: J: %.3f ',...    % Print result
        'Vol: %.3f Total time: %.3f \n\n'],J,f,Ttot+toc(T));
    TO = struct('w',wPhys,'N',N,'f',f,'J',J); figure; Ttot=Ttot+toc(T);    % Store output, create new dehomogenisation figure and accumulate time
end
%% ================================================== POST DEHOMOGENISATION
T = tic; t=tiledlayout('flow',"TileSpacing","tight","Padding","tight");    % Create tiled layout for dehomogenisation plots
fprintf(['Dehomogenisation to single-scale structure, ', ...
    'with minimum length-scale: %.3f...\n'],dMin); nexttile(t);
rhoPhys = phasorDehomogenise(deHomGrid,wMin,wPhys,N,alignItr,align);       % Call phasor-based dehomogenisation
rhoPhys(pasRhoPhys) = 1.0;                                                 % Correct passive solid elements
fp = mean(rhoPhys(:));                                                     % Compute volume fraction and
fErr = fp/f;                                                               %    volume fraction error
imagesc(1-rhoPhys); colormap(gray); axis off; daspect([1,1,1]);clim([0,1]);% Plot dehomogenisation density field
fprintf(' Vol: %.3f err: %.2f%% Time: %.3f (total: %.3f)\n',...            % Print results
    fp,(fErr-1)*100,toc(T)+[0,Ttot]); drawnow; Ttot=Ttot+toc(T);           % Accumulate time
% ----------------------------------------------- DEHOMOGENISATION ANALYSIS
if eval
    fprintf('\nAnalysing dehomogenisation result...\n'); T=tic;
    [~,edofMat,iK,jK,U,F,freeDofs,Kp] = prepFEA(...                        % Get edofmat, triplets and free dofs for dehomogenised design
        deHomGrid.fsize(2),deHomGrid.fsize(1),deHomGrid.fscale);
    C = getIsotropicElasticity(E,EMin,nu,rhoPhys);                         % Get isotropic elastic properties
    U = doFEA(F,U,C,iK,jK,freeDofs,KE0,Kp);                                % Call FE-analysis to get displacements
    Jp = F'*U;                                                             %    and global compliance of dehomogenised design
    jErr = Jp/J;                                                           % Compliance error
    strain = reshape(BMatrixMid*U(edofMat'),3,1,[]);                       % Compute element strian
    stress = pagemtimes(reshape(permute(C(:,:,[1,2,3,2,4,5,3,5,6]),...     % Compute element stress
        [3,1,2]),3,3,[]),strain);
    W0 = 0.5*squeeze(sum(stress.*strain,1));                               % Compute strain-enegy-density
    W0 = real(log10(W0)); W0(rhoPhys(:) == 0) = -inf; ax = nexttile(t);    % Scale strain-enegy-density
    imagesc(reshape(W0,deHomGrid.fsize(1:2))); clim(max(W0)-[5,0]);        % Plot strain-enegy-density
    c=colorbar; c.Label.String='log10(W)'; axis off; daspect([1,1,1]);
    colormap(ax,[1,1,1;turbo(numel(unique(W0)))]); drawnow;
    fprintf([' J: %.3f err: %.2f%% wt err: %.2f%% Time: %.3f (total: ',... % Print results
        ' %.3f)\n'],Jp,(jErr-1)*100,(jErr*fErr-1)*100,toc(T)+[0,Ttot]);
    fprintf([' Evalulated on %dx%d grid (x%d scaled), results subjecte',...
        'd to h-conv. effects\n'],deHomGrid.fsize([2,1,2])./[1,1,nelX]);
end
end
%% ================================ FUNCTION TO PREPARE FEA OF BRIDGE MODEL
function [passiveElms,edofMat,iK,jK,U,F,freeDofs,Kp] =prepFEA(nelX,nelY,sc)
nodeNumbers = reshape(1:(1+nelX)*(1+nelY),1+nelY,1+nelX);                  % Get edofmat, stiffness matrix triplets, and initialise U
edofMat = repmat(reshape(2*nodeNumbers(1:end-1,1:end-1)+1,[],1),1,8)+...
    repmat([0,1,2*nelY+[2,3,0,1],-2,-1],nelX*nelY,1);
iK = reshape(kron(edofMat,ones(8,1))',[],1);
jK = reshape(kron(edofMat,ones(1,8))',[],1);
numDof = 2*(nelY+1)*(nelX+1); allDofs = 1:numDof; U = zeros(numDof,1);
% ----------------------------------------------- DEFINE LOADS AND SUPPORTS
[px,py] = deal(ceil(nelX/(15*sc))*sc,ceil(nelY/(30*sc))*sc);               % Passive block sizes (coarse scaled)
forceElmsB1 = floor(nelY*nelX*0.5) - (0:(ceil(0.5*px)-1))*nelY;            % Elements subjected to traction
forceElmsB2 = nelX*nelY - forceElmsB1 + nelY;
forceElms = union(forceElmsB1,forceElmsB2);
iF = edofMat(forceElms,:);                                                 % Nodes subjected to traction
sF = repmat(-1*[0 1 0 1 0 0 0 0],[numel(forceElms),1]);                    % Integrate traction over element surface
F = sparse(iF(:),ones(size(iF(:))),sF(:),numDof,1); F = F/sum(abs(F(:)));  % Build RHS force vector with unit magnitude
fixedDofsB1 = 2*px*(nelY+1) + (2*(nelY+1):2*(nelY+1):2*(nelY+1)*(px+1));   % Distributed Y-dofs at first support (BC1)
fixedDofsB2 = 2*(nelY+1)*(nelX + 1) - fixedDofsB1 + 2*(nelY+1);            % Distributed Y-dofs at second support (BC2)
[iB1,jB1] = meshgrid(fixedDofsB1); [iB2,jB2] = meshgrid(fixedDofsB2);      % Distributed Y-roller support (multipoint constraint)
Kp = sparse([iB1(:);iB2(:)],[jB1(:);jB2(:)],1e4,numDof,numDof);            %    with penalty method (penalty value = 1e4)
fixedDofs = fixedDofsB1(1)-1;                                              % Fix single X-dof
freeDofs = setdiff(allDofs,fixedDofs);                                     % Get free dofs
% ------------------- SET PASSIVE ELEMENTS AT LOAD AND ROLLER SUPPORT POINT
pasBC1 = px*nelY + (nelY:nelY:nelY*px); pasBC2 = nelX*nelY - pasBC1 + nelY;% BC1 and BC2 passive elements in X direction
passiveElms = reshape([forceElms, pasBC1, pasBC2] - (0:(py-1))',[],1);     % Grid passive elements from BC and traction
end
%% ================ FUNCTION TO TRANSFORM SINGLE-SCALE WIDTH TO MULTI-SCALE
function [mu,dmudw] = getMuFromW(rho,w)
[mu,dmudw] = deal(w, repmat(w,[1,1,2]));                                   % Initialise
% ------------------------------------------------------------------ GET MU
wSum = w(:,1)+w(:,2);                                                      % Sum of w
mu(:,1) = w(:,1)./wSum.*rho;                                               % Compute mu1 from w1 and w2
mu(:,2) = w(:,2).*rho./(wSum.*(1 - mu(:,1)));                              % Compute mu2 from mu1, w1 and w2
% ------------------------------------------------ GET SENSITIVITIES WRT. W
dmudw(:,1,1)=w(:,1).*((1-w(:,2)).*wSum-rho)./(wSum.^2)+rho./wSum;          % Compute dmu1dw1, sensitivities of mu1 wrt. w1
dmudw(:,1,2)=w(:,1).*((1-w(:,1)).*wSum-rho)./(wSum.^2);                    % Compute dmu1dw2, sensitivities of mu1 wrt. w2
dmudw(:,2,1)=-(((-1+w(:,2)).*wSum+rho).*(1-mu(:,1))-rho.*dmudw(:,1,1).*... % Compute dmu2dw1, sensitivities of mu2 wrt. w1
    wSum).*w(:,2)./(wSum.^ 2.*(1-mu(:,1)).^2);
dmudw(:,2,2)=-(((-1+w(:,1)).*wSum+rho).*(1-mu(:,1))-rho.*dmudw(:,1,2).*... % Compute dmu2dw2, sensitivities of mu2 wrt. w2
    wSum).*w(:,2)./(wSum.^ 2.*(1-mu(:,1)).^2)+rho./(wSum.*(1-mu(:,1)));
[mu(wSum < realmin,:), dmudw(wSum < realmin,:,:)] = deal(0);               % Correct zero width for numerical stability
end
%% ============== FUNCTION TO BUILD RANK-2 MATERIAL MODEL AND SENSITIVITIES
function [CT,dCTdmu1,dCTdmu2,dCTda] = getConstitutiveRank2(a,mu,E,nu,EMin)
a = reshape(a,1,1,[]);                                                     % Initialise and reshape to 3D array
[mu1,mu2] = deal(reshape(mu(:,1),1,1,[]), reshape(mu(:,2),1,1,[]));
[c,s,s2] = deal(cos(a), sin(a), sin(2*a));
[C,dCdmu1,dCdmu2] = deal(zeros(3,3,numel(a)));
T2D = @(A) reshape(permute(A,[3,1,2]),[],9);                               % Function to reshape result 3D to 2D arrays
TRI = @(A) reshape(A(:,[1,2,3,5,6,9]),1,[],6);                             % Function to extract lower triangular matrix components in 3D array
% ------------------------------------------- GET PROPERTIES IN LOCAL FRAME
denominator = (1 - mu2+mu1.*mu2*(1-nu^2));                                 % Orthogonal rank-2 constitutive matrix
C(1,1,:) = E./denominator.*(mu1)+EMin/(1-nu^2);                            % C11
[C(1,2,:),C(2,1,:)] = deal(E./denominator.*(mu2.*mu1*nu)+nu*EMin/(1-nu^2));% C12 and C21
C(2,2,:) = E./denominator.*(mu2.*(1-mu2+mu2.*mu1))+EMin/(1-nu^2);          % C22
C(3,3,:) = EMin/(1-nu^2)*(1-nu)/2;                                         % C33
% ----------------------------------------------- TRANSFORM TO GLOBAL FRAME
T = [c.^2 s.^2 c.*s; s.^2 c.^2 -c.*s; -2*c.*s 2*c.*s c.^2-s.^2];           % 2D rotation-transformation matrix
Tt = pagetranspose(T);                                                     % 2D-wise transpose of T
CT = TRI(T2D(pagemtimes(pagemtimes(Tt,C),T)));                             % Transformed constitutive matrix, Ct = T'*C*T;
% ---------------------------------------- GET SENSITIVITIES IN LOCAL FRAME
dtmp=E./(denominator.^2);                                                  % Intermediate sensitivity term
dCdmu1(1,1,:) = dtmp.*(1-mu2);                                             % Compute dC11dmu1, sensitivity of C11 wrt. mu1
[dCdmu1(1,2,:),dCdmu1(2,1,:)] = deal(dtmp.*mu2*nu.*(1-mu2));               % Compute dC12dmu1 and dC21dmu1, sensitivity of C12 and C21 wrt. mu1
dCdmu1(2,2,:) = dtmp.*(mu2.^2*nu^2.*(1-mu2));                              % Compute dC22dmu1, sensitivity of C22 wrt. mu1
dCdmu2(1,1,:) = dtmp.*(mu1.*(1-mu1*(1-nu^2)));                             % Compute dC11dmu2, sensitivity of C11 wrt. mu2
[dCdmu2(1,2,:),dCdmu2(2,1,:)] = deal(dtmp.*mu1*nu);                        % Compute dC12dmu2 and dC21dmu2, sensitivity of C12 and C21 wrt. mu2
dCdmu2(2,2,:) = dtmp.*(1-2*mu2+2*mu1.*mu2+(mu2).^2-2*(mu2).^2.*mu1+...     % Compute dC22dmu2, sensitivity of C22 wrt. mu2
    (mu1).^2.*(mu2).^2+(mu2).^2*nu^2.*mu1-(mu2).^2.*(mu1).^2*nu^2);
% ----------------------------------------------- TRANSFORM TO GLOBAL FRAME
dT =  [-s2 s2 2*c.^2-1; s2 -s2 -2*c.^2+1; -4*c.^2+2 4*c.^2-2 -2*s2];       % Sensitivities of transformation matrix wrt. a
dCTdmu1 = TRI(T2D(pagemtimes(pagemtimes(Tt,dCdmu1),T)));                   % Transformed constitutive sensitivity matrix, dCtdmu1 = T'*dCdmu1*T
dCTdmu2 = TRI(T2D(pagemtimes(pagemtimes(Tt,dCdmu2),T)));                   % Transformed constitutive sensitivity matrix, dCtdmu2 = T'*dCdmu2*T
dCTda =  TRI(T2D(pagemtimes(pagemtimes(dT,'transpose',C,'none'),T) ...     % Transformed constitutive sensitivity matrix,
                 + pagemtimes(pagemtimes(Tt,C),dT)));                      %    dCtda = dTda'*C*T + T'*C*dTda
end
%% ======================= FUNCTION TO BUILD SOLID ISOTROPIC MATERIAL MODEL
function C = getIsotropicElasticity(E,EMin,nu,rho)
C = zeros(1,numel(rho),6);                                                 % Initialise SIMP constitutive matrix, only unique components considered
ESimp = EMin + (E-EMin)*rho(:);                                            % SIMP material model
[C(:,:,1),C(:,:,4)] = deal(ESimp/(1-nu^2));                                % C11 and C22
[C(:,:,2),C(:,:,6)] = deal(nu*ESimp/(1-nu^2), ESimp./(2*nu + 2));          % C12 and C33
end
%% ===================== FUNCTION TO PERFORM LINEAR FINITE ELEMENT ANALYSIS
function U = doFEA(F,U,C,iK,jK,freeDofs,KE0,Kp)
K = sparse(iK,jK,reshape(sum(pagemtimes(KE0,C),3),[],1));                  % Create stiffness matrix from triplet values, sK = sum(KE0i*Ci,i=1..6)
K = K + Kp*max(diag(K));                                                   % Apply scaled penalty term
U(freeDofs)=decomposition(K(freeDofs,freeDofs),'chol','lower')\F(freeDofs);% Direct solve linear systems
end
%% ========== FUNCTION TO BUILD STRUCT OF SIZES, PARAMETERS, AND INDEX SETS
function grid = prepPhasor(nelX,nelY,dMin,wmin)
% ------------------------------------------- SCALING OF INTERMEDIATE GRIDS
Lx = 1; Ly = nelY/nelX; if nelY > nelX; Ly = 1; Lx = Ly*nelX/nelY; end     % Map domain to bounding box w. max edge length 1
[h_c,grid.h_c] = deal(Lx/nelX);                                            % Coarse grid element edge-length within the mapped domain
grid.omega = wmin/(dMin*h_c);                                              % Set phasor fixed signal frequency from minimum feature size (dMin)
grid.fpx = 3;                                                              % Pixel mapping feature size
[sc,grid.scale] = deal(ceil(grid.omega*Lx/(0.1*nelX)));                    % Minimal upscaling factor to capture details when sampling
[fsc,grid.ifscale] = deal(max(ceil(2*grid.omega*Lx/(wmin*nelX)),sc));      % Minimal upscaling factor to capture details for branch closure
mup = ceil((grid.omega*grid.fpx/wmin)/(fsc*max(nelX,nelY)));               % Upscaling factor to ensure sufficient resolution after thresholding
[grid.fscale,grid.Lx,grid.Ly]=deal(fsc*mup,Lx,Ly);                         % Store computed sizes in struct
grid.csize = [nelY, nelX, nelY*nelX];                                      % Coarse grid sizes
grid.isize = [sc*grid.csize(1:2), sc^2*grid.csize(3)];                     % First intermediate grid sizes
grid.ifsize = [fsc*grid.csize(1:2), fsc^2*grid.csize(3)];                  % Second intermediate grid sizes
grid.fsize = [mup*grid.ifsize(1:2), mup^2*grid.ifsize(3)];                 % Final grid sizes
% -------------------------------------------------- FILTERS AND INDEX SETS
kx = [0.0585,0,-0.0585;0.0965,0,-0.0965;0.0585,0,-0.0585];                 % 2-D Kernel for Gaussian derivative filter
[grid.kx,grid.ky,grid.kxy,grid.ksum]=deal(kx,kx',kx*kx',sum(kx*kx','all'));% Store derivative kernels
grid.filter1 = setUpIndexes(1,grid.csize,1,1);                             % 3x3 filter-size index sets
grid.filter4 = setUpIndexes(4,grid.csize,1,1);                             % 6x6 filter-size index sets
% --------------------------- STRUCTURED GRIDS AT THE DIFFERENT RESOLUTIONS
[h_i,h_if,h_f] = deal(h_c/sc,h_c/fsc,h_c/fsc/mup);                         % Element sizes on first and second intermediate, and final grid
[x_c,y_c] = deal(h_c/2:h_c:Lx-h_c/2, h_c/2:h_c:Ly-h_c/2);                  % Phasor kernel coordinates on coarse grid
[grid.x_i, grid.y_i]  = deal(h_i/2:h_i:Lx-h_i/2, h_i/2:h_i:Ly-h_i/2);      % Sampling coordinates on first intermediate grid
[grid.x_if,grid.y_if] = deal(h_if/2:h_if:Lx-h_if/2, h_if/2:h_if:Ly-h_if/2);% Sampling coordinates on second intermediate grid
[grid.y_f, grid.x_f]  = deal(h_f/2:h_f:Ly-h_f/2, h_f/2:h_f:Lx-h_f/2);      % Coordinates on final resolution grid
[grid.Y_c, grid.X_c]  = ndgrid(y_c,x_c);                                   % Coarse grid
[grid.Y_i, grid.X_i]  = ndgrid(grid.y_i,grid.x_i);                         % First intermediate grid
[grid.Y_if,grid.X_if] = ndgrid(grid.y_if,grid.x_if);                       % Second intermediate grid
% -------------------------------------------------- INITIALISE INTERPOLANT
grid.intp_c = griddedInterpolant(grid.Y_c,grid.X_c,grid.Y_c*0,'linear');   % Initialise interpolant on coarse grid
grid.intp_i = griddedInterpolant(grid.Y_i,grid.X_i,grid.Y_i*0,'linear');   % Initialise interpolant on first intermediate grid
grid.intp_if = griddedInterpolant(...                                      % Initialise interpolant on second intermediate grid
    {grid.y_if',grid.x_if'},grid.Y_if*0,'linear');
% ---------------------------------------------- PHASOR PARAMETER SELECTION
beta = -log(1e-1)/(2*grid.h_c^2);                                          % Kernel and sampling filter bandwidth
grid.R = max((1.1/grid.omega)^2,4*grid.h_c^2);                             % Kernel neighbourhood radius
grid.omega_dsc = 1/(h_i*grid.omega);                                       % Discretised wave-length on first intermediate grid
grid.sigma = 1/(2*grid.omega^2);                                             % Gaussian wave-length radius squared
[grid.rx, grid.ry, grid.rrx, grid.rry] = deal(0.5,2,0.5,2);                % Alignment and Sampling anisotropy
grid.cutoff = -log(1e-3)/beta;                                             % Signal response distance cut-off
grid.bpainv = 1/(beta+0.9*beta); grid.abpainv = 0.9*beta*grid.bpainv;      % Sampling filter intermediate computations of the combination
grid.babpainv = beta*grid.abpainv;                                         %    of the phasor kernel and sampling filter bandwidhts
% ------------------------------ INDEX SETS FOR BRANCH CONNECTION PROCEDURE
pd = ceil(grid.omega_dsc/4)-1;                                             % Pad-size for bifurcation point search
grid.bsearch = setUpIndexes(pd,grid.isize,0,1);                            % Index-sets for bifurcation point search based on
[dy,dx] = ndgrid((-pd:pd).^2/(pd+1).^2); grid.bsearch.circle = dx+dy <= 1; %    1/4-wavelength circular filter kernel
wsl = grid.omega_dsc*fsc/sc;                                               % Half discretised wavelength (wsl)
pd = ceil(wsl/2)-1;                                                        % Pad-size from wsl for disconnections
grid.discomp = setUpIndexes(pd,grid.ifsize,0,0);                           % Index-sets for degree of disconnection computations
[dy,dx] = ndgrid((-pd:pd).^2/(pd+1).^2); grid.discomp.circle = dx+dy <= 1; %    corresponding to a circular neighbourhood with the radius of wsl
grid.clspatch = setUpIndexes(ceil(wsl*1.5),grid.ifsize,0,1);               % Index-sets for patch-wise localised branch closure
end
%% ========================== FUNCTION TO BUILD STRUCT OF FILTER INDEX SETS
function subgrid = setUpIndexes(pd,nc,widmat,wnopdMat)
idpad = reshape(1:(nc(2)+2*pd)*(nc(1)+2*pd),nc(1)+2*pd,nc(2)+2*pd);        % Enumerate padded matrix indexes with pad-size pd
idofVec = reshape(idpad(1:end-2*pd,1:end-2*pd),nc(3),1);                   % Index of first element in each neighbourhood
if widmat
    subgrid.idofMat = repmat(idofVec,1,(2*pd+1)^2)+...                     % Build matrix of row-wise index-sets
        repmat(reshape((0:2*pd)'+((0:2*pd)*(nc(1)+pd*2)),1,[]),nc(3),1);   %    in padded array for each element in
end
if wnopdMat
    subgrid.nopdMat = zeros(nc(1)+2*pd,nc(2)+2*pd);                        % Adjust enumeration to original matrix indexes
    subgrid.nopdMat(pd+1:end-pd,pd+1:end-pd) = reshape(1:nc(3),nc(1:2));
end
[subgrid.pd, subgrid.idpad, subgrid.idofVec] = deal(pd, idpad, idofVec);   % Store computed sets to output struct
end
%% ============== FUNCTION TO PRECOMPUTE KERNEL NEIGHBORHOODS AND DISTANCES
function [align]=prepPhasorKernels(grid)
sample(grid.csize(3)) = struct('lnrs',[],'ldist',[],'xnrs',[],'xdist',[]); % Initialise structure to store precomputations on kernel basis
align.ids = setUpIndexes(ceil(1/(min(grid.rrx,grid.rry))),grid.csize,1,1); % Construct index set for neighbour search radius
[icy,icx] = meshgrid(1:grid.csize(2),1:grid.csize(1));                     % Index mesh-grid for signal impact candidate extraction
swx =ceil(sqrt(2*grid.cutoff)/2/min(grid.rx,grid.ry)*max(grid.isize(1:2)));% Cut-off wrt. discrete sampling point radius
% ------------------------- PERFORM PRECOMPUTATIONS PER COARSE-GRID ELEMENT
for j=1:grid.csize(3)
    ms = round(([icx(j) icy(j)]-0.5)*grid.scale);                          % Potentially influenced sampling points for kernel j
    idlx = max(ms(1)-swx,1):min(ms(1)+swx+1,grid.isize(1));                % Local patch of indexes x-direction
    idly = max(ms(2)-swx,1):min(ms(2)+swx+1,grid.isize(2));                % Local patch of indexes y-direction
    [sample(j).idXY,id] = deal(reshape((idly-1)*grid.isize(1)+idlx',[],1));% Index-set of full patch
    sample(j).xdist=[grid.X_i(id)-grid.X_c(j),grid.Y_i(id)-grid.Y_c(j)];   % Distance from kernel j to potential sampling points
end
[align.idN,align.sample]=deal(align.ids.nopdMat(align.ids.idofMat),sample);% Extract global indexes of potential kernel neighbours and save precomputed neighbourhoods
end
%% ======================= MAIN PROCEDURE FOR PHASOR-BASED DEHOMOGENISATION
function rho = phasorDehomogenise(grid,wMin,w,N,alignItr,align)
nLayer = size(w,2);                                                        % Number of lamination layers
N=normalize(N,2,'norm'); [Nx,Ny]=deal(squeeze(N(:,1,:)),squeeze(N(:,2,:)));% Kernel normalized orientations    
actKernels = (w >= wMin*0.95) & (w < 1-1e-3);                              % Active kernels per layer
% ---------------------------- INTERPOLATION AND SMOOTHING OF COARSE FIELDS
indicator = getIndicator(wMin,w,grid,nLayer);                              % Get layer indicator fields at relevant resolutions
[Nx,Ny] = filterVectorField(grid,nLayer,actKernels,Nx,Ny);                 % Filter orientation layers to reduce discontinuous jumps
[iNx, iNy, ifNx, ifNy, grid] = interpOrientations(grid,Nx,Ny,nLayer);      % Upscale lamination orientations to relevant resolutions
[w_i,w_if,w_f] = interpThickness(grid,nLayer,w,wMin,indicator.coarse);     % Upscale relative thickness' to the intermediate grid and fine grid
[domain,boundary,align] = addBoundary(grid,w_f,align,wMin,indicator,Nx,Ny);% Add phasor-boundary to smoothed structural domain
% ------------------------------------------------ PHASOR NOISE METHODOLOGY
align = phaseAlignment(grid,Nx,-Ny,alignItr,actKernels,align,nLayer,w);    % Align phase shifts of active kernels
non_solid = (w_i < 1-1e-3).*indicator.intermed;                            % Binary indicator field of intermediate thickness to reduce sampling space
phasor_i = phasorSampling(grid,Nx,-Ny,iNx,-iNy,actKernels,...              % Sample aligned phasor kernels on first intermediate grid
    non_solid>=0.01,align,nLayer);
[phasor_if,triangular_if] = deal(zeros(grid.ifsize(3),nLayer));            % Initialise phasor and triangular wave field
intp_i = grid.intp_i; intp_i.Method = 'cubic';                             % Upscale sampled complex phasor field to second intermediate grid
for r=1:nLayer                                                             %    from first intermediate grid using cubic interpolation
    intp_i.Values(:) = phasor_i(:,r);                                      % Interpolate layer-wise
    phasor_if(:,r) = reshape(intp_i({grid.y_if,grid.x_if}),[],1);
end
% ------------------------------------------------- STRUCTURAL IMPROVEMENTS
for r=1:nLayer                                                             % Branch connection on second intermediate grid w. inputs:
    inp1 = reshape(phasor_i(:,r),grid.isize(1:2));                         %    Sampled complex phasor field
    inp2 = reshape(non_solid(:,r)>=0.5,grid.isize(1:2));                   %    Region of interest is intermediate thicknesses
    inp3 = reshape(phasor_if(:,r),grid.ifsize(1:2));                       %    Upscaled complex wave
    inp4 = reshape(w_if(:,r),grid.ifsize(1:2));                            %    Second intermediate thickness
    inp5 = reshape(ifNx(:,r),grid.ifsize(1:2));                            %    Orientation x-component
    inp6 = reshape(ifNy(:,r),grid.ifsize(1:2));                            %    Orientation y-component
    triangular_if(:,r) = closeBranches(grid,inp1,inp2,inp3,inp4,inp5,inp6);% Obtain triangular wave field
end
% ------------------------------------------------- COMBINE FINAL STRUCTURE
rho = reshape(max(w_f,[],2) >= 0.99,grid.fsize(1:2));                      % Starting from only fully solid regions
for r=1:nLayer                                                             % Layerwise:
    grid.intp_if.Values(:) = (triangular_if(:,r)+1)/2;                     % Normalised triangular wave to interpolate
    triangular_if_f = grid.intp_if({grid.y_f',grid.x_f'});                 % Linear interpolation for upscaling of triangular wave
    thresholded =  triangular_if_f(:) >= 1-w_f(:,r);                       % Threshold wrt. thickness
    rho(:) = rho(:)+thresholded;                                           % Add layer-contribution to structure
end
rho = removeIslands2D(min(rho.*domain+boundary,1),0.1);                                         % Concentrate to smoothed domain and add boundary
end
%% ================== FUNCTION TO DETECT SINGULAR POINTS AND CLOSE BRANCHES
function triangular_if = closeBranches(grid,phasor_i,bregion,phasor_if,...
    w_if,ifNx,ifNy)
ws = grid.omega_dsc*grid.ifscale/grid.scale;                               % Relative periodicity
[gamma_x,gamma_y,numBranch] = locateBranchingPoints(phasor_i,bregion,grid);% Locate bifurcations
atan_if = atan2(imag(phasor_if),real(phasor_if));                          % Real-valued sawtooth wave-field from upscaled phasor field
sine_if = sin(atan_if);                                                    % Translations to phasor sine-wave
triangular_if = 2/pi*asin(sine_if);                                        %    and triangular field
if numBranch == 0; triangular_if=triangular_if(:); return; end             % No bifurcations found; return
% -------------------- COMPUTE DEGREE OF CONNECTION AT EACH BRANCHING POINT
gamma_x = max(min(gamma_x,grid.ifsize(1)),1);                              % Ensure mapped bifurcation points within bounding box
gamma_y = max(min(gamma_y,grid.ifsize(2)),1);
id_gamma = (gamma_y-1)*grid.ifsize(1)+gamma_x;                             % 1D index of branching points on intermediate mesh
circle = grid.discomp.circle; idofVec = grid.discomp.idofVec;              % Extract index-sets for degree of disconnection area
pd = grid.discomp.pd; sif_pad = zeros(grid.ifsize(1:2)+2*pd);              % Pad triangular field to allow for indexing
sif_pad(pd+1:end-pd,pd+1:end-pd) = triangular_if;
idofMat = repmat(idofVec(id_gamma),1,(2*pd+1)^2)+repmat(...                % Degree of connection (condeg) computed as average sinewave
reshape((0:2*pd)'+((0:2*pd)*(grid.ifsize(1)+pd*2)),1,[]),numBranch,1);     %    value within the circular filter kernel
condeg = (sum(sif_pad(idofMat).*circle(:)',2)/sum(circle(:))+1)/2;
% ------------------------- DETERMINE CLOSURE DIRECTIONS AND CONTROL POINTS
gamma_vx = ifNx((gamma_y-1)*grid.ifsize(1)+gamma_x);                       % Extract orientations at branching points for move directions
gamma_vy = ifNy((gamma_y-1)*grid.ifsize(1)+gamma_x);
[move1x,move1y] = deal(round(-ws/3*gamma_vy), round(ws/3*gamma_vx));       % Move magnitude and direction for determining closure-direction
move2 = (1-condeg)*ws/3;                                                   % Move magnitude to new branch centre weighted by degree of disconnection
[move2x,move2y] = deal(round(move2.*gamma_vy), round(-move2.*gamma_vx));   % Move direction and magnitude from original branching point to new centre
[move3x,move3y] = deal(2*move2x, 2*move2y);                                % Move magnitude to get control points on the other side of each branching point
gamma_xp = min(max(gamma_x+move1x,1),grid.ifsize(1));                      % Coordinates representing first candidate move direction
gamma_yp = min(max(gamma_y+move1y,1),grid.ifsize(2));
gamma_xm = min(max(gamma_x-move1x,1),grid.ifsize(1));                      % Coordinates representing second candidate move direction
gamma_ym = min(max(gamma_y-move1y,1),grid.ifsize(2));
sa = sign(triangular_if((gamma_ym-1)*grid.ifsize(1)+gamma_xm)-...          % Choosing move direction based on largest density
    triangular_if((gamma_yp-1)*grid.ifsize(1)+gamma_xp));
gamma_xc = min(max(gamma_x+sa.*move2x,1),grid.ifsize(1));                  % Perform move to obtain new centre in chosen direction
gamma_yc = min(max(gamma_y+sa.*move2y,1),grid.ifsize(2));
gamma_xo = min(max(gamma_x+sa.*move3x,1),grid.ifsize(1));                  % Perform double move to obtain opposing control point in chosen direction
gamma_yo = min(max(gamma_y+sa.*move3y,1),grid.ifsize(2));
% ---------------- COMPUTE ORIENTED DISTANCES LOCALISED NEAR BRANCH CENTRES
idB = (gamma_yc-1)*grid.ifsize(1)+gamma_xc;                                % Branch point centre linearised indexes
pd = grid.clspatch.pd; idofVec = grid.clspatch.idofVec(idB);               % Extract precomputed filter indexes only for branching points
idofMat = repmat(idofVec,1,(2*pd+1)^2)+repmat(reshape((0:2*pd)'+...        % Create index-set for neighbourhood about branching point
    ((0:2*pd)*(grid.ifsize(1)+pd*2)),1,[]),numBranch,1);
tmpd = zeros(grid.ifsize(1:2)+2*pd);                                       % Temporary placeholder for zero-padded arrays
[pdy,pdx] = deal(pd+1:grid.ifsize(1)+pd, pd+1:grid.ifsize(2)+pd);          % Indexes of interior of padded array
tmpd(pdy,pdx) = grid.X_if; xdist = tmpd(idofMat)-grid.X_if(idB);           % Pad x-coordinate grid and compute localised x-distances
tmpd(pdy,pdx) = grid.Y_if; ydist = tmpd(idofMat)-grid.Y_if(idB);           % Pad y-coordinate grid and compute localised y-distances
tmpd(pdy,pdx) = ifNx; gamma_vx = tmpd(idofMat);                            % Pad orientation x-component and extract localised values
tmpd(pdy,pdx) = ifNy; gamma_vy = tmpd(idofMat);                            % Pad orientation y-component and extract localised values
xdirdist = -xdist.*gamma_vy-ydist.*gamma_vx;                               % Map distances based on orientations to prepare anisotropy
ydirdist = xdist.*gamma_vx-ydist.*gamma_vy;
% ----------------- PERFORM LOCALLY WEIGHTED PHASE SHIFTS TO BRANCH CLOSURE
wx = 1/(2*pi); sqsigma = sqrt(grid.sigma/2); sigma_inv = 1/grid.sigma;     % Precompute Gaussian weights
idofMat(:) = grid.clspatch.nopdMat(idofMat(:));                            % Extract non-padded indexes
for K=1:numBranch                                                          % For each branching point:
    idl = idofMat(K,:)>0; idL = idofMat(K,idl)';                           %    Extract non-padded indexes of neighbours
    xd = xdirdist(K,idl)'; yd = ydirdist(K,idl)';                          %    Extract oriented distances within patch
    Pihat = exp(-2*sigma_inv*(wx*xd.^2+yd.^2).*(1-0.5*triangular_if(idL)));%    Compute modified Gaussian
    Pihat2 = Pihat.*Pihat; Pihat = 3*Pihat2-2*Pihat2.*Pihat;               %    Apply smoothstep for Gaussian cutoff
    pihat = pi*Pihat.*(1-(2/pi*asin(sine_if(idL))+1)/2);                   %    Use Gaussian and local densities to weight phase shift
    sine_if(idL) = max(sine_if(idL),max(sin(atan_if(idL)-pihat),...        %    Perform weighted localised phase shift in both directions
        sin(atan_if(idL)+pihat)));                                         %    and extract union to ensure solidified branch
end
% -------------------- PINCH PROCEDURE TO REDUCE MATERIAL AND IMPROVE SHAPE
triangular_if = 2/pi*asin(sine_if);                                        % Triangular wave-transformation with branch closures
[kx,ky] = deal([1,-1;1,-1;1,-1],[1,1,1;-1,-1,-1]);                         % Derivative filter weight-matrices
[nlxo,nlyo,idofMat,idy,idx] = deal(0,0,[],[],[]);
for K=1:numBranch
    [ix,iy,ix2,iy2]=deal(gamma_xc(K),gamma_yc(K),gamma_xo(K),gamma_yo(K)); % Branch control points
    xlims = max(ix-ceil(ws*2.1)+1,1):min(ix+ceil(ws*2.1),grid.ifsize(1));  % Localised patch surrounding branch
    ylims = max(iy-ceil(ws*2.1)+1,1):min(iy+ceil(ws*2.1),grid.ifsize(2));
    [nlx,nly] = deal(length(ylims)-2,length(xlims)-2);                     % Patch-size assuming patch is padded with size 1
    if nlxo ~= nlx || nly ~= nlyo
        [ny,nx] = deal(nlx+2,nly+2); [idy,idx] = meshgrid(1:ny,1:nx);      % Create index-set used for pinch interpolation procedure
        nodenrs = reshape(1:(nlx+2)*(nly+2),nly+2,nlx+2);                  % Set-up of index sets for localised derivative filter
        idofVec = reshape(nodenrs(1:end-2,1:end-2),nlx*nly,1);
        idofMat = repmat(idofVec,1,9)+...
            repmat([0 1 2 nly+2+[0,1,2] 2*(nly+2)+[0,1,2]],nlx*nly,1);
        [nlxo,nlyo] = deal(nlx,nly);
    end
    gamma_vx = ifNx(xlims,ylims); gamma_vy = -ifNy(xlims,ylims);           % Local orientations in patch
    X_patch = grid.X_if(xlims,ylims); Y_patch = grid.Y_if(xlims,ylims);    % Global coordinates of local patch
    xdist = X_patch-grid.X_if(ix,iy); ydist = Y_patch-grid.Y_if(ix,iy);    % Distance to centre within patch
    xd = xdist.*gamma_vy-ydist.*gamma_vx;                                  % Oriented distances to branch centre
    yd = xdist.*gamma_vx+ydist.*gamma_vy;
    xdist2 =X_patch-grid.X_if(ix2,iy2); ydist2 =Y_patch-grid.Y_if(ix2,iy2);% Distance to control point within patch
    xd2 = xdist2.*gamma_vy-ydist2.*gamma_vx;                               % Oriented distances to opposing control point
    yd2 = xdist2.*gamma_vx+ydist2.*gamma_vy;
    vx = zeros(nly+2,nlx+2); vy = vx;                                      % Initialise pinch-vectors
    w_patch = w_if(xlims,ylims); c_weight = 1./((2-w_patch)*sqsigma);      % Local thicknesses and magnitude adjusting weight
    for k = 1:3                                                            % Incremental pinch procedure in 3 steps
        step = min((1-condeg(K))/2*(k-1),1);                               % Degree of connection dependent step-size
        ixx = round(step*ix2+(1-step)*ix)-min(xlims)+1;                    % Move from centre to control-point based on incremental step-size
        iyy = round(step*iy2+(1-step)*iy)-min(ylims)+1;
        xdist_k = xdist*(1-step)+xdist2*step;                              % Compute weighted average of distances between centre
        ydist_k = ydist*(1-step)+ydist2*step;                              %    and control point based on step-size
        xdx3 = xdist_k.*gamma_vy(ixx,iyy)-ydist_k.*gamma_vx(ixx,iyy);
        ydx3 = xdist_k.*gamma_vx(ixx,iyy)+ydist_k.*gamma_vy(ixx,iyy);
        Pi_pinch_k = (exp(-0.5*sigma_inv*(wx*xdx3.^2+ydx3.^2)));           % Form anisotropic gaussian about current step-point
        vx_k = sum(Pi_pinch_k(idofMat(:,[1:3,7:9])).*kx(:)',2);            % Compute derivaties by filtering
        vy_k = sum(Pi_pinch_k(idofMat(:,[1,3,4,6,7,9])).*ky(:)',2);
        vx(2:end-1,2:end-1) = reshape(vx_k,[nly,nlx]);
        vy(2:end-1,2:end-1) = reshape(vy_k,[nly,nlx]);
        vlen = max([abs(vx(:));abs(vy(:))]); vx = vx/vlen; vy = vy/vlen;   % Normalise derivate magnitudes
        step = step.^2;                                                    % Choose a smaller step-size and use to compute average distances
        xdm = (step).*xd2+(1-step).*xd; ydm = step.*yd2+(1-step).*yd;      %    between current pinch-centre and original branch centre
        Pi_local_k = exp((-sigma_inv*(1.5*wx*(1+(k-1)/2)*xdm.^2+ydm.^2)... % Modified Gaussian with impact radius reduced incrementally
            -c_weight.*abs(ydm)));
        pw = ws./3.*Pi_local_k.*(1-w_patch);                               % Compute pinch magnitude weights in pixel-count
        triangular_if(xlims,ylims)=pinchPatch(...                          % Perform localised pinch operation by linear interpolation
            triangular_if(xlims,ylims),pw.*vx,pw.*vy,idy,idx,ny,nx);
    end
end
triangular_if = triangular_if(:);
end
%% ============ FUNCTION TO LOCATE BRANCH-POINTS AT FIRST INTERMEDIATE GRID
function [gamma_x,gamma_y,numBranch] = locateBranchingPoints(...
    phasor_i,bregion,grid)
[gamma_x, gamma_y, numBranch] = deal([],[],0);                             % Initialise
if max(abs(phasor_i(:))) < 1e-6; return; end                               % No branch points found
[circle, pd] = deal(grid.bsearch.circle, grid.bsearch.pd);                 % Index sets to search for local maxima in singular field
singular = exp(-abs(phasor_i)/max(abs(phasor_i(:))));                      % Enhanced singularity field
singularp = zeros(grid.isize(1:2)+pd*2);                                   % Setup zero-padded
singularp(pd+1:end-pd,pd+1:end-pd)=singular;                               %    singularity field
[idofVec, nodenrsl] = deal(grid.bsearch.idofVec, grid.bsearch.nopdMat);
candidates = bregion.*singular > 0.95; ncsum = sum(candidates(:));         % Find candidate local minima for bifurcation points
idofMat = repmat(idofVec(candidates),1,(2*pd+1)^2)+...                     % Padded indexes within candidates
    repmat(reshape((0:2*pd)'+((0:2*pd)*(grid.isize(1)+pd*2)),1,[]),ncsum,1);
local_max = zeros(grid.isize(1:2));
local_max(candidates)=max(singularp(idofMat).*circle(:)',[],2) == ...      % Find candidates that are true local minima
    singularp(idofMat(:,pd*(2*pd+1)+pd+1));
idofMat(:) = nodenrsl(idofMat(:));                                         % Translate local minima indexes to global
idmax = idofMat(local_max(candidates)==1,pd*(2*pd+1)+pd+1);
numBranch = length(idmax);                                                 % Number of branching points found
gamma_y =floor(idmax/grid.isize(1)); gamma_x =idmax-grid.isize(1)*gamma_y; % Map found points to discrete location on second intermediate mesh
msc = grid.ifscale/grid.scale;
[gamma_x,gamma_y] =deal(floor((gamma_x-0.5)*msc),floor((gamma_y+0.5)*msc));
end
%% ========================== FUNCTION TO PERFORM AN INCREMENTAL PINCH-STEP
function patch = pinchPatch(patch,dby,dbx,ly,lx,ny,nx)                     % Function to update the pixel-wise values of 'old_patch'
rx = lx+(dbx); ry = ly+(dby); lx = floor(rx); ly = floor(ry);              %    based on the linearly interpolated value obtained by
xyzero = lx >= 1 & lx <= nx-1 & ly >= 1 & ly <= ny-1;                      %    moving 'dbx' in the x-direction and 'dby' in the
rx = rx(xyzero); lx = lx(xyzero); dlx = rx-lx; ux = lx+1; dux = ux-rx;     %    y-direction within the 'old_patch' image patch
ry = ry(xyzero); ly = ly(xyzero); dly = ry-ly; uy = ly+1; duy = uy-ry;
alxly = dux.*patch((ly-1)*nx+lx); alxuy = dux.*patch((uy-1)*nx+lx);
auxly = dlx.*patch((ly-1)*nx+ux); auxuy = dlx.*patch((uy-1)*nx+ux);
patch(xyzero) = duy.*(alxly+auxly)+dly.*(alxuy+auxuy);
end
%% ========================== FUNCTION TO CREATE VARYING THICKNESS BOUNDARY
function [indicator_f,shellwave,align] = addBoundary(grid,w,align,wmin,...
    indicator,Nx,Ny)
% -------------------- FILTER INDICATOR FIELS AND OBTAIN BOUNDARY GRADIENTS
idsfr = reshape(max(indicator.coarse,[],2),grid.csize(1:2));               % Approximate filtered indicator field w. replication padding
ids = idsfr > (1-wmin)*max(idsfr(:));                                      % Corresponding binary indicator field
[kx,ky,kxy,ksum,pd] = deal(grid.kx, grid.ky, grid.kxy, grid.ksum, 1);      % Apply Gaussian filter to binary indicator field w. zero padding
idofMat = grid.filter1.idofMat;
idsp = zeros(grid.csize(1:2)+2*pd); idsp(pd+1:end-pd,pd+1:end-pd) = ids;
idsf0 = reshape(sum(idsp(idofMat).*kxy(:)'/ksum,2),grid.csize(1:2));
idsp(:) = 0;idsp(pd+1:end-pd,pd+1:end-pd) = idsf0.*ids;                    % Apply Gaussian derivative filter to filtered field
vx = reshape(sum(idsp(idofMat).*kx(:)',2),grid.csize(1:2));                %    with binary indicator cut-off to maintain sharp
vy = reshape(sum(idsp(idofMat).*ky(:)',2),grid.csize(1:2));                %    corners and small holes in boundary orientation
[vdir,vmag] = deal(atan2(vy,vx), vx.^2+vy.^2);                             % Boundary orientation direction and magnitude
[potential,bvec] = deal(vmag > 1e-3, [cos(vdir(:)),sin(vdir(:))]);         % Indexes of potential boundary kernels and orientation vectors
% ------------------------------------- ALIGN SIMILAR BOUNDARY ORIENTATIONS
pd = 4; idofMat2 = grid.filter4.idofMat(potential(:),:);                   % Index sets for alignment neighbourhood
nodenrs2 = grid.filter4.nopdMat; nodenrs2(pd+1:end-pd,pd+1:end-pd) = ...
    nodenrs2(pd+1:end-pd,pd+1:end-pd).*potential;
idofMat2(:) = nodenrs2(idofMat2(:)); idbb = idofMat2(:,pd*(2*pd+1)+pd+1);
for j=1:length(idbb)
    v = idbb(j); lnrs = idofMat2(j,:); lnrs(lnrs<1)=[];                    % Extract non-padded neighbour indexes of vector v
    dotprod = (bvec(v,1).*bvec(lnrs,1)+bvec(v,2).*bvec(lnrs,2))';          % Dot-product between vector v and neighbours
    fct = 1-min(abs(dotprod));
    avg_angle = sum(abs(dotprod).*exp(1i*(vdir(lnrs)+(dotprod<0)*pi)));    % Weighted average of neighbouring angles
    avg_angle = fct*exp(1i*vdir(v))+(1-fct)*avg_angle;
    vdir(v) = atan2(imag(avg_angle),real(avg_angle));                      % Update boundary orientation angle
    [bvec(v,2),bvec(v,1)] = deal(sin(vdir(v)), cos(vdir(v)));              % Extract orientation vector
end
% --------------------------------- UPSCALE ORIENTATIONS IN BOUNDARY REGION
grid.intp_c.Values = vmag.*exp(1i*vdir);                                   % Linear interpolant of boundary orientations
vec_a = grid.intp_c({grid.y_i,grid.x_i});                                  % The complex angle upscaled to the first intermediate grid
[bregion,vec_a]  = deal(abs(vec_a) > 1e-4, atan2(imag(vec_a),real(vec_a)));% Boundary region indicator and real-valued angle for upscaled orientations
% ------------------------------------ SAMPLE BOUNDARY ON INTERMEDIATE MESH
idbnd = potential & vmag >= 0.25*max(vmag(:));                             % Active boundary kernels within boundary region
new_omg = min(1./(8*grid.h_c),grid.omega/2); omg_rat = new_omg/grid.omega; % Adjusted periodicity to ensure staircase artefact cut-off
[grid.omega, grid.cutoff] = deal(new_omg, 0.5/new_omg);
grid.bpainv = 2*grid.bpainv; grid.abpainv = 0.5*grid.abpainv;              % Sampling filter intermediate computations of the combination
grid.babpainv = 0.5*grid.babpainv;
atmp = align; atmp.pshift = pi*((1-idsfr(:))*0.5+1/3);                     % Set predetermined phase shift to push bounary into domain
bphasor = phasorSampling(grid,bvec(:,1),bvec(:,2),cos(vec_a(:)),...        % Sample boundary phasor field directly without phase alignment
    sin(vec_a(:)),idbnd(:),bregion(:),atmp,1);
% --------------------- PHASE ALIGNMENT PRE-COMPUTATIONS FROM BOUNDARY INFO
ddot = Nx.*bvec(:,1)+Ny.*bvec(:,2); adot = abs(ddot); nd = max(adot,[],2); % Degree of layer-wise orientations and boundary kernel alignments
[align.idbnd,align.bdist] = deal(repmat(idbnd(:),1,size(w,2)),0*adot);     % Initialise boundary information arrays
align.idbnd = align.idbnd & abs(nd-adot)==0 & adot>0.95*max(adot(:));      % Layer-wise indicator of boundary aligned kernels
for r = 1:size(w,2)                                                        % Compute distance to boundary aligned kernels
    align.bdist(:,r) = min([(grid.X_c(:)-grid.X_c(align.idbnd(:,r))').^2+...
        (grid.Y_c(:)-grid.Y_c(align.idbnd(:,r))').^2,0*idbnd(:)+1],[],2);
end
% -------------------- CONSTRUCT INDICATOR CUT-FIELD BASED ON BOUNDARY WAVE
indicator_i = max(min(max(indicator.intermed,[],2),1),0);                  % Approximate filtered indicator w. replication padding
bsawtooth = atan2(imag(bphasor),real(bphasor));                            % The real-valued conversion of the sampled phasor field
positive = (sin(bsawtooth)>1e-1).*sqrt(indicator_i.*(1-indicator_i));      % Positive region of phasor sine-wave to limit span of wave
cutfi = 10*positive.*sin(bsawtooth+pi/2)+(indicator_i);                    % Indicator cut-field along boundary
Inter_i = grid.intp_i; Inter_i.Values(:) = cutfi(:);                       % Linear interpolant of cut-field
indicator_f = Inter_i({grid.y_f,grid.x_f});                                % Fine resolution indicator cut field
outline = (indicator_f > 0.99);                                            % Restrain domain boundary to structural body%
outline(grid.fpx+1:end-grid.fpx,grid.fpx+1:end-grid.fpx) = 0;              % Set domain boundary layer to minimum feature size
indicator_f = indicator_f > 0.5;                                           % Smoothed binary fine resolution indicator field
% ------------------------------------- UPSCALE BOUNDARY WAVE AND THRESHOLD
Inter_i.Values(:) = bphasor(:);                                            % Linear interpolant of complex phasor field
bphasor_f = Inter_i({grid.y_f,grid.x_f});                                  % Interpolated fine resolution phasor field
shellwave = 2/pi*asin(sin(atan2(imag(bphasor_f),real(bphasor_f))));        % Triangular boundary wave conversion
wmaxs = max(w,[],2); shellth = 2*omg_rat*min(max(wmaxs(:),wmin),0.99);     % Scaled boundary thickness for thresholding
shellwave(:) = (indicator_f(:)).*max(shellwave(:),outline(:)) >= 1-shellth;% Threshold boundary with indicator cut-off and ensure minimum feature size at boundary
indicator_f = removeIslands2D(indicator_f-shellwave,0.1);                  % Remove potential cut-off artefacts
end
%% ============================ FUNCTION TO GET LAYER-WISE INDICATOR FIELDS
function indicator = getIndicator(wmin,w,grid,nLayer)
% ------------------------------------ FILTER COARSE LAYER INDICATOR FIELDS
indicator=struct('coarse',zeros(grid.csize(3),nLayer),...                  % Initialise struct to store indicator arrays
    'intermed',zeros(grid.isize(3),nLayer),'fine',[]);                     %    for coarse and intermediate grids
[pd,idofMat] = deal(grid.csize(1:2), grid.filter1.idofMat);                % Set up filter
for r = 1:nLayer
    idr = reshape(w(:,r),grid.csize(1:2)) >= wmin*0.95;                    % Binary coarse layer indicator
    idr1=idr([1,1:pd(1),pd(1)],[1,1:pd(2),pd(2)]);                         % Padded indicator
    [grid.intp_c.Values(:),indicator.coarse(:,r)] = ...                    % Filtered layer indicator
        deal(sum(idr1(idofMat)*1/9,2));
    indicator.intermed(:,r)=reshape(grid.intp_c({grid.y_i,grid.x_i}),[],1);% Interpolated intermediate filtered indicator field
end
end
%% ================================ FUNCTION TO UPSCALE KERNEL ORIENTATIONS
function [iNx,iNy,ifNx,ifNy,grid] = interpOrientations(grid,Nx,Ny,nLayer)
[iNx,iNy] = deal(zeros(grid.isize(3),nLayer));                             % Initialise arrays for intermediate resolution
[ifNx,ifNy] = deal(zeros(grid.ifsize(3),nLayer));                          %    orientation components
for r=1:nLayer
    grid.intp_c.Values(:) = Ny(:,r);                                       % Linear interpolation of n_2 orientation component
    iNy(:,r) = reshape(grid.intp_c({grid.y_i,grid.x_i}),[],1);             %    to first intermediate grid
    ifNy(:,r) = reshape(grid.intp_c({grid.y_if,grid.x_if}),[],1);          %    to second intermediate grid
    grid.intp_c.Values(:) = Nx(:,r);                                       % Linear interpolation of n_1 orientation component
    iNx(:,r) = reshape(grid.intp_c({grid.y_i,grid.x_i}),[],1);             %    to first intermediate grid
    ifNx(:,r) = reshape(grid.intp_c({grid.y_if,grid.x_if}),[],1);          %    to second intermediate grid
    e_veci = iNx(:,r).^2+iNy(:,r).^2 <= 0.975;                             % Locations of inconsistently interpolated orientations first intermediate grid
    e_vecif = ifNx(:,r).^2+ifNy(:,r).^2 <= 0.975;                          % Locations of inconsistently interpolated orientations second intermediate grid
    interp = grid.intp_c; interp.Method = 'nearest';                       % Correct inconsistent seamlines by nearest neighbour interpolation
    iNx(e_veci,r) = interp(grid.Y_i(e_veci),grid.X_i(e_veci));             %    for component n_1 on first intermediate grid
    ifNx(e_vecif,r) = interp(grid.Y_if(e_vecif),grid.X_if(e_vecif));       %    for component n_1 on second intermediate grid
    interp.Values(:) = Ny(:,r);
    iNy(e_veci,r) = interp(grid.Y_i(e_veci),grid.X_i(e_veci));             %    for component n_2 on first intermediate grid
    ifNy(e_vecif,r) = interp(grid.Y_if(e_vecif),grid.X_if(e_vecif));       %    for component n_2 on second intermediate grid
end
end
%% =========================== FUNCTION TO INTERPOLATE RELATIVE THICKNESSES
function [w_i,w_if,w_f] = interpThickness(grid,nLayer,w,wmin,indicator_c)
w_i = zeros(grid.isize(3),nLayer);                                         % Initialise arrays for layer-wise interpolated fields
w_if = zeros(grid.ifsize(3),nLayer); w_f = zeros(grid.fsize(3),nLayer);
for r = 1:nLayer
    grid.intp_c.Values(:) = (indicator_c(:,r)>=0.25).*max(w(:,r),wmin);    % Values to interpolate
    w_i(:,r) = reshape(grid.intp_c({grid.y_i,grid.x_i}),[],1);             % Linear interpolation for smoothly varying thickness to first intermediate grid
    w_if(:,r) = reshape(grid.intp_c({grid.y_if,grid.x_if}),[],1);          % Linear interpolation for smoothly varying thickness to second intermediate grid
    w_f(:,r) = reshape(grid.intp_c({grid.y_f,grid.x_f}),[],1);             % Linear interpolation for smoothly varying thickness at final resolution
end
end
%% ================ FUNCTION TO FILTER VECTOR FIELDS IN CASE OF LAYER-JUMPS
function [Nx,Ny] = filterVectorField(grid,nLayer,actKernels,Nx,Ny)
pd = grid.filter1.pd; idofMat = grid.filter1.idofMat;                      % Extract pre-computed filter kernels
[npy,npx] = deal(pd+1:grid.csize(1)+pd, pd+1:grid.csize(2)+pd);            % Interior indexes of padded array for filtering
Nxp = zeros(grid.csize(1:2)+2*pd); [Nyp,activep] = deal(Nxp);              % Initialise arrays to store padded arrays during ordering
active_nrm = zeros(grid.csize(3),1);                                       % Initialise array to store neighbourhood information
for r = 1:nLayer                                                           % Layer-wise ordering of orientations to allow for Rank-N applications
    Nxp(npy,npx) = reshape(Nx(:,r),grid.csize(1:2));                       % Zero-padded orientation fields
    Nyp(npy,npx) = reshape(Ny(:,r),grid.csize(1:2));
    fdot = mean(Nx(:,r).*Nxp(idofMat)+Ny(:,r).*Nyp(idofMat),2);            % Filtered dot-product of orientations wrt. to filter
    fixed = fdot > 0.975 & actKernels(:,r);                                % Fix kernels with minimal orientation variation within neighbourhood
    [it, ch, nfix] = deal(0, 1, sum(fixed));                               % Initialise convergence criteria
    while it < 50 && ch > 0
        ida = find(~fixed);                                                % Indexes of non-fixed kernels
        activep(npy,npx) = reshape(fixed(:),grid.csize(1:2));              % Padded array of fixed kernels in current layer
        active_nrm(ida) = max(sum(activep(idofMat(ida,:)),2),1);           % Normalising factor of number of fixed kernels within filter kernel radius
        for k = ida'                                                       % For each active kernel
            Nxp(npy,npx) = reshape(Nx(:,r),grid.csize(1:2));               % Zero-padded orientation fields
            Nyp(npy,npx) = reshape(Ny(:,r),grid.csize(1:2));
            [Nxp1,Nyp1] = deal(Nxp(idofMat(k,:)),Nyp(idofMat(k,:)));       % Mean of active orientation components wrt. filter
            cand=[Nx(k,r),Ny(k,r);-Nx(k,r),-Ny(k,r);...                    % Candidate orientations wrt. to rotation invariance
                -Ny(k,r),Nx(k,r);Ny(k,r),-Nx(k,r)];                        %    and layer-jumping
            fdot = sum(repmat(activep(idofMat(k,:)),4,1).*...              % Filtered dot-product of orientations wrt. to the candidates
                (cand(:,1).*Nxp1+cand(:,2).*Nyp1+1)/2,2)/active_nrm(k);
            upid = find(max(fdot)==fdot,1,'first');                        % Choose best fit candidate
            [Nx(k,r),Ny(k,r)] = deal(cand(upid,1),cand(upid,2));           % Assign chosen candidate
            fixed(k) = fdot(upid)>0.975;                                   % Fix kernel if minimal orientation variation within neighbourhood achieved
        end
        [it, ch] = deal(it+1,sum(fixed)-nfix);                             % If change recorded convergence is not reached
    end
end
end
%% ====================== FUNCTION TO PERFORM PHASOR-KERNEL PHASE-ALIGNMENT
function aCand = phaseAlignment(grid,dirx,diry,alignItr,active,aCand,...
    nLayer,w)
% --------------------- PRECOMPUTE KERNEL DISTANCES AND NEIGHBOURHOOD SIZES
idM = 1:grid.csize(3); idK = 1:nLayer;                                     % Kernel and layer index sets
idofMat = aCand.ids.idofMat; idN = aCand.idN; pd = aCand.ids.pd;           % Extract sets of potential kernel neighbours
[nopad,nmax] = deal(aCand.ids.nopdMat(idofMat),(2*pd+1)^2);                % Linearised indexes wo. padding and maximum neighbourhood size
weight = zeros(nmax,grid.csize(3),nLayer); orient = weight;                % Array to store inter-neighbourhood weights and orientation relation
indx = false(nmax,grid.csize(3),nLayer);                                   % Array to store indexes of active neighbours
xcl = zeros(grid.csize(1:2)+2*pd); ycl = xcl; pact=xcl;                    % Initialise padded arrays
xcl(pd+1:end-pd,pd+1:end-pd) = grid.X_c;                                   % Pad coordinate arrays
ycl(pd+1:end-pd,pd+1:end-pd) = grid.Y_c;
[xldista,yldista]=deal(grid.X_c(:)-xcl(idofMat),grid.Y_c(:)-ycl(idofMat)); % Compute candidate neighbourhood distances
row_ind = repmat(idM',1,nmax);                                             % Centre kernel index array
for k = idK
    idn = grid.rrx.*(xldista.*diry(:,k)-yldista.*dirx(:,k)).^2+...         % Select neighbours within anisotropic radius
        grid.rry.*(xldista.*dirx(:,k)+yldista.*diry(:,k)).^2 <= grid.R;
    pact(pd+1:end-pd,pd+1:end-pd) = reshape(active(:,k),grid.csize(1:2));
    idn = idn & pact(idofMat);                                             % Indexes of active kernels within neighbourhood radius
    select = nopad(idn); kernel = row_ind(idn);                            % Selected neihgbours and corresponding centre kernel
    dirxn = dirx(select,k); diryn = diry(select,k);                        % Extract neighbour orientations
    ndot = dirxn.*dirx(kernel,k)+diryn.*diry(kernel,k);                    % Dot-product between centre kernel and selected neighbours
    nsgn = zeros(grid.csize(3),nmax); nsgn(idn) = sign(ndot);              % Sign of dot-product for pi-rotation correction
    nwght = zeros(grid.csize(3),nmax);                                     % Precompute neighbour weight distribution
    gdist = exp(-1/grid.R*(xldista(idn).^2+yldista(idn).^2));
    nwght(idn) = gdist.*exp(1i*(nsgn(idn).*(2*pi*grid.omega*...
        (dirxn.*(xldista(idn))+diryn.*(yldista(idn))))));
    weight(:,:,k) = nwght.'; orient(:,:,k) = nsgn.'; indx(:,:,k) = idn.';  % Save to layer-wise arrays
end
% --------------------------------------------------------- PHASE ALIGNMENT
pshift = 0*aCand.idbnd; pshift(aCand.idbnd(:)) = -pi/2;                    % Extract potential initial phase shift and initialise shift at boundary aligned kernels
for k = 1:nLayer
    [orientk,idnk,weightk] = deal(orient(:,:,k),indx(:,:,k),weight(:,:,k));% Extracting precomputed arrays for layer k
    aind = idM(active(:,k));                                               % Indexes of active kernels in layer
    [~,ordid] = sortrows([aCand.bdist(aind,k),w(aind,k)],[1,-2]);          % Get alignment order according to boundary distance and thickness
    for iter = 1:alignItr                                                  % Perform phase alignment for layer k
        for j = aind(ordid)
            idnj = idnk(:,j); wj = weightk(idnj,j);                        % Extracting precomputed weights for kernel j
            ortj = orientk(idnj,j); ortj2 = pi*(1-ortj)/2;                 % Extracting orientation correction information
            idNj = idN(j,idnj);                                            % Global linearised indexes of kernel neighbours
            avg_shift = sum(wj.*exp(1i*(ortj.*pshift(idNj,k)+ortj2)));     % Compute complex weigted average of neighbour contributions
            pshift(j,k) = atan2(imag(avg_shift),real(avg_shift));          % Get real-valued phase shift for kernel j from average
        end
    end
end
aCand.pshift = pshift;                                                     % Save phase shift for sampling
end
%% ======================================== FUNCTION TO SAMPLE PHASOR FIELD
function phasor = phasorSampling(grid,Dx,Dy,Dxx,Dxy,active,indicator_i,...
                                 align,nLayer)
phasor = zeros(grid.isize(3),nLayer); phasor(~indicator_i) = -4*pi*1i+1;   % Array to store layer-wise phasor fields initialised to ensure unsampled points not singular
max_active = max(active,[],2);                                             % Cover set for active kernels across layers
[idM,idK] = deal(1:grid.csize(3), 1:nLayer);                               % Kernel and layer index sets
[sCand,pshift] = deal(align.sample, align.pshift);                         % Extract candidate sampling regions and kernel phase shifts
for j = idM(max_active)                                                    % Only consider active kernels
    idXY = sCand(j).idXY;                                                  % Linearised global indexes of candidate sampling points
    [xdist,ydist] = deal(sCand(j).xdist(:,1), sCand(j).xdist(:,2));        % Centre to candidate sample-point distances
    for k=idK(active(j,:))
        [xdv,ydv] = deal(xdist*Dx(j,k), ydist*Dy(j,k));                    % Anisotropic distance measure and extract impacted candidates
        s = grid.rx*(xdist*Dy(j,k)-ydist*Dx(j,k)).^2+grid.ry*(xdv+ydv).^2;
        idcutoff = s < grid.cutoff; idxy = idXY(idcutoff(:));              % Remove sampling points from index sets by distance cut-off
        [xdv,ydv,s] = deal(xdv(idcutoff), ydv(idcutoff), s(idcutoff));
        d_jx = grid.omega*(Dx(j,k)-Dxx(idxy,k));                           % Orientation similarity sampling filter components
        d_jy = grid.omega*(Dy(j,k)-Dxy(idxy,k));
        d_jsq = (d_jx.^2+d_jy.^2);
        wfilt_j = exp(-s*grid.babpainv-pi^2*d_jsq*grid.bpainv + ...        % Kernel contribution w. integrated sampling filter
            2*1i*(pi*grid.omega*(xdv+ydv)+(d_jx.*xdist(idcutoff) + ...
            d_jy.*ydist(idcutoff))*grid.abpainv)+1i*pshift(j,k));
        phasor(idxy,k) = phasor(idxy,k)+(wfilt_j);                         % Add response for impacted candidates to complex phasor field
    end
end
end
%% ==================== FUNCTION FOR REMOVING SMALL DISCONNECTED COMPONENTS
function rho = removeIslands2D(rho,rad)
CC = bwconncomp(rho,4); numPixels = cellfun(@numel,CC.PixelIdxList);       % Connected component labeling
rho(vertcat((CC.PixelIdxList{numPixels < rad*sum(rho,'all')}))) = 0;       % Remove components with relatively few elements
end
%% ============= FUNCTION TO GET PRE-INTEGRATED STIFFNESS MATRIX COMPONENTS
function [KE0,BMatrixMid] = preIntegrateKE()
A = [2,1;1,2]*0.5; B = [0,0;0,1]*1/3; T=[0,1;1,0];                         % Compressed pre-integrated stiffness matrix components
KE11= kron([2,-1;-1,2]*1/6,[1,0,-1,0;0,0,0,0;-1,0,1,0;0,0,0,0]);
KE12= kron([1,-1;-1,1]*0.25,[0,1,0,1;1,0,-1,0;0,-1,0,-1;1,0,-1,0]);
KE13=-kron(T*1/6,[0,-1,0,1;-1,0,1,0;0,1,0,-1;1,0,-1,0])+...
    kron([1,-1;-1,1]*0.5,[3/2,1,0,-1;1,0,-1,0;0,-1,-3/2,1;-1,0,1,0]*2/3);
KE22=[kron(A,B),-kron(T*A,B);-kron(T*A,B),kron(A,B)];
KE23=[kron(A,T*1/3),-kron(T*A,T*1/3);-kron(T*A,T*1/3),kron(A,T*1/3)]+...
    kron([-1,1;1,-1]*0.5,[0,0,0,0;0,-1,0,0;0,0,0,0;0,0,0,1]);
D=kron(A,eye(2)*1/3);D2=kron(T*A,eye(2)*1/3);E=kron(T*0.5,[0,0;0,1]);
KE33=[D-E,-D2+E;-D2+E,D-E]+kron([1,-1;-1,1]*0.25,...
    [0,1,0,-1;1,0,1,0;0,1,0,-1;-1,0,-1,0]);
KE0 = cat(3,KE11(:),KE12(:),KE13(:),KE22(:),KE23(:),KE33(:));
BMatrixMid = 1/8*[-1,0,1,0,1,0,-1,0;0,-1,0,-1,0,1,0,1;-1,-1,-1,1,1,1,1,-1];% Element centre stress-displacement matrix
end
%% ===================================================== OC UPDATE FUNCTION
% The OC routine is from the code and paper by Ferrari, F., Sigmund, O. &
% Guest, J.K. Topology optimization with linearized buckling criteria in
% 250 lines of Matlab. Struct Multidisc Optim 63, 3045-3066 (2021).
function [x,as,lmid] = ocUpdate(loop,xT,xU,xL,dg0,g1,dg1,ocPar,xOld,...
    xOld1,as,beta,restartAs)
if (loop < 2.5 || restartAs == 1); as = xT+[-0.5,0.5].*(xU-xL)./(beta+1);
else; tmp = (xT-xOld).*(xOld-xOld1); gm = ones(length(xT),1);
    [gm(tmp>0), gm(tmp<0)] = deal(ocPar(2),ocPar(1));
    as = xT+gm.*[-(xOld-as(:,1)),(as(:,2)-xOld)];
end
xL = max( 0.9*as(:,1)+0.1*xT,xL); xU = min( 0.9*as(:,2)+0.1*xT,xU);
p0_0=(dg0>0).*dg0;q0_0=(dg0<0).*dg0;p1_0=(dg1>0).*dg1;q1_0=(dg1<0).*dg1;
[p0,q0] = deal(p0_0.*(as(:,2)-xT).^2,-q0_0.*(xT-as(:,1)).^2);
[p1,q1] = deal(p1_0.*(as(:,2)-xT).^2,-q1_0.*(xT-as(:,1)).^2);
primalProj = @(lm) min(xU,max(xL,(sqrt(p0+lm*p1).*as(:,1)+ ...
    sqrt(q0+lm*q1).*as(:,2))./(sqrt(p0+lm*p1)+sqrt(q0+lm*q1))));
psiDual = @(lm) g1 - ( (as(:,2)-xT)'*p1_0 - (xT-as(:,1))'*q1_0 ) + ...
    sum(p1./(max(as(:,2)-primalProj(lm),1e-12)) + ...
    q1./(max(primalProj(lm)-as(:,1),1e-12)));
lmUp = 1e6; x = xT; lmid = -1;
if psiDual(0)*psiDual(lmUp) < 0; lmid = fzero(psiDual,[0,lmUp]);
    x = primalProj( lmid );
elseif psiDual(0) < 0; lmid = 0; x = primalProj(lmid);
elseif psiDual(lmUp) > 0; lmid = lmUp; x = primalProj(lmid);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by R.V.Woldseth, O.Sigmund and P.D.L.Jensen %
%  TopOpt Group, Department of Civil and Mechanical Engineering,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please send your comments to: pdlj@dtu.dk                                %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "An 808 line phasor-based dehomogenisation Matlab code for multi-scale   %
%  topology optimization, R.V.Woldseth, O.Sigmund and P.D.L.Jensen, 2024"  %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the website: http://www.topopt.dtu.dk                    %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserve all rights but do not guarantee that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%