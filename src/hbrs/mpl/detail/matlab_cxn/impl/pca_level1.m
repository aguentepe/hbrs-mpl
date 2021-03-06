% Principal component analysis (PCA)
%
% Uses Singular value decomposition (SVD) algorithm.
% Ready for code generation using MATLAB Coder.
%
% References:
%  ${Matlab_ROOT_DIR}/toolbox/stats/stats/pca.m

function [coeff,score,latent,mu] = pca_level1(A, Economy, Center)
    coder.varsize('A', 'x');
    
    x = A; % make copy of A in order to keep A constant / untouched

    [m,n] = size(x);
    
    vWeights = ones(1,m,'like',x);
    vVariableWeights = ones(1,n,'like',x);
    
    DOF=m-Center;
    
    if Center
        % NOTE: Use local copy of 'wnanmean' because MATLAB CODER throws:
        %  ??? The function 'wnanmean' is not supported for standalone code generation.
        %  See the documentation for coder.extrinsic to learn how you can use this
        %  function in simulation.
        
        %mu = classreg.learning.internal.wnanmean(x, vWeights);
        % NOTE: If MATLAB cannot find wnanmean, then you have to add its parent folder your search path: 
        %       ${Matlab_ROOT_DIR}/toolbox/stats/classreg/+classreg/+learning/+internal/
        mu = wnanmean(x, vWeights);
    else
        mu = zeros(1,n,'like',x);
    end
    
    if Center
        x = bsxfun(@minus,x,mu);
    end
    
    OmegaSqrt = sqrt(vWeights);
    PhiSqrt = sqrt(vVariableWeights);
    
    x = bsxfun(@times, x, OmegaSqrt');
    x = bsxfun(@times, x, PhiSqrt);
    
    if Economy
        [U,S,coeff] = svd(x, 'econ');
    else
        [U,S,coeff] = svd(x, 0);
    end
    
    U = bsxfun(@times, U, 1./OmegaSqrt');
    coeff = bsxfun(@times, coeff, 1./PhiSqrt');
    S = diag(S);
    
    score =  bsxfun(@times,U,S');
    latent = S.^2./DOF;
    
    if DOF < n
        if Economy
            coeff(:, DOF+1:end) = [];
            score(:, DOF+1:end) = [];
            latent(DOF+1:end, :) = [];
        else
            %score(:, DOF+1:n) = 0; 
            % NOTE: this increases the size of score from m*m to m*n if m<n 
            %       and Economy=false
            % NOTE: MATLAB Coder generates bounds-violating C code for the
            %       line above so we have to workaroud:
            score_ = zeros([m,n]);
            score_(:,1:DOF) = score(:,1:DOF);
            score = score_;

            %latent(DOF+1:n, 1) = 0;
            latent_ = zeros(size(latent));
            latent_(1:DOF, 1) = latent(1:DOF, 1);
            latent = latent_;
        end
    end
    
    [~,maxind] = max(abs(coeff), [], 1);
    [d1, d2] = size(coeff);
    colsign = sign(coeff(maxind + (0:d1:(d2-1)*d1)));
    coeff = bsxfun(@times, coeff, colsign);
    score = bsxfun(@times, score, colsign);
end
