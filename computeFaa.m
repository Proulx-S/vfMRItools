function [Faa,FaaLwr,FaaUpr] = computeFaa(dDoD,dVoV,method,alpha,nBoot)
    % faa from dV/V vs dD/D. With >1 output, also a confidence interval on Faa
    % (mapped from the slope CI), for the slope-fit methods.
    %   method : 'aproxSlope0' (default) | 'aproxSlope' | 'exact' | 'aprox'
    %   alpha  : CI significance (default 0.05 -> 95% CI)        [slope methods only]
    %   nBoot  : CI route (slope methods only):
    %            0   -> parametric CI (slopeCIparam: confint, from the fit covariance;
    %                   ~normal residuals; fast -- reuses the fit)            [default]
    %            >0  -> nonparametric paired bootstrap (slopeCIboot: bootci) with nBoot
    %                   resamples
    % Outputs:
    %   Faa            : faa estimate (scalar for slope methods; elementwise for exact/aprox)
    %   FaaLwr, FaaUpr : lower/upper CI bound on Faa ([] for exact/aprox). Faa = 1/2 -
    %                    1/4*slope is decreasing in slope, so the slope UPPER bound maps
    %                    to FaaLwr and the slope LOWER bound to FaaUpr.
    if ~exist('method','var') || isempty(method); method = 'aproxSlope0'; end
    if ~exist('alpha','var')  || isempty(alpha);  alpha  = 0.05;          end
    if ~exist('nBoot','var')  || isempty(nBoot);  nBoot  = 0   ;          end
    FaaLwr = []; FaaUpr = [];
    switch method
        case 'exact'
            % equation 1.7
            Faa = ( (1+dDoD).^2-dVoV-1 ) ./ ( (1-dVoV).*((dDoD+1).^4-1) );
        case 'aprox'
            % equation 1.8
            Faa = 1/2 - 1/4.*(dVoV./dDoD);
        case 'aproxSlope'
            f     = fit(dDoD(:),dVoV(:),'poly1');
            slope = f.p1;
            if nBoot==0
                [slopeLwr,slopeUpr] = slopeCIparam(f,alpha);
            else
                [slopeLwr,slopeUpr] = slopeCIboot(dDoD,dVoV,f,false,alpha,nBoot);
            end
            Faa    = 1/2 - 1/4*slope;
            FaaLwr = 1/2 - 1/4*slopeUpr;   % slope upper bound -> Faa lower bound (sign flip)
            FaaUpr = 1/2 - 1/4*slopeLwr;   % slope lower bound -> Faa upper bound
        case 'aproxSlope0'
            % intercept fixed to 0: fit through the origin (p2 constrained to 0)
            f     = fit(dDoD(:),dVoV(:),'poly1',fitoptions('poly1','Lower',[-Inf 0],'Upper',[Inf 0]));
            slope = f.p1;
            if nBoot==0
                [slopeLwr,slopeUpr] = slopeCIparam(f,alpha);
            else
                [slopeLwr,slopeUpr] = slopeCIboot(dDoD,dVoV,f,true,alpha,nBoot);
            end
            Faa    = 1/2 - 1/4*slope;
            FaaLwr = 1/2 - 1/4*slopeUpr;   % slope upper bound -> Faa lower bound (sign flip)
            FaaUpr = 1/2 - 1/4*slopeLwr;   % slope lower bound -> Faa upper bound
        otherwise
            error('computeFaa:badMethod','method must be ''exact'', ''aprox'', ''aproxSlope'' or ''aproxSlope0''');
    end
end

function [slopeLwr,slopeUpr] = slopeCIparam(f,alpha)
    % Parametric CI on the poly1 slope (p1), from the fit covariance (confint;
    % assumes ~normal residuals). Cheap -- reuses the already-computed fit.
    c = confint(f,1-alpha);     % cols = coeffs (p1,p2); rows = [lo; hi]
    slopeLwr = c(1,1);
    slopeUpr = c(2,1);          % the p2 column is NaN when the intercept is constrained
end

function [slopeLwr,slopeUpr] = slopeCIboot(x,y,~,throughOrigin,alpha,nBoot)
    % Nonparametric CI on the poly1 slope by paired (x,y) bootstrap (bootci), with a
    % closed-form slope so each resample is cheap (no re-fit). The 3rd arg (the fit
    % object f) is accepted for a call signature parallel to the fit-based path, but
    % unused -- the slope is recomputed per resample.
    x = x(:); y = y(:); ok = isfinite(x) & isfinite(y); x = x(ok); y = y(ok);
    if throughOrigin
        bfun = @(xx,yy) (xx.'*yy)./(xx.'*xx);     % OLS slope through origin
    else
        bfun = @(xx,yy) slopeFree(xx,yy);         % OLS slope, free intercept
    end
    c = sort(bootci(nBoot,{bfun,x,y},'Alpha',alpha));   % [lo; hi]
    slopeLwr = c(1);
    slopeUpr = c(2);
end

function b = slopeFree(x,y)
    % OLS slope with a free intercept (centered formula)
    x = x(:); y = y(:); mx = mean(x); my = mean(y);
    b = sum((x-mx).*(y-my)) ./ sum((x-mx).^2);
end
