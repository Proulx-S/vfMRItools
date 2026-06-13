function [Faa,xintrcp,yintrcp] = computeFaa(dDoD,dVoV,method)
    if ~exist('method','var') || isempty(method); method = 'exact'; end
    switch method
        case 'exact'
            % equation 1.7
            Faa = ( (1+dDoD).^2-dVoV-1 ) ./ ( (1-dVoV).*((dDoD+1).^4-1) );
        case 'aprox'
            % equation 1.8
            Faa = 1/2 - 1/4.*(dVoV./dDoD);
        case 'aproxSlope'
            f    = fit(dDoD(:),dVoV(:),'poly1');
            Faa  = 1/2 - 1/4*f.p1;
            yintrcp = f.p2;          % dV/V at dD/D=0 (fit y-intercept)
            xintrcp = -f.p2./f.p1;   % dD/D at dV/V=0 (fit x-intercept)
        case 'aproxSlope0'
            % intercept fixed to 0: fit through the origin (p2 constrained to 0)
            f    = fit(dDoD(:),dVoV(:),'poly1',fitoptions('poly1','Lower',[-Inf 0],'Upper',[Inf 0]));
            Faa  = 1/2 - 1/4*f.p1;
            yintrcp = f.p2;          % dV/V at dD/D=0 (fit y-intercept) = 0 by construction
            xintrcp = -f.p2./f.p1;   % dD/D at dV/V=0 (fit x-intercept) = 0 by construction
        otherwise
            error('computeFaa:badMethod','method must be ''exact'' or ''aprox''');
    end
