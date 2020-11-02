function b= sinc_corr(ntaps,fmax,fs);
    if fmax >= .499*fs
        error('fmax must be less than .499*fs')
    end
    %    
    npts= 8;                % number of freq points for goal function
    k= 0:npts-1;
    f= k*fmax/(npts-1);                          % Hz
    hsinc= sin(pi*f/fs)./(pi*f/fs + eps);        % sinx/x response
    hsinc(1)= 1;
    h_goal= 1./hsinc;                            % goal function
    %    
    % least-squares FIR design
    ff(1:2:2*npts-1)= 2*f/fs;                  % f/fnyquist  vector of freq pairs
    ff(2:2:2*npts)= 2*f/fs + .001;
    a(1:2:2*npts-1)= h_goal;                   % vector of amplitude goal pairs
    a(2:2:2*npts)= h_goal;
    %    
    b= firls(ntaps-1,ff,a);                    % corrector coeffs