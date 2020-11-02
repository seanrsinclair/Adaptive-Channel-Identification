function [ output ] = band_pass_filter( signal, t, flow, fhigh, fsample )
    
    %filters signal sampled at fsample through 8th order Butterworth BPF
    %with band flow to fhigh
    tnorm=(fsample)*t;%normalixe t, not required, just for testing and visualization
    fn1=2*flow/fsample;%normalize flow
    fn2=2*fhigh/fsample;%normalixe fhigh
    [b, a]=butter(5, [fn1 fn2]); %design filter
    
    h=fvtool(b,a);%visualize filter

    output=filter(b,a,signal); %pass data through filter
    % plot(tnorm, output)
    %figure
    %plot(t, output) %plot filtered data

end

