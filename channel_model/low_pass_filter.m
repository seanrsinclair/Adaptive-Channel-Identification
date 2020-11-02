function [rbi,rbq]=low_pass_filter(ri,rq,fs,fhigh)
    

    fcut=fhigh/fs;%normalize fhigh
        [b, a]=butter(2,fcut); %design filter
    
        h=fvtool(b,a);%visualize filter

    
    rbi = filter(b,a,ri); 
    rbq = filter(b,a,rq); 
end