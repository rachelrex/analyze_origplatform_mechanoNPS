function closestpeak=findclosestpeak(peaklist,manpt)
    %peaklist= a list of the locations of the local minima or maxima 
    %manpt= a point manually chosen by the user on the graph 
    %closestpeak= the minima/maxima closest to the manually chosen point
    
    %this function allows the user to pick a point on the plot close to the
    %local minima/maxima they want to select without having to be super
    %sure that they are actually clicking on the exact minima/maxima in
    %order to select it
    [numrow,~]=size(peaklist);
    distfrompt=ones(numrow,1).*manpt;
    distfrompt=abs(distfrompt-peaklist);
    [~,closestpeakind]=min(distfrompt);
    closestpeak=peaklist(closestpeakind,1);
    
    
    