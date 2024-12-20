function [threshpts]=classicthresholding(yasls3,upperlimit,threshpts,i,k)
%i is the index of the data - can't be less than 2 or the last elseif will
%throw an error sometimes! (because i-1) 
%k is the index of the threshpts
%upperlimit is the line that yasls3 must pass to create a threshpts

[lengthpost,~]=size(yasls3); 
%if upperlimit is just one number (instead of a matrix), just make it a
%matrix with the same length as yasls3, and all the same number
[row,col]=size(upperlimit);
if row==1 && col==1
   upperlimitmat=ones(lengthpost,1);
   upperlimitmat=upperlimitmat.*upperlimit;
   upperlimit=upperlimitmat;
end
while i<lengthpost-2 
    if yasls3(i)==upperlimit(i) % if a data point== the threshold
        threshpts(k,1)=i;
%         xline(a2,i,'k');
        %determine whether trending down or up? 
        if yasls3(i-1)>upperlimit(i) && yasls3(i+1)<upperlimit(i+1) %trending down 
            threshpts(k,2)=0;
        elseif yasls3(i-1)<upperlimit(i) && yasls3(i+1)>upperlimit(i+1) %trending up 
            threshpts(k,2)=1;
        else
            threshpts(k,2)=3;
        end
            
        k=k+1;
    elseif yasls3(i)>upperlimit(i) && yasls3(i+1)<upperlimit(i+1) && yasls3(i+2)<upperlimit(i+2) % if btwn two points going down
        threshpts(k,1)=i;
        threshpts(k,2)=0;
%         xline(a2,i,'k');
        k=k+1;
    elseif yasls3(i)<upperlimit(i) && yasls3(i+1)>upperlimit(i+1) && yasls3(i-1)<upperlimit(i-1) % if btwn two points going up 
        threshpts(k,1)=i+1;
        threshpts(k,2)=1;
%         xline(a2,i+1,'k');
        k=k+1;
    end
    
    i=i+1;
end
