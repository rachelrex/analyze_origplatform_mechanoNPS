function [truestartn,trueendn,pulsebottom2_mov,line_handles,mindiff,maxdiff,midpk]=findcorners_newway(ydetrend,yasdet,yas3det,ydiff2,yasdiff,xaxispost,vstepr,hstep,mininterval,zonenum,a4,a5,a,b,min_max_diff,known_pks_vals)

line_handles=gobjects(16,1);
%input variables: ydetrend, xaxispost, mindiff, maxdiff, minasls, vstepr, mininterval, horizline,
%a4

%output: truestart (calculated pulse start), trueend (calculated pulse
%end), pulsebottom_mov (2 column array with col1= xaxis indices and
%col2=moving average of "pulse bottom") 

%when does ydetrend cross the horizline? 
%for start look between (mindiff-vstepr) and minasls (which is the middle
%of the pulse) -- when ydetrend is<=horizline
%for end look between minasls and (maxdiff+vstepr)-- when ydetrend is>=
%horizline for the rest of the region 

%firststart=when ydetrend first crosses the horizline
%firstend=when ydetrend last crosses the horizline 
%"pulse bottom"= all points between firststart and firstend
 
%find horizline
 %make horizontal line
 
%if min_max_diff has a nonzero in the first column it means that the
%yasdiff valley was incorrect, if it is in the second column the peak was
%incorrect
               
[~,minasls]=min(yasdet(a+3:b-3));
minasls=minasls+a+3-1;
if min_max_diff(1,1)~=0 
    [~,locvas]=findpeaks(yasdiff(a:minasls).*-0.0001,'SortStr','descend');
    locvas=locvas+a-1;
    if min_max_diff(1,1)==-1
        mindiff=locvas(2,1);
    else
        %find the peak closest to the manually chosen point 
        closestpeak=findclosestpeak(locvas,min_max_diff(1,1));
        mindiff=closestpeak;
    end
else
    [~,mindiff]=min(yasdiff(a:minasls));
    mindiff=mindiff+a-1;
end
if min_max_diff(1,2)~=0
    [~,locpkas]=findpeaks(yasdiff(minasls:b).*0.0001,'SortStr','descend');
    locpkas=locpkas+minasls-1;
    if min_max_diff(1,2)==-1
        maxdiff=locpkas(2,1);
    else
        %find the peak closest to the manually chosen point 
        closestpeak=findclosestpeak(locpkas,min_max_diff(1,2));
        maxdiff=closestpeak;
    end
else
    [~,maxdiff]=max(yasdiff(minasls:b));
    maxdiff=maxdiff+minasls-1;
end



middlesize=(round((maxdiff-mindiff)/2))+mindiff;
minasls=middlesize;% just to work with later notation

line_handles(15,1)=xline(a4,minasls,'r');

%ensure that hstep does not exceed the threshpts
if minasls-hstep<a+vstepr || minasls+hstep>b-vstepr
    if b-vstepr-minasls > minasls-a+vstepr
        hstep=minasls-a+vstepr;
    else 
        hstep=b-vstepr-minasls;
    end
end


horizregion=yasdet(minasls-hstep:minasls+hstep); %was ydetrend, why did I change it? 
meanregion=mean(horizregion);
xlabel(a4,strcat('hstep=',num2str(hstep),' zone number=',num2str(zonenum)));


horizline=meanregion;
line_handles(16,1)=yline(a4,horizline,'r');
%end find horizontal line

%find first estimation of start 
firstregion=ydetrend(mindiff-vstepr:minasls);
firststart_found=0;
firststart=1;
n=1;
while firststart_found==0 && n<length(firstregion)
    if firstregion(n,1)<=horizline
        firststart_found=1;
        firststart=n;
    else
        n=n+1;
    end
    
end

firststart=firststart+(mindiff-vstepr)-1;

line_handles(1,1)=xline(a4,firststart,'Color','#77AC30');

%find first estimation of end 
firstregionend=ydetrend(minasls:maxdiff+vstepr);

firstend=1;
n=length(firstregionend);
firstend_found=0;
while  n>1 && firstend_found==0
    if firstregionend(n,1)<=horizline
        firstend_found=1;
        firstend=n;
    else
        n=n-1;
    end
end

firstend=firstend+minasls-1;
line_handles(2,1)=xline(a4,firstend,'Color','#77AC30');

% calculate second start and ends
pulsebottom=ydetrend(firststart:firstend);


pulsebottom_mov=smoothdata(pulsebottom,'movmean',mininterval/2); %was just mininterval

line_handles(3,1)=plot(a4,xaxispost(firststart:firstend),pulsebottom_mov,'LineWidth',1.5);




%try something for special case 
%find midpoint for diagonal calc
[~,locvas]=findpeaks(yasdiff(a:minasls).*-0.0001,'SortStr','descend');
locvas=locvas+a-1;
if locvas(1,1) ~=mindiff %was just locvas, but I think that was a mistake 
    'flag yasdiff valley'
    locvas(1,1)=mindiff;
end
[val,locv]=findpeaks(ydiff2(a:minasls).*-0.0001,'SortStr','descend');
locv=locv+a-1;
val=val*-1;
if length(locv)>=2
    for i=1:2
    plot(a5,locv(i,1),val(i,1),'x')
    end
end

if known_pks_vals(1,1)==0 
    i=1;
    while abs(locvas(1,1)-locv(i,1))>vstepr*3
        i=i+1;
    end
    midval=locv(i,1);

elseif  known_pks_vals(1,1)==-1
    [numlocv,~]=size(locv);
    k=1;
    i=1;
    first=0;
    while i<=numlocv && k<3
        if abs(locvas(1,1)-locv(i,1))<=vstepr*3
            k=k+1;
            first=i;
        end
        i=i+1;
    end
    if k<3
        midval=locv(first,1);
        'neverfound'
    else
        midval=locv(i-1,1); 
    end
else
    %find the peak closest to the manually chosen point 
    closestpeak=findclosestpeak(locv,known_pks_vals(1,1));
    midval=closestpeak;
end


%line_handles(4,1)=xline(a4,midval+(vstepr/2),'r');
%line_handles(5,1)=xline(a4,midval-(vstepr/2),'r');
p1v=polyfit(xaxispost(midval-(vstepr/2):midval+(vstepr/2)),ydetrend(midval-(vstepr/2):midval+(vstepr/2)),1);
x1v=xaxispost(midval-100:midval+(vstepr/2));
y1v=polyval(p1v,x1v);
line_handles(6,1)=plot(a4,x1v,y1v,'c','LineWidth',3.0);
firststartn=(horizline-p1v(2))/p1v(1);


[~,locpkas]=findpeaks(yasdiff(minasls:b).*0.0001,'SortStr','descend');
locpkas=locpkas+minasls-1;
if locpkas(1,1)~=maxdiff
    'flag yasdiff peak'
    locpkas(1,1)=maxdiff;
end
[pks,loc]=findpeaks(ydiff2(minasls:b).*0.0001,'SortStr','descend');
loc=loc+minasls-1;
if length(loc)>=4
    for i=1:4
    plot(a5,loc(i,1),pks(i,1),'x')
    end
end

if known_pks_vals(1,2)==0 
    %find midpoint for diagonal calc
    i=1;
    while abs(locpkas(1,1)-loc(i,1))>vstepr*3
        i=i+1;
    end
    midpk=loc(i,1);
elseif known_pks_vals(1,2)==-1
    [numloc,~]=size(loc);
    k=1;
    i=1;
    first=0;
    while i<=numloc && k<3
        if abs(locpkas(1,1)-loc(i,1))<=vstepr*3
            k=k+1;
            first=i;
        end
        i=i+1;
    end
    if k<3
        midpk=loc(first,1);
        'neverfound'
    else
        midpk=loc(i-1,1); 
    end
else
    closestpeak=findclosestpeak(loc,known_pks_vals(1,2));
    midpk=closestpeak;
end


%line_handles(7,1)=xline(a4,midpk+(vstepr/2),'r');
%line_handles(8,1)=xline(a4,midpk-(vstepr/2),'r');
p1=polyfit(xaxispost(midpk-(vstepr/2):midpk+(vstepr/2)),ydetrend(midpk-(vstepr/2):midpk+(vstepr/2)),1);
x1=xaxispost(midpk-100:midpk+(vstepr/2));
y1=polyval(p1,x1);
line_handles(9,1)=plot(a4,x1,y1,'c','LineWidth',3.0);
firstendn=(horizline-p1(2))/p1(1);

%need to round firststartn and firstendn
firststartn=floor(firststartn);
firstendn=ceil(firstendn);
%line_handles(10,1)=xline(a4,firststartn,'b','LineWidth',1.2);
%line_handles(11,1)=xline(a4,firstendn,'b','LineWidth',1.2);
if firststartn>firststart
    'flag start'
end
if firstendn<firstend
    'flag end'
end
if firststartn>firststart+3
    'error in start?'
end
if firstendn<firstend-3
    'error in end?'
end


pulsebottom2=ydetrend(firststartn:firstendn);
pulsebottom2_mov=smoothdata(pulsebottom2,'movmean',mininterval/2); %was just mininterval
%expand pulsebottom2_mov? 
% pulsebottom2_mov_beginning=mean(pulsebottom2_mov(1:3,1));
% pulsebottom2_mov_end=mean(pulsebottom2_mov(firstendn-firststartn+1-3:firstendn-firststartn+1,1));
% padding=zeros(vstepr,1);
% paddingbeginning=padding+pulsebottom2_mov_beginning;
% paddingend=padding+pulsebottom2_mov_end;
% pulsebottom2_mov=[paddingbeginning;pulsebottom2_mov;paddingend];
% plot(a4,xaxispost(firststartn-vstepr:firstendn+vstepr),pulsebottom2_mov,'LineWidth',1.5)
%end expand

line_handles(12,1)=plot(a4,xaxispost(firststartn:firstendn),pulsebottom2_mov,'LineWidth',1.5);

%find where the first moving average pt = the diagnoal line val 
truestartn=(pulsebottom2_mov(1,1)-p1v(2))/p1v(1);
truestartn=ceil(truestartn); %make smaller region
line_handles(13,1)=xline(a4,truestartn,'LineWidth',3,'Color','#7E2F8E');
trueendn=(pulsebottom2_mov(firstendn-firststartn+1,1)-p1(2))/p1(1);
trueendn=floor(trueendn); %make smaller region
line_handles(14,1)=xline(a4,trueendn,'LineWidth',3,'Color','#7E2F8E');

%add lines for a5 
line_handles(19,1)=xline(a5,mindiff,'m','LineWidth',1.5);    
line_handles(20,1)=xline(a5,maxdiff,'m','LineWidth',1.5);
line_handles(21,1)=xline(a5,midpk,'k','LineWidth',1.5);    
line_handles(22,1)=xline(a5,midval,'k','LineWidth',1.5);    

% or where it crosses the diagonal lines? are these points the same? 
%find truestartn
% movstart=firststartn-vstepr; 
% i=1;
% if movstart<x1v(1,1)
%     i=x1v(1,1)-movstart+1;
%     'movstart is less than x1v'
% end
% 
% while pulsebottom2_mov(i,1)<






    
    
    