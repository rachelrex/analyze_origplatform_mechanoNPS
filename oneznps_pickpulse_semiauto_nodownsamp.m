function [pulses,threshpts,fig_handles]=oneznps_pickpulse_semiauto_nodownsamp(ydetrend,yasls,yas2det,dmin,stderror,yas3det,plot_thresh,skipto,use_old_fig,fig_handles)
% These dimensions MUST be changed to reflect your specific device 
L=5150; %microns, overall length of device (including entrance node and exit node) 
szlength = 800; %microns, length of sizing channel 
sqlength = 800; %microns, length of the squeeze channel 
Deff= 29.4632;

filesaveto='emergencyrecovered';
toobig=18; %SHOULD BE REORGANIZED, the size of the cell where its so big you dont want to analyze the squeezes 
%use_old_fig is 0 or 1, depending on if you want to use the figures you already have up from a previous round, so you dont have to move the windows around  
%skipto lets you start where you left off - set to 0 if you want to start
%from beginning 

%These may need to be changed depending on your geometery and flow speed
mininterval=100; %was 200(USED IN fixthreshptserrormz) the minimum number of data points the current must cross the threshold for, before it can be counted as a possible pulse 
eventlength=4000; %(USED IN aretheytogether) approximate length (in data points) of complete cell event (all 3 zones), if this is too big you will miss events, if too small you'll catch same event more than once 
wp=eventlength+2200; %was just +200 (for vizualizing entire events) 

%dmin=minimum size you're interested in, will be used to set the "limit" of
%detection for thresholding
%minexpected= the deltR for the threshold, can be plot on the ydetrended
%graph
minexpected=(((dmin^3)/((Deff^2)*L))*(1/(1-(0.8*((dmin/Deff)^3)))))*yasls;
limit=0.9*minexpected;



[ma,~]=size(yasls);
xaxis=linspace(1,ma,ma)';

%% find where the detrended signal crosses the threshold (called limit, determined by dmin) 

%thresholding function originally written for when current signal crosses the threshold, current decreases with particle (current drop), while
%resistance increases, so to use the same thresholding function with
%resistance we *-1
yas2det_for_thresh=-1*yas2det; 
limit_for_thresh=-1*limit;
j=100+wp; 
k=1;
threshpts=zeros(2,2); %everytime the signal crosses the threshold line
threshptsc=zeros(2,2); %corrected- everytime the signal crosses the threshold line, and stays above it for at least mininterval (ensures its a true pulse, not noise) 

threshptszone=classicthresholding(yas2det_for_thresh(:,1),limit_for_thresh(:,1),threshpts(:,1),j,k);
[m,~]=size(threshptszone);
threshptszonec=fixthreshptserrormz(threshptszone,mininterval,k,m);
[mc,~]=size(threshptszonec);
threshpts(1:m,1:2)=threshptszone;
threshptsc(1:mc,1:2)=threshptszonec;


%threshpts in this code just picks out hte beginning and end of events in
%each zone- which means sizing, and sq+recov


%plot threshpts
%if plot_thresh=1--> plot all the lines that where threshptsc start
%(beginning of pulse, not end) 
%if plot_thresh=2 --> plot threshptsc and threshpts (not corrected), you
%would do this if you're worried that you're missing cell events with the
%dmin you set 
if use_old_fig==0
    fig_handles=gobjects(3,1);
    fa=figure('Name','Data Set');
    fig_handles(1,1)=fa;
else 
    fa=fig_handles(1,1);
end
sz=5;
linecolor=["#D95319";"#7E2F8E";"#77AC30"];
linew=0.75;
figure(fa)
tl=tiledlayout('flow');
ax1=nexttile;
hold(ax1,'on')


    scatter(xaxis,ydetrend(:,1),sz,'MarkerEdgeColor',linecolor(1,1))
    plot(xaxis,ydetrend(:,1),'Color',linecolor(1,1))

plot(xaxis,minexpected(:,1),'y');
plot(xaxis,limit(:,1),'c');

% %consider uncommenting this 
% if plot_thresh==2
%     ax2=nexttile;
%     hold(ax2,'on')
%     for i=1:3
%         scatter(xaxis,ydetrend(:,i),sz,'MarkerEdgeColor',linecolor(i,1))
%         plot(xaxis,ydetrend(:,i),'Color',linecolor(i,1))
%     end
%     plot(xaxis,minexpected(:,1),'y');
%     plot(xaxis,limit(:,1),'c');
% end
% 
% 
% if plot_thresh==1 || plot_thresh==2
%     for i=1:length(threshptsc)
%         for j=1:3
%             xline(ax1,threshptsc(i,j*2-1),'Color',linecolor(j,1),'LineWidth',linew)
%     %         xline(ax1,threshptsc(i,j*2),'Color',linecolor(j,1),'LineWidth',linew)
%         end
%     end
% end
% if plot_thresh==2
%     for i=1:length(threshpts)
%         for j=1:3
%             xline(ax2,threshpts(i,j*2-1),'Color',linecolor(j,1),'LineWidth',linew)
%     %         xline(ax2,threshpts(i,j*2),'Color',linecolor(j,1),'LineWidth',linew)
%         end
%     end
% end

%% prep for manually choosing 

[mc,~]=size(threshptsc);
threshptsc(mc+1,:)=zeros(1,length(threshptsc(1,:)));
mc=mc+1;

events=zeros(1,6);
threshz1=threshptsc(:,1:2);
events(1,1)=threshz1(1,1);
events(1,2)=threshz1(2,1);
stopit=0;
i=3;%identifies which row of thresh we are on 
j=2;%identifies which row of events we are on

while i<=mc && stopit==0 
   if threshz1(i,1)==0  %at the end of the thresh there will just be zeros, so we know we've come to the end
       stopit=1;
   else
       if threshz1(i,1)-events(j-1,1)<eventlength %if there is a pulse within the "eventlength" of an identified pulse then we ignore the second pulse 
           i=i+2;
       else % if the next pulse is after the last events entire "eventlength" then it is the start of a new event
           events(j,1)=threshz1(i,1);
           events(j,2)=threshz1(i+1,1);
           i=i+2;
           j=j+1;
       end
   end
end




[me,~]=size(events);
pulses=zeros(me,37); %for each column of pulses: 
%1: startsz 2:endsz 3:startsq1 4:endsq1 5:startsq2 6:ends2 7:startsq3
%8:endsq3 9:avgdeltR sizing 10: avg baselineR sizing 11:size 12:recovstart1
%13:recovend1 14:recovstart2 15:recovend2 16:recovstart3 17:recovend3 last
%3 are recovmagcheck
recovmagcheck=ones(me,3);
axesbounds=zeros(me,8);

%% semi-auto choosing 
if use_old_fig==0
    fb=figure('Name','Entire Event'); 
    fig_handles(2,1)=fb;
else
    fb=fig_handles(2,1);
end
figure(fb)
tlb=tiledlayout('flow');
axb1=nexttile; %it will show the entire event with the threshold lines
axar=nexttile;

%for window of interest
if use_old_fig==0
    f4= figure('Name','Window of Interest');
    fig_handles(3,1)=f4;
else
    f4=fig_handles(3,1);
end
figure(f4)
a4=subplot(3,1,1);
a5=subplot(3,1,2);
a6=subplot(3,1,3);

i=1;
p=1; %number of recorded pulses we're at 

%just make it manual now -- make sure endz doesnt exceed array bounds,
%yasdet2 may need to be adjusted to follow more closely, maek sure you can
%adjust the xlimits of z plots when manually picking 
while i<=me 
        cla(axb1,'reset')
        cla(axar,'reset')
        title(axb1,num2str(p))
         startev=events(i,1)-100;
         endev=startev+wp;
         if endev>ma
             endev=ma;
         end
         hold(axb1,'on')
         line_handles=gobjects(3,3);
         for j=1:1 %was 3
             %axb1
            scatter(axb1,xaxis(startev:endev,1),ydetrend(startev:endev,j),sz,'MarkerEdgeColor',linecolor(j,1))
            plot(axb1,xaxis(startev:endev,1),ydetrend(startev:endev,j),'Color',linecolor(j,1))
            
            
            %zone tiles
            if j==1
                startz=startev;
            else
                startz=events(i,j*2-2);
            end
            if j==3 
                endz=events(i,j*2)+100;
                if endz<endev
                    endz=endev;
                end
            else
                endz=events(i,j*2+2);
            end
            
            if startz==0
                startz=startev;
                'start was zero'
            end
            if endz==0 || endz==100
                endz=endev;
                'end was zero'
            end
            if endz>ma
                endz=ma;
            end
            
            hold(axar,'on')
            scatter(axar,xaxis(startz:endz,1),ydetrend(startz:endz,j),sz,'MarkerEdgeColor',linecolor(j,1))
            plot(axar,xaxis(startz:endz,1),ydetrend(startz:endz,j),'Color',linecolor(j,1))
            plot(axar,xaxis(startz:endz,1),yas2det(startz:endz,j),'b')
            
            
            if endz<startz
                'endz is less than startz'
            end
            
            %plot all thresholds in axar 
            %find threshpts
             axesbounds(p,1:2)=xlim(axb1(1,1)); %saves 
                axesbounds(p,3:4)=xlim(axar);
        
                    threshptszonec=threshptsc(:,1:2);
                    i1=1;
                    found=0;
                    while found==0 && i1<=length(threshptszonec)
                        if axesbounds(p,1*2+1)<=threshptszonec(i1,1)
                            found=1;
                        else
                            i1=i1+1;
                        end
                    end
                    threshrownumstart=i1;

                    found=0;
                    while found==0 && i1<=length(threshptszonec)
                        if axesbounds(p,1*2+2)<=threshptszonec(i1,1)
                            found=1;

                        elseif threshptszonec(i1,1)==0
                            found=1;
                        else
                            i1=i1+1;
                        end
                    end
                    threshrownumend=i1-1;

                    threshptseventzone=threshptszonec(threshrownumstart:threshrownumend,:);
                    [mc,~]=size(threshptseventzone);
                    %threshpts for this cell event in this zone found 
             %we've found all the threspts in axar        
                    
            if events(i,j*2-1)~=0
                xline(axb1,events(i,j*2-1),'Color',linecolor(j,1),'LineWidth',linew)
                line_handles(j,1)=xline(axar,events(i,j*2-1),'Color',linecolor(j,1),'LineWidth',linew);
                line_handles(j,3)=plot(axar,xaxis(startz:endz,1),limit(startz:endz,j),'k');
            end
            if events(i,j*2)~=0
                xline(axb1,events(i,j*2),'Color',linecolor(j,1),'LineWidth',linew)
                line_handles(j,2)=xline(axar,events(i,j*2),'Color',linecolor(j,1),'LineWidth',linew);
            end
             %now plot all thresh lines in axar
                     for x=1:mc
                         line_handles(x+1,1)=xline(axar,threshptseventzone(x,1),'k','LineWidth',linew);
                     end

         end
         plot(axb1,xaxis(startev:endev,1),minexpected(startev:endev,1),'y');
         plot(axb1,xaxis(startev:endev,1),limit(startev:endev,1),'c');
         % manually pick points
         
         if skipto==0
            choice=input('Type 3 to skip to later on in the data run, press enter to keep event, type 1 to discard event, type 0 to stop analyzing and exit');
         else
             choice=3;
         end
         
         if choice ==3
             if skipto==0
                skipto= input('type the approximate index number you wish to skip to');
             end
             %debug for accidental inputs
             while isnumeric(skipto)==0 || isempty(skipto)==1
                 skipto= input('type the approximate index number you wish to skip to');
             end
             while skipto>events(me,1) || skipto<=events(i,1)
                 skipto= input('your number is either too large or too small, try again');
             end
             %debugging over
             
             foundit=0;
             while events(i,1)<skipto && foundit==0
                 i=i+1; 
                 if events(i,1)>=skipto
                     foundit=1;
                     i=i-1;
                 end
             end
             skipto=0;
             
             
         elseif isempty(choice) %choice ~=1 && choice ~=0
             adjustq=1;
             while isempty(adjustq)==0
                 adjustq=input('Do you need to adjust the window? If no press enter, if yes type 1');
                 if isempty(adjustq)==0
                     adjustq=1;
                    currentsize=xlim(axar);
                    startnow=currentsize(1,1); endnow=currentsize(1,2);  
                    windowresize= input('Add 1000 points to left=1,Add 2000 points to right=2,Add 500 points to both sides=3');
                    switch windowresize
                        case 1
                            startnow=startnow-1000;
                        case 2 
                            endnow=endnow+2000;
                        case 3
                            startnow=startnow-500;
                            endnow=endnow+500;
                    end
                    cla(axar)
                    hold(axar,'on')
                    scatter(axar,xaxis(startnow:endnow,1),ydetrend(startnow:endnow,adjustq),sz,'MarkerEdgeColor',linecolor(adjustq,1))
                    plot(axar,xaxis(startnow:endnow,1),ydetrend(startnow:endnow,adjustq),'Color',linecolor(adjustq,1))
                    plot(axar,xaxis(startnow:endnow,1),yas2det(startnow:endnow,adjustq),'b')
                    if events(i,adjustq*2-1)~=0
                        line_handles(adjustq,1)=xline(axar,events(i,adjustq*2-1),'Color',linecolor(adjustq,1),'LineWidth',linew);
                        line_handles(adjustq,3)=plot(axar,xaxis(startnow:endnow,1),limit(startnow:endnow,adjustq),'k');
                    end
                    if events(i,adjustq*2)~=0
                        line_handles(adjustq,2)=xline(axar,events(i,adjustq*2),'Color',linecolor(adjustq,1),'LineWidth',linew);
                    end
                 end
             end
            % record the correct axes bounds for all 4 graphs *MAYBE BAD*
            axesbounds(p,1:2)=xlim(axb1(1,1)); %saves 
            for x=2:4
                axesbounds(p,x*2-1:x*2)=xlim(axar);
            end 
            
            
            %redo threshpts for larger pulse if needed
            threshredone=0;
            threshgood=0;
            while threshgood==0
                redo=input('Do we need to redo threshpts?, type 1=yes, type 2=no, to make a note that a recovery magnitude is unreliable type 5');
                if isempty(redo)==1
                    threshgood=1;
                elseif redo==1
                    delete(line_handles) %delete the threshpoint lines on all the individual zone graphs 
                    line_handles=gobjects(3,3);
                    newdmin=input(strcat('Current dmin=',num2str(dmin),'type whatever you want the new one to be'));
                    new_minexpected=(((newdmin^3)/((Deff^2)*L))*(1/(1-(0.8*((newdmin/Deff)^3)))))*yasls;
                    new_minexp_for_thresh=-1*new_minexpected;
                    %recalculate threshpts with "new_minexpected" as the
                    %threshold
                    new_threshpts=zeros(2,2);
                    new_threshptsc=zeros(2,2);
                    for j1=1:1
                     start_thresh=axesbounds(p,j1*2+1)-1;
                     end_thresh=axesbounds(p,j1*2+2);
                     new_threshptszone=zeros(2,2);
                     line_handles(1,j1)=plot(axar,xaxis(start_thresh:end_thresh,1),new_minexpected(start_thresh:end_thresh,j1),'k');
                     new_threshptszone=classicthresholding(yas2det_for_thresh(start_thresh:end_thresh,j1),new_minexp_for_thresh(start_thresh:end_thresh,j1),new_threshptszone,2,1);
                     [m,~]=size(new_threshptszone);
                     new_threshptszonec=fixthreshptserrormz(new_threshptszone,mininterval,1,m);
                     [mc,~]=size(new_threshptszonec);
%                      new_threshpts(1:m,j1*2-1)=new_threshptszone(1:m,1)+start_thresh-1;
%                      new_threshpts(1:m,j1*2)=new_threshptszone(1:m,2);
                     if new_threshptszonec(1,1)==0
                         'no threshpts found'
                     else
                         new_threshptsc(1:mc,j1*2-1)=new_threshptszonec(1:mc,1)+start_thresh-1;
                         new_threshptsc(1:mc,j1*2)=new_threshptszonec(1:mc,2);
                     end 
                     %now plot new thresh lines
                     for x=1:mc
                         line_handles(x+1,j1)=xline(axar,new_threshptsc(x,j1*2-1),'Color',linecolor(j1,1),'LineWidth',linew);
                     end
                    end
                    threshredone=1;
                    %etc. basically recalculate threshpts, then replot,and ask
                    %again if its good now - until we get what we want 
                    %find threshpts within the axes bounds for each zone with a
                    %new threshold
                elseif redo==5
                    recovq=input('is it the first, second or third recovery thats bad?, type 1,12,123,2,23,etc.');
                    if isempty(recovq)
                        threshgood=1;
                    else
                       %determine if all recov mag are unreliable 
                       switch recovq
                           case 123
                               recovmagcheck(p,1:3)=[0 0 0];
                           case 12
                               recovmagcheck(p,1:2)=[0 0];
                           case 13
                               recovmagcheck(p,1)=0;
                               recovmagcheck(p,3)=0;
                           case 23
                               recovmagcheck(p,2:3)=[0 0];
                           case 1
                               recovmagcheck(p,1)=0;
                           case 2
                               recovmagcheck(p,2)=0;
                           case 3
                               recovmagcheck(p,3)=0;
                       end
                    end
                else
                    threshgood=1;
                end
            end
            
            
            %do semi-auto corner picking
            cornernums=zeros(3,24);
            cellsize=0; %for each new cell event reset cellsize to 0 because it hasnt been calculated yet 
            
            for j1=1:1 %was 1:3 but now that its just mechano with one zone i only need to look at it once 
                %find the threshpts for this cell event in this zone
                %if threshredone=1, then they are already known 
                %eventnum=i
                zonenum=1;
                %if cell size is bigger than 15 microns no need to analyze
                %anymore, because I definitely can't include it
                if cellsize>=toobig
                    break
                end
                if threshredone==0
                    %find threshpts if threshredone=0
                    threshptszonec=threshptsc(:,1:2);
                    i1=1;
                    found=0;
                    while found==0 && i1<=length(threshptszonec)
                        if axesbounds(p,1*2+1)<=threshptszonec(i1,1)
                            found=1;
                        else
                            i1=i1+1;
                        end
                    end
                    threshrownumstart=i1;

                    found=0;
                    while found==0 && i1<=length(threshptszonec)
                        if axesbounds(p,1*2+2)<=threshptszonec(i1,1)
                            found=1;

                        elseif threshptszonec(i1,1)==0
                            found=1;
                        else
                            i1=i1+1;
                        end
                    end
                    threshrownumend=i1-1;

                    threshptseventzone=threshptszonec(threshrownumstart:threshrownumend,:);
                    %threshpts for this cell event in this zone found 
                else
                    threshptseventzone=new_threshptsc(:,1*2-1:1*2);
                    [m,~]=size(threshptseventzone);
                    x=1;
                    while x<=m
                        if threshptseventzone(x,1)==0
                            threshptseventzone=threshptseventzone(1:x-1,:);
                            x=m+1;
                        else
                            x=x+1;
                        end
                       
                    end
                    [m,~]=size(threshptseventzone);
                    if threshptseventzone(m,2)==0
                        threshptseventzone=threshptseventzone(1:m-1,:);
                    end
                end
                % make inverted resistance signal 
                ydetrendcorner=ydetrend(:,zonenum).*-1;
                yasdetcorner=yas2det(:,zonenum).*-1;
                yas3detcorner=yas3det(:,zonenum).*-1;

                ydiffcorner=[0;diff(ydetrendcorner)];
                ydiffcorner=(ydiffcorner.^3);

                yasdiffcorner=diff(yasdetcorner);
                yasdiffcorner=[0;yasdiffcorner];
                yasdiffcorner=(yasdiffcorner.^3);
                
                yas3diffcorner=diff(yas3detcorner);
                yas3diffcorner=[0;yas3diffcorner];
                yas3diffcorner=(yas3diffcorner.^3);

                limitcorner=limit(:,zonenum).*-1;
                minexpectedcorner=minexpected(:,zonenum).*-1;

                % run program 
                
                [cornernums(j1,:),cellsize] = mzonezone_corner_finder_nodownsamp(threshptseventzone,ydetrendcorner,yasdiffcorner,yasdetcorner,yas3detcorner,ydiffcorner,limitcorner,minexpectedcorner,dmin,stderror(1,zonenum),wp,j1,a4,a5,a6,f4,yasls,mininterval,yas3diffcorner,Deff,L,cellsize,axar);
            end
           
            pulses(p,1:4)=cornernums(1,1:4);
            
            pulses(p,5:24)=cornernums(1,5:24);
            
            
             %record resistance values for sizing pulse 
             pulses(p,25)=mean(ydetrend(pulses(p,1):pulses(p,2),1));
             pulses(p,26)=mean(yasls(pulses(p,1):pulses(p,2),1));
             pulses(p,27)=calcsize(Deff,L,pulses(p,25)/pulses(p,26));
             pulses(p,28)=pulses(p,25)/pulses(p,26); %deltaRR
             
             %for 28-37 it will be the deltaRR of each recovery subpulse 
             if pulses(p,5)>0
               %if the recovery was recorded   
                 for rcol=1:10
                     deltaRr=mean(ydetrend(pulses(p,3+(rcol*2)):pulses(p,4+(rcol*2)),1));
                     baseRr=mean(yasls(pulses(p,3+(rcol*2)):pulses(p,4+(rcol*2)),1));
                     pulses(p,28+rcol)=deltaRr/baseRr;
                 end
             end
             
             %recovmagcheck
             %pulses(p,18:20)=recovmagcheck(p,1:3);
             p=p+1;
         elseif choice == 0
             i=me+1;
         end
         
     
          
         
         
         i=i+1;
         hold(axb1,'off')
 save(filesaveto,'pulses')        
%add debug point at end  
end


