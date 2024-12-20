function [cornernums,cellsize] = mzonezone_corner_finder_nodownsamp(threshpts,ydetrend,yasdiff,yasdet,yas3det,ydiff2,upperlimit,expectedpore,dinput,stderror,wp,zonenum,a4,a5,a6,f4,yasls,mininterval,y3diff,Deff,L,cellsize,axar)
upperlimitnz=upperlimit+yasls;
ynz=ydetrend+yasls;
expectedporenz=expectedpore+yasls;
dinterest=10;
dinterest_line=((((dinterest^3)/((Deff^2)*L))*(1/(1-(0.8*((dinterest/Deff)^3)))))*yasls).*-1;
dinterest_line=dinterest_line+yasls;
yasdetnz=yasdet+yasls;
sz=5;
[rowthresh,~]=size(threshpts);
[lengthpost,~]=size(ydetrend);
xaxispost = linspace(1,lengthpost,lengthpost)';
xtitle='nope';
ytitle=num2str(cellsize); 
wp=round(wp/3); %wp/num zones 
toobig=30; %SHOULD BE REORGANIZED, the size of the cell where its so big you dont want to analyze the squeezes 
%% look at the windows 
% f4= figure('Name','Window of Interest');
% a4=subplot(2,1,1);
% a5=subplot(2,1,2);
hstep= 60; %horizontal line step size (og 10)15
hstepog=hstep;
vstep=3;%3
vstepr=16; %rising and falling lines step size (og 4)
steps=[hstep hstepog vstep vstepr];
cornercontext=zeros(1,23);
cornerindex=zeros(1,23);
cornerdiff=zeros(1,23);
guesscheck=zeros(1,4); %col1= startguess, col2=startactual, col3=endguess, col4=endactual
guesscheck(:,:)=-1;
i=1;
useranswer=1;
k=1; % the number of corner *sets* logged

%corner matrices, each row is different pulse segment (sizing, cont, recov 1-10), set of 13 rows can be from 1 cell 
%col: 1= pulse type; 2=actual start pulse; 3=actual endpulse; 4:13= the
%points surrounding the starting predicted corner (5 on each side, if it falls
%between two points); 14:23= the points surrounding the ending predicted
%corner 

while i<rowthresh && useranswer==1 && threshpts(i+1,1)+wp<lengthpost
    if threshpts(i,2)==1
        'bad start, skipping one theshpt'
        i=i+1;
    end
    a=threshpts(i,1)-vstepr; %vstepr used to be 3 
    b=threshpts(i+1,1)+vstepr; %vstepr used to be 3
    cla(a4)
    reset(a4)
    cla(a5)
    reset(a5)
    cla(a6)
    reset(a6)
        plot(a4,xaxispost(a:b),ydetrend(a:b),'r');
        plot(a5,xaxispost(a:b),ydetrend(a:b),'r');
        plot(a6,xaxispost(a-wp:b+wp),ynz(a-wp:b+wp),'r');
        linkaxes([a4,a5],'x');
        
        hold(a4,'on')
        hold(a5,'on')
        hold(a6,'on')
        yline(a5,0,'k');
        plot(a6,xaxispost(a-wp:b+wp),yasls(a-wp:b+wp),'k')
        plot(a6,xaxispost(a-wp:b+wp),upperlimitnz(a-wp:b+wp),'c')
        plot(a6,xaxispost(a-wp:b+wp),expectedporenz(a-wp:b+wp),'k')
        plot(a6,xaxispost(a-wp:b+wp),yasdetnz(a-wp:b+wp),'b','LineWidth',2)
       % plot(a6,xaxispost(a-wp:b+wp),dinterest_line(a-wp:b+wp),'Color',"#77AC30",'LineWidth',1.5)  
%         plot(a6,xaxispost(a-wp:b+wp),maxlimit(a-wp:b+wp)-yasls(a-wp:b+wp),'m')
%         plot(a6,xaxispost(a-wp:b+wp),sqlimit(a-wp:b+wp)-yasls(a-wp:b+wp),'g')
        xlabel(a6,strcat('expected pore, in black, =',num2str(dinput)));
        title(a4,strcat(' zone number=',num2str(zonenum)));
        
        scatter(a4,xaxispost(a:b),ydetrend(a:b),sz,'r');
        xlabel(a4,xtitle);
        ylabel(a4,ytitle); 
       
        scatter(a5,xaxispost(a:b),yasdiff(a:b).*(10^9),sz,'m')%was *.0001, .*.001
        plot(a5,xaxispost(a:b),yasdiff(a:b).*(10^9),'m')
        scatter(a5,xaxispost(a:b),ydiff2(a:b).*(10^8),sz,'k')%was *.0001,
        plot(a5,xaxispost(a:b),ydiff2(a:b).*(10^8),'k')
       
        
        plot(a4,xaxispost(a:b),yasdet(a:b),'b','LineWidth',2)
        
        xline(a6,threshpts(i,1),'b')
        xline(a6,threshpts(i+1,1),'b')
        
       
   
    
   
    
    
    
    useranswer0=input('empty= analyze and guess what it is, bad start=0, next=1, analyze squeeze=2, look for recovery= 3,skip the rest of pulses in this zone= 5');
    
    
    if isempty(useranswer0)==1 || useranswer0==2 || useranswer0==3
            if isempty(useranswer0)==1 %this sets the default - when its empty set answer to default 
                useranswer0=2; %current default %was =3
            end
            
            %useranswer2 is whether this is sizing or squeeze, =1 is sizing sub-pulse, =2 is squeeze sub-pulse
            if k>1 && cornercontext(k-1,1)==1 
                %this is the squeeze pulse in zone 1 
                useranswer2=2;
                %calculate diameter of cell 
                sizestart=cornercontext(k-1,2);
                sizeend=cornercontext(k-1,3);
                meandrop=-1*(mean(ydetrend(sizestart:sizeend,1)));
                meanbase=mean(yasls(sizestart:sizeend,1));
                cellsize=calcsize(Deff,L,meandrop/meanbase);
                if cellsize>=toobig
                    'cell was too big' 
                    break
                end
                ylabel(a4,num2str(cellsize)); 
            else
                useranswer2=1;
            end
            
            switch useranswer2
                case 1 %sizing pulse
                    cornercontext(k,1)=1;
                    cornerindex(k,1)=1;
                    cornerdiff(k,2)=1;
                    hstep=hstep+40;
                case 2 %squeeze/ recovery pulses
                    cornercontext(k,1)=2;
                    cornerindex(k,1)=2;
                    cornerdiff(k,2)=2;
                    hstep=hstep+10;
            end
            if k<3 || cornercontext(k,1)==1 || cornercontext(k,1)==2
            %it is a sizing or squeeze pulse 
                
                min_max_diff=[0 0]; 
                % col1=min, col2=max, when =0 are correct or haven't been checked
                % yet, when=-1 the program has guessed once after being
                % wrong, when anything else= user had to manually pick 
                pks_vals_okay=[0 0]; %when=0 the end or start wasnt chosen correctly, when=1, it was chosen correctly
                known_pks_vals=[0 0];% when =0 are correct or haven't been checked
          
                while pks_vals_okay(1,1)==0 || pks_vals_okay(1,2)==0
                    [truestart,trueend,pulsebottom_mov,line_handles,mindiff,maxdiff,midpk]=findcorners_newway(ydetrend,yasdet,yas3det, ydiff2,yasdiff,xaxispost, vstepr,hstep, mininterval,zonenum,a4,a5,a,b,min_max_diff,known_pks_vals);
                    
                    line_handles(17,1)=xline(a6,truestart,'LineWidth',1.5,'Color','#D95319');    
                    line_handles(18,1)=xline(a6,trueend,'LineWidth',1.5,'Color','#D95319');
                    guesscheck(k,1)=truestart;
                    guesscheck(k,3)=trueend;



                     %set limits so you can clearly see the pulse 
                     [mindetrend,~]=min(ydetrend(a:b));
                     [maxdetrend,~]=max(ydetrend(a:b));
                    ylim(a4,[mindetrend,maxdetrend]);
                    if useranswer2==2
                        if mindiff-a>1000
                            xlim(a4,[mindiff-200,maxdiff+200]);
                             xlim(a5,[mindiff-200,maxdiff+200]);
                        else
                            xlim(a4,[a,maxdiff+20]);
                            xlim(a5,[a,maxdiff+20]);
                        end
                    end

                    %ask user to check start
                    appropriateanswer=0;
                    while appropriateanswer==0
                        useranswer3=input('if startguess is right- press enter, if you need to fix horizline type 1, if you need to repick the "valley" used to calculate the diagonal type 6, if you decided you dont want to analyze this type 3');

                        if useranswer3==1 %horizline is wrong 
                            delete(line_handles)
                            line_handles=gobjects(2,1);
                            line_handles(1,1)=xline(a4,mindiff,'k','LineWidth',1.5);    
                            line_handles(2,1)=xline(a4,maxdiff,'k','LineWidth',1.5); 
                            whichone=0;
                            while whichone~=1 && whichone~=2
                                whichone=input('if mindiff is wrong type 1, if maxdiff is wrong type 2');
                                if isempty(whichone)==1
                                    whichone=0;
                                end
                            end
                            if whichone==2 && min_max_diff(1,2)==0
                                %this is the first time I want to repick
                                %the peak, so just let the program guess
                                %again
                                min_max_diff(1,2)=-1;
                            elseif whichone==1 && min_max_diff(1,1)==0
                                min_max_diff(1,1)=-1;
                            else
                                %not the first time it got it wrong, so the
                                %user needs to do it manually
                                d=datacursormode(f4);
                                d.Enable='on';
                                d.DisplayStyle='window';
                                startxchosen=input('click to record a point and then press enter');
                                if isempty(startxchosen)==1
                                    vals=getCursorInfo(d);
                                    min_max_diff(1,whichone)=vals.Position(1,1); %the point for the correct min or max of yasdiff is recorded in min_max_diff in either col1 or col2 depending on whether it was a peak or valley
                                end
                                d.Enable='off';   
                            end
                            appropriateanswer=1;
                            delete(line_handles)
                            truestart
                            trueend
                        elseif useranswer3==6
                            if known_pks_vals(1,1)==0
                                %this is the first time I want to repick
                                %the diag valley, so just let the program guess
                                %again
                                known_pks_vals(1,1)=-1;
                            else
                                d=datacursormode(f4);
                                 d.Enable='on';
                                 d.DisplayStyle='window';
                                 startxchosen=input('click to record a point and then press enter');
                                if isempty(startxchosen)==1
                                    vals=getCursorInfo(d);
                                    known_pks_vals(1,1)=vals.Position(1,1);
                                    guesscheck(k,2)=-2;
                                end
                                d.Enable='off';   
                                cornercontext(k,2)=useranswer3;
                                cornerindex(k,2)=useranswer3;
                            end
                            appropriateanswer=1;
                            delete(line_handles)
                            truestart
                            trueend
                        elseif isempty(useranswer3)==1 %useranswer3 was empty and therefore the guess was correct
                            if known_pks_vals(1,1)==0
                                guesscheck(k,2)=-1;
                            end
                            cornercontext(k,2)=guesscheck(k,1);
                            cornerindex(k,2)=guesscheck(k,1);
                            appropriateanswer=1;
                            pks_vals_okay(1,1)=1;
                            truestart
                            trueend
                        elseif useranswer3==3
                            %accident, i dont want to analyze this pulse 
                            appropriateanswer=1;
                            guesscheck(k,2)=7;
                            guesscheck(k,4)=7;
                            
                        end
                    end
                    if useranswer3==3 
                        break
                    end
                    if pks_vals_okay(1,1)==0
                        continue
                    end

                    %ask user to check end
                    appropriateanswer=0;
                    while appropriateanswer==0
                        useranswer4=input('if endguess is right- press enter, if you need to repick the "peak" used to calculate the diagonal type 6');
                        if useranswer4==6
                            if known_pks_vals(1,2)==0
                                %this is the first time I want to repick
                                %the diag valley, so just let the program guess
                                %again
                                known_pks_vals(1,2)=-1;
                            else
                                d=datacursormode(f4);
                                 d.Enable='on';
                                 d.DisplayStyle='window';
                                 startxchosen=input('click to record a point and then press enter');
                                if isempty(startxchosen)==1
                                    vals=getCursorInfo(d);
                                    known_pks_vals(1,2)=vals.Position(1,1);
                                    guesscheck(k,4)=-2;
                                end
                                d.Enable='off'; 
                                cornercontext(k,3)=useranswer4;
                                cornerindex(k,3)=useranswer4;
                            end
                            appropriateanswer=1;
                            delete(line_handles)
                            truestart
                            trueend
                        elseif isempty(useranswer4)==1 %useranswer4 was empty and therefore the guess was correct
                            if known_pks_vals(1,1)==0
                                guesscheck(k,4)=-1;
                            end
                            cornercontext(k,3)=guesscheck(k,3);
                            cornerindex(k,3)=guesscheck(k,3);
                            appropriateanswer=1;
                            pks_vals_okay(1,2)=1;
                            truestart
                            trueend
                        end
                    end
                    
                end
                % end pt gathering
                if isempty(useranswer3) || useranswer3~=3
                    k=k+1;
                    hstep=hstepog;
                end
                
                % ANALYZE RECOVERY
                if useranswer0==3 && useranswer2==2
                   %we need to have analyzed the sizing otherwise we cant accurately set the threshold line 
                   % cornercontext(k-1,1)=2, k>2, 
                   %find the mean drop for sizing
                   sizestart=cornercontext(1,2); %i think 1 will work 
                   sizeend=cornercontext(1,3);
                   meandrop=mean(ydetrend(sizestart:sizeend,1));
                   
                   %reregion--give region for thresholding that starts midway through
                   %the squeeze pulse -- then remove that first threshold
                   %crossing thats going up (this way we know we're not
                   %missing the actual first crossing) 
                   sqstart=cornercontext(2,2); %i think 2 will work 
                   sqend=cornercontext(2,3);
                   reregstart=round((sqend-sqstart)/2)+sqstart;
                   %reregionend should be the end of the window 
                   currentsize=xlim(axar);
                   reregend=currentsize(1,2);
                   j=reregstart;
                    threshptsrecov=zeros(2,2); %everytime the signal crosses the threshold line

                    % if the threshold isn't at a good level (not detecting
                    % a recov segment or detecting 2 when it is one, come
                    % back here to reset) 
                     %try 70% of mean drop for threshold 
                      percentofsize=0.7;
                    thresholdgood=0;
                    squeezethresh=0;
                    while thresholdgood==0
                        j=reregstart;
                    threshptsrecov=zeros(2,2); %everytime the signal crosses the threshold line

                        threshold=percentofsize*meandrop; 
                        recovthresh=classicthresholding(yasdet(1:reregend,1),threshold,threshptsrecov(:,1),j,1);
                        if recovthresh(1,2)==1
                            squeezethresh=recovthresh(1,1);
                            recovthresh(1,:)=[];
                        end


    %                     figure
    %                     plot(xaxispost(reregstart:reregend),ydetrend(reregstart:reregend))
    %                     hold on
    %                     plot(xaxispost(reregstart:reregend),yasdet(reregstart:reregend))
    %                     yline(threshold);
    %                     yline(meandrop,'r');
    %                     for j=1:length(recovthresh(:,1))
    %                         xline(recovthresh(j,1))
    %                     end

                        a=reregstart;
                        b=recovthresh(length(recovthresh(:,1)),1)+vstepr;
                        cla(a4)
                        reset(a4)
                        cla(a5)
                        reset(a5)
                        cla(a6)
                        reset(a6)
                            plot(a4,xaxispost(a:b),ydetrend(a:b),'r');
                            plot(a5,xaxispost(a:b),ydetrend(a:b),'r');
                            plot(a6,xaxispost(a-wp:b+wp),ynz(a-wp:b+wp),'r');
                            linkaxes([a4,a5],'x');

                            hold(a4,'on')
                            hold(a5,'on')
                            hold(a6,'on')
                            yline(a5,0,'k');
                            plot(a6,xaxispost(a-wp:b+wp),yasls(a-wp:b+wp),'k')
                            plot(a6,xaxispost(a-wp:b+wp),upperlimitnz(a-wp:b+wp),'c')
                            plot(a6,xaxispost(a-wp:b+wp),expectedporenz(a-wp:b+wp),'k')
                            plot(a6,xaxispost(a-wp:b+wp),yasdetnz(a-wp:b+wp),'b','LineWidth',2)
                            plot(a6,xaxispost(a-wp:b+wp),dinterest_line(a-wp:b+wp),'Color',"#77AC30",'LineWidth',1.5)  
                        %         plot(a6,xaxispost(a-wp:b+wp),maxlimit(a-wp:b+wp)-yasls(a-wp:b+wp),'m')
                        %         plot(a6,xaxispost(a-wp:b+wp),sqlimit(a-wp:b+wp)-yasls(a-wp:b+wp),'g')
                            xlabel(a6,strcat('expected pore, in black, =',num2str(dinput)));
                            title(a4,strcat(' zone number=',num2str(zonenum)));

                            scatter(a4,xaxispost(a:b),ydetrend(a:b),sz,'r');
                            xlabel(a4,xtitle);
                            ylabel(a4,ytitle); 

                            scatter(a5,xaxispost(a:b),yasdiff(a:b).*(10^9),sz,'m')%was *.0001, .*.001
                            plot(a5,xaxispost(a:b),yasdiff(a:b).*(10^9),'m')
                            scatter(a5,xaxispost(a:b),ydiff2(a:b).*(10^8),sz,'k')%was *.0001,
                            plot(a5,xaxispost(a:b),ydiff2(a:b).*(10^8),'k')


                            plot(a4,xaxispost(a:b),yasdet(a:b),'b','LineWidth',2)
                            ylim(a4,[threshold*2 -1*(threshold/2)]);
                            ylim(a5,[threshold*2 -1*threshold]);
                        for j=1:length(recovthresh(:,1))
                            xline(a6,recovthresh(j,1),'b')    
                        end
                        wrongans=1;
                        while wrongans==1
                        newthresh=input(strcat('Current threshold for recov segments is', num2str(percentofsize*100),'% of the sizing pulse, if its good, press enter, if not type the new percent of size'));
                            if isempty(newthresh)==1
                                wrongans=0;
                                thresholdgood=1;
                            elseif newthresh>=1 && newthresh<=100
                                wrongans=0;
                                thresholdgood=0;
                                percentofsize=newthresh/100;
                            else
                                'please input as a percentage eg. 80'
                            end   
                        end
                    end
                    cornercontext(k,1)=3;
                    j=1;    
                    while j<=19
                        %determine cornercontext column
                        col=j+1;
                        if j==1
                            if squeezethresh==0
                                ra=recovthresh(j,1)-vstepr; %vstepr used to be 3 
                            else
                                ra=squeezethresh;
                            end
                                
                        else
                            ra=recovthresh(j-1,1);
                        end
                        if j==19
                            rb=recovthresh(j+1,1)+vstepr; %vstepr used to be 3 
                        else
                            rb=recovthresh(j+2,1);
                        end
                         min_max_diff=[0 0]; 
                                    % col1=min, col2=max, when =0 are correct or haven't been checked
                                    % yet, when=-1 the program has guessed once after being
                                    % wrong, when anything else= user had to manually pick 
                                    pks_vals_okay=[0 0]; %when=0 the end or start wasnt chosen correctly, when=1, it was chosen correctly
                                    known_pks_vals=[0 0];% when =0 are correct or haven't been checked

                                    while pks_vals_okay(1,1)==0 || pks_vals_okay(1,2)==0
                        [truestart,trueend,pulsebottom_mov,line_handles,mindiff,maxdiff,midpk]=findcorners_newway(ydetrend,yasdet,yas3det, ydiff2,yasdiff,xaxispost, vstepr,hstep, mininterval,zonenum,a4,a5,ra,rb,min_max_diff,known_pks_vals);
                        line_handles(17,1)=xline(a6,truestart,'LineWidth',1.5,'Color','#D95319');    
                        line_handles(18,1)=xline(a6,trueend,'LineWidth',1.5,'Color','#D95319');
                              
                         guesscheck(k,1)=truestart;
                         guesscheck(k,3)=trueend;
                        %ask user to check start
                                        appropriateanswer=0;
                                        while appropriateanswer==0
                                            useranswer3=input('if startguess is right- press enter, if you need to fix horizline type 1, if you need to repick the "valley" used to calculate the diagonal type 6, if you decided you dont want to analyze this type 3');

                                            if useranswer3==1 %horizline is wrong 
                                                delete(line_handles)
                                                line_handles=gobjects(2,1);
                                                line_handles(1,1)=xline(a4,mindiff,'k','LineWidth',1.5);    
                                                line_handles(2,1)=xline(a4,maxdiff,'k','LineWidth',1.5); 
                                                whichone=0;
                                                while whichone~=1 && whichone~=2
                                                    whichone=input('if mindiff is wrong type 1, if maxdiff is wrong type 2');
                                                    if isempty(whichone)==1
                                                        whichone=0;
                                                    end
                                                end
                                                if whichone==2 && min_max_diff(1,2)==0
                                                    %this is the first time I want to repick
                                                    %the peak, so just let the program guess
                                                    %again
                                                    min_max_diff(1,2)=-1;
                                                elseif whichone==1 && min_max_diff(1,1)==0
                                                    min_max_diff(1,1)=-1;
                                                else
                                                    %not the first time it got it wrong, so the
                                                    %user needs to do it manually
                                                    d=datacursormode(f4);
                                                    d.Enable='on';
                                                    d.DisplayStyle='window';
                                                    startxchosen=input('click to record a point and then press enter');
                                                    if isempty(startxchosen)==1
                                                        vals=getCursorInfo(d);
                                                        min_max_diff(1,whichone)=vals.Position(1,1); %the point for the correct min or max of yasdiff is recorded in min_max_diff in either col1 or col2 depending on whether it was a peak or valley
                                                    end
                                                    d.Enable='off';   
                                                end
                                                appropriateanswer=1;
                                                delete(line_handles)
                                                truestart
                                                trueend
                                            elseif useranswer3==6
                                                if known_pks_vals(1,1)==0
                                                    %this is the first time I want to repick
                                                    %the diag valley, so just let the program guess
                                                    %again
                                                    known_pks_vals(1,1)=-1;
                                                else
                                                    d=datacursormode(f4);
                                                     d.Enable='on';
                                                     d.DisplayStyle='window';
                                                     startxchosen=input('click to record a point and then press enter');
                                                    if isempty(startxchosen)==1
                                                        vals=getCursorInfo(d);
                                                        known_pks_vals(1,1)=vals.Position(1,1);
                                                        guesscheck(k,2)=-2;
                                                    end
                                                    d.Enable='off';   
                                                    cornercontext(k,col)=useranswer3;
                                                    cornerindex(k,col)=useranswer3;
                                                end
                                                appropriateanswer=1;
                                                delete(line_handles)
                                                truestart
                                                trueend
                                            elseif isempty(useranswer3)==1 %useranswer3 was empty and therefore the guess was correct
                                                if known_pks_vals(1,1)==0
                                                    guesscheck(k,2)=-1;
                                                end
                                                cornercontext(k,col)=guesscheck(k,1);
                                                cornerindex(k,col)=guesscheck(k,1);
                                                appropriateanswer=1;
                                                pks_vals_okay(1,1)=1;
                                                truestart
                                                trueend
                                            elseif useranswer3==3
                                                %accident, i dont want to analyze this pulse 
                                                appropriateanswer=1;
                                                guesscheck(k,2)=7;
                                                guesscheck(k,4)=7;

                                            end
                                        end
                                        if useranswer3==3 
                                            break
                                        end
                                        if pks_vals_okay(1,1)==0
                                            continue
                                        end

                                        %ask user to check end
                                        appropriateanswer=0;
                                        while appropriateanswer==0
                                            useranswer4=input('if endguess is right- press enter, if you need to repick the "peak" used to calculate the diagonal type 6');
                                            if useranswer4==6
                                                if known_pks_vals(1,2)==0
                                                    %this is the first time I want to repick
                                                    %the diag valley, so just let the program guess
                                                    %again
                                                    known_pks_vals(1,2)=-1;
                                                else
                                                    d=datacursormode(f4);
                                                     d.Enable='on';
                                                     d.DisplayStyle='window';
                                                     startxchosen=input('click to record a point and then press enter');
                                                    if isempty(startxchosen)==1
                                                        vals=getCursorInfo(d);
                                                        known_pks_vals(1,2)=vals.Position(1,1);
                                                        guesscheck(k,4)=-2;
                                                    end
                                                    d.Enable='off'; 
                                                    cornercontext(k,col+1)=useranswer4;
                                                    cornerindex(k,col+1)=useranswer4;
                                                end
                                                appropriateanswer=1;
                                                delete(line_handles)
                                                truestart
                                                trueend
                                            elseif isempty(useranswer4)==1 %useranswer4 was empty and therefore the guess was correct
                                                if known_pks_vals(1,1)==0
                                                    guesscheck(k,4)=-1;
                                                end
                                                cornercontext(k,col+1)=guesscheck(k,3);
                                                cornerindex(k,col+1)=guesscheck(k,3);
                                                appropriateanswer=1;
                                                pks_vals_okay(1,2)=1;
                                                truestart
                                                trueend
                                            end
                                        end
                                    end

                        delete(line_handles)
                        j=j+2;
                    end
                end
                
                i=i+2;
            end
    elseif useranswer0==0
        i=i+1;
    elseif useranswer0==1 
        i=i+2;
%     elseif useranswer0==3
%         i=i+(numpul*2)-2; %skips the squeeze and recovery pulses
    elseif useranswer0==4
        i=i+(numpul*2); %skips the entire cell event 
    elseif useranswer0==5
        break
    elseif useranswer0==6
        i=i-4;
    else
        useranswer=0;
    end
    
end

[rguess,~]=size(guesscheck);
cornernums=zeros(1,24);
for i=1:rguess
    if cornercontext(i,1)==1 %is sizing
        %[startindex,endindex]=getindicesmz(cornerindex(i,:),guesscheck(i,:));
        if guesscheck(i,2)~=7
            cornernums(1,1)=guesscheck(i,1);
            cornernums(1,2)=guesscheck(i,3); 
        end
    elseif cornercontext(i,1)==2 %its squeeze
        %[startsqindex,endsqindex]=getindicesmz(cornerindex(i,:),guesscheck(i,:));
        if guesscheck(i,2)~=7
            cornernums(1,3)=guesscheck(i,1);
            cornernums(1,4)=guesscheck(i,3);
        end
    elseif cornercontext(i,1)==3 %itsrecov
        for j=1:20 
        cornernums(1,4+j)=cornercontext(i,j+1);
        end
        
    end
    
end
%done    