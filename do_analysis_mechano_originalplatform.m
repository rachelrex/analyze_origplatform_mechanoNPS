
%this code runs your analysis script 

if exist('dmin','var') ~=1
    dmin=input('dmin is not saved in your workspace, type what you want dmin to be right now');
end
if exist('plot_thresh','var') ~=1
    plot_thresh=0;
end

% run if no figures are up (automatically starts from the beginning) 

[pulses,threshpts,fig_handles]=oneznps_pickpulse_semiauto_nodownsamp(ydetrend,-1*yasls,yas2det,dmin,stderror,yas3det,plot_thresh,0,0,0);

%% run if your figures are NOT already up- but you want to start from where you left off 
if exist('dmin','var') ~=1
    dmin=input('dmin is not saved in your workspace, type what you want dmin to be right now');
end
if exist('plot_thresh','var') ~=1
    plot_thresh=0;
end
%if you copied and pasted the pulses as a new variable there needs to be at
%least one row of zeros at the end - so lets add one
numpulses=length(pulses(:,1));
pulses(numpulses+1,:)=zeros(1,length(pulses(1,:)));
numpulses=numpulses+1;
i=2;
while i<=numpulses
    if pulses(i,1)==0
        i=i-1;
        skipto=pulses(i,1);
        break
    else
        i=i+1;
    end
end
[pulses,threshpts,fig_handles]=oneznps_pickpulse_semiauto_nodownsamp(-1*ydetrend,yasls,-1*yas2det,dmin,stderror,-1*yas3det,plot_thresh,skipto,0,0);


%% run if your figures are already up and you want to start from the beginning 
use_old_fig=1;
skipto=0;
if exist('plot_thresh','var') ~=1
    plot_thresh=0;
end


[pulses,threshpts,fig_handles]=oneznps_pickpulse_semiauto_nodownsamp(ydetrend,-1*yasls,yas2det,dmin,stderror,yas3det,plot_thresh,skipto,use_old_fig,fig_handles);

%% run if your figures are already up and you want to start from where you left off
use_old_fig=1;
i=2;
while i<=length(pulses(:,1))
    if pulses(i,1)==0
        i=i-1;
        skipto=pulses(i,1);
        break
    else
        i=i+1;
    end
end

[pulses,threshpts,fig_handles]=oneznps_pickpulse_semiauto_nodownsamp(ydetrend,-1*yasls,yas2det,dmin,stderror,yas3det,plot_thresh,skipto,use_old_fig,fig_handles);


%% find fig handles if figures are up, but handles aren't saved
%this might happen if it was your first run of the session and it ended in an error 

fa=findobj('type','figure','Name','Data Set');
fb=findobj('type','figure','Name','Entire Event');
f4=findobj('type','figure','Name','Window of Interest');
fig_handles=[fa;fb;f4];
