%% 1. Rename + Reformat Data
%downsample to match multizone data 
%original samplerate is 50,000, mz is 10,000, downsample by 5
%downsampling is essential for performance of ASLS
datanew=-data';
R=downsample(datanew,5); %this is a current signal, but since the rest 
% of this code was built for resistance signals in multizone
% i'm going to keep the resistance notation

%% OPTIONAL raw data plotting
figure
plot(R)

%% 2. Filter Resistance Signals --------------------------------------------
fpass = 50; % passband frequency of the lowpass filter (in Hz)
window = 31; % averaging window width (in units of samples) 
window2= 71;
polynomial=3;
polynomial2=1;

numsgo=10;
Rfalt1= sgolayfilt(lowpass(R, fpass, sampleRate),polynomial,window);  
Rfalt2= sgolayfilt(lowpass(R, fpass, sampleRate),polynomial,window2);  
Rfalt3= sgolayfilt(lowpass(R, fpass, sampleRate),polynomial2,window2);  

for i=numsgo
    Rfalt1= sgolayfilt(Rfalt1,polynomial,window); 
    Rfalt2= sgolayfilt(Rfalt2,polynomial,window2); 
    Rfalt3= sgolayfilt(Rfalt3,polynomial2,window2); 
end
Rlow=lowpass(R,50,sampleRate);
Rlow=Rlow(window2:end-window2,:);

Rfalt1=Rfalt1(window2:end-window2,:);
Rfalt2=Rfalt2(window2:end-window2,:);
Rfalt3=Rfalt3(window2:end-window2,:);

f1=figure;
hold on

plot(Rlow(:,1),'LineWidth',2.0)
plot(Rfalt1(:,1),'LineWidth',2.0)
plot(Rfalt2(:,1),'LineWidth',2.0)

Rf=Rfalt1;

%% OPTIONAL plotting for visualization



f1=figure;
hold off
 plot(Rlow(:,1),'LineWidth',2.0)
hold on
plot(Rfalt1(:,1),'LineWidth',2.0)
plot(Rfalt2(:,1),'LineWidth',2.0)
plot(Rfalt3(:,1),'LineWidth',2.0)



Rf=Rfalt1;


%% 3. Get stderror (baseline noise)

 %pick region for stdev calculation 
 stderror=zeros(2,1);
   
        d=datacursormode(f1);
          d.Enable='on';
         d.DisplayStyle='window';
         startstdi=input('click where you want region to start, then press enter');
         if isempty(startstdi)==1
             vals=getCursorInfo(d);
             startstd=vals.Position(1,1);
         end
         d.Enable='off';
         d.Enable='on';
         endstdi=input('click where you want region to end, then press enter');
         if isempty(endstdi)==1
             vals=getCursorInfo(d);
             endstd=vals.Position(1,1);
         end
         
            stderror(1,1)=std(Rf(startstd:endstd,1));
        

%% 4. and 8. Fit baseline
%for asls if not downsampled 
lambda1 = 1e12; %larger means smoother background 
p=0; % less than 0.5 means negative is strongly supressed %0.02
maxiter=20; % maximum iterations 
noisez1=stderror(1,1);


aslsparam=struct('lambda', lambda1, 'p', p, 'max_iter', maxiter, 'noise_margin', noisez1);

       %yasls=-1*(ASLS2(Rf*-1,aslsparam));
       yasls=ASLS2(Rf,aslsparam);

%% OPTIONAL (recommended) plot yasls results (check asls fit) 
figure

plot(yasls)




    hold on
   
    plot(Rf)
    plot(Rfalt3) 


    
 %% 5. and 9. subtract baseline_replot
f1=figure;
hold on

    ydetrend=Rf-yasls;
    yas2det=Rfalt2-yasls;
    yas3det=Rfalt3-yasls;
    plot(ydetrend)
    plot(yas2det)
    plot(yas3det,'LineWidth',2)

%stop here at step 9. and save your pre-processed data 

%% 6. calculate noise again- use same region as before
%pick region for stdev calculation 

     
            stderror(2,1)=std(ydetrend(startstd:endstd,1));
   
        
         stderror


%% 6.b.(optional instead of 6.) calculate noise again - pick new region 
%pick region for stdev calculation 

       d=datacursormode(f1);
          d.Enable='on';
         d.DisplayStyle='window';
         startstdi=input('click where you want region to start, then press enter');
         if isempty(startstdi)==1
             vals=getCursorInfo(d);
             startstd=vals.Position(1,1);
         end
         d.Enable='off';
        

         d.Enable='on';
         endstdi=input('click where you want region to end, then press enter');
         if isempty(endstdi)==1
             vals=getCursorInfo(d);
             endstd=vals.Position(1,1);
         end
          std4errorplot=zeros(endstd-startstd+1,3);
       
            stderror(2,1)=std(ydetrend(startstd:endstd,1));
            std4errorplot(:,1)=ydetrend(startstd:endstd,1);
    



%% 7. reset stderror

stderror(1,1)=stderror(2,1)


   

