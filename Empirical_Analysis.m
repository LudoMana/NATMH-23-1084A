%% Script to analyse empirical EP data

%% clean and data divide it into substages

%%reorganise and clean BOLD

load cnt_ts_scale1.mat 
ts_ctrl=ts_all_scale1_mat;
load ep_ts_scale1.mat % all in folder, alphabethical order 
ts_ep=ts_all_scale1_mat;

%clean from NaNs
ts_ctrl(102,:,:)=[]; 
ts_ctrl(:,[15,48,52,55,56,58,77,110,114,117,118,120,127],:)=[];
ts_ep(:,[15,48,52,55,56,58,77,110,114,117,118,120,127],:)=[];

load STAGES_vec.mat % all in folder, alphabethical order 
%1=stage2;2=stage3a;3=stage3b;4=stage3c;5=stage4; 
stages(stages==999)=NaN;
stages = stages(~(isnan(stages)));

for i=1:5
index_stage{i}=find(stages==i);
Bold_EP_stages{i}=ts_ep(index_stage{i},:,:);
end

Bold_CNT=ts_ctrl;
Bold_EP=ts_ep;
save 'cnt_ts_scale1_clearNaN' Bold_CNT
save 'ep_ts_scale1_clearNaN' Bold_EP
save 'all_STAGES_ts_scale1_clearNaN' Bold_EP_stages

%%reorganise and clean SC

load clean_ep_sc_scale1
sc_ep=connMats;

for i=1:5
index_stage{i}=find(stages==i);
SC_EP_stages{i}=sc_ep(:,:,index_stage{i});
end
save 'all_STAGES_clean_ep_sc_scale1' SC_EP_stages


%clean NaNs EP
for cond=1:5
SC_EP_stages{cond}(:,[15,48,52,55,56,58,77,110,114,117,118,120,127],:)=[];%areas with lots of NaNs
SC_EP_stages{cond}([15,48,52,55,56,58,77,110,114,117,118,120,127],:,:)=[];%areas with lots of NaNs
end
clear connMats

load clean_cnt_sc_scale1.mat
SC_CNT=connMats;

%clean NaNs CNT
SC_CNT(:,:,102)=[]; %only subject missing pericalcarine area
SC_CNT(:,[15,48,52,55,56,58,77,110,114,117,118,120,127],:)=[];%areas with lots of NaNs
SC_CNT([15,48,52,55,56,58,77,110,114,117,118,120,127],:,:)=[];%areas with lots of NaNs


save 'all_STAGES_clean_ep_sc_scale1_clearNaN' SC_EP_stages
save 'clean_cnt_sc_scale1_clearNaN' SC_CNT


%%reorganise and clean SYMPTOMS
load SYMPTOMS_vec.mat % all in folder, alphabethical order
%column PANSSPOS	PANSSNEG	PANSSGEN	PANSSTOTAL GAF

symptoms(symptoms==999)=NaN;
b=isnan(symptoms);
symptoms = symptoms(~(b(:,1)),:);
symptoms(symptoms==0)=NaN;

BOLD_ctrl{1}=(ts_ctrl);
ALL_bold_ts=horzcat(Bold_EP_stages,BOLD_ctrl);

PANS=symptoms(:,1:4);
%PANS = PANS(:,~all(isnan(PANS)));
GAF=symptoms(:,5);
%GAF = GAF(:,~all(isnan(GAF)));

for i=1:5
PANS_EP_stages{i}=PANS(index_stage{i},:);
GAF_EP_stages{i}=GAF(index_stage{i},:);
end

save 'all_stages_SYMPTOMS' PANS_EP_stages GAF_EP_stages

%% filtering + analyisis (FC,dFC, integr/segr, metastability)

%Basic filtering parameters
delta=2; %TR
flp = .04;     % lowpass frequency of filter
fhi = .07;    % highpass
k = 2;                  % 2nd order butterworth filter
fnq = 1/(2*delta);       % Nyquist frequency
Wn = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2] = butter(k,Wn);   % construct the filter

for cond=1:6
    cond
    ts=ALL_bold_ts{cond};
    NSUB(cond) = size(ts,1);
    N = size(ts,2);
    for nsub = 1:NSUB(cond)
    nsub
    xs = double(squeeze(ts(nsub,:,:)));
    FC_emp{cond}(nsub,:,:)=corrcoef(xs'); %%FC matrix
    Tmax = size(xs,2); %%%The time is calculated inside the subjects loop 
                       %%%in case any subject has different time duration.
    T = 1:Tmax;
    timeseriedata = zeros(N,Tmax);%length(xs)
        for seed = 1:N
           x = demean(detrend(xs(seed,:)));
          timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
          Xanalytic = hilbert(demean(timeseriedata(seed,:))); %Hilbert transformation    
          Phases(seed,:) = angle(Xanalytic); %%% calculating phase of each ROI for each signal
                                             %%% which will use to compute the
                                             %%% phase synchronization measures
        end
        
    T = 1:Tmax;
    sync = zeros(1, Tmax);
        for t = T
            ku = sum(complex(cos(Phases(:,t)),sin(Phases(:,t))))/N; %kurmoto param
            ku_all{cond}(nsub,t)=ku;
            sync(t) = abs(ku);
        end
    fluct(nsub) = std(sync(:)); 
    sync_all{cond}(nsub,:)=sync; %syncronization
    %Phase-interaction matrix: %At each time point, the phase difference between two regions was calculated: 
    
     for t = T
      iFC=zeros(N);   
      Isubdiag = find(tril(ones(N),-1));
      for i = 1:N
        for j = 1:N
         iFC(i,j)=cos(Phases(i,t)-Phases(j,t));  %just for iFC trill
         dM(i,j,t) = cos(adif(Phases(i,t),Phases(j,t)));  % computes dynamic matrix/ dFC
         iFCtril = iFC(Isubdiag); 
        end
      end 
      iFCutri(:,t)=iFCtril;
     end
     t
    dM_all_m{cond}(nsub,:,:)=squeeze(mean(dM,3)); %super BIG, saves mean dFC across time points for each sub and cond     
    metastability{cond}(nsub,:)=std(sync); 
    
    niFCutritmax = size(iFCutri,2);
    count        = 1;     % count for dFC
        for t1 = 1:niFCutritmax-2
            p1 = mean(iFCutri(:,t1:t1+2));   % "smoothing" - averaging windows of three timepoints (not eigs but all the iFC conns) (TODO: flexible parameter?)
            for t2 = t1+1:niFCutritmax-2
                p2                     = mean(iFCutri(:,t2:t2+2));       
                fcd_emp{cond}(nsub,count) = dot(p1,p2)/norm(p1)/norm(p2);   % cosine similarity -  product of unit-vectors (spatial direction and magnitudes)
                count                  = count+1; 
            end
        end
        
        %load(sprintf('surrogate_UWS_%02d',nsub))
    %mmdM=squeeze(mean(mdM,1));

    %Integration is calculated: 

    %for t= T
    %cc = mean(dM,3)-mmdM; % Correct with the mean matrices calculated with the surrogates
    cc = mean(dM,3);
    %cc = dM(:,:,t);
    cc = cc-eye(N);
    pp = 1;
    PR = 0:0.01:0.99;
    cs=zeros(1,length(PR));
        for p = PR
        A = (cc)>p;
        [~, csize] = get_components(A);
        cs(pp) = max(csize);
        pp = pp+1;
        end
    integ{cond}(nsub) = sum(cs)*0.01/N; %integration
    end
    
% The modularity (as a measure of segregation) is calculated in the mean matrix and corrected with the 
% bined matrix given by the surrogate and imposing a threhsold of the
% 99% percentile
    meandM=mean(dM,3);
    
end

save 'data_empirical_analysis'

%the important variables are:

%FC_emp = pearson functional connectivity
%dM_all_m = phase coherence functional connectivity, also called dynamical %(dFC) or istantaneous (iFC)

