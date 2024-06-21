#!/bin/bash
#SBATCH --job-name=CNT_omo
#SBATCH --mail-type=END
#SBATCH --array=0-299
#SBATCH --mail-user=ludovica.mana@upf.edu
#SBATCH --output=/home/lmana/EP/hopf/Job%A_%a.out
#SBATCH --error=/home/lmana/EP/hopf/Job%A_%a.err
#SBATCH --mem-per-cpu=3000

ml MATLAB
matlab -singleCompThread -nojvm -nodisplay<<-EOF

iD = str2double(getenv('SLURM_ARRAY_TASK_ID'));
iD = iD + 1;


%% Homogeneous model
%This model based on Stuard-Landau oscillators, the local dynamics of single 
%brain regions using the normal form of a Hopf bifurcation.
%The dynamics of the N brain regions were coupled through the connectivity 
%matrix, which was given by the connectome of healthy subjects (C).
%coupling among areas is given by the SC-> by the g scaling parameter.
%For the homogegenous case,  local dynamical properties are fixed to be at
%the bifurcation point (a=0)

%This version of the model is prepared for group level, taking into account 
%all the subjects (the SC is a template).

%% Load data and basic variables

load('cnt_ts_scale1_clearNaN.mat') %timeseries of each group Nsubjects*nodes*time
ts=ts_ctrl;
NSUB = size(ts,1);
Nodes = size(ts,2);

load('clean_cnt_sc_scale1')
sc_healthy_cnt=connMats;

sc_healthy_cnt(:,:,102)=[]; 
sc_healthy_cnt(:,[15,48,52,55,56,58,77,110,114,117,118,120,127],:)=[];%areas with lots of NaNs
sc_healthy_cnt([15,48,52,55,56,58,77,110,114,117,118,120,127],:,:)=[];%areas with lots of NaNs

sc_healthy=mean(sc_healthy_cnt,3);

Tmax=size(ts,3);
wG = 0:0.05:2.5;%this is de g vector i.e,           
ldata = NSUB; %Number of subjets.
nn=1;


totalH2=symmetrization_ts_lau_cNaN(ts,Nodes,ldata);
C_old=squeeze((sc_healthy(1:Nodes,1:Nodes)));
C=symmetrization_sc_lau_cNaN(C_old,Nodes);
tseries=totalH2;


%% Details of the model

Cfg.simulID = 'Psychosis';
Cfg.opt_a.nIters = 200;%Number of iterations for initial optimization
Cfg.opt_a.updateStrength = 0.1;%Update strength for a optimization
Cfg.opt_a.abortCrit = 0.1; %maximally 2.5% error
Cfg.opt_a.gref = 2;%This is reference G, for the initial optimization only. 
                   % It should be close to the optimal G obtained by the
                   % homogeneous model
Cfg.TRsec = 1.9;
if ~isfield(Cfg, 'simulID'), Cfg.simulID = 'Unknown'; else end;

if ~isfield(Cfg, 'TRsec'), Cfg.TRsec = 2; else end;

if ~isfield(Cfg, 'opt_a'), Cfg.opt_a = []; else end;
if ~isfield(Cfg.opt_a, 'nIters'), Cfg.opt_a.nIters = 100; else end;
if ~isfield(Cfg.opt_a, 'updateStrength'), Cfg.opt_a.updateStrength = 0.1; else end
if ~isfield(Cfg.opt_a, 'abortCrit'), Cfg.opt_a.abortCrit = 0.1; else end
if ~isfield(Cfg, 'plots'), Cfg.plots = []; else end;
if ~isfield(Cfg.plots, 'showOptimization'), Cfg.plots.showOptimization = 0; else end;
if ~isfield(Cfg.plots, 'makeVideo'), Cfg.plots.makeVideo = 0; else end;


rng('shuffle');
nNodes = length(C);
nSubs = ldata; %Number of subjects
si = 1:ldata;
nWeights=numel(wG); % In the case that we only calculated for the optimal g obtaiend in the homo model.
               %If vary g -> nWeights=numel(wG);

fprintf(1, 'Fitting models for %d subjects and %d different weights\n', nSubs, nWeights);

FC_simul = zeros(nNodes, nNodes, nWeights);
fitting = zeros(1, nWeights);
meta = zeros(1, nWeights);
ksP = zeros(1, nWeights);
Phases = zeros(nNodes, Tmax, nSubs, nWeights);
bifpar = zeros(nWeights, nNodes);

%--------------------------------------------------------------------------
%CALCULATE FUNCTIONAL CONNECTIVITY MATRIX
%--------------------------------------------------------------------------
r = zeros(nNodes, nNodes, nSubs);
ts = zeros(nNodes, Tmax, nSubs);

for i = 1:nSubs;
    ts(:,:,i) = tseries{si(i)};
    r(:,:,i) = corrcoef(ts(:,:,i)');
end

FC_emp=nanmean(r,3);
C=C/max(max(C))*0.2;% Scale the SC matrix

%--------------------------------------------------------------------------
%COMPUTE POWER SPECTRA FOR
%NARROWLY FILTERED DATA WITH LOW BANDPASS (0.04 to 0.07 Hz)
%WIDELY FILTERED DATA (0.04 Hz to justBelowNyquistFrequency)
%--------------------------------------------------------------------------
TT=Tmax;
Ts = TT*Cfg.TRsec;
freq = (0:TT/2-1)/Ts;
[~, idxMinFreq] = min(abs(freq-0.04));
[~, idxMaxFreq] = min(abs(freq-0.07));
nFreqs = length(freq);

delt = 1.9;                                   % sampling interval
fnq = 1/(2*delt);                           % Nyquist frequency
k = 2;                                      % 2nd order butterworth filter

%WIDE BANDPASS
flp = .04;                                  % lowpass frequency of filter
fhi = fnq-0.001;%.249;                      % highpass needs to be limited by Nyquist frequency, which in turn depends on TR
Wn = [flp/fnq fhi/fnq];                     % butterworth bandpass non-dimensional frequency
[bfilt_wide, afilt_wide] = butter(k,Wn);    % construct the filter

%NARROW LOW BANDPASS
flp = .04;                                  % lowpass frequency of filter
fhi = .07;                                  % highpass
Wn=[flp/fnq fhi/fnq];                       % butterworth bandpass non-dimensional frequency
[bfilt_narrow,afilt_narrow] = butter(k,Wn); % construct the filter


for seed=1:Nodes
    
    for idxSub=1:nSubs %[1:96 98:99 101:nSubs]
        signaldata = double(tseries{si(idxSub)});
        x=detrend(demean(signaldata(seed,:)));
        
        ts_filt_narrow =zscore(filtfilt(bfilt_narrow,afilt_narrow,x));
        pw_filt_narrow = abs(fft(ts_filt_narrow));
        PowSpect_filt_narrow(:,seed,idxSub) = pw_filt_narrow(1:floor(TT/2)).^2/(TT/2);
        
        ts_filt_wide =zscore(filtfilt(bfilt_wide,afilt_wide,x));
        pw_filt_wide = abs(fft(ts_filt_wide));
        PowSpect_filt_wide(:,seed,idxSub) = pw_filt_wide(1:floor(TT/2)).^2/(TT/2);
    end
    
end

Power_Areas_filt_narrow_unsmoothed = nanmean(PowSpect_filt_narrow,3);
Power_Areas_filt_wide_unsmoothed = nanmean(PowSpect_filt_wide,3);
Power_Areas_filt_narrow_smoothed = zeros(nFreqs, nNodes);
Power_Areas_filt_wide_smoothed = zeros(nFreqs, nNodes);
vsig = zeros(1, nNodes);

for seed=1:Nodes
    Power_Areas_filt_narrow_smoothed(:,seed)=gaussfilt(freq,Power_Areas_filt_narrow_unsmoothed(:,seed)',0.01);
    Power_Areas_filt_wide_smoothed(:,seed)=gaussfilt(freq,Power_Areas_filt_wide_unsmoothed(:,seed)',0.01);
    %relative power in frequencies of interest (.04 - .07 Hz) with respect
    %to entire power of bandpass-filtered data (.04 - just_below_nyquist)
    vsig(seed) =...
        sum(Power_Areas_filt_wide_smoothed(idxMinFreq:idxMaxFreq,seed))/sum(Power_Areas_filt_wide_smoothed(:,seed));
end

vmax=max(vsig);
vmin=min(vsig);

%a-minimization seems to only work if we use the indices for frequency of
%maximal power from the narrowband-smoothed data
[~, idxFreqOfMaxPwr]=max(Power_Areas_filt_narrow_smoothed);
f_diff = freq(idxFreqOfMaxPwr);

%FOR EACH AREA AND TIMEPOINT COMPUTE THE INSTANTANEOUS PHASE IN THE RANGE
%OF .04 TO .07 Hz
PhasesD = zeros(nNodes, Tmax, nSubs);
for idxSub=1:nSubs %[1:96 98:99 101:nSubs]
    signaldata=double(tseries{si(idxSub)});
    for seed=1:Nodes
        x = demean(detrend(signaldata(seed,:)));
        xFilt = filtfilt(bfilt_narrow,afilt_narrow,x);    % zero phase filter the data
        Xanalytic = hilbert(demean(xFilt));
        PhasesD(seed,:,idxSub) = angle(Xanalytic);
    end 
end

%f_diff  previously computed frequency with maximal power (of narrowly filtered data) by area
omega = repmat(2*pi*f_diff',1,2); %angular velocity of each oscillator
omega(:,1) = -omega(:,1);


%% FROM HERE ON SIMULATIONS AND FITTING


%1- Optimization algorithm for obtaining approximated bifurcation parameters
%before running the model with the g. Here we optimize for g ref (Cfg.opt_a.gref) 

dt = 0.1;
sig = 0.04;
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step
a = repmat(-0.00005*ones(Nodes,1),1,2); %This vector contains the starting bif par. 
wC = Cfg.opt_a.gref*C;% This is for starting, with a fixed g.
sumC = repmat(sum(wC,2),1,2);
trackminm1 = zeros(Cfg.opt_a.nIters, nWeights); %for tracking the minimization (good for debugging)

minm=100;
bestIter = 1; 
t1=0;
current_diff = 1000000;
iter_thres=100000;
learning_rate=0.1;
diff_thres=0.01;
Cfg.opt_a.abortCrit=2*std(vsig); %This value depends on the data and the amplitude. 
                         %It highly depends on the preprocessing of the
                         %data. 
                         %can be standarized as Cfg.opt_a.abortCrit=mean(visg)+2*std(vsig)
bestIter = 1; 

%Not needed in the homo case
% for iter = 1:100
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%MODEL%%%%%%%%%%%%%%%%%
%     z = 0.1*ones(nNodes,2); % --> x = z(:,1), y = z(:,2)
%     nn=0;
%     k=1;
%     for t=1:dt:1000
%         suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
%         zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
%         z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
%         z_all(k,:)=z(:,1);
%         k=k+1;
%     end
% 
%     for t=1:dt:Tmax*Cfg.TRsec
%         suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
%         zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
%         z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
%           z_all(k,:)=z(:,1);
%         k=k+1;
%         if mod(t,Cfg.TRsec)==0
%             nn=nn+1;
%             xs(nn,:)=z(:,1)';
%         end
%     end
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%
%     %Calculate spectral properties of the simulated data
%     
%     vsigs = zeros(1, Nodes);
%     for seed=1:Nodes
% 
%         x=detrend(demean(xs(1:nn,seed)'));
%         ts_filt_wide =zscore(filtfilt(bfilt_wide,afilt_wide,x));
% 
%         TT=length(x);
%         Ts = TT*Cfg.TRsec;
%         freq = (0:TT/2-1)/Ts;
%         [~, idxMinFreqS]=min(abs(freq-0.04));
%         [~, idxMaxFreqS]=min(abs(freq-0.07));
% 
%         pw_filt_wide = abs(fft(ts_filt_wide));
%         Pow1 = pw_filt_wide(1:floor(TT/2)).^2/(TT/2);
%         Pow=gaussfilt(freq,Pow1,0.01);
% 
%         vsigs(seed)=sum(Pow(idxMinFreqS:idxMaxFreqS))/sum(Pow);
% 
%     end
% 
%     %Compare the spectral properties of the simulated and empirical data
%     %and depending on teh comparison update the bif par vector.
%     
%     vsmin=min(vsigs);
%     vsmax=max(vsigs);
%     bb=(vmax-vmin)/(vsmax-vsmin);
%     aa=vmin-bb*vsmin;  %% adaptation of local bif parameters
%     vsigs=aa+bb*vsigs;
%     minm1=max(abs(vsig-vsigs)./vsig);
%     trackminm1(iter, Cfg.opt_a.gref) = minm1;%this is for tracking of minm1
%     vsig_all(iter,:)=vsigs; %For tracking the spectral properties
% 
%     if minm1<minm
%         minm=minm1;
%         a1=a;
%         bestIter = iter;
%         best_vsigs = vsigs; 
%     end
% 
%     %--------------------------------------------------------------
%     %FEEDBACK
%     %--------------------------------------------------------------
%     if Cfg.plots.showOptimization
%         showOptimPlot(h_track_opt, idx_g, we, iter, a, a1, vsig, vsigs, bestIter, best_vsigs, trackminm1, Cfg)
%     end
%     fprintf(1, 'iter: %03d, minm1: %5.3f\n', iter, minm1);
%     %--------------------------------------------------------------
% 
%     %CRITERION REACHED?
%     if minm<Cfg.opt_a.abortCrit 
%         break;
%     end
% 
%     %UPDATE a VALUES FOR NEXT ITER
%     fun = vsigs-vsig;
%     if ~any(isnan(vsigs))
%         a(:,1)=a(:,1)+0.2*(1-vsigs./vsig)';
%         a(:,2)=a(:,2)+0.2*(1-vsigs./vsig)';
%         a_all(iter,:,:)=a;
%     else
%         warning('There are NaNs in the power spectra. Probably a-values too strong.');
%         a = a1;
%     end
% end
% 
% a=a1; %Bif par vector for the optimal step
% 
% 
% %Linearize for all gs 
% aLinear=aLin(a(:,1)',Cfg.opt_a.gref,wG);


%SIMULATE FOR EACH g USING LINEARIZED A VALUES:
%In this case we only simulate for the optimal g obtained from the
%homogeneous model so the second for loop is not neccessary but I leave it
%in case it is useful in other cases.

trial=iD
for idx_g = 1:nWeights
we = wG(idx_g); %In case defining a wG with values of g%
%load('optimal_g')
%we=g;

wC = we*C;
xs = zeros(3000/2,nNodes);
sumC = repmat(sum(wC,2),1,2);

fprintf(1, '-----------------------------------------\n');
fprintf(1, 'G(%d/%d) = %5.3f\n', idx_g, numel(wG), we);
fprintf(1, '-----------------------------------------\n');

av=aLinear(idx_g,:);
a=repmat(av',1,2);

minm=100;
bestIter = 1; 

%Optimize again but now for the optimal g and based on the traning of the
%bif par obtained before.

% for iter = 1:100
%     z = 0.1*ones(nNodes,2); % --> x = z(:,1), y = z(:,2)
%     nn=0;
%     for t=1:dt:1000
%         suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
%         zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
%         z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
%     end
% 
%     for t=1:dt:Tmax*Cfg.TRsec
%         suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
%         zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
%         z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
% 
%         if mod(t,2)==0
%             nn=nn+1;
%             xs(nn,:)=z(:,1)';
%         end
%     end
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%
%     vsigs = zeros(1, nNodes);
%     for seed=1:nNodes
% 
%         x=detrend(demean(xs(1:nn,seed)'));
%         ts_filt_wide =zscore(filtfilt(bfilt_wide,afilt_wide,x));
% 
%         TT=length(x);
%         Ts = TT*Cfg.TRsec;
%         freq = (0:TT/2-1)/Ts;
%         [~, idxMinFreqS]=min(abs(freq-0.04));
%         [~, idxMaxFreqS]=min(abs(freq-0.07));
% 
%         pw_filt_wide = abs(fft(ts_filt_wide));
%         Pow1 = pw_filt_wide(1:floor(TT/2)).^2/(TT/2);
%         Pow=gaussfilt(freq,Pow1,0.01);
% 
%         vsigs(seed)=sum(Pow(idxMinFreqS:idxMaxFreqS))/sum(Pow);
% 
%     end
% 
%     vsmin=min(vsigs);
%     vsmax=max(vsigs);
%     bb=(vmax-vmin)/(vsmax-vsmin);
%     aa=vmin-bb*vsmin;  %% adaptation of local bif parameters
%     vsigs=aa+bb*vsigs;
%     minm1=max(abs(vsig-vsigs)./vsig);
%     trackminm1(iter, idx_g) = minm1;
% 
%     if minm1<minm
%         minm=minm1;
%         a1=a;
%         bestIter = iter; 
%         best_vsigs = vsigs;
%     end
% 
%     %--------------------------------------------------------------
%     %FEEDBACK
%     %--------------------------------------------------------------
%     if Cfg.plots.showOptimization
%         showOptimPlot(h_track_opt, idx_g, we, iter, a, a1, vsig, vsigs, bestIter, best_vsigs, trackminm1, Cfg)
%     end
%     fprintf(1, 'iter: %03d, minm1: %5.3f\n', iter, minm1);
%     %--------------------------------------------------------------
% 
%     %CRITERION REACHED?
%     if minm<Cfg.opt_a.abortCrit %default is 0.1
%         break;
%     end
% 
%     %UPDATE a VALUES FOR NEXT ITER
%     if ~any(isnan(vsigs))
%         a(:,1)=a(:,1)+0.2*(1-vsigs./vsig)';
%         a(:,2)=a(:,2)+0.2*(1-vsigs./vsig)';
%     else
%         warning('There are NaNs in the power spectra. Probably a-values too strong.');
%         a = a1;
%     end
% end

a=zeros(115,2);  % a fixed at bifurcation point a=0 for all areas
bifpar(idx_g,:)=a(:,1)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Final simulation with every g and all zeros vector of a.
fprintf(1, 'SIMULATING OPTIMIZED MODEL.\n');

    xs=zeros(Tmax*1,nNodes);
    z = 0.1*ones(nNodes,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    for t=1:dt:3000 
        suma = wC*z- sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
    end

    for t=1:dt:Tmax*1.9
        suma = wC*z- sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);

        if mod(t,1.9)==0
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
fprintf(1, 'COMPUTING MODEL FIT.\n');

FC_simul(:, :, idx_g) = corrcoef(xs(1:nn,:)); 
cc=corrcoef(squareform(tril(FC_emp,-1)),squareform(tril(FC_simul(:, :, idx_g),-1)));%atanh(FC...
fitting(idx_g)=cc(2);

%%%%%%%%%%%%%%%%%%

for seed=1:nNodes
    ts_simul = detrend(demean(xs(1:Tmax,seed)'));
    ts_simul_filt_narrow = filtfilt(bfilt_narrow,afilt_narrow,ts_simul);
    Xanalytic = hilbert(ts_simul_filt_narrow);
    Phases(seed,:,1, idx_g) = angle(Xanalytic);
end

xs_all(idx_g,:,:)=xs;
cc2=corrcoef(squareform(tril(C,-1)),squareform(tril(FC_simul(:, :, idx_g),-1)));
SC_FC_g(idx_g)=cc2(2);

%Fitting of the FCD by means of the Kolmogorov-Smirnov distance:
%Compare to all subjects
k=1;
pcS(:,1)=patternCons30(Phases(:,:,1, idx_g),nNodes,Tmax);
load('pcD_cnt')
% for i = 1:nSubs %[1:96 98:99 101:nSubs]
%      pcD(:,i)=patternCons30(PhasesD(:,:,i),nNodes,Tmax);
%     k=k+1
% end

 [~, ~, ksP(idx_g)]=kstest2(pcS(:),pcD(:));

fprintf(1, 'DONE.\n\n');
end 

save(sprintf('model_homo_cnt_trial_%02d',trial),'ksP','pcS','meta','bifpar','FC_simul','xs_all')
EOF
