% here we compare values of iFC (or dFC or phase coherence FC)in the different suboups of interest
% we reorganise the data 
    % origininal oranization : 1=stage2; 2=stage3a; 3=stage3b; 4=stage3c;  5=ctrl;
    %actual oranisation: REM=stage 2, REM_REL=stage 3 REMITTINg(stage 3b+3c), nonREM=stage 3 non-REMITTINg (stage3a),  ctrl
    %NB: at the beinning I used diff labels than the actuals: REL=relapsing, REM=remitting 
%we compare global and local measures of iFC:
    %global iFC, delta iFC matrix, strength, diversity, mean delta strength 
    
load data_empirical_analysis.mat     
%% global measures

%global iFC 

for cond=1:5
    for sub=1:NSUB(cond)
        global_FC{cond}(sub)=squeeze(mean(FC_emp{cond}(sub,:,:),'all'));
    end
end

REM_gFC=global_FC{1} %stage 2
REM_REL_gFC=horzcat(global_FC{3},global_FC{4}); %stage 3 REMITTINg
nonREM_gFC=horzcat(global_FC{2}); %stage 3 nonREMITTINg

mean_CNT=mean(global_FC{5});err_CNT=std(global_FC{5})/size(global_FC{5},2);
mean_REM_gFC=mean(REM_gFC);err_REM_gFC=std(REM_gFC)/size(REM_gFC,2);
mean_REM_REL_gFC=mean(REM_REL_gFC);err_REM_REL_gFC=std(REM_REL_gFC)/size(REM_REL_gFC,2);;
mean_nonREM_gFC=mean(nonREM_gFC);err_nonREM_gFC=std(nonREM_gFC)/size(nonREM_gFC,2);

all_EP_gFC=horzcat(global_FC{1},global_FC{2},global_FC{3},global_FC{4});

condA=global_FC{5}; %ctrl 
condB=REM_REL_gFC; 
condC=nonREM_gFC;
condD=all_EP_gFC;
condE=REM_gFC;

[p_gFC(1),H]=ranksum(condA,condB); %CTRL vs REM_REL
[p_gFC(2),H]=ranksum(condA,condC);%CTRL vs nonREM
[p_gFC(3),H]=ranksum(condB,condC); %REM vs nonREM

[p_gFC(6),H]=ranksum(condA,condD);%CTRL vs EP2
[p_gFC(5),H]=ranksum(condA,condC); %CTRL vs EP3
[p_gFC(4),H]=ranksum(condA,condD);%CTRL vs allEP

figure()
C=[condA condB condC] %change cond to compare
grp = [zeros(1,size(condA,2)),ones(1,size(condB,2)),ones(1,size(condC,2))*2];
boxplot(C,grp)
yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
set(gca, 'Xtick', 1:3);
xt = get(gca, 'XTick');
hold on
    if p_gFC(1)<(0.05) %p_gFC<(0.05/n_mult_comp)
        plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k')
    end
    if p_gFC(2)<(0.05) %p_gFC<(0.05/n_mult_comp)
        plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k')
    end
    if p_gFC(3)<(0.05) %p_gFC<(0.05/n_mult_comp)
        plot(xt([1 3]), [1 1]*max(yt)*1.3, '-k',  mean(xt([2 3])), max(yt)*1.35, '*k')
    end
  
title('Global FC')

%% LOCAL MEASAURES: delta FC matrix  
% S = readtable('AREAS_1.xlsx');
% StructNames=table2array(S);
% StructNames=char(StructNames);

FC_REM=FC_emp{1}; % stage 2
FC_REM_REL=vertcat(FC_emp{3},FC_emp{4}); %stage 3 REM
FC_nonREM=vertcat(FC_emp{2}); %%stage 3 nonREM

FC_mean_ctrl=squeeze(mean(FC_emp{5})); %mean over subjects to do the diff
FC_mean_REM=squeeze(mean(FC_REM));
FC_mean_REM_REL=squeeze(mean(FC_REM_REL));
FC_mean_nonREM=squeeze(mean(FC_nonREM));

% delta FC matrix (EP-ctrl)

diff_FC1=FC_mean_REM-FC_mean_ctrl;
diff_FC2=FC_mean_REM_REL-FC_mean_ctrl;
diff_FC3=FC_mean_nonREM-FC_mean_ctrl;

diff_FC_masked1=diff_FC1;
diff_FC_masked2=diff_FC2;
diff_FC_masked3=diff_FC3;

for n=1:N
   for p=1:N
       a1=squeeze(FC_emp{5}(:,n,p));
       b1=squeeze(FC_REM(:,n,p));
       [p_FC1,H]=ranksum(a1,b1);
       b2=squeeze(FC_REM_REL(:,n,p));
       b3=squeeze(FC_nonREM(:,n,p));
       [p_FC2,H]=ranksum(a1,b2);
       [p_FC3,H]=ranksum(a1,b3);
       if p_FC1>=0.01
           diff_FC_masked1(n,p)=NaN; %mask to show only values significantly diff between cond
       end
        if p_FC2>=0.01
           diff_FC_masked2(n,p)=NaN;
        end
       if p_FC3>=0.01
           diff_FC_masked3(n,p)=NaN;
       end
    P_FC1(n,p)=p_FC1;
    P_FC2(n,p)=p_FC2;
    P_FC3(n,p)=p_FC3;
   end
end

figure();  %UNCOMMENT LINES TO PLOT COND OF INTEREST 
%li=max(diff_FC_masked1,[],'all');
%li=max(diff_FC_masked2,[],'all');%change
li=max(diff_FC_masked3,[],'all');%change
CLIM=[-0.12,li];
%f=imagesc(diff_FC_masked1); %f=imagesc(diff_FC_masked1,CLIM);
f=imagesc(diff_FC_masked2);
%f=imagesc(diff_FC_masked3);
%set(f,'AlphaData', ~isnan(diff_FC_masked1))
set(f,'AlphaData', ~isnan(diff_FC_masked2))
% set(f,'AlphaData', ~isnan(diff_FC_masked3))
axis square
%yticklabels(StructNames);
set(gca, 'Xtick',0:10:115,'fontsize',24);
set(gca, 'Ytick',0:10:115,'fontsize',24);
xlabel('Brain areas', 'FontSize',28); ylabel('Brain areas','FontSize',28);
%title('Group2 - Ctrl','FontSize',30)
title('Group3REM_REL - Ctrl','FontSize',30)
%title('Group3nonREM - Ctrl','FontSize',30)

% MEAN DIFFERENCE of FC upptr (mean over differences between paired areas in two cond: global diff,EP-ctrl)

A=triu(diff_FC1); vect_FC_diff1_R=A(:);
vect_FC_diff1_R=vect_FC_diff1_R(~vect_FC_diff1_R==0)';
B=triu(diff_FC2); vect_FC_diff2_RR=B(:);
vect_FC_diff2_RR=vect_FC_diff2_RR(~vect_FC_diff2_RR==0)';
C=triu(diff_FC3); vect_FC_diff3_nR=C(:);
vect_FC_diff3_nR=vect_FC_diff3_nR(~vect_FC_diff3_nR==0)';

Rmean=mean(vect_FC_diff1_R);
RRmean=mean(vect_FC_diff2_RR);
nRmean=mean(vect_FC_diff3_nR);
R_std=std(vect_FC_diff1_R);
RR_std=std(vect_FC_diff2_RR);
nR_std=std(vect_FC_diff3_nR);
err_R=R_std/size(vect_FC_diff1_R,2);
err_RR=RR_std/size(vect_FC_diff2_RR,2);
err_nR=nR_std/size(vect_FC_diff3_nR,2);
conf_99R=5.15.*err_R;
conf_99RR=5.15.*err_RR;
conf_99nR=5.15.*err_nR;

cohend_R=abs(Rmean)/R_std;
cohend_RR=abs(nRmean)/RR_std;
cohend_nR=abs(nRmean)/nR_std;

H0=zeros(1,size(vect_FC_diff1_R,2));

[P_diff_pair_R(1),~]=ranksum(vect_FC_diff1_R,H0);
[P_diff_pair_R(2),~]=ranksum(vect_FC_diff2_RR,H0);
[P_diff_pair_R(3),~]=ranksum(vect_FC_diff3_nR,H0);

bonf_P_diff_R=P_diff_pair_R/6786;
err=[conf_99R,conf_99RR,conf_99nR];
err=[conf_99R,conf_99RR,conf_99nR];

figure()
C=[Rmean,RRmean,nRmean];
X=[1:3];
plot(C,X,'o','MarkerSize',10,'MarkerFaceColor',[0.00,0.45,0.74]);
set(gca, 'Ytick',0:4,'fontsize',24);
xlim([-0.03 0.03])
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
%ax=get(gca, 'axis');
Y=zeros(length(yt));
hold on 
plot(Y,yt,'--r')
%axis square
yticklabels({[],'EP3 nonREM', 'EP3 REM' ,'EP2'})
title('Mean delta FC edges','fontsize',30)

%% LOCAL MEASURES: strenght & diversity FC  

for cond=1:5
    for sub=1:NSUB(cond)
        for area=1:N
         strenght_FC{cond}(sub,area)=sum(FC_emp{cond}(sub,area,:),3);
         diversity_FC{cond}(sub,area)=std(FC_emp{cond}(sub,area,:));
        end
    end
end

REM_SFC=strenght_FC{1}; %%stage 2
REM_REL_SFC=vertcat(strenght_FC{3},strenght_FC{4});%stage 3 REM
nonREM_SFC=vertcat(strenght_FC{2}); %stage 3 nonREM

condA=strenght_FC{5}; %ctrl 
condB=REM_SFC; 
condC=REM_REL_SFC;
condD=nonREM_SFC; 

x = 1 : N;
str_m1=mean(condA);
str_m2=mean(condB);
str_m3=mean(condC);
str_m4=mean(condD);
for seed=1:115
str_s1(seed)=std(condA(:,seed));
str_s2(seed)=std(condB(:,seed));
str_s3(seed)=std(condC(:,seed));
str_s4(seed)=std(condD(:,seed));
std_err1(seed)=str_s1(seed)/sqrt(size(condA,1));
std_err2(seed)=str_s2(seed)/sqrt(size(condB,1));
std_err3(seed)=str_s3(seed)/sqrt(size(condC,1));
std_err4(seed)=str_s4(seed)/sqrt(size(condD,1));
end
[Areas_sorted,indx] = sort(str_m1,'ascend'); % sorted from max to min

for n=1:N
    [P_C_R(n),H]=ranksum(condA(:,n),condB(:,n));
    [P_C_RR(n),H]=ranksum(condA(:,n),condC(:,n));
    [P_C_nR(n),H]=ranksum(condA(:,n),condD(:,n));
end
    
a=P_C_R(indx);astR=find(a<0.05);
b=P_C_RR(indx);astRR=find(b<0.05);
c=P_C_nR(indx);astnR=find(c<0.05);

bonf_P_C_R=P_C_R.*N;bonf_P_C_RR=P_C_RR.*N;  bonf_P_C_nR=P_C_nR.*N;
a_bonf=find(bonf_P_C_R<0.05);b_bonf=find(bonf_P_C_RR<0.05);c_bonf=find(bonf_P_C_nR<0.05);

sig_P_C_R=NaN(1,N); sR1=find(P_C_R<0.05);sig_P_C_R(sR1)=1;
sig_P_C_RR=NaN(1,N); sR2=find(P_C_RR<0.05);sig_P_C_RR(sR2)=1;
sig_P_C_nR=NaN(1,N); sR3=find(P_C_nR<0.05);sig_P_C_nR(sR3)=1;

figure(); 
X=[0:2];
% change str_m & std_err to plot cond of interest : 1=CTRL, 2=REM,3=REMREL, 4=nonREM
inBetween = [str_m1(indx)-std_err1(indx), fliplr(str_m1(indx)+std_err1(indx))]; %cond  1 Conf Interv
x2 = [x, fliplr(x)];
ppp=fill(x2, inBetween, 'b'); ppp.FaceAlpha=0.2; ppp.EdgeColor='none';
hold on
plot(x, str_m1(indx), 'b');%cond  1 mean
hold on
inBetween = [str_m3(indx)-std_err3(indx), fliplr(str_m3(indx)+std_err3(indx))]; %cond  2 Conf Interv
x2 = [x, fliplr(x)];
ppp=fill(x2, inBetween, 'r'); ppp.FaceAlpha=0.2; ppp.EdgeColor='none';
hold on
plot(x, str_m3(indx), 'r'); %cond  2 mean
hold on
plot(x,sig_P_C_RR(indx)+max(str_m1).*1.05,'*k','MarkerIndices',astnR);
set(gca, 'Xtick',1:2:N,'fontsize',22,'XTickLabelRotation',90);
xlim([0,115]);
aa=indx';
%xticks(1:N); xticklabels(StructNames(aa,:)); ax.XTickLabelRotation = 90;
%title('strenght CTRL, REM')
title('strenght CTRL, REM_REL')
%title('strenght CTRL, nonREM')
%title('Strength FC')

% MEAN DIFFERENCE of STRENgTH  (mean over diff between paired areas in 2cond, EP-ctrl)

H0=zeros(1,N);
diff_C_R=str_m2-str_m1;
diff_C_RR=str_m3-str_m1;
diff_C_nR=str_m4-str_m1;

diff_C_R_er=(str_m2-str_m1)./std_err1;
diff_C_RR_er=(str_m3-str_m1)./std_err1;
diff_C_nR_er=(str_m4-str_m1)./std_err1;

[P_diff_R_er(1),H]=ranksum(diff_C_R_er,H0); 
[P_diff_R_er(2),H]=ranksum(diff_C_RR_er,H0);
[P_diff_R_er(3),H]=ranksum(diff_C_nR_er,H0);
    
bonf_P_diff_R=P_diff_R_er.*3;

Rmean=mean(diff_C_R);
RRmean=mean(diff_C_RR);
nRmean=mean(diff_C_nR);

R_std=std(diff_C_R);
RR_std=std(diff_C_RR);
nR_std=std(diff_C_nR);
err_R=R_std/size(diff_C_R,2);
err_RR=RR_std/size(diff_C_RR,2);
err_nR=nR_std/size(diff_C_nR,2);
conf_99R=5.15.*err_R;
conf_99RR=5.15.*err_RR;
conf_99nR=5.15.*err_nR;

cohend_R=abs(Rmean)/R_std;
cohend_RR=abs(nRmean)/RR_std;
cohend_nR=abs(nRmean)/nR_std;


figure();
X=[1:3]; 
C=[ mean(diff_C_nR) mean(diff_C_RR) mean(diff_C_R)]; %vert
plot(C,X,'o','MarkerSize',10,'MarkerFaceColor',[0.00,0.45,0.74]);
%axis([xlim    min(yt)  ceil(max(yt)*1.2)])
set(gca, 'Ytick',0:4,'fontsize',24);
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
Y=zeros(length(yt));
hold on
plot(Y,yt,'--r')%vert
yticklabels({[],'EP3 nonREM', 'EP3 REM' ,'EP2'})
title('Mean delta cond-ctrl strength','fontsize',30)

save 'FC_comparisons'
