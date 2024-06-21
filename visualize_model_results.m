%% Results CNT
for trial=1:300
    load(sprintf('model_hetero_cnt_trial_%02d',trial))
    ksP_all_cnt(trial,:)=ksP; %optimal fit for each g with opt a
    bifpar_all_cnt(trial,:,:)=bifpar;
    FC_all_cnt(:,:,:,trial)=FC_simul;
    SC_FC_g_all_cnt(:,trial)=SC_FC_g;
    %fitting_all_cnt(:,trial)=fitting;
end

ksP_mean=mean(ksP_all_cnt);
ksP_median=median(ksP_all_cnt);

figure();
plot(ksP_mean)
figure();
plot(ksP_median)

%Bif par reorder again important to recover the original parcellation
mbifpar_cnt=squeeze(mean(bifpar_all_cnt,1));
mbifpar_cnt_ord=[mbifpar_cnt(:,1:2:112) mbifpar_cnt(:,2:2:112) mbifpar_cnt(:,113:115)];
bifpar_cnt_ord=permute(bifpar_all_cnt,[3,1,2]);
bifpar_cnt_ord=[bifpar_cnt_ord(1:2:112,:,:); bifpar_cnt_ord(2:2:112,:,:); bifpar_cnt_ord(113:115,:,:)];
bifpar_cnt_ord=permute(bifpar_cnt_ord,[2,3,1]);

[min_g,opt_g_cnt]=min(ksP_mean);
opt_mbifpar_cnt=mbifpar_cnt_ord(opt_g_cnt,:); %select optimal g!!
opt_bifpar_cnt=squeeze(bifpar_cnt_ord(:,opt_g_cnt,:));
%mean_opt_g=[];
WE = 0:0.05:2.5;
for i=1:300
    [v,g_OPT_cnt(i)]=min(ksP_all_cnt(i,:));
    val1=g_OPT_cnt(i); g_OPT_cnt_val(i)=WE(val1);
end

mean_g_OPT_cnt=mean(g_OPT_cnt_val);std_g_OPT_cnt=std(g_OPT_cnt_val);err_g_OPT_cnt=std_g_OPT_cnt./sqrt(size(g_OPT_cnt_val,2));

figure();
for g=1:51
h=histfit(ksP_all_cnt(:,g),30,'kernel');
h(1).FaceColor = 'none';h(1).EdgeColor = 'none';h(2).LineWidth = 0.5;h(2).Color = [0.65,0.65,0.65];
hold on
end
h=histfit(ksP_all_cnt(:,22),30,'kernel');
h(1).FaceColor = 'none';h(1).EdgeColor = 'none';h(2).LineWidth = 1.2;h(2).Color = 'r';

%opt_g_cnt=29;
FC_opt_cnt=squeeze(FC_all_cnt(:,:,opt_g_cnt,:));
FCmean_opt_cnt=squeeze(mean(FC_opt_cnt,3));

clear bifpar ksP FC_simul
save model_unified_SAVEFC_cnt.mat

figure();plot(mean(fitting_all_cnt,2))
figure();plot(median(fitting_all_cnt,2))

% KS distance plot
hopfModel.ksdist_cnt_median   = nanmedian(fitting_all_cnt');
hopfModel.ksdist_cnt_mean   = mean(fitting_all_cnt');
hopfModel.ksdist_cnt_25qntile = quantile(fitting_all_cnt',0.025);
hopfModel.ksdist_cnt_75qntile = quantile(fitting_all_cnt',0.975);
hopfModel.ksdist_cnt_25qntile(isnan(hopfModel.ksdist_cnt_25qntile)) = 0
hopfModel.ksdist_cnt_75qntile(isnan(hopfModel.ksdist_cnt_75qntile)) = 0
%[KSdist_cnt idx_ksdist_cnt] = min(hopfModel.ksdist_cnt_median);

figure()
xVec = 1:size(hopfModel.ksdist_cnt_mean,2);
xConfnoresp= [xVec xVec(end:-1:1)];
yConfnoresp = [(hopfModel.ksdist_cnt_25qntile) (hopfModel.ksdist_cnt_75qntile(end:-1:1))];
p4=fill(xConfnoresp,yConfnoresp,'blue');
p4.FaceColor = [0.8 0.8 1];p4.EdgeColor = 'none';p4.FaceAlpha = 0.5;
hold on
plot(hopfModel.ksdist_cnt_median,'b-','LineWidth',2)
%plot(hopfModel.ksdist_cnt_mean,'b-','LineWidth',2)


%% Results EP
clear
%Ep2

for trial=1:300
    load(sprintf('model_hetero_EP2_trial_%02d',trial))
    ksP_all_ep2(trial,:)=ksP;
    bifpar_all_ep2(trial,:,:)=bifpar;
    FC_all_ep2(:,:,:,trial)=FC_simul;
    SC_FC_g_all_ep2(:,trial)=SC_FC_g;
    fitting_all_ep2(:,trial)=fitting;
end

ksP_mean=mean(ksP_all_ep2);

figure();
plot(ksP_mean)

mbifpar_ep2=squeeze(mean(bifpar_all_ep2,1));
mbifpar_ep2_ord=[mbifpar_ep2(:,1:2:112) mbifpar_ep2(:,2:2:112) mbifpar_ep2(:,113:115)];
bifpar_ep2_ord=permute(bifpar_all_ep2,[3,1,2]);
bifpar_ep2_ord=[bifpar_ep2_ord(1:2:112,:,:); bifpar_ep2_ord(2:2:112,:,:); bifpar_ep2_ord(113:115,:,:)];
bifpar_ep2_ord=permute(bifpar_ep2_ord,[2,3,1]);

[min_g,opt_g_ep2]=min(ksP_mean);
opt_mbifpar_ep2=mbifpar_ep2_ord(opt_g_ep2,:); %select optimal g!!
opt_bifpar_ep2=squeeze(bifpar_ep2_ord(:,opt_g_ep2,:));

WE = 0:0.05:2.5;
for i=1:300
    [v,g_OPT_ep2(i)]=min(ksP_all_ep2(i,:));
    val1=g_OPT_ep2(i); g_OPT_ep2_val(i)=WE(val1);
end
mean_g_OPT_ep2=mean(g_OPT_ep2_val);std_g_OPT_ep2=std(g_OPT_ep2_val);err_g_OPT_ep2=std_g_OPT_ep2./sqrt(size(g_OPT_ep2_val,2));

%opt_g_ep2=28;
FC_opt_ep2=squeeze(FC_all_ep2(:,:,opt_g_ep2,:));
FCmean_opt_ep2=squeeze(mean(FC_opt_ep2,3));
clear bifpar ksP FC_simul
%save model_unified_SAVEFC_ep2_50.mat

figure();
for g=1:51
h=histfit(ksP_all_ep2(:,g),30,'kernel');
h(1).FaceColor = 'none';h(1).EdgeColor = 'none';h(2).LineWidth = 0.5;h(2).Color = [0.65,0.65,0.65];
hold on
end
h=histfit(ksP_all_ep2(:,24),30,'kernel');
h(1).FaceColor = 'none';h(1).EdgeColor = 'none';h(2).LineWidth = 1.2;h(2).Color = 'r';

%%
%Ep3rem
clear

for trial=1:300
    load(sprintf('model_hetero_EP3rem_trial_%02d',trial))
    ksP_all_ep3rem(trial,:)=ksP;
    bifpar_all_ep3rem(trial,:,:)=bifpar;
    FC_all_ep3r(:,:,:,trial)=FC_simul;
    SC_FC_g_all_ep3R(:,trial)=SC_FC_g;
    fitting_all_ep3R(:,trial)=fitting;
end

ksP_mean=mean(ksP_all_ep3rem);

figure();
plot(ksP_mean)

mbifpar_ep3rem=squeeze(mean(bifpar_all_ep3rem,1));
mbifpar_ep3rem_ord=[mbifpar_ep3rem(:,1:2:112) mbifpar_ep3rem(:,2:2:112) mbifpar_ep3rem(:,113:115)];
bifpar_ep3rem_ord=permute(bifpar_all_ep3rem,[3,1,2]);
bifpar_ep3rem_ord=[bifpar_ep3rem_ord(1:2:112,:,:); bifpar_ep3rem_ord(2:2:112,:,:); bifpar_ep3rem_ord(113:115,:,:)];
bifpar_ep3rem_ord=permute(bifpar_ep3rem_ord,[2,3,1]);

[min_g,opt_g_ep3rem]=min(ksP_mean);
opt_mbifpar_e3rem=mbifpar_ep3rem_ord(opt_g_ep3rem,:); %select optimal g!!
opt_bifpar_ep3rem=squeeze(bifpar_ep3rem_ord(:,opt_g_ep3rem,:));

WE = 0:0.05:2.5;
for i=1:300
    [v,g_OPT_ep3r(i)]=min(ksP_all_ep3rem(i,:));
    val1=g_OPT_ep3r(i); g_OPT_ep3r_val(i)=WE(val1);
end
mean_g_OPT_ep3r=mean(g_OPT_ep3r_val);std_g_OPT_ep3r=std(g_OPT_ep3r_val);err_g_OPT_ep3r=std_g_OPT_ep3r./sqrt(size(g_OPT_ep3r_val,2));

%opt_g_ep3rem=21;
FC_opt_ep3r=squeeze(FC_all_ep3r(:,:,opt_g_ep3rem,:));
FCmean_opt_ep3r=squeeze(mean(FC_opt_ep3r,3));

clear bifpar ksP FC_simul
save model_unified_SAVEFC_ep3rem.mat

figure();
for g=1:51
h=histfit(ksP_all_ep3rem(:,g),30,'kernel');
h(1).FaceColor = 'none';h(1).EdgeColor = 'none';h(2).LineWidth = 0.5;h(2).Color = [0.65,0.65,0.65];
hold on
end
h=histfit(ksP_all_ep3rem(:,18),30,'kernel');
h(1).FaceColor = 'none';h(1).EdgeColor = 'none';h(2).LineWidth = 1.2;h(2).Color = 'r';

%%
%Ep3Nrem
clear

for trial=1:300
    load(sprintf('model_hetero_EP3Nrem_trial_%02d',trial))
    ksP_all_ep3Nrem(trial,:)=ksP;
    bifpar_all_ep3Nrem(trial,:,:)=bifpar;
    FC_all_ep3Nr(:,:,:,trial)=FC_simul;
    SC_FC_g_all_ep3NR(:,trial)=SC_FC_g;
    fitting_all_ep3NR(:,trial)=fitting;
end

ksP_mean=mean(ksP_all_ep3Nrem);

figure();
plot(ksP_mean)

mbifpar_ep3Nrem=squeeze(mean(bifpar_all_ep3Nrem,1));
mbifpar_ep3Nrem_ord=[mbifpar_ep3Nrem(:,1:2:112) mbifpar_ep3Nrem(:,2:2:112) mbifpar_ep3Nrem(:,113:115)];
bifpar_ep3Nrem_ord=permute(bifpar_all_ep3Nrem,[3,1,2]);
bifpar_ep3Nrem_ord=[bifpar_ep3Nrem_ord(1:2:112,:,:); bifpar_ep3Nrem_ord(2:2:112,:,:); bifpar_ep3Nrem_ord(113:115,:,:)];
bifpar_ep3Nrem_ord=permute(bifpar_ep3Nrem_ord,[2,3,1]);

[min_g,opt_g_ep3Nrem]=min(ksP_mean);
opt_mbifpar_e3Nrem=mbifpar_ep3Nrem_ord(opt_g_ep3Nrem,:); %select optimal g!!
opt_bifpar_ep3Nrem=squeeze(bifpar_ep3Nrem_ord(:,opt_g_ep3Nrem,:));

WE = 0:0.05:2.5;
for i=1:300
    [v,g_OPT_ep3Nr(i)]=min(ksP_all_ep3Nrem(i,:));
    val1=g_OPT_ep3Nr(i); g_OPT_ep3Nr_val(i)=WE(val1);
end
mean_g_OPT_ep3Nr=mean(g_OPT_ep3Nr_val);std_g_OPT_ep3Nr=std(g_OPT_ep3Nr_val);err_g_OPT_ep3Nr=std_g_OPT_ep3Nr./sqrt(size(g_OPT_ep3Nr_val,2));

%opt_g_ep3Nrem=16;
FC_opt_ep3Nr=squeeze(FC_all_ep3Nr(:,:,opt_g_ep3Nrem,:));
FCmean_opt_ep3Nr=squeeze(mean(FC_opt_ep3Nr,3));

clear bifpar ksP 
save model_unified_SAVEFC_ep3Nrem.mat

figure();
for g=1:51
h=histfit(ksP_all_ep3Nrem(:,g),30,'kernel');
h(1).FaceColor = 'none';h(1).EdgeColor = 'none';h(2).LineWidth = 0.5;h(2).Color = [0.65,0.65,0.65];
hold on
end
h=histfit(ksP_all_ep3Nrem(:,14),30,'kernel');
h(1).FaceColor = 'none';h(1).EdgeColor = 'none';h(2).LineWidth = 1.2;h(2).Color = 'r';


%% homo model

for trial=1:300
    load(sprintf('model_homo_cnt_trial_%02d',trial))
    ksP_all_cnt_omo(trial,:)=ksP; %optimal fit for each g with opt a
    bifpar_all_cnt_omo(trial,:,:)=bifpar;
    FC_all_cnt_omo(:,:,:,trial)=FC_simul;
    fitting_all_cnt_omo(:,trial)=fitting;
end

ksP_mean=mean(ksP_all_cnt_omo);

figure();
plot(ksP_mean)

mbifpar_cnt_omo=squeeze(mean(bifpar_all_cnt_omo,1));
mbifpar_cnt_omo_ord=[mbifpar_cnt_omo(:,1:2:112) mbifpar_cnt_omo(:,2:2:112) mbifpar_cnt_omo(:,113:115)];
bifpar_cnt_omo_ord=permute(bifpar_all_cnt_omo,[3,1,2]);
bifpar_cnt_omo_ord=[bifpar_cnt_omo_ord(1:2:112,:,:); bifpar_cnt_omo_ord(2:2:112,:,:); bifpar_cnt_omo_ord(113:115,:,:)];
bifpar_cnt_omo_ord=permute(bifpar_cnt_omo_ord,[2,3,1]);

[min_g,opt_g_cnt_omo]=min(ksP_mean);
opt_mbifpar_cnt_omo=mbifpar_cnt_omo_ord(opt_g_cnt_omo,:); %select optimal g!!
opt_bifpar_cnt_omo=squeeze(bifpar_cnt_omo_ord(:,opt_g_cnt_omo,:));
%mean_opt_g=[];

WE = 0:0.05:2.5;
for i=1:50
    [v,g_OPT_cnt_omo(i)]=min(ksP_all_cnt_omo(i,:));
    val1=g_OPT_cnt_omo(i); g_OPT_cnt_val_omo(i)=WE(val1);
end

mean_g_OPT_cnt_omo=mean(g_OPT_cnt_val_omo);std_g_OPT_cnt_omo=std(g_OPT_cnt_val_omo);err_g_OPT_cnt_omo=std_g_OPT_cnt_omo./sqrt(size(g_OPT_cnt_val_omo,2));

%opt_g_cnt_omo=40;
FC_opt_cnt_omo=squeeze(FC_all_cnt_omo(:,:,opt_g_cnt_omo,:));
FCmean_opt_cnt_omo=squeeze(mean(FC_opt_cnt_omo,3));

clear bifpar ksP FC_simul

save model_unified_SAVEFC_cnt_HOMO.mat

figure();plot(mean(fitting_all_cnt_omo,2));

% FC distance plot
hopfModel.ksdist_cnt_median   = nanmedian(fitting_all_cnt_omo');
hopfModel.ksdist_cnt_mean   = mean(fitting_all_cnt_omo');
hopfModel.ksdist_cnt_25qntile = quantile(fitting_all_cnt_omo',0.025);
hopfModel.ksdist_cnt_75qntile = quantile(fitting_all_cnt_omo',0.975);
hopfModel.ksdist_cnt_25qntile(isnan(hopfModel.ksdist_cnt_25qntile)) = 0
hopfModel.ksdist_cnt_75qntile(isnan(hopfModel.ksdist_cnt_75qntile)) = 0
%[KSdist_cnt idx_ksdist_cnt] = min(hopfModel.ksdist_cnt_median);

figure()
xVec = 1:size(hopfModel.ksdist_cnt_mean,2);
xConfnoresp= [xVec xVec(end:-1:1)];
yConfnoresp = [(hopfModel.ksdist_cnt_25qntile) (hopfModel.ksdist_cnt_75qntile(end:-1:1))];
p4=fill(xConfnoresp,yConfnoresp,'blue');
p4.FaceColor = [0.8 0.8 1];p4.EdgeColor = 'none';p4.FaceAlpha = 0.5;
hold on
plot(hopfModel.ksdist_cnt_median,'b-','LineWidth',2)
%plot(hopfModel.ksdist_cnt_mean,'b-','LineWidth',2)

%% FIGURES
load model_unified_SAVEFC_cnt;
load model_unified_SAVEFC_cnt_HOMO;
load model_unified_SAVEFC_ep2;
load model_unified_SAVEFC_ep3Nrem;
load model_unified_SAVEFC_ep3rem;


FCmean_opt_cnt_ord=[FCmean_opt_cnt(:,1:2:112) FCmean_opt_cnt(:,2:2:112) FCmean_opt_cnt(:,113:115)];
FCmean_opt_cnt_ord=[FCmean_opt_cnt_ord(1:2:112,:)' FCmean_opt_cnt_ord(2:2:112,:)' FCmean_opt_cnt_ord(113:115,:)'];

FCmean_opt_cnt_omo_ord=[FCmean_opt_cnt_omo(:,1:2:112) FCmean_opt_cnt_omo(:,2:2:112) FCmean_opt_cnt_omo(:,113:115)];
FCmean_opt_cnt_omo_ord=[FCmean_opt_cnt_omo_ord(1:2:112,:)' FCmean_opt_cnt_omo_ord(2:2:112,:)' FCmean_opt_cnt_omo_ord(113:115,:)'];

FCmean_opt_ep2_ord=[FCmean_opt_ep2(:,1:2:112) FCmean_opt_ep2(:,2:2:112) FCmean_opt_ep2(:,113:115)];
FCmean_opt_ep2_ord=[FCmean_opt_ep2_ord(1:2:112,:)' FCmean_opt_ep2_ord(2:2:112,:)' FCmean_opt_ep2_ord(113:115,:)'];

FCmean_opt_ep3r_ord=[FCmean_opt_ep3r(:,1:2:112) FCmean_opt_ep3r(:,2:2:112) FCmean_opt_ep3r(:,113:115)];
FCmean_opt_ep3r_ord=[FCmean_opt_ep3r_ord(1:2:112,:)' FCmean_opt_ep3r_ord(2:2:112,:)' FCmean_opt_ep3r_ord(113:115,:)'];

FCmean_opt_ep3Nr_ord=[FCmean_opt_ep3Nr(:,1:2:112) FCmean_opt_ep3Nr(:,2:2:112) FCmean_opt_ep3Nr(:,113:115)];
FCmean_opt_ep3Nr_ord=[FCmean_opt_ep3Nr_ord(1:2:112,:)' FCmean_opt_ep3Nr_ord(2:2:112,:)' FCmean_opt_ep3Nr_ord(113:115,:)'];


figure();
subplot(1,2,1)
imagesc(FCmean_opt_cnt); axis square;set(gca,'fontsize', 14); 
subplot(1,2,2)
imagesc(FCmean_opt_cnt_ord); axis square;set(gca,'fontsize', 14); 

clim=[-0.06,0.42];
figure();
subplot(2,3,1)
% imagesc(FCmean_opt_cnt_omo_ord,clim); axis square;set(gca,'fontsize', 14); 
% title 'FC simul CNT homo'
subplot(2,3,2)
imagesc(FCmean_opt_cnt_ord,clim); axis square;set(gca,'fontsize', 14)
title 'FC simul CNT'
subplot(2,3,4)
imagesc(FCmean_opt_ep2_ord,clim); axis square;set(gca,'fontsize', 14)
title 'FC simul EP2'
subplot(2,3,5)
imagesc(FCmean_opt_ep3r_ord,clim); axis square;set(gca,'fontsize', 14)
title 'FC simul EP3R'
subplot(2,3,6)
imagesc(FCmean_opt_ep3Nr_ord,clim); axis square; set(gca,'fontsize', 14)
title 'FC simul EP3NR'

load FC_emp_allCOND
FC_emp_cnt=squeeze(mean(FC_emp{6}));
FC_emp_ep2=squeeze(mean(FC_emp{1}));
FC_emp_ep3rem=squeeze(mean(vertcat(FC_emp{3},FC_emp{4})));
FC_emp_ep3Nrem=squeeze(mean(vertcat(FC_emp{2},FC_emp{5})));

clim=[-0.06,0.5];
figure();
subplot(2,3,1)
imagesc(FC_emp_cnt,clim); axis square;set(gca,'fontsize', 14); 
title 'FC emp CNT homo'
subplot(2,3,4)
imagesc(FC_emp_ep2,clim); axis square;set(gca,'fontsize', 14)
title 'FC emp EP2'
subplot(2,3,5)
imagesc(FC_emp_ep3rem,clim); axis square;set(gca,'fontsize', 14)
title 'FC emp EP3R'
subplot(2,3,6)
imagesc(FC_emp_ep3Nrem,clim); axis square; set(gca,'fontsize', 14)
title 'FC emp EP3NR'

global_FC_cnt_omo=mean(FCmean_opt_cnt_omo,'all');
global_FC_cnt=mean(FCmean_opt_cnt,'all');
global_FC_ep2=mean(FCmean_opt_ep2,'all');
global_FC_ep3r=mean(FCmean_opt_ep3r,'all');
global_FC_ep3Nr=mean(FCmean_opt_ep3Nr,'all');

global_FC_cnt_omo=median(FCmean_opt_cnt_omo,'all');
global_FC_cnt=median(FCmean_opt_cnt,'all');
global_FC_ep2=median(FCmean_opt_ep2,'all');
global_FC_ep3r=median(FCmean_opt_ep3r,'all');
global_FC_ep3Nr=median(FCmean_opt_ep3Nr,'all');

figure();histogram(FCmean_opt_cnt);

global_FC_cnt_omo=mean(FCmean_opt_cnt_omo,[1,2]);
global_FC_cnt=squeeze(mean(FC_opt_cnt,[1,2]));
global_FC_ep2=squeeze(mean(FC_opt_ep2,[1,2]));
global_FC_ep3r=squeeze(mean(FC_opt_ep3r,[1,2]));
global_FC_ep3Nr=squeeze(mean(FC_opt_ep3Nr,[1,2]));

P_FC_ep2=ranksum(global_FC_cnt,global_FC_ep2);bonf_P_FC_ep2=P_FC_ep2.*3
P_FC_ep3r=ranksum(global_FC_cnt,global_FC_ep3r);bonf_P_FC_ep3r=P_FC_ep3r.*3
P_FC_ep3Nr=ranksum(global_FC_cnt,global_FC_ep3Nr);bonf_P_FC_ep3Nr=P_FC_ep3Nr.*3

figure();
C=[global_FC_cnt global_FC_ep2 global_FC_ep3r global_FC_ep3Nr] %change cond to compare
%grp = [zeros(1,size(FCmean_opt_cnt,2)),ones(1,size(FCmean_opt_ep2,2)),ones(1,size(FCmean_opt_ep3r,2))*2,ones(1,size(FCmean_opt_ep3Nr,2))*3];
boxplot(C)
yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
set(gca, 'Xtick', 1:4);
xt = get(gca, 'XTick');


N=115;
%strent
for iter=1:300
        for area=1:N
         strenght_FC{1}(iter,area)=sum(FC_opt_cnt(area,:,iter),2);
         strenght_FC{2}(iter,area)=sum(FC_opt_ep2(area,:,iter),2);
         strenght_FC{3}(iter,area)=sum(FC_opt_ep3r(area,:,iter),2);
         strenght_FC{4}(iter,area)=sum(FC_opt_ep3Nr(area,:,iter),2);
        end
end

str_m1=mean(strenght_FC{1});
str_m2=mean(strenght_FC{2});
str_m3=mean(strenght_FC{3});
str_m4=mean(strenght_FC{4});

for seed=1:115
str_s1(seed)=std(strenght_FC{1}(:,seed));
str_s2(seed)=std(strenght_FC{2}(:,seed));
str_s3(seed)=std(strenght_FC{3}(:,seed));
str_s4(seed)=std(strenght_FC{4}(:,seed));
std_err1(seed)=str_s1(seed)/sqrt(50);
std_err2(seed)=str_s2(seed)/sqrt(50);
std_err3(seed)=str_s3(seed)/sqrt(50);
std_err4(seed)=str_s4(seed)/sqrt(50);
end

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
title('Mean ? cond-ctrl strength','fontsize',30)

figure()
boxplot([g_OPT_cnt_val' g_OPT_ep2_val' g_OPT_ep3r_val' g_OPT_ep3Nr_val'])

% fit FC

%plot opt A ordered optA each cond
[Areas_sorted_opta,indx_optacnt] = sort(opt_mbifpar_cnt,'ascend');
[Areas_sorted_opta,indx_opta2] = sort(opt_mbifpar_ep2,'ascend');
[Areas_sorted_opta,indx_opta3R] = sort(opt_mbifpar_e3rem,'ascend');
[Areas_sorted_opta,indx_opta3Nr] = sort(opt_mbifpar_e3Nrem,'ascend');

figure()
plot(opt_mbifpar_cnt(indx_optacnt))
hold on
plot(opt_mbifpar_ep2(indx_opta2))
hold on
plot(opt_mbifpar_e3rem(indx_opta3R))
hold on
plot(opt_mbifpar_e3Nrem(indx_opta3Nr))
hold on
plot(zeros(115,1),'k')
set(gca,'fontsize',16,'xlim',[0,115]);
ylabel('Opt A');
xlabel('Ordered optA'); 

%area below curve
for iter = 1:300
    area_cnt(iter,:)=sum(abs(opt_bifpar_cnt(iter,:)));
    area_ep2(iter,:)=sum(abs(opt_bifpar_ep2(iter,:)));
    area_ep3R(iter,:)=sum(abs(opt_bifpar_ep3rem(iter,:)));
    %area_ep3NR(iter,:)=sum(abs(opt_bifpar_ep3NR(iter,:)));
end

figure()
boxplot([area_cnt,area_ep2,area_ep3R])