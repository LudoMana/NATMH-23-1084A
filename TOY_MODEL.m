%% TOY MODEL
%% Define the network
%The network will consist on 10 nodes and the first node will be the only node
%connected with the rest of the nodes (one hub).
nNodes=10;
C=zeros(nNodes,nNodes);

for n=1:3
C(n,1:3)=0.12;
end

C(4,:)=0.2;
C(:,4)=0.2;

for n=5:10
C(n,5:10)=0.12;
end

%plot te structural connectivity
figure
imagesc(C)
xlabel('Node')
ylabel('Node')
title('SC')
hold on
colormap('summer')


%% BUild the model

tic

% simulation parameters
dt = 0.01; % time step (in seconds)
Tmax=1000;
TRsec=2;
nSubs=1;
tvect=1:dt:Tmax*TRsec*nSubs; % vector of time stamps (in seconds)
Nt=length(tvect); % number of time points
xs=zeros(Tmax*nSubs,nNodes); % prepare variable for BOLD

% G value
G=0.3;

% omega
f=0.05*ones(nNodes,1); % frequency, i.e., number of cycles per unit time
omega = 2*pi*f; % angular velocity

% bifurcation parameter a
a_val=0;
a = a_val*ones(nNodes,1); % same for all nodes

% gC and strength
gC = G*C;
sumgC = sum(gC,2);

% noise
sig = 0.04;
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step

% Bifurcation parameters to study
a_hub_vals=0:-0.5:-3;
N_a_hub=length(a_hub_vals);

% external stimulation parameters
f_ext=0.05;
omega_ext=2*pi*f_ext;
A_ext=0.2;
z_ext=A_ext*exp(i*omega_ext*tvect); % apply function for all values of t
C_ext=0.2;
M_max=9;


% Prepare variables to store
FC=zeros(N_a_hub, M_max, nNodes, nNodes);

%% Run the simulation
for aidx=1:N_a_hub
    
    %value of the hub
    a_hub=a_hub_vals(aidx);
    a(4)=a_hub;

    % Initialize z and let it evolve under the effect of the network
    z_ini=0.1*complex(ones(nNodes,1),ones(nNodes,1)); % initialize z
    z=z_ini;
    for t=1:dt:3000 % Leaving 3000 seconds for the model to stabilize
        input = gC*z - sumgC.*z; % input from other nodes
        z = z + dt * ( a.*z - (abs(z).^2).*z + 1i*omega.*z + input) + dsig*complex(randn(nNodes,1),randn(nNodes,1));
    end
    z_ini=z; % Take the final z as the initial value for the next step

    for M=1:M_max
        z=z_ini;
        tp=1; % time steps for simulation
        tpbold=1; % time steps for BOLD
        for t=1:dt:Tmax*TRsec*nSubs
            input = gC*z - sumgC.*z;
            input(4) = input(4) + (M-1)*C_ext*(z_ext(tp)-z(4)); % add external input
            z = z + dt * ( a.*z - (abs(z).^2).*z + 1i*omega.*z + input) + dsig*complex(randn(nNodes,1),randn(nNodes,1));
            tp=tp+1;
            if abs(mod(t,TRsec))<dt/2
                xs(tpbold,:)=real(z)'; % save BOLD
                tpbold=tpbold+1; % Get index ready for next time step
            end
        end

        FC(aidx,M,:, :) = corrcoef(xs(1:end,:)); 
        
    end
end

toc

%% FIGURES

%subplots of the diff between FC matrix of perturbed - not perturbed network for each value of bifurcation parameter  

figure()
for iii=1:N_a_hub
    subplot(3,3,iii)
    clim=[0 1]
    imagesc(squeeze(FC(iii,9,:,:))-squeeze(FC(iii,1,:,:)),clim)
    colorbar()
end

%plot global FC for each value of bifurcation parameter (rows) at increasing levels of perturbation intensity (columns)
figure()
imagesc(mean(mean(FC,4),3))
















