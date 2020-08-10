clear all

load('ArchEfficiency_Simulations_workspace')
%load('data/ArchEfficiency_Simulations_workspace')

%% Promoters of interest
query = {'flp21';'npr4';'sra11';'odr2b';'dop2'};   %strong slow down

% # of unique neurons
uniqueNeurons = cell(5,1);
allNames = {allProm.name};
allQNeurons = {};   %all queried neurons
allQProm = {};  %all queried promoters
[allProm.pheno] = deal(0);
for i = 1:length(query)
    [C,ia,ib] = intersect(query{i},allNames);
    allProm(ib).pheno = 1;
    tempNeurons = allProm(ib).neurons';
    allQNeurons = [allQNeurons; tempNeurons];  %all neurons (including duplicates)
    allQProm = [allQProm; repmat({query{i}},length(tempNeurons),1)];
end
allQNeurons = [allQNeurons allQProm];

uniqueNeurons = unique(allNeurons(:,1));


% frequency of neurons
for i = 1:length(uniqueNeurons(:,1))
    tempIsMember = ismember(allQNeurons(:,1),uniqueNeurons{i,1});
    uniqueNeurons(i,3) = {allQNeurons(tempIsMember,2)};
    uniqueNeurons{i,4} = length(uniqueNeurons{i,3});
    tempIsMemberAll = ismember(allNeurons(:,1),uniqueNeurons{i,1});
    uniqueNeurons{i,5} = allNeurons(tempIsMemberAll,2);
    uniqueNeurons{i,6} = length(uniqueNeurons{i,5});
end

% class of unique neurons
for i = 1:length(uniqueNeurons(:,1))
    if ismember(uniqueNeurons{i,1},sensory)
        uniqueNeurons{i,2} = 'sensory';
    elseif ismember(uniqueNeurons{i,1},inter)
        uniqueNeurons{i,2} = 'inter';
    elseif ismember(uniqueNeurons{i,1},motor)
        uniqueNeurons{i,2} = 'motor';
    else
        uniqueNeurons{i,2} = 'unknown';
    end
end

% organize into struct
neurons = struct;
for i=1:length(uniqueNeurons(:,1))
    neurons(i).name = uniqueNeurons{i,1};
    neurons(i).class = uniqueNeurons{i,2};  %class: sensory,inter,motor,unknown
    neurons(i).promQ = uniqueNeurons{i,3};  %queried promoters expressing neurons
    neurons(i).nRepQ = uniqueNeurons{i,4};   %# of promoters expressing neuron
    neurons(i).promAll = uniqueNeurons{i,5};   %promoters (all screen) expressing neuron
    neurons(i).nRepAll = uniqueNeurons{i,6};    %number of promoters in screen expressing neuron
end

% class of the unique neurons
interInd = find(ismember({uniqueNeurons{:,2}},'inter'));
% allInter = {allUniqueNeurons{interInd,:}};
motorInd = find(ismember({uniqueNeurons{:,2}},'motor'));
sensoryInd = find(ismember({uniqueNeurons{:,2}},'sensory'));
unknownInd = find(ismember({uniqueNeurons{:,2}},'unknown'));
% distr. of neurons per promoter

% 
%% Find Sparse Solution
% Generate Neuron in Promoter Sampling Matrix
uniqueNeurons = unique(allNeurons(:,1));
uniqueNeurons = uniqueNeurons([interInd,sensoryInd,motorInd]);

A = zeros(length(allProm),length(uniqueNeurons));   % original measurement matrix;
for i=1:length(allProm)
    for j=1:length(allProm(i).neurons)
        ind = find(strcmp(allProm(i).neurons{j}, uniqueNeurons), 1);
        A(i,ind)=1;
    end
end

% Phenotype Vector
y = zeros(length(allProm),1);
y([allProm.pheno]==1) = PercentageMeanSpeedChange([1,6,21,23,27]);
y=-y;


%% CHANGE MEASUREMENT MATRIX "A" with Error

lam=21; %lambda choosen for publication in Fig1
B_Arch = zeros(1000,88);
B_MSE = zeros(1000,1);
B_ArchCell = cell(1000,1);
B_ArchCell_FitInfo = cell(1000,1);
M_Cell = cell(1000,1);
LessEffArch=[0.0 0.1 0.2 0.3 0.4 0.5];

for m=1:1000
            
    clear B FitInfo
    
    % Generate a Measurement Matrix with errors
    Anew = zeros(length(allProm),length(uniqueNeurons));
    count=0;
    for i=1:length(allProm)
        for j=1:length(allProm(i).neurons)
            ind = find(strcmp(allProm(i).neurons{j}, uniqueNeurons), 1);
            if rand>0.1  % 10 percent error
                Anew(i,ind)=1;
            else
                Anew(i,ind)=LessEffArch(randi(6));
                count=count+1;
            end
        end
    end
    
    M_Cell{m}=Anew;

    %Lasso
    [B,FitInfo] = lasso(Anew,y,'Alpha',1,'RelTol',1e-4, 'Lambda', lambda);


    B_Arch(m,:) = B(:,lam);
    B_ArchCell{m} = B;
    B_MSE(m) =  FitInfo.MSE(lam);
    B_ArchCell_FitInfo{m} = FitInfo;    
        
end


%% PLOT
figure('units','normalized','outerposition',[0 0 1 1])
imagesc(B_Arch,[-20 20])
colorbar

xTicks = [1:length(B(:,lam))];
xTickLabels = uniqueNeurons;
ax = gca;
set(ax,'XTick',xTicks);
set(ax,'XTickLabel',xTickLabels,'FontSize',12);
ax.XTickLabelRotation = -90;
xlim([0.5 length(uniqueNeurons)+0.5]);

set(gca,'ytick',[])

if ~exist('./results', 'dir')
    mkdir('./results');
end

saveas(gcf, './results/Arch_efficiency_1000correptedMatrices.png');

%% BOX PLOT
figure('units','normalized','outerposition',[0 0 1 1])
boxplot_pwhisker(B_Arch,{},10,90);   %You need boxplot_pwhisker

xTicks = [1:88];
xTickLabels = uniqueNeurons;
ax = gca;
set(ax,'FontSize',13);
set(ax,'XTick',xTicks);
set(ax,'XTickLabel',xTickLabels);
ax.XTickLabelRotation = -90;
xlim([0 length(uniqueNeurons)+1]);
ylim([-20 80])

if ~exist('./results', 'dir')
    mkdir('./results');
end

saveas(gcf, './results/Robustness_of_the_solution_to_errors.png');

