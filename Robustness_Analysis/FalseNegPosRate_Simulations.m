function FalseNegPosRate_Simulations(dataset)

  %% False Negative / False Positive Rate Simulations
  % 
  load(dataset)
  %% load workspace for necessary variables
  %load('data/simulation_workspace_M27')   %27 promoters
  %load('data/simulation_workspace_M32')  %32 promoters
  %load('data/simulation_workspace_M38')  %38 promoters
  
  %% Create necessaary variables
  uniqueNeurons = unique(allNeurons(:,1));
  
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
  
  
  % class of the unique neurons
  interInd = find(ismember({uniqueNeurons{:,2}},'inter'));
  % allInter = {allUniqueNeurons{interInd,:}};
  motorInd = find(ismember({uniqueNeurons{:,2}},'motor'));
  sensoryInd = find(ismember({uniqueNeurons{:,2}},'sensory'));
  unknownInd = find(ismember({uniqueNeurons{:,2}},'unknown'));
  
  %% Generate Measurement Matrix
  uniqueNeurons = unique(allNeurons(:,1));
  uniqueNeurons = uniqueNeurons([interInd,sensoryInd,motorInd]);
  
  A = zeros(length(allProm),length(uniqueNeurons));   %measurement matrix;
  for i=1:length(allProm)
      for j=1:length(allProm(i).neurons)
          ind = find(strcmp(allProm(i).neurons{j}, uniqueNeurons), 1);
          A(i,ind)=1;
      end
  end
  
  %% SIMULATIONS
  
  %% Lasso Solution
  N=1000;
  results_MatlabLasso = cell(N,11);
  results_MatlabLasso_inter = cell(N,11);
  results_y = cell(N,11);
  for k=[1:11]
      for indk=1:N % see top
          % Generate Phenotype Measurements
          T = randperm(size(A,2),k);  %randomly pick k neurons as key
          xSol = zeros(size(A,2),1);
          xSol(T)=rand(length(T),1);
          for i=1:length(T)
              if xSol(T(i))>=.25
                  xSol(T(i))=30+70*rand;
              else
                  xSol(T(i))=-100+70*rand;
              end
          end
  
          y = A*xSol; %generate data
  %         y(y<0)=-1; y(y>0)=1; %For nonlinear phenotype vector  %UNCOMMENT for nonlinear phenotype vector
          
          %Inference - Lasso Solution 
          x = lasso(A,y,'Alpha',1,'RelTol',1e-4,'Lambda',lambda); %EDIT lambda as lambdaBinary for nonlinear phenotype vector
                
          results_MatlabLasso{indk,k}=x;
          x(length(interInd)+1:end,:)=0;
          results_MatlabLasso_inter{indk,k}=x;
          results_y{indk,k}=y;
          resultsInter{indk,k}=xSol;
          resultsInter{indk,k}(length(interInd)+1:end,:) = 0;
          if mod(indk,50)==0
              [k indk]
          end
      end
  end
  
  
  %% Calculate FALSE NEG & POS
  falseNegMeans=zeros(1,11);
  falseNegStds=zeros(1,11);
  falsePosMeans=zeros(1,11);
  falsePosStds=zeros(1,11);
  g=21;
  for key=[1:11]
  % %% Calc overall coherence/incoherence of actual vs. results
  coherence = {};
  actual = {};
  measured = {};
  
  query = key;
  for i=1:size(resultsInter,1)
      
      temp = resultsInter{i,query};
      temp2 = results_MatlabLasso_inter{i,query}(:,g);
      
      actualSol = abs(temp)>.001;    
      measuredSol = abs(temp2)>.001;
  
      coherence{i}=actualSol.*measuredSol;
      actual{i}=actualSol;
      measured{i}=measuredSol;
  end
  
  temp = horzcat(coherence{:});
  tempMeasured = horzcat(measured{:});
  tempActual = horzcat(actual{:});
  
  % %% Calc false positives/negatives
  falsePos = [];
  falseNeg = [];
  for i=1:size(tempActual,2)
      actualInd = find(tempActual(:,i));
      measuredInd = find(tempMeasured(:,i));
      hits = intersect(actualInd,measuredInd);
      fNeg = setdiff(actualInd,hits);
      fPos = setdiff(measuredInd,hits);
      falsePos(i)=length(fPos);
      falseNeg(i)=length(fNeg);
  end
  
  %All Neg Pos
  falseNegAll(key,:)=falseNeg;
  falsePosAll(key,:)=falsePos;
  
  %Neg
  falseNegMeans(key)=mean(falseNeg);
  falseNegStds(key)=std(falseNeg);
  falseNegMedians(key)=median(falseNeg);
  %Pos
  falsePosMeans(key)=mean(falsePos);
  falsePosStds(key)=std(falsePos);
  falsePosMedians(key)=median(falsePos);
  
  end
  
  %% Plot False Neg/Pos interneurons
  
  errorbar(1:11, falsePosMeans, falsePosStds)
  hold on
  errorbar(1:11, falseNegMeans, falseNegStds,'r')
  grid on
  xlim([0 12])
  ylim([-2 12])
  xlabel('# of key neurons')
  ylabel('# of false pos/neg interneurons')
  legend('False Pos', 'False Neg')

if ~exist('./results', 'dir')
    mkdir('./results');
end
saveas(gcf, ['./results/FalseNegPos_Rate_', dataset, '.png']);
end