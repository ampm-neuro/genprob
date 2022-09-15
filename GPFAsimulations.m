% testing GPFA and HighDim methods [1,BCS] for data simulated according to
% section "Simulation with error floor" in [BCS]
% [1] https://users.ece.cmu.edu/~byronyu/software.shtml
% [BCS] Byron, M.Y., Cunningham, J.P., Santhanam, G., et. al, 2009. 
% Gaussian-process factor analysis for low-dimensional single-trial analysis 
% of neural population activity. In Advances in neural information processing systems (pp. 1881-1888).

%% 1 GENERATE SIMULATED DATA

clear; close all; clc;
nUnits      = 61;  % dimensionality of observational space
nDims       = 5;   % dimensionality of (original) state space
nTrials     = 56;
binWidth    = 20;
trialLength = 150; % trial length in ms (1 bin = 1 ms)

x          = zeros( nDims,  nTrials*trialLength ); % original "neural" trajectories
firingRate = zeros( nUnits, nTrials*trialLength ); % "observables" obtained by linear mapping
spikes     = zeros( nUnits, nTrials*trialLength*binWidth );

% generate sinusoids with different frequency
A = [3 7 17 5 2]; % set sinusoids frequency
B = [1 1 1 1 1];  % set sinusoids amplitude
for iTrial = 1:nTrials
  for iDim = 1:nDims
    x( iDim, ( iTrial - 1 )*trialLength+1:iTrial*trialLength ) = ...
                      B( iDim )*sin( A( iDim )*pi/180*( 1:trialLength ) );
  end
end

% generate diagonal covariance matrix
R = zeros( nUnits, nUnits ); % covariance matrix where diagonal elements are 
                             % independent noise variances of each neuron
for i = 1:nUnits
    R( i, i ) = 0.5;
end
multiplier = 1;
C = multiplier*rand( nUnits, nDims );
d = rand( nUnits, 1 );

% generate observables with linear Gaussian relationship
% via multivariate normal distribution 
for t = 1:nTrials*trialLength
  firingRate( :, t ) = mvnrnd( C*x( :, t ) + d, R );
end
% make firing rates always positive
firingRate = firingRate + max( max( firingRate ) );

% randomly generate 0 and 1 in time bins so that their sum
% correspond to firing rate value in this time bin
for iUnit = 1:nUnits
  for iPos = 1:nTrials*trialLength
    if ( firingRate( iUnit, iPos ) > 0 )
      spikesIndices = randi( [1 binWidth], 1, ceil( firingRate( iUnit, iPos ) ) );
      spikesIndices = spikesIndices + ( iPos - 1 )*binWidth;
    end
    spikesIndices = sort( spikesIndices );
    spikes( iUnit, spikesIndices ) = 1;
  end
end

%% 2 PLOT SIMULATED DATA
figure; 
subplot( 2, 1, 1 );
trialToPlot = 4;
for iPlot = 1:nDims
  plot( x( iPlot, ( trialToPlot - 1 )*trialLength+1:trialToPlot*trialLength ), 'LineWidth', 2 ); hold on;
end
title( 'Fourth trial of original (underlying) data for all dimensions' );

% firing rate
subplot( 2, 1, 2 );
unitsToPlot = 10;
lengthPlot  = trialLength;
trialToPlot = 4;
for iPlot = 1:unitsToPlot
  plot( firingRate( iPlot, ( trialToPlot - 1 )*trialLength+1:trialToPlot*trialLength ), 'LineWidth', 2 ); hold on;
end
title( 'Fourth trial of observations for first ten units' );

%% 3 PREPARING DATA IN PROPER FORMAT FOR GPFA
GPFA_dat( nTrials ).trialId = []; % space preallocation
GPFA_dat( nTrials ).spikes  = []; % space preallocation
for iTrial = 1:nTrials
  GPFA_dat( iTrial ).spikes = zeros( nUnits, trialLength*binWidth );
  GPFA_dat( iTrial ).spikes = spikes( :, ( iTrial - 1 )*trialLength*binWidth + 1:iTrial*trialLength*binWidth );
  disp( [ 'Trial ' num2str( iTrial ) ] );
  GPFA_dat( iTrial ).trialId = iTrial;
end
disp('FINISHED');

%% 4 LAUNCH GPFA ON SIMULATED DATA

runIdx = 1;
% this just for the first run when no results are precomputed
if ( runIdx == 1 )
  try
    mkdir 'mat_results/run001'
  catch
  end
end

method = 'gpfa'; % 'gpfa' 'fa' 'ppca' 'pca'
% Select number of latent dimensions
xDim   = nDims; % number of latent dimensions (the optimal dimensionality 
                % should be found using cross-validation (Section 2) below.

% If using a two-stage method ('fa', 'ppca', or 'pca'), select
% standard deviation (in msec) of Gaussian smoothing kernel. 
kernSD = 40; % NOTE: The optimal kernel width should be found using 
             %       cross-validation (Section 2) below.
            
result = neuralTraj( runIdx, GPFA_dat, 'binWidth', binWidth, 'method', method, 'xDim', xDim );

% Orthonormalize neural trajectories (pp.621-622 of Yu et al., J Neurophysiol, 2009.)
[ estParams, seqTrain ] = postprocess( result, 'kernSD', kernSD );

% Plot neural trajectories in 3D space
plot3D( seqTrain, 'xorth', 'dimsToPlot', 1:3 );

% Plot each dimension of neural trajectories versus time
plotEachDimVsTime( seqTrain, 'xorth', result.binWidth );

%% 5 FULL CROSS-VALIDATION to find
%  - optimal state dimensionality for all methods
%  - optimal smoothing kernel width for two-stage methods
% ========================================================
dat = GPFA_dat;
% Select number of cross-validation folds
numFolds = 4;
runIdx   = 1; 

% Perform cross-validation for different state dimensionalities.
% Results are saved in mat_results/runXXX/, where XXX is runIdx.
for xDim = 2:1:8
  neuralTraj( runIdx, dat, 'method',  'pca', 'binWidth', binWidth, 'xDim', xDim, 'numFolds', numFolds );
  neuralTraj( runIdx, dat, 'method', 'ppca', 'binWidth', binWidth, 'xDim', xDim, 'numFolds', numFolds );
  neuralTraj( runIdx, dat, 'method',   'fa', 'binWidth', binWidth, 'xDim', xDim, 'numFolds', numFolds );
  neuralTraj( runIdx, dat, 'method', 'gpfa', 'binWidth', binWidth, 'xDim', xDim, 'numFolds', numFolds );
end
fprintf( '\n' );

% Plot prediction error versus state dimensionality.
% Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
kernSD = 40; % select kernSD for two-stage methods
plotPredErrorVsDim( runIdx, kernSD );

% Plot prediction error versus kernelSD.
% Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
xDim = nDims; % select state dimensionality
plotPredErrorVsKernSD( runIdx, xDim );

%% 6 PREPARING DATA FOR DataHigh TOOLBOX (https://users.ece.cmu.edu/~byronyu/software/DataHigh/datahigh.html)

nTrials = 56;
% preallocating array of structures for speed
dataForDataHigh( nTrials ).data      = [];
dataForDataHigh( nTrials ).condition = [];
for iTrial = 1:nTrials
  dataForDataHigh( iTrial ).data = GPFA_dat( iTrial ).spikes;
  if( iTrial < 28 )
    dataForDataHigh( iTrial ).condition = '1';
  else
    dataForDataHigh( iTrial ).condition = '2';
  end
end
DataHigh( dataForDataHigh, 'DimReduce' );

%% 7 EXAMPLE FOR USING PREVIOUSLY-FIT MODEL PARAMETERS FOR EXTRACTING NEURAL TRAJECTORIES
result = load( 'mat_results/run001/gpfa_xDim05.mat' );
seqNew = getTrajNewTrials( result, GPFA_dat );
plotEachDimVsTime( seqNew, 'xorth', result.binWidth );