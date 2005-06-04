function PosteriorOddsTable = CompareModels(ModelNames,ModelPriors,type)
% 05-30-2005
%
% type is a string  = LaplaceApproximation
%                   = ModifiedHarmonicMean
% ModelPriors is a m*1 column vector
% ModelNames is m*1 cell array

global oo_

NumberOfModels = size(ModelNames,1);
MarginalLogDensity = zeros(NumberOfModels,1);

% Get the estimates of the (logged) marginal densities
for i = 1:NumberOfModels
    oo_ = load(['ModelNames(i) '_results.mat' ],'oo_');
    try
        eval(['MarginalLogDensity(i) = oo_.MarginalDensity.' type ';']) 
    catch
        if strcmpi(type,'LaplaceApproximation')
            disp(['CompareModels :: I cant''t find the Laplace approximation associated to model ' ModelNames(i,:)])
            return
        elseif strcmpi(type,'ModifiedHarmonicMean')
            disp(['CompareModels :: I cant''t find the modified harmonic mean estimate associated to model ' ModelNames(i,:)])
            return
        end
    end
end

MarginalDensity = exp(MarginalLogDensity);
ConstantOfIntegration  = ModelPriors'*MarginalDensity;
PosteriorProbabilities = ModelPriors.*MarginalDensity / ConstantOfIntegration;
PosteriorOddsTable = PosteriorProbabilities.*PosteriorProbabilities.^(-1);

% Now I display the posterior probabilities:
if NumberOfModels == 2
    disp(' ')
    disp(['Posterior odd ('ModelNames(1) '/' ModelNames(2) ') =  ' num2str(PosteriorOddsTable(1,2))])
    disp(' ')
else
    disp(' ')
    disp(' Posterior probabilities:')
    for i=1:NumberOfModels
            disp([ 'Model ' int2str(i) ' (' ModelNames(i) ') = ' num2str(PosteriorProbabilities(i))])
    end
    disp(' ')
end