function [matrixCombi] = fun_latin_hypercube(paramCombo, copies)
    %range of search
    %[0.0 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5]
    paramRange = paramCombo;
    %parameters (number)
    paramsList = { paramRange, paramRange, paramRange, paramRange, paramRange}; %// input data: cell array of vectors

    n = numel(paramsList); %// number of vectors
    combinations = cell(1,n); %// pre-define to generate comma-separated list
    [combinations{end:-1:1}] = ndgrid(paramsList{end:-1:1}); %// the reverse order in these two
    %// comma-separated lists is needed to produce the rows of the result matrix in
    %// lexicographical order 
    combinations = cat(n+1, combinations{:}); %// concat the n n-dim arrays along dimension n+1
    combinations = reshape(combinations,[],n);
    %specify how many copies of the combinations you want
    bigMatCombi = repmat(combinations, copies, 1);
    %create the multiplier matrix
    [mrows, ncolumns] = size(bigMatCombi);
    %create matrix multiplier for mx1 matrix with normally distributed
    %random numbers from normal distribution
    randMultiplier = randn(mrows,1);
    %random numbers from negbinom distribution
    %randMultiplier = transpose(nbinrnd(0.2,0.5,1,mrows));
    %muliply matrices and take absolute values of the random
    %multiplier
    bigMatCombi = bigMatCombi.*abs(randMultiplier);
    %bigMatCombi = bigMatCombi.*randMultiplier;
    %test
    matrixCombi = bigMatCombi;

end
 