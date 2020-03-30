function [responseLog,respMean,respSTD] = runTDtaskSurprisal2020( trials, interval, subjects, runName, N, optimism)
%runTDtaskModel 
%is for running the DSHM model of the Walsh & Anderson (2011) task.
%   This function calls TDtaskUOG which does all of the real work.
%   TDtaskUOG implements a single trial of the experiment.
%
%   trials: number of trials of the experiment each participant performs
%   interval: number of trials that data is aggregated across
%   subjects: number of simulated participants
%   runName: unique name for this run, for labelling data files

numOfIntervals = ceil(trials/interval);
responseLog = zeros(numOfIntervals,3,subjects);
include     = true(numOfIntervals,3,subjects);

if nargin < 5
    N = 256
end
if nargin < 6
    opitmism = 30
end

for s = 1:subjects
    
    i = 1; % current interval
    
    for t = 1:trials
        if t == 1
            [memVec,actVec,placeholder,left,conVect,numCorr,numTotal] = TDtaskSurprisal2020(N,optimism);
        else
            [memVec,actVec,placeholder,left,conVect,numCorr,numTotal] = TDtaskSurprisal2020(N,optimism,memVec,actVec,placeholder,left,conVect,numCorr,numTotal);
        end
        if mod(t,interval) == 0
            for c=1:3
                if numTotal(c) > 0
                    responseLog(i,c,s) = numCorr(c) / numTotal(c);
                else
                    % if there are no examples of this choice
                    % in this interval, due not include this interval
                    % in data aggregation (mean, std) for this choice
                    responseLog(i,c,s) = 0;
                    include(i,c,s) = false;
                end
            end
            i = i + 1;
            numCorr = zeros(1,3);
            numTotal = zeros(1,3);
        end
    end
    
    % write virtual subject data to file. File will be runName#.csv
    %   where # is the subject number, i.e., s
    %csvwrite(strcat(runName,num2str(s)), responseLog(:,:,s));
end

respMean = zeros(numOfIntervals,3);
respSTD  = zeros(numOfIntervals,3);
for v = 1:3
    for i = 1:numOfIntervals
        respMean(i,v) = mean(responseLog(i,v,include(i,v,:)));
        respSTD(i,v)  = std(responseLog(i,v,include(i,v,:)));
    end
end

% write the mean subject data to a file. File will be runName"Mean.csv"
csvwrite(strcat(runName,'Mean.csv'), respMean);

% write the standard deviations to a file. File will be runName"Error.csv"
csvwrite(strcat(runName,'Error.csv'), respSTD);

if numOfIntervals < 5
    % write other data to file for significance testing
    % J
    Ji      = zeros(subjects,numOfIntervals);
    J       = zeros(subjects*numOfIntervals,1);
    Ji(:,1) = reshape(responseLog(1,1,:),[subjects,1]);
    J(1:subjects) = Ji(:,1);
    csvwrite(strcat(runName,'J1.csv'),Ji(:,1));

    % T
    Ti      = zeros(subjects,numOfIntervals);
    T       = zeros(subjects*numOfIntervals,1);
    Ti(:,1) = reshape(responseLog(1,2,:),[subjects,1]);
    T(1:subjects) = Ti(:,1);
    csvwrite(strcat(runName,'T1.csv'),Ti(:,1));

    % V
    Vi      = zeros(subjects,numOfIntervals);
    V       = zeros(subjects*numOfIntervals,1);
    Vi(:,1) = reshape(responseLog(1,3,:),[subjects,1]);
    V(1:subjects) = Vi(:,1);
    csvwrite(strcat(runName,'V1.csv'),Vi(:,1));

    % Block
    block      = zeros(subjects*3,numOfIntervals);
    block(:,1) = [reshape(responseLog(1,1,:),[subjects,1]);reshape(responseLog(1,2,:),[subjects,1]);reshape(responseLog(1,3,:),[subjects,1])];
    csvwrite(strcat(runName,'block1.csv'),block(:,1));

    % all other intervals
    for i=2:numOfIntervals
        % range
        j = subjects*(i-1)+1;
        k = subjects*i;

        % J
        Ji(:,i) = reshape(responseLog(i,1,:),[subjects,1]);
        J(j:k)  = Ji(:,i);
        csvwrite(strcat(runName,'J',int2str(i),'.csv'), Ji(:,i));    

        % T
        Ti(:,i) = reshape(responseLog(i,2,:),[subjects,1]);
        T(j:k)  = Ti(:,i);
        csvwrite(strcat(runName,'T',int2str(i),'.csv'), Ti(:,i));    

        % V
        Vi(:,i) = reshape(responseLog(i,3,:),[subjects,1]);
        V(j:k)  = Vi(:,i);
        csvwrite(strcat(runName,'V',int2str(i),'.csv'), Vi(:,i));    

        % Block
        block(:,i) = [reshape(responseLog(i,1,:),[subjects,1]);reshape(responseLog(i,2,:),[subjects,1]);reshape(responseLog(i,3,:),[subjects,1])];
        csvwrite(strcat(runName,'block',int2str(i),'.csv'), block(:,i));    
    end
    csvwrite(strcat(runName,'J.csv'),J);
    csvwrite(strcat(runName,'T.csv'),T);
    csvwrite(strcat(runName,'V.csv'),V);
end

end