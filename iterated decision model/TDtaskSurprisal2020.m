function [memoryVectors,actionVectors,placeholder,left,contextVectors,numCorr,numTotal] = TDtaskSurprisal(N,OPTIMISM,memoryVectors,actionVectors,placeholder,left,contextVectors,numCorr,numTotal)
%TDtaskUOG
%   Implements a DSHM/BEAGLE model of Walsh & Anderson (2011) decision task
%   Uses uncontrained open grams (see Hannagan et al., 2011)
%   as a means of implementing strategy #1 (see below)

% Four possible strategies for the first decision:
%   1. if I press this key, then will some stuff happen after which I win?
%   (wildcard memory query)
%   2. what happens when I press this key? (tree search)
%   3. how can I get to winning? (backwards tree search)
%   4. do I feel good about this button? (context vector/utility estimate)
% 1, as implemented, is a bit similar to 4.

% Note: There are two kinds of memory vectors in BEAGLE models:
%   context vectors and order vectors
%
% While both are implemented by this code, only order vectors are used.
% This code could be made more efficient in a number of ways, including
% commenting out all of the context vector code.

VERBOSE = false; % print out what the model is doing

CONTEXT_WEIGHT = 0; % to what degree the context vectors contribute to decision making
ORDER_WEIGHT = 1; % to what degree the order vectors contribute to decision making

HIGH_REWARD = 0.8; % high reward probability
LOW_REWARD = 0.2; % low reward probability

SURPRISE_WEIGHT = 1.0;

%OPTIMISM = 30; % optimisim coefficient

PUNISH = false; % punishment: true if on; false if off (better if off)

%N = 256; % vector dimensionality
numOfRounds = 5; % number of steps in the experimental task
cLambda = 5; % = lag + 1, i.e, how many steps are kept in working memory
numOfActions = 9; % Start, R J, + -, T V, Good Bad are all the actions.

%moveSeq is the prior sequence of actions stored as row indices of the
%        actionVector and memoryVector matrices
%        e.g., moveSeq = [3, 1, 5] means that the last three actions
%        are rows 3, 1, and 5 (in chronological order) of the actionVector
%        matrix.
moveSeq = zeros(cLambda,1); 

%indices
start   = 1; % start
choiceR = 2; % R, results in cue +
choiceJ = 3; % J, results in either cue, with 50% probability
cue20   = 4; % +, cue indicating a 20% chance of reward
cue80   = 5; % -, cue indicating an 80% chance of reward
choiceT = 6; % T, 20% chance of reward following +
choiceV = 7; % V, 80% chance of reward following -
good    = 8; % good, symbol indicating reward
bad     = 9; % bad, symbol indicating no reward

% indices for numCorr which counts the number of correct decisions
% made for these 3 choices
choiceRJ = 1; % choice between R and J, J is correct
choice20 = 2; % choice between T and V after +, T is correct
choice80 = 3; % choice between T and V after -, V is correct


actionLabel{start}   = '!';
actionLabel{choiceR} = 'R';
actionLabel{choiceJ} = 'J';
actionLabel{cue20}   = '+';
actionLabel{cue80}   = '-';
actionLabel{choiceT} = 'T';
actionLabel{choiceV} = 'V';
actionLabel{good}    = 'Good';
actionLabel{bad}     = 'Bad';

actionL(start)   = '!';
actionL(choiceR) = 'R';
actionL(choiceJ) = 'J';
actionL(cue20)   = '+';
actionL(cue80)   = '-';
actionL(choiceT) = 'T';
actionL(choiceV) = 'V';
actionL(good)    = 'g';
actionL(bad)     = 'b';
actionL(10)      = ' ';

actionPairs = [start; ...
               choiceJ; ...
               choiceR; ...
               cue80; ...
               cue20; ...
               choiceV; ...
               choiceT; ...
               bad; ...
               good; ...
              ];

% initialize vectors if no arguments given
if nargin < 3
    if VERBOSE
        'initializing'
    end
    
    left  = randperm(N);
    
    numCorr  = zeros(1,3);
    numTotal = zeros(1,3);
    
    actionVectors = zeros(numOfActions,N);
    
    % create a random vector for each of the possible actions
    for i=1:numOfActions
        actionVectors(i,:) = realNorm(makeVector(N))';
    end
    actionVectors(numOfActions+1,:) = 0;
    
    % Make the symbol for bad negative good to punish model
    if PUNISH
        actionVectors(bad,:) = -1 * actionVectors(good,:);
    end
    
    % create a random vector to act as the placeholder vector
    placeholder = realNorm(makeVector(N))';
    
    % optimism
    for i=1:numOfActions
        memoryVectors(i,:) = actionVectors(i,:) + OPTIMISM*cconv(placeholder(left),actionVectors(good,:),N);
    end
    memoryVectors(good,:) = actionVectors(good,:);
    memoryVectors(bad,:) = actionVectors(bad,:);
    
    % optimism
    for i=1:numOfActions
        contextVectors(i,:) = actionVectors(good,:);
    end
    contextVectors(bad,:) = actionVectors(bad,:);

    %memStr = cell(numOfActions,1);
    %contextStr = memStr;
    
end

% inverse permutation of the "left" permutation
%   can be used for retrieval (decoding)
%unleft(left) = 1:N;

signal = start;

%predictCount = zeros(numOfActions,1);
predictRound = zeros(numOfRounds,1) + 10;
for round=1:numOfRounds % iterate on rounds of play
    
    if VERBOSE
        actionLabel(signal)
    end
    
    % push vectors back in matrix if we've run out of 'room'
    if round > cLambda
        for i=2:cLambda
            percepts(i - 1,:) = percepts(i,:);
            moveSeq(i - 1)    = moveSeq(i);
        end
        percepts(cLambda,:) = actionVectors(signal,:);
        moveSeq(cLambda)    = signal;
    else % for the first few turns, the percept matrix grows to cLambda
        percepts(round,:) = actionVectors(signal,:);
        moveSeq(round)    = signal;
    end
    
    % number of actions stored in the moveSeq
    seqLength = min(round,cLambda);
  
    % update the appropriate memory vectors (as indicated by moveSeq)
    if round == numOfRounds
        for i = 1:seqLength
            % get the information we need to update the memory vectors
            %conceptStr     = getUOGstr(actionL(moveSeq(1:seqLength)),'*',cLambda,i);
            percepts2      = percepts;
            percepts2(i,:) = placeholder;
            % surprisal weighting here
            %opposites      = actionVectors(actionPairs(moveSeq(1:seqLength)),:);
            if VERBOSE
                disp('Current sequence of actions:');
                disp(actionL(moveSeq(1:seqLength)));
                disp('Predictions:');
                disp(actionL(predictRound));
            end
            %predicts       = predictCount * opposites;
            %percepts2      = percepts2 - predicts;
            surpriseSeq                         = predictRound;
            surpriseSeq(surpriseSeq == moveSeq) = 10; % blank
            surpriseSeq(i)                      = 10; % blank / zeros
            if VERBOSE
                disp('Surprises:');
                disp(actionL(surpriseSeq));
            end
            surpriseVec    = actionVectors(surpriseSeq,:);
            concept        = hdmUOG(percepts2 - (SURPRISE_WEIGHT*surpriseVec),i,left);
            
            memoryVectors(moveSeq(i),:) = memoryVectors(moveSeq(i),:) + concept;
            %memStr(moveSeq(i)) = strcat(memStr(moveSeq(i)),'_',conceptStr);
            
            % update the context vector with the sum of percepts
            %   from percept i+1 to the last percept
            if i < seqLength
                contextVectors(moveSeq(i),:) = contextVectors(moveSeq(i),:) + sum(percepts(i+1:seqLength,:));
                %for j = i+1:seqLength
                %    contextStr(moveSeq(i)) = strcat(contextStr(moveSeq(i)),'_',actionL(moveSeq(j)));
                %end
            end
        end
    end
    
    % construct the context probe
    goodVector = actionVectors(good,:);
    %badVector  = actionVectors(bad,:);
    
    % construct the order probe query sequence
    % which consists of cLambda items
    % such that the last two items are 'placeholder' and 'good'
    % and the preceding items are the preceding cLambda - 2 moves in
    % moveSeq
    
    queryStart  = 1; %max(1,seqLength - 2);
    queryLength = round; %min(round, cLambda - 2);
    
    queryVectors(1:queryLength,:) = percepts(queryStart:seqLength,:);
    queryLabels(1:queryLength)    = actionL(moveSeq(queryStart:seqLength));
    
    p = queryLength+1;
    queryVectors(p,:) = placeholder;
    queryLabels(p) = 'x';
    % Putting good into the query yields weird results
    % such as allowing probe strings that don't contain good
    %queryVectors(queryLength+2,:) = goodVector;
    %queryLabels(queryLength+2) = actionL(good);

    % Using UOG for the probe works terribly because it
    % more heavily weights the unconditional effect of choosing T
    % and so it learns to avoid T entirely
    %probe = hdmUOG(queryVectors,p,left);

    % create the order probe by getting all combinations of the query seq
    %probe = hdmNgram(queryVectors,p,left);
    probe = getCombos(queryVectors,placeholder,cLambda-1,p,left);
    %probeStr = getCombosStr(queryLabels,'?g',cLambda-1,p);
    
    % this doesn't quite work, for reasons unclear to me
    % it fails to learn to chose T quite badly
    % possibly due to an interaction effect with how optimism is
    % implemented
    %queryVectors(p+1,:) = actionVectors(good,:);
    %probe    = queryVectors(1,:);
    %probeStr = queryLabels(1);
    %queryLabels(p)   = '?';
    %queryLabels(p+1) = actionL(good);
    %for item=2:p
    %    probe    = cconv(probe(left),queryVectors(item,:),N);
    %    probeStr = [probeStr,queryLabels(item)];
    %end
    
    % prediction probe / probe for decoding (presently unused)
    probeD = probe;

    % we're going to add *g (placeholder good) to the probe 
    % for the purpose of implementing optimism
    probe = probe + placeholder;
    %probeStr = strcat(probeStr,'g_','?g');
    probe = cconv(probe(left),goodVector,N);

    if all(probe == 0)
        error('empty probe');
    end
    
    % calculate the similarities of the probe to the memoryVectors
    orderSimilarities    = zeros(numOfActions,1);
    contextSimilarities  = zeros(numOfActions,1);
    decodeSimilarities   = zeros(numOfActions,numOfActions);
    printOut = cell(numOfActions,5);
    
    for i = 1:numOfActions
        contextSimilarities(i) = vectorCosine(goodVector, contextVectors(i,:)); %- vectorCosine(badVector, contextVectors(i,:));
        orderSimilarities(i)   = vectorCosine(probe, memoryVectors(i,:));
        decoded = ccorr(probeD(left),memoryVectors(i,:));
        
        for j = 1:numOfActions
            decodeSimilarities(i,j) = vectorCosine(decoded,actionVectors(j,:));
        end
        
        %printOut(i,1) = actionLabel(i);
        %printOut(i,2) = {orderSimilarities(i)};
        %printOut(i,3) = {contextSimilarities(i)};
        %printOut(i,4) = memStr(i);
        %printOut(i,5) = contextStr(i);
    end
    
    %if VERBOSE
    %    probeStr
    %    printOut
    %end
    
    %decodeSimilarities
    % choices available each round
    if round == 1
        firstopt = choiceR;
        lastopt  = choiceJ;
    elseif round == 2
        firstopt = cue20;
        lastopt  = cue80;
    elseif round == 3
        firstopt = choiceT;
        lastopt  = choiceV;
    elseif round == 4
        firstopt = good;
        lastopt  = bad;
    end
    
    % make a choice
    evaluations = ORDER_WEIGHT*orderSimilarities + CONTEXT_WEIGHT*contextSimilarities;
    [~,predicted] = max(evaluations(firstopt:lastopt));
    predicted = predicted + firstopt - 1;
    [~,prediction2] = max(decodeSimilarities(predicted,:));
    if round+2 <= numOfRounds
        predictRound(round+2) = prediction2;
    end
    %predictCount(prediction2) = predictCount(prediction2)+1;
    
    %if VERBOSE
    %    strcat('computer picks...',actionLabel(predicted))
    %    strcat('and predicts...',actionLabel(prediction2))
    %end
    
    
    
    % signal for next round, implements fitness function/problem structure
    if round == 1
        signal = predicted;
    elseif round == 2
        %increment total number of choices made count
        numTotal(choiceRJ) = numTotal(choiceRJ) + 1;
        
        if moveSeq(2) == choiceR
            signal = cue20;
        else
            % increment correct choice count
            numCorr(choiceRJ) = numCorr(choiceRJ) + 1;
            
            % randomly determine the cue
            if rand(1) <= 0.5
                signal = cue20;
            else
                signal = cue80;
            end
        end
    elseif round == 3
        signal = predicted;
    elseif round == 4
        signal = bad;
        if moveSeq(3) == cue20
            numTotal(choice20) = numTotal(choice20) + 1;
            if moveSeq(4) == choiceT
                % increment correct choice count
                numCorr(choice20) = numCorr(choice20) + 1;
                
                % small chance of reward
                if rand(1) <= LOW_REWARD
                    signal = good;
                end
            end
        elseif moveSeq(3) == cue80
            numTotal(choice80) = numTotal(choice80) + 1;
            if moveSeq(4) == choiceV
                % increment correct choice count
                numCorr(choice80) = numCorr(choice80) + 1;
                
                % big chance of reward
                if rand(1) <= HIGH_REWARD
                    signal = good;
                end
            end
        end
    end
end % end iteration on rounds of play 

if VERBOSE
    actionLabel(signal)
end
end % end program