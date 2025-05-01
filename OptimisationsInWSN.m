clear all;
clc;

% Simulation Parameters
simArea = 1000;
noOfNodes = 600;
transmissionRange = 25;
interferenceRange = 30;
transmissionRate = 4;
deliveryRate = 8;
initialEnergyLevel = 3;
numberOfInteractions = 300;
maliciousPercentage = 30;
noOfMaliciousNodes = ceil(maliciousPercentage / 100 * noOfNodes);

% Energy Model Parameters (Typical WSN Values)
Eelec = 50e-9;          % Energy to run transmitter/receiver circuitry (Joules/bit)
epsilon_fs = 10e-12;    % Free space model amplifier energy (Joules/bit/m^2)
epsilon_mp = 0.0013e-12; % Multi-path model amplifier energy (Joules/bit/m^4)
d0 = sqrt(epsilon_fs/epsilon_mp); % Distance threshold (meters)
packetSize = 2000;      % Packet size in bits

% ACO Parameters
alpha_aco = 1;          % Importance of pheromone
beta_aco = 2;           % Importance of heuristic (trust, energy)
rho = 0.02;              % Pheromone evaporation rate

% Node Deployment
nodeXLoc = rand(1, noOfNodes) * simArea;
nodeYLoc = rand(1, noOfNodes) * simArea;

% Initialize node properties and energy
nodeEnergy = ones(1, noOfNodes) * initialEnergyLevel;
nodeInteractionCount = zeros(1, noOfNodes);

% Initialize binary node properties
newNode = containers.Map();
for i = 1:noOfNodes
    newNode(num2str(i)) = randi([0, 1]);
end

% Initialize recommendations
recommendedNodes = containers.Map();
for i = 1:noOfNodes
    recommendedNodes(num2str(i)) = randi([1, noOfNodes], 1, noOfNodes);
end

% Calculate neighbor relationships
neighborNode = zeros(noOfNodes, noOfNodes);
for i = 1:noOfNodes
    for j = 1:noOfNodes
        if i ~= j
            distance = sqrt((nodeXLoc(i) - nodeXLoc(j))^2 + (nodeYLoc(i) - nodeYLoc(j))^2);
            if distance <= transmissionRange
                neighborNode(i, j) = 1;
            end
        end
    end
end

% Identify malicious nodes
maliciousNodes = randperm(noOfNodes, noOfMaliciousNodes);
isMalicious = zeros(1, noOfNodes);
isMalicious(maliciousNodes) = 1;

% Initialize Interaction Tracking
totalInteractions = 0;
successfulInteractions = 0;
reciprocatedInteractions = 0;
totalResponseTime = zeros(1, noOfNodes);

% Initialize trust matrices
deliveryRatio = ones(noOfNodes, noOfNodes);
compatibility = ones(noOfNodes, noOfNodes);
cooperativeness = ones(noOfNodes, noOfNodes);
trustValue = ones(noOfNodes, noOfNodes) * 3;
pheromone = ones(noOfNodes, noOfNodes); % Initial pheromone for ACO

% Apply malicious node properties
for i = maliciousNodes
    trustValue(:, i) = trustValue(:, i) * 0.7; % Lower initial trust for malicious nodes
end

% Arrays to track metrics
trustValues = zeros(1, numberOfInteractions);
deliveryRatios = zeros(1, numberOfInteractions);
successfulInteraction_mat = zeros(1, numberOfInteractions);
avgTrustHistory = zeros(1, numberOfInteractions);
movingAvgTrust = zeros(1, numberOfInteractions);
interactionHistory = zeros(noOfNodes, noOfNodes);
successHistory = zeros(noOfNodes, noOfNodes);

% Helper functions
function similarity = jaccardSimilarity(A, B)
    % For binary vectors
    intersection = sum(A & B);
    union = sum(A | B);
    if union == 0
        similarity = 0;
    else
        similarity = intersection / union;
    end
end

function DR = deliveryRatioCalc(transmissionRange, interferenceRange, deliveryRate, transmissionRate)
    DR = ((transmissionRange / interferenceRange)^2) * (deliveryRate / transmissionRate);
end

function isNewNodeFlag = isNewNode(a, newNode)
    if isKey(newNode, num2str(a)) && newNode(num2str(a)) == 1
        isNewNodeFlag = true;
    else
        isNewNodeFlag = false;
    end
end

function isRec = isRecommended(a, b, recommendedNodes)
    if isKey(recommendedNodes, num2str(a))
        recom = recommendedNodes(num2str(a));
        isRec = any(recom == b);
    else
        isRec = false; 
    end
end

% Initialize trust values based on delivery ratio
for i = 1:noOfNodes
    for j = 1:noOfNodes
        if i ~= j
            deliveryRatio(i, j) = deliveryRatioCalc(transmissionRange, interferenceRange, deliveryRate, transmissionRate);
        end
    end
end

% Calculate compatibility values based on node properties
for i = 1:noOfNodes
    for j = 1:noOfNodes
        if i ~= j
            % Binary feature vectors for similarity comparison
            featuresA = [nodeEnergy(i) > initialEnergyLevel/2, isNewNode(i, newNode), isMalicious(i) == 0];
            featuresB = [nodeEnergy(j) > initialEnergyLevel/2, isNewNode(j, newNode), isMalicious(j) == 0];
            
            compatibility(i, j) = jaccardSimilarity(featuresA, featuresB);
            % Simulate interaction between node i and node j
            totalInteractions = totalInteractions + 1;
            
            % Assume interaction results in successful transmission
            successfulInteractions = successfulInteractions + randi([0,1]);
            
            % Track response time (random value for demonstration)
            responseTime = randi([2, 4]); % Placeholder for actual response time
            totalResponseTime(i) = totalResponseTime(i) + responseTime;
            
            % Check reciprocity (assume all interactions are reciprocated for demonstration)
            reciprocatedInteractions = reciprocatedInteractions + randi([0,1]);
            
            % Calculate Success Rate
            successRate = successfulInteractions / totalInteractions;

            % Calculate Reciprocity Rate
            reciprocityRate = reciprocatedInteractions / totalInteractions;

            % Calculate Average Response Time per Node
            averageResponseTime = mean(totalResponseTime)/noOfNodes; 

            % Calculate Cooperativeness
            cooperativeness(i, j) = max(1, successRate + reciprocityRate + averageResponseTime);
        end
    end
end

% Calculate initial node trust values for clustering
nodeTrustValues = zeros(1, noOfNodes);
for node = 1:noOfNodes
    neighbors = find(neighborNode(node,:) == 1);
    if ~isempty(neighbors)
        nodeTrustValues(node) = mean(trustValue(node, neighbors));
    end
end

% =========================== CLUSTERING USING RFO ===========================
noOfClusters = round(0.05 * noOfNodes); % 5% of nodes as cluster heads
clusterHeads = zeros(1, noOfClusters);
clusters = cell(1, noOfClusters);
clusterMembership = zeros(1, noOfNodes);

% Select initial cluster heads with high trust and not malicious if possible
potentialHeads = find(isMalicious == 0);
if length(potentialHeads) >= noOfClusters
    [~, sortedIndices] = sort(nodeTrustValues(potentialHeads), 'descend');
    clusterHeads = potentialHeads(sortedIndices(1:noOfClusters));
else
    % If not enough non-malicious nodes, use the best available
    [~, sortedIndices] = sort(nodeTrustValues, 'descend');
    clusterHeads = sortedIndices(1:noOfClusters);
end

% Assign nodes to clusters based on distance
for i = 1:noOfNodes
    if ~ismember(i, clusterHeads)
        dists = zeros(1, noOfClusters);
        for k = 1:noOfClusters
            dists(k) = sqrt((nodeXLoc(i) - nodeXLoc(clusterHeads(k)))^2 + (nodeYLoc(i) - nodeYLoc(clusterHeads(k)))^2);
        end
        [~, idx] = min(dists);
        clusters{idx} = [clusters{idx}, i];
        clusterMembership(i) = idx;
    else
        % Cluster head belongs to its own cluster
        idx = find(clusterHeads == i);
        clusters{idx} = [clusters{idx}, i];
        clusterMembership(i) = idx;
    end
end

% ========================== TRUST EVALUATION FUNCTIONS ==========================
function success = TrustEvaluation(a, b, trustValue, newNode, recommendedNodes, compatibility, cooperativeness, deliveryRatio, clusterHeads, clusterMembership, isMalicious, nodeInteractionCount, neighborNode)
    % Apply malicious penalty if detected
    if isMalicious(b) == 1 && nodeInteractionCount(b) > 5
        maliciousPenalty = 0.7;
    else
        maliciousPenalty = 0;
    end
    
    % Select appropriate trust evaluation method
    if trustValue(a, b) > 0 && nodeInteractionCount(a) > 0 && nodeInteractionCount(b) > 0
        % Direct trust for known nodes
        success = TrustEvaluationDirect(a, b, trustValue, compatibility, cooperativeness, deliveryRatio, clusterMembership, maliciousPenalty, clusterHeads);
    elseif isNewNode(b, newNode)
        % Absolute trust for new nodes
        success = TrustEvaluationAbsolute(a, b, trustValue, clusterMembership, compatibility, cooperativeness, deliveryRatio, maliciousPenalty, clusterHeads);
    else
        % Recommendation-based trust for other cases
        success = TrustEvaluationRecommendation(a, b, trustValue, recommendedNodes, compatibility, cooperativeness, deliveryRatio, clusterHeads, clusterMembership, maliciousPenalty, neighborNode);
    end
end

function success = TrustEvaluationDirect(a, b, trustValue, compatibility, cooperativeness, deliveryRatio, clusterMembership, maliciousPenalty, clusterHeads)
    % Calculate direct trust components
    compatibility_a_b = compatibility(a, b);
    cooperativeness_a_b = cooperativeness(a, b);
    deliveryRatio_a_b = deliveryRatio(a, b);
    
    % Add cluster membership bonus
    clusterBonus = 0;
    if clusterMembership(a) == clusterMembership(b)
        clusterBonus = 2.5;
    end
    if b == clusterHeads(clusterMembership(b)) || trustValue(clusterHeads(clusterMembership(b)), b) > 5.0
        clusterBonus = clusterBonus + 1;
    end
    
    % Calculate final trust score
    trust_a_to_b = compatibility_a_b + cooperativeness_a_b + deliveryRatio_a_b + clusterBonus - maliciousPenalty;
    trustValue(a, b) = trust_a_to_b;
    
    % Decision threshold
    threshold = 5.0;
    success = (trust_a_to_b >= threshold);
end

function success = TrustEvaluationAbsolute(a, b, trustValue, clusterMembership, compatibility, cooperativeness, deliveryRatio, maliciousPenalty, clusterHeads)
    % For new nodes with no history
    compatibility_b = compatibility(a, b);
    cooperativeness_b = cooperativeness(a, b);
    deliveryRatio_b = deliveryRatio(a, b);

    % Cluster membership bonus
    clusterBonus = 0;
    if clusterMembership(a) == clusterMembership(b)
        clusterBonus = 2;
    end
    if b == clusterHeads(clusterMembership(b)) || trustValue(clusterHeads(clusterMembership(b)), b) > 5.0
        clusterBonus = clusterBonus + 1;
    end
    
    % Calculate trust score
    trust_score = compatibility_b + cooperativeness_b + deliveryRatio_b + clusterBonus - maliciousPenalty;
    trustValue(a, b) = trust_score;
    
    % Lower threshold for new nodes
    threshold = 4.0;
    success = (trust_score >= threshold);
end

function success = TrustEvaluationRecommendation(a, b, trustValue, recommendedNodes, compatibility, cooperativeness, deliveryRatio, clusterHeads, clusterMembership, maliciousPenalty, neighborNode)
    % Get cluster heads of both nodes
    ch_a = clusterHeads(clusterMembership(a));
    ch_b = clusterHeads(clusterMembership(b));
    
    % Start with cluster-based trust
    cluster_trust = 0;
    if ch_a ~= a && ch_a ~= b
        % Trust from cluster head
        if trustValue(ch_a, b) > 0
            cluster_trust = trustValue(ch_a, b);
        elseif ch_a == ch_b || trustValue(ch_a, ch_b) > 0
            cluster_trust = 2.0; % Default trust between cluster heads
        end
    end

    if ch_b == b || trustValue(ch_b, b) > 5.0
        cluster_trust = cluster_trust + 2.0; % Bonus for cluster head
    end
    
    % Add recommendation-based trust
    recom_trust = 0;
    count = 0;
    
    % Direct recommendations
    if isRecommended(a, b, recommendedNodes)
        recom_trust = recom_trust + 1.5;
        count = count + 1;
    end
    
    % Check recommendations from neighbors
    neighbors = find(neighborNode(a,:) == 1);
    for neighbor = neighbors
        if neighbor ~= a && neighbor ~= b && isRecommended(neighbor, b, recommendedNodes)
            recom_trust = recom_trust + 0.8 * trustValue(a, neighbor);
            count = count + 1;
        end
    end
    
    % Calculate final trust score
    if count > 0
        final_trust = cluster_trust + recom_trust / count - maliciousPenalty;
    else
        final_trust = cluster_trust - maliciousPenalty;
    end
    trustValue(a, b) = final_trust;
    
    % Threshold for recommendation-based trust
    threshold = 3.0;
    success = (final_trust >= threshold);
end

% ========================== RED FOX OPTIMIZATION ==========================
% Initial RFO for cluster heads
maxIterations = 30;
for iter = 1:maxIterations
    % Update node trust values
    for node = 1:noOfNodes
        neighbors = find(neighborNode(node,:) == 1);
        if ~isempty(neighbors)
            nodeTrustValues(node) = mean(trustValue(node, neighbors));
        end
    end
    
    % Optimize cluster heads using RFO
    for k = 1:noOfClusters
        members = clusters{k};
        if length(members) > 1
            % Calculate combined trust and energy score
            memberScores = nodeTrustValues(members) + nodeEnergy(members);
            memberScores(isMalicious(members) == 1) = memberScores(isMalicious(members) == 1) * 0.5;
            [~, bestIdx] = max(memberScores);
            bestMember = members(bestIdx);
            
            % Replace cluster head if better candidate found
            currentScore = nodeTrustValues(clusterHeads(k)) + nodeEnergy(clusterHeads(k));
            if isMalicious(clusterHeads(k)) == 1
                currentScore = currentScore * 0.5;
            end
            
            if bestMember ~= clusterHeads(k) && memberScores(bestIdx) > currentScore
                clusterHeads(k) = bestMember;
            end
        end
    end
    
    % Reassign nodes based on new cluster heads
    clusters = cell(1, noOfClusters);
    for i = 1:noOfNodes
        if ~ismember(i, clusterHeads)
            dists = zeros(1, noOfClusters);
            for k = 1:noOfClusters
                dists(k) = sqrt((nodeXLoc(i) - nodeXLoc(clusterHeads(k)))^2 + (nodeYLoc(i) - nodeYLoc(clusterHeads(k)))^2);
            end
            [~, idx] = min(dists);
            clusters{idx} = [clusters{idx}, i];
            clusterMembership(i) = idx;
        else
            idx = find(clusterHeads == i);
            clusters{idx} = [clusters{idx}, i];
            clusterMembership(i) = idx;
        end
    end
end

% ========================== RUN SIMULATION ==========================
% Helper function for ACO routing
function nextNode = selectNextNode(currentNode, neighbors, pheromone, trustValue, nodeEnergy, initialEnergyLevel, nodeXLoc, nodeYLoc, alpha_aco, beta_aco)

    if isempty(neighbors)
        nextNode = currentNode;
        return;
    end

    % Calculate desirability
    tau = pheromone(currentNode, neighbors); % Pheromone strength
    trust = trustValue(currentNode, neighbors) / 10; % Normalize trust to 0-1
    distance = sqrt((nodeXLoc(currentNode) - nodeXLoc(neighbors)).^2 + (nodeYLoc(currentNode) - nodeYLoc(neighbors)).^2);
    eta = 1 ./ (distance + 1e-6); % Inverse distance to prefer closer nodes

    energyFactor = nodeEnergy(neighbors) / initialEnergyLevel; % Healthy nodes preferred

    % New combined desirability
    prob = (tau.^1.5) .* (trust.^2.0) .* (eta.^1.0) .* (energyFactor.^1.2);

    prob = prob / sum(prob); % Normalize to probabilities

    % Roulette wheel selection
    cumProb = cumsum(prob);
    r = rand();
    nextNode = neighbors(find(cumProb >= r, 1, 'first'));
end

% Main simulation loop
for interaction = 1:numberOfInteractions
    % Select random source and destination
    nodeA = randi(noOfNodes);
    nodeB = randi(noOfNodes);
    while nodeB == nodeA
        nodeB = randi(noOfNodes);
    end
    
    currentCluster = clusterMembership(nodeA);
    destinationCluster = clusterMembership(nodeB);
    path = [];
    maxHops = 20; % to avoid infinite loop
    
    if currentCluster == destinationCluster
        % Intra-cluster routing
        currentNode = nodeA;
        path = [currentNode];
    
        while currentNode ~= nodeB && length(path) < maxHops
            neighbors = find(neighborNode(currentNode, :) == 1 & clusterMembership == currentCluster);
            nextNode = selectNextNode(currentNode, neighbors, pheromone, trustValue, nodeEnergy, initialEnergyLevel, nodeXLoc, nodeYLoc, alpha_aco, beta_aco);
            
            if isempty(nextNode) || nextNode == currentNode
                break; % No forward progress
            end
            
            path = [path, nextNode];
            currentNode = nextNode;
        end
    
    else
        % Inter-cluster routing
    
        % Phase 1: NodeA to Source Cluster Head
        currentNode = nodeA;
        sourceCH = clusterHeads(currentCluster);
        path = [currentNode];
        
        while currentNode ~= sourceCH && length(path) < maxHops
            neighbors = find(neighborNode(currentNode, :) == 1 & clusterMembership == currentCluster);
            nextNode = selectNextNode(currentNode, neighbors, pheromone, trustValue, nodeEnergy, initialEnergyLevel, nodeXLoc, nodeYLoc, alpha_aco, beta_aco);
            
            if isempty(nextNode) || nextNode == currentNode
                break;
            end
            
            path = [path, nextNode];
            currentNode = nextNode;
        end
    
        % Phase 2: Source Cluster Head to Destination Cluster Head
        destinationCH = clusterHeads(destinationCluster);
        if currentNode == sourceCH
            while currentNode ~= destinationCH && length(path) < maxHops
                clusterHeadNeighbors = find(neighborNode(currentNode, :) == 1 & ismember(1:noOfNodes, clusterHeads));
                nextNode = selectNextNode(currentNode, clusterHeadNeighbors, pheromone, trustValue, nodeEnergy, initialEnergyLevel, nodeXLoc, nodeYLoc, alpha_aco, beta_aco);
    
                if isempty(nextNode) || nextNode == currentNode
                    break;
                end
    
                path = [path, nextNode];
                currentNode = nextNode;
            end
        end
    
        % Phase 3: Destination Cluster Head to nodeB
        if currentNode == destinationCH
            while currentNode ~= nodeB && length(path) < maxHops
                neighbors = find(neighborNode(currentNode, :) == 1 & clusterMembership == destinationCluster);
                nextNode = selectNextNode(currentNode, neighbors, pheromone, trustValue, nodeEnergy, initialEnergyLevel, nodeXLoc, nodeYLoc, alpha_aco, beta_aco);
    
                if isempty(nextNode) || nextNode == currentNode
                    break;
                end
    
                path = [path, nextNode];
                currentNode = nextNode;
            end
        end
    end

    % Calculate energy consumption along the path
    for p = 1:length(path)-1
        sender = path(p);
        receiver = path(p+1);
        distance = sqrt((nodeXLoc(sender) - nodeXLoc(receiver))^2 + (nodeYLoc(sender) - nodeYLoc(receiver))^2);
        
        % Energy model for transmission/reception
        if distance < d0
            Etx = packetSize * Eelec + packetSize * epsilon_fs * distance^2;
        else
            Etx = packetSize * Eelec + packetSize * epsilon_mp * distance^4;
        end
        
        Erx = packetSize * Eelec;
        
        nodeEnergy(sender) = max(0, nodeEnergy(sender) - Etx);
        nodeEnergy(receiver) = max(0, nodeEnergy(receiver) - Erx);
    end
 
    % Evaluate trust between source and destination
    successfulInteraction_mat(interaction) = TrustEvaluation(nodeA, nodeB, trustValue, newNode, recommendedNodes, compatibility, cooperativeness, deliveryRatio, clusterHeads, clusterMembership, isMalicious, nodeInteractionCount, neighborNode);
    
    % Update interaction counters
    nodeInteractionCount(nodeA) = nodeInteractionCount(nodeA) + 1;
    nodeInteractionCount(nodeB) = nodeInteractionCount(nodeB) + 1;
    interactionHistory(nodeA, nodeB) = interactionHistory(nodeA, nodeB) + 1;
    
    if successfulInteraction_mat(interaction) == 1
        successHistory(nodeA, nodeB) = successHistory(nodeA, nodeB) + 1;
    end
    
    % Calculate success rate for the pair
    successRate = 0;
    if interactionHistory(nodeA, nodeB) > 0
        successRate = successHistory(nodeA, nodeB) / interactionHistory(nodeA, nodeB);
    end
    
    % Update network parameters
    drate = deliveryRate + randi([-1, 1]);
    trate = transmissionRate + randi([-1, 1]);
    
    % Update energy based on interaction result
    if successfulInteraction_mat(interaction) == 1
        nodeEnergy(nodeA) = min(nodeEnergy(nodeA) + 2, 10);
        nodeEnergy(nodeB) = min(nodeEnergy(nodeB) + 2, 10);
        
        % Update pheromone for ACO
        % Dynamic pheromone update after success
        pheromone = (1 - rho) * pheromone; % Evaporation
        
        pathTrustBoost = mean(mean(trustValue(path(1:end-1), path(2:end)))); % Average trust along path
        
        for p = 1:length(path)-1
            sender = path(p);
            receiver = path(p+1);
            pheromone(sender, receiver) = pheromone(sender, receiver) + 0.4 * (trustValue(sender, receiver)/10);
        end

        
        % Propagate positive trust through cluster
        clusterA = clusterMembership(nodeA);
        for member = clusters{clusterA}
            if member ~= nodeA && member ~= nodeB && isMalicious(member) == 0
                % Increase trust gradually with dampening factor
                currentTrust = trustValue(member, nodeB);
                trustValue(member, nodeB) = currentTrust + 1.2 * (1 - (currentTrust / 10));
            end
        end
    else
        nodeEnergy(nodeA) = max(nodeEnergy(nodeA) - 0.5, 0);
    end
    
    % Calculate new delivery ratio with energy consideration
    energyFactor = nodeEnergy(nodeA) / initialEnergyLevel; % Node health
    pathPenalty = 1 / (length(path) + 1); % Longer paths are penalized
    trustFactor = trustValue(nodeA, nodeB);
    
    DR = ((transmissionRange / interferenceRange)^2) * (drate / trate) * energyFactor * pathPenalty * trustFactor;

    
    % Apply exponential smoothing to delivery ratio
    alpha = 0.6; % Smoothing factor
    if interaction > 1
        DR = alpha * DR + (1 - alpha) * deliveryRatios(interaction - 1);
    end
    
    deliveryRatios(interaction) = DR;
    deliveryRatio(nodeA, nodeB) = DR;
    
    % Update trust values
    % Base Increment
    if successfulInteraction_mat(interaction)
        baseIncrement = 0.6 + 0.4 * (successRate); % Successful and stable = more bonus
    else
        baseIncrement = -0.01; % Failures punish more if failure trend
    end

    % Modify by path quality
    pathQualityBoost = min(1, 3 / length(path)); % Good if shorter path (max boost 1)
    
    % Modify by energy status
    energyBoost = nodeEnergy(nodeA) / initialEnergyLevel; % Healthier node = higher boost
    
    % Modify by trust momentum
    trustMomentum = (trustValue(nodeA, nodeB) - 5) / 5; % Positive if >5, Negative if <5
    
    % Final Trust Change Rate
    trustChangeRate = baseIncrement * (1 + 0.4 * pathQualityBoost + 0.3 * energyBoost + 0.3 * trustMomentum);
    
    % If malicious behavior suspected, cap the increase
    if isMalicious(nodeA) || isMalicious(nodeB)
        trustChangeRate = min(trustChangeRate, 0.4); % Trust rises slowly for malicious nodes
    end
    
    % Apply update with boundaries
    trustValue(nodeA, nodeB) = max(0, min(20, trustValue(nodeA, nodeB) + trustChangeRate));
    trustValues(interaction) = trustValue(nodeA, nodeB);


    % Scale increment based on history to promote stability
    if interactionHistory(nodeA, nodeB) > 5
        successTrend = successRate > 0.5; % Check if generally successful
        if successTrend
            trustChangeRate = baseIncrement * 1.5; % Adjust change rate for stable relationships
        else
            trustChangeRate = baseIncrement;
        end
    else
        trustChangeRate = baseIncrement;
    end

    neighborA = find(neighborNode(nodeA, :) == 1);
    neighborB = find(neighborNode(nodeB, :) == 1);
    
    for n = neighborA
        trustValue(n, nodeB) = min(10, trustValue(n, nodeB) + 0.2);
    end
    
    for n = neighborB
        trustValue(n, nodeA) = min(10, trustValue(n, nodeA) + 0.2);
    end

    
    % Apply cluster-based trust bonus
    if clusterMembership(nodeA) == clusterMembership(nodeB)
        trustChangeRate = trustChangeRate * 2; % Cluster bonus
    end
    
    % Apply trust update
    trustValue(nodeA, nodeB) = max(0, min(10, trustValue(nodeA, nodeB) + trustChangeRate * 3 + nodeEnergy(nodeA)/initialEnergyLevel));
    trustValues(interaction) = trustValue(nodeA, nodeB);
    
    % Update trust for nodes in the path
    if successRate > 0.7
        trustBoost = 0.5; % bigger boost
    else
        trustBoost = 0.2; % small cautious boost
    end
    for p = 1:length(path)-1
        sender = path(p);
        receiver = path(p+1);
        trustValue(sender, receiver) = min(10, trustValue(sender, receiver) + trustBoost);
    end

    alpha_smooth = 0.6; % very low alpha for no spikes in trust of the system for stability
    if interaction > 1
    trustValues(interaction) = alpha_smooth * trustValues(interaction) + (1 - alpha_smooth) * trustValues(interaction-1);
    trustValue(nodeA, nodeB) = trustValues(interaction);
    end
    
    % Calculate moving average for trust visualization
    windowSize = min(10, interaction);
    movingAvgTrust(interaction) = mean(trustValues(max(1, interaction-windowSize+1):interaction));
    
    % Periodically update cluster heads
    if mod(interaction, 20) == 0
        % Update node trust values
        for node = 1:noOfNodes
            neighbors = find(neighborNode(node,:) == 1);
            if ~isempty(neighbors)
                nodeTrustValues(node) = mean(trustValue(node, neighbors));
            end
        end
        
        % Red Fox Optimization for cluster head selection
        for k = 1:noOfClusters
            members = clusters{k};
            for i = 1:length(members)
                for j = 1:length(members)
                    if i ~= j && isMalicious(members(i)) == 0 && isMalicious(members(j)) == 0
                        trustValue(members(i), members(j)) = min(10, trustValue(members(i), members(j)) + 0.2);
                    end
                end
            end
            
            % Check if current cluster head has low energy
            if nodeEnergy(clusterHeads(k)) < 0.2 * initialEnergyLevel
                % Force cluster head replacement due to low energy
                potentialHeads = members(nodeEnergy(members) >= 0.3 * initialEnergyLevel);
                if ~isempty(potentialHeads)
                    [~, bestIdx] = max(nodeTrustValues(potentialHeads) + nodeEnergy(potentialHeads));
                    clusterHeads(k) = potentialHeads(bestIdx);
                end
            elseif length(members) > 1
                % Regular optimization based on trust and energy
                memberScores = nodeTrustValues(members) + nodeEnergy(members);
                memberScores(isMalicious(members) == 1) = memberScores(isMalicious(members) == 1) * 0.5;
                [~, bestIdx] = max(memberScores);
                bestMember = members(bestIdx);
                
                currentScore = nodeTrustValues(clusterHeads(k)) + nodeEnergy(clusterHeads(k));
                if isMalicious(clusterHeads(k)) == 1
                    currentScore = currentScore * 0.5;
                end
                
                if bestMember ~= clusterHeads(k) && memberScores(bestIdx) > currentScore
                    clusterHeads(k) = bestMember;
                end
            end
        end
        
        % Reassign nodes to clusters
        clusters = cell(1, noOfClusters);
        for i = 1:noOfNodes
            if ~ismember(i, clusterHeads)
                dists = zeros(1, noOfClusters);
                for k = 1:noOfClusters
                    dists(k) = sqrt((nodeXLoc(i) - nodeXLoc(clusterHeads(k)))^2 + (nodeYLoc(i) - nodeYLoc(clusterHeads(k)))^2);
                end
                [~, idx] = min(dists);
                clusters{idx} = [clusters{idx}, i];
                clusterMembership(i) = idx;
            else
                idx = find(clusterHeads == i);
                clusters{idx} = [clusters{idx}, i];
                clusterMembership(i) = idx;
            end
        end
    end
end

% ========================== RESULTS ANALYSIS ==========================
% Calculate final metrics
maxTrustValue = max(trustValues);
maxDeliveryRatio = max(deliveryRatios);
successfulInteractions = sum(successfulInteraction_mat == 1);
overheadRatio = (numberOfInteractions - successfulInteractions) / numberOfInteractions;

% Output results
fprintf('The number of successful interactions is: %d\n', successfulInteractions);
fprintf('Maximum trust value among all the interactions is: %.6f\n', maxTrustValue);
fprintf('Maximum delivery ratio among all the interactions is: %.6f\n', maxDeliveryRatio);
fprintf('The Overhead Ratio of the Trust Agent is: %.2f\n', overheadRatio);

% ========================== VISUALIZATIONS ==========================
% Plot interaction metrics
figure;
subplot(3, 1, 1);
plot(1:numberOfInteractions, successfulInteraction_mat, 'LineWidth', 1.5);
xlabel('Number of Interactions');
ylabel('Successful Interactions');
title('Successful Interactions vs Number of Interactions');
grid on;

subplot(3, 1, 2);
plot(1:numberOfInteractions, deliveryRatios, 'LineWidth', 1.5);
xlabel('Number of Interactions');
ylabel('Delivery Ratio');
title('Delivery Ratio vs Number of Interactions');
grid on;

subplot(3, 1, 3);
hold on;
plot(1:numberOfInteractions, trustValues, 'b', 'LineWidth', 1);
plot(1:numberOfInteractions, movingAvgTrust, 'r', 'LineWidth', 2);
xlabel('Number of Interactions');
ylabel('Trust Value');
title('Trust Value vs Number of Interactions');
legend('Individual Interaction Trust', 'Moving Average (10)');
grid on;
hold off;

% Plot final clusters
figure;
hold on;
colors = lines(noOfClusters);
for k = 1:noOfClusters
    members = clusters{k};
    if ~isempty(members)
        nonMaliciousMembers = members(isMalicious(members) == 0);
        maliciousMembers = members(isMalicious(members) == 1);
        
        % Plot non-malicious nodes
        if ~isempty(nonMaliciousMembers)
            plot(nodeXLoc(nonMaliciousMembers), nodeYLoc(nonMaliciousMembers), '.', 'Color', colors(mod(k-1,size(colors,1))+1,:), 'MarkerSize', 12);
        end
        
        % Plot malicious nodes with x marker
        if ~isempty(maliciousMembers)
            plot(nodeXLoc(maliciousMembers), nodeYLoc(maliciousMembers), 'x', 'Color', 'r', 'MarkerSize', 8);
        end
        
        % Plot cluster heads
        plot(nodeXLoc(clusterHeads(k)), nodeYLoc(clusterHeads(k)), 'o', 'MarkerSize', 12, 'MarkerFaceColor', colors(mod(k-1,size(colors,1))+1,:), 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    end
end
title('Cluster Formation using Red Fox Optimization');
xlabel('X position (m)');
ylabel('Y position (m)');
legend('Regular Nodes', 'Malicious Nodes', 'Cluster Heads');
hold off;
