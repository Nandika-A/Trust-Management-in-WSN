clear all;
clc;

% Simulation Parameters
simArea = 1000; % Simulation area size (square area in meters)
noOfNodes = 600; % Number of nodes
transmissionRange = 25; % Transmission range in meters
interferenceRange = 30; % Interference range in meters
transmissionRate = 4; % Transmission rate in Mbps
deliveryRate = 8; % Delivery rate (packets per 30 seconds)
initialEnergyLevel = 3; % Initial energy level

% Malicious Node Percentage
maliciousPercentage = 30; % Percentage of nodes that are malicious
noOfMaliciousNodes = ceil(maliciousPercentage / 100 * noOfNodes);

% Node Deployment
nodeXLoc = rand(1, noOfNodes) * simArea; % Random x-coordinates within simulation area
nodeYLoc = rand(1, noOfNodes) * simArea; % Random y-coordinates within simulation area
global newNode;
% Define a mapping called newNode
newNode = containers.Map();
for i = 1:noOfNodes
    % Initialize the values randomly as zeros and ones
    newNode(num2str(i)) = randi([0, 1]);
end

global recommendedNodes;
% Define a mapping called recommendedNodes
recommendedNodes = containers.Map();
recommendedNodes(num2str(i)) = randi([1, noOfNodes], 1, noOfNodes);

% Neighbor Node Determination
neighborNode = zeros(noOfNodes, noOfNodes);

for i = 1:noOfNodes
    for j = 1:noOfNodes
        if i ~= j
            % Calculate distance between nodes
            distance = sqrt((nodeXLoc(i) - nodeXLoc(j))^2 + (nodeYLoc(i) - nodeYLoc(j))^2);
            
            % Check if nodes are within transmission range
            if distance <= transmissionRange
                neighborNode(i, j) = 1;
            end
        end
    end
end

% Initialize Interaction Tracking
totalInteractions = 0;
successfulInteractions = 0;
reciprocatedInteractions = 0;
totalResponseTime = zeros(1, noOfNodes);

% Trust Calculation and Matrix Initialization
deliveryRatio = zeros(noOfNodes, noOfNodes);
compatibility = zeros(noOfNodes, noOfNodes);
cooperativeness = zeros(noOfNodes, noOfNodes);
recommendationList = containers.Map();

% Specify number of interactions
numberOfInteractions = 100;
trustValues = zeros(1, numberOfInteractions);
deliveryRatios = zeros(1, numberOfInteractions);
successfulInteraction_mat = zeros(1, numberOfInteractions);
trustValue = zeros(noOfNodes, noOfNodes);

% Function to calculate Jaccards Similarity
function similarity = jaccardSimilarity(A, B)
    intersection = sum(A & B); % Number of common elements
    union = sum(A | B); % Number of unique elements
    similarity = intersection / union;
end

function DR = deliveryRatioCalc(transmissionRange, interferenceRange, deliveryRate, transmissionRate)
    DR = ((transmissionRange / interferenceRange)^2) * (deliveryRate / transmissionRate);
end

function isNew = isNewNode(a)
    if newNode(a) == 1
        isNew = 1;
    else
        isNew = 0;
    end
end

% Function to check if b is recommended by a
function isRec = isRecommended(a, b)
    % Function to check if b is recommended by a based on a recommendation list
    
    % Check if node b is in the recommendation list from a
    if isKey(recommendedNodes, a)
        isRec = ismember(b, recommendedNodes(a));
    else
        isRec = false; % No recommendations from a
    end
end

for i = 1:noOfNodes
    for j = 1:noOfNodes
        if i ~= j && neighborNode(i, j) == 1
            % Calculate Delivery Ratio (DR)
            deliveryRatio(i, j) = deliveryRatioCalc(transmissionRange, interferenceRange, deliveryRate, transmissionRate);
            
            % Calculate Compatibility using Jaccards Similarity
            nodeA = [transmissionRange, interferenceRange, transmissionRate, deliveryRate, initialEnergyLevel];
            nodeB = [transmissionRange, interferenceRange, transmissionRate, deliveryRate, initialEnergyLevel];
            compatibility(i, j) = jaccardSimilarity(nodeA, nodeB);
            
            % Simulate interaction between node i and node j
            totalInteractions = totalInteractions + 1;
            
            % Assume interaction results in successful transmission
            successfulInteractions = successfulInteractions + 1;
            
            % Track response time (random value for demonstration)
            responseTime = randi([0, 5]); % Placeholder for actual response time
            totalResponseTime(i) = totalResponseTime(i) + responseTime;
            
            % Check reciprocity (assume all interactions are reciprocated for demonstration)
            reciprocatedInteractions = reciprocatedInteractions + 1;
            
            % Calculate Success Rate
            successRate = successfulInteractions / totalInteractions;

            % Calculate Reciprocity Rate
            reciprocityRate = reciprocatedInteractions / totalInteractions;

            % Calculate Average Response Time per Node
            averageResponseTime = sum(totalResponseTime) / noOfNodes;

            % Calculate Cooperativeness
            cooperativeness(i, j) = successRate + reciprocityRate + averageResponseTime;
        end
    end
end

% Algorithm 1: Direct Observations-Based Trust Evaluation
function success = TrustEvaluation(a, b, trustValue)
    % a is the trustor, b is the trustee

    % Step 2: Identify trustee (b)
    bid = b;

    % Step 3: check whether a and b have any previous trust values

        if trustValue(a, b) == 0
          success = TrustEvaluationAbsolute(a, b, trustValue); % Proceed to absolute trust evaluation
        end 

    % if there are previous trust values, go within algorithm
        % Step 5: Calculate aggregated trust from a to b
        trust_a_to_b = 0;

       
        % Calculate Compatibility of a towards b
        compatibility_a_b = compatibility(a, b);

        % Calculate Cooperativeness of a towards b
        cooperativeness_a_b = cooperativeness(a, b);

        % Calculate Delivery ratio of a towards b
        deliveryRatio_a_b = deliveryRatio(a, b);

        % Sum up individual trust parameters to get total trust
        trust_a_to_b = compatibility_a_b + cooperativeness_a_b + deliveryRatio_a_b;
        
        % Step 7: Threshold comparison
        threshold = 10; % Trust threshold
        if trust_a_to_b >= threshold
            success = 1;
            disp('Provide Services'); % Provide services if trust threshold is met
        else
            success = 0;
            disp('Decline'); % Decline service request if trust threshold is not met
        end
    end


% Algorithm 2: Absolute Trust Formulation of Parameters
function success = TrustEvaluationAbsolute(a, b, trustValue)
    % a is the trustor, b is the trustee
    
    % Step 2: Identify trustee (b)
    bid = b;
    
    % Step 3: Check if b is a new node
    if isNewNode(b)
        % Step 3: Evaluate Trust
        compatibility_b = compatibility(a, b); % Compatibility of b towards a
        cooperativeness_b = cooperativeness(a, b); % Cooperativeness of b towards a
        deliveryRatio_b = deliveryRatio(a, b); % Delivery ratio of b towards a
        
        summation_Tform_a_b = compatibility_b + cooperativeness_b + deliveryRatio_b;
        trustValue(a, b) = summation_Tform_a_b; % Store trust value
        
        % Step 4: Threshold comparison
        threshold = 10;
        if summation_Tform_a_b >= threshold
            success = 1;
            disp('Provide Services'); % Provide services if trust threshold is met
        else
            success = 0;
            disp('Decline'); % Decline service request if trust threshold is not met
        end
    else
        success = TrustEvaluationRecommendation(a, b, trustValue); % Proceed to recommendation-based trust evaluation
    end
end

% Algorithm 3: Recommendation-Based Indirect Trust Evaluation
function success = TrustEvaluationRecommendation(a, b, trustValue)
    % a is the trustor, b is the trustee
    
    % Step 2: Identify trustee (b)
    bid = b;
    
    % Step 3: Check recommendations from a to b (not implemented here)
    if isRecommended(a, b)
        % Step 4: Get recommended nodes from a to b (not implemented here)
        recommended = recommendedNodes(b);
        
        % Step 5: Calculate aggregated trust from recommended nodes to b
        recom_form_b_kth = 0;
        
        for k = recommended
            % Calculate Compatibility of k towards b
            compatibility_k_b = compatibility(k, b);
            
            % Calculate Cooperativeness of k towards b
            cooperativeness_k_b = cooperativeness(k, b);
            
            % Calculate Delivery ratio of k towards b
            deliveryRatio_k_b = deliveryRatio(k, b);
            
            % Sum up individual trust parameters for each recommended node
            recom_form_b_kth = recom_form_b_kth + compatibility_k_b + cooperativeness_k_b + deliveryRatio_k_b;
        end
        
        % Step 6: Data Update
        Tupdate_a_b = recom_form_b_kth;
        trustValue(a, b) = Tupdate_a_b;
        
        % Step 7: TB certification (allocate certificate)
        TBcert_allocate(a, b); % Allocate trust certificate
        
        % Step 8: Threshold comparison
        threshold = 10; % Trust threshold
        if Tupdate_a_b >= threshold
            success = 1;
            disp('Provide Services'); % Provide services if trust threshold is met
        else
            success = 0;
            disp('Decline'); % Decline service request if trust threshold is not met
        end
    else
        disp('Decline'); % Decline if no valid recommendations
    end
end

energy = initialEnergyLevel;
% Loop for simulating interactions
for interaction = 1:numberOfInteractions
    % Select random pair of nodes for interaction
    nodeA = randi([1, noOfNodes]);
    nodeB = randi([1, noOfNodes]);
    
    % Check if nodes are neighbors and not the same node
    if neighborNode(nodeA, nodeB) == 1 && nodeA ~= nodeB

        % Simulate the algorithm
        success = TrustEvaluation(nodeA, nodeB, trustValue);

        % Trust updation based on algorithm results
        drate = deliveryRate + randi([-1, 1]);
        trate = transmissionRate + randi([-1, 1]);
        A = [transmissionRange, interferenceRange, trate, drate, energy];
        B = [transmissionRange, interferenceRange, transmissionRate, deliveryRate, initialEnergyLevel];
        compatibility = jaccardSimilarity(A, B);
        DR = ((transmissionRange / interferenceRange)^2) * (drate / trate);
        rt = randi([0, 5]);

        if(success == 1)
             energy = energy + 1;
             cooperativeness = 1 + rt;
        else
            energy = energy - 1;
            cooperativeness = rt;
        end

        % Populate the matrices to plot
        deliveryRatios(interaction) = DR + energy;
        trustValues(interaction) = DR + compatibility + energy + cooperativeness;
        successfulInteraction_mat(interaction) = success;
    end
end

% Plot trustValue against number of interactions
subplot(3, 1, 3);
plot(1:numberOfInteractions, trustValues);
xlabel('Number of Interactions');
ylabel('Trust Value');
title('Trust Value vs Number of Interactions');

% Plot deliveryRatio against number of interactions
subplot(3, 1, 2);
plot(1:numberOfInteractions, deliveryRatios);
xlabel('Number of Interactions');
ylabel('Delivery Ratio');
title('Delivery Ratio vs Number of Interactions');

% Plot successfulInteractions against the number of interactions
subplot(3, 1, 1);
plot(1:numberOfInteractions, successfulInteraction_mat);
xlabel('Number of Interactions');
ylabel('Successful Interactions');
title('Successful Interactions vs Number of Interactions');