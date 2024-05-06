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

% Define a mapping called newNode
newNode = containers.Map();
for i = 1:noOfNodes
    % Initialize the values randomly as zeros and ones
    newNode(num2str(i)) = randi([0, 1]);
end

% Define a mapping called recommendedNodes
recommendedNodes = containers.Map();
for i = 1:noOfNodes
    recommendedNodes(num2str(i)) = randi([1, noOfNodes], 1, noOfNodes);
end

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
                
                % Calculate Delivery Ratio (DR)
                deliveryRatio(i, j) = deliveryRatioCalc(transmissionRange, interferenceRange, deliveryRate, transmissionRate, i, j);
            end
        end
    end
end

% Initialize Interaction Tracking
totalInteractions = 0;
successfulInteractions = 0;
reciprocatedInteractions = 0;
success = 0;
totalResponseTime = zeros(1, noOfNodes);

% Trust Calculation and Matrix Initialization
compatibility = zeros(noOfNodes, noOfNodes);
cooperativeness = zeros(noOfNodes, noOfNodes);
trustValue = zeros(noOfNodes, noOfNodes);

for i = 1:noOfNodes
    for j = 1:noOfNodes
        if i ~= j && neighborNode(i, j) == 1
            % Calculate Compatibility using Jaccards Similarity
            nodeA = [transmissionRange, interferenceRange, transmissionRate, deliveryRate, initialEnergyLevel];
            nodeB = [transmissionRange, interferenceRange, transmissionRate, deliveryRate, initialEnergyLevel];
            compatibility(i, j) = jaccardSimilarity(nodeA, nodeB);
            
            % Simulate interaction between node i and node j
            totalInteractions = totalInteractions + 1;
            
            % Assume interaction results in successful transmission
            successfulInteractions = successfulInteractions + 1;
            
            % Track response time (random value for demonstration)
            responseTime = rand(); % Placeholder for actual response time
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

% Function to calculate Jaccards Similarity
function similarity = jaccardSimilarity(A, B)
    intersection = sum(A & B); % Number of common elements
    union = sum(A | B); % Number of unique elements
    similarity = intersection / union;
end

% Function to calculate Delivery Ratio
function DR = deliveryRatioCalc(transmissionRange, interferenceRange, deliveryRate, transmissionRate, i, j)
    DR = ((transmissionRange / interferenceRange)^2) * (deliveryRate / transmissionRate);
end

% Algorithm 1: Direct Observations-Based Trust Evaluation
function TrustEvaluation(a, b)
    % a is the trustor, b is the trustee

    % Step 2: Identify trustee (b)
    bid = b;

    % Step 3: check whether a and b have any previous trust values
    if trustValue(a, b) == 0
        TrustEvaluationAbsolute(a, b); % Proceed to absolute trust evaluation
    end 

    % if there are previous trust values, go within algorithm
    % Step 5: Calculate aggregated trust from a to b
    trust_a_to_b = 0;


    % Calculate Compatibility of a towards b
    compatibility_a_b = compatibility(a, b);

    % Calculate Cooperativeness of a towards b
    cooperativeness_a_b = cooperativeness(a, b);

    % Sum up individual trust parameters to get total trust
    trust_a_to_b = compatibility_a_b + cooperativeness_a_b;

    % Step 7: Threshold comparison
    threshold = 5; % Trust threshold
    if trust_a_to_b >= threshold
        success = 1;
        disp('Provide Services'); % Provide services if trust threshold is met
    else
        success = 0;
        disp('Decline'); % Decline service request if trust threshold is not met
    end
end


% Algorithm 2: Absolute Trust Formulation of Parameters
function TrustEvaluationAbsolute(a, b)
    % a is the trustor, b is the trustee

    % Step 2: Identify trustee (b)
    bid = b;

    % Step 3: Evaluate Trust
    compatibility_b = compatibility(a, b); % Compatibility of b towards a
    cooperativeness_b = cooperativeness(a, b); % Cooperativeness of b towards a

    summation_Tform_a_b = compatibility_b + cooperativeness_b;
    trustValue(a, b) = summation_Tform_a_b; % Store trust value

    % Step 4: Threshold comparison
    threshold = 5;
    if summation_Tform_a_b >= threshold
        success = 1;
        disp('Provide Services'); % Provide services if trust threshold is met
    else
        success = 0;
        disp('Decline'); % Decline service request if trust threshold is not met
    end
end

% Specify number of interactions
numberOfInteractions = 100;
trustValues = zeros(1, numberOfInteractions);
deliveryRatios = zeros(1, numberOfInteractions);
successfulInteractions = zeros(1, numberOfInteractions);

% Loop for simulating interactions
for interaction = 1:numberOfInteractions
    % Select random pair of nodes for interaction
    nodeA = randi([1, noOfNodes]);
    nodeB = randi([1, noOfNodes]);
    
    % Check if nodes are neighbors and not the same node
    if neighborNode(nodeA, nodeB) == 1 && nodeA ~= nodeB
        success = 0;
        % Simulate interaction between nodeA and nodeB
        TrustEvaluation(nodeA, nodeB);
        deliveryRate_new = deliveryRate + randi([-2, 2]); % Random delivery rate change
        deliveryRatio(nodeA, nodeB) = deliveryRatioCalc(transmissionRange, interferenceRange, deliveryRate_new, transmissionRate)
                
        % Store values in array
        trustValues(interaction) = trustValue(nodeA, nodeB);
        deliveryRatios(interaction) = deliveryRatio(nodeA, nodeB);
        successfulInteractions(interaction) = success;
    end
end

% Plot trustValue against number of interactions
// subplot(3, 1, 1);
plot(1:numberOfInteractions, trustValues);
xlabel('Number of Interactions');
ylabel('Trust Value');
title('Trust Value vs Number of Interactions');

% Plot deliveryRatio against number of interactions
// subplot(3, 1, 2);
plot(1:numberOfInteractions, deliveryRatios);
xlabel('Number of Interactions');
ylabel('Delivery Ratio');
title('Delivery Ratio vs Number of Interactions');

% Plot successfulInteractions against the number of interactions
// subplot(3, 1, 3);
plot(1:numberOfInteractions, successfulInteractions);
xlabel('Number of Interactions');
ylabel('Successful Interactions');
title('Successful Interactions vs Number of Interactions');