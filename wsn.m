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

for i = 1:noOfNodes
    for j = 1:noOfNodes
        if i ~= j && neighborNode(i, j) == 1
            % Calculate Delivery Ratio (DR)
            deliveryRatio(i, j) = ((transmissionRange / interferenceRange)^2) * (deliveryRate / transmissionRate);
            
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
            cooperativeness(i, j) = successRate + reciprocityRate + avgResponseTime;
        end
    end
end

% Function to calculate Jaccards Similarity
function similarity = jaccardSimilarity(A, B)
    intersection = sum(A & B); % Number of common elements
    union = sum(A | B); % Number of unique elements
    similarity = intersection / union;
end


% Algorithm 1: Direct Observations-Based Trust Evaluation
function TrustEvaluation(a, b)
    % a is the trustor, b is the trustee

    % Step 2: Identify trustee (b)
    bid = b;

    % Step 3: check whether a and b have any previous trust values

        if compatibility(a, b) == 0 && cooperativeness(a, b) && deliveryRatio(a, b)
          TrustEvaluationAbsolute(a, b); % Proceed to absolute trust evaluation
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
        
        
       % step 6 : formulation of trust value
       form_trust_a_to_b = 0;


        % Step 7: Threshold comparison
        threshold = 5; % Trust threshold
        if form_trust_a_to_b >= threshold
            disp('Provide Services'); % Provide services if trust threshold is met
        else
            disp('Decline'); % Decline service request if trust threshold is not met
        end
    end


% Algorithm 2: Absolute Trust Formulation of Parameters
function TrustEvaluationAbsolute(a, b)
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
        
        % Step 4: Threshold comparison
        threshold = 5;
        if summation_Tform_a_b >= threshold
            disp('Provide Services'); % Provide services if trust threshold is met
        else
            disp('Decline'); % Decline service request if trust threshold is not met
        end
    else
        TrustEvaluationRecommendation(a, b); % Proceed to recommendation-based trust evaluation
    end
end

% Algorithm 3: Recommendation-Based Indirect Trust Evaluation
function TrustEvaluationRecommendation(a, b)
    % a is the trustor, b is the trustee
    
    % Step 2: Identify trustee (b)
    bid = b;
    
    % Step 3: Check recommendations from a to b (not implemented here)
    if isRecommended(a, b)
        % Step 4: Get recommended nodes from a to b (not implemented here)
        recommendedNodes = getRecommendations(a, b);
        
        % Step 5: Calculate aggregated trust from recommended nodes to b
        recom_form_b_kth = 0;
        
        for k = recommendedNodes
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
        
        % Step 7: TB certification (allocate certificate)
        TBcert_allocate(a, b); % Allocate trust certificate
        
        % Step 8: Threshold comparison
        threshold = 5; % Trust threshold
        if Tupdate_a_b >= threshold
            disp('Provide Services'); % Provide services if trust threshold is met
        else
            disp('Decline'); % Decline service request if trust threshold is not met
        end
    else
        disp('Decline'); % Decline if no valid recommendations
    end
end

% Function to check if b is recommended by a (placeholder)
function isRec = isRecommended(a, b)
    % Placeholder function to check if b is recommended by a
    % Implement logic to check recommendations (not implemented here)
    isRec = true; % Placeholder for demonstration
    % Implement actual logic to check recommendations
end

% Function to get recommended nodes from a to b (placeholder)
function recommendedNodes = getRecommendations(a, b)
    % Placeholder function to retrieve recommended nodes from a to b
    % Implement logic to get recommended nodes (not implemented here)
    recommendedNodes = [1, 2, 3]; % Placeholder for recommended nodes
    % Implement actual logic to get recommended nodes
end