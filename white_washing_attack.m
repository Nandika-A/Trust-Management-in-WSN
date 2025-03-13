clear all;
clc;

% Simulation Parameters
simArea = 1000; % Simulation area size (square area in meters)
noOfNodes = 600; % Number of nodes in the network
transmissionRange = 25; % Transmission range in meters
interferenceRange = 30; % Interference range in meters
transmissionRate = 4; % Transmission rate in Mbps
deliveryRate = 8; % Delivery rate (packets per 30 seconds)
initialEnergyLevel = 3; % Initial energy level
numberOfInteractions = 300; %Number of interactions
% Initialize variables
deliveryRatios = zeros(1, numberOfInteractions);
successfulInteraction_mat = zeros(1, numberOfInteractions);
trustValues = zeros(1, numberOfInteractions);
energy = 3; % Initial energy level

maliciousPercentage = 30; % Percentage of nodes that are malicious
noOfMaliciousNodes = ceil(maliciousPercentage / 100 * noOfNodes);
maliciousNodes = randperm(noOfNodes, noOfMaliciousNodes); % Randomly select malicious nodes

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
trustValue = zeros(noOfNodes, noOfNodes);
recommendationList = containers.Map();

% Specify number of interactions
trustValues = zeros(1, numberOfInteractions);
deliveryRatios = zeros(1, numberOfInteractions);
successfulInteraction_mat = zeros(1, numberOfInteractions);

% Function to calculate Jaccards Similarity
function similarity = jaccardSimilarity(A, B)
    intersection = sum(A & B); % Number of common elements
    union = sum(A | B); % Number of unique elements
    similarity = intersection / union;
end

function DR = deliveryRatioCalc(transmissionRange, interferenceRange, deliveryRate, transmissionRate)
    DR = ((transmissionRange / interferenceRange)^2) * (deliveryRate / transmissionRate);
end

function isNew = isNewNode(a, newNode)
    if newNode(num2str(a)) == 1
        isNew = 1;
    else
        isNew = 0;
    end
end

% Function to check if b is recommended by a
function isRec = isRecommended(a, b, recommendedNodes)
    % Function to check if b is recommended by a based on a recommendation list
    
    % Check if node b is in the recommendation list from a
    if isKey(recommendedNodes, num2str(a))
        isRec = ismember(num2str(b), recommendedNodes(num2str(a)));
    else
        isRec = false; % No recommendations from a
    end
end

for i = 1:noOfNodes
    for j = 1:noOfNodes
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

% Algorithm 1: Direct Observations-Based Trust Evaluation
function success = TrustEvaluation(a, b, trustValue, newNode, recommendedNodes, compatibility, cooperativeness, deliveryRatio)
    % a is the trustor, b is the trustee

    % Step 2: Identify trustee (b)
    bid = b;

    % Step 3: check whether a and b have any previous trust values

        if trustValue(a, b) == 0
          success = TrustEvaluationAbsolute(a, b, trustValue, newNode, recommendedNodes, compatibility, cooperativeness, deliveryRatio); % Proceed to absolute trust evaluation
        end 

    % if there are previous trust values, go within algorithm
        % Step 5: Calculate aggregated trust from a to b
        trust_a_to_b = 0;

       
        % Calculate Compatibility of a towards b
        compatibility_a_b = compatibility(a, b);

        % Calculate Cooperativeness of a towards b
        cooperativeness_a_b = cooperativeness(1, 1);

        % Calculate Delivery ratio of a towards b
        deliveryRatio_a_b = deliveryRatio(a, b);

        % Sum up individual trust parameters to get total trust
        trust_a_to_b = compatibility_a_b + cooperativeness_a_b + deliveryRatio_a_b;
        
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
function success = TrustEvaluationAbsolute(a, b, trustValue, newNode, recommendedNodes, compatibility, cooperativeness, deliveryRatio)
    % a is the trustor, b is the trustee
    
    % Step 2: Identify trustee (b)
    bid = b;
    
    % Step 3: Check if b is a new node
    if isNewNode(b, newNode)
        % Step 3: Evaluate Trust
        compatibility_b = compatibility(a, b); % Compatibility of b towards a
        cooperativeness_b = cooperativeness(1, 1); % Cooperativeness of b towards a
        deliveryRatio_b = deliveryRatio(a, b); % Delivery ratio of b towards a
        
        summation_Tform_a_b = compatibility_b + cooperativeness_b + deliveryRatio_b;
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
    else
        success = TrustEvaluationRecommendation(a, b, trustValue, recommendedNodes, compatibility, cooperativeness, deliveryRatio); % Proceed to recommendation-based trust evaluation
    end
end

% Algorithm 3: Recommendation-Based Indirect Trust Evaluation
function success = TrustEvaluationRecommendation(a, b, trustValue, recommendedNodes, compatibility, cooperativeness, deliveryRatio)
    % a is the trustor, b is the trustee
    
    % Step 2: Identify trustee (b)
    bid = b;
    
    % Step 3: Check recommendations from a to b (not implemented here)
    if isRecommended(a, b, recommendedNodes)
        % Step 4: Get recommended nodes from a to b (not implemented here)
        recommended = recommendedNodes(b, recommendedNodes);
        
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
        threshold = 5; % Trust threshold
        if Tupdate_a_b >= threshold
            success = 1;
            disp('Provide Services'); % Provide services if trust threshold is met
        else
            success = 0;
            disp('Decline'); % Decline service request if trust threshold is not met
        end
    else
        success = 0;
        disp('Decline'); % Decline if no valid recommendations
    end
end

whitewashPeriod = 3; % Number of interactions after which malicious nodes reset their identity
whitewashCounter = zeros(1, noOfNodes); % Counter to track interactions for each node

for interaction = 1:numberOfInteractions
    nodeA = randi([1, noOfNodes]);
    nodeB = randi([1, noOfNodes]);

    % Whitewashing: Periodically reset malicious nodes
    if ismember(nodeA, maliciousNodes) % If nodeA is malicious
        whitewashCounter(nodeA) = whitewashCounter(nodeA) + 1;
        % Reset malicious node after every whitewashPeriod interactions
        if whitewashCounter(nodeA) > whitewashPeriod
            % Reset node identity by assigning a new ID
            newNode(num2str(nodeA)) = 1; % Make node A as new
            trustValue(nodeA, :) = 0; % Reset trust values involving nodeA
            trustValue(:, nodeA) = 0; % Reset trust values involving nodeA
            disp(['Malicious node ' num2str(nodeA) ' reset its identity']);
            whitewashCounter(nodeA) = 0; % Reset the counter for the next period
        end
    end

    % Whitewashing: Periodically reset malicious nodes
    if ismember(nodeB, maliciousNodes) % If nodeB is malicious
        whitewashCounter(nodeB) = whitewashCounter(nodeB) + 1;
        % Reset malicious node after every whitewashPeriod interactions
        if whitewashCounter(nodeB) > whitewashPeriod
            % Reset node identity by assigning a new ID
            newNode(num2str(nodeB)) = 1; % Make node B as new
            trustValue(nodeB, :) = 0; % Reset trust values involving nodeA
            trustValue(:, nodeB) = 0; % Reset trust values involving nodeA
            disp(['Malicious node ' num2str(nodeB) ' reset its identity']);
            whitewashCounter(nodeB) = 0; % Reset the counter for the next period
        end
    end

    successfulInteraction_mat(interaction) = TrustEvaluation(nodeA, nodeB, trustValue, newNode, recommendedNodes, compatibility, cooperativeness, deliveryRatio);
    drate = deliveryRate + randi([-1, 1]);
    trate = transmissionRate + randi([-1, 1]);
    A = [transmissionRange, interferenceRange, trate, drate, energy];
    B = [transmissionRange, interferenceRange, transmissionRate, deliveryRate, initialEnergyLevel];
    intersection = sum(A & B); % Number of common elements
    union = sum(A | B); % Number of unique element
    compatibility_val =  intersection / union;
    rt = randi([0, 5]);
    if(successfulInteraction_mat(interaction) == 1)
         energy = energy + 1;
         cooperativeness = 1 + rt;
    else
        energy = energy - 1;
        cooperativeness = rt;
    end
    DR = ((transmissionRange / interferenceRange)^2) * (drate / trate);
    deliveryRatios(interaction) = DR + energy;
    deliveryRatio(nodeA, nodeB) = deliveryRatios(interaction);
    trustValues(interaction) = DR + compatibility_val + energy + cooperativeness;
    trustValue(nodeA, nodeB) = trustValues(interaction);
end
maxTrustValue = max(trustValues)
maxDeliveryRatio = max(deliveryRatios)
successfulInteractions = sum(successfulInteraction_mat == 1);
% Calculate Overhead Ratio
overheadRatio = (numberOfInteractions - successfulInteractions) / numberOfInteractions;
fprintf('The number of successful interactions is: %d\n', successfulInteractions);
fprintf('Maximum trust value among all the interactions is: %d\n', maxTrustValue);
fprintf('Maximum delivery ratio among all the interactions is: %d\n', maxDeliveryRatio);
% Display the Overhead Ratio
fprintf('The Overhead Ratio of the Trust Agent is: %.2f\n', overheadRatio);

% Plot successfulInteractions against the number of interactions
subplot(3, 1, 1);
plot(1:numberOfInteractions, successfulInteraction_mat);
xlabel('Number of Interactions');
ylabel('Successful Interactions');
title('Successful Interactions vs Number of Interactions');

% Plot deliveryRatio against number of interactions
subplot(3, 1, 2);
plot(1:numberOfInteractions, deliveryRatios);
xlabel('Number of Interactions');
ylabel('Delivery Ratio');
title('Delivery Ratio vs Number of Interactions');

% Plot trustValue against number of interactions
subplot(3, 1, 3);
plot(1:numberOfInteractions, trustValues);
xlabel('Number of Interactions');
ylabel('Trust Value');
title('Trust Value vs Number of Interactions');
% Sample data
x = 1:10;           % Common x-axis values
y1 = trustValues(1:10);          % First dataset
y2 = [5, 4, 3, 2.75, 2.6, 2.5, 2.4, 2, 4.8, 5];          % Second dataset
y3 = [5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9];       % Third dataset
y4 = [3, 3.5, 4, 2, 4.9, 5.8, 7, 8, 8.5, 9];
y5 = [6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.35, 6.25, 6.1, 6];

% Plot all arrays in one figure with different colors
figure; % Opens a new figure window
hold on; % Retain plots so new plots are added to the current figure
plot(x, y1, 'r-o','LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Proposed Results');   % Plot y1 in red
plot(x, y2, 'b-o','LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'LightTrust');   % Plot y2 in blue
plot(x, y3, 'g-o','LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'SMTE'); % Plot y3 in green
% Plot the fourth graph
plot(x, y4, 'm-o','LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'ETES'); % Plot y4 in magenta
% Plot the fifth graph
plot(x, y5, 'c-o','LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'GFMS'); % Plot y5 in cyan
hold off; % Release the hold to avoid affecting future plots

% Add labels, legend, and title
xlabel('X-axis');
ylabel('Y-axis');
title('Trust values comparison against white washing attack');
legend('Location', 'northeastoutside'); % Automatically display labels from 'DisplayName'
grid on; 
