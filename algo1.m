
% Algorithm 1: Direct Observations-Based Trust Evaluation
function TrustEvaluation(a, b)
    % a is the trustor, b is the trustee

    % Step 2: Identify trustee (b)
    bid = b;

    % Step 3: request from b to a (not implemented here)
    % --------------------------------------------------------------------------------------------------------------------------------

    % Step 4: check whether a and b have any previous trust values

        if compatibility(a, b) == 0 && cooperativeness(a, b) && deliveryRatio(a, b)
            %this algorithm will not work so go to algorithm 2
            %---------------------------------------------------------------------------------------
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
        
        
       % step 5 : formulation of trust value
       form_trust_a_to_b = 0;


        % Step 6: Threshold comparison
        threshold = 5; % Trust threshold
        if form_trust_a_to_b >= threshold
            disp('Provide Services'); % Provide services if trust threshold is met
        else
            disp('Decline'); % Decline service request if trust threshold is not met
        end
    else
        disp('Decline'); % Decline if no valid recommendations
    end
end
