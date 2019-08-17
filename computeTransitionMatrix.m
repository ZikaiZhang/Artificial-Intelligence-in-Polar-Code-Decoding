function Transition_matrix = computeTransitionMatrix(numberOfIterationForOneTimeWhileLoop, EbN0, N)
%computeTransitionMatrix
    % Parse the Markow_chain_data.txt file and get the value for type
    [type, number, value] = textread('Markov_chain_data.txt', '%s%d%d', 'delimiter', ',');
    
    % get the number of differennt EbN0 and store it in variable
    % differentEbN0
    differentEbN0 = length(EbN0);
    
    % create a 3D all zero matrix for storing the number of changes from 2 to 1
    ratioMatrixFor1 = zeros(differentEbN0, numberOfIterationForOneTimeWhileLoop, N - 1);
    
    % create a 3D all zero matrix for storing the number of changes from 1 to 2
    ratioMatrixFor2 = zeros(differentEbN0, numberOfIterationForOneTimeWhileLoop, N - 1);
    
    % create a 2D all zero matrix for storing the transition probability from 2
    % to 1
    averRatioMatrixFor1Temp = zeros(differentEbN0, numberOfIterationForOneTimeWhileLoop);
    
    % create a 2D all zero matrix for storing the transition probability from 1
    % to 2
    averRatioMatrixFor2Temp = zeros(differentEbN0, numberOfIterationForOneTimeWhileLoop);
    
    % create an 1D all zero matrix for storing the transition probability from 2
    % to 1
    averRatioMatrixFor1 = zeros(differentEbN0);
    
    % create an 1D all zero matrix for storing the transition probability from 1
    % to 2
    averRatioMatrixFor2 = zeros(differentEbN0);
    
    % create an 3D all zero matrix for storing the transition probability
    Transition_matrix = zeros(differentEbN0, 2, 2);
    
    for index = 1 : differentEbN0
        
        for i = 1 : numberOfIterationForOneTimeWhileLoop
            
            % set an identifier indicating some form of appearance of 'D' or 'G'
            DorGAppeared = false;
            
            for j = 1 : N - 1
                
                % now is the value of the current type
                now = char(type((index - 1) * numberOfIterationForOneTimeWhileLoop * (N - 1) + (i - 1) * (N - 1) + j));

                % perform some tasks if the list size is changed from 1 to
                % 2
                % if the type value 'D' or 'G' is encountered
                if now == 'D' || now == 'G'
                    
                    % set the identifier to true
                    DorGAppeared = true;
                    
                end
                
                % if 'E' is being visited and there is a 'D' or 'G'
                % appeared (some special form) 
                if now == 'E' && DorGAppeared == true
                    
                    % set the corresponding element in the 3D matrix to be
                    % 1
                    ratioMatrixFor1(index, i, j) = 1;
                    
                    % set the identifier to be false
                    DorGAppeared = false;
                    
                end

            end
            
            % set an identifier indicating some form of appearance of 'E'
            EAppeared = false;
            
            % set an identifier indicating the appearance of 'D' or 'G'
            DorGAppeared = false;
            
            for k = 1 : N - 1
                
                % now is the value of the current type
                now = char(type((index - 1) * numberOfIterationForOneTimeWhileLoop * (N - 1) + (i - 1) * (N - 1) + j));
                
                % perform some tasks if the list size is changed from 2 to
                % 1
                % if the type value 'E' is encountered
                if now == 'E'
                    
                    % set the identifier to be true
                    EAppeared = true;
                    
                end
                
                % if 'D' or 'G' is being visited and there is an 'E'
                % appeared(some special form)
                if (now == 'D' || now == 'G') && (EAppeared == true)
                    
                    % set the corresponding element in the 3D ratio matrix
                    % to be 1
                    ratioMatrixFor2(index, i, j) = 1;
                    
                    % set the identifier to be false
                    EAppeared = false;
                    
                    % set the identifier to be true
                    DorGAppeared = true;
                    
                end
                
                % if 'D' or 'G' is encountered for the first time in the
                % innermost iteration without a former appearance of E
                if (now == 'D' || now == 'G') && (EAppeared == false) && (DorGAppeared == false)
                    
                    % set the corresponding element in the 3D matrix to be
                    % 1
                    ratioMatrixFor2(index, i, j) = 1;
                    
                end
                
            end
            
            % compute the temparory ratio matrix
            averRatioMatrixFor1Temp(index, i) = sum(ratioMatrixFor1(index, i, 1:N - 1)) / (N - 1);
            
            % compute the temparory ratio matrix
            averRatioMatrixFor2Temp(index, i) = sum(ratioMatrixFor2(index, i, 1:N - 1)) / (N - 1);
            
        end
        
        % compute the transition probability matrix
        averRatioMatrixFor1(index) = sum(averRatioMatrixFor1Temp(index, 1:numberOfIterationForOneTimeWhileLoop)) / numberOfIterationForOneTimeWhileLoop;
        
        % compute the transition probability matrix
        averRatioMatrixFor2(index) = sum(averRatioMatrixFor2Temp(index, 1:numberOfIterationForOneTimeWhileLoop)) / numberOfIterationForOneTimeWhileLoop;
        
        % compute the final transition probability matrix
        Transition_matrix(index, 1, 1) = 1 - averRatioMatrixFor2(index);
        
        Transition_matrix(index, 1, 2) = averRatioMatrixFor2(index);
        
        Transition_matrix(index, 2, 1) = averRatioMatrixFor1(index);
        
        Transition_matrix(index, 2, 2) = 1 - averRatioMatrixFor1(index);    
    
    end
    
end

