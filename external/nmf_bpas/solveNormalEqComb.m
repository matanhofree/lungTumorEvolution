function [ Z,iter, pinvCounter ] = solveNormalEqComb( AtA,AtB,PassSet )
% Solve normal equations using combinatorial grouping.
% Although this function was originally adopted from the code of
% "M. H. Van Benthem and M. R. Keenan, J. Chemometrics 2004; 18: 441-450",
% important modifications were made to fix bugs.
%
% Modified by Jingu Kim (jingu@cc.gatech.edu)
%             School of Computational Science and Engineering,
%             Georgia Institute of Technology
%
% Last updated Aug-12-2009

    iter = 0;
    nearSingSolver = exist('lsqminnorm','file') == 2;
    pinvCounter = 0;
    
    if (nargin ==2) || isempty(PassSet) || all(PassSet(:))
        Z = AtA\AtB;
        iter = iter + 1;
    else
        Z = zeros(size(AtB));
        [n,k1] = size(PassSet);

        %% Fixed on Aug-12-2009
        if k1==1
            if rcond(full(AtA(PassSet,PassSet))) < 1e-5
                if pinvCounter < 1
                    fprintf('Warn: near singular using pinv (k1)\n');                
                end
                pinvCounter = pinvCounter + 1;
                Z(PassSet)=pinv(full(AtA(PassSet,PassSet)))*AtB(PassSet); 
            else
                Z(PassSet)=AtA(PassSet,PassSet)\AtB(PassSet); 
            end
        else
            %% Fixed on Aug-12-2009
            % The following bug was identified by investigating a bug report by Hanseung Lee.
            [sortedPassSet,sortIx] = sortrows(PassSet');
            breaks = any(diff(sortedPassSet)');
            breakIx = [0 find(breaks) k1];
            % codedPassSet = 2.^(n-1:-1:0)*PassSet;
            % [sortedPassSet,sortIx] = sort(codedPassSet);
            % breaks = diff(sortedPassSet);
            % breakIx = [0 find(breaks) k1];

            for k=1:length(breakIx)-1
                cols = sortIx(breakIx(k)+1:breakIx(k+1));
                vars = PassSet(:,sortIx(breakIx(k)+1));
                if rcond(full(AtA(vars,vars))) < 1e-5
                    if pinvCounter < 1
                        fprintf('Warn: near singular using pinv\n');                
                    end
                    pinvCounter = pinvCounter + 1;

                    if nearSingSolver
                        Z(vars,cols) = lsqminnorm(full(AtA(vars,vars)),AtB(vars,cols));
                    else
                        Z(vars,cols) = pinv(full(AtA(vars,vars)))*AtB(vars,cols);
                    end                    
                    % lsqminnorm -- might work better
                else
                    Z(vars,cols) = AtA(vars,vars)\AtB(vars,cols);
                end
                iter = iter + 1;
            end
        end
    end

    if pinvCounter > 1
        fprintf('Warn: near singular observed (%d)\n',pinvCounter);                
    end

end
