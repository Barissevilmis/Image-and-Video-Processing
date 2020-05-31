function J = weber(L1,L2,Lb)

% Filled in according to the given coordinates on lab assignment sheet
% Matrix by 600 x 600
    J = ones(600,600);
    J = Lb .* J;
    J(221:380,221:300) = L1;
    J(221:380,301:380) = L2;
    
end