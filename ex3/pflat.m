%% Dana Joffe 312129240
function output = pflat(matrix)
    output = matrix ./ matrix(end, :); 
end