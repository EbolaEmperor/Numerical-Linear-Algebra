n = 1024;
fileID = fopen('output.txt', 'w');

%%user_code_here

h = 1 / (n + 1);
b = (h * h) * ones(n * n, 1);
x = solveLaplacian(n, b);
A = genDiff2(n);

if norm(A * x - b) > 1e-6
    fprintf(fileID, 'Too large error!\n');
    fclose(fileID);
    return
end

fprintf(fileID, 'The North remembers\n');
fclose(fileID);