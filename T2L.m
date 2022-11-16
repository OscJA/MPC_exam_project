function T2L(Col, Row, Mat, Path)
n = size(Mat, 1);
m = size(Mat, 2);

cs = [];
for i=1:m-1
    cs = [cs, ' c'];
end

fileID = fopen(Path, 'w');
fprintf(fileID, ['\\begin{tabular}{|c|', cs, ' |}\n']);
fprintf(fileID, "\\hline\n");

tex = ' ';
for i=1:m-1
    tex = [tex, '& $', Col(i), '$'];
end
tex = strjoin([tex, '\\\\\n']);

fprintf(fileID, tex);
fprintf(fileID, "\\hline\n");

for j=1:m
    tex = ['$', Row(j), '$'];
    for i=1:n
        tex = [tex, '& $', Mat(i, j), '$'];
    end
    tex = strjoin([tex, '\\\\\n']);
    fprintf(fileID, tex);
    fprintf(fileID, '\\hline\n');
end
fprintf(fileID, '\\end{tabular}\n');

