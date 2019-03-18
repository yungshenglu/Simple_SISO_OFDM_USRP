function v = read_complex_binary (filename, start, count)
%% usage: read_complex_binary (filename, [count])
%%
%% open filename and return the contents as a column vector, 
%% treating them as 32 bit complex numbers
%%
m = nargchk(1, 3, nargin);

if m
    usage (m);
end

if nargin < 2
    start = 0;
    count = Inf;
end

f = fopen(filename, 'rb');
if f >= 0
    fseek(f, start, -1);
    t = fread(f, [2, count], 'float');
    fclose(f);
    v = t(1, :) + t(2, :) * 1i;
    [r, c] = size(v);
end