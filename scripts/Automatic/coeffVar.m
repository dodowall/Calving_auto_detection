function cv = coeffVar(v)
vstd = std(v);
vmean = mean(v);
cv = vstd./vmean;
end