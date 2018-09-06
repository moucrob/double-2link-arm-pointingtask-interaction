function [i1, j1, i2, j2] = bijecNtoN2xN2(i, nbRows,nbColumns)
nbPos = nbRows*nbColumns;

k1 = ceil(i/(nbPos-1)); %k1 being the index into posVec(k1,:)
%and revealing which "part" of the long request list we are.
%(and the index of the removed row into posVec for establish k2)

%rows are divided into 24 subsets of 23 doublets of doublets.
indexes = 1:24 ; indexes(k1) = [];
k2tmp = i - (k1-1)*(nbPos-1);

if k2tmp < k1
    k2 = k2tmp;
else %k2tmp >= k1
    k2 = k2tmp + 1;
end

%Then, we know how posVec reshaped posMat:
%Beginning from the extreme top left of the matrix,
%and browsing from left to right.
i1 = ceil(k1/nbColumns); %the row of the realWorldPos1
j1 = k1 - (i1-1)*nbColumns;
i2 = ceil(k2/nbColumns);
j2 = k2 - (i2-1)*nbColumns;

end