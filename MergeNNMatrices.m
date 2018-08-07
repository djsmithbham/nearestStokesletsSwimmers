%MERGENNMATRICES Merges two nearest-neighbour matrices
function NN = MergeNNMatrices(NNa,NNb)

% Keep only top third section of NN matrices
[mA,nA] = size(NNa);
mA = mA/3; nA = nA/3;
NNa = NNa(1:mA,1:nA);

[mB,nB] = size(NNb);
mB = mB/3; nB = nB/3;
NNb = NNb(1:mB,1:nB);

% Get data from matrices
[iA,jA] = find(NNa);
[iB,jB] = find(NNb);

% Create new NN matrix
NN = sparse([iA(:);mA + iB(:)],[jA(:);nA + jB(:)], ...
    ones(length([iA(:);iB(:)]),1),mA+mB,nA+nB);

NN = kron(speye(3),NN);

end