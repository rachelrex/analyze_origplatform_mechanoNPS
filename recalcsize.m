function newsize=recalcsize(matrix,deff,L)

[numrows,~]=size(matrix);
newsize=zeros(numrows,1);
for i=1:numrows
    newsize(i,1)=calcsize(deff,L,matrix(i,1));
end
