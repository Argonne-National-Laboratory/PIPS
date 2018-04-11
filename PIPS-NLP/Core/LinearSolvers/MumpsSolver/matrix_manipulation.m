A = load('matrix.txt');
n = A(1,1); nnz=A(1,2);
A = A(2:end,:);

n
nnz


if(nnz~=size(A,1)); error('nnz does not match'); end
im = max(A(:,1));
if(im<0 || im>=n); error('size does not match i indexes'); end
jm = max(A(:,2));
if(jm<0 || jm>=n); error('size does not match j indexes'); end

i=A(:,1)+1; j=A(:,2)+1; els=A(:,3);
M= sparse(i,j,els);
M=0.5*(M+M');
spy(M)
M\ones(n,1)