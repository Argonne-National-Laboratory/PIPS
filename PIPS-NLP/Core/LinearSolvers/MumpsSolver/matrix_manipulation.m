clear all

A = load('matrix.txt');
n = A(1,1); nonz=A(1,2);
A = A(2:end,:);

[n nonz]

close all

if(nonz~=size(A,1)); error('nnz does not match'); end
im = max(A(:,1));
if(im<0 || im>=n); error('size does not match i indexes'); end
jm = max(A(:,2));
if(jm<0 || jm>=n); error('size does not match j indexes'); end

i=A(:,1)+1; j=A(:,2)+1; els=A(:,3);
M0= sparse(i,j,els);
M=(tril(M0,-1)+M0');
figure;spy(M)
tic;
%x = M\ones(n,1);
[L,D,P,S] = ldl(M, 0.01, 'lower', 'vector');
figure; spy(L)
nnz(L)
toc

tic;
[L,D,P] = ldl(M, 'vector');
figure; spy(M(P,P))
toc

tic;
[L,D,P] = ldl(M);
figure; spy(L)
toc

tic;
x=M\ones(n,1);
toc;
fprintf('solve error %g\n', norm(M*x-ones(n,1)));
fprintf(' vec %g %g\n', x(1), x(2));

figure;spy(P'*M*P)


P = load('ordering.txt');
if(length(P)~=n) error('size of the ordering does not match n\n'); end
%I = speye(n,n);
%P=I(P,:);
figure;spy(M(P,P)); 
aaa=0;



%figure;spy(M0)
    
%L = chol(M);

