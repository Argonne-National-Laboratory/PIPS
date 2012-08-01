function hsl_ma86_test()
%
% Unit tests for hsl_ma86 matlab interface
%
fails = 0;

fprintf('Testing Poisson(2) real:\n')
A = gallery('poisson', 2);
fails = fails + test_with_matrix(A);

fprintf('Testing toy example real:\n')
A = sparse ([1 1 1 2 2 3 3 3 4 4], [2 3 4  1 3  1 2 3  1 4], [1.1 2.2 3.3, 1.1 4.4, 2.2 4.4 5.5, 3.3 6.6]);
fails = fails + test_with_matrix(A);

fprintf('Testing toy example complex Hermitian:\n')
A = sparse ([1 1 1 2 2 3 3 3 4 4], [2 3 4  1 3  1 2 3  1 4], [1.1-i 2.2-i 3.3-i, 1.1+i 4.4-i, 2.2+i 4.4+i 5.5, 3.3+i 6.6]);
fails = fails + test_with_matrix(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(fails == 0)
   fprintf('Test OK.\n')
else
   fprintf('Failed %i tests.\n', fails)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fails = test_with_matrix(A)
% Run through all tests with a specified matrix
fails = 0;

control.nemin = 8;
if(isreal(A))
   x = rand(size(A,1));
   x2 = rand(size(A,1));
else
   x = rand(size(A,1)) + rand(size(A,1))*i;
   x2 = rand(size(A,1)) + rand(size(A,1))*i;
   if(A(1,2) == conj(A(2,1)))
      control.hermitian = true;
   end
end
b = A*x;
b2 = A*x2;

fprintf('   - mc68 ordering, separate calls\n')
[handleA, info] = hsl_ma86_factor(A, control);
soln = hsl_ma86_solve(handleA, b);
res = norm(A*soln - b, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
hsl_ma86_destroy(handleA);

fprintf('   - symamd ordering, separate calls\n')
handleA = hsl_ma86_factor(A, control, symamd(A));
soln = hsl_ma86_solve(handleA, b);
res = norm(A*soln - b, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
hsl_ma86_destroy(handleA);

fprintf('   - symamd ordering, backslash with handle\n')
[soln, info, handleA] = hsl_ma86_backslash(A, b, control, symamd(A));
res = norm(A*soln - b, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end

fprintf('   - subsequent solve\n')
[soln, info] = hsl_ma86_solve(handleA, b2, control);
res = norm(A*soln - b2, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b2, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
hsl_ma86_destroy(handleA);

fprintf('   - simple backslash\n')
soln = hsl_ma86_backslash(A, b, control);
res = norm(A*soln - b, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
