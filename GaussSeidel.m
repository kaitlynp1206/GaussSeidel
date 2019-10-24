%Gauss Seidel Algorithm
%x parameter = initial guesses
%lambda = relaxation parameter
function x = GaussSeidel(A, B, x, lambda, convergence)
    error = 100;
    iterations = 0;
    %Iterate through values until error (et) meets convergence requirements
    while(error > convergence)
        prev = x;
        for i = 1:size(A,1)
            x(i) = B(i);
            for j = 1:size(A,1)
                %check for diff   
                if(j ~= i)
                    x(i) = x(i)-(A(i,j)*x(j));
                end
            end
            x(i) = lambda*x(i)/A(i,i)+(1-lambda)*prev(i);
        end
        error = norm((x-prev)/x);
        iterations = iterations + 1;
    end
    fprintf("relaxation parameter=%.2f results in %d iterations of Gauss Seidel\n", lambda, iterations);
end