function result = func(lambda, X, Y)
    result = 0
    
    m = size(X, 2)
    
    first = -1/2 * lambda.' * (Y.' * Y) * (X.' * X )* lambda
    first = first + sum(lambda)
end