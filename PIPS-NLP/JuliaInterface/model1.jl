include("PipsNlp.jl")

using PipsNlp

function eval_f(x)
    x1 = x[1]
    x2 = x[2]
    return (1.0 - x1)^2 + 100 *(x2-x1^2)^2
end

function eval_g(x,g)
    x1 = x[1]
    x2 = x[2]
    
    g[1] = x1+x2
end

function eval_grad_f(x, grad_f)
    x1 = x[1]
    x2 = x[2]

    grad_f[1] = (-2 * (1 - x1) + 100 * (2 * -(2x1) * (x2 - x1 ^ 2)))
    grad_f[2] = (100 * (2 * (x2 - x1 ^ 2)))  
end

function eval_jac_g(x, mode, irow, kcol, values)
    println("julia - eval_jac_g - ",mode)
    @show irow
    @show kcol
    
    x1 = x[1]
    x2 = x[2]

    if(mode == :Structure)
        val = ones(Float64,length(Ijac))
        mat = sparse(Ijac,Jjac,val)
        @show mat
        array_copy(mat.rowval,irow)
        array_copy(mat.colptr,kcol)
        @show irow, kcol
    elseif(mode == :Values)
        values[1] = 1
        values[2] = 1
    end
end

function array_copy(src,dest)
    assert(length(src)==length(dest))
    for i in 1:length(src)
        dest[i] = src[i]-1
    end    
end

function eval_h(x,mode,irow,kcol,obj_factor,lambda, values)
    x1 = x[1]
    x2 = x[2]
    
    if(mode == :Structure)
        val = ones(Float64,length(Ihess))
        mat = sparse(Ihess,Jhess,val)
        array_copy(mat.rowval,irow)
        array_copy(mat.colptr,kcol)
    else
        values[1] = obj_factor*(2 + 100 * (-4 * (x2 - x1 ^ 2) + 2 * -(2x1) * -(2x1)))   
        values[2] = obj_factor*(100 * (2 * -(2x1)))
        values[3] = obj_factor*200
    end        
end

Ijac = [1,1]
Jjac = [1,2]
Ihess = [1,1,2]
Jhess = [1,2,2]
n = 2
m = 1
x_L = fill(typemin(Float64),n)
x_U = fill(typemax(Float64),n)
g_L = fill(100.0,m)
g_U = fill(100.0,m)

prob = createProblem(n,m,
    x_L,x_U,
    g_L,g_U,
    length(Ijac),length(Ihess),
    eval_f,eval_g,eval_grad_f,eval_jac_g, eval_h)



ret = solveProblem(prob)

println(ret)

ret = freeProblem(prob)

println(ret)



