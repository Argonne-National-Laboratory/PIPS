
module PipsNlp

#update this to the path of the PIPS-NLP serial share library. 
libpipsnlp=Libdl.dlopen("/home/fqiang/workspace/PIPS/build_pips/PIPS-NLP/libpipsnlp.so")

type ProblemData
    
end

type UserData

end

export createProblem, solveProblem, freeProblem

type PipsNlpProblem
    ref::Ptr{Void}
    n::Int
    m::Int
    x::Vector{Float64}
    g::Vector{Float64}
    obj_val::Float64
    status::Int

    # Callbacks
    eval_f::Function
    eval_g::Function
    eval_grad_f::Function
    eval_jac_g::Function
    eval_h  # Can be nothing
    
    #jac , hess
    nzJac::Int
    nzHess::Int

    # For MathProgBase
    sense::Symbol

    
    function PipsNlpProblem(ref::Ptr{Void}, n, m, eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h, nzJac, nzHess)
        prob = new(ref, n, m, zeros(Float64, n), zeros(Float64, m), 0.0, 0,
        eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h, nzJac, nzHess,
        :Min)
        
        # Free the internal PipsNlpProblem structure when
        # the Julia IpoptProblem instance goes out of scope
        finalizer(prob, freeProblem)
        # Return the object we just made
        prob
    end
end


###########################################################################
# Callback wrappers
###########################################################################
# Objective (eval_f)
function eval_f_wrapper(x_ptr::Ptr{Float64}, obj_ptr::Ptr{Float64}, user_data::Ptr{Void})
    println(" julia - eval_f_wrapper " ); 
    # Extract Julia the problem from the pointer
    prob = unsafe_pointer_to_objref(user_data)::PipsNlpProblem
    # Calculate the new objective
    new_obj = convert(Float64, prob.eval_f(pointer_to_array(x_ptr, prob.n)))::Float64
    # Fill out the pointer
    unsafe_store!(obj_ptr, new_obj)
    # Done
    return Int32(1)
end

# Constraints (eval_g)
function eval_g_wrapper(x_ptr::Ptr{Float64}, g_ptr::Ptr{Float64}, user_data::Ptr{Void})
    println(" julia - eval_g_wrapper " ); 
    # Extract Julia the problem from the pointer
    prob = unsafe_pointer_to_objref(user_data)::PipsNlpProblem
    # Calculate the new constraint values
    new_g = pointer_to_array(g_ptr, prob.m)
    prob.eval_g(pointer_to_array(x_ptr, prob.n), new_g)
    # Done
    return Int32(1)
end

# Objective gradient (eval_grad_f)
function eval_grad_f_wrapper(x_ptr::Ptr{Float64}, grad_f_ptr::Ptr{Float64}, user_data::Ptr{Void})
    println(" julia -  eval_grad_f_wrapper " );    
    # Extract Julia the problem from the pointer
    prob = unsafe_pointer_to_objref(user_data)::PipsNlpProblem
    # Calculate the gradient
    new_grad_f = pointer_to_array(grad_f_ptr, Int(prob.n))
    prob.eval_grad_f(pointer_to_array(x_ptr, Int(prob.n)), new_grad_f)
    if prob.sense == :Max
        new_grad_f *= -1.0
    end
    # Done
    return Int32(1)
end

# Jacobian (eval_jac_g)
function eval_jac_g_wrapper(x_ptr::Ptr{Float64}, values_ptr::Ptr{Float64}, iRow::Ptr{Cint}, jCol::Ptr{Cint},  user_data::Ptr{Void})
    println(" julia -  eval_jac_g_wrapper " );
    # Extract Julia the problem from the pointer  
    #@show user_data  
    prob = unsafe_pointer_to_objref(user_data)::PipsNlpProblem
    #@show prob
    # Determine mode
    mode = (values_ptr == C_NULL) ? (:Structure) : (:Values)
    x = pointer_to_array(x_ptr, prob.n)
    irows = pointer_to_array(iRow, Int(prob.nzJac))
    kcols = pointer_to_array(jCol, Int(prob.n+1))
    values = pointer_to_array(values_ptr, Int(prob.nzJac))
    prob.eval_jac_g(x, mode, irows, kcols, values)
    # Done
    return Int32(1)
end

# Hessian
function eval_h_wrapper(x_ptr::Ptr{Float64}, lambda_ptr::Ptr{Float64}, values_ptr::Ptr{Float64}, iRow::Ptr{Cint}, jCol::Ptr{Cint}, user_data::Ptr{Void})
    println(" julia - eval_h_wrapper " ); 
    # Extract Julia the problem from the pointer
    prob = unsafe_pointer_to_objref(user_data)::PipsNlpProblem
    # Did the user specify a Hessian
    if prob.eval_h === nothing
        # No Hessian provided
        return Int32(0)
    else
        # Determine mode
        mode = (values_ptr == C_NULL) ? (:Structure) : (:Values)
        x = pointer_to_array(x_ptr, prob.n)
        lambda = pointer_to_array(lambda_ptr, prob.m)
        @show prob.nzHess
        @show prob.n
        irows = pointer_to_array(iRow, Int(prob.nzHess))
        kcols = pointer_to_array(jCol, Int(prob.n+1))
        values = pointer_to_array(values_ptr, Int(prob.nzHess))
        obj_factor = 1.0
        if prob.sense == :Max
            obj_factor *= -1.0
        end
        prob.eval_h(x, mode, irows, kcols, obj_factor, lambda, values)
        # Done
        return Int32(1)
    end
end

###########################################################################
# C function wrappers
###########################################################################
function createProblem(n::Int,m::Int,
    x_L::Vector{Float64},x_U::Vector{Float64},
    g_L::Vector{Float64},g_U::Vector{Float64},
    nzJac::Int, nzHess::Int,
    eval_f, eval_g,eval_grad_f,eval_jac_g,eval_h = nothing)

    @assert n == length(x_L) == length(x_U)
    @assert m == length(g_L) == length(g_U)
    eval_f_cb = cfunction(eval_f_wrapper,Cint, (Ptr{Float64}, Ptr{Float64}, Ptr{Void}) )
    eval_g_cb = cfunction(eval_g_wrapper,Cint, (Ptr{Float64}, Ptr{Float64}, Ptr{Void}) )
    eval_grad_f_cb = cfunction(eval_grad_f_wrapper, Cint, (Ptr{Float64}, Ptr{Float64}, Ptr{Void}) )
    eval_jac_g_cb = cfunction(eval_jac_g_wrapper, Cint, (Ptr{Float64}, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Void}))
    eval_h_cb = cfunction(eval_h_wrapper, Cint, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Void}))
    
    ret = ccall((:CreatePipsNlpProblem,:libpipsnlp),Ptr{Void},
            (Cint, Cint,
            Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64},
            Cint, Cint,
            Ptr{Void},Ptr{Void},
            Ptr{Void},Ptr{Void}, Ptr{Void}),
            n,m,
            x_L,x_U,
            g_L,g_U,
            nzJac, nzHess,
            eval_f_cb, eval_g_cb,
            eval_grad_f_cb, eval_jac_g_cb, eval_h_cb
            )
    println(" ccall CreatePipsNlpProblem done ")
    @show ret   
    
    if ret == C_NULL
        error("PIPS-NLP: Failed to construct problem.")
    else
        return(PipsNlpProblem(ret, n, m, eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h, nzJac, nzHess))
    end
end

function solveProblem(prob::PipsNlpProblem)
    @show "solveProblem"    
    @show prob
    
    final_objval = [0.0]
    ret = ccall((:PipsNlpSolve, :libpipsnlp), Cint, 
            (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Any),
            prob.ref, final_objval, prob.x, prob)
    prob.obj_val = final_objval[1]
    prob.status = Int(ret)

    return prob.status
end

function freeProblem(prob::PipsNlpProblem)
    @show "freeProblem"
    ret = ccall((:FreePipsNlpProblem, :libpipsnlp),
            Void, (Ptr{Void},),
            prob.ref)
    @show ret
    return ret
end

end



