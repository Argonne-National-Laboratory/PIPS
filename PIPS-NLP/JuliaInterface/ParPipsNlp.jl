
module ParPipsNlp

import MPI

libparpipsnlp=Libdl.dlopen("/home/fqiang/workspace/PIPS/build_pips/PIPS-NLP/libparpipsnlp.so")

type PipsNlpProblemStruct
    ref::Ptr{Void}
    comm::MPI.Comm
    n::Int
    m::Int
    nnodes::Int
    # Callbacks
    str_init_x0::Function
    str_prob_info::Function
    str_eval_f::Function
    str_eval_g::Function
    str_eval_grad_f::Function
    str_eval_jac_g::Function
    str_eval_h::Function

    # For MathProgBase
    sense::Symbol
	status::Int
    
    rowmap::Dict{Int,Int}
    colmap::Dict{Int,Int}
    eq_rowmap::Dict{Int,Int}
    inq_rowmap::Dict{Int,Int}
    
    function PipsNlpProblemStruct(comm, n, m, nnodes, str_init_x0, str_prob_info, str_eval_f, str_eval_g, str_eval_grad_f, str_eval_jac_g, str_eval_h)
        prob = new(C_NULL, comm, n, m, nnodes,
        str_init_x0, str_prob_info, str_eval_f, str_eval_g, str_eval_grad_f, str_eval_jac_g, str_eval_h,
        :Min
        ,0
        ,Dict{Int,Int}()
        ,Dict{Int,Int}()
        ,Dict{Int,Int}()
       	,Dict{Int,Int}()
        )
        # Free the internal PipsNlpProblem structure when
        # the Julia IpoptProblem instance goes out of scope
        finalizer(prob, freeProblemStruct)
        # Return the object we just made
        return prob
    end
end

immutable CallBackData
	prob::Ptr{Void}
	row_node_id::Cint
    col_node_id::Cint
end

export createProblemStruct, solveProblemStruct, freeProblemStruct

###########################################################################
# Callback wrappers
###########################################################################

function str_init_x0_wrapper(x0_ptr::Ptr{Float64}, cbd::Ptr{CallBackData})
	println(" julia - str_init_x0_wrapper ")
    @show cbd
    data = unsafe_load(cbd)
    # @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    @show prob
    # data = unsafe_pointer_to_objref(cbd)::CallBackData
    # out = Array(Ptr{CallBackData},1)
    rowid = data.row_node_id
    colid = data.col_node_id
    assert(rowid == colid)
    n0 = prob.colmap[colid]
    x0 = pointer_to_array(x0_ptr, n0)
    prob.str_init_x0(colid,x0)

    return Int32(1)
end

# prob info (prob_info)
function str_prob_info_wrapper(n_ptr::Ptr{Cint}, col_lb_ptr::Ptr{Float64}, col_ub_ptr::Ptr{Float64}, m_ptr::Ptr{Cint}, row_lb_ptr::Ptr{Float64}, row_ub_ptr::Ptr{Float64}, cbd::Ptr{CallBackData})
    println(" julia - str_prob_info_wrapper ")
    @show cbd
    data = unsafe_load(cbd)
    # @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    @show prob
    # data = unsafe_pointer_to_objref(cbd)::CallBackData
    # out = Array(Ptr{CallBackData},1)
    rowid = data.row_node_id
    colid = data.col_node_id
    assert(rowid == colid)
	
	mode = (col_lb_ptr == C_NULL) ? (:Structure) : (:Values)
	if(mode==:Structure)
		col_lb = pointer_to_array(col_lb_ptr,0)
		col_ub = pointer_to_array(col_ub_ptr,0)
		row_lb = pointer_to_array(row_lb_ptr,0)
		row_ub = pointer_to_array(row_ub_ptr,0)
		(n,m) = prob.str_prob_info(colid,mode,col_lb,col_ub,row_lb,row_ub)
		unsafe_store!(n_ptr,convert(Cint,n)::Cint)
		unsafe_store!(m_ptr,convert(Cint,m)::Cint)
		prob.rowmap[colid] = m
		prob.colmap[colid] = n
	else
		n = unsafe_load(n_ptr)
		m = unsafe_load(m_ptr)
		col_lb = pointer_to_array(col_lb_ptr,n)
		col_ub = pointer_to_array(col_ub_ptr,n)
		row_lb = pointer_to_array(row_lb_ptr,m)
		row_ub = pointer_to_array(row_ub_ptr,m)
		prob.str_prob_info(colid,mode,col_lb,col_ub,row_lb,row_ub)

		neq = 0
		nineq = 0
		for i = 1:length(row_lb)
			if row_lb[i] == row_ub[i]
				neq += 1
			else
				nineq += 1
			end
		end
		assert(neq+nineq == length(row_lb) == m)
		prob.eq_rowmap[colid] = neq
		prob.inq_rowmap[colid] = nineq
	end
	return Int32(1)
end
# Objective (eval_f)
function str_eval_f_wrapper(x0_ptr::Ptr{Float64}, x1_ptr::Ptr{Float64}, obj_ptr::Ptr{Float64}, cbd::Ptr{CallBackData})
    println(" julia - eval_f_wrapper " ); 
    data = unsafe_load(cbd)
    @show data
    # @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    rowid = data.row_node_id
    colid = data.col_node_id
    assert(rowid == colid)
    n0 = prob.colmap[0]
    n1 = prob.colmap[colid]
    # Calculate the new objective
    x0 = pointer_to_array(x0_ptr, n0)
    x1 = pointer_to_array(x1_ptr, n1)
    new_obj = convert(Float64, prob.str_eval_f(colid,x0,x1))::Float64
    # Fill out the pointer
    unsafe_store!(obj_ptr, new_obj)
    # Done
    return Int32(1)
end

# Constraints (eval_g)
function str_eval_g_wrapper(x0_ptr::Ptr{Float64}, x1_ptr::Ptr{Float64}, eq_g_ptr::Ptr{Float64}, inq_g_ptr::Ptr{Float64}, cbd::Ptr{CallBackData})
    println(" julia - eval_g_wrapper " ); 
    data = unsafe_load(cbd)
    # @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    rowid = data.row_node_id
    colid = data.col_node_id
    assert(rowid == colid)
    n0 = prob.colmap[0]
    n1 = prob.colmap[colid]
    x0 = pointer_to_array(x0_ptr, n0)
    x1 = pointer_to_array(x1_ptr, n1)
    # Calculate the new constraint values
    neq = prob.eq_rowmap[rowid]
    nineq = prob.inq_rowmap[rowid]
    new_eq_g = pointer_to_array(eq_g_ptr,neq)
    new_inq_g = pointer_to_array(inq_g_ptr, nineq)
    prob.str_eval_g(colid,x0,x1,new_eq_g,new_inq_g)
    # Done
    return Int32(1)
end

# Objective gradient (eval_grad_f)
function str_eval_grad_f_wrapper(x0_ptr::Ptr{Float64}, x1_ptr::Ptr{Float64}, grad_f_ptr::Ptr{Float64}, cbd::Ptr{CallBackData})
    println(" julia -  eval_grad_f_wrapper  - $(x0_ptr) ,$(x1_ptr)");    
    # Extract Julia the problem from the pointer
    @show cbd
    data = unsafe_load(cbd)
    @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    rowid = data.row_node_id
    colid = data.col_node_id
    n0 = prob.colmap[0]
    n1 = prob.colmap[rowid]
    @show n0,n1
    x0 = pointer_to_array(x0_ptr, n0)
    x1 = pointer_to_array(x1_ptr, n1)
    # Calculate the gradient
    grad_len = prob.colmap[colid]
    new_grad_f = pointer_to_array(grad_f_ptr, grad_len)
    prob.str_eval_grad_f(rowid,colid,x0,x1,new_grad_f)
    if prob.sense == :Max
        new_grad_f *= -1.0
    end
    # Done
    return Int32(1)
end

# Jacobian (eval_jac_g)
function str_eval_jac_g_wrapper(x0_ptr::Ptr{Float64}, x1_ptr::Ptr{Float64}, 
	e_nz_ptr::Ptr{Cint}, e_values_ptr::Ptr{Float64}, e_row_ptr::Ptr{Cint}, e_col_ptr::Ptr{Cint}, 
	i_nz_ptr::Ptr{Cint}, i_values_ptr::Ptr{Float64}, i_row_ptr::Ptr{Cint}, i_col_ptr::Ptr{Cint},  
	cbd::Ptr{CallBackData}
	)
    println(" julia -  eval_jac_g_wrapper " );
    # Extract Julia the problem from the pointer  
    data = unsafe_load(cbd)
    @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    rowid = data.row_node_id
    colid = data.col_node_id
    n0 = prob.colmap[0]
    n1 = prob.colmap[colid]
    x0 = pointer_to_array(x0_ptr, n0)
    x1 = pointer_to_array(x1_ptr, n1)
    nrow = prob.rowmap[rowid]
    ncol = prob.colmap[colid]
    #@show prob
    # Determine mode
    mode = (e_values_ptr == C_NULL && i_values_ptr == C_NULL) ? (:Structure) : (:Values)
    if(mode == :Structure)
    	e_values = pointer_to_array(e_values_ptr,0)
		e_colptr = pointer_to_array(e_col_ptr,0)
		e_rowidx = pointer_to_array(e_row_ptr,0)
		i_values = pointer_to_array(i_values_ptr,0)
		i_colptr = pointer_to_array(i_col_ptr,0)
		i_rowidx = pointer_to_array(i_row_ptr,0)
		(e_nz,i_nz) = prob.str_eval_jac_g(rowid,colid,x0,x1,mode,e_rowidx,e_colptr,e_values,i_rowidx,i_colptr,i_values)
		unsafe_store!(e_nz_ptr,convert(Cint,e_nz)::Cint)
		unsafe_store!(i_nz_ptr,convert(Cint,i_nz)::Cint)
		@show "structure - ",(e_nz,i_nz)
    else
    	e_nz = unsafe_load(e_nz_ptr)
    	e_values = pointer_to_array(e_values_ptr,e_nz)
    	e_rowidx = pointer_to_array(e_row_ptr, e_nz)
    	e_colptr = pointer_to_array(e_col_ptr, ncol+1)
    	i_nz = unsafe_load(i_nz_ptr)
    	@show "values - ",(e_nz,i_nz)
    	i_values = pointer_to_array(i_values_ptr,i_nz)
    	i_rowidx = pointer_to_array(i_row_ptr, i_nz)
    	i_colptr = pointer_to_array(i_col_ptr, ncol+1)
    	prob.str_eval_jac_g(rowid,colid,x0,x1,mode,e_rowidx,e_colptr,e_values,i_rowidx,i_colptr,i_values)
    end
    # Done
    return Int32(1)
end

# Hessian
function str_eval_h_wrapper(x0_ptr::Ptr{Float64}, x1_ptr::Ptr{Float64}, lambda_ptr::Ptr{Float64}, nz_ptr::Ptr{Cint}, values_ptr::Ptr{Float64}, row_ptr::Ptr{Cint}, col_ptr::Ptr{Cint}, cbd::Ptr{CallBackData})
    println(" julia - eval_h_wrapper " ); 
    # Extract Julia the problem from the pointer
    data = unsafe_load(cbd)
    @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    rowid = data.row_node_id
    colid = data.col_node_id
    
    high = max(rowid,colid)
    low  = min(rowid,colid)
    n0 = prob.colmap[0]
    n1 = prob.colmap[high]
    x0 = pointer_to_array(x0_ptr, n0)
    x1 = pointer_to_array(x1_ptr, n1)
    ncol = prob.colmap[low]
    g0 = prob.rowmap[high]
    lambda = pointer_to_array(lambda_ptr, g0)
    obj_factor = 1.0
    if prob.sense == :Max
        obj_factor *= -1.0
    end
    # Did the user specify a Hessian
    mode = (values_ptr == C_NULL) ? (:Structure) : (:Values)
    if(mode == :Structure)
    	values = pointer_to_array(values_ptr,0)
		colptr = pointer_to_array(col_ptr,0)
		rowidx = pointer_to_array(row_ptr,0)
		nz = prob.str_eval_h(rowid,colid,x0,x1,obj_factor,lambda,mode,rowidx,colptr,values)
		unsafe_store!(nz_ptr,convert(Cint,nz)::Cint)
		@show "structure - ", nz
    else
    	nz = unsafe_load(nz_ptr)
    	values = pointer_to_array(values_ptr, nz)
    	rowidx = pointer_to_array(row_ptr, nz)
    	colptr = pointer_to_array(col_ptr, ncol+1)
    	@show "value - ", nz
    	prob.str_eval_h(rowid,colid,x0,x1,obj_factor,lambda,mode,rowidx,colptr,values)
    end
    # Done
    return Int32(1)
end

###########################################################################
# C function wrappers
###########################################################################
function createProblemStruct(comm::MPI.Comm, nnodes::Int, n::Int,m::Int, 
    str_init_x0, str_prob_info, str_eval_f, str_eval_g, str_eval_grad_f, str_eval_jac_g, str_eval_h)
	println(" createProblemStruct  -- julia")
	str_init_x0_cb = cfunction(str_init_x0_wrapper, Cint, (Ptr{Float64}, Ptr{CallBackData}) )
    str_prob_info_cb = cfunction(str_prob_info_wrapper, Cint, (Ptr{Cint}, Ptr{Float64}, Ptr{Float64}, Ptr{Cint}, Ptr{Float64}, Ptr{Float64}, Ptr{CallBackData}) )
    str_eval_f_cb = cfunction(str_eval_f_wrapper,Cint, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{CallBackData}) )
    str_eval_g_cb = cfunction(str_eval_g_wrapper,Cint, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{CallBackData}) )
    str_eval_grad_f_cb = cfunction(str_eval_grad_f_wrapper, Cint, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{CallBackData}) )
    str_eval_jac_g_cb = cfunction(str_eval_jac_g_wrapper, Cint, (Ptr{Float64}, Ptr{Float64}, 
    	Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, 
    	Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, 
    	Ptr{CallBackData}))
    str_eval_h_cb = cfunction(str_eval_h_wrapper, Cint, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{CallBackData}))
    
    println(" callback created ")
    prob = PipsNlpProblemStruct(comm, nnodes,n, m, str_init_x0, str_prob_info, str_eval_f, str_eval_g, str_eval_grad_f, str_eval_jac_g, str_eval_h)
    @show prob
    ret = ccall((:CreatePipsNlpProblemStruct,:libparpipsnlp),Ptr{Void},
            (MPI.Comm, 
            Cint, 
            Cint,
            Cint,
            Ptr{Void}, 
            Ptr{Void}, Ptr{Void}, Ptr{Void},
            Ptr{Void}, Ptr{Void}, Ptr{Void}
            ,Any
            ),
            comm, 
            nnodes,
            n,
            m, 
            str_init_x0_cb,
            str_prob_info_cb,
            str_eval_f_cb, 
            str_eval_g_cb,
            str_eval_grad_f_cb, 
            str_eval_jac_g_cb, 
            str_eval_h_cb,
            prob
            )
    println(" ccall CreatePipsNlpProblemStruct done ")
    @show ret   
    
    if ret == C_NULL
        error("PIPS-NLP: Failed to construct problem.")
    else
        prob.ref = ret
    end
    @show prob
    println("end createProblemStruct - julia")
    return prob
end

function solveProblemStruct(prob::PipsNlpProblemStruct)
	println("solveProblemStruct - julia")
    @show prob
    
    ret = ccall((:PipsNlpSolveStruct, :libparpipsnlp), Cint, 
            (Ptr{Void},),
            prob.ref)
    
    prob.status = Int(ret)

    return prob.status
end

function freeProblemStruct(prob::PipsNlpProblemStruct)
    @show "freeProblemStruct"
    ret = ccall((:FreePipsNlpProblemStruct, :libparpipsnlp),
            Void, (Ptr{Void},),
            prob.ref)
    @show ret
    return ret
end

# include("PipsNlpSolverInterface.jl")

end



