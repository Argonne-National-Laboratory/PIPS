module ParPipsNlp

import MPI

try
  sharedLib=ENV["PIPS_NLP_PAR_SHARED_LIB"]
  #explicitly check if the file exists (because dlopen sometimes does not throw an error for invalid filenames, resulting in a seg fault)
  if(!isfile(sharedLib))
    error(string("The specified shared library ([", sharedLib, "]) does not exist"))
  end  
  const libparpipsnlp=Libdl.dlopen(get(ENV,"PIPS_NLP_PAR_SHARED_LIB",""))
catch 
  warn("Could not load PIPS-NLP shared library. Make sure the ENV variable 'PIPS_NLP_PAR_SHARED_LIB' points to its location, usually in the PIPS repo at PIPS/build_pips/PIPS-NLP/libparpipsnlp.so")
  rethrow()
end


# Struct Model interface
abstract ModelInterface

#######################

type FakeModel <: ModelInterface
    sense::Symbol
    status::Int
    nscen::Int
    rowmap::Dict{Int,Int}
    colmap::Dict{Int,Int}
    eq_rowmap::Dict{Int,Int}
    inq_rowmap::Dict{Int,Int}

    get_num_scen::Function
    get_sense::Function

    get_status::Function
    get_num_rows::Function
    get_num_cols::Function
    get_num_eq_cons::Function
    get_num_ineq_cons::Function

    set_status::Function
    set_num_rows::Function
    set_num_cols::Function
    set_num_eq_cons::Function
    set_num_ineq_cons::Function

    str_init_x0::Function
    str_prob_info::Function
    str_eval_f::Function
    str_eval_g::Function
    str_eval_grad_f::Function
    str_eval_jac_g::Function
    str_eval_h::Function
           

    function FakeModel(sense::Symbol,status::Int,nscen::Int, str_init_x0, str_prob_info, str_eval_f, str_eval_g,str_eval_grad_f,str_eval_jac_g, str_eval_h)
        instance = new(sense,status,nscen,Dict{Int,Int}(),Dict{Int,Int}(),Dict{Int,Int}(),Dict{Int,Int}())
        instance.str_init_x0 = str_init_x0
        instance.str_prob_info = str_prob_info 
        instance.str_eval_f = str_eval_f
        instance.str_eval_g = str_eval_g
        instance.str_eval_grad_f = str_eval_grad_f
        instance.str_eval_jac_g = str_eval_jac_g
        instance.str_eval_h = str_eval_h


        instance.get_num_scen = function()
            return instance.nscen
        end
        instance.get_sense = function()
            return instance.sense
        end
        instance.get_status = function()
            return instance.status
        end
        instance.get_num_rows = function(id::Integer)
            return instance.rowmap[id]
        end
        instance.get_num_cols = function(id::Integer)
            return instance.colmap[id]
        end
        instance.get_num_eq_cons = function(id::Integer)
            return instance.eq_rowmap[id]
        end
        instance.get_num_ineq_cons = function(id::Integer)
            return instance.inq_rowmap[id]
        end
        instance.set_status = function(s::Integer)
            instance.status = s
        end
        instance.set_num_rows = function(id::Integer, v::Integer)
            return instance.rowmap[id] = v
        end
        instance.set_num_cols = function(id::Integer, v::Integer)
            return instance.colmap[id] = v
        end
        instance.set_num_eq_cons = function(id::Integer, v::Integer)
            return instance.eq_rowmap[id] = v
        end
        instance.set_num_ineq_cons = function(id::Integer, v::Integer)
            return instance.inq_rowmap[id] = v
        end
        return instance
    end
end
 


#######################


type PipsNlpProblemStruct
    ref::Ptr{Void}
    model::ModelInterface
    comm::MPI.Comm
    
    function PipsNlpProblemStruct(comm, model)
        prob = new(C_NULL, model, comm)
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

export  ModelInterface, FakeModel,
        createProblemStruct, solveProblemStruct, freeProblemStruct

###########################################################################
# Callback wrappers
###########################################################################

function str_init_x0_wrapper(x0_ptr::Ptr{Float64}, cbd::Ptr{CallBackData})
	# @show " julia - str_init_x0_wrapper "
    # @show cbd
    data = unsafe_load(cbd)
    # @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    # @show prob
    # data = unsafe_pointer_to_objref(cbd)::CallBackData
    # out = Array(Ptr{CallBackData},1)
    rowid = data.row_node_id
    colid = data.col_node_id
    assert(rowid == colid)
    n0 = prob.model.get_num_cols(colid)
    x0 = pointer_to_array(x0_ptr, n0)
    prob.model.str_init_x0(colid,x0)

    return Int32(1)
end

# prob info (prob_info)
function str_prob_info_wrapper(n_ptr::Ptr{Cint}, col_lb_ptr::Ptr{Float64}, col_ub_ptr::Ptr{Float64}, m_ptr::Ptr{Cint}, row_lb_ptr::Ptr{Float64}, row_ub_ptr::Ptr{Float64}, cbd::Ptr{CallBackData})
    # @show " julia - str_prob_info_wrapper "
    # @show cbd
    data = unsafe_load(cbd)
    # @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    # @show prob
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
		(n,m) = prob.model.str_prob_info(colid,mode,col_lb,col_ub,row_lb,row_ub)
		unsafe_store!(n_ptr,convert(Cint,n)::Cint)
		unsafe_store!(m_ptr,convert(Cint,m)::Cint)
        # @show typeof(colid), typeof(m)
		prob.model.set_num_rows(colid, m)
		prob.model.set_num_cols(colid, n)
	else
		n = unsafe_load(n_ptr)
		m = unsafe_load(m_ptr)
		col_lb = pointer_to_array(col_lb_ptr,n)
		col_ub = pointer_to_array(col_ub_ptr,n)
		row_lb = pointer_to_array(row_lb_ptr,m)
		row_ub = pointer_to_array(row_ub_ptr,m)
		prob.model.str_prob_info(colid,mode,col_lb,col_ub,row_lb,row_ub)

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
		prob.model.set_num_eq_cons(colid,neq)
		prob.model.set_num_ineq_cons(colid,nineq) 
	end
	return Int32(1)
end
# Objective (eval_f)
function str_eval_f_wrapper(x0_ptr::Ptr{Float64}, x1_ptr::Ptr{Float64}, obj_ptr::Ptr{Float64}, cbd::Ptr{CallBackData})
    # @show " julia - eval_f_wrapper "
    data = unsafe_load(cbd)
    # @show data
    # @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    rowid = data.row_node_id
    colid = data.col_node_id
    assert(rowid == colid)
    n0 = prob.model.get_num_cols(0)
    n1 = prob.model.get_num_cols(colid)
    # Calculate the new objective
    x0 = pointer_to_array(x0_ptr, n0)
    x1 = pointer_to_array(x1_ptr, n1)
    new_obj = convert(Float64, prob.model.str_eval_f(colid,x0,x1))::Float64
    # Fill out the pointer
    unsafe_store!(obj_ptr, new_obj)
    # Done
    return Int32(1)
end

# Constraints (eval_g)
function str_eval_g_wrapper(x0_ptr::Ptr{Float64}, x1_ptr::Ptr{Float64}, eq_g_ptr::Ptr{Float64}, inq_g_ptr::Ptr{Float64}, cbd::Ptr{CallBackData})
    # @show " julia - eval_g_wrapper " 
    data = unsafe_load(cbd)
    # @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    rowid = data.row_node_id
    colid = data.col_node_id
    assert(rowid == colid)
    n0 = prob.model.get_num_cols(0)
    n1 = prob.model.get_num_cols(colid)
    x0 = pointer_to_array(x0_ptr, n0)
    x1 = pointer_to_array(x1_ptr, n1)
    # Calculate the new constraint values
    neq = prob.model.get_num_eq_cons(rowid)
    nineq = prob.model.get_num_ineq_cons(rowid)
    new_eq_g = pointer_to_array(eq_g_ptr,neq)
    new_inq_g = pointer_to_array(inq_g_ptr, nineq)
    prob.model.str_eval_g(colid,x0,x1,new_eq_g,new_inq_g)
    # Done
    return Int32(1)
end

# Objective gradient (eval_grad_f)
function str_eval_grad_f_wrapper(x0_ptr::Ptr{Float64}, x1_ptr::Ptr{Float64}, grad_f_ptr::Ptr{Float64}, cbd::Ptr{CallBackData})
    # @show " julia -  eval_grad_f_wrapper "  
    # Extract Julia the problem from the pointer
    # @show cbd
    data = unsafe_load(cbd)
    # @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    rowid = data.row_node_id
    colid = data.col_node_id
    n0 = prob.model.get_num_cols(0)
    n1 = prob.model.get_num_cols(rowid)
    # @show n0,n1
    x0 = pointer_to_array(x0_ptr, n0)
    x1 = pointer_to_array(x1_ptr, n1)
    # Calculate the gradient
    grad_len = prob.model.get_num_cols(colid)
    new_grad_f = pointer_to_array(grad_f_ptr, grad_len)
    prob.model.str_eval_grad_f(rowid,colid,x0,x1,new_grad_f)
    if prob.model.get_sense() == :Max
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
    # @show " julia -  eval_jac_g_wrapper " 
    # Extract Julia the problem from the pointer  
    data = unsafe_load(cbd)
    # @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    rowid = data.row_node_id
    colid = data.col_node_id
    n0 = prob.model.get_num_cols(0)
    n1 = prob.model.get_num_cols(rowid) #we can do this because of 2-level and no linking constraint
    # @show n0, n1 
    x0 = pointer_to_array(x0_ptr, n0)
    x1 = pointer_to_array(x1_ptr, n1)
    # @show x0
    # @show x1 
    nrow = prob.model.get_num_rows(rowid) 
    ncol = prob.model.get_num_cols(colid) 
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
        (e_nz,i_nz) = prob.model.str_eval_jac_g(rowid,colid,x0,x1,mode,e_rowidx,e_colptr,e_values,i_rowidx,i_colptr,i_values)
		unsafe_store!(e_nz_ptr,convert(Cint,e_nz)::Cint)
		unsafe_store!(i_nz_ptr,convert(Cint,i_nz)::Cint)
		# @show "structure - ",(e_nz,i_nz)
    else
    	e_nz = unsafe_load(e_nz_ptr)
    	e_values = pointer_to_array(e_values_ptr,e_nz)
    	e_rowidx = pointer_to_array(e_row_ptr, e_nz)
    	e_colptr = pointer_to_array(e_col_ptr, ncol+1)
    	i_nz = unsafe_load(i_nz_ptr)
    	# @show "values - ",(e_nz,i_nz), ncol
    	i_values = pointer_to_array(i_values_ptr,i_nz)
    	i_rowidx = pointer_to_array(i_row_ptr, i_nz)
    	i_colptr = pointer_to_array(i_col_ptr, ncol+1)
        # @show x0
        # @show x1 
    	prob.model.str_eval_jac_g(rowid,colid,x0,x1,mode,e_rowidx,e_colptr,e_values,i_rowidx,i_colptr,i_values)
    end
    # Done
    return Int32(1)
end

# Hessian
function str_eval_h_wrapper(x0_ptr::Ptr{Float64}, x1_ptr::Ptr{Float64}, lambda_ptr::Ptr{Float64}, nz_ptr::Ptr{Cint}, values_ptr::Ptr{Float64}, row_ptr::Ptr{Cint}, col_ptr::Ptr{Cint}, cbd::Ptr{CallBackData})
    # @show " julia - eval_h_wrapper " 
    # Extract Julia the problem from the pointer
    data = unsafe_load(cbd)
    # @show data
    userdata = data.prob
    prob = unsafe_pointer_to_objref(userdata)::PipsNlpProblemStruct
    rowid = data.row_node_id
    colid = data.col_node_id
    
    high = max(rowid,colid)
    low  = min(rowid,colid)
    n0 = prob.model.get_num_cols(0) 
    n1 = prob.model.get_num_cols(high)
    x0 = pointer_to_array(x0_ptr, n0)
    x1 = pointer_to_array(x1_ptr, n1)
    # @show x0
    # @show x1
    ncol = prob.model.get_num_cols(low)
    g0 = prob.model.get_num_rows(high) 
    # @show g0
    # @show ncol
    lambda = pointer_to_array(lambda_ptr, g0)
    obj_factor = 1.0
    if prob.model.get_sense() == :Max
        obj_factor *= -1.0
    end
    # Did the user specify a Hessian
    mode = (values_ptr == C_NULL) ? (:Structure) : (:Values)
    if(mode == :Structure)
    	values = pointer_to_array(values_ptr,0)
		colptr = pointer_to_array(col_ptr,0)
		rowidx = pointer_to_array(row_ptr,0)
		nz = prob.model.str_eval_h(rowid,colid,x0,x1,obj_factor,lambda,mode,rowidx,colptr,values)
		unsafe_store!(nz_ptr,convert(Cint,nz)::Cint)
		# @show "structure - ", nz
    else
    	nz = unsafe_load(nz_ptr)
    	values = pointer_to_array(values_ptr, nz)
    	rowidx = pointer_to_array(row_ptr, nz)
    	colptr = pointer_to_array(col_ptr, ncol+1)
    	# @show "value - ", nz
    	prob.model.str_eval_h(rowid,colid,x0,x1,obj_factor,lambda,mode,rowidx,colptr,values)
    end
    # Done
    return Int32(1)
end

###########################################################################
# C function wrappers
###########################################################################
function createProblemStruct(comm::MPI.Comm, model::ModelInterface)
	# println(" createProblemStruct  -- julia")
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
    
    # println(" callback created ")
    prob = PipsNlpProblemStruct(comm, model)
    # @show prob
    ret = ccall((:CreatePipsNlpProblemStruct,:libparpipsnlp),Ptr{Void},
            (MPI.Comm, 
            Cint, 
            Ptr{Void}, 
            Ptr{Void}, Ptr{Void}, Ptr{Void},
            Ptr{Void}, Ptr{Void}, Ptr{Void}
            ,Any
            ),
            comm, 
            model.get_num_scen(),
            str_init_x0_cb,
            str_prob_info_cb,
            str_eval_f_cb, 
            str_eval_g_cb,
            str_eval_grad_f_cb, 
            str_eval_jac_g_cb, 
            str_eval_h_cb,
            prob
            )
    # println(" ccall CreatePipsNlpProblemStruct done ")
    # @show ret   
    
    if ret == C_NULL
        error("PIPS-NLP: Failed to construct problem.")
    else
        prob.ref = ret
    end
    # @show prob
    # println("end createProblemStruct - julia")
    return prob
end

function solveProblemStruct(prob::PipsNlpProblemStruct)
	# println("solveProblemStruct - julia")
    # @show prob
    
    ret = ccall((:PipsNlpSolveStruct,:libparpipsnlp), Cint, 
            (Ptr{Void},),
            prob.ref)
    
    prob.model.set_status(Int(ret))

    return prob.model.get_status()
end

function freeProblemStruct(prob::PipsNlpProblemStruct)
    # @show "freeProblemStruct"
    ret = ccall((:FreePipsNlpProblemStruct,:libparpipsnlp),
            Void, (Ptr{Void},),
            prob.ref)
    # @show ret
    return ret
end

# include("PipsNlpSolverInterface.jl")

end



