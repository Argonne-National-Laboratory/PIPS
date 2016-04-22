import MPI

#min x1^2 + x3^2 + x4^2 + x5^2 + x6^2
# st.
#     x3 + x4 = 200
#     x5 + x6 = 300
#     x1 = x3
#     x1 = x5
# x free variables

include("../../JuliaInterface/ParPipsNlp.jl")

using ParPipsNlp

function str_init_x0(nodeid, x0)
    if(nodeid == 0)
        x0[1] = 0.0
    elseif nodeid == 1
        x0[1] = 0.0
        x0[2] = 0.0
    elseif nodeid == 2
        x0[1] = 0.0
        x0[2] = 0.0
    end
end

function str_prob_info(nodeid,mode,col_lb,col_ub,row_lb,row_ub)
    println("str_prob_info")             
    if(mode == :Values)
        #setting col lb ub, row lb ub
        if(nodeid==0)
            assert(length(row_lb) == length(row_ub))
            assert(length(col_lb) == length(col_ub))
            fill!(col_lb,-Inf)
            fill!(col_ub, Inf)
        elseif (nodeid == 1 )
            assert(length(row_lb) == length(row_ub))
            assert(length(col_lb) == length(col_ub))
            fill!(col_lb, -Inf)
            fill!(col_ub, Inf)
        	row_lb[1] = 200
            row_ub[1] = 200
            row_lb[2] =  0
            row_ub[2] = 0
        elseif (nodeid == 2)
            assert(length(row_lb) == length(row_ub))
            assert(length(col_lb) == length(col_ub))
            fill!(col_lb, -Inf)
            fill!(col_ub, Inf)
            row_lb[1] = 300
            row_ub[1] = 300
            row_lb[2] =  0
            row_ub[2] = 0
        else
            assert(false)
        end  
	    @show nodeid
	    @show mode
        @show col_lb
        @show col_ub
        @show row_lb
        @show row_ub  
    end

    if(nodeid==0)
        return (1,0)
    elseif(nodeid == 1 || nodeid == 2)
        return (2,2)
    else
        assert(false)
    end
end

function str_eval_f(nodeid,x0,x1)
  println("str_eval_f")
  @show nodeid
  @show x0
  @show x1
    fval = 0.0
    if nodeid ==  0
        x11 = x0[1]
        fval = x11*x11 
    elseif nodeid == 1
        x33 = x1[1]
        x44 = x1[2]
        fval = x33*x33 + x44*x44
    elseif nodeid == 2
        x55 = x1[1]
        x66 = x1[1]
        fval = x55*x55 + x66*x66
    else
        assert(false)
    end
    @show fval
    return fval
end

function str_eval_g(nodeid,x0,x1,new_eq_g, new_inq_g)
     println("str_eval_g")
    @show nodeid
    @show x0
    @show x1
    x11 = x0[1]
    println("eval_g  ", nodeid)
    println(x0)
    println(x1)
    println(new_eq_g)
    println(new_inq_g)
    if nodeid == 0
        assert(x0 == x1)
        assert(length(new_inq_g) == 0 )
        assert(length(new_eq_g) == 0)
    elseif nodeid == 1
        x33 = x1[1]
        x44 = x1[2]
    	assert(length(new_inq_g) == 0 )
        assert(length(new_eq_g) == 2)
    	new_eq_g[1] = x33+x44
        new_eq_g[2] = x11-x33
    elseif nodeid == 2
        x55 = x1[1]
        x66 = x1[2]
        assert(length(new_inq_g) == 0 )
        assert(length(new_eq_g) == 2)
    	new_eq_g[1] = x55+x66
        new_eq_g[2] = x11-x55
    else
        assert(false)
    end
  	@show new_eq_g
  	@show new_inq_g
end


function str_eval_grad_f(rowid,colid,x0,x1,new_grad_f)
    println("str_eval_grad_f")   
    @show rowid,colid
    @show x0
    @show x1
    x11 = x0[1]
    if (0,0) == (rowid, colid)
        assert(length(new_grad_f) == 1)
        assert(x1 == x0)
        new_grad_f[1] = 2.0*x11
     elseif (1,1) == (rowid, colid)
        assert(length(new_grad_f) == 2)
        x33 = x1[1]
        x44 = x1[2]
        new_grad_f[1] = 2.0*x33
        new_grad_f[2] = 2.0*x44
    elseif (2,2) == (rowid, colid)
        assert(length(new_grad_f) == 2)
        x55 = x1[1]
        x66 = x1[2]
        new_grad_f[1] = 2.0*x55
        new_grad_f[2] = 2.0*x66
    elseif (1,0) == (rowid, colid)
        assert(length(new_grad_f) == 1)
        new_grad_f[1] = 0.0
    elseif (2,0) == (rowid, colid)
        assert(length(new_grad_f) == 1)
        new_grad_f[1] = 0.0
    else
        assert(false)
    end
    @show new_grad_f
end


function str_eval_jac_g(rowid,colid,x0,x1,mode,e_rowidx,e_colptr,e_values,i_rowidx,i_colptr,i_values)
    println("str_eval_jac_g")    
    @show rowid,colid
    @show x0
    @show x1
    @show mode
    if(mode == :Structure)
        if (rowid,colid) == (0,0)
            return (0,0)
        elseif (rowid, colid) == (1,1)
            return (3,0)
        elseif (rowid, colid) == (2,2)
            return (3,0)
        elseif (rowid, colid) == (1,0)
            return (1,0)
        elseif (rowid, colid) == (2,0)
            return (1,0)
        else
            assert(false)
        end
    else
        x11 = x0[1]
        if (rowid,colid) == (0,0)   
        elseif (rowid, colid) == (1,1)
            assert(length(e_rowidx) == length(e_values) == 3)
            assert(length(i_rowidx) == length(i_values) == 0)
            e_rowidx[1] = 1
            e_rowidx[2] = 2
        	e_rowidx[3] = 1
            e_colptr[1] = 1
            e_colptr[2] = 3
            e_colptr[3] = 4
            e_values[1] = 1.0
            e_values[2] = -1.0
        	e_values[3] = 1.0
        elseif (rowid, colid) == (2,2)
            assert(length(e_rowidx) == length(e_values) == 3)
            assert(length(i_rowidx) == length(i_values) == 0)
            e_rowidx[1] = 1
            e_rowidx[2] = 2
        	e_rowidx[3] = 1
            e_colptr[1] = 1
            e_colptr[2] = 3
            e_colptr[3] = 4
            e_values[1] = 1.0
            e_values[2] = -1.0
        	e_values[3] = 1.0
        elseif (rowid, colid) == (1,0)
            assert(length(e_rowidx) == length(e_values) == 1)
            assert(length(i_rowidx) == length(i_values) == 0)
            e_rowidx[1] = 2
            e_colptr[1] = 1
            e_colptr[2] = 2
            e_values[1] = 1.0
        elseif (rowid, colid) == (2,0)
            assert(length(e_rowidx) == length(e_values) == 1)
            assert(length(i_rowidx) == length(i_values) == 0)
            e_rowidx[1] = 2
            e_colptr[1] = 1
            e_colptr[2] = 2
            e_values[1] = 1.0
        else
            assert(false)
        end
        convert_to_c_idx(i_rowidx)
        convert_to_c_idx(i_colptr)
        convert_to_c_idx(e_rowidx)
        convert_to_c_idx(e_colptr)
	    @show e_rowidx
	    @show e_colptr
	    @show e_values
	    @show i_rowidx
	    @show i_colptr
	    @show i_values
    end
end

function str_eval_h(rowid,colid,x0,x1,obj_factor,lambda,mode,rowidx,colptr,values)
    println("str_eval_h")
    @show rowid,colid
    @show x0
    @show x1
    @show mode   
    @show obj_factor
    if(mode == :Structure)
        if(rowid,colid) == (0,0)
            return 1
        elseif (rowid,colid) == (1,1)
            return 2
        elseif (rowid,colid) == (2,2)
            return 2
        elseif (rowid,colid) == (0,1)
            return 1
        elseif (rowid,colid) == (0,2)
            return 1
        elseif (rowid,colid) == (1,0)
            return 0
        elseif (rowid,colid) == (2,0)
            return 0
        else
            assert(false)
        end
    else
        x11 = x0[1]
        # @show (rowid,colid)
        if(rowid,colid) == (0,0)
            # @show "diag 0,0 "
            rowidx[1] = 1
            colptr[1] = 1
            colptr[2] = 2
            values[1] = 2.0 * obj_factor
        elseif (rowid,colid) == (1,1)
            # @show "1,1"
            rowidx[1] = 1
            rowidx[2] = 2
            colptr[1] = 1
            colptr[2] = 2
            colptr[3] = 3
            values[1] = 2.0 * obj_factor
            values[2] = 2.0 * obj_factor
        elseif (rowid,colid) == (2,2)
            # @show "2,2"
            rowidx[1] = 1
            rowidx[2] = 2
            colptr[1] = 1
            colptr[2] = 2
            colptr[3] = 3
            values[1] = 2.0 * obj_factor
            values[2] = 2.0 * obj_factor
        elseif (rowid,colid) == (0,1) #the diagnal contribution to 1st stage
            # @show "diag  contr - 1, 0 "
            rowidx[1] = 1
            colptr[1] = 1
            colptr[2] = 2
            values[1] = 0.0
        elseif (rowid,colid) == (0,2) #the diagnal contribution  to 1st stage
            # @show "diag contr - 2, 0 "
            rowidx[1] = 1
            colptr[1] = 1
            colptr[2] = 2
            values[1] = 0.0
        elseif (rowid,colid) == (1,0) #liking border
            # @show "linking border - 0 ,1 "
            # colptr[1] = 1
            # colptr[2] = 2
            # colptr[3] = 3
            # rowidx[1] = 1
            # rowidx[2] = 1
            # values[1] = 1.0 * obj_factor + 1.0 * lambda[1]
            # values[2] = 1.0 * obj_factor
        elseif (rowid,colid) == (2,0) #liking border
            # @show "linking border - ", (rowid,colid)
            # colptr[1] = 1
            # colptr[2] = 2
            # colptr[3] = 3
            # rowidx[1] = 1
            # rowidx[2] = 1
            # values[1] = 1.0 * obj_factor + 1.0 * lambda[1]
            # values[2] = 1.0 * obj_factor
        else
            assert(false)
        end
        convert_to_c_idx(colptr)
        convert_to_c_idx(rowidx)
        @show colptr
        @show rowidx
        @show values
    end
end

function convert_to_c_idx(indicies)
    for i in 1:length(indicies)
        indicies[i] = indicies[i] - 1
    end
end

function create()
    comm = MPI.COMM_WORLD
    println("[$(MPI.Comm_rank(comm))/$(MPI.Comm_size(comm))] create problem ")
    # @show comm
    model = FakeModel(:Min,0, 2,
        str_init_x0, str_prob_info, str_eval_f, str_eval_g, str_eval_grad_f, str_eval_jac_g, str_eval_h)
    prob = createProblemStruct(comm, model) 
    # println("end create problem ")
    return prob
end

function solve(prob)
    # println("Solve problem")
    # @show prob
    ret = solveProblemStruct(prob)
    # println("end solve problem")
end

function main()
    MPI.Init()

    prob = create()
    


    solve(prob)


    MPI.Finalize()
end

main()
