import MPI

# min x1^2 + x2^2 +x1x2 + x3^2 + x4^2 +x3x4 + x5^2 + x6^2 +x5x6
# st.
#     x1 + x2 = 100
#     0< x1   + x3 + x4                                 < 500
#     0<   x2           + x5 + x6                       < 500
#     0< x1                       + x7 + x8             < 500
#     0<   .                                ...         < 500
#     0<   .                                      ...   < 500
# x free variables

include("./ParPipsNlp.jl")

const scen = 100
const nnodes = scen + 1

using ParPipsNlp

function str_init_x0(nodeid, x0)
    assert(length(x0) ==2)
    assert(nodeid<nnodes)
    x0[1] = 1.0
    x0[2] = 1.0
end

function str_prob_info(nodeid,mode,col_lb,col_ub,row_lb,row_ub)
    assert(nodeid<nnodes)
    if(mode == :Values)
        #setting col lb ub, row lb ub
        if(nodeid==0)
            assert(length(row_lb) == length(row_ub))
            assert(length(col_lb) == length(col_ub))
            fill!(col_lb,-Inf)
            fill!(col_ub, Inf)
            fill!(row_lb, 100.0)
            fill!(row_ub, 100.0)
        else
            assert(length(row_lb) == length(row_ub))
            assert(length(col_lb) == length(col_ub))
            fill!(col_lb, -Inf)
            fill!(col_ub, Inf)
            fill!(row_lb, 0)
            fill!(row_ub, 500)
        end  
    end
    # @show col_lb
    # @show col_ub
    # @show row_lb
    # @show row_ub  
    return (2,1)
end

function str_eval_f(nodeid,x0,x1)
    # @show nodeid
    # @show x0
    # @show x1
    assert(nodeid<nnodes)
    fval = 0.0
    if nodeid ==  0
        assert(x0==x1)
        fval = x1[1]*x1[2]
    else
        fval = x1[1]*x1[2]
    end
   
    return fval
end

function str_eval_g(nodeid,x0,x1,new_eq_g, new_inq_g)
    # @show nodeid
    assert(nodeid<nnodes)
    x11 = x0[1]
    x22 = x0[2] 
    x33 = x1[1]
    x44 = x1[2]   
    if nodeid == 0
        assert(x0 == x1)
        assert(length(new_inq_g) == 0 )
        assert(length(new_eq_g) == 1)
        new_eq_g[1] = x11 + x22
    elseif nodeid%2==0
        assert(length(new_inq_g) == 1 )
        assert(length(new_eq_g) == 0)
        new_inq_g[1] = x22 + x33 + x44
    else #if nodeid%2 == 1
        assert(length(new_inq_g) == 1 )
        assert(length(new_eq_g) == 0)
        new_inq_g[1] = x11 + x33 + x44
    end
end

function str_eval_grad_f(rowid,colid,x0,x1,new_grad_f)
    # @show rowid,colid
    # @show x0
    # @show x1
    assert(rowid<nnodes && colid<nnodes)  
    x11 = x0[1]
    x22 = x0[2]
    x33 = x1[1]
    x44 = x1[2] 
    if rowid == colid
        assert(length(new_grad_f) == 2)
        if(rowid == 0) assert(x1 == x0) end
        new_grad_f[1] = 2.0*x33
        new_grad_f[2] = 2.0*x44
    else 
    	assert(colid==0&&rowid!=0)
        assert(length(new_grad_f) == 2)
        new_grad_f[1] = 0.0
        new_grad_f[2] = 0.0
    end
end

function str_eval_jac_g(rowid,colid,x0,x1,mode,e_rowidx,e_colptr,e_values,i_rowidx,i_colptr,i_values)
    assert(rowid<nnodes && colid<nnodes)
    if(mode == :Structure)
        if (rowid,colid) == (0,0)
            return (2,0)
        elseif (rowid == colid) 
            return (0,2)
        else
            return (0,1)
        end
    else
        x11 = x0[1]
        x22 = x0[2]
        x33 = x1[1]
        x44 = x1[2]
        if (rowid,colid) == (0,0)
            assert(length(e_rowidx) == length(e_values) == 2)
            assert(length(i_rowidx) == length(i_values) == 0)
            e_rowidx[1] = 1
            e_rowidx[2] = 1
            e_colptr[1] = 1
            e_colptr[2] = 2
            e_colptr[3] = 3
            e_values[1] = 1.0
            e_values[2] = 1.0
        elseif rowid == colid
            assert(length(e_rowidx) == length(e_values) == 0)
            assert(length(i_rowidx) == length(i_values) == 2)
            i_rowidx[1] = 1
            i_rowidx[2] = 1
            i_colptr[1] = 1
            i_colptr[2] = 2
            i_colptr[3] = 3
            i_values[1] = 1.0
            i_values[2] = 1.0
        elseif colid%2==0
            assert(length(e_rowidx) == length(e_values) == 0)
            assert(length(i_rowidx) == length(i_values) == 1)
            i_rowidx[1] = 1
            i_colptr[1] = 1
            i_colptr[2] = 2
            i_colptr[3] = 2
            i_values[1] = 1.0
        elseif colid%2==1
            assert(length(e_rowidx) == length(e_values) == 0)
            assert(length(i_rowidx) == length(i_values) == 1)
            i_rowidx[1] = 1
            i_colptr[1] = 1
            i_colptr[2] = 1
            i_colptr[3] = 2
            i_values[1] = 1.0
        end
        convert_to_c_idx(i_rowidx)
        convert_to_c_idx(i_colptr)
        convert_to_c_idx(e_rowidx)
        convert_to_c_idx(e_colptr)
    end
end

function str_eval_h(rowid,colid,x0,x1,obj_factor,lambda,mode,rowidx,colptr,values)
    assert(rowid<nnodes && colid<nnodes)
    # @show obj_factor
    if(mode == :Structure)
        if rowid == colid   #diagonal block
            return 2
        else
            if(colid == 0)  #root diagonal contribution
                return 2
            else
                return 0   #linker border
            end
        end
    else
        assert(length(lambda) == 1)
        x11 = x0[1]
        x22 = x0[2]
        # @show (rowid,colid)
        if rowid==colid			#diagonal block
            # @show "diag 0,0 "
            colptr[1] = 1
            colptr[2] = 2
            colptr[3] = 3
            rowidx[1] = 1
            rowidx[2] = 2
            values[1] = 2.0 * obj_factor
            values[2] = 2.0 * obj_factor
        else
            if(colid == 0) 		#root diagonal contribution
                colptr[1] = 1
                colptr[2] = 2
                colptr[3] = 3
                rowidx[1] = 1
                rowidx[2] = 2
                values[1] = 0.0
                values[2] = 0.0
            end
        end    
        
        convert_to_c_idx(colptr)
        convert_to_c_idx(rowidx)
        # @show colptr
        # @show rowidx
        # @show values
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
    prob = createProblemStruct(comm,
        scen,  #number scen
        str_init_x0, str_prob_info, str_eval_f, str_eval_g, str_eval_grad_f,
        str_eval_jac_g, str_eval_h) 
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
