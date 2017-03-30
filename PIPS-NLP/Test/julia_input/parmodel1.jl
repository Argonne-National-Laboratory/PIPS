import MPI

# min (x1+x2)^2+(x1+x2)*x3 + (x1+x2)*x4
# st.
#     x1 * x2 = 10
#     x2^2 + x3*x1 < 5
#     x2^2 + x4*x1 < 6
# x1, x2 , x3, x4 free variables

#include("../../JuliaInterface/ParPipsNlp.jl")
include("/Users/frank/.julia/v0.4/StructJuMPSolverInterface/src/pips_parallel_cfunc.jl")

using ParPipsNlp

function str_init_x0(nodeid, x0)
    if(nodeid == 0)
        assert(length(x0) ==2)
        x0[1] = 1.0
        x0[2] = 1.0
    elseif nodeid == 1
        assert(length(x0) ==1)
        x0[1] = 1.0
    elseif nodeid == 2
        assert(length(x0) ==1)
        x0[1] = 1.0
    end
end

function str_prob_info(nodeid,mode,col_lb,col_ub,row_lb,row_ub)

    if(mode == :Values)
        #setting col lb ub, row lb ub
        if(nodeid==0)
            fill!(col_lb,-Inf)
            fill!(col_ub, Inf)
            fill!(row_lb, 10.0)
            fill!(row_ub, 10.0)
        elseif (nodeid == 1 )
            fill!(col_lb, -Inf)
            fill!(col_ub, Inf)
            fill!(row_lb, -Inf)
            fill!(row_ub, 5)
        elseif (nodeid == 2)
            fill!(col_lb, -Inf)
            fill!(col_ub, Inf)
            fill!(row_lb, -Inf)
            fill!(row_ub, 6)
        else
            assert(false)
        end  

        @show col_lb
        @show col_ub
        @show row_lb
        @show row_ub  
    end

    if(nodeid==0)
        return (2,1)
    elseif(nodeid == 1 || nodeid == 2)
        return (1,1)
    else
        assert(false)
    end
end

function str_eval_f(nodeid,x0,x1)
    @show nodeid
    @show x0
    @show x1
    fval = 0.0
    if nodeid ==  0
        x01 = x0[1]
        x02 = x0[2]
        fval = (x01+x02)^2
    elseif nodeid == 1
        x01 = x0[1]
        x02 = x0[2]
        x13 = x1[1]
        fval = (x01+x02)*x13
    elseif nodeid == 2
        x01 = x0[1]
        x02 = x0[2]
        x24 = x1[1]
        fval = (x01+x02)*x24
    else
        assert(false)
    end
    return fval
end

function str_eval_g(nodeid,x0,x1,new_eq_g, new_inq_g)
    if nodeid == 0
        assert(x0 == x1)
        assert(length(new_inq_g) == 0 )
        assert(length(new_eq_g) == 1)
        x01 = x0[1]
        x02 = x0[2]
        new_eq_g[1] = x01 * x02
    elseif nodeid == 1
        x01 = x0[1]
        x02 = x0[2]
        x13 = x1[1]
        assert(length(new_inq_g) == 1 )
        assert(length(new_eq_g) == 0)
        new_inq_g[1] = x02^2 + x13*x01
    elseif nodeid == 2
        x01 = x0[1]
        x02 = x0[2]
        x24 = x1[1]
        assert(length(new_inq_g) == 1 )
        assert(length(new_eq_g) == 0)
        new_inq_g[1] = x02^2 + x24*x01
    else
        assert(false)
    end
end

function str_eval_grad_f(rowid,colid,x0,x1,new_grad_f)
    @show rowid
    @show colid
    @show x0
    @show x1
    if (0,0) == (rowid, colid)
        assert(length(new_grad_f) == 2)
        assert(x1 == x0)
        x01 = x0[1]
        x02 = x0[2]
        new_grad_f[1] = 2*(x01 + x02)
        new_grad_f[2] = 2*(x01 + x02)
    elseif (1,0) == (rowid, colid)
        assert(length(new_grad_f) == 2)
        x01 = x0[1]
        x02 = x0[2]
        x13 = x1[1]
        new_grad_f[1] = x13
        new_grad_f[2] = x13
    elseif (1,1) == (rowid, colid)
        assert(length(new_grad_f) == 1)
        x01 = x0[1]
        x02 = x0[2]
        x13 = x1[1]
        new_grad_f[1] = x01+x02
    elseif (2,0) == (rowid, colid)
        assert(length(new_grad_f) == 2)
        x01 = x0[1]
        x02 = x0[2]
        x24 = x1[1]
        new_grad_f[1] = x24
        new_grad_f[2] = x24
    elseif (2,2) == (rowid, colid)
        assert(length(new_grad_f) == 1)
        x01 = x0[1]
        x02 = x0[2]
        x24 = x1[1]
        new_grad_f[1] = x01+x02
    else
        assert(false)
    end
end

function str_eval_jac_g(rowid,colid,x0,x1,mode,e_rowidx,e_colptr,e_values,i_rowidx,i_colptr,i_values)
    if(mode == :Structure)
        if (rowid,colid) == (0,0)
            return (2,0)
        elseif (rowid, colid) == (1,1)
            return (0,1)
        elseif (rowid, colid) == (2,2)
            return (0,1)
        elseif (rowid, colid) == (1,0)
            return (0,2)
        elseif (rowid, colid) == (2,0)
            return (0,2)
        else
            assert(false)
        end
    else
        x01 = x0[1]
        x02 = x0[2]
        if (rowid,colid) == (0,0)
            assert(x1==x0)
            assert(length(e_rowidx) == length(e_values) == 2)
            assert(length(i_rowidx) == length(i_values) == 0)
            e_rowidx[1] = 1
            e_rowidx[2] = 1
            e_colptr[1] = 1
            e_colptr[2] = 2
            e_colptr[3] = 3
            e_values[1] = x02
            e_values[2] = x01
        elseif (rowid, colid) == (1,1)
            assert(length(e_rowidx) == length(e_values) == 0)
            assert(length(i_rowidx) == length(i_values) == 1)
            x13 = x1[1]
            i_rowidx[1] = 1
            i_colptr[1] = 1
            i_colptr[2] = 2
            i_values[1] = x01
        elseif (rowid, colid) == (2,2)
            assert(length(e_rowidx) == length(e_values) == 0)
            assert(length(i_rowidx) == length(i_values) == 1)
            x24 = x1[1]
            i_rowidx[1] = 1
            i_colptr[1] = 1
            i_colptr[2] = 2
            i_values[1] = x01
        elseif (rowid, colid) == (1,0)
            assert(length(e_rowidx) == length(e_values) == 0)
            assert(length(i_rowidx) == length(i_values) == 2)
            x13 = x1[1]
            i_rowidx[1] = 1
            i_rowidx[2] = 1
            i_colptr[1] = 1
            i_colptr[2] = 2
            i_colptr[3] = 3
            i_values[1] = x13
            i_values[2] = 2.0*x02
        elseif (rowid, colid) == (2,0)
            assert(length(e_rowidx) == length(e_values) == 0)
            assert(length(i_rowidx) == length(i_values) == 2)
            x24 = x1[1]
            i_rowidx[1] = 1
            i_rowidx[2] = 1
            i_colptr[1] = 1
            i_colptr[2] = 2
            i_colptr[3] = 3
            i_values[1] = x24
            i_values[2] = 2.0*x02
        else
            assert(false)
        end
        convert_to_c_idx(i_rowidx)
        convert_to_c_idx(i_colptr)
        convert_to_c_idx(e_rowidx)
        convert_to_c_idx(e_colptr)
    end
end

function str_eval_h(rowid,colid,x0,x1,obj_factor,lambda,mode,rowidx,colptr,values)
    if(mode == :Structure)
        if(rowid,colid) == (0,0)
            return 3
        elseif (rowid,colid) == (1,1)
            return 0
        elseif (rowid,colid) == (2,2)
            return 0
        elseif (rowid,colid) == (1,0)
            return 3
        elseif (rowid,colid) == (2,0)
            return 3
        elseif (rowid,colid) == (0,1)
            return 2
        elseif (rowid,colid) == (0,2)
            return 2
        else
            assert(false)
        end
    else
        assert(length(lambda) == 1)
        @show (rowid,colid)
        @show colptr
        @show rowidx
        @show values
        if(rowid,colid) == (0,0)
            @show "diag 0,0 "
            rowidx[1] = 1
            rowidx[2] = 2
            rowidx[3] = 2
            colptr[1] = 1
            colptr[2] = 3
            colptr[3] = 4
            values[1] = 2.0 * obj_factor
            values[2] = 2.0 * obj_factor + 1.0 * lambda[1]
            values[3] = 2.0 * obj_factor
        elseif (rowid,colid) == (1,1)
            @show "1,1"
        elseif (rowid,colid) == (2,2)
            @show "2,2"
        elseif (rowid,colid) == (1,0) #the diagnal contribution to 1st stage
            @show "diag  contr - 1, 0 "
            rowidx[1] = 1
            rowidx[2] = 2
            rowidx[3] = 2
            colptr[1] = 1
            colptr[2] = 3
            colptr[3] = 4
            values[1] = 0.0
            values[2] = 0.0
            values[3] = 2.0 * lambda[1]
        elseif (rowid,colid) == (2,0) #the diagnal contribution  to 1st stage
            @show "diag contr - 2, 0 "
            rowidx[1] = 1
            rowidx[2] = 2
            rowidx[3] = 2
            colptr[1] = 1
            colptr[2] = 3
            colptr[3] = 4
            values[1] = 0.0
            values[2] = 0.0
            values[3] = 2.0 * lambda[1]
        elseif (rowid,colid) == (0,1) #liking border
            @show "linking border - 0 ,1 "
            colptr[1] = 1
            colptr[2] = 2
            colptr[3] = 3
            rowidx[1] = 1
            rowidx[2] = 1
            values[1] = 1.0 * obj_factor + 1.0 * lambda[1]
            values[2] = 1.0 * obj_factor  
        elseif (rowid,colid) == (0,2) #liking border
            @show "linking border - ", (rowid,colid)
            colptr[1] = 1
            colptr[2] = 2
            colptr[3] = 3
            rowidx[1] = 1
            rowidx[2] = 1
            values[1] = 1.0 * obj_factor + 1.0 * lambda[1]
            values[2] = 1.0 * obj_factor  
        else
            assert(false)
        end
        convert_to_c_idx(colptr)
        convert_to_c_idx(rowidx)
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
    @show comm
    model = FakeModel(:Min,0, 2,
        str_init_x0, str_prob_info, str_eval_f, str_eval_g, str_eval_grad_f, str_eval_jac_g, str_eval_h)
    prob = createProblemStruct(comm, model) 
    println("end create problem ")
    return prob
end

function solve(prob)
    println("Solve problem")
    @show prob
    ret = solveProblemStruct(prob)
    println("end solve problem")
end

function main()
    MPI.Init()

    prob = create()
    

    solve(prob)


    MPI.Finalize()
end

#main()
