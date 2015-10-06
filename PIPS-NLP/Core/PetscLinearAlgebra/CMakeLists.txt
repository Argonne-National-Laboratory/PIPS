
if(HAVE_PETSC)
  add_library(petsc_solver
    PetscIterativeSolver_Schur.C
  )
  add_library(petsc_solver_sparse
    PetscIterativeSolver_Sparse.C
  )
endif(HAVE_PETSC)


