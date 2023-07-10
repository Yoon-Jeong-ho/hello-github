#pragma once

#include "HigherOrderBeamLib.h"
#include "slepceps.h"
#include <chrono>

using namespace std::chrono;

static char help[] = "Standard symmetric eigenproblem corresponding to the Laplacian operator in 1 dimension.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";

tuple<VectorXd, vector<VectorXd>, int, int> EigSolver(SparseMatrix<double> KS, SparseMatrix<double> MS, int freqn, int tdof, int nnz)
{
    Mat         A, B;
    EPS         eps;
    ST          st;
    KSP         ksp;
    PC          pc;
    PetscReal   error1, tol = 1e-30, re, im;
    PetscScalar kr, ki;
    Vec         xr, xi;
    VectorXd    eval, freq;
    vector<VectorXd>    evec;
    PetscInt n = KS.rows(), nev = freqn, maxit = 1e6, nconv = nev;
    //SlepcInitialize(&argc,&argv,(char*)0,help);
    SlepcInitialize((int*)0, (char***)0, (char*)0, help);

    MatCreateSeqAIJ(MPI_COMM_WORLD,n,n,nnz,0,&A);
    //PetscCall(MatCreateAIJ(MPI_COMM_WORLD,n,n,PETSC_DETERMINE,PETSC_DETERMINE,1000,PETSC_NULL,1000,PETSC_NULL,&A));
    MatSetFromOptions(A);
    MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
    MatSetUp(A);
    MatCreateSeqAIJ(MPI_COMM_WORLD,n,n,nnz,0,&B);
    //PetscCall(MatCreateAIJ(MPI_COMM_WORLD,n,n,PETSC_DETERMINE,PETSC_DETERMINE,1000,PETSC_NULL,1000,PETSC_NULL,&B));
    MatSetFromOptions(B);
    MatSetOption(B,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
    MatSetUp(B);

    for (int i = 0; i < KS.outerSize(); ++i)
        for (SparseMatrix<double>::InnerIterator it(KS, i); it; ++it)
        {
            MatSetValue(A, it.row(), it.col(), it.value(), INSERT_VALUES);
        }

    for (int i = 0; i < MS.outerSize(); ++i)
        for (SparseMatrix<double>::InnerIterator jt(MS, i); jt; ++jt)
        {
            MatSetValue(B, jt.row(), jt.col(), jt.value(), INSERT_VALUES);
        }

    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);

    MatCreateVecs(A,NULL,&xr);
    MatCreateVecs(A,NULL,&xi);

    cout << "Creating Eigen Solver" << endl;

    auto start = high_resolution_clock::now();
    EPSCreate(MPI_COMM_WORLD,&eps);

    EPSSetOperators(eps,A,B);
    EPSSetProblemType(eps,EPS_GHEP);

    EPSSetType(eps,EPSKRYLOVSCHUR);
    EPSGetST(eps,&st);
    STSetType(st,STSINVERT);
    STSetShift(st,0);
    STGetKSP(st,&ksp);
    KSPSetType(ksp,KSPPREONLY);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCCHOLESKY);
    PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
    EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE);
    EPSSetTarget(eps,0.0);
    EPSSetTolerances(eps,tol,maxit);
    EPSSetDimensions(eps, nev, nev * 2, tdof / 2);
    EPSSetFromOptions(eps);

    cout << endl << "Solving the Eigensystem:" << endl << endl;

    EPSSolve(eps);

    PetscPrintf(MPI_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev);
    EPSGetTolerances(eps,&tol,&maxit);

    EPSGetConverged(eps,&nconv);

    PetscInt *ei;
    PetscScalar *eveci;
    PetscMalloc1(tdof, &ei);
    PetscMalloc1(tdof, &eveci);

    for (int i = 0; i < tdof; i++)
        ei[i] = i;
    if (nconv>0)
    {
        eval.setZero(nconv);

        for (int i = 0; i < nconv; i++)
        {
            EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);
            EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error1);

            #if defined(PETSC_USE_COMPLEX)
                re = PetscRealPart(kr);
                im = PetscImaginaryPart(kr);
            #else
                re = kr;
                im = ki;
            #endif

            eval(i) = re;
            VecGetValues(xr, tdof, ei, eveci);
            Map<VectorXd> evecj(eveci, tdof);
            evec.push_back(evecj);
        }
        PetscPrintf(MPI_COMM_WORLD,"\n");
    }

    PetscFree(ei);
    PetscFree(eveci);

    EPSDestroy(&eps);
    MatDestroy(&A);
    MatDestroy(&B);
    VecDestroy(&xr);
    VecDestroy(&xi);
    SlepcFinalize();

    freq = eval.cwiseAbs().cwiseSqrt()/(2*M_PI);

    cout << "Frequencies: \n" << freq << endl << endl;

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Modal Frequencies computed. EigenSolver Time:  " << duration.count() << endl << endl;

    if (freq.size() < nev)
        nev = freq.size();

    return make_tuple(freq, evec, nev, nconv);
}
