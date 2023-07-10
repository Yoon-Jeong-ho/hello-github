#pragma once

#include "HigherOrderBeamLib.h"
#include "Inertia.h"

tuple<vector<Matrix<double50, -1, -1>>, Matrix<double50, -1, -1>, Matrix<double50, -1, -1>, vector<int>> probtune(vector<Matrix<double50, -1, -1>> coef, Matrix<double50, -1, -1> Q, Array<double50, -1, 1> Le, int t)
{
    Index Nf = coef[0].cols(), NC = coef[0].rows();
    int NE = coef.size();
    vector<Matrix<double50, -1, -1>> Coef = coef;
    vector<int> z, c2;
    ColPivHouseholderQR<Matrix<double50, -1, -1>> qr(Q);
    qr.compute(Q);
    Matrix<double50, -1, -1> P = qr.colsPermutation(), Qi, H, H1, Pi, S, d;
    Tensor<double50, 3> coef2;
    ArrayXi c, C, h;

    c.setZero(P.cols());
    for (int i = 0;i < P.cols(); i++)
    {
        for (int j = 0; j < P.rows(); j++)
        {
            if (P(j, i) == 1)
                c(i) = j;
        }
    }
    if (t != 0)
    {
        for (int i = t; i < t * 2;i++)
        {
            for (int j = 0; j < c.size();j++)
            {
                if (i == c(j))
                {
                    c2.push_back(j);
                }
            }
        }
    }
    for (int i = 0;i < NE;i++)
    {
        coef[i] = Coef[i](all, c);
    }
    
    Pi = Inertia(coef, coef, Le);

    Matrix<double50, -1, -1> evec;
    Array<double50, -1, 1> eval;
    if (Q.rows() == 0)
    {
        EigenSolver<Matrix<double50, -1, -1>> eig1;
        eig1.compute(Pi);
        eval = eig1.eigenvalues().real();
        evec = eig1.eigenvectors().real();
    }
    else
    {
        Matrix<double50, -1, -1> QR = qr.matrixR();
        Matrix<double50, -1, -1> R = QR.triangularView<Upper>();
        Matrix<double50, -1, -1> R1 = R.block(0, 0, R.rows(), R.rows());
        Matrix<double50, -1, -1> R2 = R.block(0, R.rows(), R.rows(), R.cols() - R.rows());
        Matrix<double50, -1, -1> Z, Z1, Z2;
        Z1 = -R1.inverse() * R2;
        Z2.setIdentity(R2.cols(), R2.cols());
        Z.setZero(Z1.rows() + Z2.rows(), Z1.cols());
        Z.block(0, 0, Z1.rows(), Z1.cols()) = Z1;
        Z.block(Z1.rows(), 0, Z2.rows(), Z2.cols()) = Z2;

        Matrix<double50, -1, -1> M1 = Z.transpose() * Pi * Z;
        Matrix<double50, -1, -1> M2 = Z.transpose() * Z;

        GeneralizedEigenSolver<Matrix<double50, -1, -1>> eig1;
        eig1.compute(M1, M2);
        eval = eig1.eigenvalues().real();
        Matrix<double50, -1, -1> x1 = eig1.eigenvectors().real();

        evec = Z * x1;
    }
    vector<int> rng1, rng2, rng, redundant, rnd1;

    for (int i = 0;i < eval.size();i++)
    {
        if (isfinite(eval(i)) && eval(i) > 0)
            rng1.push_back(i);
    }
    for (int i = 0;i < eval.size();i++)
    {
        if (eval.abs()(i) > eval.abs()(rng1).maxCoeff() * 1e-15)
            rng2.push_back(i);
    }
    for (int i = 0;i < rng1.size();i++)
    {
        for (int j = 0;j < rng2.size();j++)
        {
            if (rng1[i] == rng2[j])
                rng.push_back(rng2[j]);
        }
    }

    for (int i = 0;i < rng.size();i++)
    {
        for (int j = 0;j < redundant.size();j++)
        {
            if (i == redundant[j])
                rnd1.push_back(j);
        }

        if (rnd1.size() == 0 && i < rng.size())
        {
            double50 di = eval(rng[i]);
            for (int j = i + 1;j < rng.size(); j++)
            {
                double50 dj = eval(rng[j]);
                if (abs((di - dj)) < (abs(di) * 1e-20))
                    redundant.push_back(j);
            }
        }
    }

    vector<int> effrngv;
    VectorXi effective, rngv, redundv;
    rngv.setZero(rng.size());
    redundv.setZero(redundant.size());
    effective.setZero(rng.size());

    for (int j = 0; j < rng.size();j++)
        rngv(j) = rng[j];

    for (int j = 0; j < redundant.size();j++)
        redundv(j) = redundant[j];

    sort(rngv.begin(), rngv.end());
    sort(redundv.begin(), redundv.end());
    VectorXi::iterator it;

    it = set_difference(rngv.begin(), rngv.end(), redundv.begin(), redundv.end(), effective.begin());
    effective.conservativeResize(it - effective.begin());

    for (int i = 0; i < effective.size(); i++)
        for (int j = 0; j < rngv.size(); j++)
        {
            if (effective(i) == rngv[j])
                effrngv.push_back(j);
        }

    Matrix<double50, -1, -1> evec1 = evec(all, rngv(effrngv));
    Array<double50, -1, 1> eval1 = eval(rngv(effrngv));

    return make_tuple(coef, Pi, evec1, c2);
}

Matrix<double50, -1, -1> eig(Matrix<double50, -1, -1> P, Matrix<double50, -1, -1> S)
{
    Matrix<double50, -1, -1> evec;
    Array<double50, -1, 1> eval;
    EigenSolver<Matrix<double50, -1, -1>> eig1;
    eig1.compute(P);
    eval = eig1.eigenvalues().real();
    evec = eig1.eigenvectors().real();

    vector<int> rng1, rng2, rng, redundant, rnd1;

    for (int i = 0;i < eval.size();i++)
    {
        if (isfinite(eval(i)) && eval(i) > 0)
            rng1.push_back(i);
    }
    for (int i = 0;i < eval.size();i++)
    {
        if (eval.abs()(i) > eval.abs()(rng1).maxCoeff() * 1e-10)
            rng2.push_back(i);
    }
    for (int i = 0;i < rng1.size();i++)
    {
        for (int j = 0;j < rng2.size();j++)
        {
            if (rng1[i] == rng2[j])
                rng.push_back(rng2[j]);
        }
    }

    for (int i = 0;i < rng.size();i++)
    {
        for (int j = 0;j < redundant.size();j++)
        {
            if (i == redundant[j])
                rnd1.push_back(j);
        }

        if (rnd1.size() == 0 && i < rng.size())
        {
            double50 di = eval(rng[i]);
            for (int j = i + 1;j < rng.size(); j++)
            {
                double50 dj = eval(rng[j]);
                if (abs((di - dj)) < (abs(di) * 1e-20))
                    redundant.push_back(j);
            }
        }
    }

    vector<int> effrngv;
    VectorXi effective, rngv, redundv;
    rngv.setZero(rng.size());
    redundv.setZero(redundant.size());
    effective.setZero(rng.size());

    for (int j = 0; j < rng.size();j++)
        rngv(j) = rng[j];

    for (int j = 0; j < redundant.size();j++)
        redundv(j) = redundant[j];

    sort(rngv.begin(), rngv.end());
    sort(redundv.begin(), redundv.end());
    VectorXi::iterator it;

    it = set_difference(rngv.begin(), rngv.end(), redundv.begin(), redundv.end(), effective.begin());
    effective.conservativeResize(it - effective.begin());

    for (int i = 0; i < effective.size(); i++)
        for (int j = 0; j < rngv.size(); j++)
        {
            if (effective(i) == rngv[j])
                effrngv.push_back(j);
        }

    Matrix<double50, -1, -1> evec1 = evec(all, rngv(effrngv));
    Array<double50, -1, 1> eval1 = eval(rngv(effrngv));

    return evec1;
}