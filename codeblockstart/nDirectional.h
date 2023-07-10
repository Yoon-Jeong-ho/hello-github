#pragma once

#include "HigherOrderBeamLib.h"
#include "Edge.h"
#include "Continuity.h"

class Equation
{
public:
    Matrix<double50, -1, -1> h1;
    ArrayXi h2;
    Matrix<double50, -1, -1> h(MatrixXi cscc, double50 s)
    {
        if (cscc.rows() == 2)
        {
            h1.setZero(1, 3);
            h1(0, 0) = 1;
            h1(0, 1) = s * s;
            h1(0, 2) = s * s * s;
        }
        else
        {
            h1.setZero(1, 4);
            h1(0, 0) = 1;
            h1(0, 1) = s;
            h1(0, 2) = s * s;
            h1(0, 3) = s * s * s;
        }
        return h1;
    }
    Matrix<double50, -1, -1> dh(MatrixXi cscc, double50 s)
    {
        if (cscc.rows() == 2)
        {
            h1.setZero(1, 3);
            h1(0, 0) = 0;
            h1(0, 1) = 2 * s;
            h1(0, 2) = 3 * s * s;
        }
        else
        {
            h1.setZero(1, 4);
            h1(0, 0) = 0;
            h1(0, 1) = 1;
            h1(0, 2) = 2 * s;
            h1(0, 3) = 3 * s * s;
        }
        return h1;
    }
    Matrix<double50, -1, -1> ddh(MatrixXi cscc, double50 s)
    {
        if (cscc.rows() == 2)
        {
            h1.setZero(1, 3);
            h1(0, 0) = 0;
            h1(0, 1) = 2;
            h1(0, 2) = 6 * s;
        }
        else
        {
            h1.setZero(1, 4);
            h1(0, 0) = 0;
            h1(0, 1) = 0;
            h1(0, 2) = 2;
            h1(0, 3) = 6 * s;
        }
        return h1;
    }
    ArrayXi rng(ArrayXi ne, int e)
    {
        int c = 0;
        for (int i = 0;i <= e; i++)
        {
            c += ne(i);
        }
        h2.setLinSpaced(ne(e), c - ne(e), c - 1);
        return h2;
    }
    ArrayXi hrng(ArrayXi ne, int e)
    {
        return h2.setLinSpaced(ne(e), 0, ne(e) - 1);
    }
};

vector<Matrix<double50, -1, -1>> nDirectional(Matrix<double50, -1, -1> csc, MatrixXi cscc, vector<Matrix<double50, -1, -1>> coef_s, vector<Matrix<double50, -1, -1>> coef_n)
{
    int NE = coef_s.size(), Nc = csc.rows(), ord = coef_s[0].rows() - 1, NC = ord + 1;
    Index NMs = coef_s[0].cols(), NMn = coef_n[0].cols();

    Array<double50, -1, 1> alpha, Le, s;
    ArrayXi mi, nc(Nc, 1), ne(NE, 1);
    Edge edge; Continuity cont;
    alpha = edge.Anglemp(csc, cscc);
    Le = edge.Lengthmp(csc, cscc);
    vector<vector<int>> cont_edg;
    vector<Matrix<double50, -1, -1>> cont_Rxy, coef;
    vector<Array<double50, -1, 1>> cont_sc;
    vector<int> e, nrow = { 1,3,4 };
    Equation hf;

    Matrix<double50, -1, -1> C, V, Uxy, Ci, hval, Vi, c0, c1, c2, vej(2, 1), D;
    coef.resize(NE);

    for (int m = NMn; m < NMs; m++)
    {
        int countUxy = 0;
        mi.setLinSpaced(NE, m, m + (NE - 1) * NMs);
        for (int i = 0; i < NE; i++)
            coef[i] = coef_s[i](all, m);
        tie(cont_edg, cont_Rxy, cont_sc) = cont.Edges(csc, cscc, coef, 1);
        nc.setZero(), ne.setConstant(4);
        for (int i = 0; i < Nc;i++)
        {
            e = cont_edg[i];
            if (e.size() == 1)
            {
                nc(i) = 1;
                ne(e) = 3;
            }
            else
                nc(i) = e.size() * 2;
        }
        C.resize(nc.sum(), nc.sum());
        V.resize(nc.sum(), 1);
        C.setZero();
        V.setZero();
        for (int i = 0; i < Nc; i++)
        {
            e = cont_edg[i];
            s = cont_sc[i];

            if (e.size() == 1)
            {
                Ci.setZero(1, C.rows());
                hval = hf.ddh(cscc, s(0));
                Ci(0, hf.rng(ne, e[0])) = hval(0, hf.hrng(ne, e[0]));
            }
            else
            {
                Uxy = cont_Rxy[countUxy];
                countUxy++;
                Vi.resize(nc(i), 1);
                c0.resize(e.size(), C.rows());
                Vi.setZero();
                c0.setZero();
                for (int j = 0; j < e.size(); j++)
                {
                    hval = hf.h(cscc, s(j));
                    c0(j, hf.rng(ne, e[j])) = hval(0, hf.hrng(ne, e[j]));
                    vej(0, 0) = sin(alpha(e[j]));
                    vej(1, 0) = -cos(alpha(e[j]));
                    Vi(j, 0) = vej.cwiseProduct(Uxy).sum();
                }
                c1.resize(e.size() - 1, C.rows());
                c1.setZero();
                for (int j = 0; j < e.size() - 1;j++)
                {
                    for (int k = 0; k < 2;k++)
                    {
                        hval = hf.dh(cscc, s(static_cast<Eigen::Index>(j) + k));
                        c1(j, hf.rng(ne, e[static_cast<std::vector<int, std::allocator<int>>::size_type>(j) + k])) = hval(0, hf.hrng(ne, e[static_cast<std::vector<int, std::allocator<int>>::size_type>(j) + k])) * pow(-1, k + 1);
                    }
                }
                c2.resize(1, C.rows());
                c2.setZero();
                for (int j = 0; j < e.size();j++)
                {
                    hval = hf.ddh(cscc, s(j));
                    if (s(j) == 0)
                        c2(0, hf.rng(ne, e[j])) = hval(0, hf.hrng(ne, e[j])) * (-1);
                    else
                        c2(0, hf.rng(ne, e[j])) = hval(0, hf.hrng(ne, e[j]));
                }
                Ci.resize(c0.rows() + c1.rows() + c2.rows(), c2.cols());
                Ci.block(0, 0, c0.rows(), c0.cols()) = c0;
                Ci.block(c0.rows(), 0, c1.rows(), c0.cols()) = c1;
                Ci.block(c0.rows() + c1.rows(), 0, c2.rows(), c0.cols()) = c2;

                V(hf.rng(nc, i), 0) = Vi;
            }
            C(hf.rng(nc, i), all) = Ci;
        }
        D = C.inverse() * V;

        if (cscc.rows() == 2)
        {
            for (int j = 0;j < NE; j++)
            {
                coef_n[j].conservativeResize(NC, NMs);
                coef_n[j](nrow, m) = D(hf.rng(ne, j), 0);
            }
        }
        else
        {
            for (int j = 0; j < NE; j++)
            {
                coef_n[j].conservativeResize(NC, NMs);
                coef_n[j].block(0, m, hf.rng(ne, j).size(), 1) = D(hf.rng(ne, j), 0);
            }
        }
    }

    return coef_n;
}