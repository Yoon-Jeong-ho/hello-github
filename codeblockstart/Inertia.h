#pragma once

#include "HigherOrderBeamLib.h"
#include "Calculus.h"

Matrix<double50, -1, -1> shift(Matrix<double50, -1, -1> m1, int n)
{
    int row1 = m1.rows(), col1 = m1.cols();
    ArrayXi r1;
    r1.setLinSpaced(static_cast<Eigen::Index>(row1) - n, 0, row1 - 1 - n);
    Matrix<double50, -1, -1> m2;
    m2.setZero(row1, col1);
    m2.block(n, 0, row1 - n, col1) = m1.block(0, 0, row1 - n, col1);

    return m2;
}

Matrix<double50, -1, -1> Inertia(vector<Matrix<double50, -1, -1>> coef1, vector<Matrix<double50, -1, -1>> coef2, Array<double50, -1, 1> Le)
{
    int NE = Le.rows();
    Index NC = coef1[0].rows() * 2, NM1 = coef1[0].cols(), NM2 = coef2[0].cols();

    Matrix<double50, -1, -1> c2;
    vector<Matrix<double50, -1, -1>> c3, c3i; 
    Matrix<double50, -1, -1> J(NM1, NM2), val(1, NM2), s;
    Vector<double50, -1> p, c1; 
    double50 c1k, Li;
    Calculus integ;

    c1.setZero(NC, 1), c2.setZero(NC, NM2), c3.resize(1), c3[0].setZero(NC, NM2);
    int n{};
    J.setZero();

    p.setLinSpaced(NC, 0, NC - 1);

    for (int i = 0;i < NM1;i++)
    {
        val.setZero();
        for (int e = 0;e < NE;e++)
        {
            Li = Le(e);
            s.setZero(1, NC);
            c1.block(0, 0, NC / 2, 1) = coef1[e].block(0, i, NC / 2, 1);
            c2.block(0, 0, NC / 2, NM2) = coef2[e];
            for (int j = 0;j < NC;j++)
            {
                s(j) = pow(Li, p(j));
                if (c1(j) != 0)
                {
                    n = j + 1;
                }
            }
            c3[0].setZero();
            if (n != 0)
            {
                for (int k = 1;k <= n;k++)
                {
                    c1k = c1(k - 1);
                    c3[0] += c1k * shift(c2, k - 1);
                }
                c3i = integ.Intcoef(c3);
                val += s * c3i[0];
            }
        }
        J.block(i, 0, 1, NM2) = val;
    }

    return J;
}
