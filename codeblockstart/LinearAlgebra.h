#pragma once

#include "HigherOrderBeamLib.h"

ArrayXd axistrans(MatrixXd v1, Vector<double, 3> org, Vector<double, 3> z, Vector<double, 3> x, int p)
{
    int t = 0;
    Vector<double, 3> y, v2;
    if (v1.cols() == 3)
    {
        v1 = v1.transpose();
        t = 1;
    }

    z = z / z.norm();
    x = x / x.norm();
    y = z.cross(x);
    y = y / y.norm();

    Matrix3d T(3, 3);
    T.setZero();
    T << x.transpose(),
        y.transpose(),
        z.transpose();

    if (p == 0)
        v2 = T.inverse() * v1 + org;
    else
        v2 = T * (v1 - org);

    if (t == 1)
        v2 = v2.transpose();

    return v2;
}

Vector<double, 3> normalize(Vector<double, 3> v1)
{
    Vector<double, 3> v2 = v1 / v1.norm();
    return v2;
}

tuple<double, double> intersection(Array3d u0, Vector3d u1, Array3d v0, Vector3d v1)
{
    double a, b;
    int tol = 8;
    vector<int> zz, nn;
    MatrixXd w1(3, 2);
    w1.setZero();
    w1(all, 0) = u1;
    w1(all, 1) = -v1;
    Vector3d w0 = v0 - u0;
    VectorXd w00, ab;
    MatrixXd ab1;

    Array3d w11 = w1.cwiseAbs().rowwise().sum();
    for (int i = 0; i < w11.size(); i++)
        if (w11(i) == 0)
            zz.push_back(i);

    w00.setZero(zz.size());
    for (int i = 0; i < zz.size(); i++)
        w00(i) = round(w0(zz[i], 0) * 1e8) / 1e8;
    if (zz.size() == 1)
    {
        if (w00(0) == 0)
        {
            for (int i = 0; i < 3; i++)
                if (i != zz[0])
                    nn.push_back(i);
            ab = (w0(nn).transpose() * w1(nn, all).transpose().inverse()).transpose();
            a = ab(0, 0), b = ab(1, 0);
        }
        else
        {
            a = INT_MAX;
            b = INT_MAX;
        }
    }
    else if (zz.size() == 2)
    {
        if (w00(0) == 0 && w00(1) == 0)
        {
            a = 0;
            b = 0;
        }
        else
        {
            a = INT_MAX;
            b = INT_MAX;
        }
    }
    else
    {
        vector<ArrayXd> Ab;
        Array4i row;
        row << 0, 1, 2, 0;
        Array2i j;
        MatrixXd w1_i, w11_i;
        Vector2d w0_i;
        w0_i.setZero();
        for (int i = 0; i < 3; i++)
        {
            j.setLinSpaced(2, row(i), row(i + 1));
            w1_i.setZero(2, 2);
            w11_i.setZero(2, 2);
            w1_i = w1(j, all);
            for (int k = 0; k < w1_i.rows(); k++)
            {
                for (int l = 0; l < w1_i.cols(); l++)
                    w11_i(k, l) = round(w1_i(k, l) * 1e8) / 1e8;
            }
            w0_i = w0(j, 0);
            FullPivLU<Matrix2d> LU(w11_i);
            if (LU.rank() == 2)
            {
                ab1 = (w0_i.transpose() * w1_i.transpose().inverse()).transpose();
                Ab.push_back(ab1);
            }
        }

        MatrixXd AB(2, Ab.size()), AB1, ABi;
        AB.setZero();
        for (int i = 0; i < Ab.size(); i++)
            AB(all, i) = Ab[i];


        if (AB.size() == 0)
        {
            a = INT_MAX;
            b = INT_MAX;
        }
        else
        {
            ArrayXd issame;
            issame.setZero(AB.cols());
            AB1.setZero(AB(all, 0).rows(), 1);
            for (int i = 0; i < AB(all, 0).size(); i++)
                AB1(i, 0) = round(AB(i, 0) * 100000000) / 100000000;
            for (int i = 0; i < AB.cols(); i++)
            {
                ABi.setZero(AB.rows(), 1);
                for (int j = 0; j < AB.rows(); j++)
                    ABi(j, 0) = round(AB(j, i) * 100000000) / 100000000;
                if (AB1(0, 0) == ABi(0, 0) && AB1(1, 0) == ABi(1, 0))
                    issame(i) = 1;
            }
            if (issame.minCoeff() == 1)
            {
                a = AB(0, 0);
                b = AB(1, 0);
            }
            else
            {
                a = INT_MAX;
                b = INT_MAX;
            }
        }
    }

    return make_tuple(a, b);
}