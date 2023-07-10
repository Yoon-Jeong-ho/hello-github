#pragma once
#define _USE_MATH_DEFINES

#include "HigherOrderBeamLib.h"

tuple<ArrayXd, ArrayXd> LGQ(float N)
{
    double a = -1, b = 1;
    double Nn = N - 1;
    double N1 = Nn + 1, N2 = Nn + 2;
    ArrayXd xu, x0;
    xu.setLinSpaced(N1, -1, 1);
    x0.setLinSpaced(Nn + 1, 0, Nn);

    MatrixXd L, Lp, y, y0, one, w;
    y = cos((2 * x0 + 1) * M_PI / (2 * Nn + 2)) + (0.27 / N1) * sin(M_PI * xu * Nn / N2);
    L.setZero(N1, N2);
    one.setOnes(L.rows(), 1);
    y0 = y0.setOnes(y.rows(), 1) * 2;

    while ((y - y0).cwiseAbs().maxCoeff() > DBL_EPSILON)
    {
        L.block(0, 0, L.rows(), 1).setOnes();
        L.block(0, 1, L.rows(), 1) = y;

        for (int k = 1;k < N1;k++)
        {
            L.block(0, k + 1, L.rows(), 1) = ((2 * k + 1) * y.cwiseProduct(L.block(0, k, L.rows(), 1)) - k * L.block(0, k - 1, L.rows(), 1)) / (k + 1);
        }
        Lp = N2 * (L.block(0, N1 - 1, L.rows(), 1) - y.cwiseProduct(L.block(0, N2 - 1, L.rows(), 1))).cwiseProduct((one - y.cwiseProduct(y)).cwiseInverse());

        y0 = y;
        y = y0 - L.block(0, N2 - 1, L.rows(), 1).cwiseProduct(Lp.cwiseInverse());
    }
    w = (b - a) * one.cwiseProduct(((one - y.cwiseProduct(y)).cwiseProduct(Lp.cwiseProduct(Lp))).cwiseInverse()) * pow(N2 / N1, 2);

    return make_tuple(y.reverse(), w);

}