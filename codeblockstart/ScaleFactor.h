#pragma once

#include "HigherOrderBeamLib.h"

Matrix<double50, -1, -1> ScaleFactor(Array<double50, -1, 1> Le, Matrix<double50, -1, -1> coef, Matrix<double50, -1, -1> M)
{
    int Nx = coef.cols();
    Matrix<double50, -1, -1> sf(Nx, Nx), cofT, cof, recip;
    sf.setZero();
    for (int i = 0;i < Nx;i++)
    {
        cofT = coef(all, i).transpose();
        cof = coef(all, i);
        recip = cofT * M * cof;
        sf(i, i) = sqrt(Le.sum() / recip(0, 0));
    }

    return sf;
}
