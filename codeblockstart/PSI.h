#pragma once

#include "HigherOrderBeamLib.h"
#include "Edge.h"

MatrixXd PSI(vector<MatrixXd> coef_s, vector<MatrixXd> coef_n, vector<MatrixXd> coef_z, int e, ArrayXd s, int d)
{
    int NMi = coef_s[0].cols(), NMo = coef_z[0].cols();
    MatrixXd fs, fn, fz;
    Edge edge;

    fs = edge.Continuity2(coef_s, e, s, d);
    fn = edge.Continuity2(coef_n, e, s, d);
    fz = edge.Continuity2(coef_z, e, s, d);

    MatrixXd fval;
    fval.setZero(3, (NMi + NMo) * s.rows());
    for (int i = 0; i < s.rows(); i++)
    {
        fval.block(0, (NMi + NMo) * i, 1, NMi) = fs.block(i, 0, 1, NMi);
        fval.block(1, (NMi + NMo) * i, 1, NMi) = fn.block(i, 0, 1, NMi);
        fval.block(2, (NMi + NMo) * i + NMi, 1, NMo) = fz.block(i, 0, 1, NMo);
    }

    return fval;
}
