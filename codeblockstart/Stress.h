#pragma once

#include "HigherOrderBeamLib.h"
#include "Edge.h"
#include "PSI.h"
#include "Hermite.h"

MatrixXd stress3D(MatrixXd csc, MatrixXi cscc, vector<MatrixXd> coef_s, vector<MatrixXd> coef_n, vector<MatrixXd> coef_z, double n, VectorXd D, VectorXd zcoord, double sm, double E, double nu)
{
    MatrixXd stress(0, 0);
    ArrayXd alpha, Le, Nes;
    int Nmode = coef_s[0].cols() + coef_z[0].cols(), Nnz, Nez, Ndof = Nmode * 2;
    double E1;
    Edge edge;

    alpha = edge.Angle(csc, cscc);
    Le = edge.Length(csc, cscc);

    Nnz = zcoord.size();
    Nez = Nnz - 1;
    Nes = ceil(Le / sm);

    Matrix3d C;
    if (Nmode == 6)
        C << E, 0, 0,
        0, E, 0,
        0, 0, E / 2 / (1 + nu);
    else
    {
        E1 = E / (1 - pow(nu, 2));
        C << 1, nu, 0,
            nu, 1, 0,
            0, 0, (1 - nu) / 2;
        C *= E1;
    }

    MatrixXd d;
    d.setZero(Ndof * 2, Nez);
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
    VectorXi rng;
    for (int i = 0; i < Nez; i++)
    {
        rng.setLinSpaced(Ndof * 2, i * Ndof, i * Ndof + Ndof * 2 - 1);
        d(all, i) = D(rng);
    }

    MatrixXd L1(2, 3), L2(2, 3), L3(2, 3), L4(3, 2), L5(3, 2), La, Lb, Lc, Ld, Lf;
    L1 << 1, 0, 0,
        0, 0, 1;
    L2 << 0, 1, 0,
        0, 0, 0;
    L3 << 0, 0, 0,
        0, 1, 0;
    L4 << 1, 0,
        0, 0,
        0, 1;
    L5 << 0, 0,
        0, 1,
        1, 0;

    La = L4 * L1, Lb = L5 * L1, Lc = L4 * L2, Ld = L4 * L3 + L5 * L2, Lf = L5 * L3;

    int Ns, col, si = 0;
    VectorXd scoord, dj;
    MatrixXd F, dF, ddF, strs_e, f, df, ddf, strain;
    double Li, k;
    Hermite herm;

    for (int e = 0; e < Le.size(); e++)
    {
        Ns = Nes(e) + 1;
        scoord.setLinSpaced(Ns, 0, Le(e));
        F = PSI(coef_s, coef_n, coef_z, e, scoord, 0);
        dF = PSI(coef_s, coef_n, coef_z, e, scoord, 1);
        ddF = PSI(coef_s, coef_n, coef_z, e, scoord, 2);
        strs_e.setZero(3, Nnz * Ns);
        for (int i = 0; i < Ns; i++)
        {
            f = F.block(0, Nmode * i, 3, Nmode);
            df = dF.block(0, Nmode * i, 3, Nmode);
            ddf = ddF.block(0, Nmode * i, 3, Nmode);
            for (int j = 0; j < Nnz; j++)
            {
                if (j == 0)
                {
                    Li = zcoord(1) - zcoord(0);
                    dj = d(all, 0);
                    k = -1;
                }
                else
                {
                    Li = zcoord(j) - zcoord(j - 1);
                    dj = d(all, j - 1);
                    k = 1;
                }
                col = i * Nnz + j;
                strain = La * df * herm.N(Li, k, Nmode) + Lb * f * herm.dN(Li, k, Nmode) - n * Lc * ddf * herm.N(Li, k, Nmode) - n * Ld * df * herm.dN(Li, k, Nmode) - n * Lf * f * herm.ddN(Li, k, Nmode);
                strs_e(all, col) = C * strain * dj;
            }
        }
        si = stress.rows();
        stress.conservativeResize(strs_e.cols() + si, 3);
        stress.block(si, 0, strs_e.cols(), strs_e.rows()) = strs_e.transpose();
    }

    return stress;
}