#pragma once

#include "HigherOrderBeamLib.h"
#include "BeamNodeDof.h"
#include "LinearAlgebra.h"
#include "Edge.h"
#include "PSI.h"

VectorXd Force(MatrixXd FNode, MatrixXd FEdge, VectorXi Nmode, VectorXi Nels, VectorXi sbdof, vector<vector<MatrixXd>> coef_s, vector<vector<MatrixXd>> coef_n, vector<vector<MatrixXd>> coef_z, RowVectorXi csnum_s, MatrixXd pnt_s, MatrixXi pnt_con_s, MatrixXd xdir_s, vector<MatrixXd> cscd, vector<MatrixXi> cscc, vector<VectorXd> zcoord, int ks, MatrixXi ShellNdof, vector<Matrix3d> shellTransform)
{
    VectorXd Force3d, F; F.setZero(ks);

    VectorXi forcedofi, dofe;
    Vector3i Forcedof;
    MatrixXd spnt;
    Vector3d sdir, sdir1, sdir2, psd1, psd2;
    int bi, ni, cf, di, dofix, ei;
    ArrayXd six; six.setZero(1);
    MatrixXd bpnt;
    Vector3d bzdir, bxdir;
    Vector3d pnt_void, Fdir; pnt_void.setZero();
    VectorXd sFdir;

    if (FNode.rows() != 0)
    {
        for (int i = 0; i < FNode.rows(); i++)
        {
            bi = FNode(i, 1);
            if (bi == -1)
            {
                ni = FNode(i, 2);
                di = FNode(i, 3);
                for (int j = 0; j < ShellNdof.rows(); j++)
                {
                    if (ni == ShellNdof(j, 1))
                    {
                        sFdir.setZero(6);
                        sFdir(di - 1) = FNode(i, 4);
                        forcedofi = ShellNdof(j, { 2,3,4,5,6,7 });
                        break;
                    }
                }
                F(forcedofi) = sFdir;
            }
            else
            {
                ni = FNode(i, 2);
                cf = coef_s[csnum_s(bi)][0].cols();
                di = FNode(i, 3);
                forcedofi = BeamNodeDof(Nmode, Nels, sbdof, bi, ni);
                if (FNode(i, 0) == 0)
                {
                    if (di == 1)
                    {
                        Fdir = { 1,0,0 };
                        bpnt = pnt_s(pnt_con_s(bi, all), all);
                        bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                        bxdir = xdir_s(bi, all);

                        Force3d = axistrans(Fdir, pnt_void, bzdir, bxdir, 1) * FNode(i, 4);
                        Forcedof = { forcedofi(0),forcedofi(1),forcedofi(cf) };
                    }
                    else if (di == 2)
                    {
                        Fdir = { 0,1,0 };
                        bpnt = pnt_s(pnt_con_s(bi, all), all);
                        bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                        bxdir = xdir_s(bi, all);

                        Force3d = axistrans(Fdir, pnt_void, bzdir, bxdir, 1) * FNode(i, 4);
                        Forcedof = { forcedofi(0),forcedofi(1),forcedofi(cf) };
                    }
                    else if (di == 3)
                    {
                        Fdir = { 0,0,1 };
                        bpnt = pnt_s(pnt_con_s(bi, all), all);
                        bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                        bxdir = xdir_s(bi, all);

                        Force3d = axistrans(Fdir, pnt_void, bzdir, bxdir, 1) * FNode(i, 4);
                        Forcedof = { forcedofi(0),forcedofi(1),forcedofi(cf) };
                    }
                    else if (di == 4)
                    {
                        Fdir = { 1,0,0 };
                        bpnt = pnt_s(pnt_con_s(bi, all), all);
                        bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                        bxdir = xdir_s(bi, all);

                        Force3d = axistrans(Fdir, pnt_void, bzdir, bxdir, 1) * FNode(i, 4);
                        Forcedof = { forcedofi(cf + 1),forcedofi(cf + 2),forcedofi(2) };
                    }
                    else if (di == 5)
                    {
                        Fdir = { 0,1,0 };
                        bpnt = pnt_s(pnt_con_s(bi, all), all);
                        bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                        bxdir = xdir_s(bi, all);

                        Force3d = axistrans(Fdir, pnt_void, bzdir, bxdir, 1) * FNode(i, 4);
                        Forcedof = { forcedofi(cf + 1),forcedofi(cf + 2),forcedofi(2) };
                    }
                    else if (di == 6)
                    {
                        Fdir = { 0,0,1 };
                        bpnt = pnt_s(pnt_con_s(bi, all), all);
                        bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                        bxdir = xdir_s(bi, all);

                        Force3d = axistrans(Fdir, pnt_void, bzdir, bxdir, 1) * FNode(i, 4);
                        Forcedof = { forcedofi(cf + 1),forcedofi(cf + 2),forcedofi(2) };
                    }
                    F(Forcedof) += Force3d;
                }
                else if (FNode(i, 0) == 1)
                {
                    dofix = forcedofi(cf + di + 2);
                    F(dofix) = FNode(i, 4);
                }
                else if (FNode(i, 0) == 2)
                {
                    dofix = forcedofi(di + 2);
                    F(dofix) = FNode(i, 4);
                }
            }

        }
    }

    Edge len;
    ArrayXd eL;

    if (FEdge.rows() != 0)
    {
        for (int i = 0; i < FEdge.rows(); i++)
        {
            bi = FEdge(i, 0);
            ni = FEdge(i, 1);
            cf = coef_s[csnum_s(bi)][0].cols();
            di = FEdge(i, 4);
            ei = FEdge(i, 2) - 1;
            eL = len.Length(cscd[csnum_s(bi)], cscc[csnum_s(bi)]);
            six[0] = FEdge(i, 3) * eL(ei);
            forcedofi = BeamNodeDof(Nmode, Nels, sbdof, bi, ni);
            dofe = forcedofi(VectorXi::LinSpaced(Nmode(bi), 0, Nmode(bi) - 1));
            if (di == 1)
            {
                Fdir = { 1,0,0 };
                bpnt = pnt_s(pnt_con_s(bi, all), all);
                bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                bxdir = xdir_s(bi, all);

                spnt = cscd[csnum_s(bi)](cscc[csnum_s(bi)](ei, all), all);
                sdir1.setZero(), sdir2.setZero();
                sdir1({ 0,1 }) = spnt(0, all);
                sdir1(2) = zcoord[bi](ni);
                sdir2({ 0,1 }) = spnt(1, all);
                sdir2(2) = zcoord[bi](ni);
                psd1 = axistrans(sdir1, pnt_void, bzdir, bxdir, 0);
                psd2 = axistrans(sdir2, pnt_void, bzdir, bxdir, 0);
                sdir = normalize(psd2 - psd1);

                Force3d = axistrans(Fdir, pnt_void, bzdir, sdir, 1) * FEdge(i, 5);
                F(dofe) += (Force3d.transpose() * PSI(coef_s[csnum_s(bi)], coef_n[csnum_s(bi)], coef_z[csnum_s(bi)], ei, six, 0)).transpose();
            }
            else if (di == 2)
            {
                Fdir = { 0,1,0 };
                bpnt = pnt_s(pnt_con_s(bi, all), all);
                bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                bxdir = xdir_s(bi, all);

                spnt = cscd[csnum_s(bi)](cscc[csnum_s(bi)](ei, all), all);
                sdir1.setZero(), sdir2.setZero();
                sdir1({ 0,1 }) = spnt(0, all);
                sdir1(2) = zcoord[bi](ni);
                sdir2({ 0,1 }) = spnt(1, all);
                sdir2(2) = zcoord[bi](ni);
                psd1 = axistrans(sdir1, pnt_void, bzdir, bxdir, 0);
                psd2 = axistrans(sdir2, pnt_void, bzdir, bxdir, 0);
                sdir = normalize(psd2 - psd1);

                Force3d = axistrans(Fdir, pnt_void, bzdir, sdir, 1) * FEdge(i, 5);
                F(dofe) += (Force3d.transpose() * PSI(coef_s[csnum_s(bi)], coef_n[csnum_s(bi)], coef_z[csnum_s(bi)], ei, six, 0)).transpose();
            }
            else if (di == 3)
            {
                Fdir = { 0,0,1 };
                bpnt = pnt_s(pnt_con_s(bi, all), all);
                bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                bxdir = xdir_s(bi, all);

                spnt = cscd[csnum_s(bi)](cscc[csnum_s(bi)](ei, all), all);
                sdir1.setZero(), sdir2.setZero();
                sdir1({ 0,1 }) = spnt(0, all);
                sdir1(2) = zcoord[bi](ni);
                sdir2({ 0,1 }) = spnt(1, all);
                sdir2(2) = zcoord[bi](ni);
                psd1 = axistrans(sdir1, pnt_void, bzdir, bxdir, 0);
                psd2 = axistrans(sdir2, pnt_void, bzdir, bxdir, 0);
                sdir = normalize(psd2 - psd1);

                Force3d = axistrans(Fdir, pnt_void, bzdir, sdir, 1) * FEdge(i, 5);
                F(dofe) += (Force3d.transpose() * PSI(coef_s[csnum_s(bi)], coef_n[csnum_s(bi)], coef_z[csnum_s(bi)], ei, six, 0)).transpose();
            }
            else if (di == 4 || di == 5 || di == 6)
            {
                cout << "Error: Only point force in x, y or z directions can be applied on Edge!!" << endl << endl;
                throw std::exception();
            }
        }
    }

    return F;
}
