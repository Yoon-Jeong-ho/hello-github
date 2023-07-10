#pragma once

#include "HigherOrderBeamLib.h"
#include "BeamNodeDof.h"
#include "LinearAlgebra.h"

tuple< vector<VectorXi>, vector<MatrixXd>, int> BeamBoundCs(MatrixXd BoundCs, VectorXi Nmode, VectorXi Nels, VectorXi sbdof, vector<vector<MatrixXd>> coef_s, RowVectorXi csnum_s, MatrixXd pnt_s, MatrixXi pnt_con_s, MatrixXd xdir_s, VectorXi Mset)
{
    vector<VectorXi> FixDofs;
    VectorXi fix1, fixi;
    vector<int> fix2;
    MatrixXd fixedXY, bpnt;
    vector<MatrixXd> fixed_xy;
    Vector3d bzdir, bxdir, Udir, fj;
    Vector3d pnt_void; pnt_void.setZero();
    ArrayXd Fix3d;
    int bi, ni, bn, fixsize = 0, cf;

    for (int i = 0; i < BoundCs.rows(); i++)
    {
        bi = BoundCs(i, 1);
        ni = BoundCs(i, 2);
        cf = coef_s[csnum_s(bi)][0].cols();
        fix1 = BeamNodeDof(Nmode, Nels, sbdof, bi, ni);
        fix2.clear();
        if (BoundCs(i, 0) == 0)
        {
            bn = BoundCs(i, 3);
            while (bn > 0) { fix2.push_back((bn % 10) - 1); bn = bn / 10; }
            sort(fix2.begin(), fix2.end());
            if (fix2.size() == 6)
            {
                fixi = fix1;
                FixDofs.push_back(fixi);
                fixsize += fixi.size();

            }
            else if (fix2.size() != 6)
            {
                fixi.setZero(fix2.size());
                for (int j = 0; j < fix2.size(); j++)
                {
                    fixedXY.setZero(3, 2);
                    if (fix2[j] == 0)
                    {
                        Udir = { 1,0,0 };
                        bpnt = pnt_s(pnt_con_s(bi, all), all);
                        bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                        bxdir = xdir_s(bi, all);

                        Fix3d = axistrans(Udir, pnt_void, bzdir, bxdir, 1);
                        fj = { double(fix1(0)),double(fix1(1)),double(fix1(cf)) };
                        fixedXY(all, 0) = fj;
                        fixedXY(all, 1) = Fix3d;
                    }
                    else if (fix2[j] == 1)
                    {
                        Udir = { 0,1,0 };
                        bpnt = pnt_s(pnt_con_s(bi, all), all);
                        bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                        bxdir = xdir_s(bi, all);

                        Fix3d = axistrans(Udir, pnt_void, bzdir, bxdir, 1);
                        fj = { double(fix1(0)),double(fix1(1)),double(fix1(cf)) };
                        fixedXY(all, 0) = fj;
                        fixedXY(all, 1) = Fix3d;
                    }
                    else if (fix2[j] == 2)
                    {
                        Udir = { 0,0,1 };
                        bpnt = pnt_s(pnt_con_s(bi, all), all);
                        bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                        bxdir = xdir_s(bi, all);

                        Fix3d = axistrans(Udir, pnt_void, bzdir, bxdir, 1);
                        fj = { double(fix1(0)),double(fix1(1)),double(fix1(cf)) };
                        fixedXY(all, 0) = fj;
                        fixedXY(all, 1) = Fix3d;
                    }
                    else if (fix2[j] == 3)
                    {
                        Udir = { 1,0,0 };
                        bpnt = pnt_s(pnt_con_s(bi, all), all);
                        bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                        bxdir = xdir_s(bi, all);

                        Fix3d = axistrans(Udir, pnt_void, bzdir, bxdir, 1);
                        fj = { double(fix1(cf + 1)),double(fix1(cf + 2)),double(fix1(2)) };
                        fixedXY(all, 0) = fj;
                        fixedXY(all, 1) = Fix3d;
                    }
                    else if (fix2[j] == 4)
                    {
                        Udir = { 0,1,0 };
                        bpnt = pnt_s(pnt_con_s(bi, all), all);
                        bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                        bxdir = xdir_s(bi, all);

                        Fix3d = axistrans(Udir, pnt_void, bzdir, bxdir, 1);
                        fj = { double(fix1(cf + 1)),double(fix1(cf + 2)),double(fix1(2)) };
                        fixedXY(all, 0) = fj;
                        fixedXY(all, 1) = Fix3d;
                    }
                    else if (fix2[j] == 5)
                    {
                        Udir = { 0,0,1 };
                        bpnt = pnt_s(pnt_con_s(bi, all), all);
                        bzdir = normalize(bpnt(1, all) - bpnt(0, all));
                        bxdir = xdir_s(bi, all);

                        Fix3d = axistrans(Udir, pnt_void, bzdir, bxdir, 1);
                        fj = { double(fix1(cf + 1)),double(fix1(cf + 2)),double(fix1(2)) };
                        fixedXY(all, 0) = fj;
                        fixedXY(all, 1) = Fix3d;
                    }
                    fixed_xy.push_back(fixedXY);
                }
            }
        }
        else if (BoundCs(i, 0) == 1)
        {
            if (BoundCs(i, 3) >= Mset(csnum_s(bi)))
            {
                fixi = fix1(VectorXi::LinSpaced(Nmode(bi) - cf - 3, cf + 3, Nmode(bi)));
                FixDofs.push_back(fixi);
                fixsize += fixi.size();
            }
        }
        else if (BoundCs(i, 0) == 2)
        {
            if (BoundCs(i, 3) >= Mset(csnum_s(bi)))
            {
                fixi = fix1(VectorXi::LinSpaced(cf - 3, 3, cf));
                FixDofs.push_back(fixi);
                fixsize += fixi.size();
            }
        }
    }

    return make_tuple(FixDofs, fixed_xy, fixsize);
}