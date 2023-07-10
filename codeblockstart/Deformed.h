#pragma once

#include "HigherOrderBeamLib.h"
#include "Edge.h"
#include "PSI.h"

tuple<MatrixXi, MatrixXd, MatrixXd, double> deformed(MatrixXd csc, MatrixXi cscc, vector<MatrixXd> coef_s, vector<MatrixXd> coef_n, vector<MatrixXd> coef_z, VectorXd D, VectorXd zcoord, double sm)
{
    ArrayXd alpha, Le, Nes;
    Edge edge;
    int Nmode = coef_s[0].cols() + coef_z[0].cols(), Nnz, Nez;
    double dim = max(csc(all, 0).maxCoeff() - csc(all, 0).minCoeff(), csc(all, 1).maxCoeff() - csc(all, 1).minCoeff());
    MatrixXd di;
    vector<int> ji;
    double scw;

    alpha = edge.Angle(csc, cscc);
    Le = edge.Length(csc, cscc);

    class Transform
    {
    public:
        Matrix3d T1(double a)
        {
            Matrix3d T2;
            T2 << cos(a), sin(a), 0,
                sin(a), -cos(a), 0,
                0, 0, 1;

            return T2;
        }
    };

    Nnz = zcoord.size();
    Nez = Nnz - 1;
    Nes = ceil(Le / sm);
    di.setZero(Nnz, Nmode);
    for (int j = 0; j < Nmode; j++)
    {
        ji.clear();
        for (int k = j; k < D.size() - 1; k = k + Nmode * 2)
            ji.push_back(k);
        di(all, j) = D(ji);
    }

    MatrixXi elnode;
    MatrixXd nodcoord, noddisp, F;
    int see, sne, k, ens = 0, ncs, nds, maxno;
    double ae;
    Vector2i ec;
    Vector2d cnre;
    Vector3d cnre0;
    MatrixXi en;
    MatrixXd nc, nd, f, nct, nctt, ndt;
    VectorXd scoord, row;
    Transform T;
    vector<MatrixXi> envec;

    for (int e = 0; e < Le.size(); e++)
    {
        see = Nes(e);
        sne = see + 1;
        ae = alpha(e);
        ec = cscc(e, all).transpose();
        cnre = csc(ec(0), all);

        // Element Node Connectivity
        en.setZero(see * Nez, 4);
        for (int i = 0; i < see; i++)
        {
            for (int j = 0; j < Nez; j++)
            {
                k = i * Nez + j;
                en(k, 0) = k + i;
                en(k, 1) = k + 1 + i;
                en(k, 2) = k + 1 + Nnz + i;
                en(k, 3) = k + Nnz + i;
            }
        }
        if (elnode.size() == 0)
            maxno = 0;
        else
            maxno = elnode.maxCoeff() + 1;
        ens = elnode.rows();
        //envec.push_back(en);
        elnode.conservativeResize(ens + en.rows(), en.cols());
        elnode.block(ens, 0, en.rows(), en.cols()) = en + MatrixXi::Ones(en.rows(), en.cols()) * maxno;

        // Node Coordinate and Displacement
        scoord.setLinSpaced(sne, 0, Le(e));
        F = PSI(coef_s, coef_n, coef_z, e, scoord, 0);
        nc.setZero(Nnz * sne, 3);
        nd.setZero(Nnz * sne, 3);

        for (int i = 0; i < sne; i++)
        {
            f = F.block(0, Nmode * i, 3, Nmode);
            row.setLinSpaced(Nnz, i * Nnz, (i + 1) * Nnz - 1);
            nc(row, 0) = scoord(i) * VectorXd::Ones(Nnz, 1);
            nc(row, 1).setZero();
            nc(row, 2) = zcoord;
            nd(row, all) = di * f.transpose();
        }

        nct = nc * T.T1(ae).transpose();
        cnre0 << cnre,
            0;
        nctt = nct + VectorXd::Ones(Nnz * sne) * cnre0.transpose();
        ncs = nodcoord.rows();
        nodcoord.conservativeResize(nodcoord.rows() + nctt.rows(), nctt.cols());
        nodcoord.block(ncs, 0, nctt.rows(), nctt.cols()) = nctt;
        ndt = nd * T.T1(ae).transpose();
        nds = noddisp.rows();
        noddisp.conservativeResize(noddisp.rows() + ndt.rows(), ndt.cols());
        noddisp.block(nds, 0, ndt.rows(), ndt.cols()) = ndt;
    }

    scw = dim / (5 * noddisp.cwiseAbs().maxCoeff());

    return make_tuple(elnode, nodcoord, noddisp, scw);
}