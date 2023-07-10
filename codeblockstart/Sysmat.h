#pragma once

#include "HigherOrderBeamLib.h"
#include "Edge.h"
#include "Hermite.h"
#include "PSI.h"
#include "LGQ.h"

tuple<vector<MatrixXi>, vector<MatrixXi>, vector<MatrixXi>, vector<MatrixXi>, vector<MatrixXd>, vector<MatrixXd>, VectorXd, VectorXi> sysmat(vector<MatrixXd> csc, vector<MatrixXi> cscc, vector<vector<MatrixXd>> coef_s, vector<vector<MatrixXd>> coef_n, vector<vector<MatrixXd>> coef_z, MatrixXd pnt, MatrixXi pnt_con, ArrayXi csnum, MatrixXd besz, ArrayXd t, MatrixXd enr, vector<VectorXd> zcoord, VectorXi Nel)
{
    double tol1 = 1e-8, tol2 = 8;
    int Nb = pnt_con.rows();
    vector<MatrixXi> Ki, Kj, Mi, Mj;
    vector<MatrixXd> Kall, Mall;
    VectorXd L(Nb);
    VectorXi Nmode(Nb);
    MatrixXd pnti, csci, KK, MM, MM1;
    double Li;
    int cs, rck = 0, bdof = 0;
    MatrixXi cscci, ki, kj, mi, mj;
    vector<MatrixXd> coef_si, coef_ni, coef_zi;
    VectorXd zcoordi, Le;
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

    for (int i = 0; i < Nb; i++)
    {
        cout << "BEAM NO.:  " << i << endl;
        pnti = pnt(pnt_con(i, all), all);
        Li = (pnti(1, all) - pnti(0, all)).norm();
        L(i) = Li;
        cs = csnum(i);
        csci = csc[cs];
        cscci = cscc[cs];
        coef_si = coef_s[cs];
        coef_ni = coef_n[cs];
        coef_zi = coef_z[cs];
        Nmode(i) = coef_si[0].cols() + coef_zi[0].cols();
        zcoordi = zcoord[i];
        Le = zcoordi.block(1, 0, zcoordi.size() - 1, 1) - zcoordi.block(0, 0, zcoordi.size() - 1, 1);

        tie(ki, kj, KK) = stiffness(csci, cscci, coef_si, coef_ni, coef_zi, t(i), Li, enr(i, 0), enr(i, 1), Le, INT_MAX);
        tie(mi, mj, MM1) = Mass(csci, cscci, coef_si, coef_ni, coef_zi, t(i), Li, enr(i, 0), enr(i, 1), Le, INT_MAX);
        MM = MM1 * enr(i, 2);

        rck += KK.rows();

        Kall.push_back(KK);
        Mall.push_back(MM);

        if (i == 0)
        {
            Ki.push_back(ki);
            Kj.push_back(kj);
            Mi.push_back(mi);
            Mj.push_back(mj);
        }
        else
        {
            Ki.push_back(ki + MatrixXi::Ones(ki.rows(), ki.cols()) * bdof);
            Kj.push_back(kj + MatrixXi::Ones(kj.rows(), kj.cols()) * bdof);
            Mi.push_back(mi + MatrixXi::Ones(mi.rows(), mi.cols()) * bdof);
            Mj.push_back(mj + MatrixXi::Ones(mj.rows(), mj.cols()) * bdof);
        }

        bdof += (Nel(i) + 1) * Nmode(i) * 2;

    }
    cout << endl;

    return make_tuple(Ki, Kj, Mi, Mj, Kall, Mall, L, Nmode);
}

tuple<MatrixXi, MatrixXi, MatrixXd> stiffness(MatrixXd csc, MatrixXi cscc, vector<MatrixXd> coef_s, vector<MatrixXd> coef_n, vector<MatrixXd> coef_z, double t, double L, double E, double nu, VectorXd Li1, int Neli)
{
    //std::cout << "Stiffness Matrix Computation Starts:" << endl << endl;
    int NE = cscc.rows(), tol = 2;
    int NMs = coef_s[0].cols(), NMz = coef_z[0].cols(), Nel, Nmode = NMs + NMz, Ndof = Nmode * 2;
    double E1 = E / (1 - pow(nu, 2));
    ArrayXd alpha, Le, Lei;
    VectorXd Li;
    MatrixXd Di(3, 3), D(3, 3);
    Edge edge;
    Order Ord;

    if (Neli != INT_MAX)
    {
        Nel = Neli;
        Li = Lei.Ones(Nel, 1) * L / Nel;
    }
    else
    {
        Nel = Li1.size();
        Li = Li1;
    }
    int tdof = Ndof * (Nel + 1);

    if (Nmode == 6)
    {
        Di << E, 0, 0,
            0, E, 0,
            0, 0, E / (2 * (1 - nu));
    }
    else
    {
        Di << 1, nu, 0,
            nu, 1, 0,
            0, 0, (1 - nu) / 2;
    }
    D = E1 * Di;
    alpha = edge.Angle(csc, cscc);
    Le = edge.Length(csc, cscc);

    MatrixXi ord(NE, NMs * 2 + NMz);
    ord.setZero();
    ord.block(0, 0, NE, NMs) = Ord.order(coef_s);
    ord.block(0, NMs, NE, NMs) = Ord.order(coef_n);
    ord.block(0, NMs * 2, NE, NMz) = Ord.order(coef_z);

    float degs = ord.maxCoeff() * 2, degz = 3 * 2;
    float Ngps = ceil((degs + 1) / 2) + tol;
    float Ngpz = ceil((degz + 1) / 2) + tol;

    ArrayXd GPs, WFs, GPz, WFz;
    std::tie(GPs, WFs) = LGQ(Ngps);
    std::tie(GPz, WFz) = LGQ(Ngpz);

    MatrixXi elnode(Nel, 2), one1;
    for (int i = 0;i < Nel;i++)
        elnode.block(i, 0, 1, elnode.cols()) << i + 1, i + 2;
    one1.setOnes(elnode.rows(), 1);

    MatrixXi eldof(Nel, Ndof * 2);
    eldof.setZero();
    for (int i = 0; i < Ndof; i++)
    {
        eldof.block(0, i, Nel, 1) = elnode.block(0, 0, Nel, 1) * Ndof - one1 * (Ndof - i);
        eldof.block(0, i + Ndof, Nel, 1) = elnode.block(0, 1, Nel, 1) * Ndof - one1 * (Ndof - i);
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

    MatrixXd D1, D2, D3, D4, D5, D6, D7, D8, D9;
    D1 = La.transpose() * D * La;
    D2 = La.transpose() * D * Lb;
    D3 = Lb.transpose() * D * Lb;
    D4 = Lc.transpose() * D * Lc;
    D5 = Lc.transpose() * D * Ld;
    D6 = Ld.transpose() * D * Ld;
    D7 = Lc.transpose() * D * Lf;
    D8 = Ld.transpose() * D * Lf;
    D9 = Lf.transpose() * D * Lf;

    MatrixXd H1, H2, H3, H4, H5, H6, H7, H8, H9;
    H1.setZero(Nmode, Nmode);
    H2.setZero(Nmode, Nmode);
    H3.setZero(Nmode, Nmode);
    H4.setZero(Nmode, Nmode);
    H5.setZero(Nmode, Nmode);
    H6.setZero(Nmode, Nmode);
    H7.setZero(Nmode, Nmode);
    H8.setZero(Nmode, Nmode);
    H9.setZero(Nmode, Nmode);

    double tt = pow(t, 3) / 12, jcob, wf;
    ArrayXd s;
    MatrixXd F, dF, ddF, f, df, ddf;

    for (int i = 0; i < NE; i++)
    {
        jcob = Le(i) / 2;
        s = jcob * (GPs + 1);

        F = PSI(coef_s, coef_n, coef_z, i, s, 0);
        dF = PSI(coef_s, coef_n, coef_z, i, s, 1);
        ddF = PSI(coef_s, coef_n, coef_z, i, s, 2);

        for (int j = 0; j < Ngps; j++)
        {
            wf = WFs(j);
            f = F.block(0, (NMs + NMz) * j, 3, NMs + NMz);
            df = dF.block(0, (NMs + NMz) * j, 3, NMs + NMz);
            ddf = ddF.block(0, (NMs + NMz) * j, 3, NMs + NMz);

            H1 += df.transpose() * D1 * df * jcob * wf * t;
            H2 += df.transpose() * D2 * f * jcob * wf * t;
            H3 += f.transpose() * D3 * f * jcob * wf * t;
            H4 += ddf.transpose() * D4 * ddf * jcob * wf * tt;
            H5 += ddf.transpose() * D5 * df * jcob * wf * tt;
            H6 += df.transpose() * D6 * df * jcob * wf * tt;
            H7 += ddf.transpose() * D7 * f * jcob * wf * tt;
            H8 += df.transpose() * D8 * f * jcob * wf * tt;
            H9 += f.transpose() * D9 * f * jcob * wf * tt;
        }
    }
    vector<double> Lij;
    for (int i = 0; i < Li.size(); i++)
        Lij.push_back(Li(i));

    auto end = Lij.end();
    for (auto it = Lij.begin(); it != end; ++it)
        end = remove(it + 1, end, *it);

    Lij.erase(end, Lij.end());

    VectorXd Ltype, Const;
    Ltype.setZero(Nel);
    MatrixXd k1, k2, k3, k4, k5, k6, k7, k8, k9;
    vector<MatrixXd> ke;
    double Lm;
    double ksi;
    Hermite shN;
    vector<int> Lij_m;

    for (int m = 0; m < Lij.size(); m++)
    {
        Lm = Lij[m];
        Lij_m.clear();
        for (int ii = 0; ii < Li.size(); ii++)
            if (Li(ii) == Lm)
                Lij_m.push_back(ii);

        Ltype(Lij_m) = Const.setConstant(Lij_m.size(), m);

        k1.setZero(Ndof * 2, Ndof * 2);
        k2.setZero(Ndof * 2, Ndof * 2);
        k3.setZero(Ndof * 2, Ndof * 2);
        k4.setZero(Ndof * 2, Ndof * 2);
        k5.setZero(Ndof * 2, Ndof * 2);
        k6.setZero(Ndof * 2, Ndof * 2);
        k7.setZero(Ndof * 2, Ndof * 2);
        k8.setZero(Ndof * 2, Ndof * 2);
        k9.setZero(Ndof * 2, Ndof * 2);

        jcob = Lm / 2;
        for (int i = 0; i < Ngpz; i++)
        {
            ksi = GPz(i);
            wf = WFz(i);

            k1 += shN.N(Lm, ksi, Nmode).transpose() * H1 * shN.N(Lm, ksi, Nmode) * jcob * wf;
            k2 += shN.N(Lm, ksi, Nmode).transpose() * H2 * shN.dN(Lm, ksi, Nmode) * jcob * wf;
            k3 += shN.dN(Lm, ksi, Nmode).transpose() * H3 * shN.dN(Lm, ksi, Nmode) * jcob * wf;
            k4 += shN.N(Lm, ksi, Nmode).transpose() * H4 * shN.N(Lm, ksi, Nmode) * jcob * wf;
            k5 += shN.N(Lm, ksi, Nmode).transpose() * H5 * shN.dN(Lm, ksi, Nmode) * jcob * wf;
            k6 += shN.dN(Lm, ksi, Nmode).transpose() * H6 * shN.dN(Lm, ksi, Nmode) * jcob * wf;
            k7 += shN.N(Lm, ksi, Nmode).transpose() * H7 * shN.ddN(Lm, ksi, Nmode) * jcob * wf;
            k8 += shN.dN(Lm, ksi, Nmode).transpose() * H8 * shN.ddN(Lm, ksi, Nmode) * jcob * wf;
            k9 += shN.ddN(Lm, ksi, Nmode).transpose() * H9 * shN.ddN(Lm, ksi, Nmode) * jcob * wf;
        }
        ke.push_back(k1 + k2 + k2.transpose() + k3 + k4 + k5 + k5.transpose() + k6 + k7 + k7.transpose() + k8 + k8.transpose() + k9);
    }

    VectorXd Ke(4 * Ndof * Ndof);
    MatrixXi edofs(2 * Ndof, Nel), Ki(2 * Ndof * 2 * Ndof, Nel), Kj(2 * Ndof * 2 * Ndof, Nel);
    MatrixXd Kall(2 * Ndof * 2 * Ndof, Nel);
    Ke.setZero(), edofs = eldof.transpose(), Ki.setZero(), Kj.setZero(), Kall.setZero();
    vector<MatrixXd> Kke;
    VectorXd en;
    for (int k = 0; k < Nel; k++)
    {
        for (int i = 0; i < 2 * Ndof;i++)
        {
            Ke.block(2 * Ndof * i, 0, 2 * Ndof, 1) = ke[Ltype(k)].block(0, i, 2 * Ndof, 1);
        }
        Kke.push_back(Ke);
    }
    for (int i = 0; i < Nel; i++)
        Kall.block(0, i, 2 * Ndof * 2 * Ndof, 1) = Kke[i];
    for (int i = 0; i < 2 * Ndof; i++)
    {
        Ki.block(2 * Ndof * i, 0, 2 * Ndof, Nel) = edofs;
        for (int j = 0; j < 2 * Ndof;j++)
            Kj.block(2 * Ndof * i + j, 0, 1, Nel) = edofs.block(i, 0, 1, Nel);
    }

    //std::cout << "Stiffness Matrix Computed!" << endl << endl;

    return make_tuple(Ki, Kj, Kall);
}

tuple<MatrixXi, MatrixXi, MatrixXd> Mass(MatrixXd csc, MatrixXi cscc, vector<MatrixXd> coef_s, vector<MatrixXd> coef_n, vector<MatrixXd> coef_z, double t, double L, double E, double nu, VectorXd Li1, int Neli)
{
    //std::cout << "Mass Matrix Computation Starts:" << endl << endl;

    int NE = cscc.rows(), tol = 2;
    int NMs = coef_s[0].cols(), NMz = coef_z[0].cols(), Nel, Nmode = NMs + NMz, Ndof = Nmode * 2;
    double E1 = E / (1 - pow(nu, 2));
    ArrayXd alpha, Le, Lei;
    VectorXd Li;
    MatrixXd Di(3, 3), D(3, 3);
    Edge edge;
    Order Ord;

    if (Neli != INT_MAX)
    {
        Nel = Neli;
        Li = Lei.Ones(Nel, 1) * L / Nel;
    }
    else
    {
        Nel = Li1.size();
        Li = Li1;
    }

    int tdof = Ndof * (Nel + 1);
    Di << 1, nu, 0,
        nu, 1, 0,
        0, 0, (1 - nu) / 2;
    D = E1 * Di;

    alpha = edge.Angle(csc, cscc);
    Le = edge.Length(csc, cscc);

    MatrixXi ord(NE, NMs * 2 + NMz);
    ord.setZero();
    ord.block(0, 0, NE, NMs) = Ord.order(coef_s);
    ord.block(0, NMs, NE, NMs) = Ord.order(coef_n);
    ord.block(0, NMs * 2, NE, NMz) = Ord.order(coef_z);


    float degs = ord.maxCoeff() * 2, degz = 3 * 2;
    float Ngps = ceil((degs + 1) / 2) + tol;
    float Ngpz = ceil((degz + 1) / 2) + tol;

    ArrayXd GPs, WFs, GPz, WFz;
    std::tie(GPs, WFs) = LGQ(Ngps);
    std::tie(GPz, WFz) = LGQ(Ngpz);

    MatrixXi elnode(Nel, 2), one1;
    for (int i = 0;i < Nel;i++)
        elnode.block(i, 0, 1, elnode.cols()) << i + 1, i + 2;
    one1.setOnes(elnode.rows(), 1);

    MatrixXi eldof(Nel, Ndof * 2);
    eldof.setZero();
    for (int i = 0; i < Ndof; i++)
    {
        eldof.block(0, i, Nel, 1) = elnode.block(0, 0, Nel, 1) * Ndof - one1 * (Ndof - i);
        eldof.block(0, i + Ndof, Nel, 1) = elnode.block(0, 1, Nel, 1) * Ndof - one1 * (Ndof - i);
    }

    MatrixXd Lm1(3, 3), Lm2(3, 3);
    Lm1 << 0, 1, 0,
        0, 0, 0,
        0, 0, 0;
    Lm2 << 0, 0, 0,
        0, 0, 0,
        0, 1, 0;

    MatrixXd Hm1, Hm2, Hm3, Hm4;
    Hm1.setZero(Nmode, Nmode);
    Hm2.setZero(Nmode, Nmode);
    Hm3.setZero(Nmode, Nmode);
    Hm4.setZero(Nmode, Nmode);

    double tt = pow(t, 3) / 12, jcob, wf;
    ArrayXd s;
    MatrixXd F, dF, f, df;

    for (int i = 0; i < NE; i++)
    {
        jcob = Le(i) / 2;
        s = jcob * (GPs + 1);

        F = PSI(coef_s, coef_n, coef_z, i, s, 0);
        dF = PSI(coef_s, coef_n, coef_z, i, s, 1);

        for (int j = 0; j < Ngps; j++)
        {
            wf = WFs(j);
            f = F.block(0, (NMs + NMz) * j, 3, NMs + NMz);
            df = dF.block(0, (NMs + NMz) * j, 3, NMs + NMz);


            Hm1 += f.transpose() * f * jcob * wf * t;
            Hm2 += df.transpose() * Lm1.transpose() * Lm1 * df * jcob * wf * tt;
            Hm3 += df.transpose() * Lm1.transpose() * Lm2 * f * jcob * wf * tt;
            Hm4 += f.transpose() * Lm2.transpose() * Lm2 * f * jcob * wf * tt;
        }
    }

    vector<double> Lij;
    for (int i = 0; i < Li.size(); i++)
        Lij.push_back(Li(i));

    auto end = Lij.end();
    for (auto it = Lij.begin(); it != end; ++it)
        end = remove(it + 1, end, *it);

    Lij.erase(end, Lij.end());

    VectorXd Ltype, Const;
    Ltype.setZero(Nel);
    vector<MatrixXd> me;
    double Lm;
    double ksi;
    Hermite shN;
    vector<int> Lij_m;
    MatrixXd m1, m2, m3, m4;

    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
    for (int m = 0; m < Lij.size(); m++)
    {
        Lm = Lij[m];
        Lij_m.clear();
        for (int ii = 0; ii < Li.size(); ii++)
            if (Li(ii) == Lm)
                Lij_m.push_back(ii);

        Ltype(Lij_m) = Const.setConstant(Lij_m.size(), m);
        m1.setZero(Ndof * 2, Ndof * 2);
        m2.setZero(Ndof * 2, Ndof * 2);
        m3.setZero(Ndof * 2, Ndof * 2);
        m4.setZero(Ndof * 2, Ndof * 2);

        jcob = Lm / 2;
        for (int i = 0; i < Ngpz; i++)
        {
            ksi = GPz(i);
            wf = WFz(i);

            m1 += shN.N(Lm, ksi, Nmode).transpose() * Hm1 * shN.N(Lm, ksi, Nmode) * jcob * wf;
            m2 += shN.N(Lm, ksi, Nmode).transpose() * Hm2 * shN.N(Lm, ksi, Nmode) * jcob * wf;
            m3 += shN.N(Lm, ksi, Nmode).transpose() * Hm3 * shN.dN(Lm, ksi, Nmode) * jcob * wf;
            m4 += shN.dN(Lm, ksi, Nmode).transpose() * Hm4 * shN.dN(Lm, ksi, Nmode) * jcob * wf;
        }
        me.push_back(m1 + m2 + m3 + m3.transpose() + m4);
    }

    VectorXd Me(4 * Ndof * Ndof);
    MatrixXi edofs(2 * Ndof, Nel), Mi(2 * Ndof * 2 * Ndof, Nel), Mj(2 * Ndof * 2 * Ndof, Nel);
    MatrixXd MM((Nel + 1) * Ndof, (Nel + 1) * Ndof), Mall(2 * Ndof * 2 * Ndof, Nel);
    Me.setZero(), edofs = eldof.transpose(), Mi.setZero(), Mj.setZero(), Mall.setZero();
    vector<MatrixXd> Mme;
    MM.setZero();
    VectorXd en;
    for (int k = 0; k < Nel; k++)
    {
        for (int i = 0; i < 2 * Ndof;i++)
        {
            Me.block(2 * Ndof * i, 0, 2 * Ndof, 1) = me[Ltype(k)].block(0, i, 2 * Ndof, 1);
        }
        Mme.push_back(Me);
    }
    for (int i = 0; i < Nel; i++)
        Mall.block(0, i, 2 * Ndof * 2 * Ndof, 1) = Mme[i];
    for (int i = 0; i < 2 * Ndof; i++)
    {
        Mi.block(2 * Ndof * i, 0, 2 * Ndof, Nel) = edofs;
        for (int j = 0; j < 2 * Ndof; j++)
            Mj.block(2 * Ndof * i + j, 0, 1, Nel) = edofs.block(i, 0, 1, Nel);
    }

    //std::cout << "Mass Matrix Computed!" << endl << endl;

    return make_tuple(Mi, Mj, Mall);
}
