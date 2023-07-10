#pragma once

#include "HigherOrderBeamLib.h"
#include "LinearAlgebra.h"
#include "Stress.h"

class NodalValues
{
public:
    tuple<MatrixXd, MatrixXi, MatrixXd, MatrixXd> NodalDispAndStress(vector<VectorXd> D,vector<MatrixXd> cscd, vector<MatrixXi> cscc, vector<vector<MatrixXd>> coef_s, vector<vector<MatrixXd>> coef_n, vector<vector<MatrixXd>> coef_z, vector<VectorXd> zcoord, MatrixXd enr_s, RowVectorXi csnum_s, MatrixXd pnt_s, MatrixXi pnt_con_s, MatrixXd xdir_s)
    {
        int cs;
        double sm = 0.01;
        vector<MatrixXi> elnodeF;
        vector <MatrixXd> nodcoordF, noddispF;
        MatrixXi elnode, elnodei, Eni;
        MatrixXd nodcoord, nodcoordi, noddisp, noddispi, pci, Nci, Ndi, Sti, stress, stressbeam;
        Vector3d zdiri, xdiri;
        VectorXd scl;
        double scli, n;
        scl.setZero(pnt_con_s.rows());
        int els, ncs, nds, si;

        for (int i = 0; i < pnt_con_s.rows(); i++)
        {
            n = 0;
            cs = csnum_s(i);
            pci = pnt_s(pnt_con_s(i, all), all);
            zdiri = normalize(pci(1, all) - pci(0, all));
            xdiri = normalize(xdir_s(i, all));
            tie(elnodei, nodcoordi, noddispi, scli) = deformed(cscd[cs], cscc[cs], coef_s[cs], coef_n[cs], coef_z[cs], D[i], zcoord[i], sm);
            stressbeam = stress3D(cscd[cs], cscc[cs], coef_s[cs], coef_n[cs], coef_z[cs], n, D[i], zcoord[i], sm, enr_s(i, 0), enr_s(i, 1));
            Eni.setZero(elnodei.rows(), 5);
            Eni(all, 0) = VectorXi::Ones(elnodei.rows()) * i;
            Eni.block(0, 1, elnodei.rows(), elnodei.cols()) = elnodei + MatrixXi::Ones(elnodei.rows(), elnodei.cols()) * nodcoord.rows();
            Nci.setZero(nodcoordi.rows(), 3);
            for (int j = 0; j < nodcoordi.rows(); j++)
                Nci(j, all) = axistrans(nodcoordi(j, all), pci(0, all), zdiri, xdiri, 0).transpose();
            Ndi.setZero(noddispi.rows(), 3);
            for (int j = 0; j < noddispi.rows(); j++)
                Ndi(j, all) = axistrans(noddispi(j, all), Vector3d::Zero(), zdiri, xdiri, 0).transpose();
            Sti.setZero(stressbeam.rows(), 3);
            for (int j = 0; j < stressbeam.rows(); j++)
                Sti(j, all) = axistrans(stressbeam(j, all), Vector3d::Zero(), zdiri, xdiri, 0).transpose();
            scl(i) = scli;
            els = elnode.rows(), ncs = nodcoord.rows(), nds = noddisp.rows(), si = stress.rows();
            stress.conservativeResize(Sti.rows() + si, 3);
            stress.block(si, 0, Sti.rows(), Sti.cols()) = Sti;
            elnode.conservativeResize(els + Eni.rows(), 5);
            nodcoord.conservativeResize(ncs + Nci.rows(), 3);
            noddisp.conservativeResize(nds + Ndi.rows(), 3);
            elnode.block(els, 0, Eni.rows(), Eni.cols()) = Eni;
            nodcoord.block(ncs, 0, Nci.rows(), Nci.cols()) = Nci;
            noddisp.block(nds, 0, Ndi.rows(), Ndi.cols()) = Ndi;
        }

        cout << "Deformed Shape and Stress Calculated" << endl << endl;

        return make_tuple(nodcoord, elnode, noddisp, stress);
    }

    tuple<MatrixXd, MatrixXi, MatrixXd> ModalFreqDisp(vector<VectorXd> D,vector<MatrixXd> cscd, vector<MatrixXi> cscc, vector<vector<MatrixXd>> coef_s, vector<vector<MatrixXd>> coef_n, vector<vector<MatrixXd>> coef_z, vector<VectorXd> zcoord, MatrixXd enr_s, RowVectorXi csnum_s, MatrixXd pnt_s, MatrixXi pnt_con_s, MatrixXd xdir_s)
    {
        int cs;
        double sm = 0.01;
        vector<MatrixXi> elnodeF;
        vector <MatrixXd> nodcoordF, noddispF;
        MatrixXi elnode, elnodei, Eni;
        MatrixXd nodcoord, nodcoordi, noddisp, noddispi, pci, Nci, Ndi, stress, stressbeam;
        Vector3d zdiri, xdiri;
        VectorXd scl;
        double scli;
        int els, ncs, nds, si;

        scl.setZero(pnt_con_s.rows());
        for (int i = 0; i < pnt_con_s.rows(); i++)
        {
            //cout << "i: " << i << endl << endl;
            elnodei.setZero(), nodcoordi.setZero(), noddispi.setZero();
            cs = csnum_s(i);
            pci = pnt_s(pnt_con_s(i, all), all);
            zdiri = normalize(pci(1, all) - pci(0, all));
            xdiri = normalize(xdir_s(i, all));
            tie(elnodei, nodcoordi, noddispi, scli) = deformed(cscd[cs], cscc[cs], coef_s[cs], coef_n[cs], coef_z[cs], D[i], zcoord[i], sm);
            Eni.setZero(elnodei.rows(), 5);
            Eni(all, 0) = VectorXi::Ones(elnodei.rows()) * i;
            Eni.block(0, 1, elnodei.rows(), elnodei.cols()) = elnodei + MatrixXi::Ones(elnodei.rows(), elnodei.cols()) * nodcoord.rows();
            Nci.setZero(nodcoordi.rows(), 3);
            for (int k = 0; k < nodcoordi.rows(); k++)
                Nci(k, all) = axistrans(nodcoordi(k, all), pci(0, all), zdiri, xdiri, 0).transpose();
            Ndi.setZero(noddispi.rows(), 3);
            for (int k = 0; k < noddispi.rows(); k++)
                Ndi(k, all) = axistrans(noddispi(k, all), Vector3d::Zero(), zdiri, xdiri, 0).transpose();
            scl(i) = scli;
            els = elnode.rows(), ncs = nodcoord.rows(), nds = noddisp.rows();
            elnode.conservativeResize(els + Eni.rows(), 5);
            nodcoord.conservativeResize(ncs + Nci.rows(), 3);
            noddisp.conservativeResize(nds + Ndi.rows(), 3);
            elnode.block(els, 0, Eni.rows(), Eni.cols()) = Eni;
            nodcoord.block(ncs, 0, Nci.rows(), Nci.cols()) = Nci;
            noddisp.block(nds, 0, Ndi.rows(), Ndi.cols()) = Ndi;
        }

        return make_tuple(nodcoord, elnode, noddisp);
    }

    tuple<MatrixXd, MatrixXd, MatrixXd> ShellDisp(MatrixXi ShellNdof, MatrixXd Ds, MatrixXd Shetable, vector<MatrixXd> shellcoord, MatrixXd scoord)
    {
        MatrixXd disp, stress, coords;

        disp = scoord(all, { 1,2,3 }) + Ds(all, { 0,1,2 });

        /*std::cout << "Hello waiter" << std::endl;
        std::chrono::seconds dura(300);
        std::this_thread::sleep_for( dura );
        std::cout << "Waited 5s\n";*/

        return make_tuple(disp, stress, coords);
    }
};
