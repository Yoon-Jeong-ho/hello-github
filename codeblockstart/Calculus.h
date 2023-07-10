#pragma once

#include "Order.h"

class Calculus
{
public:
    vector<Matrix<double50, -1, -1>> Intcoef(vector<Matrix<double50, -1, -1>> coef)
    {
        int NE = coef.size();
        Index NM = coef[0].cols(), NC = coef[0].rows();
        Array<double50, -1, 1> sp, spo;
        sp.setOnes(NC, 1);
        for (int i = 1;i < NC;i++)
        {
            sp(i, 0) = sp(i, 0) / i;
        }
        Matrix<double50, -1, -1> sp1(NC, NM);
        for (int i = 0;i < NM;i++)
        {
            sp1(all, i) = sp;
        }

        Matrix<double50, -1, -1> Icoef;
        Icoef.setZero(NC, NM);

        for (int i = 0; i < NE; i++)
        {
            Icoef.block(1, 0, NC - 1, NM) = coef[i].block(0, 0, NC - 1, NM).cwiseProduct(sp1.block(1, 0, NC - 1, NM));
            coef[i] = Icoef;
        }

        return coef;
    }
    vector<Matrix<double50, -1, -1>> Intcoef(vector<Matrix<double50, -1, -1>> coef, int n)
    {
        /*int NE = coef.size();
        Index NM = coef[0].cols(), NC = coef[0].rows();
        Order order;
        MatrixXi mo = order.ordermp(coef);
        for (int i = 0; i < mo.size(); i++)
            mo(i) = max(mo(i), 1);

        cout << mo << endl << endl;
        cout << endl;
        Matrix<double50, -1, 1> sp, spo;
        sp.setOnes(NC, 1);
        for (int i = 1;i < NC + 1;i++)
        {
            sp(i - 1, 0) = sp(i - 1, 0) / i;
        }
        Matrix<double50, -1, -1> Icoef1 = MatrixCast(coef, NC, NM * NE);
        Matrix<double50, -1, -1> Icoef2, Icoef3;
        Icoef2.setZero(NC, 1), Icoef3.setZero(NC, static_cast<Eigen::Index>(NM) * NE);
        int mm;
        for (int m = 0; m < NM; m++)
        {
            mm = mo(m);
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < NE; i++)
                {
                    Icoef2.block(1, 0, NC - 1, 1) = (mm + j + 1) * Icoef1.block(0, static_cast<Eigen::Index>(i) * NM + m, NC - 1, 1).cwiseProduct(sp.block(0, 0, NC - 1, 1));
                    Icoef1.block(0, static_cast<Eigen::Index>(i) * NM + m, NC, 1) = Icoef2;
                }
            }
        }

        Tensor<double50, 3> Icoef = TensorCast(Icoef1, NC, NM, NE); // Cast Eigen::Matrix --> Eigen::Tensor*/

        return coef;
    }
    vector<Matrix<double50, -1, -1>> diffcoefmp(vector<Matrix<double50, -1, -1>> coef, int d)
    {
        vector<Matrix<double50, -1, -1>> Dcoef;
        if (d > 0)
        {
            int NE = coef.size();
            Index NM = coef[0].cols();
            Index NC = coef[0].rows();
            Array<double50, -1, 1> spow;
            spow.setLinSpaced(NC, 0, NC - 1);
            Matrix<double50, -1, -1> sp(NC, NM);
            for (int i = 0;i < NM;i++)
            {
                sp(all, i) = spow;
            }

            vector<Matrix<double50, -1, -1>> dcoef1 = coef, dcoef3 = coef;
            Matrix<double50, -1, -1> dcoef2;
            dcoef2.setZero(NC, NM);

            for (int i = 0; i < d;i++)
            {
                for (int j = 0; j < NE;j++)
                {
                    dcoef2.block(0, 0, NC - 1, NM) = dcoef1[j].block(1, 0, NC - 1, NM).cwiseProduct(sp.block(1, 0, NC - 1, NM));
                    dcoef1[j].block(0, 0, NC, NM) = dcoef2;
                    dcoef3[j] = dcoef2; // Check this loop
                }
            }
            Dcoef = dcoef3;
        }
        else
            Dcoef = coef;

        return Dcoef;
    }
    vector<MatrixXd> diffcoef(vector<MatrixXd> coef, int d)
    {
        vector<MatrixXd> Dcoef;
        if (d > 0)
        {
            int NE = coef.size();
            Index NM = coef[0].cols();
            Index NC = coef[0].rows();
            ArrayXd spow;
            spow.setLinSpaced(NC, 0, NC - 1);
            MatrixXd sp(NC, NM);
            for (int i = 0;i < NM;i++)
            {
                sp(all, i) = spow;
            }

            vector<MatrixXd> dcoef1 = coef, dcoef3 = coef;
            MatrixXd dcoef2;
            dcoef2.setZero(NC, NM);

            for (int i = 0; i < d;i++)
            {
                for (int j = 0; j < NE;j++)
                {
                    dcoef2.block(0, 0, NC - 1, NM) = dcoef1[j].block(1, 0, NC - 1, NM).cwiseProduct(sp.block(1, 0, NC - 1, NM));
                    dcoef1[j].block(0, 0, NC, NM) = dcoef2;
                    dcoef3[j] = dcoef2; // Check this loop
                }
            }
            Dcoef = dcoef3;
        }
        else
            Dcoef = coef;

        return Dcoef;
    }
};