#pragma once

#include "HigherOrderBeamLib.h"

class Order
{
public:
    MatrixXi ordermp(vector<Matrix<double50, -1, -1>> coef)
    {
        int NE = coef.size();
        Index NM = coef[0].cols(), NC = coef[0].rows();
        Array<double50, -1, 1> f;
        MatrixXi ord(NE, NM);
        vector<int> ord_ij;
        ord.setZero();
        for (int i = 0; i < NE; i++)
        {
            for (int j = 0; j < NM; j++)
            {
                ord_ij.clear();
                f = coef[i].block(0, j, NC, 1);
                for (int k = 0; k < f.rows(); k++)
                {
                    if (f(k) != 0)
                        ord_ij.push_back(k);
                }
                if (ord_ij.size() > 0)
                    ord(i, j) = ord_ij[ord_ij.size() - 1];
                else
                    ord(i, j) = 0;
            }
        }

        return ord;
    }
    MatrixXi order(vector<MatrixXd> coef)
    {
        int NE = coef.size();
        Index NM = coef[0].cols(), NC = coef[0].rows();
        ArrayXd f;
        MatrixXi ord(NE, NM);
        vector<int> ord_ij;
        ord.setZero();
        for (int i = 0; i < NE; i++)
        {
            for (int j = 0; j < NM; j++)
            {
                ord_ij.clear();
                f = coef[i].block(0, j, NC, 1);
                for (int k = 0; k < f.rows(); k++)
                {
                    if (f(k) != 0)
                        ord_ij.push_back(k);
                }
                if (ord_ij.size() > 0)
                    ord(i, j) = ord_ij[ord_ij.size() - 1];
                else
                    ord(i, j) = 0;
            }
        }

        return ord;
    }
};