#pragma once

#include "HigherOrderBeamLib.h"
#include "Edge.h"

class Continuity
{
public:
    tuple<vector<vector<int>>, vector<Matrix<double50, -1, -1>>, vector<Array<double50, -1, 1>>> Edges(Matrix<double50, -1, -1>  csc, MatrixXi cscc, vector<Matrix<double50, -1, -1>> coef, int t)
    {
        Index NC = csc.rows(), NM = coef[0].cols();
        Array<double50, -1, 1> alpha, Le, sc;
        VectorXi rng_r(2), edgv, rng_s;
        Edge edge;
        alpha = edge.Anglemp(csc, cscc);
        Le = edge.Lengthmp(csc, cscc);
        Matrix<double50, -1, -1> cont_cond{}, Z(Le.rows(), 1), f1, f2, f3, C, T(2, 2), Rxyf(2, NM), Rxy, Uxyf(2, NM), Uxy, temp1(1, 2), c0, c1, c2;
        vector<vector<int>> cont_edg;
        vector<Matrix<double50, -1, -1>> cont_Rxy;
        vector<Array<double50, -1, 1>> cont_sc;
        Z.setZero();
        vector<int> edg;
        vector<double50> ss;
        int count1 = 0, count2 = 0, count4;

        for (int i = 0; i < NC;i++)
        {
            edg.clear();
            ss.clear();
            for (int j = 0;j < 2;j++)
            {
                for (int k = 0; k < cscc.rows();k++)
                {
                    if (cscc(k, j) == i)
                    {
                        edg.push_back(k);
                    }
                }
                if (j == 0)
                {
                    for (int k = 0; k < edg.size();k++)
                    {
                        ss.push_back(Z(edg[k]));
                    }
                    count4 = ss.size();
                }
                else if (j == 1)
                    for (int k = count4; k < edg.size();k++)
                    {
                        ss.push_back(Le(edg[k]));
                    }
            }
            sc.setZero(ss.size());
            for (int k = 0; k < ss.size(); k++)
                sc(k) = ss[k];
            cont_edg.push_back(edg);
            cont_sc.push_back(sc);
            if (t == 0) // OUT OF PLANE MODES
            {
                if (edg.size() > 1)
                {
                    count1 += edg.size() - 1;
                    cont_cond.conservativeResize(count1, NM);
                    for (int j = 0;j < edg.size() - 1;j++)
                    {
                        f1 = edge.Continuity(coef, edg[j], sc(j), 0);
                        f2 = edge.Continuity(coef, edg[static_cast<std::vector<int, std::allocator<int>>::size_type>(j) + 1], sc(static_cast<Eigen::Index>(j) + 1), 0);
                        cont_cond(count2, all) = f1 - f2;
                        count2++;
                    }
                }
                // Corner Angle Rotation
                if (edg.size() > 2)
                {
                    if (abs(sin(alpha(edg[0]) - alpha(edg[1]))) < 1E-10)
                        rng_r << 0, 2;
                    else
                        rng_r << 0, 1;

                    f1 = edge.Continuity(coef, edg[rng_r(0)], sc(rng_r(0)), 1);
                    f2 = edge.Continuity(coef, edg[rng_r(1)], sc(rng_r(1)), 1);
                    T << cos(alpha(edg[rng_r(1)])), -cos(alpha(edg[rng_r(0)])),
                        sin(alpha(edg[rng_r(1)])), -sin(alpha(edg[rng_r(0)]));
                    T = T / (sin(alpha(edg[rng_r(0)]) - alpha(edg[rng_r(1)])));
                    Rxyf(0, all) = f1, Rxyf(1, all) = f2;
                    Rxy = T * Rxyf;

                    cont_Rxy.push_back(Rxy);
                }
                else if (edg.size() == 2)
                {
                    f1 = edge.Continuity(coef, edg[0], sc(0), 1);
                    f2 = edge.Continuity(coef, edg[1], sc(1), 1);
                    T(0, 0) = cos(alpha(edg[1]));
                    T(0, 1) = -cos(alpha(edg[0]));
                    T(1, 0) = sin(alpha(edg[1]));
                    T(1, 1) = -sin(alpha(edg[0]));
                    T = T / (sin(alpha(edg[0]) - alpha(edg[1])));
                    Rxyf(0, all) = f1, Rxyf(1, all) = f2;
                    Rxy = T * Rxyf;


                    cont_Rxy.push_back(Rxy);
                }
            }
            else if (t == 1) // IN-PLANE MODES
            {
                if (edg.size() > 2)
                {
                    if (abs(sin(alpha(edg[0]) - alpha(edg[1]))) < 1E-10)
                        rng_r << 0, 2;
                    else
                        rng_r << 0, 1;

                    edgv.setZero(edg.size());
                    rng_s.setZero(edgv.size());

                    for (int j = 0; j < edg.size();j++)
                        edgv(j) = j;

                    sort(edgv.begin(), edgv.end());
                    VectorXi::iterator it;

                    it = set_difference(edgv.begin(), edgv.end(), rng_r.begin(), rng_r.end(), rng_s.begin());
                    rng_s.conservativeResize(it - rng_s.begin());

                    f1 = edge.Continuity(coef, edg[rng_r(0)], sc(rng_r(0)), 0);
                    f2 = edge.Continuity(coef, edg[rng_r(1)], sc(rng_r(1)), 0);
                    T(0, 0) = -sin(alpha(edg[rng_r(1)]));
                    T(0, 1) = sin(alpha(edg[rng_r(0)]));
                    T(1, 0) = cos(alpha(edg[rng_r(1)]));
                    T(1, 1) = -cos(alpha(edg[rng_r(0)]));
                    T = T / (sin(alpha(edg[rng_r(0)]) - alpha(edg[rng_r(1)])));
                    Uxyf(0, all) = f1, Uxyf(1, all) = f2;
                    Uxy = T * Uxyf;

                    cont_Rxy.push_back(Uxy);
                    count1 += rng_s.size();
                    cont_cond.conservativeResize(count1, NM);

                    for (int j = 0;j < rng_s.size(); j++)
                    {
                        f3 = edge.Continuity(coef, edg[rng_s(j)], sc(rng_s(j)), 0);
                        temp1 << cos(alpha(edg[rng_s(j)])), sin(alpha(edg[rng_s(j)]));
                        cont_cond(count2, all) = f3 - temp1 * Uxy;
                        count2++;
                    }
                }
                else if (edg.size() == 2)
                {
                    f1 = edge.Continuity(coef, edg[0], sc(0), 0);
                    f2 = edge.Continuity(coef, edg[1], sc(1), 0);
                    T(0, 0) = -sin(alpha(edg[1]));
                    T(0, 1) = sin(alpha(edg[0]));
                    T(1, 0) = cos(alpha(edg[1]));
                    T(1, 1) = -cos(alpha(edg[0]));
                    T = T / (sin(alpha(edg[0]) - alpha(edg[1])));
                    Uxyf(0, all) = f1, Uxyf(1, all) = f2;
                    Uxy = T * Uxyf;

                    cont_Rxy.push_back(Uxy);
                }
            }
            else if (t == 2)// WALL BENDING MODES
            {
                if (edg.size() == 1)
                {
                    C = edge.Continuity(coef, edg[0], sc(0), 2);
                }
                else
                {
                    C.setZero(2 * edg.size(), NM);
                    c0.setZero(edg.size(), NM);
                    for (int j = 0;j < edg.size();j++)
                    {
                        c0(j, all) = edge.Continuity(coef, edg[j], sc(j), 0);
                    }

                    c1.setZero(edg.size() - 1, NM);
                    for (int j = 0;j < edg.size() - 1;j++)
                    {
                        for (int k = 0;k < 2;k++)
                            c1(j, all) += edge.Continuity(coef, edg[j + k], sc(j + k), 1) * pow(-1, k);
                    }

                    c2.setZero(1, NM);
                    for (int j = 0;j < edg.size();j++)
                    {
                        if (sc(j) == 0)
                            c2 += edge.Continuity(coef, edg[j], sc[j], 2) * (-1);
                        else
                            c2 += edge.Continuity(coef, edg[j], sc[j], 2);
                    }
                    C.block(0, 0, edg.size(), NM) = c0;
                    C.block(edg.size(), 0, edg.size() - 1, NM) = c1;
                    C.block(edg.size() + edg.size() - 1, 0, 1, NM) = c2;
                }
                count1 += C.rows();
                cont_cond.conservativeResize(count1, NM);
                cont_cond.block(count1 - C.rows(), 0, C.rows(), NM) = C;
            }
        }

        return make_tuple(cont_edg, cont_Rxy, cont_sc);
    }

    Matrix<double50, -1, -1> Constraint(Matrix<double50, -1, -1>  csc, MatrixXi cscc, vector<Matrix<double50, -1, -1>> coef, int t)
    {
        Index NC = csc.rows(), NM = coef[0].cols();
        Array<double50, -1, 1> alpha, Le, sc;
        VectorXi rng_r(2), edgv, rng_s;
        Edge edge;
        alpha = edge.Anglemp(csc, cscc);
        Le = edge.Lengthmp(csc, cscc);
        Matrix<double50, -1, -1> cont_cond{}, Z(Le.rows(), 1), f1, f2, f3, C, T(2, 2), Rxyf(2, NM), Rxy, Uxyf(2, NM), Uxy, temp1(1, 2), c0, c1, c2;
        vector<vector<int>> cont_edg;
        vector<Matrix<double50, -1, -1>> cont_Rxy;
        vector<Array<double50, -1, 1>> cont_sc;
        Z.setZero();
        vector<int> edg;
        vector<double50> ss;
        int count1 = 0, count2 = 0, count4;

        for (int i = 0; i < NC;i++)
        {
            edg.clear();
            ss.clear();
            for (int j = 0;j < 2;j++)
            {
                for (int k = 0; k < cscc.rows();k++)
                {
                    if (cscc(k, j) == i)
                    {
                        edg.push_back(k);
                    }
                }
                if (j == 0)
                {
                    for (int k = 0; k < edg.size();k++)
                    {
                        ss.push_back(Z(edg[k]));
                    }
                    count4 = ss.size();
                }
                else if (j == 1)
                    for (int k = count4; k < edg.size();k++)
                    {
                        ss.push_back(Le(edg[k]));
                    }
            }
            sc.setZero(ss.size());
            for (int k = 0; k < ss.size(); k++)
                sc(k) = ss[k];
            cont_edg.push_back(edg);
            cont_sc.push_back(sc);
            if (t == 0) // OUT OF PLANE MODES
            {
                if (edg.size() > 1)
                {
                    count1 += edg.size() - 1;
                    cont_cond.conservativeResize(count1, NM);
                    for (int j = 0;j < edg.size() - 1;j++)
                    {
                        f1 = edge.Continuity(coef, edg[j], sc(j), 0);
                        f2 = edge.Continuity(coef, edg[static_cast<std::vector<int, std::allocator<int>>::size_type>(j) + 1], sc(static_cast<Eigen::Index>(j) + 1), 0);
                        cont_cond(count2, all) = f1 - f2;
                        count2++;
                    }
                }
                // Corner Angle Rotation
                if (edg.size() > 2)
                {
                    if (abs(sin(alpha(edg[0]) - alpha(edg[1]))) < 1E-10)
                        rng_r << 0, 2;
                    else
                        rng_r << 0, 1;

                    f1 = edge.Continuity(coef, edg[rng_r(0)], sc(rng_r(0)), 1);
                    f2 = edge.Continuity(coef, edg[rng_r(1)], sc(rng_r(1)), 1);
                    T << cos(alpha(edg[rng_r(1)])), -cos(alpha(edg[rng_r(0)])),
                        sin(alpha(edg[rng_r(1)])), -sin(alpha(edg[rng_r(0)]));
                    T = T / (sin(alpha(edg[rng_r(0)]) - alpha(edg[rng_r(1)])));
                    Rxyf(0, all) = f1, Rxyf(1, all) = f2;
                    Rxy = T * Rxyf;

                    cont_Rxy.push_back(Rxy);
                }
                else if (edg.size() == 2)
                {
                    f1 = edge.Continuity(coef, edg[0], sc(0), 1);
                    f2 = edge.Continuity(coef, edg[1], sc(1), 1);
                    T(0, 0) = cos(alpha(edg[1]));
                    T(0, 1) = -cos(alpha(edg[0]));
                    T(1, 0) = sin(alpha(edg[1]));
                    T(1, 1) = -sin(alpha(edg[0]));
                    T = T / (sin(alpha(edg[0]) - alpha(edg[1])));
                    Rxyf(0, all) = f1, Rxyf(1, all) = f2;
                    Rxy = T * Rxyf;

                    cont_Rxy.push_back(Rxy);
                }
            }
            else if (t == 1) // IN-PLANE MODES
            {
                if (edg.size() > 2)
                {
                    if (abs(sin(alpha(edg[0]) - alpha(edg[1]))) < 1E-10)
                        rng_r << 0, 2;
                    else
                        rng_r << 0, 1;

                    edgv.setZero(edg.size());
                    rng_s.setZero(edgv.size());

                    for (int j = 0; j < edg.size();j++)
                        edgv(j) = j;

                    sort(edgv.begin(), edgv.end());
                    VectorXi::iterator it;

                    it = set_difference(edgv.begin(), edgv.end(), rng_r.begin(), rng_r.end(), rng_s.begin());
                    rng_s.conservativeResize(it - rng_s.begin());

                    f1 = edge.Continuity(coef, edg[rng_r(0)], sc(rng_r(0)), 0);
                    f2 = edge.Continuity(coef, edg[rng_r(1)], sc(rng_r(1)), 0);
                    T(0, 0) = -sin(alpha(edg[rng_r(1)]));
                    T(0, 1) = sin(alpha(edg[rng_r(0)]));
                    T(1, 0) = cos(alpha(edg[rng_r(1)]));
                    T(1, 1) = -cos(alpha(edg[rng_r(0)]));
                    T = T / (sin(alpha(edg[rng_r(0)]) - alpha(edg[rng_r(1)])));
                    Uxyf(0, all) = f1, Uxyf(1, all) = f2;
                    Uxy = T * Uxyf;

                    cont_Rxy.push_back(Uxy);
                    count1 += rng_s.size();
                    cont_cond.conservativeResize(count1, NM);

                    for (int j = 0;j < rng_s.size(); j++)
                    {
                        f3 = edge.Continuity(coef, edg[rng_s(j)], sc(rng_s(j)), 0);
                        temp1 << cos(alpha(edg[rng_s(j)])), sin(alpha(edg[rng_s(j)]));
                        cont_cond(count2, all) = f3 - temp1 * Uxy;
                        count2++;
                    }
                }
                else if (edg.size() == 2)
                {
                    f1 = edge.Continuity(coef, edg[0], sc(0), 0);
                    f2 = edge.Continuity(coef, edg[1], sc(1), 0);
                    T(0, 0) = -sin(alpha(edg[1]));
                    T(0, 1) = sin(alpha(edg[0]));
                    T(1, 0) = cos(alpha(edg[1]));
                    T(1, 1) = -cos(alpha(edg[0]));
                    T = T / (sin(alpha(edg[0]) - alpha(edg[1])));
                    Uxyf(0, all) = f1, Uxyf(1, all) = f2;
                    Uxy = T * Uxyf;

                    cont_Rxy.push_back(Uxy);
                }
            }
            else if (t == 2)// WALL BENDING MODES
            {
                if (edg.size() == 1)
                {
                    C = edge.Continuity(coef, edg[0], sc(0), 2);
                }
                else
                {
                    C.setZero(2 * edg.size(), NM);
                    c0.setZero(edg.size(), NM);
                    for (int j = 0;j < edg.size();j++)
                    {
                        c0(j, all) = edge.Continuity(coef, edg[j], sc(j), 0);
                    }

                    c1.setZero(edg.size() - 1, NM);
                    for (int j = 0;j < edg.size() - 1;j++)
                    {
                        for (int k = 0;k < 2;k++)
                            c1(j, all) += edge.Continuity(coef, edg[j + k], sc(j + k), 1) * pow(-1, k);
                    }

                    c2.setZero(1, NM);
                    for (int j = 0;j < edg.size();j++)
                    {
                        if (sc(j) == 0)
                            c2 += edge.Continuity(coef, edg[j], sc[j], 2) * (-1);
                        else
                            c2 += edge.Continuity(coef, edg[j], sc[j], 2);
                    }
                    C.block(0, 0, edg.size(), NM) = c0;
                    C.block(edg.size(), 0, edg.size() - 1, NM) = c1;
                    C.block(edg.size() + edg.size() - 1, 0, 1, NM) = c2;
                }
                count1 += C.rows();
                cont_cond.conservativeResize(count1, NM);
                cont_cond.block(count1 - C.rows(), 0, C.rows(), NM) = C;
            }
        }

        return cont_cond;
    }
};

