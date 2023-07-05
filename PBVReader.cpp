// PBVReader.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <iomanip>
#include <thread>
#include <chrono>
#include "Center.h"
#include "Edge.h"

using namespace std;
using namespace Eigen;
using namespace boost::multiprecision;
typedef cpp_dec_float_50 double50;

class Extract
{
public:
    int CrossSection(string str)
    {
        stringstream ss;
        string temp;
        ss << str;
        int found = 0;
        while (found == 0) {
            ss >> temp;
            if (temp.length() > 0){
            if (stringstream(temp) >> found) {}
            temp = "";
        }}
        return found;
    }

    RowVector2d EdgeLandT(string str)
    {
        stringstream ss;
        string temp;
        /* Storing the whole string into string stream */
        ss << str;

        /* Running loop till the end of the stream */
        RowVector2d found; found.setZero();
        int si = 0;
        while (!ss.eof()) {

            /* extracting word by word from stream */
            ss >> temp;

            /* Checking the given word is integer or not */
            if (temp.length() > 0){
            if (stringstream(temp) >> found(si)) {
                //cout << found << " ";
                si++;
            }}
            /* To save from space at the end of string */
            temp = "";
        }

        return found;
    }

    RowVector4d Nodes(string str)
    {
        stringstream ss;
        string temp;
        ss << str;
        RowVector4d found; found.setZero();
        int si = 0;
        while (!ss.eof()) {
            ss >> temp;
            if (temp.length() > 0){
            if (stringstream(temp) >> found(si)) {
                si++;
            }}
            temp = "";
        }
        return found;
    }

    RowVector3i CSMandM(string str)
    {
        stringstream ss;
        string temp;
        ss << str;
        RowVector3i found; found.setZero();
        int si = 0;
        while (!ss.eof()) {
            ss >> temp;
            if (temp.length() > 0){
            if (stringstream(temp) >> found(si)) {
                si++;
            }}
            temp = "";
        }
        return found;
    }

    RowVector3d BCs(string str)
    {
        stringstream ss;
        string temp;
        ss << str;
        RowVector3d found; found.setZero();
        int si = 0;
        while (!ss.eof()) {
            ss >> temp;
            if (temp.length() > 0){
            if (stringstream(temp) >> found(si)) {
                si++;
            }}
            temp = "";
        }
        return found;
    }

    RowVector4d Material(string str)
    {
        stringstream ss;
        string temp;
        ss << str;
        RowVector4d found; found.setZero();
        int si = 0;
        while (!ss.eof()) {
            ss >> temp;
            if (temp.length() > 0){
            if (stringstream(temp) >> found(si)) {
                si++;
            }}
            temp = "";
        }
        return found;
    }

    RowVectorXd UserCross(string str)
    {
        stringstream ss;
        string temp;
        ss << str;
        RowVectorXd found; found.setZero(6);
        int si = 0;
        while (!ss.eof()) {
            ss >> temp;
            if (temp.length() > 0){
            if (stringstream(temp) >> found(si)) {
                si++;
            }}
            temp = "";
        }
        return found;
    }

    RowVectorXd Elements(string str)
    {
        stringstream ss;
        string temp;
        ss << str;
        RowVectorXd found; found.setZero(7);
        int si = 0;
        while (!ss.eof()) {
            ss >> temp;
            if (temp.length() > 0){
            if (stringstream(temp) >> found(si)) {
                si++;
            }}
            temp = "";
        }
        return found;
    }

    RowVectorXd Joint(string str)
    {
        stringstream ss;
        string temp;
        ss << str;
        RowVectorXd found; found.setZero(6);
        int si = 0;
        while (!ss.eof()) {
            ss >> temp;
            if (temp.length() > 0){
            if (stringstream(temp) >> found(si)) {
                si++;
            }}
            temp = "";
        }
        return found;
    }

    RowVectorXd Forces(string str, int n)
    {
        stringstream ss;
        string temp;
        ss << str;
        RowVectorXd found; found.setZero(n);
        int si = 0;
        while (!ss.eof()) {
            ss >> temp;
            if (temp.length() > 0){
            if (stringstream(temp) >> found(si)) {
                si++;
            }}
            temp = "";
        }
        return found;
    }
};

class Cross_Section
{
public:
    Vector<double50, -1> EdgeLength(Matrix<double50, -1, -1> csc, MatrixXi cscc)
    {
        Vector<double50, -1> Le;
        VectorXi e;
        Le.setZero(cscc.rows());
        for (int i = 0; i < cscc.rows(); i++)
        {
            e = cscc(i, all);
            Le(i) = sqrt(pow((csc(e(1), 1) - csc(e(0), 1)), 2) + pow((csc(e(1), 0) - csc(e(0), 0)), 2));
        }
        return Le;
    }
    Vector<double50, -1> EdgeAngle(Matrix<double50, -1, -1> csc, MatrixXi cscc)
    {
        Vector<double50, -1> alpha;
        VectorXi e;
        alpha.setZero(cscc.rows());
        for (int i = 0; i < cscc.rows(); i++)
        {
            e = cscc(i, all);
            alpha(i) = atan2(csc(e(1), 1) - csc(e(0), 1), csc(e(1), 0) - csc(e(0), 0));
        }
        return alpha;
    }
};

void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

    matrix.conservativeResize(numRows, numCols);
}

tuple< vector<Matrix<double50, -1, -1>>, vector<MatrixXi>, RowVectorXi, MatrixXd, MatrixXd, MatrixXi, MatrixXd, VectorXd, MatrixXi, MatrixXd, MatrixXd, MatrixXd, MatrixXd, vector<VectorXd>, VectorXi, VectorXi, string, int, ArrayXi, RowVector3d, string, MatrixXd, MatrixXd, MatrixXd, MatrixXd, vector<MatrixXd>, MatrixXd> inputReader(string txtFile)
{
    ifstream inFile;
    inFile.open(txtFile);

    if (!inFile) {
        cerr << "Unable to open file datafile.txt";
        exit(1);   // call system to stop
    }

    //cout << inFile.rdbuf() << endl << endl;

    string tp, analysis;
    int si = 0, csnumi = 0, freq = 0;
    vector<int> sec, sectype, bcs, fcs, nw1;
    vector<RowVector<double50, -1>> dim;
    vector<RowVectorXd> Elements, Joints, UCross1, BoundCs, Forces, PShell, ShellElements;
    vector<RowVector3i> CSMandM;
    vector<RowVector2d> shelljoint;
    vector<Matrix<double50, -1, -1>> UserCross;
    Matrix<double50, -1, -1> UCross2;
    RowVector2d csdim;
    RowVector3d csdimch, FRF;
    RowVectorXd jointcheck;
    Extract dimensions;
    double csthick;
    vector<double> usthick;
    RowVector<double50, -1> csType(5), csTypei(6);
    vector<RowVector4d> Nodes, Materials;
    MatrixXd HO_Nodes, HO_Elements, HO_Joints, HO_Materials, HO_BCs, FEdge, FNode, HO_PShell, Shell, ShellJoints;
    MatrixXi ModenMat;
    ArrayXi NW1;

    while (getline(inFile, tp))
    {
        if (tp.substr(0,11) == "*HO_SECTION")
        {
            //sec.push_back(si); // Number of Cross Sections and their places in the test file
            getline(inFile, tp); si++;
            getline(inFile, tp); si++;
            CSMandM.push_back(dimensions.CSMandM(tp));
            if (tp.find("OPEN") != string::npos)
            {
                sectype.push_back(0); // 0 means open Cross Section
                getline(inFile, tp); si++;
                csdim(0) = dimensions.EdgeLandT(tp)(0); // Extracting Edge Length and Thickness
                getline(inFile, tp); si++;
                csdim(1) = dimensions.EdgeLandT(tp)(0);
                csthick = dimensions.EdgeLandT(tp)(1);
                csType << csnumi + 1, sectype[csnumi], csdim, csthick; // 결정된 모양에 따른 숫자 넘겨주기
                dim.push_back(csType);
                nw1.push_back(2);
            }
            else if (tp.find("RECT") != string::npos)
            {
                sectype.push_back(1); // 1 means closed cross section
                getline(inFile, tp); si++;
                csdim(0) = dimensions.EdgeLandT(tp)(0); // Extracting Edge Length and Thickness
                getline(inFile, tp); si++;
                csdim(1) = dimensions.EdgeLandT(tp)(0);
                csthick = dimensions.EdgeLandT(tp)(1);
                csType << csnumi + 1, sectype[csnumi], csdim, csthick;
                dim.push_back(csType);
                nw1.push_back(1);
            }
            else if (tp.find("USER") != string::npos)
            {
                sectype.push_back(2); // 2 means User Defined cross section
                UCross1.clear();
                getline(inFile, tp); si++;
                while (tp.length() > 1)
                {
                    UCross1.push_back(dimensions.UserCross(tp));
                    getline(inFile, tp); si++;
                }
                UCross2.setZero(UCross1.size(), 4);
                for (int i = 0; i < UCross1.size(); i++)
                    UCross2(i, all) = UCross1[i]({ 0,1,2,3 });
                UserCross.push_back(UCross2);
                usthick.push_back(UCross1[0](4));
                nw1.push_back(1);
            }
            else if (tp.find("CHAMFER") != string::npos)
            {
                sectype.push_back(3); // 3 means CHAMFERED cross section
                getline(inFile, tp); si++;
                csdimch(0) = dimensions.EdgeLandT(tp)(0); // Extracting Edge Length and Thickness
                getline(inFile, tp); si++;
                csdimch(1) = dimensions.EdgeLandT(tp)(0);
                getline(inFile, tp); si++;
                csdimch(2) = dimensions.EdgeLandT(tp)(0);
                csthick = dimensions.EdgeLandT(tp)(1);
                csTypei << csnumi + 1, sectype[csnumi], csdimch, csthick;
                dim.push_back(csTypei);
                nw1.push_back(1);
            }
            csnumi++;
        }
        else if (tp.substr(0,12) == "*HO_ANALYSIS")
        {
            getline(inFile, tp); si++;
            analysis = tp;
        }
        else if (tp.substr(0,12) == "*HO_MATERIAL")
        {
            getline(inFile, tp); si++;
            getline(inFile, tp); si++;
            // Length should be checked with > 1 not != 1
            //while (tp.length() != 1)
            while (tp.length() > 1)
            {
                Materials.push_back(dimensions.Material(tp));
                getline(inFile, tp); si++;
            }
        }
        else if (tp.substr(0,10) == "*HO_PSHELL")
        {
            getline(inFile, tp); si++;
            getline(inFile, tp); si++;
            while (tp.empty() == 0)
            {
                PShell.push_back(dimensions.BCs(tp));
                getline(inFile, tp); si++;
            }

        }
        else if (tp.substr(0,8) == "*HO_NODE")
        {
            getline(inFile, tp); si++;
            while (tp.length() > 1)
            {
                Nodes.push_back(dimensions.Nodes(tp));
                getline(inFile, tp); si++;
            }
        }
        else if (tp.substr(0,8) == "*HO_BEAM")
        {
            getline(inFile, tp); si++;
            while (tp.length() > 1)
            {
                Elements.push_back(dimensions.Elements(tp));
                getline(inFile, tp); si++;
            }
        }
        else if (tp.substr(0,9) == "*HO_SHELL")
        {
            getline(inFile, tp); si++;
            while (tp.empty() == 0)
            {
                ShellElements.push_back(dimensions.Joint(tp));
                getline(inFile, tp); si++;
            }
        }
        else if (tp.substr(0,9) == "*HO_JOINT")
        {
            getline(inFile, tp); si++;
            getline(inFile, tp); si++;
            while (tp.length() > 1)
            {
                jointcheck = dimensions.Joint(tp);
                if (jointcheck({ 2,3,4,5 }) == RowVector4d::Zero())
                    shelljoint.push_back(jointcheck({ 0,1 }));
                else
                    Joints.push_back(jointcheck);
                getline(inFile, tp); si++;
                if (inFile.eof() != 0)
                {
                    jointcheck = dimensions.Joint(tp);
                    if (jointcheck({ 2,3,4,5 }) == RowVector4d::Zero())
                        shelljoint.push_back(jointcheck({ 0,1 }));
                    else
                        Joints.push_back(jointcheck);
                    tp = "";
                }
            }
        }
        else if (tp.substr(0,7) == "*HO_SPC")
        {
            getline(inFile, tp); si++;
            while (tp.length() > 1)
            {
                if (tp.find("U") != string::npos)
                    bcs.push_back(0);                       // 0 for Rigid Body Modes
                else if (tp.find("W") != string::npos)
                    bcs.push_back(1);                       // 1 for Warping Modes
                else if (tp.find("X") != string::npos)
                    bcs.push_back(2);                       // 2 for Distortion Modes

                BoundCs.push_back(dimensions.BCs(tp));
                getline(inFile, tp); si++;
            }

        }
        else if (tp.substr(0,9) == "*HO_FORCE")
        {
            getline(inFile, tp); si++;
            while (tp.length() > 1)
            {
                if (tp.find("F") != string::npos)
                {
                    fcs.push_back(0);                       // 0 for point force
                    Forces.push_back(dimensions.Forces(tp, 3));
                }
                else if (tp.find("B") != string::npos)
                {
                    fcs.push_back(1);                       // 1 for Warping force
                    Forces.push_back(dimensions.Forces(tp, 3));
                }
                else if (tp.find("Q") != string::npos)
                {
                    fcs.push_back(2);                       // 2 for Distortion force
                    Forces.push_back(dimensions.Forces(tp, 3));
                }
                else
                {
                    Forces.push_back(dimensions.Forces(tp, 5));
                }
                getline(inFile, tp); si++;
            }
        }
        else if (tp.substr(0,8) == "*HO_EIGR")
        {
            getline(inFile, tp); si++;
            freq = dimensions.CSMandM(tp)(2);
        }
        else if (tp.substr(0,8) == "*HO_FREQ")
        {
            getline(inFile, tp); si++;
            FRF = dimensions.BCs(tp);
        }
        si++;
    }
    inFile.close();

    ModenMat.setZero(CSMandM.size(), 3); // 색션에 대한 특징
    for (int i = 0; i < CSMandM.size(); i++)
    {
        ModenMat(i, { 0,2 }) = CSMandM[i]({ 0,2 });
        ModenMat(i, 1) = CSMandM[i](1) - 1;
    }

    Shell.setZero(ShellElements.size(), 7);
    for (int i = 0; i < ShellElements.size(); i++)
    {   
        Shell(i, 0) = i;
        Shell.block(i, 1, 1, 6) = ShellElements[i];
    }

    HO_Nodes.setZero(Nodes.size(), 4);
    for (int i = 0; i < Nodes.size(); i++)
        HO_Nodes(i, all) =  Nodes[i];

    HO_Elements.setZero(Elements.size(), 7);
    for (int i = 0; i < Elements.size(); i++)
        HO_Elements(i, all) = Elements[i]; // 한꺼번에 모아주기

    HO_Joints.setZero(Joints.size(), 6);
    for (int i = 0; i < Joints.size(); i++)
        HO_Joints(i, all) = Joints[i];

    ShellJoints.setZero(shelljoint.size(), 2);
    for (int i = 0; i < shelljoint.size(); i++)
        ShellJoints(i, all) = shelljoint[i];

    NW1.setZero(nw1.size());
    for (int i = 0; i < nw1.size(); i++)
        NW1(i) = nw1[i]; // 이게 무슨 타입인지 저장

    for (int i = 0; i < HO_Elements.rows(); i++)
    {
        for (int j = 4; j < 7; j++)
            if (abs(HO_Elements(i, j)) < 1e-9) // 너무 작은 숫자 0으로 만들어줌
                HO_Elements(i, j) = 0;
        HO_Elements(i, { 4,5,6 }) = HO_Elements(i, { 4,5,6 }) / HO_Elements(i, { 4,5,6 }).norm(); 
    }
    // Calculating Cross Section Coordinates

    int icp1 = 0, icp2 = 0, jp;
    double50 unitcheck;
    RowVector<double50, 2> csp1, csp2, csp3;
    Matrix<double50, -1, -1> UCross3;
    vector <RowVector<double50, 2>> csp;
    VectorXi cvec;
    vector<double50> Thickness;
    Matrix<double50, -1, -1> csc1;
    MatrixXi cscc1;
    vector<Matrix<double50, -1  , -1>> csc;
    vector<MatrixXi> cscc;
    string units;

    if (dim.size() > 0)
        unitcheck = abs(dim[0](2));
    else if (UserCross.size() > 0)
        unitcheck = abs(UserCross[0](0,0));
    else
        unitcheck = HO_Nodes(all, {1,2,3}).maxCoeff() / 10;

    if(unitcheck > 1)
    {
        units = "millimeters"; // 왜 미리미터 일까?
        for (int se = 0; se < sectype.size(); se++)
        {
            if (sectype[se] == 0 || sectype[se] == 1)
            {
                csc1.setZero(4, 2);
                csc1 << -dim[icp1](2) / 2 / 1000, -dim[icp1](3) / 2 / 1000,
                    dim[icp1](2) / 2 / 1000, -dim[icp1](3) / 2 / 1000,
                    dim[icp1](2) / 2 / 1000, dim[icp1](3) / 2 / 1000,
                    -dim[icp1](2) / 2 / 1000, dim[icp1](3) / 2 / 1000;

                csc.push_back(csc1); // 좌표 저장
                if (dim[icp1](1) == 0) // open
                {
                    cscc1.setZero(3, 2);
                    cscc1 << 0, 1,
                        1, 2,
                        2, 3;
                }
                else if (dim[icp1](1) == 1) // rect
                {
                    cscc1.setZero(4, 2);
                    cscc1 << 0, 1,
                        1, 2,
                        2, 3,
                        3, 0;
                }
                cscc.push_back(cscc1); // 순서 저장
                Thickness.push_back(dim[icp1](4) / 1000);// 두께 저장
                icp1++;
            }
            else if (sectype[se] == 2) //user
            {
                UCross3 = UserCross[icp2];
                csp.clear();
                for (int j = 0; j < UCross3.rows(); j++)
                {
                    csp1 = UCross3(j, { 0,1 });
                    for (int k = 0; k < 2; k++)
                    {
                        if (csp.size() == 0)
                            csp.push_back(csp1);
                        else
                        {
                            jp = 0;
                            for (int m = 0; m < csp.size(); m++)
                                if (csp1 == csp[m])
                                    jp++;
                            if (jp == 0)
                                csp.push_back(csp1);
                        }
                        csp1 = UCross3(j, { 2,3 });
                    }
                }
                csc1.setZero(csp.size(), 2);
                cscc1.setZero(UCross3.rows(), 2);
                cvec.setLinSpaced(csp.size(), 0, csp.size() - 1);
                for (int i = 0; i < cvec.size(); i++)
                {
                    for (int j = 0; j < UCross3.rows(); j++)
                    {
                        if (csp[cvec(i)] == UCross3(j, { 0,1 }))
                            cscc1(j, 0) = cvec(i);
                        if (csp[cvec(i)] == UCross3(j, { 2,3 }))
                            cscc1(j, 1) = cvec(i);
                    }
                    csc1(i, all) = csp[i] / 1000;
                }
                csc1(all, 0) = -csc1(all, 0);
                csc.push_back(csc1);
                cscc.push_back(cscc1);
                Thickness.push_back(usthick[icp2] / 1000);
                icp2++;
            }
            else if (sectype[se] == 3) // chamfer
            {
                csc1.setZero(8, 2);
                csc1 << (-dim[icp1](2) / 2 + dim[icp1](4)) / 1000, -dim[icp1](3) / 2 / 1000,
                    (dim[icp1](2) / 2 - dim[icp1](4)) / 1000, -dim[icp1](3) / 2 / 1000,
                    dim[icp1](2) / 2 / 1000, (-dim[icp1](3) / 2 + dim[icp1](4)) / 1000,
                    dim[icp1](2) / 2 / 1000, (dim[icp1](3) / 2 - dim[icp1](4)) / 1000,
                    (dim[icp1](2) / 2 - dim[icp1](4)) / 1000, dim[icp1](3) / 2 / 1000,
                    (-dim[icp1](2) / 2 + dim[icp1](4)) / 1000, dim[icp1](3) / 2 / 1000,
                    -dim[icp1](2) / 2 / 1000, (dim[icp1](3) / 2 - dim[icp1](4)) / 1000,
                    -dim[icp1](2) / 2 / 1000, (-dim[icp1](3) / 2 + dim[icp1](4)) / 1000;

                csc.push_back(csc1);
                cscc1.setZero(8, 2);
                cscc1 << 0, 1,
                    1, 2,
                    2, 3,
                    3, 4,
                    4, 5,
                    5, 6,
                    6, 7,
                    7, 0;
                cscc.push_back(cscc1);
                Thickness.push_back(dim[icp1](5) / 1000);
                icp1++;
            }
        }
    }
    else
    {
        units = "meters";
        for (int se = 0; se < sectype.size(); se++)
        {
            if (sectype[se] == 0 || sectype[se] == 1)
            {
                csc1.setZero(4, 2);
                csc1 << -dim[icp1](2) / 2, -dim[icp1](3) / 2,
                    dim[icp1](2) / 2, -dim[icp1](3) / 2,
                    dim[icp1](2) / 2, dim[icp1](3) / 2,
                    -dim[icp1](2) / 2, dim[icp1](3) / 2;

                csc.push_back(csc1);
                if (dim[icp1](1) == 0)
                {
                    cscc1.setZero(3, 2);
                    cscc1 << 0, 1,
                        1, 2,
                        2, 3;
                }
                else if (dim[icp1](1) == 1)
                {
                    cscc1.setZero(4, 2);
                    cscc1 << 0, 1,
                        1, 2,
                        2, 3,
                        3, 0;
                }
                cscc.push_back(cscc1);
                Thickness.push_back(dim[icp1](4));
                icp1++;
            }
            else if (sectype[se] == 2)
            {
                UCross3 = UserCross[icp2];
                csp.clear();
                for (int j = 0; j < UCross3.rows(); j++)
                {
                    csp1 = UCross3(j, { 0,1 });
                    for (int k = 0; k < 2; k++)
                    {
                        if (csp.size() == 0)
                            csp.push_back(csp1);
                        else
                        {
                            jp = 0;
                            for (int m = 0; m < csp.size(); m++)
                                if (csp1 == csp[m])
                                    jp++;
                            if (jp == 0)
                                csp.push_back(csp1);
                        }
                        csp1 = UCross3(j, { 2,3 });
                    }
                }
                csc1.setZero(csp.size(), 2);
                cscc1.setZero(UCross3.rows(), 2);
                cvec.setLinSpaced(csp.size(), 0, csp.size() - 1);
                for (int i = 0; i < cvec.size(); i++)
                {
                    for (int j = 0; j < UCross3.rows(); j++)
                    {
                        if (csp[cvec(i)] == UCross3(j, { 0,1 }))
                            cscc1(j, 0) = cvec(i);
                        if (csp[cvec(i)] == UCross3(j, { 2,3 }))
                            cscc1(j, 1) = cvec(i);
                    }
                    csc1(i, all) = csp[i];
                }
                csc1(all, 0) = -csc1(all, 0);
                csc.push_back(csc1);
                cscc.push_back(cscc1);
                Thickness.push_back(usthick[icp2]);
                icp2++;
            }
            else if (sectype[se] == 3)
            {
                csc1.setZero(8, 2);
                csc1 << (-dim[icp1](2) / 2 + dim[icp1](4)), -dim[icp1](3) / 2,
                    (dim[icp1](2) / 2 - dim[icp1](4)), -dim[icp1](3) / 2,
                    dim[icp1](2) / 2, (-dim[icp1](3) / 2 + dim[icp1](4)),
                    dim[icp1](2) / 2, (dim[icp1](3) / 2 - dim[icp1](4)),
                    (dim[icp1](2) / 2 - dim[icp1](4)), dim[icp1](3) / 2,
                    (-dim[icp1](2) / 2 + dim[icp1](4)), dim[icp1](3) / 2,
                    -dim[icp1](2) / 2, (dim[icp1](3) / 2 - dim[icp1](4)),
                    -dim[icp1](2) / 2, (-dim[icp1](3) / 2 + dim[icp1](4));

                csc.push_back(csc1);
                cscc1.setZero(8, 2);
                cscc1 << 0, 1,
                    1, 2,
                    2, 3,
                    3, 4,
                    4, 5,
                    5, 6,
                    6, 7,
                    7, 0;
                cscc.push_back(cscc1);
                Thickness.push_back(dim[icp1](5));
                icp1++;
            }
        }
    }
    Center coord;
    Edge edge;
    Array<double50, -1, 1> alpha, Le, cR;
    double50 beta;
    for (int i = 0; i < csc.size(); i++)
    {
        if (sectype[i] == 2) // user
        {
            alpha = edge.Anglemp(csc[i], cscc[i]);
            Le = edge.Lengthmp(csc[i], cscc[i]);
            cR = coord.Rot.Coordinates(csc[i], cscc[i], alpha, Le);
            if (cR(0) != 0 && cR(1) != 0)
            {
                csc[i](all, 0) -= cR(0) * Vector<double50, -1>::Ones(csc[i].rows());
                csc[i](all, 1) -= cR(1) * Vector<double50, -1>::Ones(csc[i].rows());
            }
        }
    }

    HO_Materials.setZero(Materials.size(), 4);
    for (int i = 0; i < Materials.size(); i++)// 물질 값 변화
    {
        if (units == "millimeters")
        {
            HO_Materials(i, { 0,2 }) = Materials[i]({ 0,2 });
            HO_Materials(i, 1) = Materials[i](1) * 1e6;
            HO_Materials(i, 3) = Materials[i](3) * 1e12;
        }
        else if (units == "meters")
        {
            HO_Materials(i, { 0,2 }) = Materials[i]({ 0,2 });
            HO_Materials(i, 1) = Materials[i](1);
            HO_Materials(i, 3) = Materials[i](3);
        }
    }

    vector<MatrixXd> shellcoord;
    vector<RowVectorXd> scoordi;
    MatrixXd scoord;
    MatrixXd sncoord;

    for (int i = 0; i < Shell.rows(); i++)
    {
        sncoord.setZero(3, 4);
        if (units == "millimeters")
        {
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < HO_Nodes.rows(); k++)
                    if (Shell(i, j + 3) == HO_Nodes(k, 0))
                    {
                        sncoord(all, j) = HO_Nodes(k, { 1,2,3 }).transpose() / 1000;
                        break;
                    }
            }
        }
        else if (units == "meters")
        {
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < HO_Nodes.rows(); k++)
                    if (Shell(i, j + 3) == HO_Nodes(k, 0))
                    {
                        sncoord(all, j) = HO_Nodes(k, { 1,2,3 }).transpose();
                        break;
                    }
            }
        }
        shellcoord.push_back(sncoord);
    }

    if (Shell.rows() > 0)
    {
        vector<int> vec1, vec2, ndtable;
        for (int i = 3; i < 7; i++)
            vec1.push_back(Shell(0, i));

        for (int i = 1; i < Shell.rows(); i++)
        {
            for (int j = 3; j < 7; j++)
                vec2.push_back(Shell(i, j));

            ndtable.resize(vec1.size() + vec2.size());

            std::sort(vec1.begin(), vec1.end()); // 정렬 이 라이브러리는 ?
            std::sort(vec2.begin(), vec2.end());

            auto it = std::set_union(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), ndtable.begin()); // ? 합집합

            ndtable.resize(it - ndtable.begin());

            vec1 = ndtable;
            vec2.clear();
        }

        for (int i = 0; i < ndtable.size(); i++)
        {
            for (int j = 0; j < HO_Nodes.rows(); j++)
                if (ndtable[i] == HO_Nodes(j, 0))
                {
                    scoordi.push_back(HO_Nodes(j, all));
                    break;
                }
        }

        scoord.setZero(scoordi.size(), 4);
        if (units == "millimeters")
        {
            for (int i = 0; i < scoordi.size(); i++)
            {
                scoord(i, 0) = scoordi[i](0);
                scoord(i, { 1,2,3 }) = scoordi[i]({ 1,2,3 }) / 1000;
            }
        }
        if (units == "meters")
        {
            for (int i = 0; i < scoordi.size(); i++)
            {
                scoord(i, 0) = scoordi[i](0);
                scoord(i, { 1,2,3 }) = scoordi[i]({ 1,2,3 });
            }
        }
    }

    VectorXd thick(Thickness.size());
    thick.setZero();
    for (int i = 0; i < Thickness.size(); i++)
        thick(i) = double(Thickness[i]);

    vector<int> iN; // 누가 몇번 들어가는지 확인하여 시작과 끝을 알 수 있다. 
    MatrixXi Node_Place(HO_Nodes.rows(), 2); // 여기에 몇번 들어가는지 저장
    for (int i = 0; i < HO_Nodes.rows(); i++)
    {
        iN.clear();
        for (int j = 0; j < HO_Elements.rows(); j++)
        {
            if (HO_Nodes(i, 0) == HO_Elements(j, 2))
                iN.push_back(j);
            if (HO_Nodes(i, 0) == HO_Elements(j, 3))
                iN.push_back(j);
        }
        Node_Place(i, 0) = HO_Nodes(i, 0);
        Node_Place(i, 1) = iN.size();
    }
    vector<int> nE;
    for (int i = 0; i < Node_Place.rows(); i++)
        if (Node_Place(i, 1) == 1)
        {
            for (int j = 0; j < HO_Elements.rows(); j++)
                if (Node_Place(i, 0) == HO_Elements(j, 2))
                    nE.push_back(i);

        }

    MatrixXi NodeEnds = Node_Place(nE, all);

    int kN, jN, iNj = 0;
    bool check1;
    vector<RowVectorXd> belem;
    MatrixXd beam1, pnt1(2, 3);
    vector<MatrixXd> beam, beamcoord;
    vector<int> csnum1;
    MatrixXd besz(HO_Joints.rows(), 8);
    besz.setZero();

    while (iNj < NodeEnds.rows())
    {
        kN = NodeEnds(iNj, 0);
        check1 = true;
        belem.clear();
        pnt1.setZero();
        while (check1)
        {
            jN = kN;
            for (int i = 0; i < HO_Elements.rows(); i++)
            {
                for (int j = 0; j < HO_Joints.rows(); j++)
                {
                    if (kN == HO_Joints(j, 0))
                    {
                        besz(j, 0) = beam.size();
                        besz(j, 1) = HO_Joints(j, 1) - 1;
                    }
                    if (kN == HO_Joints(j, 3))
                    {
                        besz(j, 4) = beam.size();
                        besz(j, 5) = HO_Joints(j, 4) - 1;
                    }
                }
                if (kN == HO_Elements(i, 2))
                {
                    belem.push_back(HO_Elements(i, all));
                    kN = HO_Elements(i, 3);
                    break;
                }
            }
            if (kN == jN)
            {
                check1 = false;
                beam1.setZero(belem.size(), 7);
                for (int k = 0; k < belem.size(); k++)
                    beam1(k, all) = belem[k];
                if (beam1.rows() != 0)
                {
                    beam.push_back(beam1);
                    csnum1.push_back(beam1(0, 1));
                }
                for (int i = 0; i < HO_Nodes.rows(); i++)
                {
                    if (beam1(0, 2) == HO_Nodes(i, 0))
                        pnt1(0, all) = HO_Nodes(i, { 1, 2, 3 });
                    if (beam1(beam1.rows() - 1, 3) == HO_Nodes(i, 0))
                        pnt1(1, all) = HO_Nodes(i, { 1, 2, 3 });
                }
                beamcoord.push_back(pnt1);
            }
        }
        iNj++;
    }

    RowVectorXi csnum(csnum1.size());
    for (int i = 0; i < csnum1.size(); i++)
    {
        for (int j = 0; j < CSMandM.size(); j++)
            if (csnum1[i] == CSMandM[j](0))
                csnum(i) = j;
    }

    Cross_Section csLA;
    int b1, b2, e1, e2;
    RowVectorXd p11, p12, p21, p22;

    for (int i = 0; i < besz.rows(); i++)
    {
        b1 = besz(i, 0), b2 = besz(i, 4);
        e1 = besz(i, 1), e2 = besz(i, 5);
        besz(i, 2) = double(csLA.EdgeLength(csc[csnum(b1)], cscc[csnum(b1)])(e1)) * HO_Joints(i, 2);
        besz(i, 6) = double(csLA.EdgeLength(csc[csnum(b2)], cscc[csnum(b2)])(e2)) * HO_Joints(i, 5);

        p11 = beamcoord[b1](0, all);
        for (int j = 0; j < HO_Nodes.rows(); j++)
            if (HO_Joints(i, 0) == HO_Nodes(j, 0))
                p12 = HO_Nodes(j, { 1, 2, 3 });

        if(units == "millimeters")
            besz(i, 3) = ((p12 - p11) / 1000).norm();
        else
            besz(i, 3) = (p12 - p11).norm();


        p21 = beamcoord[b2](0, all);
        for (int j = 0; j < HO_Nodes.rows(); j++)
            if (HO_Joints(i, 3) == HO_Nodes(j, 0))
                p22 = HO_Nodes(j, { 1, 2, 3 });

        if(units == "millimeters")
            besz(i, 7) = ((p22 - p21) / 1000).norm();
        else
            besz(i, 7) = (p22 - p21).norm();
    }

    HO_PShell.setZero(PShell.size(), 3);
    for (int i = 0; i < PShell.size(); i++)
    {
        HO_PShell(i, {0, 1}) = PShell[i]({0, 1});
        if (units == "millimeters")
            HO_PShell(i, 2) = PShell[i](2) / 1000;
        if (units == "meters")
            HO_PShell(i, 2) = PShell[i](2);
    }

    MatrixXd xdir(beam.size(), 3), pnts(beam.size() * 2, 3), bbns;
    MatrixXi pnt_con_s(beam.size(), 2);
    xdir.setZero(), pnts.setZero();
    RowVector3d ydir, zdir, zi;
    vector<RowVector3d> Bbns;
    VectorXd zcoordi;
    vector<VectorXd> zcoord;
    double nd1, nd2;

    for (int i = 0; i < beam.size(); i++)
    {
        ydir = beam[i](0, { 4, 5, 6 });
        nd1 = beam[i](0, 2);
        nd2 = beam[i](beam[i].rows() - 1, 3);
        for (int j = 0; j < HO_Nodes.rows(); j++)
        {
            if (nd1 == HO_Nodes(j, 0))
                p11 = HO_Nodes(j, { 1, 2, 3 });
            if (nd2 == HO_Nodes(j, 0))
                p12 = HO_Nodes(j, { 1, 2, 3 });
        }

        if(units == "millimeters")
        {
            zi = (p12 - p11) / 1000;
            zdir = zi / zi.norm();
            xdir(i, all) = ydir.cross(zdir);
            pnts({ 2 * i, 2 * i + 1 }, all) = beamcoord[i] / 1000;
        }
        else
        {
            zi = (p12 - p11);
            zdir = zi / zi.norm();
            xdir(i, all) = ydir.cross(zdir);
            pnts({ 2 * i, 2 * i + 1 }, all) = beamcoord[i];
        }

        for (int k = 0; k < ShellJoints.rows(); k++)
            if (nd1 == ShellJoints(k, 0))
                Bbns.push_back({ double(i), 0, ShellJoints(k, 1) });
            else if (nd1 == ShellJoints(k, 1))
                Bbns.push_back({ double(i), 0, ShellJoints(k, 0) });

        zcoordi.setZero(beam[i].rows() + 1);
        for (int j = 0; j < beam[i].rows(); j++)
        {
            nd2 = beam[i](j, 3);
            for (int k = 0; k < HO_Nodes.rows(); k++)
                if (nd2 == HO_Nodes(k, 0))
                    p12 = HO_Nodes(k, { 1, 2, 3 });

            for (int k = 0; k < ShellJoints.rows(); k++)
                if (nd2 == ShellJoints(k, 0))
                    Bbns.push_back({ double(i), double(j + 1), ShellJoints(k, 1) });
                else if (nd2 == ShellJoints(k, 1))
                    Bbns.push_back({ double(i), double(j + 1), ShellJoints(k, 0) });

            if(units == "millimeters")
            {
                zi = (p12 - p11) / 1000;
                zcoordi(j + 1) = sqrt(pow(zi(0), 2) + pow(zi(1), 2) + pow(zi(2), 2));
            }
            else
            {
                zi = (p12 - p11);
                zcoordi(j + 1) = sqrt(pow(zi(0), 2) + pow(zi(1), 2) + pow(zi(2), 2));
            }
        }
        zcoord.push_back(zcoordi);
    }

    bbns.setZero(Bbns.size(), 3);
    for (int i = 0; i < Bbns.size(); i++)
        bbns(i, all) = Bbns[i];

    vector<RowVector2d> bni, bci;
    vector<int> sni, bi, shi;
    int bcci = 0;

    for (int i = 0; i < Forces.size(); i++)
    {
        if (beam.size() > 0)
            for (int j = 0; j < beam.size(); j++)
            {
                for (int k = 0; k < beam[j].rows(); k++)
                {
                    if (Forces[i](0) == beam[j](k, 2))
                    {
                        bni.push_back({ j,k });
                        break;
                    }
                    else if (Forces[i](0) == beam[j](k, 3))
                    {
                        bni.push_back({ j,k + 1 });
                        break;
                    }
                }
            }
        else
            bni.push_back({ -1,Forces[i](0) });

    }

    for (int i = 0; i < BoundCs.size(); i++)
    {
        bcci = bci.size();
        for (int j = 0; j < beam.size(); j++)
        {
            for (int k = 0; k < beam[j].rows(); k++)
            {
                if (BoundCs[i](0) == beam[j](k, 2))
                {
                    bci.push_back({ j,k });
                    break;
                }
                else if (BoundCs[i](0) == beam[j](k, 3))
                {
                    bci.push_back({ j,k + 1 });
                    break;
                }
            }
        }

        if (bcci == bci.size())
        {
            sni.push_back(BoundCs[i](0));
            shi.push_back(i);
        }
    }


    FNode.setZero(fcs.size(), 5);
    FEdge.setZero(Forces.size() - fcs.size(), 6);
    int fni = 0, fei = 0;
    for (int i = 0; i < Forces.size(); i++)
    {
        cout << Forces[i] << endl << endl;
        if (Forces[i].size() == 3)
        {
            FNode(fni, 0) = fcs[fni];
            FNode(fni, { 1,2 }) = bni[i];
            FNode(fni, { 3,4 }) = Forces[i]({ 1,2 });
            fni++;
        }
        else if (Forces[i].size() == 5)
        {
            FEdge(fei, { 0,1 }) = bni[i];
            FEdge(fei, { 2,3,4,5 }) = Forces[i]({ 1,2,3,4 });
            fei++;
        }
    }

    HO_BCs.setZero(bci.size(), 5);
    for (int i = 0; i < bci.size(); i++)
    {
        HO_BCs(i, 0) = bcs[i];
        HO_BCs(i, { 1,2 }) = bci[i];
        HO_BCs(i, { 3,4 }) = BoundCs[i]({ 1,2 });
    }

    MatrixXd HO_SBCs;
    HO_SBCs.setZero(shi.size(), 4);
    for (int i = 0; i < shi.size(); i++)
    {
        HO_SBCs(i, 0) = bcs[shi[i]];
        HO_SBCs(i, 1) = sni[i];
        HO_SBCs(i, { 2,3 }) = BoundCs[shi[i]]({ 1,2 });
    }

    VectorXi Nel, Nn; Nel.setZero(zcoord.size());
    for (int i = 0; i < zcoord.size(); i++)
        Nel(i) = zcoord[i].size() - 1;
    Nn = Nel + VectorXi::Ones(Nel.size());

    vector<RowVector3d> pnts1;
    for (int i = 0; i < pnts.rows(); i++)
        pnts1.push_back(pnts(i, all));

    std::sort(pnts1.begin(), pnts1.end(), [](Eigen::RowVectorXd const& t1, Eigen::RowVectorXd const& t2) { return t1(0) < t2(0); });
    auto it = std::unique(pnts1.begin(), pnts1.end());
    pnts1.resize(std::distance(pnts1.begin(), it));

    MatrixXd pnt_s(pnts1.size(), 3);
    for (int i = 0; i < pnts1.size(); i++)
        pnt_s(i, all) = pnts1[i];

    if(units == "millimeters")
    {
        for (int i = 0; i < beamcoord.size(); i++)
            for (int j = 0; j < pnt_s.rows(); j++)
            {
                if (beamcoord[i](0, all) / 1000 == pnt_s(j, all))
                    pnt_con_s(i, 0) = j;
                if (beamcoord[i](1, all) / 1000 == pnt_s(j, all))
                    pnt_con_s(i, 1) = j;
            }
    }
    else
    {
        for (int i = 0; i < beamcoord.size(); i++)
            for (int j = 0; j < pnt_s.rows(); j++)
            {
                if (beamcoord[i](0, all) == pnt_s(j, all))
                    pnt_con_s(i, 0) = j;
                if (beamcoord[i](1, all) == pnt_s(j, all))
                    pnt_con_s(i, 1) = j;
            }
    }

    return make_tuple(csc, cscc, csnum, xdir, pnt_s, pnt_con_s, besz, thick, ModenMat, HO_Materials, FNode, FEdge, HO_BCs, zcoord, Nel, Nn, analysis, freq, NW1, FRF, units, HO_SBCs, HO_PShell, Shell, bbns, shellcoord, scoord);
}
