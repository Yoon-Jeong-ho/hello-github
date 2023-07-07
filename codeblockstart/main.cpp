#include<fstream>
#include<iostream>
#include<string>
#include<cmath>
#include <vector>
#include <sstream>
#include <Eigen/Core>
#include<Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <iomanip>
#include <thread>
#include <chrono>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace Eigen;
using namespace std;
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
            if (temp.length() > 0) {
                if (stringstream(temp) >> found) {}
                temp = "";
            }
        }
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
            if (temp.length() > 0) {
                if (stringstream(temp) >> found(si)) {
                    //cout << found << " ";
                    si++;
                }
            }
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
            if (temp.length() > 0) {
                if (stringstream(temp) >> found(si)) {
                    si++;
                }
            }
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
            if (temp.length() > 0) {
                if (stringstream(temp) >> found(si)) {
                    si++;
                }
            }
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
            if (temp.length() > 0) {
                if (stringstream(temp) >> found(si)) {
                    si++;
                }
            }
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
            if (temp.length() > 0) {
                if (stringstream(temp) >> found(si)) {
                    si++;
                }
            }
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
            if (temp.length() > 0) {
                if (stringstream(temp) >> found(si)) {
                    si++;
                }
            }
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
            if (temp.length() > 0) {
                if (stringstream(temp) >> found(si)) {
                    si++;
                }
            }
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
            if (temp.length() > 0) {
                if (stringstream(temp) >> found(si)) {
                    si++;
                }
            }
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
            if (temp.length() > 0) {
                if (stringstream(temp) >> found(si)) {
                    si++;
                }
            }
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

int main() {
	ifstream in("Tjoint_diap_xdir2.hob");
	string str, analysis, model, material, sectiontype,tp;
	int intmp;
	char tmp, jointtype, spctype[2], f;
	vector<MatrixXd> eigvectors;
	VectorXd material_properties(4), sectionlength(2), sectionthick(2), force(3);
	MatrixXd eigvec(32, 7), nodecoordination(46, 4), storebeam(44, 7), joint(6, 6), spc(2, 3), coordinate(4, 2), seq(4, 2);
	eigvec.setZero();
	vector<RowVector3i> CSMandM;
    vector<RowVector4d> Nodes, Materials;
    MatrixXd HO_Nodes, HO_Elements, HO_Joints, HO_Materials, HO_BCs, FEdge, FNode, HO_PShell, Shell, ShellJoints;
    vector<RowVectorXd> Elements, Joints, UCross1, BoundCs, Forces, PShell, ShellElements;
    RowVectorXd jointcheck;
    RowVector2d csdim;
    RowVector3d csdimch, FRF;
	Extract dimensions;
    double csthick;
    vector<double> usthick;
    int si = 0, csnumi = 0, freq = 0;
    vector<int> sec, sectype, bcs, fcs, nw1;
    vector<RowVector<double50, -1>> dim;
    vector<Matrix<double50, -1, -1>> UserCross;
    RowVector<double50, -1> csType(5), csTypei(6);
    Matrix<double50, -1, -1> UCross2;
    MatrixXi ModenMat;
    ArrayXi NW1;
    vector<RowVector2d> shelljoint;
	material_properties.setZero();
    while (getline(in, str))
    {
        if (str.substr(0, 11) == "*HO_SECTION")
        {
            //sec.push_back(si); // Number of Cross Sections and their places in the test file
            getline(in, tp); si++;
            getline(in, tp); si++;
            CSMandM.push_back(dimensions.CSMandM(tp));
            if (tp.find("OPEN") != string::npos)
            {
                sectype.push_back(0); // 0 means open Cross Section
                getline(in, tp); si++;
                csdim(0) = dimensions.EdgeLandT(tp)(0); // Extracting Edge Length and Thickness
                getline(in, tp); si++;
                csdim(1) = dimensions.EdgeLandT(tp)(0);
                csthick = dimensions.EdgeLandT(tp)(1);
                csType << csnumi + 1, sectype[csnumi], csdim, csthick; // 결정된 모양에 따른 숫자 넘겨주기
                dim.push_back(csType);
                nw1.push_back(2);
            }
            else if (tp.find("RECT") != string::npos)
            {
                sectype.push_back(1); // 1 means closed cross section
                getline(in, tp); si++;
                csdim(0) = dimensions.EdgeLandT(tp)(0); // Extracting Edge Length and Thickness
                getline(in, tp); si++;
                csdim(1) = dimensions.EdgeLandT(tp)(0);
                csthick = dimensions.EdgeLandT(tp)(1);
                csType << csnumi + 1, sectype[csnumi], csdim, csthick;
                dim.push_back(csType);
                nw1.push_back(1);
                cout <<dim[0]<<endl<<endl;
                cout << csType;
            }
            else if (tp.find("USER") != string::npos)
            {
                sectype.push_back(2); // 2 means User Defined cross section
                UCross1.clear();
                getline(in, tp); si++;
                while (tp.length() > 1)
                {
                    UCross1.push_back(dimensions.UserCross(tp));
                    getline(in, tp); si++;
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
                getline(in, tp); si++;
                csdimch(0) = dimensions.EdgeLandT(tp)(0); // Extracting Edge Length and Thickness
                getline(in, tp); si++;
                csdimch(1) = dimensions.EdgeLandT(tp)(0);
                getline(in, tp); si++;
                csdimch(2) = dimensions.EdgeLandT(tp)(0);
                csthick = dimensions.EdgeLandT(tp)(1);
                csTypei << csnumi + 1, sectype[csnumi], csdimch, csthick;
                dim.push_back(csTypei);
                nw1.push_back(1);
            }
            csnumi++;
        }
        else if (str.substr(0,12) == "*HO_ANALYSIS")
        {
            getline(in, tp); si++;
            analysis = tp;
        }
        else if (str == "*HO_JOINT")
        {
            getline(in, str);
            getline(in, str);
            while (str.length() > 1)
            {
                jointcheck = dimensions.Joint(str);


                if (jointcheck({ 2,3,4,5 }) == RowVector4d::Zero()) {
                    shelljoint.push_back(jointcheck({ 0,1 }));
                }

                else {
                    Joints.push_back(jointcheck);
                }
                getline(in, str);
                if (in.eof() != 0)
                {
                   jointcheck = dimensions.Joint(str);
                    if (jointcheck({ 2,3,4,5 }) == RowVector4d::Zero())
                        shelljoint.push_back(jointcheck({ 0,1 }));
                    else
                        Joints.push_back(jointcheck);
                    str = "";
                }
            }
        }
        else if (str.substr(0, 12) == "*HO_MATERIAL")
        {
            getline(in, tp); si++;
            getline(in, tp); si++;
            // Length should be checked with > 1 not != 1
            //while (tp.length() != 1)
            while (tp.length() > 1)
            {
                Materials.push_back(dimensions.Material(tp));
                getline(in, tp); si++;
            }
        }
        else if (str.substr(0, 8) == "*HO_NODE")
        {
            getline(in, str);
            while (str.length() > 1)
            {
                Nodes.push_back(dimensions.Nodes(str));
                getline(in, str);
            }
        }
        else if (str.substr(0, 8) == "*HO_BEAM")
        {
            getline(in, str);
            while (str.length() > 1)
            {
                Elements.push_back(dimensions.Elements(str));
                getline(in, str);
            }
        }
        else if (str.substr(0, 7) == "*HO_SPC")
        {
            getline(in, tp); si++;
            while (tp.length() > 1)
            {
                if (tp.find("U") != string::npos)
                    bcs.push_back(0);                       // 0 for Rigid Body Modes
                else if (tp.find("W") != string::npos)
                    bcs.push_back(1);                       // 1 for Warping Modes
                else if (tp.find("X") != string::npos)
                    bcs.push_back(2);                       // 2 for Distortion Modes

                BoundCs.push_back(dimensions.BCs(tp));
                getline(in, tp); si++;
            }

        }
        else if (str.substr(0, 9) == "*HO_FORCE")
        {
            getline(in, tp); si++;
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
                getline(in, tp); si++;
            }
        }
    }
    HO_Nodes.setZero(Nodes.size(), 4);

    ModenMat.setZero(CSMandM.size(), 3);
    for (int i = 0; i < CSMandM.size(); i++) // 섹션에 대한 특징 저장
    {
        ModenMat(i, { 0,2 }) = CSMandM[i]({ 0,2 });
        ModenMat(i, 1) = CSMandM[i](1) - 1;

    }


    for (int i = 0; i < Nodes.size(); i++) {
        HO_Nodes(i, all) = Nodes[i];

    }


    HO_Elements.setZero(Elements.size(), 7);
    for (int i = 0; i < Elements.size(); i++)
        HO_Elements(i, all) = Elements[i];

    HO_Joints.setZero(Joints.size(), 6);
    for (int i = 0; i < Joints.size(); i++)
        HO_Joints(i, all) = Joints[i];

    NW1.setZero(nw1.size());
    for (int i = 0; i < nw1.size(); i++)
        NW1(i) = nw1[i]; // 이게 무슨 타입인지 저장


    for (int i = 0; i < HO_Elements.rows(); i++)
    {
        for (int j = 4; j < 7; j++) {
            if (abs(HO_Elements(i, j)) < 1e-9) {
                HO_Elements(i, j) = 0;
            }
        }
        HO_Elements(i, { 4,5,6 }) = HO_Elements(i, { 4,5,6 }) / HO_Elements(i, { 4,5,6 }).norm(); // 너무 작은 숫자 0으로 만들어줌
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
    vector<Matrix<double50, -1, -1>> csc;
    vector<MatrixXi> cscc;
    string units;

    if (dim.size() > 0)
        unitcheck = abs(dim[0](2));
    cout <<"unitcheck"<<unitcheck<<endl<<endl;

    if (unitcheck > 1)
    {
        units = "millimeters";
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


    HO_Materials.setZero(Materials.size(), 4);
    for (int i = 0; i < Materials.size(); i++)
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

    VectorXd thick(Thickness.size());
    thick.setZero();
    for (int i = 0; i < Thickness.size(); i++)
        thick(i) = double(Thickness[i]);
    vector<int> iN;
    MatrixXi Node_Place(HO_Nodes.rows(), 2);
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
        besz(i, 2) = double(csLA.EdgeLength(csc[csnum(b1)], cscc[csnum(b1)])(e1)) * HO_Joints(i, 2);  //csc = 모양에 따른 좌표
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
        nd1 = beam[i](0, 2); //
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

            for (int k = 0; k < ShellJoints.rows(); k++) // 이게 날아갔네
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
        else // 이것도 날라가쎈
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
    for(int i = 0 ; i <beamcoord.size();++i)

    in.close();


    for(int i = 0 ; i <csc.size() ; ++i)
    cout << "csc["<<i<<"] =\n" <<csc[i] << endl<<endl;//
    for(int i = 0 ; i <cscc.size() ; ++i)
    cout << "cscc["<<i<<"] =\n"<<cscc[i]<<endl<<endl;
    cout<<"csnum =\n"<< csnum<<endl<<endl;
    cout<<"xdir=\n"<< xdir<<"\n\n";
    cout<< "pnt_s =\n"<< pnt_s<<endl<<endl;
    cout<< "pnt_con_s =\n"<<pnt_con_s<<endl<<endl;
    cout<< "besz =\n"<<besz<<endl<<endl;
    cout<< "thick =\n"<<thick<<endl<<endl;
    cout<< "ModenMat =\n"<<ModenMat<<endl<<endl;
    cout<< "HO_Materials =\n"<<HO_Materials<<endl<<endl;
    cout<< "FNode =\n"<<FNode<<endl<<endl;
    cout<< "FEdge =\n"<<FEdge<<endl<<endl;
    cout<< "HO_BCs =\n"<<HO_BCs<<endl<<endl;
    for(int i = 0 ; i <zcoord.size() ; ++i)
    cout << "zcoord["<<i<<"] =\n" <<zcoord[i]<<endl<<endl;
    cout<< "Nel =\n"<<Nel<<endl<<endl;
    cout<< "Nn =\n"<<Nn<<endl<<endl;
    cout<< "analysis =\n"<<analysis<<endl<<endl;
    cout<< "freq =\n"<<freq<<endl<<endl;
    cout<< "NW1 =\n"<<NW1<<endl<<endl;
    cout<< "FRF =\n"<<FRF<<endl<<endl;
    cout<< "units =\n"<<units<<endl<<endl;


	return 0;
}



int main4() {
	ifstream in("Tjoint_diap_xdir2.hob");
	string str, analysis, model, material, sectiontype;
	int intmp;
	char tmp, jointtype, spctype[2], f;
	vector<MatrixXd> eigvectors;
	VectorXd material_properties(4), sectionlength(2), sectionthick(2), force(3);
	MatrixXd eigvec(32, 7), nodecoordination(46, 4), storebeam(44, 7), joint(6, 6), spc(2, 3), coordinate(4, 2), seq(4, 2);
	eigvec.setZero();
	material_properties.setZero();

	if (in.is_open()) {
		while (!in.eof()) {
			getline(in, str);
			if (str == "*HO_ANALYSIS") {
				getline(in, analysis);
			}
			else if (str == "*HO_MODEL") {
				getline(in, model);
			}
			else if (str == "*HO_MATERIAL") {
				getline(in, material);
				for (int i = 0; i < 3; i++) {
					in >> material_properties(i);
					in >> tmp;
				}
				in >> material_properties(3);
			}
			else if (str == "*HO_SECTION") {
				getline(in, str);
				in >> sectiontype;
				getline(in, str);
				in >> sectionlength(0);
				in >> tmp;
				in >> sectionthick(0);
				in >> sectionlength(1);
				in >> tmp;
				in >> sectionthick(1);
				coordinate(0, 0) = sectionlength(0) * 0.5;
				coordinate(0, 1) = sectionlength(1) * 0.5 * (-1);
				coordinate(1, 0) = sectionlength(0) * 0.5;
				coordinate(1, 1) = sectionlength(1) * 0.5;
				coordinate(2, 0) = sectionlength(0) * 0.5 * (-1);
				coordinate(2, 1) = sectionlength(1) * 0.5;
				coordinate(3, 0) = sectionlength(0) * 0.5 * (-1);
				coordinate(3, 1) = sectionlength(1) * 0.5 * (-1);
				seq << 1, 2,
					2, 3,
					3, 4,
					4, 1;

			}
			else if (str == "*HO_NODE") {
				for (int i = 0; i < 46; i++) {
					for (int j = 0; j < 3; j++) {
						in >> nodecoordination(i, j);
						in >> tmp;
					}
					in >> nodecoordination(i, 3);
				}
			}
			else if (str == "*HO_BEAM") {
				for (int i = 0; i < 44; i++) {
					for (int j = 0; j < 6; j++) {
						in >> storebeam(i, j);
						in >> tmp;
					}
					in >> storebeam(i, 6);
				}
			}
			else if (str == "*HO_JOINT") {
				in >> jointtype;
				getline(in, str);
				for (int i = 0; i < 6; i++) {
					for (int j = 0; j < 5; j++) {
						in >> joint(i, j);
						in >> tmp;
					}
					in >> joint(i, 5);
				}
			}
			else if (str == "*HO_SPC") {
				for (int i = 0; i < 2; i++) {
					in >> spc(i, 0);
					in >> tmp;
					in >> spctype[i];
					in >> tmp;
					in >> spc(i, 1);
					in >> tmp;
					in >> spc(i, 2);
				}
			}
			else if (str == "*HO_FORCE") {
				in >> force(0);
				in >> tmp;
				in >> f;
				in >> tmp;
				in >> force(1);
				in >> tmp;
				in >> force(2);
			}

		}

	}
	else {
		cout << "file open fail";
		return 1;
	}

	int start[2] = {-1,-1}, count = 0;
	for (int i = 0; i < 44; i++) {
		for (int j = 0; j < 44; j++) {
			if (storebeam(i, 2) == storebeam(j, 3))break;
			if (j == 43) start[0] = i;
		}
		if (start[0] != -1) break;
	}

	for (int i = start[0] + 1; i < 44; i++) {
		for (int j = 0; j < 44; j++) {
			if (storebeam(i, 2) == storebeam(j, 3))break;
			if (j == 43) start[1] = i;
		}
		if (start[1] != -1) break;
	}

	vector<RowVectorXd> beams;
	MatrixXd beami;
	vector<MatrixXd> Beam;
	intmp = 0;
	for (int k = 0; k < 2; k++)
	{
		beams.clear();
		beams.push_back(storebeam(start[k], all));
		for (int i = 0; i < storebeam.rows(); i++)
		{
			for (int j = 0; j < storebeam.rows(); j++)
			{
				if (storebeam(start[k],3) == storebeam(j, 2))
				{
					start[k] = j;
					beams.push_back(storebeam(j, all));
					break;
				}

			}
		}
		beami.setZero(beams.size(), 7);
		for (int i = 0; i < beams.size(); i++)
				beami(i, all) = beams[i];
		Beam.push_back(beami);
	}

	cout << "This is the " << analysis << " Analysis of " << model << " Model." << endl;
	cout << "element " << material_properties(0) << " " << material << "(E = " << material_properties(1) << ", nu = " << material_properties(2) << ", rho = " << material_properties(3) << ")" << endl;
	cout << "Section Type = " << sectiontype << " sectionlength = " << sectionlength.transpose() << " sectionthick = " << sectionthick.transpose() << endl;
	cout << " coordinate is = \n" << coordinate << endl << "seq is = \n" << seq << endl;
	cout << nodecoordination << endl;
	cout << "joint type is = " << jointtype << endl << joint << endl;
	cout << "spc = \n" << spc << endl;
	cout << "force = " << force.transpose() << endl;
	for (int i = 0; i < 2; i++)cout << Beam[i] << endl << endl;
	in.close();
	return 0;
}
