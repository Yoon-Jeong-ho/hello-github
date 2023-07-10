#ifndef HIGHERORDERBEAMLIB_H_INCLUDED
#define HIGHERORDERBEAMLIB_H_INCLUDED

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>

using namespace Eigen;
using namespace std;
using namespace boost::multiprecision;
using namespace boost::math;
typedef cpp_dec_float_50 double50;

/*
	Reads Input File
*/
tuple< vector<Matrix<cpp_dec_float_50, -1, -1>>, vector<MatrixXi>, RowVectorXi, MatrixXd, MatrixXd, MatrixXi, MatrixXd, VectorXd, MatrixXi, MatrixXd, MatrixXd, MatrixXd, MatrixXd, vector<VectorXd>, VectorXi, VectorXi, string, int, ArrayXi, RowVector3d, string> inputReader(string txtFile);

/*
	Derives Cross-Section Modes
	Requirements:
		csc: Coordinates (Matrix double)
		cscc: Coordinates Connectivity
*/
tuple<vector<MatrixXd>, vector<MatrixXd>, vector<MatrixXd>, MatrixXi> Modegen(Matrix<cpp_dec_float_50, -1, -1>, MatrixXi, int, MatrixXi);

/*
	Derives Rigid Body Modes of the Cross Section
	Requirements:
		NE: Number of Edges (int)
		csc: Coordinates (Matrix double)
		alpha: Edge Angles (Array double)
		betaT: Torsion Angle (double)
		cT: Torsion Center (Array double)
		cR: Rotation Center (Array double)
		betaR: Rotation Angle (double)
*/
tuple<vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>> Rigid(Index, Matrix<cpp_dec_float_50, -1, -1>, MatrixXi, Array<cpp_dec_float_50, -1, 1>, cpp_dec_float_50, Array<cpp_dec_float_50, -1, 1>, Array<cpp_dec_float_50, -1, 1>, cpp_dec_float_50);

tuple<vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, ArrayXi> LW_InXD(Matrix<cpp_dec_float_50, -1, -1>, MatrixXi, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, Array<cpp_dec_float_50, -1, 1>);

/*
	Returns S-Directional Shape Functions of the Extensional DIstortion modes of the Cross Section
	Requirements:
		csc: Coordinates (Matrix double)
		cscc: A matrix of type int containing cooridnate connections
		coef_x: 3 Coefficient Matrices
		Le: An array of type double containing Edge Lengths
*/
tuple<vector<Matrix<double50, -1, -1>>, ArrayXi> ExtDist(Matrix<cpp_dec_float_50, -1, -1>, MatrixXi, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, Array<cpp_dec_float_50, -1, 1>);

/*
	Returns Wall Bending Modes for the Cross Section
*/
tuple<vector<Matrix<double50, -1, -1>>, ArrayXi> WallBend(Matrix<cpp_dec_float_50, -1, -1>, MatrixXi, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, Array<cpp_dec_float_50, -1, 1>, int);

/*
	Returns Non-Linear Warping Modes for the Cross Section
*/
tuple<vector<Matrix<double50, -1, -1>>, ArrayXi> NLW(Matrix<cpp_dec_float_50, -1, -1>, MatrixXi, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, vector<Matrix<double50, -1, -1>>, Array<cpp_dec_float_50, -1, 1>);

/*
	Saves Tensors in a csv file
*/
void saveData(string fileName, Tensor<cpp_dec_float_50, 3> tensor);

/*
	Calculates Stiffness Matrix
*/
tuple<MatrixXi, MatrixXi, MatrixXd> stiffness(Matrix<double, -1, -1>, MatrixXi, vector<MatrixXd>, vector<MatrixXd>, vector<MatrixXd>, double, double, double, double, VectorXd, int);

/*
	Calculates Mass Matrix
*/
tuple<MatrixXi, MatrixXi, MatrixXd> Mass(Matrix<double, -1, -1>, MatrixXi, vector<MatrixXd>, vector<MatrixXd>, vector<MatrixXd>, double, double, double, double, VectorXd, int);

vector<Array<double, 8, 1>> jntcon(vector<MatrixXd>, vector<MatrixXi>, MatrixXd, MatrixXi, MatrixXd, VectorXi, MatrixXi);

tuple<MatrixXd, MatrixXd> jntdisp(MatrixXd, MatrixXi, vector<MatrixXd>, vector<MatrixXd>, vector<MatrixXd>, int, double, Vector3d, int);

MatrixXd bendrot(MatrixXd, MatrixXi, vector<MatrixXd>, vector<MatrixXd>, vector<MatrixXd>, Vector2i);

tuple<MatrixXi, MatrixXd, MatrixXd, double> deformed(MatrixXd, MatrixXi, vector<MatrixXd>, vector<MatrixXd>, vector<MatrixXd>, VectorXd, VectorXd, double);

MatrixXd stress3D(MatrixXd, MatrixXi, Tensor<double, 3>, Tensor<double, 3>, Tensor<double, 3>, double, VectorXd, VectorXd, double, double, double);

tuple<double, double> c2es(MatrixXd, MatrixXi, int);
#endif // !HIGHERORDERBEAMLIB_H_INCLUDED
