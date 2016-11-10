/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #1:
 * svmmath.cpp
 *		a routine for intersecting >2 lines (for vanishing point
 *		computation);
 *		routines for computing the homography for the reference
 *		plane and arbitrary polygons
 **************************************************************/

#pragma warning(disable : 4996)

#include "svmmath.h"
#include "jacob.h"
#include "vec.h"
#include <cstring>
#include <cstdio>
#include <assert.h>
#include <iostream>

#include "Eigen/Core"
#include "MinEig.h"

using namespace Eigen;
using namespace std;

//
// TODO 1: BestFitIntersect()
//		Given lines, the list of 3 or more lines to be intersected,
//		find the best fit intersection point.
//		See http://www-2.cs.cmu.edu/~ph/869/www/notes/vanishing.txt.
//
SVMPoint BestFitIntersect(const std::list<SVMLine> &lines, int imgWidth, int imgHeight)
{
  // check
  if (lines.size() < 2)
    {
      fprintf(stderr, "Not enough lines to compute the best fit.");
      abort();
    }

  SVMPoint bestfit;
  list<SVMLine>::const_iterator iter;

  // To accumulate stuff

  // Dynamic rows, fixed cols.
  typedef Matrix<double, Dynamic, 3, RowMajor> Matrix3;
  // int numLines = (int) lines.size();
  // Matrix3 A = Matrix3::Zero(numLines, 3);
  Matrix3 M = Matrix3::Zero(3, 3);
  Matrix3 tmp = Matrix3::Zero(3, 3);

  // Transformation for numerical stability

  // Note: iterate through the lines list as follows:
  //		for (iter = lines.begin(); iter != lines.end(); iter++) {
  //			...iter is the pointer to the current line...
  //		}
  // Note: Function to find eigenvector with smallest eigenvalue is MinEig(A, eval, evec)
  //
  /******** BEGIN TODO ********/
  double inf = 1e-10;
  for (iter = lines.begin(); iter != lines.end(); iter++){
    // cout<<"Start lines iterations"<<endl;
    //endpoint1
    double x1 = iter->pnt1->u;
    double y1 = iter->pnt1->v;
    double w1 = iter->pnt1->w;

    //endpoint2
    double x2 = iter->pnt2->u;
    double y2 = iter->pnt2->v;
    double w2 = iter->pnt2->w;

    Vec3d e1 = Vec3d(x1, y1, w1);
    Vec3d e2 = Vec3d(x2, y2, w2);
    Vec3d l = cross(e1, e2);
    
    tmp << l[0]*l[0], l[0]*l[1], l[0]*l[2], 
           l[0]*l[1], l[1]*l[1], l[1]*l[2],
           l[0]*l[2], l[1]*l[2], l[2]*l[2];

    cout << "tmp:" <<tmp<<endl;
    M += tmp;
}
 cout<<"Finish lines iterations"<<endl;
  
   cout<<"M: "<<M<<endl;
  double eval;
  double *evec;
  MinEig(M, eval, evec);
  
  // cout<<eval<<endl;
  if(evec[2]<inf){
    bestfit.u = evec[0]/evec[2];
    bestfit.v = evec[1]/evec[2];
  //bestfit.w = evec[2];
  }else{
    bestfit.u = evec[0];
    bestfit.v = evec[1];
    bestfit.w = evec[2];
  }
  

  printf("TODO: svmmath.cpp:61\n"); 
  fl_message("TODO: svmmath.cpp:61\n");

	/******** END TODO ********/

 return bestfit;
}


//
// TODO 2: ConvertToPlaneCoordinate()
//		Given a plane defined by points, converts their coordinates into
//		a plane coordinate of your choise.
//              See the pdf titled "Homography from Polygon in R^3 to Image Plane",
//              whose link can be found from the project page.
//
//      The final divisors you apply to the u and v coordinates should be saved uScale and vScale
//
void ConvertToPlaneCoordinate(const vector<SVMPoint>& points, vector<Vec3d>& basisPts, double &uScale, double &vScale)
{
  int numPoints = points.size();
  typedef Vec4<double> Vec4d;
  /******** BEGIN TODO ********/
  SVMPoint R = points[0];
  SVMPoint P = points[1];
  Vec4d r(R.X, R.Y, R.Z, R.W);
  Vec4d p(P.X, P.Y, P.Z, P.W);
  Vec4d ex = p - r;
  ex.normalize();
  Vec4d pr = p - r;
  Vec4d q, qr;


  double cos;
  double bestcos = 1.0; // the worst case: angle = 0

  for(int i = 2; i < numPoints; i++){
    Vec4d q_t(points[i].X, points[i].Y, points[i].Z, points[i].W);
    Vec4d qr_t = q_t - r;
    qr_t.normalize();
    
    cos = qr_t * ex;
    cos = abs(cos); // sign is not important
    if(cos < bestcos){
      bestcos = cos;
      q = q_t;
      qr = qr_t;
    }
  }

  Vec4d ey = qr - (qr * ex) * ex;
  ey.normalize();
  double inf = 1e10;
  double uMin = inf, uMax = 0.0, vMin = inf, vMax = 0.0;
  for(int i=0; i< numPoints; i++){
    Vec4d a(points[i].X, points[i].Y, points[i].Z, points[i].W);
    Vec4d ar = a - r;
    double bx = ar * ex;
    double by = ar * ey;
    Vec3d coor(bx, by, 1);

    if(coor[0] > uMax){
      uMax = coor[0];
    }else if(coor[0] < uMin){
      uMin = coor[0];
    }

    if(coor[1] > vMax){
      vMax = coor[1];
    }else if(coor[1] < vMin){
      vMin = coor[1];
    }
    basisPts.push_back(coor);
    uScale = uMax - uMin;
    vScale = vMax - vMax;

  }
  
printf("TODO: svmmath.cpp:101\n"); 
fl_message("TODO: svmmath.cpp:101\n");

	/******** END TODO ********/
}



//
// TODO 3: ComputeHomography()
//		Computes the homography H from the plane specified by "points" to the image plane,
//		and its inverse Hinv.
//		If the plane is the reference plane (isRefPlane == true), don't convert the
//		coordinate system to the plane. Only do this for polygon patches where
//		texture mapping is necessary.
//		Coordinate system conversion is to be implemented in a separate routine
//		ConvertToPlaneCoordinate.
//		For more detailed explaination, see the pdf titled
//              "Homography from Polygon in R^3 to Image Plane", whose link can be found from
//              the project page.
//
void ComputeHomography(CTransform3x3 &H, CTransform3x3 &Hinv, const vector<SVMPoint> &points, vector<Vec3d> &basisPts, bool isRefPlane)
{
  int i;
  int numPoints = (int) points.size();
  assert( numPoints >= 4 );

  basisPts.clear();
  if (isRefPlane) // reference plane
    {
      for (i=0; i < numPoints; i++)
        {
          Vec3d tmp = Vec3d(points[i].X, points[i].Y, points[i].W); // was Z, not W
          basisPts.push_back(tmp);
        }
    }
  else // arbitrary polygon
    {
      double uScale, vScale; // unused in this function
      ConvertToPlaneCoordinate(points, basisPts, uScale, vScale);
    }

  // A: 2n x 9 matrix where n is the number of points on the plane
  //    as discussed in lecture
  int numRows = 2 * numPoints;
  const int numCols = 9;

  typedef Matrix<double, Dynamic, 9, RowMajor> MatrixType;
  MatrixType A = MatrixType::Zero(numRows, numCols);

  /******** BEGIN TODO ********/
  

  for(i=0; i < numPoints; i++){
    // cout<<"Start constructing matrix A..."<<endl;
    SVMPoint tmp = points[i];

    A.block(2*i,0,2,9) << basisPts[i][0], basisPts[i][1], 1, 0, 0, 0, -tmp.u * basisPts[i][0], -tmp.u * basisPts[i][1], -tmp.u,
                          0, 0, 0, basisPts[i][0], basisPts[i][1], 1, -tmp.v * basisPts[i][0], -tmp.v * basisPts[i][1], -tmp.v;
  }

// cout<<"Finish constructing matrix A..."<<endl;
// cout<<"A: "<<A<<endl;

printf("TODO: svmmath.cpp:187\n"); 
fl_message("TODO: svmmath.cpp:187\n");


 double eval, h[9];
 MinEig(A, eval, h);

 H[0][0] = h[0];
 H[0][1] = h[1];
 H[0][2] = h[2];

 H[1][0] = h[3];
 H[1][1] = h[4];
 H[1][2] = h[5];

 H[2][0] = h[6];
 H[2][1] = h[7];
 H[2][2] = h[8];

 /******** END TODO ********/

 // compute inverse of H
 if (H.Determinant() == 0)
   fl_alert("Computed homography matrix is uninvertible \n");
 else
   Hinv = H.Inverse();

 int ii;
 printf("\nH=[\n");
 for (ii=0; ii<3; ii++)
   printf("%e\t%e\t%e;\n", H[ii][0]/H[2][2], H[ii][1]/H[2][2], H[ii][2]/H[2][2]);
 printf("]\nHinv=[\n");

 for (ii=0; ii<3; ii++)
   printf("%e\t%e\t%e;\n", Hinv[ii][0]/Hinv[2][2], Hinv[ii][1]/Hinv[2][2], Hinv[ii][2]/Hinv[2][2]);

 printf("]\n\n");
}
