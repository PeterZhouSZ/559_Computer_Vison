/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #2:
 * ImgView.inl (included from ImgView.cpp)
 *		contains routines for computing the 3D position of points
 ***************************************************************/

//
// TODO 4: sameXY()
//		Computes the 3D position of newPoint using knownPoint
//		that has the same X and Y coordinate, i.e. is directly
//		below or above newPoint.
//		See lecture slide on measuring heights.
//
// HINT1: make sure to dehomogenize points when necessary
// HINT2: there is a degeneracy that you should look out for involving points already in line with the reference
// HINT3: make sure to get the sign of the result right, i.e. whether it is above or below ground
void ImgView::sameXY()
{
  if (pntSelStack.size() < 2)
    {
      fl_alert("Not enough points on the stack.");
      return;
    }
  
  SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
  SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];
  
  if( !knownPoint.known() )
    {
      fl_alert("Can't compute relative values for unknown point.");
      return;
    }
  
  if( refPointOffPlane == NULL )
    {
      fl_alert("Need to specify the reference height first.");
      return;
    }
  
  /******** BEGIN TODO ********/
  
  // See the lecture note on measuring heights
  // using a known point directly below the new point.
  
  // printf("sameXY() to be implemented!\n");
  
  Mat3d HM(H[0][0], H[0][1], H[0][2],
           H[1][0], H[1][1], H[1][2],
           H[2][0], H[2][1], H[2][2]);

  newPoint.X = knownPoint.X;
  newPoint.Y = knownPoint.Y;
  newPoint.W = 1;

  //Initialize points
  Vec3d vx(xVanish.u/xVanish.w, xVanish.v/xVanish.w, 1);
  Vec3d vy(yVanish.u/yVanish.w, yVanish.v/yVanish.w, 1);
  Vec3d vz(zVanish.u/zVanish.w, zVanish.v/zVanish.w, 1);

  Vec3d knownPt(knownPoint.X, knownPoint.Y, 1);
  Vec3d refPt(refPointOffPlane->X, refPointOffPlane->Y, 1);  
  Vec3d b = HM * refPt;
  Vec3d b0 = HM * knownPt; 
  Vec3d t0(newPoint.u/newPoint.w, newPoint.v/newPoint.w, 1);
  Vec3d r(refPointOffPlane->u/refPointOffPlane->w, refPointOffPlane->v/refPointOffPlane->w, 1);
  
  //Initialize intersections of lines
  Vec3d v = cross(cross(b, b0), cross(vx, vy));

  Vec3d t;
  //Consider involving points already in line with the reference
  if(v.iszero()){
    t = t0 - b0 + b;
  }else{
    t = cross(cross(v, t0), cross(r, b));
  }

  t /= t[2];
  b /= b[2];
  b0 /= b0[2];

  double h = sqrt((t - b)*(t - b)*(vz - r)*(vz - r) / ((r - b)*(r - b)) / ((vz - t)*(vz - t)))*referenceHeight;
  Vec3d t0b0 = b0 - t0 ;
  Vec3d vzb = b - vz;
  
  // cout<<"t0:"<<t0[0]<<t0[1]<<endl;
  // cout<<"b0:"<<b0[0]<<b0[1]<<endl;
  // cout<<"vz:"<<vz[0]<<vz[1]<<endl;
  // cout<<"b:"<<b[0]<<b[1]<<endl;
  // cout<<"t0b0:"<<t0b0[0]<<t0b0[1]<<endl;
  // cout<<"vzb:"<<vzb[0]<<vzb[1]<<endl;
  // cout<<t0b0*vzb<<endl;

  if(t0b0 * vzb >= 0){
    newPoint.Z = h;
  }
  else{
    newPoint.Z = -h;
  }

  
printf("TODO: ImgView.inl:49\n"); 
fl_message("TODO: ImgView.inl:49\n");

	/******** END TODO ********/
 
 newPoint.known(true);
 
 printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );
 
 redraw();
}



//
// TODO 5: sameZPlane()
//		Compute the 3D position of newPoint using knownPoint
//		that lies on the same plane and whose 3D position is known.
//		See the man on the box lecture slide.
//		If newPoint is on the reference plane (Z==0), use homography (this->H, or simply H) directly.
//
// HINT: For this function, you will only need to use the three vanishing points and the reference homography 
//       (in addition to the known 3D location of knownPoint, and the 2D location of newPoint)
void ImgView::sameZPlane()
{
  if (pntSelStack.size() < 2)
    {
      fl_alert("Not enough points on the stack.");
      return;
    }
  
  SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
  SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];
  
  if( !knownPoint.known() )
    {
      fl_alert("Can't compute relative values for unknown point.");
      return;
    }
  
  /******** BEGIN TODO ********/

Mat3d HM(H[0][0], H[0][1], H[0][2],
         H[1][0], H[1][1], H[1][2],
         H[2][0], H[2][1], H[2][2]);

Mat3d HinvM(Hinv[0][0], Hinv[0][1], Hinv[0][2],
            Hinv[1][0], Hinv[1][1], Hinv[1][2], 
            Hinv[2][0], Hinv[2][1], Hinv[2][2]);

Vec3d t1(newPoint.u/newPoint.w, newPoint.v/newPoint.w, 1);
Vec3d b1(newPoint.u/newPoint.w, newPoint.v/newPoint.w, 1);

Vec3d t0(newPoint.u/newPoint.w, newPoint.v/newPoint.w, 1);
Vec3d m0(knownPoint.u/knownPoint.w, knownPoint.v/knownPoint.w, 1);

Vec3d vx(xVanish.u/xVanish.w, xVanish.v/xVanish.w, 1);
Vec3d vy(yVanish.u/yVanish.w, yVanish.v/yVanish.w, 1);
Vec3d vz(zVanish.u/zVanish.w, zVanish.v/zVanish.w, 1);

if(knownPoint.Z != 0)
{
  Vec3d knownPt(knownPoint.X, knownPoint.Y, 1);
  Vec3d v = cross(cross(t1, m0), cross(vx, vy));
  Vec3d b0 = HM * knownPt;
  
  if(v.iszero()){
    b1 = cross(b0, cross(t1, vz));
  }
  else{
    b1 = cross(cross(b0, v), cross(t1, vz));
  }
}

b1 /= b1[2];
t0 = HinvM * b1;

newPoint.X = t0[0] / t0[2];
newPoint.Y = t0[1] / t0[2];
newPoint.Z = knownPoint.Z;
newPoint.W = 1;

printf("TODO: ImgView.inl:142\n"); 
fl_message("TODO: ImgView.inl:142\n");

	/******** END TODO ********/
 
 newPoint.known(true);
 
 printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );
 
 redraw();
}

