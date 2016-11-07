
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: eigfaces.cpp                                                                         //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"
using namespace std;

EigFaces::EigFaces()
:
Faces()
{
	//empty
}

EigFaces::EigFaces(int count, int width, int height)
:
Faces(count, width, height)
{
	//empty
}

void EigFaces::projectFace(const Face& face, Vector& coefficients) const
{
	if (face.getWidth()!=width || face.getHeight()!=height) {
		throw Error("Project: Face to project has different dimensions");
	}

	coefficients.resize(getSize());
	// ----------- TODO #2: compute the coefficients for the face and store in coefficients.
	int num = getSize(); 
	for(int i=0; i<num; i++){
		Face subface;
		Face average_face=(*this).getAverage();
		face.sub(average_face, subface);
		coefficients[i]=(*this)[i].dot(subface);
	}
}

void EigFaces::constructFace(const Vector& coefficients, Face& result) const
{	
	// ----------- TODO #3: construct a face given the coefficients
	int num = getSize();
	result=(*this).getAverage();
	for(int i=0; i<num; i++){
		Vector tmp(vector_size);
		tmp = (*this)[i];
		tmp *= coefficients[i];
		result += tmp;
	}
}

bool EigFaces::isFace(const Face& face, double max_reconstructed_mse, double& mse) const
{
	// ----------- TODO #4: Determine if an image is a face and return true if it is. Return the actual
	// MSE you calculated for the determination in mse
	// Be sure to test this method out with some face images and some non face images
	// to verify it is working correctly.
	Vector co;
	projectFace(face, co);
	Face before = face;
	Face after;
	constructFace(co, after);
	mse = after.mse(before);
	if(mse > max_reconstructed_mse)
		return false;
	else
		return true;
}

bool EigFaces::verifyFace(const Face& face, const Vector& user_coefficients, double max_coefficients_mse, double& mse) const
{
	// ----------- TODO #5 : Determine if face is the same user give the user's coefficients.
	// return the MSE you calculated for the determination in mse.
	Vector co;
	projectFace(face, co);
	mse = co.mse(user_coefficients);

	if(mse > max_coefficients_mse)
		return false;
	else 
		return true;

}

void EigFaces::recognizeFace(const Face& face, Users& users) const
{
	// ----------- TODO #6: Sort the users by closeness of match to the face
	Vector co;
	projectFace(face, co);
	for(int i=0; i<users.getSize(); i++){
		users[i].setMse(co.mse(users[i]));
	}
	users.sort();
}

bool compare_error(const FacePosition& first, const FacePosition& second){

	return(first.error<second.error);

}


void EigFaces::findFace(const Image& img, double min_scale, double max_scale, double step, int n, bool crop, Image& result) const
{
	// ----------- TODO #7: Find the faces in Image. Search image scales from min_scale to max_scale inclusive,
	// stepping by step in between. Find the best n faces that do not overlap each other. If crop is true,
	// n is one and you should return the cropped original img in result. The result must be identical
	// to the original besides being cropped. It cannot be scaled and it must be full color. If crop is
	// false, draw green boxes (use r=100, g=255, b=100) around the n faces found. The result must be
	// identical to the original image except for the addition of the boxes.
	
	std::list<FacePosition> lsFP; // keep a list of the top best face positions
	
	for(double k = min_scale; k<=max_scale; k += step){  
		int imgW = img.getWidth();
		int imgH = img.getHeight();
		int imgC = img.getColors();
		int samW = imgW * k; // down-sample size
		int samH = imgH * k;

		Image sampleImg(samW, samH, imgC); //  
		img.resample(sampleImg); // 

		Face f(width, height);

		for(int j=0; j< samH - height; j++){
			for(int i=0; i<samW - width; i++){
				//cout<<"lsFP size is : "<<lsFP.size()<<endl;
				//get face in (i,j), subimage size must be the same as face size
				f.subimage(i, i+width-1, j, j+height-1, sampleImg, false);

				double mse;
				isFace(f, 800, mse);
				//cout<<"start checking isFace..."<<endl;
				Face normF(width, height);
				f.sub((*this).getAverage(), normF);
				double error = mse*normF.mag()/f.var();
				// double r = f.pixel(i+width/2,j+height/2,0);
				// double g= f.pixel(i+width/2,j+height/2,1);
				// double b = f.pixel(i+width/2,j+height/2,2);


				if(isFace(f, 800, mse)){
					//cout<<"This face is valid: mse="<<mse<<endl;

					FacePosition fp;
					fp.x = i;
					fp.y = j;
					fp.scale = k;
					fp.mse = mse;


					Face normF(width, height); // face img after normarization
					f.sub((*this).getAverage(), normF); // normarize face f

					// To prevent the alg from getting fooled by low-texture areas or areas close to face space 
					// but far from the facial mean. The approach is to multiply the MSE by the distance of the 
					// face from the face mean and then dividing by the variance of the face
					fp.error = mse * normF.mag()/f.var(); //mag(): distance of face to average_face
					//cout<<"fp.error:" <<fp.error<<endl;
					// Position is found
					
					// First check if overlap
					std::list<FacePosition>::iterator del = lsFP.end(); 
					bool isOverlap = false; // flag
					bool isReplace = false;
					//cout<<"start finding overlap..."<<endl;

					for(std::list<FacePosition>::iterator t=lsFP.begin(); t!=lsFP.end();t++){
						FacePosition compareFP= *t;
						// The following coordinates are all based on oringinal image size!
						//cout<<"lsFP size"<<lsFP.size()<<endl;
						// cout<<"stuck in here"<<lsFP.size()<<i<<j<<endl;
						
						if(!isReplace){
							int compW = width / compareFP.scale;
							int compH = height / compareFP.scale;
							int compX0 = compareFP.x/compareFP.scale;
							int compY0 = compareFP.y/compareFP.scale;
							int compX1 = compX0 + compW;
							int compY1 = compY0 + compH;
							int W = width / k;
							int H = height / k;
							int X0 = i / k;
							int Y0 = j / k;
							int X1 = X0 + W;
							int Y1 = Y0 + H;
							// Handle the situation when there is already a face position in the list 
							// that overlaps the one you just found.
							int overlap = 20;

							if( (abs(compX0-X0)<=overlap) && (abs(compY0-Y0)<=overlap) && (abs(compX1-X1)<=overlap) && (abs(compY1-Y1)<=overlap) ) 
							{
								//cout<<"find overlap!!!"<<endl;
								isOverlap = true;
								//choose the one with the better face and throw the other away
								if(compareFP.error>fp.error){
								// if(compareFP.error<mse){
									//cout<<"will be replaced!!!"<<endl;
									del = t;
									isReplace = true;
								}
								else{
									//cout<<"will not be replaced!!!"<<endl;
									
								}
							}
						}

					}// END of finding overlap

					//cout<<"quit finding overlap..."<<endl;

					// for(std::list<FacePosition>::iterator it=lsFP.begin();it!=lsFP.end();it++){
					// 	std::cout << "subFace from["<<(int)(*it).x<<","<<(int)(*it).y<<"] has mse: " << (*it).error << std::endl;
					// }

					// Check if Position is better than the last position in the list
					if(!isOverlap){
						//cout<<"start comparing with the last position..."<<endl;
						if(lsFP.size()>=n){
							std::list<FacePosition>::iterator lastFP= --lsFP.end();
							if(((*lastFP).error>fp.error)){
							// if((*lastFP).error<mse){
								lsFP.pop_back(); // delete the last element of the list
								lsFP.push_back(fp);
								cout<<"error is: "<<fp.error<<endl;
								cout<<"mse is: "<<mse<<endl;
								

							}
						}else{
							lsFP.push_back(fp);
							cout<<"error is: "<<fp.error<<endl;
							cout<<"mse is: "<<mse<<endl;
							

						}
						//cout<<"quit comparing with the last position..."<<endl;
					}// END of comparing with the last element
					
					
					if(isReplace){
						//cout<<"start checking replacement..."<<endl;

						lsFP.erase(del); // delete the overlapped one
						//lsFP.insert(del, fp); // insert the new one
						lsFP.push_back(fp);
						cout<<"error is: "<<fp.error<<endl;
						cout<<"mse is: "<<mse<<endl;
						

						//for(std::list<FacePosition>::iterator it=lsFP.begin();it!=lsFP.end();it++){
						// std::cout << "subFace from["<<(int)(*it).x<<","<<(int)(*it).y<<"] has mse: " << (*it).error << std::endl;
						// }
						//cout<<"quit checking replacement..."<<endl;
					}
					
				}// END of a face in (i, j)

				//cout<<"handle a face in ["<< i << ", "<< j<<"]"<<std::endl;
				lsFP.sort(compare_error); // sort list after a face

			}
		}// END of all (i, j)

	}// END of all scale

	//cout<<"Loop through all scale..."<<std::endl;

	// for(std::list<FacePosition>::iterator it=lsFP.begin();it!=lsFP.end();it++){
	// 	std::cout << "subFace from["<<(int)((*it).x/(*it).scale)<<","<<(int)((*it).y/(*it).scale)<<"] has mse: " << (*it).error << std::endl;
	// }

	if(crop){
		FacePosition cropFP = *lsFP.begin();
		int x0 = cropFP.x/cropFP.scale;
		int x1 = (cropFP.x+width)/cropFP.scale;
		int y0 = cropFP.y/cropFP.scale;
		int y1 = (cropFP.y+height)/cropFP.scale;
		img.crop(x0, y0, x1, y1, result);
		std::cout<<"scale is "<<cropFP.scale<<std::endl;
		std::cout<<"origin x is "<<x0<<std::endl;
		std::cout<<"origin y is "<<y0<<std::endl;
		std::cout<<"width is "<<width<<std::endl;
		std::cout<<"height is "<<height<<std::endl;
		std::cout<<"crop x is "<<x1<<std::endl;
		std::cout<<"crop y is "<<y1<<std::endl;
	}else{
		result = img;
		for(std::list<FacePosition>::iterator t=lsFP.begin(); t!=lsFP.end(); t++){
			FacePosition fp= *t;
			int x0 = fp.x/fp.scale;
			int x1 = (fp.x+width)/fp.scale;
			int y0 = fp.y/fp.scale;
			int y1 = (fp.y+height)/fp.scale;
			result.line(x0, y0, x1, y0, 100, 255, 100);
			result.line(x1, y0, x1, y1, 100, 255, 100);
			result.line(x0, y1, x1, y1, 100, 255, 100);
			result.line(x0, y0, x0, y1, 100, 255, 100);
		}
		
	}

}

void EigFaces::morphFaces(const Face& face1, const Face& face2, double distance, Face& result) const
{
	// TODO (extra credit): MORPH along *distance* fraction of the vector from face1 to face2 by
	// interpolating between the coefficients for the two faces and reconstructing the result.
	// For example, distance 0.0 will approximate the first, while distance 1.0 will approximate the second.
	// Negative distances are ok two.

}

const Face& EigFaces::getAverage() const
{
	return average_face;
}

void EigFaces::setAverage(const Face& average)
{
	average_face=average;
}



