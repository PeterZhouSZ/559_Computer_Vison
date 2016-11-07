#include <assert.h>
#include <math.h>
#include <FL/Fl.H>
#include <FL/Fl_Image.H>
#include "features.h"
#include "ImageLib/FileIO.h"
#include <iostream>
#include <vector>
#include <algorithm>
#define PI 3.14159265358979323846

// Compute features of an image.
bool computeFeatures(CFloatImage &image, FeatureSet &features, int featureType) {
	// TODO: Instead of calling dummyComputeFeatures, write your own
	// feature computation routines and call them here.
	switch (featureType) {
	case 1:
		dummyComputeFeatures(image, features);
		break;
	case 2:
		ComputeHarrisFeatures(image, features);
		break;
	case 3: 
		ComputeHarrisFeatures_MOPS(image, features);
		break;
	case 4: 
		
		break;
	default:
		return false;
	}

	// This is just to make sure the IDs are assigned in order, because
	// the ID gets used to index into the feature array.
	for (unsigned int i=0; i<features.size(); i++) {
		features[i].id = i+1;
	}

	return true;
}

// Perform a query on the database.  This simply runs matchFeatures on
// each image in the database, and returns the feature set of the best
// matching image.
bool performQuery(const FeatureSet &f, const ImageDatabase &db, int &bestIndex, vector<FeatureMatch> &bestMatches, double &bestScore, int matchType) {
	// Here's a nice low number.
	bestScore = -1e100;

	vector<FeatureMatch> tempMatches;
	double tempScore;

	for (unsigned int i=0; i<db.size(); i++) {
		if (!matchFeatures(f, db[i].features, tempMatches, tempScore, matchType)) {
			return false;
		}

		if (tempScore > bestScore) {
			bestIndex = i;
			bestScore = tempScore;
			bestMatches = tempMatches;
		}
	}

	return true;
}

// Match one feature set with another.
bool matchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore, int matchType) {
	// TODO: We have given you the ssd matching function, you must write your own
	// feature matching function for the ratio test.
	
	printf("\nMatching features.......\n");

	switch (matchType) {
	case 1:
		ssdMatchFeatures(f1, f2, matches, totalScore);
		return true;
	case 2:
		ratioMatchFeatures(f1, f2, matches, totalScore);
		return true;
	default:
		return false;
	}
}

// Evaluate a match using a ground truth homography.  This computes the
// average SSD distance between the matched feature points and
// the actual transformed positions.
double evaluateMatch(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9]) {
	double d = 0;
	int n = 0;

	double xNew;
	double yNew;

    unsigned int num_matches = matches.size();
	for (unsigned int i=0; i<num_matches; i++) {
		int id1 = matches[i].id1;
        int id2 = matches[i].id2;
        applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);
		d += sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
		n++;
	}	

	return d / n;
}

void addRocData(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9],vector<bool> &isMatch,double threshold,double &maxD) {
	double d = 0;

	double xNew;
	double yNew;

    unsigned int num_matches = matches.size();
	for (unsigned int i=0; i<num_matches; i++) {
		int id1 = matches[i].id1;
        int id2 = matches[i].id2;
		applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);

		// Ignore unmatched points.  There might be a better way to
		// handle this.
		d = sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
		if (d<=threshold)
		{
			isMatch.push_back(1);
		}
		else
		{
			isMatch.push_back(0);
		}

		if (matches[i].score>maxD)
			maxD=matches[i].score;
	}	
}

vector<ROCPoint> computeRocCurve(vector<FeatureMatch> &matches,vector<bool> &isMatch,vector<double> &thresholds)
{
	vector<ROCPoint> dataPoints;

	for (int i=0; i < (int)thresholds.size();i++)
	{
		//printf("Checking threshold: %lf.\r\n",thresholds[i]);
		int tp=0;
		int actualCorrect=0;
		int fp=0;
		int actualError=0;
		int total=0;

        int num_matches = (int) matches.size();
		for (int j=0;j < num_matches;j++)
		{
			if (isMatch[j])
			{
				actualCorrect++;
				if (matches[j].score<thresholds[i])
				{
					tp++;
				}
			}
			else
			{
				actualError++;
				if (matches[j].score<thresholds[i])
				{
					fp++;
				}
            }
			
			total++;
		}

		ROCPoint newPoint;
		//printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);
		newPoint.trueRate=(double(tp)/max(actualCorrect,1));
		newPoint.falseRate=(double(fp)/actualError);
		//printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);

		dataPoints.push_back(newPoint);
	}

	return dataPoints;
}


// Compute silly example features.  This doesn't do anything
// meaningful.
void dummyComputeFeatures(CFloatImage &image, FeatureSet &features) {
	CShape sh = image.Shape();
	Feature f;

	for (int y=0; y<sh.height; y++) {
		for (int x=0; x<sh.width; x++) {
			double r = image.Pixel(x,y,0);
			double g = image.Pixel(x,y,1);
			double b = image.Pixel(x,y,2);

			if ((int)(255*(r+g+b)+0.5) % 100  == 1) {
				// If the pixel satisfies this meaningless criterion,
				// make it a feature.
				
				f.type = 1;
				f.id += 1;
				f.x = x;
				f.y = y;

				f.data.resize(1);
				f.data[0] = r + g + b;

				features.push_back(f);
			}
		}
	}
}

void ComputeHarrisFeatures(CFloatImage &image, FeatureSet &features)
{
	//Create grayscale image used for Harris detection
	CFloatImage grayImage=ConvertToGray(image);


	//Create image to store Harris values
	CFloatImage harrisImage(image.Shape().width,image.Shape().height,1);

	//Create image to store local maximum harris values as 1, other pixels 0
	CByteImage harrisMaxImage(image.Shape().width,image.Shape().height,1);

	
	//compute Harris values puts harris values at each pixel position in harrisImage. 
	//You'll need to implement this function.
    computeHarrisValues(grayImage, harrisImage);
	
	// Threshold the harris image and compute local maxima.  You'll need to implement this function.
	computeLocalMaxima(harrisImage,harrisMaxImage);

    // Prints out the harris image for debugging purposes
	CByteImage tmp(harrisImage.Shape());
	convertToByteImage(harrisImage, tmp);
    WriteFile(tmp, "harris.tga");

	// Prints out the harrisMaxImage for debugging purposes
	//CByteImage tmp_max(harrisMaxImage.Shape());
    WriteFile(harrisMaxImage, "harrisMaxImage.tga");
    
    CShape s = harrisMaxImage.Shape();
    int w = s.width;
    int h = s.height;
	
	CFloatImage orientation(w, h, 1);
	computeOrientation(grayImage, orientation);

	
	// TO DO--------------------------------------------------------------------
	//Loop through feature points in harrisMaxImage and create feature descriptor 
	//for each point above a threshold

    for (int y=0;y<h;y++) {
		for (int x=0;x<w;x++) {
		
			// Skip over non-maxima
            if (harrisMaxImage.Pixel(x, y, 0) == 0)
                continue;

            //TO DO---------------------------------------------------------------------
		    // Fill in feature with descriptor data here. 
            Feature f;
			f.type = 2;
			f.id += 1;
			f.x = x;
			f.y = y;
			f.angleRadians=orientation.Pixel(x, y, 0);

            // Add the feature to the list of features
            features.push_back(f);
        }
	}


	ComputeSimpleDescriptors(grayImage, features);
}

void ComputeHarrisFeatures_MOPS(CFloatImage &image, FeatureSet &features)
{
	//Create grayscale image used for Harris detection
	CFloatImage grayImage=ConvertToGray(image);


	//Create image to store Harris values
	CFloatImage harrisImage(image.Shape().width,image.Shape().height,1);

	//Create image to store local maximum harris values as 1, other pixels 0
	CByteImage harrisMaxImage(image.Shape().width,image.Shape().height,1);

	
	//compute Harris values puts harris values at each pixel position in harrisImage. 
	//You'll need to implement this function.
    computeHarrisValues(grayImage, harrisImage);
	
	// Threshold the harris image and compute local maxima.  You'll need to implement this function.
	computeLocalMaxima(harrisImage,harrisMaxImage);

	/*for(int y=0; y<harrisImage.Shape().height; y++)
	{
		for(int x=0; x<harrisImage.Shape().width; x++)
		{
			harrisImage.Pixel(x,y,0)=10*harrisImage.Pixel(x,y,0);
		}
	}*/

    // Prints out the harris image for debugging purposes
	CByteImage tmp(harrisImage.Shape());
	convertToByteImage(harrisImage, tmp);
    WriteFile(tmp, "harris.tga");

	// Prints out the harrisMaxImage for debugging purposes
	//CByteImage tmp_max(harrisMaxImage.Shape());
    WriteFile(harrisMaxImage, "harrisMaxImage.tga");
    
	//Get the orientation for image gradient
	CFloatImage orientation(harrisMaxImage.Shape().width,harrisMaxImage.Shape().height,1);

	computeOrientation(grayImage, orientation);
	//Feature Center
	
	// TO DO--------------------------------------------------------------------
	//Loop through feature points in harrisMaxImage and create feature descriptor 
	//for each point above a threshold

    for (int y=0;y<harrisMaxImage.Shape().height;y++) {
		for (int x=0;x<harrisMaxImage.Shape().width;x++) {
		
			// Skip over non-maxima
            if (harrisMaxImage.Pixel(x, y, 0) == 0)
                continue;

            //TO DO---------------------------------------------------------------------
		    // Fill in feature with descriptor data here. 
            Feature f;
			f.type = 2;
			f.id += 1;
			f.x = x;
			f.y = y;
			f.angleRadians=orientation.Pixel(x,y,0);

			//cout<<"Feature id for level 1 is"<<f.id<<endl;

            // Add the feature to the list of features
            features.push_back(f);
        }
	}
	ComputeMOPSDescriptors(grayImage, features);
}

//TO DO---------------------------------------------------------------------
//Loop through the image to compute the harris corner values as described in class
// srcImage:  grayscale of original image
// harrisImage:  populate the harris values per pixel in this image
void computeHarrisValues(CFloatImage &srcImage, CFloatImage &harrisImage)
{
    CShape s = srcImage.Shape();
	int w = s.width;
    int h = s.height;

    //Ix, Iy
	CFloatImage px = CFloatImage(s);
	CFloatImage py = CFloatImage(s);

	CFloatImage pxX = CFloatImage(s);
	CFloatImage pyY = CFloatImage(s);
	CFloatImage pxY = CFloatImage(s);

	//Below alg is the outline of a basic feature detection alg mentioned in book 4.1

    //1. Compute the horizontal and vertical derivatives of the srcImage Ix and Iy by convolving the original image with
    //   derivatives of Gaussians (According to paper, choose sigma = 1)
	Convolve(srcImage, px, ConvolveKernel_SobelX);
	Convolve(srcImage, py, ConvolveKernel_SobelY);

	//2. Compute the three images corresponding to the outer products of these gradients

	for (int y = 0; y < h; y++) 
	{	for (int x = 0; x < w; x++) 
		{
			double Ix = px.Pixel(x, y, 0);
			double Iy = py.Pixel(x, y, 0);

			pxX.Pixel(x, y, 0) = Ix * Ix;
        	pxY.Pixel(x, y, 0) = Ix * Iy;
        	pyY.Pixel(x, y, 0) = Iy * Iy;
		}
	}
	//3. Convolve each of these images with a larger Gaussian

	// CFloatImage H11 = CFloatImage(s);
 //    //H11.borderMode = eBorderReflect;
 //    //H11.origin[0] = 2;
 //    //H11.origin[1] = 2;
 //    CFloatImage H12 = CFloatImage(s);
 //    //H12.borderMode = eBorderReflect;
 //    //H12.origin[0] = 2;
 //    //H12.origin[1] = 2;
 //    CFloatImage H22 = CFloatImage(s);
 //    //H22.borderMode = eBorderReflect;
 //    //H22.origin[0] = 2;
 //    //H22.origin[1] = 2;


 //    Convolve(h11, H11, ConvolveKernel_Gaussian_5x5_2);
	// Convolve(h12, H12, ConvolveKernel_Gaussian_5x5_2);
	// Convolve(h22, H22, ConvolveKernel_Gaussian_5x5_2);
	

	CFloatImage weight(5,5,1);
	for(int y=0; y<5; y++)
	{
		for(int x=0; x<5; x++)
		{
			weight.Pixel(x,y,0)=gaussian5x5[y*5+x];
		}
	}

	CFloatImage H11(srcImage.Shape());
	CFloatImage H22(srcImage.Shape());
	CFloatImage H12(srcImage.Shape());

	Convolve(pxX, H11, weight);
	Convolve(pyY, H22, weight);
	Convolve(pxY, H12, weight);

	//4. Compute a scalar interest measure using one of the formulas : harmonic mean
	
	// double maxima = 0;
	// double min  = 0;

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            
			// TODO:  Compute the harris score for 'srcImage' at this pixel and store in 'harrisImage'.  See the project
            //   page for pointers on how to do this

   //          double ph11 = H11.Pixel(x, y, 0);
			// double ph12 = H12.Pixel(x, y, 0);
			// double ph22 = H22.Pixel(x, y, 0);

			// double sum = (ph11 * ph22 - ph12 * ph12)/(ph11 + ph22);
			// harrisImage.Pixel(x, y, 0) = sum;

			// if(sum>maxima){
			// 	maxima = sum;
			// }else if(sum <min){
			// 	min = sum;
			// }

			if(H11.Pixel(x,y,0)+H22.Pixel(x,y,0)==0)
				harrisImage.Pixel(x,y,0)=0;
			else
			{
				float f=(H11.Pixel(x,y,0)*H22.Pixel(x,y,0)-H12.Pixel(x,y,0)*H12.Pixel(x,y,0))/(H11.Pixel(x,y,0)+H22.Pixel(x,y,0));
				harrisImage.Pixel(x,y,0)=f;
			}
        }
    }   
	CByteImage tmp(harrisImage.Shape());
    convertToByteImage(harrisImage, tmp);
    WriteFile(tmp, "harrisImage.tga");
}

void computeOrientation(CFloatImage &srcImage, CFloatImage &orientation)
{
	int w = srcImage.Shape().width;
    int h = srcImage.Shape().height;

	CFloatImage px(srcImage.Shape());
	CFloatImage py(srcImage.Shape());

	Convolve(srcImage, px, ConvolveKernel_SobelX);
	Convolve(srcImage, py, ConvolveKernel_SobelY);

	for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
			orientation.Pixel(x,y,0)=atan2(py.Pixel(x,y,0),px.Pixel(x,y,0));
		}
	}
}

// TO DO---------------------------------------------------------------------
// Loop through the harrisImage to threshold and compute the local maxima in a neighborhood
// srcImage:  image with Harris values
// destImage: Assign 1 to a pixel if it is above a threshold and is the local maximum in 3x3 window, 0 otherwise.
//    You'll need to find a good threshold to use.
void computeLocalMaxima(CFloatImage &srcImage,CByteImage &destImage)
{
	int w=srcImage.Shape().width;
	int h=srcImage.Shape().height;
	for (int y = 0; y < h; y++) 
	{
        for (int x = 0; x < w; x++) 
		{
			float center=srcImage.Pixel(x,y,0);
			if(center>=0.01)
				destImage.Pixel(x,y,0)=1;
			else
				destImage.Pixel(x,y,0)=0;
			for(int w_y=0; w_y<5; w_y++)
			{
				for(int w_x=0; w_x<5; w_x++)
				{
					int yy=y+w_y-2;
					int xx=x+w_x-2;
					float value;
					if(yy == y && xx == x)
					{
						continue;
					}
					//handle the critical cases
					if(yy<0 || xx<0 || yy>=h || xx>=w)
					{
						value=center;
					}
					else
					{
						value=srcImage.Pixel(xx, yy, 0);
					}
					if(center<value)
						destImage.Pixel(x,y,0)=0;
				}
			}
		}	
	}

	std::vector<double> distVect;
	CFloatImage distImage(srcImage.Shape());
	//ANMS algorithm:
	for (int y = 0; y < h; y++){
		for (int x = 0; x < w; x++){


			if(destImage.Pixel(x, y, 0)==1){

				double minpoint = 99999999;
				
				for (int yy = 0; yy < h; yy++)
				{
					for (int xx = 0; xx < w; xx++)
					{
						if(destImage.Pixel(xx, yy, 0)==1)
						{
							double out = srcImage.Pixel(x, y, 0);
							double in = srcImage.Pixel(xx, yy, 0);

							if((x != xx)&&(y!=yy)&&(out<0.9*in))
							{
								double dist = sqrt(pow((xx-x),2)+pow((yy-y),2));
								if (dist < minpoint)
								{
									minpoint = dist;

								}
							}
						}

			
					}
				}
			distVect.push_back(minpoint);
			distImage.Pixel(x, y, 0)=minpoint;
			}else{
				distImage.Pixel(x, y, 0)=-1;
			}
			
		}
	}
	CByteImage tmppp(distImage.Shape());
    convertToByteImage(distImage, tmppp);
    WriteFile(tmppp, "distImage.tga");

    int size = distVect.size();

	std::nth_element(distVect.begin(), distVect.begin()+700, distVect.end(), std::greater<double>());
	double n500th=distVect[700];
	double nmin = distVect[1];
	cout<<nmin<<std::endl;
	cout <<n500th<<std::endl;
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			double center = distImage.Pixel(x, y, 0);
			if (center >= n500th){
				destImage.Pixel(x, y, 0)=1;
			}else{
				destImage.Pixel(x, y, 0)=0;
			}
		}
	}
}

// Perform simple feature matching.  This just uses the SSD
// distance between two feature vectors, and matches a feature in the
// first image with the closest feature in the second image.  It can
// match multiple features in the first image to the same feature in
// the second image.
void ssdMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) {
	int m = f1.size();
	int n = f2.size();

	matches.resize(m);
	totalScore = 0;

	double d;
	double dBest;
	int idBest;

	for (int i=0; i<m; i++) {
		dBest = 1e100;
		idBest = 0;

		for (int j=0; j<n; j++) {
			d = distanceSSD(f1[i].data, f2[j].data);

			if (d < dBest) {
				dBest = d;
				idBest = f2[j].id;
			}
		}
        matches[i].id1 = f1[i].id;
		matches[i].id2 = idBest;
		matches[i].score = dBest;
		totalScore += matches[i].score;
	}
}

// TODO: Write this function to perform ratio feature matching.  
// This just uses the ratio of the SSD distance of the two best matches as the score
// and matches a feature in the first image with the closest feature in the second image.
// It can match multiple features in the first image to the same feature in
// the second image.  (See class notes for more information, and the sshMatchFeatures function above as a reference)
void ratioMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) {
	int m = f1.size();
	int n = f2.size();

	//matches.resize(m);
	totalScore = 0;

	double d;
	double dBest;
	double dSecond;
	int idBest;
	int idSecond;

	for (int i=0; i<m; i++) {
		dBest = 1e100;
		idBest = 0;
		dSecond = 1e100;
		idSecond = 0;

		for (int j=0; j<n; j++) {
			d = distanceSSD(f1[i].data, f2[j].data);

			if (d < dBest) {
				dSecond = dBest;
				idSecond = idBest;
				dBest = d;
				idBest = f2[j].id;
			}
			else if (d < dSecond && d >= dBest) {
				dSecond = d;
				idSecond = f2[j].id;
			}
		}

			FeatureMatch goodFeature;
			goodFeature.id1 = f1[i].id;
			goodFeature.id2 = idBest;
			goodFeature.score = dBest/dSecond;
			matches.push_back(goodFeature);
		//cout<<"match "<<matches[i].id1<<" with "<<matches[i].id2<<endl;
	}
    
}


// Convert Fl_Image to CFloatImage.
bool convertImage(const Fl_Image *image, CFloatImage &convertedImage) {
	if (image == NULL) {
		return false;
	}

	// Let's not handle indexed color images.
	if (image->count() != 1) {
		return false;
	}

	int w = image->w();
	int h = image->h();
	int d = image->d();

	// Get the image data.
	const char *const *data = image->data();

	int index = 0;

	for (int y=0; y<h; y++) {
		for (int x=0; x<w; x++) {
			if (d < 3) {
				// If there are fewer than 3 channels, just use the
				// first one for all colors.
				convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,1) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,2) = ((uchar) data[0][index]) / 255.0f;
			}
			else {
				// Otherwise, use the first 3.
				convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
				convertedImage.Pixel(x,y,1) = ((uchar) data[0][index+1]) / 255.0f;
				convertedImage.Pixel(x,y,2) = ((uchar) data[0][index+2]) / 255.0f;
			}

			index += d;
		}
	}
	
	return true;
}

// Convert CFloatImage to CByteImage.
void convertToByteImage(CFloatImage &floatImage, CByteImage &byteImage) {
	CShape sh = floatImage.Shape();

    assert(floatImage.Shape().nBands == byteImage.Shape().nBands);
	for (int y=0; y<sh.height; y++) {
		for (int x=0; x<sh.width; x++) {
			for (int c=0; c<sh.nBands; c++) {
				float value = floor(255*floatImage.Pixel(x,y,c) + 0.5f);

				if (value < byteImage.MinVal()) {
					value = byteImage.MinVal();
				}
				else if (value > byteImage.MaxVal()) {
					value = byteImage.MaxVal();
				}

				// We have to flip the image and reverse the color
				// channels to get it to come out right.  How silly!
				byteImage.Pixel(x,sh.height-y-1,sh.nBands-c-1) = (uchar) value;
			}
		}
	}
}

// Compute SSD distance between two vectors.
double distanceSSD(const vector<double> &v1, const vector<double> &v2) {
	int m = v1.size();
	int n = v2.size();

	if (m != n) {
		// Here's a big number.
		return 1e100;
	}

	double dist = 0;

	for (int i=0; i<m; i++) {
		dist += pow(v1[i]-v2[i], 2);
	}

	return sqrt(dist);
}

// Transform point by homography.
void applyHomography(double x, double y, double &xNew, double &yNew, double h[9]) {
	double d = h[6]*x + h[7]*y + h[8];

	xNew = (h[0]*x + h[1]*y + h[2]) / d;
	yNew = (h[3]*x + h[4]*y + h[5]) / d;
}

// Compute AUC given a ROC curve
double computeAUC(vector<ROCPoint> &results)
{
	double auc=0;
	double xdiff,ydiff;
	for (int i = 1; i < (int) results.size(); i++)
    {
        //fprintf(stream,"%lf\t%lf\t%lf\n",thresholdList[i],results[i].falseRate,results[i].trueRate);
		xdiff=(results[i].falseRate-results[i-1].falseRate);
		ydiff=(results[i].trueRate-results[i-1].trueRate);
		auc=auc+xdiff*results[i-1].trueRate+xdiff*ydiff/2;
    }
	cout<<"auc for this img is: "<<auc<<endl;
	return auc;
}

// Compute Simple descriptors.
void ComputeSimpleDescriptors(CFloatImage &grayImage, FeatureSet &features)
{

    const int windowSize = 5;
    CFloatImage destImage(windowSize, windowSize, 1);

    for (vector<Feature>::iterator i = features.begin(); i != features.end(); i++) {
        Feature &f = *i;
        int x = f.x;
        int y = f.y;
        f.data.resize(25);

		double sum=0.0, mean=0.0, variance=0.0;

		for(int w=0; w<5; w++)
		{
			for(int h=0; h<5; h++)
			{
				if(x+w-2==x && y+h-2==y)
					continue;
				if(x+w-2<0 || y+h-2<0 || x+w-2>=grayImage.Shape().width || y+h-2>=grayImage.Shape().height)
					sum+=grayImage.Pixel(x,y,0);
				else
					sum+=grayImage.Pixel(x+w-2,y+h-2,0);
			}
		}
		mean=sum/25;
		
		for(int w=0; w<5; w++)
		{
			for(int h=0; h<5; h++)
			{
				if(x+w-2==x && y+h-2==y)
					continue;
				if(x+w-2<0 || y+h-2<0 || x+w-2>=grayImage.Shape().width || y+h-2>=grayImage.Shape().height)
					variance+=pow(grayImage.Pixel(x,y,0)-mean,0);
				else
					variance+=pow(grayImage.Pixel(x+w-2,y+h-2,0)-mean,0);
			}
		}
        
		double stddev=sqrt(variance/25);

		for(int w=0; w<5; w++)
		{
			for(int h=0; h<5; h++)
			{
				if(x+w-2==x && y+h-2==y)
					continue;
				if(x+w-2<0 || y+h-2<0 || x+w-2>=grayImage.Shape().width || y+h-2>=grayImage.Shape().height)
					f.data[w*5+h]=(grayImage.Pixel(x,y,0)-mean)/stddev;
				else
					f.data[w*5+h]=(grayImage.Pixel(x+w-2,y+h-2,0)-mean)/stddev;	
			}
		}
    }
}

void ComputeMOPSDescriptors(CFloatImage &grayImage, FeatureSet &features)
{      
	const int windowSize = 8;
    CFloatImage destImage(windowSize, windowSize, 1);
    CFloatImage blurImg(grayImage.Shape().width, grayImage.Shape().height, 1);
    CFloatImage blurG(5, 5, 1);

	for(int y=0; y<5; y++)
	{
		for(int x=0; x<5; x++)
		{
			blurG.Pixel(x,y,0)=gaussian5x5[y*5+x];
		}
	}

    Convolve(grayImage, blurImg, blurG);
    vector<Feature>::iterator i;

    for (i = features.begin(); i != features.end(); i++) {
        Feature &f = *i;
        CTransform3x3 xform, origin, rotate, trans, scale;

        origin = origin.Translation(float(-windowSize/2), float(-windowSize/2));
        rotate = rotate.Rotation(f.angleRadians * 180.0 / PI); 
        trans = trans.Translation(f.x, f.y);
        scale[0][0] = 40/windowSize;
        scale[1][1] = 40/windowSize;
        xform = trans * rotate * scale * origin;
        WarpGlobal(blurImg, destImage, xform, eWarpInterpLinear);
        f.data.resize(windowSize * windowSize);

        // Normalize
        double mean = 0.0, variance = 0.0;

        for (int i = 0; i < windowSize; i++ ){
            for (int j = 0; j < windowSize; j++){
                mean += destImage.Pixel(j, i, 0);
            }
        }
        mean /= (windowSize * windowSize);
       
        for (int j = 0; j < windowSize; j++ ){
            for (int i = 0; i < windowSize; i++){
                variance += pow((destImage.Pixel(i, j, 0) - mean), 2);
            }
        }
        double stddev = sqrt(variance / ((windowSize * windowSize) - 1));

        for (int j = 0; j < windowSize; j++ ){
            for (int i = 0; i < windowSize; i++){
                f.data[j * windowSize + i] = (destImage.Pixel(i, j, 0) - mean)/stddev;
            }
        }
    }
}

double round(double x) { return (floor(x + 0.5)); }