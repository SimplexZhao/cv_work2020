#include"localtypes.h"
#include"opencv2/opencv.hpp"
using namespace cv;
#include<iostream>
using namespace std;
#include"getCPts.hpp"
#include"m_solver.hpp"


bool convert_type(vector<vector<RotatedRect>> ells, vector<mPoint2d>& Ipts, vector<vector<int>> code_list)
{
    Ipts.resize(ells.size()*4);
    for(int i=0;i<ells.size();++i)
    {
        for(int j=0;j<4;++j)
        {
            Ipts[i*4+j].x=ells[i][j].center.x;
            Ipts[i*4+j].y=ells[i][j].center.y;
            Ipts[i*4+j].id=code_list[i][j];
        }
    }
    return true;
}

// process the img and get image point list.
void procImg(string imgname, vector<mPoint2d>& IPt_list);

int main(int argc, char** argv)
{
    // string imgname;
    // cin>>imgname;
    if(argc!=5){cout<<"Error in loading files."<<endl;return -1;}
    string imagename_L=argv[1];//"./编程实习材料/无畸变影像-左.bmp";
    string imagename_R=argv[2];//"./编程实习材料/无畸变影像-右.bmp";
    IOPs iops;
    EOPs eop_left,eop_right;
// From srcimag to intersection main begin.
    vector<mPoint3d> CP_list;
    solver::loadIn_Params(argv[3]/*"./编程实习材料/内部参数.txt"*/, iops);
    solver::loadControlfile(argv[4]/*"./编程实习材料/控制点坐标.txt"*/, CP_list);
    cout<<"Files loaded."<<endl;
    vector<mPoint2d> IPt_list_Left, IPt_list_Right;
    procImg(imagename_L, IPt_list_Left);
    procImg(imagename_R, IPt_list_Right);
    solver::getEO(iops, eop_left, IPt_list_Left, CP_list);
    solver::nonlinear_EO(iops, eop_left, IPt_list_Left, CP_list,"left");
    solver::getEO(iops, eop_right, IPt_list_Right, CP_list);
    solver::nonlinear_EO(iops, eop_right, IPt_list_Right, CP_list,"right");
    solver::writeEOPfile(eop_left, eop_right);
    solver::DualImg_Intersection(iops, eop_left, eop_right, IPt_list_Left, IPt_list_Right,CP_list);
    cout<<"EOPs computation over, report written."<<endl;
// main end.
// Rectify begin
    solver::loadIn_Params(argv[3]/*"./编程实习材料/内部参数.txt"*/, iops);
    solver::loadEOPfile(eop_left, eop_right);
    Mat leftimg=imread(imagename_L);
    Mat rightimg=imread(imagename_R);
    Mat KCam=(Mat_<double>(3,3)<<iops.f,0,iops.x0,0,iops.f,iops.y0,0,0,1);
    double Xpi[9]; MatOperation::GetRotationMat(Xpi, 0.0, 3.14159265358979, 0.0,2);
    double RL_op[9]; MatOperation::GetRotationMat(RL_op, eop_left.phi, eop_left.omega, eop_left.kappa,2);
    double RL[9]; MatOperation::MatrixMulti(Xpi, RL_op, RL,3,3,3);
    double RLt[9]; MatOperation::Rotate(RL, RLt, 3,3);
    double RR_op[9]; MatOperation::GetRotationMat(RR_op, eop_right.phi, eop_right.omega, eop_right.kappa,2);
    double RR[9]; MatOperation::MatrixMulti(Xpi, RR_op, RR,3,3,3);
    double RRL[9]; MatOperation::MatrixMulti(RR, RLt, RRL, 3,3,3);// RRL=RR*RL^(-1);
    //double Ttmp[3]={eop_right.offX-eop_left.offX, eop_right.offY-eop_left.offY, eop_right.offZ-eop_left.offZ}; // wrong translation vector
    double Ttmp[3]={-eop_right.offX+eop_left.offX, -eop_right.offY+eop_left.offY, -eop_right.offZ+eop_left.offZ};
    double Ttmp2[3];
    Mat R, T;
    R=(Mat_<double>(3,3)<<RRL[0],RRL[1],RRL[2],RRL[3],RRL[4],RRL[5],RRL[6],RRL[7],RRL[8]);
    MatOperation::MatrixMulti(RR, Ttmp, Ttmp2,3,1,3);
    
    T=(Mat_<double>(3,1)<<Ttmp2[0], Ttmp2[1], Ttmp2[2]);
    Mat R1, R2, P1, P2, Q;
    cout<<"Rectify begin."<<endl;
    stereoRectify(KCam,Mat::zeros(1,5,CV_32F),KCam,Mat::zeros(1,5,CV_32F), leftimg.size(),R,T,R1,R2,P1,P2,Q,0);
    // cout<<P1<<endl<<P2<<endl;// 输出相机矩阵
    //计算映射
    Mat rmap[2][2];
    initUndistortRectifyMap(KCam, Mat::zeros(1,5,CV_32F), R1, P1, leftimg.size(), CV_32FC1, rmap[0][0], rmap[0][1]);
    initUndistortRectifyMap(KCam, Mat::zeros(1,5,CV_32F), R2, P2, leftimg.size(), CV_32FC1, rmap[1][0], rmap[1][1]);
    Mat imgLr, imgRr;
    remap(leftimg, imgLr, rmap[0][0], rmap[0][1], INTER_AREA);//左校正
    remap(rightimg, imgRr, rmap[1][0], rmap[1][1], INTER_AREA);//右校正
    // save epipolar image pair.
	imwrite("imgLr.jpg",imgLr);
	imwrite("imgRr.jpg",imgRr);
    cout<<"Rectify over."<<endl;
    // Mat showImage(leftimg.size().height,2*leftimg.size().width,CV_8UC3);
    // Rect rectLeft(0,0,leftimg.size().width,leftimg.size().height);
    // Rect rectRight(leftimg.size().width,0,leftimg.size().width,leftimg.size().height);
    // imgLr.copyTo(showImage(rectLeft));
    // imgRr.copyTo(showImage(rectRight));
    // imwrite("compare.bmp", showImage);
// Rectify end
// Get diaparity map begin
    // struct _stat buffer[1];
    //int flag_dispmap=_stat("./vocabulary.xml",buffer);
    // if(flag_dispmap==-1)
    // {
        cout<<"Dense matching begin."<<endl;
        Ptr<StereoSGBM> sgbm = StereoSGBM::create(-160,640,3);
        int sgbmWinSize = 3;
        int cn = imgLr.channels();
        int numDisparities = 480;
        sgbm->setPreFilterCap(63);
        // sgbm->setBlockSize(sgbmWinSize);
        sgbm->setP1(8*cn*sgbmWinSize*sgbmWinSize);
        sgbm->setP2(32*cn*sgbmWinSize*sgbmWinSize);
        sgbm->setUniquenessRatio(10);
        sgbm->setSpeckleWindowSize(100);
        sgbm->setSpeckleRange(2);
        sgbm->setDisp12MaxDiff(20);
        sgbm->setMode(StereoSGBM::MODE_SGBM);
        cv::Mat disp, disp8(imgLr.size(),CV_32FC1);
        sgbm->compute(imgLr, imgRr, disp);
        disp.convertTo(disp, CV_32FC1, 1/16.0);
        normalize(disp,disp8,1.0,0.0,NORM_MINMAX);
        Mat show;
        disp8.convertTo(show, CV_8U, 255.0);
        cv::imwrite("disp.jpg", show);
        cout<<"Dense matching over."<<endl;
        // FileStorage fs("./vocabulary.xml", FileStorage::WRITE);
        // fs<<"dispmap"<<disp; fs.release();
    // }
// Get diaparity map end
// Compute 3d points start
    Mat M1,M2;
    M1=Mat::zeros(4,4,CV_64F);
    M2=Mat::zeros(4,4,CV_64F);
    R1.convertTo(M1(Range(0, 3), Range(0, 3)),CV_64F);
    R2.convertTo(M2(Range(0, 3), Range(0, 3)),CV_64F);
    M1.ptr<float>(2)[3]=1.0; M2.ptr<float>(2)[3]=1.0;
    // load disparity file
    // FileStorage ifs(".\\vocabulary.xml", FileStorage::READ);
    // Mat disp_in;
    // ifs["dispmap"] >> disp_in;
    // ifs.release();
    // load over.
    Mat disp_in=disp;
    int count=0;
    medianBlur(disp_in, disp_in,13);
    vector<Point2f> leftpoints, rightpoints;
    for(int j=0;j<disp_in.rows;++j)
    {
        float *k=disp_in.ptr<float>(j);
        for(int i=0;i<disp_in.cols;++i)
        {
            if(k[i]>-161&&k[i]<480)
            {
                leftpoints.push_back(Point2f(i,j));
                rightpoints.push_back(Point2f(i-k[i],j));
                count++;
            }
        }
    }
    Mat structure;
    cout<<"TriangulatePoints start."<<endl;
    triangulatePoints(P1, P2, leftpoints, rightpoints, structure);
    cout<<"TriangulatePoints over."<<endl;
    cout<<"Write OBJ file start."<<endl;
    ofstream objfile("重建结果.obj");
    float* xptr=structure.ptr<float>(0);
    float* yptr=structure.ptr<float>(1);
    float* zptr=structure.ptr<float>(2);
    float* eptr=structure.ptr<float>(3);
    for(int i=0;i<structure.cols;i++)
    objfile<<"v"<<" "<<xptr[i]/eptr[i]<<" "<<yptr[i]/eptr[i]<<" "<<zptr[i]/eptr[i]<<endl;
    objfile.close();
    cout<<"Write OBJ file over."<<endl;
// Compute end and write file end.

    return 0;
}

void procImg(string imgname, vector<mPoint2d>& IPt_list)
{
    Mat input=imread(imgname,0);
    Mat edge_out(input.rows, input.cols, CV_8UC1, Scalar(0));
    Canny(input, edge_out, 128, 255);
    vector<RotatedRect> ells;
    vector<Vec2f> ls;
    getCPts::EclipseExtract(edge_out, ells, ls);
    vector<vector<RotatedRect>> Ell_group;
    getCPts::Points2Group(ells, ls, Ell_group);
    vector<vector<int>> code_list;
    getCPts::EncodePts(Ell_group, code_list);
    convert_type(Ell_group, IPt_list, code_list);
}

// SIFT match begin
    // Ptr<SIFT> sift=SIFT::create(30);
    // vector<KeyPoint> keypoints1,keypoints2;
    // sift->detect(leftimg, keypoints1);
    // sift->detect(rightimg, keypoints2);
    // Mat result1, result2;
    // drawKeypoints(leftimg, keypoints1, result1, Scalar::all(-1), DrawMatchesFlags::DEFAULT);
    // drawKeypoints(rightimg, keypoints2, result2, Scalar::all(-1), DrawMatchesFlags::DEFAULT);
    // Mat descriptors1, descriptors2;
    // sift->compute(leftimg, keypoints1, descriptors1);
    // sift->compute(rightimg, keypoints2, descriptors2);
    // Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("BruteForce");
    // vector<DMatch> matches;
    // matcher->match(descriptors1, descriptors2, matches);
    // Mat imgMatches;
    // drawMatches(leftimg, keypoints1, rightimg, keypoints2, matches, imgMatches);
    // imwrite("output.bmp", imgMatches);
// SIFT match end

    // vector<uchar> RansacStatus;
	// Mat Fundamental = findHomography(leftpoints, rightpoints, RansacStatus, LMEDS, 3);
    // Fundamental.release();
    // vector<Point2f> leftpts_filtered, rightpts_filtered;
    // cout<<endl<<"RANSAC start."<<endl;
    // for(int i=0;i<leftpoints.size();++i)
    // {
    //     if (RansacStatus[i] != 0)
	// 	{
	// 		leftpts_filtered.push_back(leftpoints[i]);
	// 		rightpts_filtered.push_back(rightpoints[i]);
	// 	}
    // }
    // cout<<endl<<"RANSAC over."<<endl;