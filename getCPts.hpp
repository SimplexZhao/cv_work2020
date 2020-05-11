#include<iostream>
#include<fstream>
using namespace std;
typedef unsigned char uchar;
#include"opencv2/opencv.hpp"
using namespace cv;

bool drawlines(vector<Vec2f> lines, Mat image)
{
    for (size_t i = 0; i < lines.size(); i++)
	{
        
		float rho = lines[i][0], theta = lines[i][1];
		Point pt1, pt2;
		double a = cos(theta), b = sin(theta);
		double x0 = a*rho, y0 = b*rho;
		pt1.x = cvRound(x0 + 3000 * (-b));  //把浮点数转化成整数
		pt1.y = cvRound(y0 + 3000 * (a));
		pt2.x = cvRound(x0 - 3000 * (-b));
		pt2.y = cvRound(y0 - 3000 * (a));
		line(image, pt1, pt2, Scalar(255), 1);
	}
    return true;
}


class getCPts
{
    public:
        getCPts();
        virtual ~getCPts();
        static bool EclipseExtract(Mat input, vector<RotatedRect>& ell, vector<Vec2f>& lines); // 提取椭圆（圆心）
        static bool Points2Group(vector<RotatedRect> ells, vector<Vec2f> lines, vector<vector<RotatedRect>>& ell_group);
        static bool EncodePts(vector<vector<RotatedRect>> ell_group, vector<vector<int>>& seq_list); // 获取控制点编号如‘0001’
        static bool ShowCodes(Mat& img, vector<vector<RotatedRect>> ell_group, vector<vector<int>> seq_list);
};


Point2f getCrossPoint(Vec2f lA, Vec2f lB)// get crosspoint of lA/lB in pole coordinates
{
    float thetaA = lA[1], thetaB = lB[1];
    float rhoA = lA[0], rhoB = lB[0];
    double k_a = -1.0/tan(thetaA), k_b = -1.0/tan(thetaB);
    double b_a= rhoA/sin(thetaA), b_b= rhoB/sin(thetaB);
    Point2f crossPoint;
    crossPoint.x = (b_b - b_a) / (k_a - k_b);
    crossPoint.y = k_a*crossPoint.x+b_a;
    return crossPoint;
}
//
Point2f ctr_pt;
bool closer(RotatedRect ell1, RotatedRect ell2)// ell1 or ell2 closer to ctr_pt
{
    if((fabs(ctr_pt.x-ell2.center.x)+fabs(ctr_pt.y-ell2.center.y))>(fabs(ctr_pt.x-ell1.center.x)+fabs(ctr_pt.y-ell1.center.y)))
    return true;// Ascend
    else return false;
};

inline double getp2l_dis(Point pt, Vec2f line)// distance from point pt to line
{
    double dist=0;
    double pt_rho = sqrt(pt.x*pt.x+pt.y*pt.y);
    double pt_theta = atan((double)pt.y/pt.x);// convert x,y to ρ,θ
    dist = fabs(pt_rho*cos(pt_theta-line[1])-line[0]);
    return dist;
}


bool getCPts::EclipseExtract(Mat input,vector<RotatedRect>& ell, vector<Vec2f>& lines)
{
    vector<vector<Point>> output;
    findContours(input, output, RETR_EXTERNAL, CHAIN_APPROX_NONE);
    // filter the contours
    // vector<RotatedRect> ell;
    vector<Point> centers;
    for(auto p=output.begin();p!=output.end();){
        if(p->size()>input.cols/8||p->size()<160)
        p=output.erase(p);
        else 
        {
            RotatedRect tmp=fitEllipse(*p);
            double ratio = tmp.size.height/tmp.size.width;
            if(ratio>3||tmp.angle<80||tmp.angle>100)p=output.erase(p);
            else
            {
                centers.push_back(tmp.center);
                ell.push_back(tmp);
                ++p;
            }
        }
    }

    Mat contour(input.size(),CV_8UC3);
    for(int i=0;i<output.size();++i){
        drawContours(contour, output, i, Scalar(255,255,255));
    }
    imwrite("contours.jpg",contour);

    Mat lineImg(input.size(),CV_8UC1, Scalar(0));
    for(int i=0;i<centers.size();++i){
        lineImg.at<uchar>(centers[i])=255;
    }
    // vector<Vec2f> lines;
    HoughLines(lineImg, lines, 10, CV_PI/300, 3);
    // remove similar lines
    contour.copyTo(lineImg);
    drawlines(lines, lineImg);
    imwrite("lines.jpg", lineImg);
    return true;
}

bool getCPts::Points2Group(vector<RotatedRect> ells, vector<Vec2f> lines, vector<vector<RotatedRect>>& ell_group)
{
    if(lines.size()<2)return false;
    double d=0;
    ell_group.resize(lines.size());
    // for(int i=0;i<lines.size();++i){
    //     ell_group[i].resize(4);
    // }
    int size=ells.size();
    for(auto ell=ells.begin();ell!=ells.end();)
    {
        for(int i=0;i<lines.size();++i)
        {
            
            d = getp2l_dis((*ell).center,lines[i]);
            if(d<10)
            {
                ell_group[i].push_back(*ell);
                ell=ells.erase(ell);
                i=-1;
            }
        }
        if(size==ells.size())ell=ells.erase(ell);
        else size=ells.size();
    }
    // remove the residual lines
    auto line=lines.begin();
    for(auto ell=ell_group.begin();ell!=ell_group.end();)
    {
        if((*ell).size()==0)
        {
            line=lines.erase(line);
            ell=ell_group.erase(ell);
        }
        else 
        {
            line++;
            ell++;
        }
    }

    ctr_pt = getCrossPoint(lines[0],lines[1]);
    // After grouping pts, sort each group according to dist from img center.
                // int dist[4]={0};
    for(int i=0;i<ell_group.size();++i)
    {
                // for(int j=0;j<4;++j){dist[j]=fabs(ctr_pt.x-ell_group[i][j].center.x)+fabs(ctr_pt.y-ell_group[i][j].center.y);}
        sort(ell_group[i].begin(),ell_group[i].end(),closer);
        
    }
    cout<<"done: "<<ell_group.size()<<" groups."<<endl;
    return true;
}

bool getCPts::EncodePts(vector<vector<RotatedRect>> ell_group, vector<vector<int>>& seq_list)
{
    seq_list.resize(ell_group.size());
    for(int i=0;i<ell_group.size();++i)
    {
        seq_list[i].resize(4);
        int flag=0;
        int a = ell_group[i][0].size.height<125?0:1;// smaller one:0
        int b = ell_group[i][1].size.height<125?0:1;// bigger one:1
        int c = ell_group[i][2].size.height<125?0:1;
        int d = ell_group[i][3].size.height<125?0:1;
        flag=((a<<3)+(b<<2)+(c<<1)+d)<<2;
        seq_list[i][0]=flag-1;
        seq_list[i][1]=flag-2;
        seq_list[i][2]=flag-3;
        seq_list[i][3]=flag-4;
    }
    return true;
}

bool getCPts::ShowCodes(Mat& img, vector<vector<RotatedRect>> ell_group, vector<vector<int>> seq_list)// 显示每个
{
    if(img.empty())return false;
    for(int i=0;i<ell_group.size();++i)
    {
        for(int j=0;j<ell_group[i].size();++j)
        {
            putText(img, to_string(seq_list[i][j]), ell_group[i][j].center, 1, 12, Scalar(255,255,255),2);
        }
    }
    return true;
}
