#include<iostream>
#include<vector>
#include<string>
#include<fstream>
using namespace std;
#include"basic_mat.hpp"

class solver
{
    public:
    solver(){};
    virtual ~solver(){};
    static bool loadControlfile(string filename, vector<mPoint3d>& CPts); // 读取控制点文件
    static bool loadIn_Params(string filename, IOPs& params);
    static bool getEO(IOPs params, EOPs& params2, vector<mPoint2d> IPt_list, vector<mPoint3d> CPt_list);
    static bool nonlinear_EO(IOPs params, EOPs& params2, vector<mPoint2d> IPt_list, vector<mPoint3d> CPt_list, string flag="");
    static void writeEOPfile(EOPs params_Left, EOPs params_Right);
    static void loadEOPfile(EOPs& params_Left, EOPs& params_Right);
    static bool DualImg_Intersection(IOPs params, EOPs& params2_L, EOPs& params2_R, vector<mPoint2d> IPt_list_L, vector<mPoint2d> IPt_list_R, vector<mPoint3d> CPt_list);
};

bool solver::loadControlfile(string filename, vector<mPoint3d>& CPts)
{
    int size;
    ifstream fp(filename);
    fp>>size;
    CPts.resize(size);
    string line_buf;
    getline(fp, line_buf);
    int i=0;
    while(!fp.eof())
    {
        fp>>i;
        fp>>CPts[i].x>>CPts[i].y>>CPts[i].z>>CPts[i].flag;
    }
    fp.close();
    return true;
}

bool solver::loadIn_Params(string filename, IOPs& params)
{
    ifstream fp(filename);
    string line_buf;
    getline(fp, line_buf);
    params.dx=params.dy=1;
    params.height=-1;
    fp>>params.x0>>params.y0>>params.f>>params.width>>params.height;
    fp.close();
    if(params.height==-1)return false;
    return true;
}

bool solver::getEO(IOPs params, EOPs& params2, vector<mPoint2d> IPt_list, vector<mPoint3d> CPt_list)
{
    int size = IPt_list.size();// Image pts' number.
    double* K=new double[size*2*8];
    double* tmpK=K;
    double h[8];// h matrix in 2D-DLT
    for(int i=0;i<8;i++)h[i]=0;// Initialize h
    double* U=new double[size*2];
    double* tmpU=U;
    for(int i=0;i<IPt_list.size();++i)
    {
        // Every image point can form two equations.
        // from u,v to x,y
        IPt_list[i].x-=params.x0;
        *tmpU = IPt_list[i].x;
        ++tmpU;
        // row fill start
        *tmpK=CPt_list[IPt_list[i].id].x; ++tmpK;
        *tmpK=CPt_list[IPt_list[i].id].y; ++tmpK;
        *tmpK=1; ++tmpK;
        *tmpK=0; ++tmpK;
        *tmpK=0; ++tmpK;
        *tmpK=0; ++tmpK;
        // *tmpK=0; ++tmpK;
        *tmpK=-CPt_list[IPt_list[i].id].x*IPt_list[i].x; ++tmpK;
        *tmpK=-CPt_list[IPt_list[i].id].y*IPt_list[i].x; ++tmpK;
        // row fill end
        IPt_list[i].y=params.y0-IPt_list[i].y;
        *tmpU=IPt_list[i].y;
        ++tmpU;
        // row fill start
        *tmpK=0; ++tmpK;
        *tmpK=0; ++tmpK;
        *tmpK=0; ++tmpK;
        //*tmpK=0; ++tmpK;
        *tmpK=CPt_list[IPt_list[i].id].x; ++tmpK;
        *tmpK=CPt_list[IPt_list[i].id].y; ++tmpK;
        *tmpK=1; ++tmpK;
        *tmpK=-CPt_list[IPt_list[i].id].x*IPt_list[i].y; ++tmpK;
        *tmpK=-CPt_list[IPt_list[i].id].y*IPt_list[i].y; ++tmpK;
        // row fill end
    }
    double* KT=new double[size*2*8];
    double KTK[64];
    double KTK2[64];
    // double v[64];// test the inversion accuracy.
    double KTU[8];
    MatOperation::Rotate(K, KT, size*2, 8);
    MatOperation::MatrixMulti(KT,K,KTK, 8, 8, size*2); // memcpy(KTK2,KTK,sizeof(double)*64); //record KTK to verify inversion
    MatOperation::MatrixMulti(KT,U,KTU,8,1,size*2);
    MatOperation::MatrixInversion(KTK, 8); // KTK^(-1)      // double* KTKKT = new double[2*9*size]; MatOperation::MatrixMulti(KTK,KTK2,v,8,8,8);// verify the inversion // MatOperation::MatrixMulti(KTK,KT,KTKKT,9,size*2,9);
    MatOperation::MatrixMulti(KTK,KTU,h,8,1,8);           // MatOperation::MatrixMulti(KTKKT,U,M,9,1,size*2); delete[] KTKKT;
    
    for(int i=0;i<6;++i)h[i]/=-params.f;
    double t3 = -(sqrt(1/(h[0]*h[0]+h[3]*h[3]+h[6]*h[6]))+sqrt(1/(h[1]*h[1]+h[4]*h[4]+h[7]*h[7])))/2.0;
    for(int i=0;i<8;++i)h[i]*=t3;
    // MatOperation::ShowMat(h,8);
    params2.omega=asin(-h[7]);
    params2.phi=asin(-h[6]/cos(params2.omega));
    params2.kappa = atan2(h[1]/cos(params2.omega),h[4]/cos(params2.omega));// better than above line.
    double t[3]={h[2],h[5],t3};
    double rt[3]={0.0};
    double R1[9]={0.0};
    double R2[9];
    MatOperation::GetRotationMat(R2, params2.phi, params2.omega, params2.kappa, 2);
    // memcpy(R1, h, sizeof(double)*8);
    // MatOperation::ShowMat(R1,3,3);
    // MatOperation::Rotate(R1,R2,3,3);
    // MatOperation::Xprod_3d(R2,R2+3,R2+6);
    MatOperation::Rotate(R2,R1,3,3);
    // double I[9]; MatOperation::MatrixMulti(R2,R1,I,3,3,3); // verification of R
    MatOperation::MatrixMulti(R1,t,rt,3,1,3);
    // Compute offset in Photogrammetry.
    params2.offX=-rt[0];
    params2.offY=-rt[1];
    params2.offZ=-rt[2];
    delete[] K;
    delete[] KT;
    delete[] U;
    return true;
}

bool solver::nonlinear_EO(IOPs params, EOPs& params2, vector<mPoint2d> IPt_list, vector<mPoint3d> CPt_list,string flag)
{
    struct Pt_Pair{
        double Ix;
        double Iy;
        double Gx;
        double Gy;
        double Gz;
    };
    int size = IPt_list.size();
	vector<Pt_Pair> group;
    group.resize(size);
    int tmpi=0;
    for(int i=0;i<IPt_list.size();++i)
    {
        group[i].Ix=IPt_list[i].x-params.x0;
        group[i].Iy=-IPt_list[i].y+params.y0;
        group[i].Gx=CPt_list[IPt_list[i].id].x;
        group[i].Gy=CPt_list[IPt_list[i].id].y;
        group[i].Gz=CPt_list[IPt_list[i].id].z;
    }
    // Init val
	// for (auto& val : group)params2.offX += val.Gx / 4.0;
	// for (auto& val : group)params2.offY += val.Gy / 4.0;
	// //for (auto& val : group)res.Z += val.Gz / 4.0;
	// for (auto& val : group) {
	// 	val.Ix /= 1000.0; val.Iy /= 1000.0;
	// }
	// params2.offZ += 15000 * params.f;
	double rotate[9];
	double Z_ba, f_XbyZ, f_YbyZ;
	Z_ba=f_XbyZ = f_YbyZ = 0;

	double* v=new double[size*2];
	double* B=new double[size*2*6];
	double* l = new double[size*2];
	double* p2B=B;
	double* Bt=new double[size*2*6];
	double Nbb[36];
	double Btl[6];
	double x[6] = { 0 };
	x[5] = 1;
	short count = 0;
	while (fabs(x[5])>5e-6||fabs(x[4])>5e-6||fabs(x[3])>5e-6||count<6) {
		// Calculate R
		rotate[0] = cos(params2.phi)*cos(params2.kappa) - sin(params2.phi)*sin(params2.omega)*sin(params2.kappa);
		rotate[1] = -cos(params2.phi)*sin(params2.kappa) - sin(params2.phi)*sin(params2.omega)*cos(params2.kappa);
		rotate[2] = -sin(params2.phi)*cos(params2.omega);
		rotate[3] = cos(params2.omega)*sin(params2.kappa);
		rotate[4] = cos(params2.omega)*cos(params2.kappa);
		rotate[5] = -sin(params2.omega);
		rotate[6] = sin(params2.phi)*cos(params2.kappa) + cos(params2.phi)*sin(params2.omega)*sin(params2.kappa);
		rotate[7] = -sin(params2.phi)*sin(params2.kappa) + cos(params2.phi)*sin(params2.omega)*cos(params2.kappa);
		rotate[8] = cos(params2.phi)*cos(params2.omega);
		// Compute L, B
		for (int i = 0; i < size; i++) {
			Z_ba = rotate[2] * (group[i].Gx - params2.offX) + rotate[5] * (group[i].Gy - params2.offY) + rotate[8] * (group[i].Gz - params2.offZ);
			f_XbyZ = params.f * (rotate[0] * (group[i].Gx - params2.offX) + rotate[3] * (group[i].Gy - params2.offY) + rotate[6] * (group[i].Gz - params2.offZ))
				/ Z_ba;
			f_YbyZ = params.f * (rotate[1] * (group[i].Gx - params2.offX) + rotate[4] * (group[i].Gy - params2.offY) + rotate[7] * (group[i].Gz - params2.offZ))
				/ Z_ba;

			l[i * 2] = group[i].Ix + f_XbyZ;//double lx
			l[i * 2 + 1] = group[i].Iy + f_YbyZ;//double ly

			*p2B = (rotate[0] * params.f + rotate[2] * group[i].Ix) / Z_ba; p2B++;//double xX
			*p2B = (rotate[3] * params.f + rotate[5] * group[i].Ix) / Z_ba; p2B++;//double xY
			*p2B = (rotate[6] * params.f + rotate[8] * group[i].Ix) / Z_ba; p2B++;//double xZ
			*p2B = group[i].Iy*sin(params2.omega) - ((group[i].Ix / params.f)*(group[i].Ix*cos(params2.kappa) - group[i].Iy*sin(params2.kappa)) + params.f*cos(params2.kappa))*cos(params2.omega);
			p2B++;		*p2B = -params.f*sin(params2.kappa) - (group[i].Ix / params.f)*(group[i].Ix*sin(params2.kappa) + group[i].Iy*cos(params2.kappa));
			p2B++;		*p2B = group[i].Iy;// partial phi, omega, kappa
			p2B++;
			*p2B = (rotate[1] * params.f + rotate[2] * group[i].Iy) / Z_ba; p2B++;//double yX
			*p2B = (rotate[4] * params.f + rotate[5] * group[i].Iy) / Z_ba; p2B++;//double yY
			*p2B = (rotate[7] * params.f + rotate[8] * group[i].Iy) / Z_ba; p2B++;//double yZ
			*p2B = -group[i].Ix*sin(params2.omega) - ((group[i].Iy / params.f)*(group[i].Ix*cos(params2.kappa) - group[i].Iy*sin(params2.kappa)) - params.f*sin(params2.kappa))*cos(params2.omega);
			p2B++;	*p2B = -params.f*cos(params2.kappa) - (group[i].Iy / params.f)*(group[i].Ix*sin(params2.kappa) + group[i].Iy*cos(params2.kappa));
			p2B++;	*p2B = -group[i].Ix;// partial phi, omega, kappa
			p2B++;
		}
		p2B = B;
		MatOperation::Rotate(B, Bt, size*2 , 6);
		MatOperation::MatrixMulti(Bt, B, Nbb, 6, 6, size*2);
		MatOperation::MatrixMulti(Bt, l, Btl, 6, 1, size*2);
		MatOperation::MatrixInversion(Nbb, 6);
		MatOperation::MatrixMulti(Nbb, Btl, x, 6, 1, 6);
		MatOperation::MatrixMulti(B, x, v, size*2, 1, 6);
        MatOperation::Subtraction(v,v,l,size*2);
		params2.offX += x[0];
		params2.offY += x[1];
		params2.offZ += x[2];
		params2.phi += x[3];
		params2.omega += x[4];
		params2.kappa += x[5];
		count++;
	}
	ofstream pf;
	pf.open("EOPAdj_"+flag+".txt",std::ofstream::trunc);
	if (!pf.is_open())return false;
	pf << "迭代次数：" << count<<endl;
	for (int i = 0; i < 9; i++) {
		pf << rotate[i];
		if ((i+1)% 3 == 0)pf << endl;
		else pf << " ";
	}

	int p=0;
	for (auto ele : group) {
		pf << "(x,y)=" << "("<<ele.Ix+v[2*p]<<"," <<ele.Iy+v[2*p+1]<<")"<<endl;
		++p;
	}
	double sigma[1];
	MatOperation::MatrixMulti(v, v, sigma, 1, 1, 8);
	*sigma = sqrt(*sigma / 2);
	pf.setf(ios::fixed);
	pf << "外方位元素：" << endl << "Xs=" << params2.offX << " 精度："<<*sigma*sqrt(Nbb[0])<< endl;
	pf << "Ys=" << params2.offY << " 精度：" << *sigma*sqrt(Nbb[7]) << endl;
	pf << "Zs=" << params2.offZ << " 精度：" << *sigma*sqrt(Nbb[14]) << endl;
	pf << "φ=" << -atan(rotate[2] / rotate[8])<< " 精度：" << *sigma*sqrt(Nbb[21]) << endl;
	pf << "ω=" << -asin(rotate[5])<< " 精度：" << *sigma*sqrt(Nbb[28]) << endl;
	pf << "κ=" << atan(rotate[3]/rotate[4])+3.14159265358979<< " 精度：" << *sigma*sqrt(Nbb[35]) << endl;
	pf << "单位权中误差：" <<*sigma<< endl;
	pf.close();

    delete[] v;
	delete[] l;
	delete[] B;
	delete[] Bt;
    return true;
}

void solver::writeEOPfile(EOPs params_Left, EOPs params_Right)
{
    ofstream ofs("EOPs.txt");
    ofs<<"φ "<<"ω "<<"κ "<<"Xs "<<"Ys "<<"Zs"<<endl;
    ofs<<params_Left.phi<<" "<<params_Left.omega<<" "<<params_Left.kappa<<" "<<params_Left.offX<<" "<<params_Left.offY<<" "<<params_Left.offZ<<endl;
    ofs<<params_Right.phi<<" "<<params_Right.omega<<" "<<params_Right.kappa<<" "<<params_Right.offX<<" "<<params_Right.offY<<" "<<params_Right.offZ<<endl;
}

void solver::loadEOPfile(EOPs& params_Left, EOPs& params_Right)
{
    ifstream ifs("EOPs.txt");
    string tmp;
    getline(ifs, tmp);
    ifs>>params_Left.phi>>params_Left.omega>>params_Left.kappa>>params_Left.offX>>params_Left.offY>>params_Left.offZ;
    ifs>>params_Right.phi>>params_Right.omega>>params_Right.kappa>>params_Right.offX>>params_Right.offY>>params_Right.offZ;
}

bool id_bigger(mPoint2d pt1, mPoint2d pt2)// ell1 or ell2 closer to ctr_pt
{
    if(pt1.id<pt2.id)return true;// id Ascend
    else return false;
};

bool solver::DualImg_Intersection(IOPs params, EOPs& params2_L, EOPs& params2_R, vector<mPoint2d> IPt_list_L, vector<mPoint2d> IPt_list_R, vector<mPoint3d> CPt_list)
{
    for(int i=0;i<IPt_list_L.size();++i)
    {
            IPt_list_L[i].x-=params.x0;
            IPt_list_L[i].y=params.y0-IPt_list_L[i].y;
            IPt_list_R[i].x-=params.x0;
            IPt_list_R[i].y=params.y0-IPt_list_R[i].y;
    }
    // Compute Rotation matrix for left and right image, k=0 and k=1 respectively.
    double a1[2], a2[2], a3[2],b1[2], b2[2], b3[2],c1[2], c2[2], c3[2];
    // left
    a1[0] = cos(params2_L.phi)*cos(params2_L.kappa) - sin(params2_L.phi)*sin(params2_L.omega)*sin(params2_L.kappa);
    a2[0] = -cos(params2_L.phi)*sin(params2_L.kappa) - sin(params2_L.phi)*sin(params2_L.omega)*cos(params2_L.kappa);
    a3[0] = -sin(params2_L.phi)*cos(params2_L.omega);
    b1[0] = cos(params2_L.omega)*sin(params2_L.kappa);
    b2[0] = cos(params2_L.omega)*cos(params2_L.kappa);
    b3[0] = -sin(params2_L.omega);
    c1[0] = sin(params2_L.phi)*cos(params2_L.kappa) + cos(params2_L.phi)*sin(params2_L.omega)*sin(params2_L.kappa);
    c2[0] = -sin(params2_L.phi)*sin(params2_L.kappa) + cos(params2_L.phi)*sin(params2_L.omega)*cos(params2_L.kappa);
    c3[0] = cos(params2_L.phi)*cos(params2_L.omega);
    // right 
    a1[1] = cos(params2_R.phi)*cos(params2_R.kappa) - sin(params2_R.phi)*sin(params2_R.omega)*sin(params2_R.kappa);
    a2[1] = -cos(params2_R.phi)*sin(params2_R.kappa) - sin(params2_R.phi)*sin(params2_R.omega)*cos(params2_R.kappa);
    a3[1] = -sin(params2_R.phi)*cos(params2_R.omega);
    b1[1] = cos(params2_R.omega)*sin(params2_R.kappa);
    b2[1] = cos(params2_R.omega)*cos(params2_R.kappa);
    b3[1] = -sin(params2_R.omega);
    c1[1] = sin(params2_R.phi)*cos(params2_R.kappa) + cos(params2_R.phi)*sin(params2_R.omega)*sin(params2_R.kappa);
    c2[1] = -sin(params2_R.phi)*sin(params2_R.kappa) + cos(params2_R.phi)*sin(params2_R.omega)*cos(params2_R.kappa);
    c3[1] = cos(params2_R.phi)*cos(params2_R.omega);

    sort(IPt_list_L.begin(), IPt_list_L.end(), id_bigger);
    sort(IPt_list_R.begin(), IPt_list_R.end(), id_bigger);

    auto ILL=IPt_list_L.begin();
    auto ILR=IPt_list_R.begin();
    while(ILL!=IPt_list_L.end()&&ILR!=IPt_list_R.end())
    {
        if(ILL->id==ILR->id)
        {++ILL;
        ++ILR;}
        else if(ILL->id>ILR->id)
         ILR=IPt_list_R.erase(ILR);
        else ILL=IPt_list_L.erase(ILL);
    }
    int minisize=0;
    if(IPt_list_L.size()<IPt_list_R.size())
        minisize=IPt_list_L.size();
    else 
        minisize=IPt_list_R.size();
    mPoint3d* Pt3D = new mPoint3d[minisize];
    
    double dP[3];
    double C[9];//B'*B
    double W[3];//B'*l
    double *dV = new double[2*2];
    double*B = new double[3*2*2];
    double *p2B = B;
    double *l = new double[2*2];
    double *p2l = l;
    for(int i=0;i<minisize;++i)
    {
        Pt3D[i].x=Pt3D[i].y=Pt3D[i].z=0;
        Pt3D[i].id=IPt_list_L[i].id;
        for (int iter = 0; iter < 10; iter++)// maximum iteration number: 10
        {
            p2B = B; p2l = l;
            //构造B,l
            {// Left image
                double Z_ba = a3[0] * (Pt3D[i].x - params2_L.offX) + b3[0] * (Pt3D[i].y -params2_L.offY) + c3[0] * (Pt3D[i].z -params2_L.offZ);
                double f_XbyZ = params.f*(a1[0] * (Pt3D[i].x -params2_L.offX) + b1[0] * (Pt3D[i].y -params2_L.offY) + c1[0] * (Pt3D[i].z -params2_L.offZ))
                    / Z_ba;
                double f_YbyZ = params.f*(a2[0] * (Pt3D[i].x -params2_L.offX) + b2[0] * (Pt3D[i].y -params2_L.offY) + c2[0] * (Pt3D[i].z -params2_L.offZ))
                    / Z_ba;

                *p2l = IPt_list_L[i].x + f_XbyZ; p2l++;//double lx
                *p2l = IPt_list_L[i].y + f_YbyZ; p2l++;//double ly
                *p2B = -(a1[0] * params.f + a3[0] * IPt_list_L[i].x) / Z_ba; p2B++;//double xX
                *p2B = -(b1[0] * params.f + b3[0] * IPt_list_L[i].x) / Z_ba; p2B++;//double xY
                *p2B = -(c1[0] * params.f + c3[0] * IPt_list_L[i].x) / Z_ba; p2B++;//double xZ
                *p2B = -(a2[0] * params.f + a3[0] * IPt_list_L[i].y) / Z_ba; p2B++;//double yX
                *p2B = -(b2[0] * params.f + b3[0] * IPt_list_L[i].y) / Z_ba; p2B++;//double yY
                *p2B = -(c2[0] * params.f + c3[0] * IPt_list_L[i].y) / Z_ba; p2B++;//double yZ
            }
            {// Right image
                double Z_ba = a3[1] * (Pt3D[i].x - params2_R.offX) + b3[1] * (Pt3D[i].y -params2_R.offY) + c3[1] * (Pt3D[i].z -params2_R.offZ);
                double f_XbyZ = params.f*(a1[1] * (Pt3D[i].x -params2_R.offX) + b1[1] * (Pt3D[i].y -params2_R.offY) + c1[1] * (Pt3D[i].z -params2_R.offZ))
                    / Z_ba;
                double f_YbyZ = params.f*(a2[1] * (Pt3D[i].x -params2_R.offX) + b2[1] * (Pt3D[i].y -params2_R.offY) + c2[1] * (Pt3D[i].z -params2_R.offZ))
                    / Z_ba;

                *p2l = IPt_list_R[i].x + f_XbyZ; p2l++;//double lx
                *p2l = IPt_list_R[i].y + f_YbyZ; p2l++;//double ly
                *p2B = -(a1[1] * params.f + a3[1] * IPt_list_R[i].x) / Z_ba; p2B++;//double xX
                *p2B = -(b1[1] * params.f + b3[1] * IPt_list_R[i].x) / Z_ba; p2B++;//double xY
                *p2B = -(c1[1] * params.f + c3[1] * IPt_list_R[i].x) / Z_ba; p2B++;//double xZ
                *p2B = -(a2[1] * params.f + a3[1] * IPt_list_R[i].y) / Z_ba; p2B++;//double yX
                *p2B = -(b2[1] * params.f + b3[1] * IPt_list_R[i].y) / Z_ba; p2B++;//double yY
                *p2B = -(c2[1] * params.f + c3[1] * IPt_list_R[i].y) / Z_ba; p2B++;//double yZ
            }            
            if ((p2B - B) != 6 * 2 || (p2l - l) != 2 * 2) {std::cout<<"Error!!"<<endl; }//赋值不完全则暂停
            double *Bt = new double[6*2];
            MatOperation::Rotate(B, Bt, 2*2, 3);
            MatOperation::MatrixMulti(Bt, B, C, 3, 3, 2*2);
            MatOperation::MatrixMulti(Bt, l, W, 3, 1, 2*2);
            delete[]Bt;
            MatOperation::MatrixInversion(C, 3);
            MatOperation::MatrixMulti(C, W, dP, 3, 1, 3);
            Pt3D[i].x += dP[0];
            Pt3D[i].y += dP[1];
            Pt3D[i].z += dP[2];
            MatOperation::MatrixMulti(B, dP, dV, 2*2, 1, 3);
            if((dP[0]+dP[1]+dP[2])<1e-10)break;
        }
        for (int count = 0; count < 2*2; count++)
        {
            dV[count] -= l[count];
        }
        // for (int num_Img = 0, i = 0; num_Img< 2; num_Img++)//改正Point2D数组
        // {
        //     IPt_list_L[num_Img].x += dV[i]; ++i;
        //     IPt_list_L[num_Img].y += dV[i]; ++i;
        // }
    }
    //std::cout << dP[0] << " " << dP[1] << " " << dP[2]<<endl;//to see if deltas converge
    ofstream ofs_result;
    ofs_result.open("intersection_result.txt");
    // precision calculation
    double rms[3], max[3], tmp[3];// rms[0]=RMSx,...rms[2]=RMSz
    max[0]=max[1]=max[2]=0;
    for(int i=0;i<minisize;++i)
    {
        ofs_result << Pt3D[i].id << " " << Pt3D[i].x << " " << Pt3D[i].y << " " << Pt3D[i].z << endl;
        tmp[0]=fabs(Pt3D[i].x-CPt_list[Pt3D[i].id].x);
        tmp[1]=fabs(Pt3D[i].y-CPt_list[Pt3D[i].id].y);
        tmp[2]=fabs(Pt3D[i].z-CPt_list[Pt3D[i].id].z);
        for(int j=0;j<3;++j)
        {
            if(tmp[j]>max[j])max[j]=tmp[j];
            rms[j]+=tmp[j]*tmp[j];
        }
    }
    ofs_result.close();// intersection result written
    for(int j=0;j<3;++j)
    {
        rms[j]/=minisize;
        rms[j]=sqrt(rms[j]);
    }
    ofstream ofsprec;
    ofsprec.open("precision_report.txt");
    ofsprec<<"RMS X: "<<rms[0]<<endl;
    ofsprec<<"RMS Y: "<<rms[1]<<endl;
    ofsprec<<"RMS Z: "<<rms[2]<<endl;
    ofsprec<<"MAX X: "<<max[0]<<endl;
    ofsprec<<"MAX Y: "<<max[1]<<endl;
    ofsprec<<"MAX Z: "<<max[2]<<endl;
    ofsprec.close();
    // calculation over.
    delete[]B; delete[]l; delete[] dV; 
    delete[] Pt3D;
    return true;
}
