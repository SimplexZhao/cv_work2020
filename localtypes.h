#ifndef _ZXB_TYPES
#define _ZXB_TYPES

struct mPoint2d
{
    double x;
    double y;
    int id;
    int flag;
};

struct mPoint3d
{
    double x;
    double y;
    double z;
    int id;
    int flag; // control(1) or check(2)
};

struct IOPs
{
    double f;
    double x0;
    double y0;
    int width;
    int height;
    double dx;
    double dy;
};

struct EOPs
{
    double phi;
    double omega;
    double kappa;
    double offX;
    double offY;
    double offZ;
    // double DLT_s;
};

struct ImgPt
{
    float x;
    float y;
};

struct Disp2D
{
    float p;
    float q;
};

#endif