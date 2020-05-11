# 计算机视觉实习——源码使用说明

## 使用方法

### 1. 编译源代码

> In basic_mat.hpp, `MatrixMulti and MatrixInversion` are to be completed by the user. Recommended to replace the methods in basic_mat.hpp with those in Eigen.

安装OpenCV，vs201x；

已经安装vs某版本，配置一个OpenCV环境生成可执行文件；或者打开命令行`x64 Native Tools Command Prompt for VS 2017`进入源代码文件夹，运行：

```shell
cl.exe /Zi /EHsc /Fe: Dirname1/main.exe ./main.cpp /I D:/Libs/opencv4.1.0/include /I D:/Libs/opencv4.1.0/include/opencv2 /link D:/Libs/opencv4.1.0/x64/vc15/lib/opencv_world410.lib
```

即cl.exe /Zi /EHsc /Fe: `可执行文件名` `主程序代码` /I `OpenCV包含路径1` /I `OpenCV包含路径2` /link `OpenCV的lib库`  
编译成功即可

### 运行程序

在main.exe所在文件夹打开命令行，运行如下命令

```shell
.\main.exe "./编程实习材料/无畸变影像-左.bmp" "./编程实习材料/无畸变影像-右.bmp" "./编程实习材料/内部参数.txt" "./编程实习材料/控制点坐标.txt"
```

四个参数分别指明：`左影像路径`，`右影像路径`，`内部参数文件路径`，`控制点坐标文件路径`

### 运行结果

在源码相同文件夹下，生成了一系列文件：

文件名|说明
--|--
contours.jpg|控制点轮廓检测结果
lines.jpg|提取共线控制点结果
EOPAdj_left.txt/EOPAdj_right.txt|左右影像后交平差结果
EOPs.txt|左右影像外参数
precision_report.txt|前交精度报告
intersection——result.txt|前交控制点坐标结果
imgLr.jpg/imgRr.jpg|左右影像核线纠正结果
disp.jpg|密集匹配视差图
重建结果.obj|三维重建结果

打开obj文件推荐使用[Meshlab](http://www.meshlab.net/)
