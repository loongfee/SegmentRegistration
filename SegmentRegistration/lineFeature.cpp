#include "lineFeature.h"
void LineFeature(IplImage* img,vector<cvline_polar>& lines_polar_list,vector<vector<double>> &LinesMSLD,vector<int>&LineNum)
{
	int i;//检测出的直线的下标
	int j;//检测出的直线上的点的下标
	int r;
	int m;//检测出的直线的平行邻域内直线的下标
	const int R=20;//平行邻域的大小
	int  M=2*R+1;
	int N;
	vector <vector<double>> Feature;//所有直线的灰度矩阵
	vector<vector<double>> LinesFeature;//所有直线的梯度描述矩阵
	vector<int>LineLength;
	for(i=0;i!=int(lines_polar_list.size());i++)
{       
	if( (lines_polar_list[i].pt1.x>=R)&&(lines_polar_list[i].pt1.x<=(img->width-R))&&
		(lines_polar_list[i].pt1.y>=R)&& (lines_polar_list[i].pt1.y<=(img->height-R))&&
		(lines_polar_list[i].pt2.x>=R)&&(lines_polar_list[i].pt2.x<=(img->width-R))&&
	    (lines_polar_list[i].pt2.y>=R)&& (lines_polar_list[i].pt2.y<=(img->height-R)))
 {
	    N=abs(floor(lines_polar_list[i].pt1.x)-floor(lines_polar_list[i].pt2.x) );//检测出的直线上点的个数
		if(N>=1)
 {
		LineNum.push_back (i);
	    LineLength.push_back(N);
		vector<fPoint> PointInLine;//检测出的直线上的点
		vector<fPoint> PointsInRegion;//存放平行邻域内各点坐标的数组
		vector<double>FeatureValue;//平行域内各点灰度值
		double x0,y0,xj,yj,dx,dy;
	    if (lines_polar_list[i].theta<PI/2&&lines_polar_list[i].theta>=0)//theta在0~π/2之间时
    {  
		x0=(lines_polar_list[i].pt1.x<lines_polar_list[i].pt2.x)?lines_polar_list[i].pt1.x:lines_polar_list[i].pt2.x;//从具有较小x、y坐标的点开始
		y0=(lines_polar_list[i].pt1.y<lines_polar_list[i].pt2.y)?lines_polar_list[i].pt1.y:lines_polar_list[i].pt2.y;
		PointInLine.push_back(fPoint(x0,y0));
		for(j=1;j!=N;j++)
		{
			
			xj=PointInLine[0].x+j;
			yj=PointInLine[0].y+j*tan(lines_polar_list[i].theta);
			PointInLine.push_back(fPoint(xj,yj));
		}	
		 for(r=R;r>=-R;r--) //获得平行邻域内各点坐标值
		{
		   for(j=0;j!=N;j++)
		   {
			   dx=PointInLine[j].x-r*sin(lines_polar_list[i].theta);
			   dy=PointInLine[j].y+r*cos(lines_polar_list[i].theta);
			   PointsInRegion.push_back(fPoint(dx,dy));
		   }
		}
	}

	else if (lines_polar_list[i].theta>-PI/2&&lines_polar_list[i].theta<0)//theta在-π/2~0之间时
	{
		x0=(lines_polar_list[i].pt1.x>lines_polar_list[i].pt2.x)?lines_polar_list[i].pt1.x:lines_polar_list[i].pt2.x;
		y0=(lines_polar_list[i].pt1.y<lines_polar_list[i].pt2.y)?lines_polar_list[i].pt1.y:lines_polar_list[i].pt2.y;
		PointInLine.push_back(fPoint(x0,y0));
		for(j=1;j!=N;j++)
		{
			xj=PointInLine[0].x-j;
			yj=PointInLine[0].y-j*tan(lines_polar_list[i].theta);
			PointInLine.push_back(fPoint(xj,yj));
		}
		 for(r=R;r>=-R;r--) //获得平行邻域内各点坐标值
		{
		   for(j=0;j!=N;j++)
		   {
			   dx=PointInLine[j].x+r*sin(lines_polar_list[i].theta);
			   dy=PointInLine[j].y-r*cos(lines_polar_list[i].theta);
			   PointsInRegion.push_back(fPoint(dx,dy));
		   }
		}
	}

	else                                                                                             //theta等于π/2时
	{
        x0=lines_polar_list[i].pt1.x;
	    y0=(lines_polar_list[i].pt1.y<lines_polar_list[i].pt2.y)?lines_polar_list[i].pt1.y:lines_polar_list[i].pt2.y;//从直线中具有较小y坐标的点开始
		PointInLine.push_back(fPoint(x0,y0));
		for(j=1;j!=N;j++)
		{
           xj=PointInLine[0].x;
		   yj=PointInLine[0].y+j;
		   PointInLine.push_back(fPoint(xj,yj));
		}	
		 for(r=R;r>=-R;r--) //获得平行邻域内各点坐标值
		{
		   for(j=0;j!=N;j++)
		   {
			   dx=PointInLine[j].x-r;
			   dy=PointInLine[j].y;
			   PointsInRegion.push_back(fPoint(dx,dy));
		   }
		}
	}

		  for(m=0;m!=M;m++)//获得每个点的灰度值
		{
		  for(j=0;j!=N;j++)
		   {
		      CvScalar s= cvGet2D(img,PointsInRegion[j +m*N].y,PointsInRegion[j +m*N].x);
		      double pix= s.val[0];
		      FeatureValue.push_back(pix);
		   }
		}
		Feature.push_back(FeatureValue);
		PointInLine.swap(vector<fPoint>());
		PointsInRegion.swap(vector<fPoint>());
		FeatureValue.swap(vector<double>());
	}
  }
	else
	{;}
}
	//构造梯度均值-标准差描述子，对灰度矩阵中的(M-1)*(N-1)求梯度
	vector <vector<double>> All_Grad_dL;//dL方向的梯度
	for(i=0;i!=Feature.size();i++)
	{
		vector<double>Each_Grad_dL;
		N=LineLength[i];
		for(m=0;m!=M-1;m++)
		{
			for(j=0;j!=N-1;j++)
			{
				Each_Grad_dL.push_back (Feature[i][j+1+m*N]-Feature[i][j+m*N]);
			}
		}
		All_Grad_dL.push_back (Each_Grad_dL);
		Each_Grad_dL.swap(vector<double>());
	}

	vector <vector<double>> All_Grad_dT;//dT方向的梯度
	for(i=0;i!=Feature.size();i++)
	{
		vector<double>Each_Grad_dT;
		N=LineLength[i];
		for(j=0;j!=N-1;j++)
		{
			for(m=0;m!=M-1;m++)
			{
				Each_Grad_dT.push_back (Feature[i][j+(m+1)*N]-Feature[i][j+m*N]);
			}
		}
		All_Grad_dT.push_back (Each_Grad_dT);
		Each_Grad_dT.swap(vector<double>());
	}
	//两个方向的梯度合并
	for(i=0;i!=Feature.size();i++)
	{   
		vector<double> Each_Grad_dLdT;
		N=LineLength[i];
		for(m=0;m!=M-1;m++)
		{
			for(j=0;j!=N-1;j++)
				{
					Each_Grad_dLdT.push_back (All_Grad_dL[i][j+m*(N-1)]);	
			    } 
	     }
		for(m=0;m!=M-1;m++)
		{
			for(j=0;j!=N-1;j++)
				{
					Each_Grad_dLdT.push_back (All_Grad_dT[i][m+j*(M-1)]);//沿dT方向计算的梯度按dL方向排列，与dL方向的梯度对应。
			    } 
	     }
		 LinesFeature.push_back (Each_Grad_dLdT);//得到梯度描述矩阵
		 Each_Grad_dLdT.swap(vector<double>());
 	}
	    All_Grad_dL.swap(vector<vector<double>>());
	    All_Grad_dT.swap(vector<vector<double>>());

    //计算梯度描述矩阵的均值和标准差
    vector<vector<double>> AllMeanGrad;
	vector<vector<double>> AllStdGrad;
	for(i=0;i!=Feature.size();i++)
	{
		vector<double> EachMeanGrad;
		vector<double> EachStd;
		N=LineLength[i];
		for(m=0;m!=2*(M-1);m++)
		{
			double SumRowGray=0.0;
			for(j=0;j!=N-1;j++)
			{
				SumRowGray+=LinesFeature[i][j+m*(N-1)];
			}
			EachMeanGrad.push_back (SumRowGray/(N-1));//得到每行的梯度均值
			double SumSqr=0.0;
			for(j=0;j!=N-1;j++)
			{
				SumSqr+=(LinesFeature[i][j+m*(N-1)]-SumRowGray/(N-1))*(LinesFeature[i][j+m*(N-1)]-SumRowGray/(N-1));
			}
		    double Std=sqrt(SumSqr);
			EachStd.push_back(Std);//得到每行的标准差
		}
		AllMeanGrad.push_back (EachMeanGrad);//所有直线的梯度均值矩阵
		AllStdGrad.push_back (EachStd);//所有直线的梯度标准差矩阵
		EachMeanGrad.swap(vector<double>());
        EachStd.swap(vector<double>());
	}
	for(i=0;i!=Feature.size();i++)
	{
		double SumMeanGray=0.0;
		double SumStd=0.0;
		for(m=0;m!=2*(M-1);m++)
		{
		    SumMeanGray+=AllMeanGrad[i][m]*AllMeanGrad[i][m];
			SumStd+=AllStdGrad[i][m]*AllStdGrad[i][m];
		}
		double LengthMean=sqrt(SumMeanGray);
		double LengthStd=sqrt(SumStd);
		for(m=0;m!=2*(M-1);m++)
		{
			AllMeanGrad[i][m]=AllMeanGrad[i][m]/ LengthMean;//归一化的梯度均值
			AllStdGrad[i][m]=AllStdGrad[i][m]/LengthStd;//归一化的标准差
		}
	}

	for(i=0;i!=Feature.size();i++)
	{
		for(m=0;m!=2*(M-1);m++)
		{
			AllMeanGrad[i].push_back (AllStdGrad[i][m]);//将标准差与梯度均值合并
		}
		LinesMSLD.push_back (AllMeanGrad[i]);//得到MSLD描述矩阵
	}
	AllMeanGrad.swap(vector<vector<double>>());
	AllStdGrad.swap(vector<vector<double>>());
	LinesFeature.swap(vector<vector<double>>());
}

