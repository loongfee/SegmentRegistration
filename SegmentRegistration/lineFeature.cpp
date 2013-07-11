#include "lineFeature.h"
void LineFeature(IplImage* img,vector<cvline_polar>& lines_polar_list,vector<vector<double>> &LinesMSLD,vector<int>&LineNum)
{
	int i;//������ֱ�ߵ��±�
	int j;//������ֱ���ϵĵ���±�
	int r;
	int m;//������ֱ�ߵ�ƽ��������ֱ�ߵ��±�
	const int R=20;//ƽ������Ĵ�С
	int  M=2*R+1;
	int N;
	vector <vector<double>> Feature;//����ֱ�ߵĻҶȾ���
	vector<vector<double>> LinesFeature;//����ֱ�ߵ��ݶ���������
	vector<int>LineLength;
	for(i=0;i!=int(lines_polar_list.size());i++)
{       
	if( (lines_polar_list[i].pt1.x>=R)&&(lines_polar_list[i].pt1.x<=(img->width-R))&&
		(lines_polar_list[i].pt1.y>=R)&& (lines_polar_list[i].pt1.y<=(img->height-R))&&
		(lines_polar_list[i].pt2.x>=R)&&(lines_polar_list[i].pt2.x<=(img->width-R))&&
	    (lines_polar_list[i].pt2.y>=R)&& (lines_polar_list[i].pt2.y<=(img->height-R)))
 {
	    N=abs(floor(lines_polar_list[i].pt1.x)-floor(lines_polar_list[i].pt2.x) );//������ֱ���ϵ�ĸ���
		if(N>=1)
 {
		LineNum.push_back (i);
	    LineLength.push_back(N);
		vector<fPoint> PointInLine;//������ֱ���ϵĵ�
		vector<fPoint> PointsInRegion;//���ƽ�������ڸ������������
		vector<double>FeatureValue;//ƽ�����ڸ���Ҷ�ֵ
		double x0,y0,xj,yj,dx,dy;
	    if (lines_polar_list[i].theta<PI/2&&lines_polar_list[i].theta>=0)//theta��0~��/2֮��ʱ
    {  
		x0=(lines_polar_list[i].pt1.x<lines_polar_list[i].pt2.x)?lines_polar_list[i].pt1.x:lines_polar_list[i].pt2.x;//�Ӿ��н�Сx��y����ĵ㿪ʼ
		y0=(lines_polar_list[i].pt1.y<lines_polar_list[i].pt2.y)?lines_polar_list[i].pt1.y:lines_polar_list[i].pt2.y;
		PointInLine.push_back(fPoint(x0,y0));
		for(j=1;j!=N;j++)
		{
			
			xj=PointInLine[0].x+j;
			yj=PointInLine[0].y+j*tan(lines_polar_list[i].theta);
			PointInLine.push_back(fPoint(xj,yj));
		}	
		 for(r=R;r>=-R;r--) //���ƽ�������ڸ�������ֵ
		{
		   for(j=0;j!=N;j++)
		   {
			   dx=PointInLine[j].x-r*sin(lines_polar_list[i].theta);
			   dy=PointInLine[j].y+r*cos(lines_polar_list[i].theta);
			   PointsInRegion.push_back(fPoint(dx,dy));
		   }
		}
	}

	else if (lines_polar_list[i].theta>-PI/2&&lines_polar_list[i].theta<0)//theta��-��/2~0֮��ʱ
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
		 for(r=R;r>=-R;r--) //���ƽ�������ڸ�������ֵ
		{
		   for(j=0;j!=N;j++)
		   {
			   dx=PointInLine[j].x+r*sin(lines_polar_list[i].theta);
			   dy=PointInLine[j].y-r*cos(lines_polar_list[i].theta);
			   PointsInRegion.push_back(fPoint(dx,dy));
		   }
		}
	}

	else                                                                                             //theta���ڦ�/2ʱ
	{
        x0=lines_polar_list[i].pt1.x;
	    y0=(lines_polar_list[i].pt1.y<lines_polar_list[i].pt2.y)?lines_polar_list[i].pt1.y:lines_polar_list[i].pt2.y;//��ֱ���о��н�Сy����ĵ㿪ʼ
		PointInLine.push_back(fPoint(x0,y0));
		for(j=1;j!=N;j++)
		{
           xj=PointInLine[0].x;
		   yj=PointInLine[0].y+j;
		   PointInLine.push_back(fPoint(xj,yj));
		}	
		 for(r=R;r>=-R;r--) //���ƽ�������ڸ�������ֵ
		{
		   for(j=0;j!=N;j++)
		   {
			   dx=PointInLine[j].x-r;
			   dy=PointInLine[j].y;
			   PointsInRegion.push_back(fPoint(dx,dy));
		   }
		}
	}

		  for(m=0;m!=M;m++)//���ÿ����ĻҶ�ֵ
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
	//�����ݶȾ�ֵ-��׼�������ӣ��ԻҶȾ����е�(M-1)*(N-1)���ݶ�
	vector <vector<double>> All_Grad_dL;//dL������ݶ�
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

	vector <vector<double>> All_Grad_dT;//dT������ݶ�
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
	//����������ݶȺϲ�
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
					Each_Grad_dLdT.push_back (All_Grad_dT[i][m+j*(M-1)]);//��dT���������ݶȰ�dL�������У���dL������ݶȶ�Ӧ��
			    } 
	     }
		 LinesFeature.push_back (Each_Grad_dLdT);//�õ��ݶ���������
		 Each_Grad_dLdT.swap(vector<double>());
 	}
	    All_Grad_dL.swap(vector<vector<double>>());
	    All_Grad_dT.swap(vector<vector<double>>());

    //�����ݶ���������ľ�ֵ�ͱ�׼��
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
			EachMeanGrad.push_back (SumRowGray/(N-1));//�õ�ÿ�е��ݶȾ�ֵ
			double SumSqr=0.0;
			for(j=0;j!=N-1;j++)
			{
				SumSqr+=(LinesFeature[i][j+m*(N-1)]-SumRowGray/(N-1))*(LinesFeature[i][j+m*(N-1)]-SumRowGray/(N-1));
			}
		    double Std=sqrt(SumSqr);
			EachStd.push_back(Std);//�õ�ÿ�еı�׼��
		}
		AllMeanGrad.push_back (EachMeanGrad);//����ֱ�ߵ��ݶȾ�ֵ����
		AllStdGrad.push_back (EachStd);//����ֱ�ߵ��ݶȱ�׼�����
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
			AllMeanGrad[i][m]=AllMeanGrad[i][m]/ LengthMean;//��һ�����ݶȾ�ֵ
			AllStdGrad[i][m]=AllStdGrad[i][m]/LengthStd;//��һ���ı�׼��
		}
	}

	for(i=0;i!=Feature.size();i++)
	{
		for(m=0;m!=2*(M-1);m++)
		{
			AllMeanGrad[i].push_back (AllStdGrad[i][m]);//����׼�����ݶȾ�ֵ�ϲ�
		}
		LinesMSLD.push_back (AllMeanGrad[i]);//�õ�MSLD��������
	}
	AllMeanGrad.swap(vector<vector<double>>());
	AllStdGrad.swap(vector<vector<double>>());
	LinesFeature.swap(vector<vector<double>>());
}

