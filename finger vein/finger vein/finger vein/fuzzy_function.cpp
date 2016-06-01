#include "stdafx.h"
#include <cv.h>
#include <highgui.h>
#include "function.h"
#include <math.h>
#include "constant.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "fuzzy_retinex.h"

void print_file(float **file_point, int height, int width, char *filename)
{
	FILE *fp = fopen(filename,"w");

	int i,j;

	for(i=0; i < height; i++){
		for(j=0; j < width; j++){
			fprintf(fp, "%f	", file_point[i][j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);

	return;

}

void membership_function1(IplImage *read_image, int a, int b, int c, float **membership_value)
{

	int i,j;
	int Xmn;
	float dividend,
		  divisor;

	for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				Xmn = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];

				if(0 <= Xmn && Xmn < a)	membership_value[i][j] = 0;
				else if(a <= Xmn && Xmn < b){
					dividend = ((float)Xmn - a) * ((float)Xmn - a);
					divisor = ((float)b - a) * ((float)c - a);
					membership_value[i][j] = dividend / divisor;
				}
				else if(b <= Xmn && Xmn < c){
					dividend = ((float)Xmn - c) * ((float)Xmn - c);
					divisor = ((float)c - b) * ((float)c - a);
					membership_value[i][j] =1.0 - (dividend / divisor);
				}
				else	membership_value[i][j] = 1.0;
			}
		}

	return;
	
}

float entropy(IplImage *read_image, float **membership_value)
{
	int i,j;
	float logA,
		  logB;
	float AAA,
		  BBB;
	float fsum;
	float entropy_value;


	fsum = 0;
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){

			logA = log10(membership_value[i][j]) / log10(2.0);
			logB = log10(1.0 - membership_value[i][j]) / log10(2.0);
				
			AAA = (-1.0) * membership_value[i][j] * logA;
			BBB = (1.0 - membership_value[i][j]) * logB;
				
			if(membership_value[i][j] == 1.0 || membership_value[i][j] == 0.0)
				fsum = fsum + 0.0;
			else
				fsum = fsum + (AAA - BBB);
		}
	}

	entropy_value = fsum / (read_image->height * read_image->width);

	return entropy_value;

}

void edge_function(IplImage *read_image, float **membership_value, float **edge_value)
{
	int i,j;
	int up, down, left, right;
	float XXX[9],YYY[9];
	float sumx,sumy;

	int cx_sobel[3][3] = {-1,	0,	1,
						  -2,	0,	2,
						  -1,	0,	1};

	int cy_sobel[3][3] = {-1,	-2,	-1,
						  0,	0,	0,
						  1,	2,	1};


	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){

			up = i-1;
			down = i+1;
			left = j-1;
			right = j+1;

			if( up < 0 )	up = 0;	
			if( down >= read_image->height ) down = read_image->height - 1;
			if( left < 0 ) left = 0;
			if( right >= read_image->width ) right = read_image->width - 1;

			XXX[0] = membership_value[up][left] * cx_sobel[0][0];
			XXX[1] = membership_value[up][j] * cx_sobel[0][1];
			XXX[2] = membership_value[up][right] * cx_sobel[0][2];

			XXX[3] = membership_value[i][left] * cx_sobel[1][0];
			XXX[4] = membership_value[i][j] * cx_sobel[1][1];
			XXX[5] = membership_value[i][right] * cx_sobel[1][2];

			XXX[6] = membership_value[down][left] * cx_sobel[2][0];
			XXX[7] = membership_value[down][j] * cx_sobel[2][1];
			XXX[8] = membership_value[down][right] * cx_sobel[2][2];



			YYY[0] = membership_value[up][left] * cy_sobel[0][0];
			YYY[1] = membership_value[up][j] * cy_sobel[0][1];
			YYY[2] = membership_value[up][right] * cy_sobel[0][2];

			YYY[3] = membership_value[i][left] * cy_sobel[1][0];
			YYY[4] = membership_value[i][j] * cy_sobel[1][1];
			YYY[5] = membership_value[i][right] * cy_sobel[1][2];

			YYY[6] = membership_value[down][left] * cy_sobel[2][0];
			YYY[7] = membership_value[down][j] * cy_sobel[2][1];
			YYY[8] = membership_value[down][right] * cy_sobel[2][2];

			sumx = XXX[0] + XXX[1] + XXX[2] + XXX[3] + XXX[4] + XXX[5] + XXX[6] + XXX[7] + XXX[8];
			sumy = YYY[0] + YYY[1] + YYY[2] + YYY[3] + YYY[4] + YYY[5] + YYY[6] + YYY[7] + YYY[8];

			edge_value[i][j] = sqrt(sumx*sumx + sumy*sumy);

			
		}
	}

	return;

}

void window_entropy(IplImage *read_image, int windowSize, float **membership_value, float **edge_value, float **window_fuzzy_entropy)
{
	int i,j;

	float fsum;
	float dividend,	//被除數(分子)
		  divisor;	//除數(分母)
	float AAA,	//數學式的前式
		  BBB;	//數學式的後式
	int subX, subY;	//遮罩中的X Y索引
	int tmpX, tmpY;	//遮罩中的X Y暫存索引

	float sum_window_entropy;

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){

			//先計算window下所有β的總和
			fsum = 0.0;
			for(subX = i - (windowSize / 2); subX < i + (windowSize / 2); subX++){
				for(subY = j - (windowSize / 2); subY < j + (windowSize / 2); subY++){
					tmpX = subX;
					tmpY = subY;
					if(tmpX < 0)
						tmpX = 0;
					else if(tmpX >= read_image->height)
						tmpX = read_image->height - 1;
					
					if(tmpY < 0)
						tmpY = 0;
					else if(tmpY >= read_image->width)
						tmpY = read_image->width - 1;

					fsum += (membership_value[tmpX][tmpY] * edge_value[tmpX][tmpY]);

				}
			}

			sum_window_entropy = 0.0;
			for(subX = i - (windowSize / 2); subX < i + (windowSize / 2); subX++){
				for(subY = j - (windowSize / 2); subY < j + (windowSize / 2); subY++){
					tmpX = subX;
					tmpY = subY;
					if(tmpX < 0)
						tmpX = 0;
					else if(tmpX >= read_image->height)
						tmpX = read_image->height - 1;
					
					if(tmpY < 0)
						tmpY = 0;
					else if(tmpY >= read_image->width)
						tmpY = read_image->width - 1;

					dividend = membership_value[tmpX][tmpY] * edge_value[tmpX][tmpY];
					divisor = fsum;

					if(dividend == 0 || divisor == 0)
						sum_window_entropy += 0;
					else{
						AAA = dividend / divisor;
						BBB = log10(AAA) / log10(2.0);
						sum_window_entropy += AAA * BBB;
					}

				}
			}
			window_fuzzy_entropy[i][j] =(-1.0) * sum_window_entropy / (log10((float)windowSize * (float)windowSize) / log10(2.0));

		}
	}
	return;
}

void find_power_value(IplImage *read_image, float *fhist, int a, int b, int c, float f, float **membership_value, float **window_fuzzy_entropy, float **power_value)
{
	int i,j;
	float fsum;

	float dividend,	//被除數(分子)
		  divisor;	//除數(分母)
	float AAA,	//數學式的前式
		  BBB;	//數學式的後式

	int Lmax, Lmin;
	float power_min, power_max;
	int gl, gh;
	float u_gl, u_gh;
	float entropy_max, entropy_min;

	//先找到Lmin和Lmax
	for(i=0; i < 256; i++){
		if(fhist[i] != 0){
			Lmin = i;
			break;
		}
	}
	for(i=255; i > 0; i--){
		if(fhist[i] != 0){
			Lmax = i;
			break;
		}
	}

	dividend = c - a;
	divisor = 2 * (Lmax - Lmin);
	power_min = dividend / divisor;
	power_max = 1.0;


	//1. 決定gl和gh
	fsum = 0.0;
	for(i = a; i <= c; i++){
		fsum += fhist[i];
		if(fsum > f){
			gl = i - 1;
			break;
		}
	}

	fsum = 0.0;
	for(i = c; i >= a; i--){
		fsum += fhist[i];
		if(fsum > f){
			gh = i + 1;
			break;
		}
	}

	
	//先算出gl的membership_value
	if(0 <= gl && gl < a)	
	u_gl = 0;
	else if(a <= gl && gl < b){
		dividend = ((float)gl - a) * ((float)gl - a);
		divisor = ((float)b - a) * ((float)c - a);
		u_gl = dividend / divisor;
	}
	else if(b <= gl && gl < c){
		dividend = ((float)gl - c) * ((float)gl - c);
		divisor = ((float)c - b) * ((float)c - a);
		u_gl =1.0 - (dividend / divisor);
	}
	else	u_gl = 1.0;

	//先算出gh的membership_value
	if(0 <= gh && gh <= a)	
		u_gh = 0;
	else if(a <= gh && gh <= b){
		dividend = ((float)gh - a) * ((float)gh - a);
		divisor = ((float)b - a) * ((float)c - a);
		u_gh = dividend / divisor;
	}
	else if(b <= gh && gh <= c){
		dividend = ((float)gh - c) * ((float)gh - c);
		divisor = ((float)c - b) * ((float)c - a);
		u_gh =1.0 - (dividend / divisor);
	}
	else	u_gh = 1.0;

	//找出entropy_max和entropy_min
	entropy_max = -999.0;
	entropy_min = 999.0;
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			if(window_fuzzy_entropy[i][j] > entropy_max)
				entropy_max = window_fuzzy_entropy[i][j];
			
			if(window_fuzzy_entropy[i][j] < entropy_min)
				entropy_min = window_fuzzy_entropy[i][j];
		}
	}

	//開始判斷
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){

			if(membership_value[i][j] < u_gl){
				dividend = (window_fuzzy_entropy[i][j] - entropy_min) * (power_max - power_min);
				divisor = entropy_max - entropy_min;

				AAA = u_gl / membership_value[i][j];
				BBB = power_min + (dividend / divisor);

				power_value[i][j] = AAA * BBB;

			}else if(u_gl <= membership_value[i][j] && membership_value[i][j] <= u_gh){
				dividend = (window_fuzzy_entropy[i][j] - entropy_min) * (power_max - power_min);
				divisor = entropy_max - entropy_min;

				AAA = power_min;
				BBB = dividend / divisor;

				power_value[i][j] = AAA + BBB;

			}else if(membership_value[i][j] > u_gh){
				dividend = (window_fuzzy_entropy[i][j] - entropy_min) * (power_max - power_min);
				divisor = entropy_max - entropy_min;

				AAA = membership_value[i][j] / u_gh;
				BBB = power_min + (dividend / divisor);

				power_value[i][j] = AAA * BBB;
			}

			
		}
	}

	return;


}

void find_avg_edge_value(IplImage *read_image, int windowSize, float **membership_value, float **edge_value, float **avg_edge_value)
{
	int i,j;
	float fsum;

	float dividend,	//被除數(分子)
		  divisor;	//除數(分母)
	int subX, subY;	//遮罩中的X Y索引
	int tmpX, tmpY;	//遮罩中的X Y暫存索引

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			
			fsum = 0.0;
			dividend = 0.0;	//分子
			divisor = 0.0;	//分母
			for(subX = i - (windowSize / 2); subX < i + (windowSize / 2); subX++){
				for(subY = j - (windowSize / 2); subY < j + (windowSize / 2); subY++){
					tmpX = subX;
					tmpY = subY;

					if(tmpX < 0)	
						tmpX = 0;
					else if(tmpX >= read_image->height)	
						tmpX = read_image->height - 1;

					if(tmpY < 0)	
						tmpY = 0;
					else if(tmpY >= read_image->width)
						tmpY = read_image->width - 1;
					
					dividend += (membership_value[tmpX][tmpY] * edge_value[tmpX][tmpY]);
					divisor += edge_value[tmpX][tmpY];
					
				}
			}
			if(divisor == 0)
				avg_edge_value[i][j] = 0;
			else
				avg_edge_value[i][j] = dividend / divisor;

		}

	}

	return;

}

void find_contrast_membership_value(IplImage *read_image, float **membership_value, float **avg_edge_value, float **contrast_mem_value)
{
	int i,j;
	float dividend,	//被除數(分子)
		  divisor;	//除數(分母)

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			dividend = abs(membership_value[i][j] - avg_edge_value[i][j]);
			divisor = abs(membership_value[i][j] + avg_edge_value[i][j]);

			contrast_mem_value[i][j] = dividend / divisor;
		}
	}

	return;

}

void find_transform_contrast_membership_value(IplImage *read_image, float **contrast_mem_value, float **power_value, float **trans_contrast_mem_value)
{
	int i,j;

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			trans_contrast_mem_value[i][j] = pow(contrast_mem_value[i][j], power_value[i][j]);
		}
	}
	return;
}

void find_modified_membership_value(IplImage *read_image, float **membership_value, float **avg_edge_value, float **trans_contrast_mem_value, float **modified_mem_value)
{
	int i,j;
	float dividend,	//被除數(分子)
		  divisor;	//除數(分母)
	float AAA,	//數學式的前式
		  BBB;	//數學式的後式

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			if(membership_value[i][j] <= avg_edge_value[i][j]){
				dividend = 1.0 - trans_contrast_mem_value[i][j];
				divisor = 1.0 + trans_contrast_mem_value[i][j];

				AAA = avg_edge_value[i][j];
				BBB = dividend / divisor;

				modified_mem_value[i][j] = AAA * BBB;
			}else{
				dividend = 1.0 + trans_contrast_mem_value[i][j];
				divisor = 1.0 - trans_contrast_mem_value[i][j];

				AAA = avg_edge_value[i][j];
				BBB = dividend / divisor;

				modified_mem_value[i][j] = AAA * BBB;
				if(modified_mem_value[i][j] > 1.0)
					modified_mem_value[i][j] = 1.0;
			}

		}
	}

	return;

}

void defuzzification(IplImage *read_image, int a, int b, int c, int Lmax, int Lmin, float **modified_mem_value, IplImage *out_image)
{
	int i,j;
	float dividend,	//被除數(分子)
		  divisor;	//除數(分母)
	float AAA,	//數學式的前式
		  BBB;	//數學式的後式

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			if(modified_mem_value[i][j] == 0){
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = Lmin;
			}else if(0 < modified_mem_value[i][j] && modified_mem_value[i][j] <= (float)(b - a) / (float)(c - a)){
				dividend = Lmax - Lmin;
				divisor = c - a;

				AAA = Lmin;
				BBB = dividend * sqrt(modified_mem_value[i][j] * (float)(b - a) * (float)(c - a)) / (float)(c - a);

				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = AAA + BBB;
			}else if((float)(b - a) / (float)(c - a) < modified_mem_value[i][j] && modified_mem_value[i][j] < 1.0){
				dividend = Lmax-Lmin;
				divisor = c - a;

				AAA = Lmin;
				BBB = (Lmax - Lmin) * (c - a - sqrt((1.0 - modified_mem_value[i][j]) * (float)(c - b) * (float)(c - a))) / (float)(c - a);

				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = AAA + BBB;
			}else{
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = Lmax;
			}
		}
	}

	return;

}


void fuzzy_contrast(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	int sum;
	int cunt;
	float fsum;
	float fmax, fmin;

	int windowSize = 16;

	/*********************************Algorithm 1宣告****************************************/
	//目的: 決定a c
	int hist[256] = {0};		//直方圖統計(pixel)
	float fhist[256] = {0.0};	//直方圖統計(機率)
	int histIndex[256] = {0};	//1:Local max	2:peak
	int left, mid, right;	//尋找local max要用的label
	int avg_height;	//所有local max的平均值
	int Lmin, Lmax;	//Lmin:此張圖的最小灰階值	Lmax:此張圖的最大灰階值
	float f1 = 0.01, 
		  f2 = 0.5;
	int g1,	//the first peak的灰階值
		gk;	//the last peak的灰階值
	int B1, B2;
	int a,b,c,bopt;

	float **membership_value = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(歸屬值)
	for(i=0; i < read_image->height; i++){
		membership_value[i] = (float *)malloc(sizeof(float) * read_image->width);
	}
	float H[256] = {0};	//fuzzy entropy


	/*********************************Algorithm 3宣告****************************************/
	//目的: 1.決定低對比強度範圍 2.判斷是否要做enhancement
	float f = 0.005;
	
	float **window_fuzzy_entropy = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(遮罩下的fuzzy entropy)
	for(i=0; i < read_image->height; i++){
		window_fuzzy_entropy[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	float **power_value = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(增強對比次方數) 1:
	for(i=0; i < read_image->height; i++){
		power_value[i] = (float *)malloc(sizeof(float) * read_image->width);
	}


	/*********************************Algorithm 2宣告****************************************/
	//目的: adaptive fuzzy contrast enhancement
	
	float **edge_value = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(邊界值)
	for(i=0; i < read_image->height; i++){
		edge_value[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	float **avg_edge_value = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(平均邊界值)
	for(i=0; i < read_image->height; i++){
		avg_edge_value[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	float **contrast_mem_value = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(對比歸屬值)
	for(i=0; i < read_image->height; i++){
		contrast_mem_value[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	float **trans_contrast_mem_value = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(轉換對比歸屬值)
	for(i=0; i < read_image->height; i++){
		trans_contrast_mem_value[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	float **modified_mem_value = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(改寫過的membership value)
	for(i=0; i < read_image->height; i++){
		modified_mem_value[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	

	
	/*******************Algorithm 1*******************************/

	//1. 直方圖統計
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			hist[ ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] ] ++;
		}

	}

	for(i=0; i < 256; i++){
		fhist[i] = (float)hist[i] / (read_image->height * read_image->width);
	}

	//先找到Lmin和Lmax
	for(i=0; i < 256; i++){
		if(hist[i] != 0){
			Lmin = i;
			break;
		}
	}
	for(i=255; i > 0; i--){
		if(hist[i] != 0){
			Lmax = i;
			break;
		}
	}

	//2. Find Local maxima
	for(i = Lmin; i <= Lmax; i++){
		left = i - 1;
		mid = i;
		right = i+1;

		if(left < 0)		left = i;
		if(right > Lmax)	right = Lmax;

		if((hist[mid] - hist[left]) > 0 && (hist[mid] - hist[right]) > 0){
			histIndex[i] = 1;
			printf("i = %d gray = %d\n", i, hist[i]);
		}
	}

	//3. Calculate average
	sum = 0;
	cunt = 0;
	for(i = Lmin; i <= Lmax; i++){
		if(histIndex[i] == 1){
			sum += hist[i];
			cunt ++;
		}
	}
	avg_height = sum / cunt;


	//4. Select peak(gl和gk)
	for(i = Lmin; i <= Lmax; i++){
		if(histIndex[i] == 1 && hist[i] > avg_height){
			g1 = i;
			break;
		}
	}
	for(i = Lmax; i >= Lmin; i--){
		if(histIndex[i] == 1 && hist[i] > avg_height){
			gk = i;
			break;
		}
	}

	//5. 決定B1、B2
	fsum = 0;
	for(i = Lmin; i <= Lmax; i++){
		fsum += fhist[i];
		if(fsum > f1){
			B1 = i-1;
			break;
		}
		
	}

	fsum = 0;
	for(i = Lmax; i >= Lmin; i--){
		fsum += fhist[i];
		if(fsum > f1){
			B2 = i+1;
			break;
		}
		
	}


	//6. 決定a,c
	a = (1.0 - f2) * (float)(g1 - Lmin) + (float)Lmin;
	if(a > B1)	
		a = B1;
	
	c = f2 * (float)(Lmax - gk) + (float)gk;
	if(c < B2)	
		c = B2;
	

	/*******************決定bopt*******************************/
	
	for(b = a+1; b < c; b++){
		
		membership_function1(read_image, a, b, c, membership_value);
		//membership function
		/*for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				Xmn = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];

				if(0 <= Xmn && Xmn <= a)	membership_value[i][j] = 0;
				else if(a <= Xmn && Xmn <= b){
					dividend = ((float)Xmn - a) * ((float)Xmn - a);
					divisor = ((float)b - a) * ((float)c - a);
					membership_value[i][j] = dividend / divisor;
				}
				else if(b <= Xmn && Xmn <= c){
					dividend = ((float)Xmn - c) * ((float)Xmn - c);
					divisor = ((float)c - b) * ((float)c - a);
					membership_value[i][j] =1.0 - (dividend / divisor);
				}
				else	membership_value[i][j] = 1.0;
			}
		}*/

		//Calculate entropy (Shannon function)

		//H[b] = entropy(read_image, membership_value);

		H[b] = entropy(read_image, membership_value);
		/*fsum = 0;
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){

				logA = log10(membership_value[i][j]) / log10(2.0);
				logB = log10(1.0 - membership_value[i][j]) / log10(2.0);
				
				AAA = (-1.0) * membership_value[i][j] * logA;
				BBB = (1.0 - membership_value[i][j]) * logB;
				
				if(membership_value[i][j] == 1.0 || membership_value[i][j] == 0.0)
					fsum = fsum + 0.0;
				else
					fsum = fsum + (AAA - BBB);
			}
		}
		H[b] = fsum / (read_image->height * read_image->width);*/
	}

	fmax = -999;
	for(i = a+1; i < c; i++){
		if(H[i] > fmax){
			fmax = H[i];
			bopt = i;
		}
	}


	/*******************Algorithm 3*******************************/

	
	//要先建立membership function
	membership_function1(read_image, a, bopt, c, membership_value);
	print_file(membership_value, read_image->height, read_image->width, "membership_value.txt");
	
	//Find dege value
	edge_function(read_image, membership_value, edge_value);
	print_file(edge_value, read_image->height, read_image->width, "edge_value.txt");

	//2. 計算window下的entropy
	window_entropy(read_image, windowSize, membership_value, edge_value, window_fuzzy_entropy);
	print_file(window_fuzzy_entropy, read_image->height, read_image->width, "window_fuzzy_entropy.txt");


	//3. 決定power value
	find_power_value(read_image, &fhist[0], a, bopt, c, f, membership_value, window_fuzzy_entropy, power_value);
	print_file(power_value, read_image->height, read_image->width, "power_value.txt");



	/*******************Algorithm 2*******************************/

	//1. construct membership(algorithm 3 有先做了)
	

	//2. Find dege value(algorithm 3 有先做了)
	//print_file(edge_value, read_image->height, read_image->width, "edge_value.txt");


	//3. 計算遮罩中的平均edge value
	find_avg_edge_value(read_image, windowSize, membership_value, edge_value, avg_edge_value);
	print_file(avg_edge_value, read_image->height, read_image->width, "avg_edge_value.txt");

	//4. Evaluate contrast membership value
	find_contrast_membership_value(read_image, membership_value, avg_edge_value, contrast_mem_value);
	print_file(contrast_mem_value, read_image->height, read_image->width, "contrast_mem_value.txt");

	//5. Transform the contrast
	find_transform_contrast_membership_value(read_image, contrast_mem_value, power_value, trans_contrast_mem_value);
	print_file(trans_contrast_mem_value, read_image->height, read_image->width, "trans_contrast_mem_value.txt");

	//6. Obtain modified membership value
	find_modified_membership_value(read_image, membership_value, avg_edge_value, trans_contrast_mem_value, modified_mem_value);
	print_file(modified_mem_value, read_image->height, read_image->width, "modified_mem_value.txt");

	//7. Defuzzification(解模糊化)
	defuzzification(read_image, a, bopt, c, Lmax, Lmin, modified_mem_value, out_image);


	return;
	


}