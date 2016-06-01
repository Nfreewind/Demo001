#include "stdafx.h"
#include <cv.h>
#include <highgui.h>
#include "function.h"
#include <math.h>
#include "constant.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "fuzzy_function.h"

void membership_function2(IplImage * read_image, int a, int b, int c, float r, float **membership_value)
{
	int i,j;
	int Xmn;
	
	float dividend,
		  divisor;
	

				
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){

			Xmn = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
			if(0 <= Xmn && Xmn < a)	
				membership_value[i][j] = 0;
			else if(a <= Xmn && Xmn < b){
				dividend = pow((float)(Xmn - a), r);
				divisor = pow((float)(b - a), (r-1)) * (float)(c - a);
				membership_value[i][j] = dividend / divisor;
			}
			else if(b <= Xmn && Xmn < c){
				dividend = pow((float)(c - Xmn), r);
				divisor = pow((float)(c - b), (r - 1)) * (float)(c - a);
				membership_value[i][j] =1.0 - (dividend / divisor);
			}
			else	membership_value[i][j] = 1.0;
	
		}
	}

	return;

}

void find_power_value_retinex(IplImage *read_image, float *fhist, int a, int b, int c, float r, float f, float **membership_value, float **window_fuzzy_entropy, float **power_value)
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
		dividend = pow((float)(gl - a),r);
		divisor = pow((float)(b - a), (r - 1)) * (float)(c - a);
		u_gl = dividend / divisor;
	}
	else if(b <= gl && gl < c){
		dividend = pow((float)(c - gl), (r));
		divisor = pow((float)(c - b), (r - 1)) * (float)(c - a);
		u_gl =1.0 - (dividend / divisor);
	}
	else	u_gl = 1.0;

	//先算出gh的membership_value
	if(0 <= gh && gh < a)	
		u_gh = 0;
	else if(a <= gh && gh < b){
		dividend = pow((float)(gh - a),r);
		divisor = pow((float)(b - a), (r - 1)) * (float)(c - a);
		u_gh = dividend / divisor;
	}
	else if(b <= gh && gh < c){
		dividend = pow((float)(c - gh), (r));
		divisor = pow((float)(c - b), (r - 1)) * (float)(c - a);
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


void defuzzification_retinex(IplImage *read_image, int a, int b, int c, float r, int Lmax, int Lmin, float **modified_mem_value, IplImage *out_image)
{
	int i,j;
	float dividend,	//被除數(分子)
		  divisor;	//除數(分母)
	float AAA,	//數學式的前式
		  BBB,	//數學式的中式
		  CCC;	//數學式的後式
	float C1, C2, C3;	//數學式後式的第1 2 3 項

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			if(modified_mem_value[i][j] == 0){
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = Lmin;

			}else if(0 < modified_mem_value[i][j] && modified_mem_value[i][j] <= (float)(b - a) / (float)(c - a)){
				dividend = Lmax - Lmin;
				divisor = c - a;
				C1 = pow((float)(b - a), (r - 1));
				C2 = (float)(c - a);
				C3 = modified_mem_value[i][j];


				AAA = Lmin;
				BBB = dividend / divisor;
				CCC = pow((C1 * C2 * C3), (1 / r));

				/*C1 = modified_mem_value[i][j];
				C2 = pow((float)(b - a), (r - 1));
				C3 = c - a;
				CCC = pow((C1 * C2 * C3), (1 / r));*/

				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = AAA + BBB * CCC;

			}else if((float)(b - a) / (float)(c - a) < modified_mem_value[i][j] && modified_mem_value[i][j] < 1.0){
				dividend = Lmax-Lmin;
				divisor = c - a;
				C1 = pow((float)(c - b), (r - 1));
				C2 = (float)(c - a);
				C3 = 1.0 - modified_mem_value[i][j];

				AAA = Lmin;
				BBB = dividend / divisor;;
				CCC = pow((C1 * C2 * C3),(1 / r));

				/*C1 = pow((float)(c - b), (r - 1));
				C2 = c - a;
				C3 = 1 - modified_mem_value[i][j];
				CCC = pow((C1 * C2 * C3), (1 / r));*/

				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = AAA + BBB * (c - a - CCC);
			}else{
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = Lmax;
			}
		}
	}

	return;

}

void retinex_contrast(IplImage *read_image, float **membership_value, float **membership_enhanced)
{
	int i,j;
	float fsum,
		  constant,
		  relation;
	int heightStep,
		widthStep;

	float a,b;
	float fmax,fmin;

	float **membership_input = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(歸屬值)
	for(i=0; i < read_image->height; i++){
		membership_input[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	float **membership_tmp = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(歸屬值)
	for(i=0; i < read_image->height; i++){
		membership_tmp[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	float **membership_output = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(歸屬值)
	for(i=0; i < read_image->height; i++){
		membership_output[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	fsum = 0.0;
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			fsum += membership_value[i][j];
		}
	}
	constant = fsum / (read_image->height * read_image->width);

	heightStep = read_image->height / 2;
	widthStep = read_image->width / 2;
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			membership_input[i][j] = constant;
			//printf("membership_value[i][j] = %f\n",membership_value[i][j]);
			//printf("membership_input[i][j] = %f\n\n",membership_input[i][j]);
		}
	}
	
	do{
		print_file(membership_input, read_image->height, read_image->width, "membership_input.txt");
		//先做水平
		for(i=0; i < read_image->height; i++){
			for(j=0; j < read_image->width - widthStep; j++){
				//relation = log10(membership_value[i][j + widthStep]) - log10(membership_value[i][j]);
				relation = log10((float)((uchar *)(read_image->imageData + i*read_image->widthStep))[j + widthStep]) - log10((float)((uchar *)(read_image->imageData + i*read_image->widthStep))[j]);
				membership_tmp[i][j] = membership_input[i][j] - relation;
				membership_tmp[i][j + widthStep] = membership_input[i][j + widthStep] + relation;
			}
		}
		print_file(membership_tmp, read_image->height, read_image->width, "membership_tmp.txt");

		//再做垂直
		for(j=0; j < read_image->width; j++){
			for(i=0; i < read_image->height - heightStep; i++){
				//relation = log10(membership_value[i + heightStep][j]) - log10(membership_value[i][j]);
				relation = log10((float)((uchar *)(read_image->imageData + (i + heightStep)*read_image->widthStep))[j]) - log10((float)((uchar *)(read_image->imageData + i*read_image->widthStep))[j]);
				membership_output[i][j] = membership_tmp[i][j] - relation;
				membership_output[i + heightStep][j] = membership_tmp[i + heightStep][j] + relation;
			}
		}
		print_file(membership_output, read_image->height, read_image->width, "membership_output.txt");

		//把輸出再當作輸入，繼續疊代
		for(i=0; i < read_image->height; i++){
			for(j=0; j < read_image->width; j++){
				membership_input[i][j] = membership_output[i][j];
			}
		}
		
		heightStep /= 2;
		widthStep /= 2;
	} while(heightStep != 0);

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			if(membership_output[i][j] < 0)
				membership_enhanced[i][j] = 0;
			else if(membership_output[i][j] > 1)
				membership_enhanced[i][j] = 1;
			else
				membership_enhanced[i][j] = membership_output[i][j];
		}
	}

	/*fmax = -999;
	fmin = 999;
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			if(membership_output[i][j] > fmax)
				fmax = membership_output[i][j];
			if(membership_output[i][j] < fmin)
				fmin = membership_output[i][j];
		}
	}
	a = (1.0 - 0.0) / (fmax - fmin);
	b = ((-1.0)*(1.0)*fmin) / (fmax - fmin);
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			membership_enhanced[i][j] = a * membership_output[i][j] + b;
		}
	}*/

	for(i=0; i < read_image->height; i++){
		free(membership_input[i]);
	}
	free(membership_input);
	
	for(i=0; i < read_image->height; i++){
		free(membership_tmp[i]);
	}
	free(membership_tmp);
	
	for(i=0; i < read_image->height; i++){
		free(membership_output[i]);
	}
	free(membership_output);

	return;
	
	
}

void fuzzy_retinex(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	int sum;
	int cunt;
	float fsum;
	float fmax, fmin;

	float f1 = 0.01, 
		  f2 = 0.5;
	float f = 0.005;
	int windowSize = 16;

	/*********************************Algorithm 1宣告****************************************/
	//目的: 決定a c
	int hist[256] = {0};		//直方圖統計(pixel)
	float fhist[256] = {0.0};	//直方圖統計(機率)
	int histIndex[256] = {0};	//1:Local max	2:peak
	int left, mid, right;	//尋找local max要用的label
	int avg_height;	//所有local max的平均值
	int Lmin, Lmax;	//Lmin:此張圖的最小灰階值	Lmax:此張圖的最大灰階值
	
	int g1,	//the first peak的灰階值
		gk;	//the last peak的灰階值
	int B1, B2;
	int a,b,c,bopt;

	float **membership_value = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(歸屬值)
	for(i=0; i < read_image->height; i++){
		membership_value[i] = (float *)malloc(sizeof(float) * read_image->width);
	}
	//int Xmn;	//讀取到的灰階值
	float H[256] = {0};	//fuzzy entropy
	/*float logA, logB;
	float dividend,	//被除數(分子)
		  divisor;	//除數(分母)
	float AAA,	//數學式的前式
		  BBB;	//數學式的後式*/

	float r;	
	float entropy_tmp;
	float entropy_max;
	float ropt;	//最佳r

	/*********************************Algorithm 3宣告****************************************/
	//目的: 1.決定低對比強度範圍 2.判斷是否要做enhancement
	
	
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

	/**************************************************/
	float **membership_enhanced = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體(改寫過的membership value)
	for(i=0; i < read_image->height; i++){
		membership_enhanced[i] = (float *)malloc(sizeof(float) * read_image->width);
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


	//4. Select peak
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

	/*******************決定ropt*******************************/
	
	entropy_max = -999.0;
	for(r = 0.1; r <= 10; r += 0.1){

		membership_function2(read_image, a, bopt, c, r, membership_value);	
		/*for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){

				Xmn = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
				if(0 <= Xmn && Xmn <= a)	
					membership_value[i][j] = 0;
				else if(a <= Xmn && Xmn <= bopt){
					dividend = pow((float)(Xmn - a), r);
					divisor = pow((float)(bopt - a), (r-1)) * (float)(c - a);
					membership_value[i][j] = dividend / divisor;
				}
				else if(bopt <= Xmn && Xmn <= c){
					dividend = pow((float)(Xmn - c), r);
					divisor = pow((float)(c - bopt), (r - 1)) * (float)(c - a);
					membership_value[i][j] =1.0 - (dividend / divisor);
				}
				else	membership_value[i][j] = 1.0;
		
			}
		}*/

		//計算fuzzy entropy
		
		entropy_tmp = entropy(read_image, membership_value);
		/*for(i=0; i < read_image->height; i++){
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
		}*/
		if(entropy_tmp > entropy_max){
			entropy_max = entropy_tmp;
			ropt = r;
		}
	}


	/*for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){

			Xmn = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
			if(0 <= Xmn && Xmn <= a)	
				membership_value[i][j] = 0;
			else if(a <= Xmn && Xmn <= bopt){
				dividend = pow((float)(Xmn - a), ropt);
				divisor = pow((float)(bopt - a), (ropt-1)) * (float)(c - a);
				membership_value[i][j] = dividend / divisor;
			}
			else if(bopt <= Xmn && Xmn <= c){
				dividend = pow((float)(Xmn - c), ropt);
				divisor = pow((float)(c - bopt), (ropt - 1)) * (float)(c - a);
				membership_value[i][j] =1.0 - (dividend / divisor);
			}
			else	membership_value[i][j] = 1.0;
	
		}
	}*/



	;
	;
	/*******************Algorithm 3*******************************/

	
	//要先建立membership function
	membership_function2(read_image, a, bopt, c, ropt, membership_value);
	print_file(membership_value, read_image->height, read_image->width, "membership_value.txt");
	
	//Find dege value
	/*edge_function(read_image, membership_value, edge_value);
	print_file(edge_value, read_image->height, read_image->width, "edge_value.txt");

	//2. 計算window下的entropy
	window_entropy(read_image, windowSize, membership_value, edge_value, window_fuzzy_entropy);
	print_file(window_fuzzy_entropy, read_image->height, read_image->width, "window_fuzzy_entropy.txt");


	//3. 決定power value
	find_power_value_retinex(read_image, &fhist[0], a, bopt, c, ropt, f, membership_value, window_fuzzy_entropy, power_value);
	//find_power_value(read_image, &fhist[0], a, bopt, c, f, membership_value, window_fuzzy_entropy, power_value);
	print_file(power_value, read_image->height, read_image->width, "power_value.txt");*/



	/*******************Algorithm 2*******************************/

	//1. construct membership(algorithm 3 有先做了)
	

	//2. Find dege value(algorithm 3 有先做了)
	//print_file(edge_value, read_image->height, read_image->width, "edge_value.txt");


	//3. 計算遮罩中的平均edge value
	/*find_avg_edge_value(read_image, windowSize, membership_value, edge_value, avg_edge_value);
	print_file(avg_edge_value, read_image->height, read_image->width, "avg_edge_value.txt");

	//4. Evaluate contrast membership value
	find_contrast_membership_value(read_image, membership_value, avg_edge_value, contrast_mem_value);
	print_file(contrast_mem_value, read_image->height, read_image->width, "contrast_mem_value.txt");

	//5. Transform the contrast
	find_transform_contrast_membership_value(read_image, contrast_mem_value, power_value, trans_contrast_mem_value);
	print_file(trans_contrast_mem_value, read_image->height, read_image->width, "trans_contrast_mem_value.txt");

	//6. Obtain modified membership value
	find_modified_membership_value(read_image, membership_value, avg_edge_value, trans_contrast_mem_value, modified_mem_value);
	print_file(modified_mem_value, read_image->height, read_image->width, "modified_mem_value.txt");*/


	retinex_contrast(read_image, membership_value, membership_enhanced);	//retinex的對比強化
	print_file(membership_enhanced, read_image->height, read_image->width, "membership_enhanced.txt");

	//7. Defuzzification(解模糊化)
	defuzzification_retinex(read_image, a, bopt, c, ropt, Lmax, Lmin, membership_enhanced, out_image);

	return;
}