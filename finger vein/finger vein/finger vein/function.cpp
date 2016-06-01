#include "stdafx.h"
#include <cv.h>
#include <highgui.h>
#include "function.h"
#include <math.h>
#include "constant.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>






int is_overflow(int a)
{
	if(a > 255)
		a = 255;
	else if(a < 0)
		a = 0;

	return a;
}

/*
complex *Alloc1DArray(int size, int item_size)
{
	/* ------------------------------------------------------------
　　　作用：取得一維陣列動態記憶體
　　　參數：size = 一維陣列大小
			item_size = 陣列元素大小
　　　傳回：記憶體空間指標，若為NULL表記憶體不夠
	  限制：size*item_size必須小於65536
　　　------------------------------------------------------------ 
	int i, j;
	complex * ptr;

	ptr = (complex *) malloc(item_size * size);
	if (ptr == NULL){
		printf("記憶體不夠, 釋放原取得之記憶體!\n");
		free(ptr);
		return NULL;
	}
	
	return ptr;
}


complex **Alloc2DArray(int xsize,int ysize,int item_size)
{
	/* ------------------------------------------------------------
　　　作用：取得二維陣列動態記憶體
　　　參數：xsize, ysize = 二維陣列大小
			item_size = 陣列元素大小
　　　傳回：記憶體空間指標，若為NULL表記憶體不夠
	  限制：xsize*item_size必須小於65536
　　　------------------------------------------------------------ 
	
	int i, j;

	complex ** ptr;

	ptr = (complex**) malloc(sizeof(complex*) * ysize);

	if (ptr == NULL){
		printf("取得記憶體失敗!\n");	
		return NULL;
	}
	
	for (i=0; i<ysize; i++){
		ptr[i] = (void*) malloc(item_size * xsize);
		if (ptr[i] == NULL){
			printf("記憶體不夠, 釋放原取得之記憶體\n"); 
			for (j=0; j<i; j++){
				free(ptr[j]);
			}
			free(ptr);
			return NULL;
		}
	}
	return ptr;
}

void Free2DArray(int xsize,int ysize, complex **ptr)
{
	/* ------------------------------------------------------------
　　作用：釋放二維陣列動態記憶體
　　參數：xsize, ysize = 二維陣列大小
　　　　　ptr = 記憶體空間指標
　　備註：xsize實際上並沒用到，但為容易辨識使用，故加上此一參數
   ------------------------------------------------------------ 
	
	int i;

	if(ptr == NULL){
		printf("空指標!!\n");
		return;
	}
	for (i=0; i<ysize; i++) free(ptr[i]);
	free(ptr);
}*/

//按比例轉換灰階值
void scal_trans(int **read_image, IplImage *out_image)
{
	int i,j;
	int gL = 999999, gH = -9999999;

	for(i=0; i < out_image->height; i++){
		for(j=0; j < out_image->width; j++){
			if( read_image[i][j] < gL )
				gL = read_image[i][j];
			if( read_image[i][j] > gH)
				gH = read_image[i][j];
		}
	}

	for(i=0; i < out_image->height; i++){
		for(j=0; j < out_image->width; j++){
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = (255 * (read_image[i][j] - gL)) / (gH - gL);
		}
	}

	return;
	
}

//顯示影像
void imshow(IplImage *read_image, char *windowName)
{
	cvNamedWindow(windowName,0);	//創建一個新視窗，命名為"showpicture1"
	cvResizeWindow(windowName, read_image->width, read_image->height);	//設定視窗大小
	cvShowImage(windowName, read_image);	//顯示圖像
	cvWaitKey(0);

	cvDestroyWindow(windowName);	//關閉視窗

	return;
}

void U8toF64(IplImage *read_image, IplImage *out_image)
{
	int i,j;

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((float *)(out_image->imageData + i*out_image->widthStep))[j] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
		}
	}

	return;
}

void mat2grat(int **read_matrix, IplImage *out_image)
{
	int i,j;
	int min, max, shift, range;

	min = 99999999;
	max = -99999999;

	for(i = 0; i < out_image->height; i++){
		for(j = 0; j < out_image->width; j++){
			if( read_matrix[i][j] < min )	min = read_matrix[i][j];
			if( read_matrix[i][j] > max )	max = read_matrix[i][j];
		}

	}

	shift = abs(min);
	range = max - min;

	for(i = 0; i < out_image->height; i++){
		for(j = 0; j < out_image->width; j++){
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = (((read_matrix[i][j] + shift) * 255) / range);
		}

	}

	return;

}

void print_int_mat_txt(int **in, int height, int width, char *fileName)
{
	int i,j;
	FILE *fp = fopen(fileName,"w");

	for(i=0; i < height; i++){
		for(j=0 ; j < width; j++){
			fprintf(fp, "%d ", in[i][j]);
		}
		fprintf(fp, "\n");
		
	}

	fclose(fp);
	return;

}

/*void print_inttxt(int *in, char *fileName)
{
	int i;
	FILE *filterPtr = fopen(fileName,"W");

	for(i=0; i < 256; i++){
		fprintf(filterPtr, "%d ", in[i]);
	}

	fclose(filterPtr);
	return;
	
}*/

void complex_div(int dim, complex **dividend, complex **divisor, complex **quot)
{
	int i,j;
	float a,b,c,d;

	for(i=0; i < dim; i++){
		for(j=0; j < dim; j++){
			a = dividend[i][j].real;
			b = dividend[i][j].image;
			c = divisor[i][j].real;
			d = divisor[i][j].image;
			quot[i][j].real = ((a*c) + (b*d)) / (c*c + d*d);
			quot[i][j].image = ((-1)*(a*d) + (b*c)) / (c*c + d*d);
		}
	}

	return;
}


//轉成單通道灰階
void gray(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	int r,g,b;
	double y;

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			b=((uchar *)(read_image->imageData + i*read_image->widthStep))[j*read_image->nChannels + 0]; // B
			g=((uchar *)(read_image->imageData + i*read_image->widthStep))[j*read_image->nChannels + 1]; // G
			r=((uchar *)(read_image->imageData + i*read_image->widthStep))[j*read_image->nChannels + 2]; // R
			y=0.114*b + 0.587*g + 0.299*r;  										//灰階轉換公式
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j]=y;
		}

	}

	return ;
}

//二值化
void binarization(IplImage *read_image, IplImage *out_image, int TH)
{
	
	int i,j;
	

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			if( ( ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] ) > TH)
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 255;
			else
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 0;
		}

	}

	

	return ;
}

void binarization2(IplImage *read_image, IplImage *out_image, int TH1, int TH2)
{
	int i,j;
	

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			if( ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] >= TH1 && ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] <= TH2 )
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 255;
			else
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 0;
		}

	}

	

	return ;
}

//位元平面
void bit_plane(IplImage *read_image, IplImage *out_image, int bitNumber)
{
	/************************************************************************************************************************

	輸入: bitNumber 為想要取出的位元平面。假設bitNumber = 3，代表取出第三位元平面(c3)。最低位元平面為bitNumber = 0(c0)。

	*************************************************************************************************************************/
	int i,j,k;
	int value;
	int BCD[8] = {0};

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			value = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
			for(k=7; k >= 0; k--){
				BCD[k] = value % 2;
				value = value / 2;
			}


			for(k=7; k >= 0; k--){
				if((7-k) == bitNumber){
					if(BCD[k] == 1)
						((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 255;
					else
						((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 0;
				}
			}

		}

	}
}

//空間解析度
void im_resize(IplImage *read_image, IplImage *out_image, int resolution)
{
	/********************************************************************************************

	輸入: resolution 為想要調整的解析度。假設resolution = 128，代表調整為128 X 128 解析度

	*********************************************************************************************/
	int i,j;
	int tmpi,tmpj;
	int maskN;

	maskN = 256 / resolution;	//maskN代表是用N乘N的遮罩下去降解析度，也就是N乘N裡面的值都要相同
	for(i=0; i < read_image->height; i = i + maskN){
		for(j=0; j < read_image->width; j = j + maskN){
			
			
			for(tmpi = i; tmpi < i+maskN; tmpi++){
				for(tmpj = j; tmpj < j+maskN; tmpj++){
					((uchar *)(out_image->imageData + tmpi*out_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					
				}
			}

		}

	}
	return;
}

//均勻量化
void quantizatioin(IplImage *read_image, IplImage *out_image, int grayNumber)
{
	
	/*****************************************************************************************************************************

	輸入: grayNumber 為想要量化的灰階數。如果grayNumber = 128，代表使用128灰階混色，大部分影像是使用256灰階混色(grayNumber = 256)

	******************************************************************************************************************************/
	int value;	//影像piexl的值
	int unit;	//一等分為多少piexl值
	int range;	//量化的區間
	int i,j,k;

	unit = 256 / grayNumber;
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			
			value = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
			range = 0;
			for(k=0; k < grayNumber; k++){

				if(range <= value && value < (range + unit) ){
					((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = is_overflow(range + unit);
					break;
				}

				range = range + unit;
			}
			

		}

	}
	return;
}

//混色(半色調)
void dither(IplImage *read_image, IplImage *out_image)
{
	IplImage *dither_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);

	/**************************
		請決定混色矩陣的維度
	***************************/
	int scale = 2;

	int ditherMatrixD[2][2] = {0, 128,
							   192, 64};

	int ditherMatrixD2[4][4] = {0, 128, 32, 160,
							  192, 64, 224, 96,
							  48, 176, 16, 144,
							  240, 112, 208, 80};
	int i,j, tmpi, tmpj;
	int ditherIndexI = 0, ditherIndexJ = 0;
	
	for(i=0; i < read_image->height; i = i + scale){
		for(j=0; j < read_image->width; j = j + scale){
			
			
			for(tmpi = i; tmpi < i + scale; tmpi++){
				for(tmpj = j; tmpj < j + scale; tmpj++){
					((uchar *)(dither_image->imageData + tmpi*dither_image->widthStep))[tmpj] = ditherMatrixD[ditherIndexI][ditherIndexJ];
					ditherIndexJ++;
					
				}
				ditherIndexI++;
				ditherIndexJ = 0;
			}
			ditherIndexI = 0;

		}

	}
	

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			
			if( ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] > ((uchar *)(dither_image->imageData + i*dither_image->widthStep))[j])
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 255;
			else
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 0;

		}

	}

	cvReleaseImage(&dither_image);
	return;
}

//量化多層灰階混色
void multi_dither(IplImage *read_image, IplImage *out_image, int grayNumber)
{
	/*****************************************************************************************************************************

	輸入: grayNumber 為想要量化的灰階數。如果grayNumber = 128，代表使用128灰階混色，大部分影像是使用256灰階混色(grayNumber = 256)

	******************************************************************************************************************************/
	IplImage *normal_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	IplImage *dist_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);

	int i,j,k;
	int tmpi, tmpj;
	int normalMatrix[2][2];
	
	int dist;
	int unit;	//一等分為多少piexl值
	int indexI = 0, indexJ = 0;
	int TH;
	int value;

	unit = 256 / grayNumber;
	normalMatrix[0][0] = 0;
	normalMatrix[0][1] = (unit + unit) / 3;
	normalMatrix[1][0] = unit;
	normalMatrix[1][1] = unit / 3;
	
	
	for(i=0; i < read_image->height; i = i + 2){
		for(j=0; j < read_image->width; j = j + 2){
			
			
			for(tmpi = i; tmpi < i + 2; tmpi++){
				for(tmpj = j; tmpj < j + 2; tmpj++){
					((uchar *)(normal_image->imageData + tmpi*normal_image->widthStep))[tmpj] = normalMatrix[indexI][indexJ];
					indexJ++;
					
				}
				indexI++;
				indexJ = 0;
			}
			indexI = 0;

		}

	}

	
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			value = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
			TH = ( ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] / unit ) * unit;

			if( (value - TH) > ((uchar *)(normal_image->imageData + i*normal_image->widthStep))[j])
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = is_overflow(TH + unit);
			else
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = is_overflow(TH);

		}

	}

	
	cvReleaseImage(&normal_image);
	cvReleaseImage(&dist_image);
	return;
}

//誤差擴散法
void error_diff(IplImage *read_image, IplImage *out_image)
{
	IplImage *error_image = cvCreateImage(cvSize(read_image->width + 2, read_image->height + 2), IPL_DEPTH_8U,	1);

	int i,j;
	int e;	//誤差矩陣灰階值和量化值的誤差

	//將誤差矩陣初始化為0
	for(i=0; i < error_image->height; i++){
		for(j=0; j < error_image->width; j++){

			((uchar *)(error_image->imageData + i*error_image->widthStep))[j] = 0;

		}

	}

	//Copy原影像(n x n)至誤差矩陣(n+2 x n+2)
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j+1] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];

		}

	}

	
	for(i=1; i < error_image->height - 1; i++){
		for(j=1; j < error_image->width - 1; j++){
			if(((uchar *)(error_image->imageData + i*error_image->widthStep))[j] > 128){	//比較誤差矩陣灰階值和128，如果大於128，代表量化值為255
				((uchar *)(out_image->imageData + (i-1)*out_image->widthStep))[j-1] = 255;	//將誤差矩陣量化後，輸出影像
				e = ((uchar *)(error_image->imageData + i*error_image->widthStep))[j] - 255;	//計算誤差矩陣灰階值和量化值的誤差
			}else{
				((uchar *)(out_image->imageData + (i-1)*out_image->widthStep))[j-1] = 0;	
				e = ((uchar *)(error_image->imageData + i*error_image->widthStep))[j];
			}

			//調整誤差矩陣
			((uchar *)(error_image->imageData + i*error_image->widthStep))[j+1] = is_overflow(((uchar *)(error_image->imageData + i*error_image->widthStep))[j+1] + (e*7)/16);
			((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j-1] = is_overflow(((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j-1] + (e*3)/16);
			((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j] = is_overflow(((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j] + (e*5)/16);
			((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j+1] = is_overflow(((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j+1] + e/16);
			
		}
	}
	cvReleaseImage(&error_image);
	return;
}

//影像加法
void im_add(IplImage *read_image, IplImage *out_image, int number)
{
	int i,j;

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = is_overflow( ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] + number );

		}

	}

	return;
}

//影像減法
void im_sub(IplImage *read_image, IplImage *out_image, int number)
{
	int i,j;

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = is_overflow( ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] - number );

		}

	}

	return;
}

//影像乘法
void im_multiply(IplImage *read_image, IplImage *out_image, int number)
{
	int i,j;

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = is_overflow( ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] * number );

		}

	}

	return;
}

//影像除法
void im_divide(IplImage *read_image, IplImage *out_image, int number)
{
	int i,j;

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = is_overflow( ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] / number );

		}

	}

	return;
}

//影像補色
void im_complement(IplImage *read_image, IplImage *out_image)
{
	int i,j;

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = is_overflow( 255 - ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] );

		}

	}

	return;
}

//統計直方圖
void calc_Hist(IplImage *read_image)
{
	FILE *calHist_fp;
	int grayNumberArray[256] = {0};
	int i,j;

	calHist_fp = fopen("out.txt", "w");

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			grayNumberArray[ ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] ] ++;
		}

	}

	for(i = 0; i < 256; i ++){
		fprintf(calHist_fp, "%d ", grayNumberArray[i]);

	}

	fclose(calHist_fp);
	return;
}

//直方圖擴展法
void hist_expand(IplImage *read_image)
{
	FILE *out_fp;
	//FILE *out2_fp;
	int grayNumberArray[256] = {0};
	int expandArray[256] = {0};
	int i,j;
	int flagMin, flagMax ;
	int expandMin = 0, expandMax = 255;

	out_fp = fopen("out.txt", "w");
	//out2_fp = fopen("out2.txt", "w");

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			grayNumberArray[ ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] ] ++;
		}

	}

	for(i=0; i < 256; i++){
		if( grayNumberArray[i] != 0){
			flagMin = i;
			break;
		}
	}

	for(i=255; i >= 0; i--){
		if( grayNumberArray[i] != 0){
			flagMax = i;
			break;
		}
	}

	for(i=flagMin; i <= flagMax; i++){
		if( grayNumberArray[i] ){
			j = ( (expandMax - expandMin) * (i - flagMin) ) / (flagMax - flagMin) + expandMin;
			expandArray[j] = grayNumberArray[i];
		}
		
	}

	

	for(i = 0; i < 256; i ++){
		fprintf(out_fp, "%d ", grayNumberArray[i]);
		//fprintf(out2_fp, "%d ", expandArray[i]);
		
	}

	fclose(out_fp);
	//fclose(out2_fp);
	return;
}

//直方圖擴展HS
void hist_stretch(IplImage *read_image, IplImage *out_image)
{
	

	int i,j;

	int Hist[256] = {0};	//統計灰階出現次數矩陣
	int LUT[256] = {0};
	int greyMin, greyMax;
	int NumPixels;
	int sum;

	NumPixels = read_image->width * read_image->height;

	//統計灰階出現次數
	greyMin = 9999;
	greyMax = -9999;
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			if(((uchar *)(read_image->imageData + i*read_image->widthStep))[j] < greyMin)
				greyMin = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];

			if(((uchar *)(read_image->imageData + i*read_image->widthStep))[j] > greyMax)
				greyMax = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];

			Hist[ ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] ] ++;
		}

	}

	

	/*********使用累積比例延展***********/
	;
	sum = 0;
	for(i=0; i < 256; i++){
		if(Hist[i]){
			sum += Hist[i];
			LUT[i] = (255 * sum) / NumPixels;
		}
	}

	/*********使用線段比例延展***********/
	;
	/*
	for(i = greyMin; i <= greyMax; i++){
		LUT[i] = (255 * (i - greyMin)) / (greyMax - greyMin);
	}*/



	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = LUT[ ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] ];
		}

	}

	float accu[256] = {0.0};

	for(i=0; i < 256; i++)
		Hist[i] = 0;

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			Hist[ ((uchar *)(out_image->imageData + i*out_image->widthStep))[j] ] ++;
		}

	}

	for(i=0; i < 256; i++){
		accu[i] = (float)Hist[i] /  (float)NumPixels;
	}

	FILE *fp = fopen("hist.txt","w");

	for(i=0; i < 256; i++){
		fprintf(fp, "%f ", accu[i]);
		
	}

	fclose(fp);

	
	
	return;
}


void sort(IplImage *read_image, xyw* buf, int low)
{
	int i, j, num, wt;
	int k;
	xyw tmp;
	int up, down, right, left;
	int A[8] = {0};

	
	num = 0;
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			if( ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] == low ){
				up = i-1;
				down = i+1;
				right = j+1;
				left = j-1;

				if(up < 0)	up = 0;
				if(down >= read_image->height)	down = read_image->height - 1;
				if(right >= read_image->width)	right = read_image->width - 1;
				if(left < 0)	left = 0;

				A[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left];
				A[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j];
				A[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right];

				A[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left];
				A[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right];

				A[5] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left];
				A[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j];
				A[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right];

				wt = 0;
				for(k = 0; k < 8; k++)
					wt += A[k];

				buf[num].w = wt;
				buf[num].x = j;
				buf[num].y = i;
				num++;
				
			}
		}
	}

	for(i = 0; i < num - 1; i++){
		for(j = i+1; j < num; j++){
			if(buf[i].w <= buf[j].w){
				tmp.w = buf[i].w;
				tmp.x = buf[i].x;
				tmp.y = buf[i].y;

				buf[i].w = buf[j].w;
				buf[i].x = buf[j].x;
				buf[i].y = buf[j].y;

				buf[j].w = tmp.w;
				buf[j].x = tmp.x;
				buf[j].y = tmp.y;

			}
		}
	}

	return;
}

void hist_equ(IplImage *read_image, IplImage *out_image)
{
	int i,j,x,y,sum;
	int delt;	//根據周圍像素的階層選擇像素數
	int low,high;	//處理階層的範圍
	int BinIncr;	//等化後，1個灰階階層的像素數
	long Hist[256] = {0};
	xyw *buf = (xyw *)malloc(sizeof(xyw) * 100000);
	IplImage *buf_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);

	BinIncr = (read_image->width * read_image->height) / 256;
	high = 255;
	low = 255;
	
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 0;
			((uchar *)(buf_image->imageData + i*buf_image->widthStep))[j] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
			Hist[ ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] ]++;
		}
	}

	for(i = 255; i > 0; i--){
		for(sum = 0; sum < BinIncr; low--)
			sum += Hist[low];

		low ++;

		delt = Hist[low] - (sum - BinIncr);
		sort(buf_image, &buf[0], low);	//按照周圍像素的灰階由高向低排列變化

		if(low < high){
			for(y = 0; y < read_image->height; y++){
				for(x = 0; x < read_image->width; x++){
					if( ((uchar *)(read_image->imageData + y*read_image->widthStep))[x] >= (low + 1) &&  ((uchar *)(read_image->imageData + y*read_image->widthStep))[x] <= high)
						((uchar *)(out_image->imageData + y*out_image->widthStep))[x] = i;
				}
			}
		}

		for(j = 0; j < delt; j++){
			((uchar *)(out_image->imageData + (buf[j].y)*out_image->widthStep))[buf[j].x] = i;
			((uchar *)(buf_image->imageData + (buf[j].y)*buf_image->widthStep))[buf[j].x] = 0;
		}
		Hist[low] = Hist[low] - delt;
		high = low;

	}


	for(i=0; i < 256; i++)
		Hist[i] = 0;

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			Hist[ ((uchar *)(out_image->imageData + i*out_image->widthStep))[j] ] ++;
		}

	}

	FILE *fp = fopen("hist.txt","w");

	for(i=0; i < 256; i++){
		fprintf(fp, "%d ", Hist[i]);
		
	}

	fclose(fp);

	return;

}

//LUT直方圖等化
void LUT_hist_eq(IplImage *read_image, IplImage *out_image)
{
	FILE *out_fp;
	

	int grayNumberArray[256] = {0};	//統計灰階出現次數陣列
	int cdfArray[256] = {0};
	int equArray[256] = {0};	//直方圖等化陣列
	int LUT[256] = {0};
	int i,j;
	int flagMin, flagMax ;	//flagMin為灰階最小值發生的地方，flagMax為灰階最大值發生的地方
	int cdfNumber = 0, cdfMin, cdfMax;	//cdfNumber為用來累積計算cdf的數字，cdfMin為累積分布函數最小值，cdfMax為累積分布函數最大值(其實也就是piexl總數)
	int changedGray;

	out_fp = fopen("out.txt", "w");
	

	//課本enginner.tif練習
	;
	/***********************************
	FILE *out_divi_fp;
	out_divi_fp = fopen("out_divi.txt", "w");
	IplImage *tmp2_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	cvNamedWindow("showpicture3",0);	//創建一個新視窗，命名為"showpicture1"
	cvResizeWindow("showpicture3", read_image->width, read_image->height);	//設定視窗大小
	cvNamedWindow("showpicture4",0);	//創建一個新視窗，命名為"showpicture1"
	cvResizeWindow("showpicture4", read_image->width, read_image->height);	//設定視窗大小
	cvNamedWindow("showpicture5",0);	//創建一個新視窗，命名為"showpicture1"
	cvResizeWindow("showpicture5", read_image->width, read_image->height);	//設定視窗大小
	int tmpCalcArray[256] = {0};

	im_divide(read_image, tmp2_image, 4);
	cvShowImage("showpicture3", read_image);	//顯示圖像
	cvShowImage("showpicture4", tmp2_image);	//顯示圖像
	cvWaitKey(0);
	//統計灰階出現次數
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			tmpCalcArray[ ((uchar *)(tmp2_image->imageData + i*tmp2_image->widthStep))[j] ] ++;
		}

	}

	for(i = 0; i < 256; i ++){
		fprintf(out_divi_fp, "%d ", tmpCalcArray[i]);
		
	}
	
	fclose(out_divi_fp);


	/************************************/

	//統計灰階出現次數
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			grayNumberArray[ ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] ] ++;
		}

	}

	//計算機率分布函數pdf
	for(i=0; i < 256; i++){
		cdfNumber = cdfNumber + grayNumberArray[i];
		cdfArray[i] = cdfNumber;
		
	}

	//找到灰階最小值發生的地方、累積分布函數最小值
	for(i=0; i < 256; i++){
		if( grayNumberArray[i] ){
			flagMin = i;
			cdfMin = cdfArray[i];
			break;
		}
	}

	//找到灰階最大值發生的地方、累積分布函數最大值
	for(i=255; i >= 0; i--){
		if( grayNumberArray[i] ){
			flagMax = i;
			cdfMax = cdfArray[i];
			break;
		}
	}

	//只對感興趣的地方做直方圖等化
	for(i = flagMin; i <= flagMax; i++){
		if(grayNumberArray[i]){
			changedGray = (int)((((float)cdfArray[i] - cdfMin) * 255) / ((float)cdfMax - cdfMin) + 0.5);
			LUT[i] = changedGray;
		}
	}

	//查詢LUT表格，找到相對應的灰階值
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = LUT[ ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] ];
		}

	}

	
	//統計查詢LUT等化後的出現次數(也就是直方圖等化)
	for(i=0; i < out_image->height; i++){
		for(j=0; j < out_image->width; j++){

			equArray[ ((uchar *)(out_image->imageData + i*out_image->widthStep))[j] ] ++;
		}

	}

	for(i = 0; i < 256; i ++){
		fprintf(out_fp, "%d ", equArray[i]);
		
	}

	fclose(out_fp);
	


	return;
}


//平均濾波器
void avg_filter(IplImage *read_image, IplImage *out_image, int dim)
{
	
	/*****************************************************************************************************************************

	輸入: shape 為指定影像邊緣部分的處理方法。-1:忽略邊緣 0:補零 1:鏡射
		  dim 為指定遮罩大小。若dim = 5，則產生5x5遮罩

	******************************************************************************************************************************/

	//int** mask;
	int sum;
	int i,j,tmpi,tmpj;
	int row,col;
	int filDexI,filDexJ;

	IplImage *copy_image = cvCreateImage(cvSize(read_image->width + (2*(dim / 2)), read_image->height + (2*(dim / 2))), IPL_DEPTH_8U,	1);

	float **avgFilter = (float **)malloc(sizeof(float *) * dim);
	for(i=0; i < dim; i++){
		avgFilter[i] = (float *)malloc(sizeof(float) * dim);
	}
	

	
	
	//初始化遮罩
	for(i=0; i < dim; i++){
		for(j=0; j < dim; j++){
			avgFilter[i][j] = (float)1 / (dim*dim);
		}
	}
	



	for(i=0; i < copy_image->height; i++){
		for(j=0; j < copy_image->width; j++){

			if(i < (dim/2) ){
				row = 0;
			}else if(i >= (copy_image->height - (dim/2)) ){
				row = read_image->height - 1;
			}else{
				row = i - (dim/2);
			}

			if( j < (dim/2) ){
				col = 0;
			}else if(j >= (copy_image->width - (dim/2)) ){
				col = read_image->width - 1;
			}else{
				col = j - (dim/2);
			}

			((uchar *)(copy_image->imageData + i*copy_image->widthStep))[j] = ((uchar *)(read_image->imageData + (row)*read_image->widthStep))[col];
		}

	}

	for(i = (dim / 2); i < copy_image->height - (dim / 2); i++){		
		for(j = (dim / 2); j < copy_image->width - (dim / 2); j++){		
			sum = 0;
			filDexI = 0;

			for(tmpi = i - (dim / 2); tmpi <= i + (dim / 2); tmpi++){	
				filDexJ = 0;
				for(tmpj = j - (dim / 2); tmpj <= j + (dim / 2); tmpj++){			
					sum = sum + ((uchar *)(copy_image->imageData + (tmpi)*copy_image->widthStep))[tmpj] * avgFilter[filDexI][filDexJ];
					filDexJ ++;

				}
				filDexI ++;
			}

			((uchar *)(out_image->imageData + (i - (dim / 2))*out_image->widthStep))[j - (dim / 2)] = is_overflow(sum);	
		}
	}


	cvReleaseImage(&copy_image);		//釋放圖像記憶體

	for (i=0; i < dim; i++) free(avgFilter[i]);
	free(avgFilter);

	return;


}

//Laplacian濾波器
void lapla_filter(IplImage *read_image, IplImage *out_image, int dim)
{
	float Lap[3][3] = {0.1667, 0.6667, 0.1667,
					   0.6667, -3.3333, 0.6667,
					   0.1667, 0.6667, 0.1667};
	int sum;
	int i,j,tmpi,tmpj;
	int maskIndexI = 0, maskIndexJ = 0;

	IplImage *tmpA_image = cvCreateImage(cvSize(read_image->width + (2*(dim / 2)), read_image->height + (2*(dim / 2))), IPL_DEPTH_8U,	1);
	

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			if(i == 0 && j == 0){	/*****左上角****/
				for(tmpi = i; tmpi < (dim / 2); tmpi++){
					for(tmpj = j; tmpj < (dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*tmpA_image->widthStep))[j];
					}
				}
			}else if(i == 0 && j == read_image->width - 1){	/*****右上角****/
				for(tmpi = 0; tmpi < (dim / 2); tmpi++){
					for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					}
				}
			}else if(i == read_image->height - 1 && j == 0){	/*****左下角****/
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					for(tmpj = 0; tmpj < (dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*tmpA_image->widthStep))[j];
					}
				}
			}else if(i == read_image->height - 1 && j == read_image->width - 1){	/*****右下角****/
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					}
				}
			}else if(i == 0){	/*****上邊****/
				for(tmpi = i; tmpi < (dim / 2); tmpi++){
					((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[j] = ((uchar *)(read_image->imageData + (0)*read_image->widthStep))[j];
				}
			}else if(j == read_image->width - 1){	/*****右邊****/
				for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
					((uchar *)(tmpA_image->imageData + i*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[read_image->width - 1];
				}
			}else if(i == read_image->height - 1){	/*****下邊****/
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[j] = ((uchar *)(read_image->imageData + (read_image->height - 1)*read_image->widthStep))[j];
				}
			}else if(j == 0){	/*****左邊****/
				for(tmpj = 0; tmpj < (dim / 2); tmpj++){
					((uchar *)(tmpA_image->imageData + i*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[0];
				}
			}else{
				((uchar *)(tmpA_image->imageData + (i + (dim / 2))*tmpA_image->widthStep))[j + (dim / 2)] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
			}
		}
	}
			

			
	for(i = (dim / 2); i < tmpA_image->height - (dim / 2); i++){
		for(j = (dim / 2); j < tmpA_image->width - (dim / 2); j++){
			sum = 0;
			
			maskIndexI = 0;
			for(tmpi = i - (dim / 2); tmpi <= i + (dim / 2); tmpi++){
				maskIndexJ = 0;
				for(tmpj = j - (dim / 2); tmpj <= j + (dim / 2); tmpj++){
					sum = sum + ((uchar *)(tmpA_image->imageData + (tmpi)*tmpA_image->widthStep))[tmpj] * Lap[maskIndexI][maskIndexJ];
					maskIndexJ++;
				}
				maskIndexI++;
				
			}
			((uchar *)(out_image->imageData + (i - (dim / 2))*out_image->widthStep))[j - (dim / 2)] = is_overflow(sum);
		}
	}

	
	cvReleaseImage(&tmpA_image);		//釋放圖像記憶體
	
	return;
}

//LoG濾波
void LoG_filter(IplImage *read_image, IplImage *out_image, int dim)
{
	float LoG[5][5] = {0.0448, 0.0468, 0.0564, 0.0468, 0.0448,
					   0.0468, 0.3167, 0.7146, 0.3167, 0.0468,
					   0.0564, 0.7146, -4.9048, 0.7146, 0.0564,
					   0.0468, 0.3167, 0.7146, 0.3167, 0.0468,
					   0.0448, 0.0468, 0.0564, 0.0468, 0.0448};
	int sum;
	int i,j,tmpi,tmpj;
	int maskIndexI = 0, maskIndexJ = 0;

	IplImage *tmpA_image = cvCreateImage(cvSize(read_image->width + (2*(dim / 2)), read_image->height + (2*(dim / 2))), IPL_DEPTH_8U,	1);
	

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			if(i == 0 && j == 0){	/*****左上角****/
				for(tmpi = i; tmpi < (dim / 2); tmpi++){
					for(tmpj = j; tmpj < (dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*tmpA_image->widthStep))[j];
					}
				}
			}else if(i == 0 && j == read_image->width - 1){	/*****右上角****/
				for(tmpi = 0; tmpi < (dim / 2); tmpi++){
					for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					}
				}
			}else if(i == read_image->height - 1 && j == 0){	/*****左下角****/
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					for(tmpj = 0; tmpj < (dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*tmpA_image->widthStep))[j];
					}
				}
			}else if(i == read_image->height - 1 && j == read_image->width - 1){	/*****右下角****/
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					}
				}
			}else if(i == 0){	/*****上邊****/
				for(tmpi = i; tmpi < (dim / 2); tmpi++){
					((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[j] = ((uchar *)(read_image->imageData + (0)*read_image->widthStep))[j];
				}
			}else if(j == read_image->width - 1){	/*****右邊****/
				for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
					((uchar *)(tmpA_image->imageData + i*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[read_image->width - 1];
				}
			}else if(i == read_image->height - 1){	/*****下邊****/
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[j] = ((uchar *)(read_image->imageData + (read_image->height - 1)*read_image->widthStep))[j];
				}
			}else if(j == 0){	/*****左邊****/
				for(tmpj = 0; tmpj < (dim / 2); tmpj++){
					((uchar *)(tmpA_image->imageData + i*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[0];
				}
			}else{
				((uchar *)(tmpA_image->imageData + (i + (dim / 2))*tmpA_image->widthStep))[j + (dim / 2)] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
			}
		}
	}
			

			
	for(i = (dim / 2); i < tmpA_image->height - (dim / 2); i++){
		for(j = (dim / 2); j < tmpA_image->width - (dim / 2); j++){
			sum = 0;
			
			maskIndexI = 0;
			for(tmpi = i - (dim / 2); tmpi <= i + (dim / 2); tmpi++){
				maskIndexJ = 0;
				for(tmpj = j - (dim / 2); tmpj <= j + (dim / 2); tmpj++){
					sum = sum + ((uchar *)(tmpA_image->imageData + (tmpi)*tmpA_image->widthStep))[tmpj] * LoG[maskIndexI][maskIndexJ];
					maskIndexJ++;
				}
				maskIndexI++;
				
			}
			((uchar *)(out_image->imageData + (i - (dim / 2))*out_image->widthStep))[j - (dim / 2)] = is_overflow(sum);
		}
	}

	
	cvReleaseImage(&tmpA_image);		//釋放圖像記憶體
	return;
}

//高通濾波器
void HPF_filter(IplImage *read_image, IplImage *out_image, int dim)
{
	int Filter[3][3] = {1, -2, 1,
					   -2, 4, -2, 
					   1, -2, -1};
	int sum;
	int i,j,tmpi,tmpj;
	int up, right, left, down;
	int A[9] = {0};
	
	
	//int **tmpB_image;

	//配置記憶體空間給比例轉換後的image
	;
	/*
	tmpB_image = (int **)malloc(sizeof(int *) * read_image->height);	
	for(i=0; i < read_image->height; i++){
		tmpB_image[i] = (int *)malloc(sizeof(int ) * read_image->width);
	}*/

	;
	//copy一張比較大的圖
	;
	/*
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			if(i == 0 && j == 0){	//*****左上角****
				for(tmpi = i; tmpi < (dim / 2); tmpi++){
					for(tmpj = j; tmpj < (dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*tmpA_image->widthStep))[j];
					}
				}
			}else if(i == 0 && j == read_image->width - 1){	//****右上角***
				for(tmpi = 0; tmpi < (dim / 2); tmpi++){
					for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					}
				}
			}else if(i == read_image->height - 1 && j == 0){	//*****左下角****
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					for(tmpj = 0; tmpj < (dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*tmpA_image->widthStep))[j];
					}
				}
			}else if(i == read_image->height - 1 && j == read_image->width - 1){	//*****右下角****
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					}
				}
			}else if(i == 0){	//*****上邊****
				for(tmpi = i; tmpi < (dim / 2); tmpi++){
					((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[j] = ((uchar *)(read_image->imageData + (0)*read_image->widthStep))[j];
				}
			}else if(j == read_image->width - 1){	//*****右邊****
				for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
					((uchar *)(tmpA_image->imageData + i*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[read_image->width - 1];
				}
			}else if(i == read_image->height - 1){	//*****下邊****
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[j] = ((uchar *)(read_image->imageData + (read_image->height - 1)*read_image->widthStep))[j];
				}
			}else if(j == 0){	*****左邊****
				for(tmpj = 0; tmpj < (dim / 2); tmpj++){
					((uchar *)(tmpA_image->imageData + i*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[0];
				}
			}else{
				((uchar *)(tmpA_image->imageData + (i + (dim / 2))*tmpA_image->widthStep))[j + (dim / 2)] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
			}
		}
	}*/
			
	
	for(i = 0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){

			up = i-1;
			down = i+1;
			right = j+1;
			left = j-1;

			if(up < 0)	up = 0;
			if(down >= read_image->height)	down = read_image->height - 1;
			if(right >= read_image->width)	right = read_image->width - 1;
			if(left < 0)	left = 0;


			A[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * Filter[0][0];
			A[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * Filter[0][1];
			A[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * Filter[0][2];

			A[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * Filter[1][0];
			A[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * Filter[1][1];
			A[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * Filter[1][2];

			A[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * Filter[2][0];
			A[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * Filter[2][1];
			A[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * Filter[2][2];

			sum = A[0] + A[1] + A[2] + A[3] + A[4] + A[5] + A[6] + A[7] + A[8] ;
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = is_overflow(sum);
		}
	}

	
	
	
	//scal_trans(tmpB_image, out_image);	//按比例轉換灰階值

	return;
}

//去銳利化
void unsharp(IplImage *read_image, IplImage *out_image, int gain, float k)
{
	/*****************************************************************************************************************************

	輸入: gain 為調整原影像的放大比率
		  k 為調整比例 k < 1

	******************************************************************************************************************************/
	float unit[3][3] = {0, 0, 0,
					    0, 1, 0,
						0, 0, 0};
	float avgFilter[3][3] = {0.1111, 0.1111, 0.1111,
							 0.1111, 0.1111, 0.1111,
							 0.1111, 0.1111, 0.1111};
	float unsharpFilter[3][3] = {0};

	int i,j;
	int up, right, left, down;
	int A[9] = {0};
	int sum = 0;

	
	

	//去銳利化影像 = 直接算出去銳利畫遮罩，用遮罩下去濾波
	;
	
	for(i=0; i < 3; i++){
		for(j=0; j < 3; j++){
			unsharpFilter[i][j] = ( gain * unit[i][j] ) - ( (1 / k) * avgFilter[i][j]);
		}
	}

	for(i = 0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){

			up = i-1;
			down = i+1;
			right = j+1;
			left = j-1;

			if(up < 0)	up = 0;
			if(down >= read_image->height)	down = read_image->height - 1;
			if(right >= read_image->width)	right = read_image->width - 1;
			if(left < 0)	left = 0;


			A[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * unsharpFilter[0][0];
			A[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * unsharpFilter[0][1];
			A[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * unsharpFilter[0][2];

			A[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * unsharpFilter[1][0];
			A[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * unsharpFilter[1][1];
			A[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * unsharpFilter[1][2];

			A[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * unsharpFilter[2][0];
			A[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * unsharpFilter[2][1];
			A[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * unsharpFilter[2][2];

			sum = A[0] + A[1] + A[2] + A[3] + A[4] + A[5] + A[6] + A[7] + A[8] ;
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = is_overflow(sum);
		}
	}
	//if(!cvSaveImage("unsharp2.bmp",out_image)) printf("Could not save: %s\n", "out.bmp");
	
	
	;
	
	//去銳利化影像 = 原始影像放大gain倍 - 平均濾波影像放大1/k倍()
	;
	/*****************************************************
	IplImage *avg_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);

	for(i = 0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){

			up = i-1;
			down = i+1;
			right = j+1;
			left = j-1;

			if(up < 0)	up = 0;
			if(down >= read_image->width)	down = read_image->width - 1;
			if(right >= read_image->width)	right = read_image->width - 1;
			if(left < 0)	left = 0;


			A[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * avgFilter[0][0];
			A[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * avgFilter[0][1];
			A[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * avgFilter[0][2];

			A[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * avgFilter[1][0];
			A[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * avgFilter[1][1];
			A[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * avgFilter[1][2];

			A[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * avgFilter[2][0];
			A[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * avgFilter[2][1];
			A[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * avgFilter[2][2];

			sum = A[0] + A[1] + A[2] + A[3] + A[4] + A[5] + A[6] + A[7] + A[8] ;
			((uchar *)(avg_image->imageData + (i)*avg_image->widthStep))[j] = is_overflow(sum);
		}
	}

	for(i = 0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = is_overflow( ( ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * gain ) - ( ((float) 1/k) *  ((uchar *)(avg_image->imageData + (i)*avg_image->widthStep))[j]) );
		}
	}

	cvNamedWindow("read",0);	//創建一個新視窗，命名為"showpicture1"
	cvResizeWindow("read", read_image->width, read_image->height);	//設定視窗大小
	cvShowImage("read", read_image);	//顯示圖像

	cvNamedWindow("avg",0);	//創建一個新視窗，命名為"showpicture1"
	cvResizeWindow("avg", read_image->width, read_image->height);	//設定視窗大小
	cvShowImage("avg", avg_image);	//顯示圖像

	cvNamedWindow("unsharp",0);	//創建一個新視窗，命名為"showpicture1"
	cvResizeWindow("unsharp", read_image->width, read_image->height);	//設定視窗大小
	cvShowImage("unsharp", out_image);	//顯示圖像

	if(!cvSaveImage("unsharp1.bmp",out_image)) printf("Could not save: %s\n", "out.bmp");

	cvWaitKey(0);

	cvDestroyWindow("read");	//關閉視窗
	cvDestroyWindow("unsharp");	//關閉視窗

	cvReleaseImage(&avg_image);		//釋放圖像記憶體

	****************************************************/
	

	return;


}

//高增幅濾波
void hb_filter(IplImage *read_image, IplImage *out_image, float A)
{
	float avgFilter[3][3] = {0.1111, 0.1111, 0.1111,
							 0.1111, 0.1111, 0.1111,
							 0.1111, 0.1111, 0.1111};
	int i,j;
	int up, right, left, down;
	float AA[9] = {0};
	float sum = 0;

	IplImage *avg_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			up = i-1;
			down = i+1;
			right = j+1;
			left = j-1;

			if(up < 0)	up = 0;
			if(down >= read_image->height)	down = read_image->height - 1;
			if(right >= read_image->width)	right = read_image->width - 1;
			if(left < 0)	left = 0;


			AA[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * avgFilter[0][0];
			AA[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * avgFilter[0][1];
			AA[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * avgFilter[0][2];

			AA[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * avgFilter[1][0];
			AA[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * avgFilter[1][1];
			AA[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * avgFilter[1][2];

			AA[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * avgFilter[2][0];
			AA[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * avgFilter[2][1];
			AA[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * avgFilter[2][2];

			sum = AA[0] + AA[1] + AA[2] + AA[3] + AA[4] + AA[5] + AA[6] + AA[7] + AA[8] ;
			((uchar *)(avg_image->imageData + (i)*avg_image->widthStep))[j] = is_overflow(sum);
		}
	}

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = is_overflow( ((A / (2*A - 1)) * ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] ) - ( ((1-A) / (2*A - 1)) * ((uchar *)(avg_image->imageData + (i)*avg_image->widthStep))[j]) );
		}
	}

	cvReleaseImage(&avg_image);		//釋放圖像記憶體
	return;
}

//最大最小濾波器
void max_min_filter(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	int m,n;
	int up, right, left, down;
	int AA[9] = {0};
	int tmp;

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			up = i-1;
			down = i+1;
			right = j+1;
			left = j-1;

			if(up < 0)	up = 0;
			if(down >= read_image->height)	down = read_image->height - 1;
			if(right >= read_image->width)	right = read_image->width - 1;
			if(left < 0)	left = 0;


			AA[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left];
			AA[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j];
			AA[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right];

			AA[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left];
			AA[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
			AA[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right];

			AA[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left];
			AA[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j];
			AA[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right];

			for(m=0; m < 9-1; m++){
				for(n=0; n < 9-1-m; n++){
					if(AA[n] > AA[n+1]){
						tmp = AA[n+1];
						AA[n+1] = AA[n];
						AA[n] = tmp;
					}
				}
			}

			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = AA[0];
		}
	}

	return;
}

//中位濾波器
void mid_filter(IplImage *read_image, IplImage *out_image, int dim)
{
	/*****************************************************************************************************************************

	輸入: shape 為指定影像邊緣部分的處理方法。-1:忽略邊緣 0:補零 1:鏡射
		  dim 為指定遮罩大小。若dim = 5，則產生5x5遮罩

	******************************************************************************************************************************/

	int tmp;
	int A,B;
	int i,j,tmpi,tmpj;
	int row, col;
	int m,n;
	int filDex;

	IplImage *copy_image = cvCreateImage(cvSize(read_image->width + (2*(dim / 2)), read_image->height + (2*(dim / 2))), IPL_DEPTH_8U,	1);

	int *midFilter = (int *)malloc(sizeof(int ) * (dim*dim));
	
	

	
	
	//初始化遮罩
	for(i=0; i < (dim*dim); i++){
		midFilter[i] = 0;
	}
	



	for(i=0; i < copy_image->height; i++){
		for(j=0; j < copy_image->width; j++){

			if(i < (dim/2) ){
				row = 0;
			}else if(i >= (copy_image->height - (dim/2)) ){
				row = read_image->height - 1;
			}else{
				row = i - (dim/2);
			}

			if( j < (dim/2) ){
				col = 0;
			}else if(j >= (copy_image->width - (dim/2)) ){
				col = read_image->width - 1;
			}else{
				col = j - (dim/2);
			}

			((uchar *)(copy_image->imageData + i*copy_image->widthStep))[j] = ((uchar *)(read_image->imageData + (row)*read_image->widthStep))[col];
		}

	}


	for(i = (dim / 2); i < copy_image->height - (dim / 2); i++){		
		for(j = (dim / 2); j < copy_image->width - (dim / 2); j++){		
			
			filDex = 0;
			for(tmpi = i - (dim / 2); tmpi <= i + (dim / 2); tmpi++){	
				for(tmpj = j - (dim / 2); tmpj <= j + (dim / 2); tmpj++){			
					midFilter[filDex] = ((uchar *)(copy_image->imageData + (tmpi)*copy_image->widthStep))[tmpj];
					filDex ++;

				}
			}

			for(m = 0; m < (dim*dim) - 1; m++){
				for(n=0; n < (dim*dim) - m - 1; n++){
					if(midFilter[n] > midFilter[n+ 1 ] ){
						tmp = midFilter[n + 1];
						midFilter[n + 1] = midFilter[n];
						midFilter[n] = tmp;
					}
				}	
			}

			if( (dim % 2) == 0 ){
				A = midFilter[((dim*dim) / 2) - 1];
				B = midFilter[(dim*dim) / 2];
				((uchar *)(out_image->imageData + (i - (dim / 2))*out_image->widthStep))[j - (dim / 2)] = (A+B) / 2;
			}else{
				((uchar *)(out_image->imageData + (i - (dim / 2))*out_image->widthStep))[j - (dim / 2)] = midFilter[(dim*dim) / 2];
			}
		}
	}

	cvReleaseImage(&copy_image);		//釋放圖像記憶體

	free(midFilter);

	return;
}

void outliers_filter(IplImage *read_image, IplImage *out_image, int dim)
{
	int i,j,tmpi,tmpj;
	int row, col;
	int sum;
	float avg;

	float D = 0.2;

	IplImage *copy_image = cvCreateImage(cvSize(read_image->width + (2*(dim / 2)), read_image->height + (2*(dim / 2))), IPL_DEPTH_8U,	1);

	for(i=0; i < copy_image->height; i++){
		for(j=0; j < copy_image->width; j++){

			if(i < (dim/2) ){
				row = 0;
			}else if(i >= (copy_image->height - (dim/2)) ){
				row = read_image->height - 1;
			}else{
				row = i - (dim/2);
			}

			if( j < (dim/2) ){
				col = 0;
			}else if(j >= (copy_image->width - (dim/2)) ){
				col = read_image->width - 1;
			}else{
				col = j - (dim/2);
			}

			((uchar *)(copy_image->imageData + i*copy_image->widthStep))[j] = ((uchar *)(read_image->imageData + (row)*read_image->widthStep))[col];
		}

	}

	for(i = (dim / 2); i < copy_image->height - (dim / 2); i++){		
		for(j = (dim / 2); j < copy_image->width - (dim / 2); j++){		
			
			sum = 0;
			for(tmpi = i - (dim / 2); tmpi <= i + (dim / 2); tmpi++){	
				for(tmpj = j - (dim / 2); tmpj <= j + (dim / 2); tmpj++){
					if(tmpi != i || tmpj != j)
						sum =sum + ((uchar *)(copy_image->imageData + (tmpi)*copy_image->widthStep))[tmpj];

				}
			}
			avg = (float)sum / (dim*dim - 1);

			if( abs( ((uchar *)(copy_image->imageData + (i)*copy_image->widthStep))[j] - avg) > D )
				 ((uchar *)(out_image->imageData + (i - (dim / 2))*out_image->widthStep))[j - (dim / 2)] = (int)avg;
			else
				((uchar *)(out_image->imageData + (i - (dim / 2))*out_image->widthStep))[j - (dim / 2)] = ((uchar *)(copy_image->imageData + (i)*copy_image->widthStep))[j];

		}
	}

	cvReleaseImage(&copy_image);		//釋放圖像記憶體

	return;
}

void adaptive_filter(IplImage *read_image, IplImage *out_image, int dim)
{
	int i,j,tmpi,tmpj;
	int row, col;

	float A;
	float sumG, avgG, varG;
	float sumF, avgF, varF;


	IplImage *copy_image = cvCreateImage(cvSize(read_image->width + (2*(dim / 2)), read_image->height + (2*(dim / 2))), IPL_DEPTH_8U,	1);

	for(i=0; i < copy_image->height; i++){
		for(j=0; j < copy_image->width; j++){

			if(i < (dim/2) ){
				row = 0;
			}else if(i >= (copy_image->height - (dim/2)) ){
				row = read_image->height - 1;
			}else{
				row = i - (dim/2);
			}

			if( j < (dim/2) ){
				col = 0;
			}else if(j >= (copy_image->width - (dim/2)) ){
				col = read_image->width - 1;
			}else{
				col = j - (dim/2);
			}

			((uchar *)(copy_image->imageData + i*copy_image->widthStep))[j] = ((uchar *)(read_image->imageData + (row)*read_image->widthStep))[col];
		}

	}

	sumG = 0;
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			sumG = sumG + ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
		}

	}
	avgG = sumG / (read_image->height * read_image->width);

	sumG = 0;
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			A = ( ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] - avgG ) * ( ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] - avgG );
			sumG = sumG + A;
		}

	}
	varG = sumG / (read_image->height * read_image->width - 1);

	////////////////////////////////////////////////////////////////////////////////////////////

	for(i = (dim / 2); i < copy_image->height - (dim / 2); i++){		
		for(j = (dim / 2); j < copy_image->width - (dim / 2); j++){		
			
			sumF = 0;
			for(tmpi = i - (dim / 2); tmpi <= i + (dim / 2); tmpi++){	
				for(tmpj = j - (dim / 2); tmpj <= j + (dim / 2); tmpj++){
					sumF =sumF + ((uchar *)(copy_image->imageData + (tmpi)*copy_image->widthStep))[tmpj];

				}
			}
			avgF = sumF / (dim*dim);

			sumF = 0;
			for(tmpi = i - (dim / 2); tmpi <= i + (dim / 2); tmpi++){	
				for(tmpj = j - (dim / 2); tmpj <= j + (dim / 2); tmpj++){
					A = ( ((uchar *)(copy_image->imageData + (tmpi)*copy_image->widthStep))[tmpj] - avgF ) * ( ((uchar *)(copy_image->imageData + (tmpi)*copy_image->widthStep))[tmpj] - avgF );
					sumF =sumF + A;

				}
			}
			varF = sumF / (dim*dim - 1);

			((uchar *)(out_image->imageData + (i - (dim / 2))*out_image->widthStep))[j - (dim / 2)] = avgF + (varF * (((uchar *)(copy_image->imageData + (i)*copy_image->widthStep))[j] - avgF)) / (varF + varG);

		}
	}

	cvReleaseImage(&copy_image);		//釋放圖像記憶體

	return;
}

void band_reject_folter(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	float distance;
	float **filter = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體;
	for(i=0; i < in_no; i++){
		filter[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	complex **read_DFT = (complex **)malloc(sizeof(complex *) * read_image->height);	//取得二維陣列動態記憶體;
	for(i=0; i < in_no; i++){
		read_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	complex **read_pass = (complex **)malloc(sizeof(complex *) * read_image->height);	//取得二維陣列動態記憶體;
	for(i=0; i < in_no; i++){
		read_pass[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}


	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			distance = sqrt( ((float)i-127)*((float)i-127) + ((float)j-127)*((float)j-127) );
			if( 47 < distance && distance < 51)
				filter[i][j] = 0.0;
			else
				filter[i][j] = 1.0;
		}
	}

	DFT2(read_image, read_DFT);	//二維傅立葉轉換
	

	for(i=0;i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			read_pass[i][j].real = read_DFT[i][j].real * filter[i][j];
			read_pass[i][j].image = read_DFT[i][j].image * filter[i][j];
		}
	}
	fftshow(256, read_DFT,"read_DFT");	//顯示頻譜圖
	IDFT2(read_pass, out_image);	//二維傅立葉「反」轉換

	for (i=0; i < read_image->width; i++) free(filter[i]);
	free(filter);

	for (i=0; i < read_image->width; i++) free(read_DFT[i]);
	free(read_DFT);

	for (i=0; i < read_image->width; i++) free(read_pass[i]);
	free(read_pass);

	return;
}

void notch_filter(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	float distance;
	float **filter = (float **)malloc(sizeof(float *) * read_image->height);	//取得二維陣列動態記憶體;
	for(i=0; i < in_no; i++){
		filter[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	complex **read_DFT = (complex **)malloc(sizeof(complex *) * read_image->height);	//取得二維陣列動態記憶體;
	for(i=0; i < in_no; i++){
		read_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	complex **read_pass = (complex **)malloc(sizeof(complex *) * read_image->height);	//取得二維陣列動態記憶體;
	for(i=0; i < in_no; i++){
		read_pass[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}


	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			if(i == 155 || i == 101 || j == 169 || j == 88)
				filter[i][j] = 0.0;
			else
				filter[i][j] = 1.0;
			
		}
	}

	DFT2(read_image, read_DFT);	//二維傅立葉轉換
	

	for(i=0;i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			read_pass[i][j].real = read_DFT[i][j].real * filter[i][j];
			read_pass[i][j].image = read_DFT[i][j].image * filter[i][j];
		}
	}
	fftshow(256, read_pass,"read_pass");	//顯示頻譜圖
	IDFT2(read_pass, out_image);	//二維傅立葉「反」轉換

	for (i=0; i < read_image->width; i++) free(filter[i]);
	free(filter);

	for (i=0; i < read_image->width; i++) free(read_DFT[i]);
	free(read_DFT);

	for (i=0; i < read_image->width; i++) free(read_pass[i]);
	free(read_pass);

	return;

}

//ROI濾波
void ROI_foiter(IplImage *read_image, IplImage *out_image, int startI, int endI, int startJ, int endJ, int filterMode)
{
	float avgMask[3][3] = {0.1111, 0.1111, 0.1111,
							 0.1111, 0.1111, 0.1111,
							 0.1111, 0.1111, 0.1111};

	float unsharpMask[3][3] = {-0.3333, -0.3333, -0.3333,
								 -0.3333, 3.6667, -0.3333,
								 -0.3333, -0.3333, -0.3333};

	float LapMask[3][3] = {0.1667,    0.6667,    0.1667,
						   0.6667,   -3.3333,    0.6667,
						   0.1667,    0.6667,    0.1667};

	float LoGMask[5][5] = {0.0448,    0.0468,    0.0564,    0.0468,    0.0448,
							 0.0468,    0.3167,    0.7146,    0.3167,    0.0468,
							 0.0564,    0.7146,   -4.9048,    0.7146,    0.0564,
							 0.0468,    0.3167,    0.7146,    0.3167,    0.0468,
							 0.0448,    0.0468,    0.0564,    0.0468,    0.0448};

	int i,j,k;
	int up, right, left, down, U2, R2, D2, L2;
	float AA[9] = {0};
	float BB[25] = {0};
	float sum = 0;

	switch(filterMode){
		case 1:	//平均濾波器

			for(i=startI; i <= endI; i++){
				for(j=startJ; j <= endJ; j++){
					up = i-1;
					down = i+1;
					right = j+1;
					left = j-1;

					if(up < startI)	up = startI;
					if(down >= endI)	down = endI - 1;
					if(right >= endJ)	right = endJ - 1;
					if(left < startJ)	left = startJ;


					AA[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * avgMask[0][0];
					AA[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * avgMask[0][1];
					AA[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * avgMask[0][2];

					AA[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * avgMask[1][0];
					AA[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * avgMask[1][1];
					AA[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * avgMask[1][2];

					AA[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * avgMask[2][0];
					AA[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * avgMask[2][1];
					AA[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * avgMask[2][2];

					sum = AA[0] + AA[1] + AA[2] + AA[3] + AA[4] + AA[5] + AA[6] + AA[7] + AA[8] ;
					((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = is_overflow(sum);
				}
			}
			break;

		case 2:
			
			for(i=startI; i <= endI; i++){
				for(j=startJ; j <= endJ; j++){
					up = i-1;
					down = i+1;
					right = j+1;
					left = j-1;

					if(up < startI)	up = startI;
					if(down >= endI)	down = endI - 1;
					if(right >= endJ)	right = endJ - 1;
					if(left < startJ)	left = startJ;


					AA[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * unsharpMask[0][0];
					AA[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * unsharpMask[0][1];
					AA[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * unsharpMask[0][2];

					AA[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * unsharpMask[1][0];
					AA[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * unsharpMask[1][1];
					AA[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * unsharpMask[1][2];

					AA[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * unsharpMask[2][0];
					AA[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * unsharpMask[2][1];
					AA[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * unsharpMask[2][2];

					sum = AA[0] + AA[1] + AA[2] + AA[3] + AA[4] + AA[5] + AA[6] + AA[7] + AA[8] ;
					((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = is_overflow(sum);
				}
			}
			break;

		case 3:
			for(i=startI; i <= endI; i++){
				for(j=startJ; j <= endJ; j++){
					up = i-1;
					down = i+1;
					right = j+1;
					left = j-1;

					if(up < startI)	up = startI;
					if(down >= endI)	down = endI - 1;
					if(right >= endJ)	right = endJ - 1;
					if(left < startJ)	left = startJ;


					AA[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * LapMask[0][0];
					AA[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * LapMask[0][1];
					AA[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * LapMask[0][2];

					AA[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * LapMask[1][0];
					AA[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * LapMask[1][1];
					AA[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * LapMask[1][2];

					AA[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * LapMask[2][0];
					AA[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * LapMask[2][1];
					AA[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * LapMask[2][2];

					sum = AA[0] + AA[1] + AA[2] + AA[3] + AA[4] + AA[5] + AA[6] + AA[7] + AA[8] ;
					((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = is_overflow(sum);
				}
			}	
			break;
		
		case 4:

			for(i=startI; i <= endI; i++){
				for(j=startJ; j <= endJ; j++){
					up = i-1;
					down = i+1;
					right = j+1;
					left = j-1;
					L2 = j - 2;
					U2 = i - 2;
					R2 = j + 2;
					D2 = i + 2;

					if(up < startI)	up = startI;
					if(down >= endI)	down = endI - 1;
					if(right >= endJ)	right = endJ - 1;
					if(left < startJ)	left = startJ;

					if(U2 < startI)	U2 = startI;
					if(D2 >= endI)	D2 = endI - 2;
					if(R2 >= endJ)	R2 = endJ - 2;
					if(L2 < startJ)	L2 = startJ;


					BB[0] = ((uchar *)(read_image->imageData + (U2)*read_image->widthStep))[L2] * LoGMask[0][0];
					BB[1] = ((uchar *)(read_image->imageData + (U2)*read_image->widthStep))[left] * LoGMask[0][1];
					BB[2] = ((uchar *)(read_image->imageData + (U2)*read_image->widthStep))[j] * LoGMask[0][2];
					BB[3] = ((uchar *)(read_image->imageData + (U2)*read_image->widthStep))[right] * LoGMask[0][3];
					BB[4] = ((uchar *)(read_image->imageData + (U2)*read_image->widthStep))[R2] * LoGMask[0][4];

					BB[5] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[L2] * LoGMask[1][0];
					BB[6] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * LoGMask[1][1];
					BB[7] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * LoGMask[1][2];
					BB[8] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * LoGMask[1][3];
					BB[9] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[R2] * LoGMask[1][4];

					BB[10] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[L2] * LoGMask[2][0];
					BB[11] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * LoGMask[2][1];
					BB[12] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * LoGMask[2][2];
					BB[13] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * LoGMask[2][3];
					BB[14] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[R2] * LoGMask[2][4];

					BB[15] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[L2] * LoGMask[3][0];
					BB[16] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * LoGMask[3][1];
					BB[17] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * LoGMask[3][2];
					BB[18] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * LoGMask[3][3];
					BB[19] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[R2] * LoGMask[3][4];

					BB[20] = ((uchar *)(read_image->imageData + (D2)*read_image->widthStep))[L2] * LoGMask[4][0];
					BB[21] = ((uchar *)(read_image->imageData + (D2)*read_image->widthStep))[left] * LoGMask[4][1];
					BB[22] = ((uchar *)(read_image->imageData + (D2)*read_image->widthStep))[j] * LoGMask[4][2];
					BB[23] = ((uchar *)(read_image->imageData + (D2)*read_image->widthStep))[right] * LoGMask[4][3];
					BB[24] = ((uchar *)(read_image->imageData + (D2)*read_image->widthStep))[R2] * LoGMask[4][4];

					for(k=0; k < 25; k++){
						sum = BB[k];
					}
					
					((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = is_overflow(sum);
				}
			}

			if(!cvSaveImage("out.bmp",out_image)) printf("Could not save: %s\n", "out.bmp");
			break;

		default:
			printf("FilterMode is failed!!\n\n");
			break;

	}

	return;

}




/*
complex in_data[data_no] = {2, 0,
							3, 0,
							4, 0,
							5, 0,
							6, 0,
							7, 0,
							8, 0,
							1, 0};

complex in_dataA[M] = {1, 0,
							  2, 0,
							  3, 0,
							  4, 0};
complex in_dataB[N] = {5, 0,
							  6, 0,
							  7, 0,
							  8, 0};




in_data[0].real = 1;	in_data[0].image = 0;
in_data[1].real = 2;	in_data[1].image = 0;
in_data[2].real = 3;	in_data[2].image = 0;
in_data[3].real = 4;	in_data[3].image = 0;
//in_data[4].real = 6;	in_data[4].image = 0;
//in_data[5].real = 7;	in_data[5].image = 0;
//in_data[6].real = 8;	in_data[6].image = 0;
//in_data[7].real = 1;	in_data[7].image = 0;*/


//印出複數陣列
void print_complex(int data_no, complex *in_data)
{
	int i;
	for(i=0; i < data_no; i++){
		printf("%f + %fi\n", in_data[i].real, in_data[i].image);
	}
	printf("\n");
	return;
}

//一維數列平移
void shift_array(int data_no, complex *in_data, complex *out_data)
{
	int i;
	int sign = 1;

	out_data[0].real = in_data[0].real;
	out_data[0].image = in_data[0].image;
	for(i=1; i < data_no; i++){
		sign = sign * (-1);
		out_data[i].real = in_data[i].real * sign;
		out_data[i].image = in_data[i].image * sign;
	}

	return;
}

//一維傅立葉轉換
void DFT1(int data_no, complex *in_data, complex *out_data)
{
	
	float angle_step;	//計算-2*PI / T
	float w_cos, w_sin;
	float angle;
	int x,u;

	

	angle_step =  (((float)-2*PI) / data_no);
	for(u = 0; u < data_no; u++){
		out_data[u].real = in_data[0].real;
		out_data[u].image = in_data[0].image;

		for(x = 1; x < data_no; x ++){
			angle = (float) u * x * angle_step;
			w_cos = cos(angle);
			w_sin = sin(angle);
			out_data[u].real = out_data[u].real + w_cos * in_data[x].real - w_sin * in_data[x].image;
			out_data[u].image = out_data[u].image +  w_cos * in_data[x].image + w_sin * in_data[x].real;
		}

	}

	return;
}



//一維「反」傅立葉轉換
void iDFT1(int data_no, complex *in_data, complex *out_data)
{
	
	float angle_step;	//計算-2*PI / T
	float w_cos, w_sin;
	float angle;
	int x,u;

	angle_step = (float) (2*PI / data_no);
	for(u = 0; u < data_no; u++){
		out_data[u].real = in_data[0].real;
		out_data[u].image = in_data[0].image;

		for(x = 1; x < data_no; x ++){
			angle = (float) u * x * angle_step;
			w_cos = cos(angle);
			w_sin = sin(angle);
			out_data[u].real = out_data[u].real + w_cos * in_data[x].real - w_sin * in_data[x].image;
			out_data[u].image = out_data[u].image +  w_cos * in_data[x].image + w_sin * in_data[x].real;
		}

		out_data[u].real = out_data[u].real / data_no;
		out_data[u].image = out_data[u].image / data_no;

	}

	return;
}

//旋積
void conv(int data1_M, complex *data1, int data2_N, complex *data2, complex *out_data)
{
	//complex tmpA[ data1_M + 2*(data2_N - 1) ];
	//complex tmpB[ data1_M + 2*(data2_N - 1) ];
	complex *tmpA = (complex *)malloc(sizeof(complex) * (data1_M + 2*(data2_N - 1)));
	complex *tmpB = (complex *)malloc(sizeof(complex) * (data1_M + 2*(data2_N - 1)));
	complex sum;
	complex tmp;
	

	int i,j,k;


	//tmpB水平翻轉，並且補零
	for(j = M + 2*(N - 1) - 1; j >= 0 ; j--){
		if( j < N ){
			tmpB[j].real = data2[N - j - 1].real;
			tmpB[j].image = data2[N - j - 1].image;
		}else{
			tmpB[j].real = 0;
			tmpB[j].image = 0;
		}
	}

	//tmpA上下都補零，中間放數字
	for(j=0; j < (M + 2*(N - 1)); j++){
		if( N-1 <= j && j <= M+N-2 ){
			tmpA[j].real = data1[j-M+1].real;
			tmpA[j].image = data1[j-M+1].image;
		}else{
			tmpA[j].real = 0;
			tmpA[j].image = 0;
		}
	}

	//矩陣相乘
	for(i=0; i < M+N-1; i++){
		sum.real = 0;
		sum.image = 0;
		for(j=0; j < (M + 2*(N - 1)); j++){
			sum.real = sum.real + ( tmpB[j].real * tmpA[j].real - tmpB[j].image * tmpA[j].image );
			sum.image = sum.image + ( tmpB[j].real * tmpA[j].image + tmpB[j].image * tmpA[j].real );
		}
		out_data[i].real = sum.real;
		out_data[i].image = sum.image;

		//做完一次乘法之後，向量向右移1格
		tmp.real = tmpB[ M + 2*(N - 1) -1].real;
		tmp.image = tmpB[ M + 2*(N - 1) -1].image;
		for(k = M + 2*(N - 1) -1; k > 0; k--){
			tmpB[k].real = tmpB[k-1].real;
			tmpB[k].image = tmpB[k-1].image;

		}
		tmpB[0].real = tmp.real;
		tmpB[0].image= tmp.image;

	}

	return;
}

//環形旋積
void cconv(int data_no, complex *data1, complex *data2, complex *out_data)
{
	complex *tmpB = (complex *)malloc(sizeof(complex) * data_no);
	complex sum;
	complex tmp;

	int i,j,k;

	//水平翻轉
	for(j = data_no-1; j >= 0 ; j--){
		tmpB[j].real = data2[data_no - j - 1].real;
		tmpB[j].image = data2[data_no - j - 1].image;
	}

	//然後tmpB向右循環
	tmp.real = tmpB[ data_no -1].real;
	tmp.image = tmpB[ data_no -1].image;
	for(k = data_no -1; k > 0; k--){
		tmpB[k].real = tmpB[k-1].real;
		tmpB[k].image = tmpB[k-1].image;

	}
	tmpB[0].real = tmp.real;
	tmpB[0].image= tmp.image;

	for(i=0; i < data_no; i++){
		out_data[i].real = 0;
		out_data[i].image = 0;
		for(j=0; j < data_no; j++){
			out_data[i].real = out_data[i].real + ( tmpB[j].real * data1[j].real - tmpB[j].image * data1[j].image );
			out_data[i].image = out_data[i].image + ( tmpB[j].real * data1[j].image + tmpB[j].image * data1[j].real ); 
		}

		//tmpB向量平移
		tmp.real = tmpB[ N -1].real;
		tmp.image = tmpB[ N -1].image;
		for(k = N -1; k > 0; k--){
			tmpB[k].real = tmpB[k-1].real;
			tmpB[k].image = tmpB[k-1].image;

		}
		tmpB[0].real = tmp.real;
		tmpB[0].image= tmp.image;
	}

	return;
	
}

//二維傅立葉轉換
void DFT2(IplImage *read_image, complex **out_data)
{
	/**************************************************************
	功能: 二維傅立葉轉換
	輸入: data_no = 輸入資料的列/行資料個數
		  read_image = 輸入影像，型態為IplImage *
	輸出: out_data = 傅立葉轉換結果，型態為complex **

	**************************************************************/
	int i,j,k;
	int sign;
	int data_no = read_image->width;

	

	/****************************************************************
	動態記憶體分配
	read_data = 複數的二維矩陣，用來存取讀進來的影像
	row		  = 複數的一維矩陣，用來暫存取列向量
	col		  = 複數的一維矩陣，用來暫存取行向量
	tmp		  = 複數的一維矩陣，用來暫存每次計算後的結果
	****************************************************************/
	complex **read_data = (complex **)malloc(sizeof(complex *) * data_no);	//取得二維陣列動態記憶體;
	for(i=0; i < data_no; i++){
		read_data[i] = (complex *)malloc(sizeof(complex) * data_no);
	}

	complex *row = (complex *)malloc(sizeof(complex) * data_no);	//取得一維陣列動態記憶體
	complex *col = (complex *)malloc(sizeof(complex) * data_no);	//取得一維陣列動態記憶體
	complex *tmp = (complex *)malloc(sizeof(complex) * data_no);	//取得一維陣列動態記憶體

	
	/****************************************************************
								先平移
	****************************************************************/
	for(i=0; i < data_no; i++){
		for(j=0; j < data_no; j++){
			if( ((i+j) % 2) == 0 )
				sign = 1;
			else
				sign = -1;
			
			read_data[i][j].real = ((float *)(read_image->imageData + i*read_image->widthStep))[j] * sign;
			read_data[i][j].image = 0;

		}
	}
	
	

	/****************************************************************
	複製read_data的列向量 至 row
	對row做傅立葉轉換 至 tmp
	複製tmp 至 out_data的列向量
	****************************************************************/
	for(i=0; i < data_no; i++){
		for(j=0; j < data_no; j++){
			row[j].real = read_data[i][j].real;
			row[j].image = read_data[i][j].image;
		}
		DFT1(data_no, &row[0], &tmp[0]);

		for(j=0; j < data_no; j++){
			out_data[i][j].real = tmp[j].real;
			out_data[i][j].image = tmp[j].image;
		}
	}

	/****************************************************************
	複製out_data的行向量 至 col
	對col做傅立葉轉換 至 tmp
	複製tmp 至 out_data的行向量
	****************************************************************/
	for(j=0; j < data_no; j++){
		for(i=0; i < data_no; i++){
			col[i].real = out_data[i][j].real;
			col[i].image = out_data[i][j].image;
		}
		DFT1(data_no, &col[0], &tmp[0]);

		for(i=0; i < data_no; i++){
			out_data[i][j].real = tmp[i].real;
			out_data[i][j].image = tmp[i].image;
		}
	}

	/******************************************
				釋放記憶體
	******************************************/
	for (i=0; i < data_no; i++) free(read_data[i]);
	free(read_data);

	free(row);
	free(col);
	free(tmp);
	return;
}

//讀取濾波器
void load_filter(FILE *datap, float **out_data)
{
	/*****************************************************************************************************************************
	功能: 讀取濾波器
	輸入: datap = 要讀取的濾波器檔名
	輸出: out_data = 濾波器，型態為float **
	******************************************************************************************************************************/
	int i,j;
	int data_no;

	if(datap==NULL) {
		printf("test file name error!!\n");
		return;
	}

	i=0;
	j=0;
	fscanf(datap, "%d", &data_no );
	for(i=0 ; i < data_no ; i++){
		for(j=0; j < data_no; j++){
			fscanf(datap, "%f", &out_data[i][j] );
		}

	}

	return;
}

//顯示頻譜圖
void fftshow(int data_no, complex **in_data, char *windowName)
{
	/*****************************************************************************************************************************
	功能: 顯示頻譜圖
	輸入: data_no = 輸入資料的列/行資料個數
		  in_data = 要顯示的頻譜圖，型態為complex **
	輸出: 無
	******************************************************************************************************************************/
	IplImage *FreqImage = cvCreateImage(cvSize(data_no, data_no), IPL_DEPTH_8U,	1);

	

	int i,j;
	float max = -999;

	/****************************************************************
	動態記憶體分配
	tmp = 浮點數的二維矩陣，用來儲存複數強度 |F(u,v)|
	****************************************************************/
	float **tmp = (float **)malloc(sizeof(float *) * data_no);	//取得二維陣列動態記憶體;
	for(i=0; i < data_no; i++){
		tmp[i] = (float *)malloc(sizeof(float) * data_no);
	}

	/****************************************************************
	將複數轉換為強度 |F(u,v)| ，同時找到最大值
	****************************************************************/
	for(i = 0; i < data_no; i++){
		for(j=0; j < data_no; j++){
			tmp[i][j] = log(1 + sqrt( in_data[i][j].real * in_data[i][j].real + in_data[i][j].image * in_data[i][j].image ));
			if(tmp[i][j] > max)
				max = tmp[i][j];
		}
	}

	/****************************************************************
	將每個 |F(u,v)| 正規化至0~255，顯示影像
	****************************************************************/
	for(i=0; i < data_no; i++){
		for(j=0; j < data_no; j++){
			((uchar *)(FreqImage->imageData + i*FreqImage->widthStep))[j] = (tmp[i][j] * 255) / max;
		}
	}

	cvNamedWindow(windowName,0);	//創建一個新視窗，命名為"FreqImage"
	cvResizeWindow(windowName, data_no, data_no);	//設定視窗大小
	cvShowImage(windowName, FreqImage);	//顯示圖像
	cvWaitKey(0);

	
	cvDestroyWindow(windowName);	//關閉視窗
	cvReleaseImage(&FreqImage);		//釋放圖像記憶體
	for (i=0; i < data_no; i++) free(tmp[i]);	//釋放二維記憶體
	free(tmp);
    
	return;
}

//二維傅立葉「反」轉換
void IDFT2(complex **read_data, IplImage *out_image)
{
	/*****************************************************************************************************************************
	功能: 二維傅立葉「反」轉換
	輸入: data_no = 輸入資料的列/行資料個數
		  read_data = 濾波後的頻譜圖，型態為complex **
	輸出: out_image = 輸出影像，型態為IplImage *

	******************************************************************************************************************************/
	int i,j,k;
	int sign;
	float max = -999;
	int data_no = out_image->width;

	/****************************************************************
	動態記憶體分配
	out_data  = 複數的二維矩陣，用來儲存反轉換後的數值
	tmpInv	  = 浮點數的二維矩陣，用來儲存正規化數值
	row		  = 複數的一維矩陣，用來暫存取列向量
	col		  = 複數的一維矩陣，用來暫存取行向量
	tmp		  = 複數的一維矩陣，用來暫存每次計算後的結果
	****************************************************************/

	complex **out_data = (complex **)malloc(sizeof(complex *) * data_no);	//取得二維陣列動態記憶體;
	for(i=0; i < data_no; i++){
		out_data[i] = (complex *)malloc(sizeof(complex) * data_no);
	}

	float **tmpInv = (float **)malloc(sizeof(float *) * data_no);	//取得二維陣列動態記憶體;
	for(i=0; i < data_no; i++){
		tmpInv[i] = (float *)malloc(sizeof(float) * data_no);
	}
	
	complex *row = (complex *)malloc(sizeof(complex) * data_no);
	complex *col = (complex *)malloc(sizeof(complex) * data_no);
	complex *tmp = (complex *)malloc(sizeof(complex) * data_no);



	/***************************************
	複製read_data的列向量 至 row
	對row做傅立葉轉換 至 tmp
	複製tmp 至 out_data的列向量
	****************************************/
	for(i=0; i < data_no; i++){
		for(j=0; j < data_no; j++){
			row[j].real = read_data[i][j].real;
			row[j].image = read_data[i][j].image;
		}
		iDFT1(data_no, &row[0], &tmp[0]);

		for(j=0; j < data_no; j++){
			out_data[i][j].real = tmp[j].real;
			out_data[i][j].image = tmp[j].image;
		}
	}

	/****************************************
	複製out_data的行向量 至 col
	對col做傅立葉轉換 至 tmp
	複製tmp 至 out_data的行向量
	****************************************/
	for(j=0; j < data_no; j++){
		for(i=0; i < data_no; i++){
			col[i].real = out_data[i][j].real;
			col[i].image = out_data[i][j].image;
		}
		iDFT1(data_no, &col[0], &tmp[0]);

		for(i=0; i < data_no; i++){
			out_data[i][j].real = tmp[i].real;
			out_data[i][j].image = tmp[i].image;
		}
	}

	/****************************************************************
			將數值轉換為原本在空間域的數值 ，同時找到最大值
	****************************************************************/
	for(i = 0; i < data_no; i++){
		for(j=0; j < data_no; j++){
			
			tmpInv[i][j] = sqrt( out_data[i][j].real * out_data[i][j].real + out_data[i][j].image * out_data[i][j].image );
			if(tmpInv[i][j] > max)
				max = tmpInv[i][j];
		}
	}

	/****************************************************************
				將每個數值正規化至0~255
	****************************************************************/
	for(i=0; i < data_no; i++){
		for(j=0; j < data_no; j++){
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = (tmpInv[i][j] * 255) / max;
		}
	}

	/****************************************************************
						釋放記憶體
	****************************************************************/
	for (i=0; i < data_no; i++) free(tmpInv[i]);
	free(tmpInv);

	for (i=0; i < data_no; i++) free(out_data[i]);
	free(out_data);

	free(row);
	free(col);
	free(tmp);
	
	return;
}

void remove_butter(IplImage *read_image, IplImage *out_image)
{

	int i,j;

	IplImage *read_image64F = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_64F,	1);	//分配記憶體給tmp_image	單通道

	float **lowpass_butterD15N2 = (float **)malloc(sizeof(float *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		lowpass_butterD15N2[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	float **lowpass_butterD40N10 = (float **)malloc(sizeof(float *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		lowpass_butterD40N10[i] = (float *)malloc(sizeof(float) * read_image->width);
	}
	

	complex **read_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		read_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	complex **read_lowpass = (complex **)malloc(sizeof(complex *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		read_lowpass[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	IplImage *noise_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//分配記憶體給tmp_image	單通道
	IplImage *noise_image64F = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_64F,	1);	//分配記憶體給tmp_image	單通道

	complex **noise_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		noise_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	complex **noise_lowpass = (complex **)malloc(sizeof(complex *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		noise_lowpass[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	FILE *lowpass_butterD15N2_ptr = fopen("butter_lowpass_D15N2_256.txt","r");
	FILE *lowpass_butterD40N10_ptr = fopen("butter_lowpass_D100N10_256.txt","r");


	load_filter(lowpass_butterD15N2_ptr, lowpass_butterD15N2);	//讀取濾波器
	load_filter(lowpass_butterD40N10_ptr, lowpass_butterD40N10);	//讀取濾波器


	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((float *)(read_image64F->imageData + (i)*read_image64F->widthStep))[j] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
		}
	}
	DFT2(read_image64F, read_DFT);	//二維傅立葉轉換
	fftshow(256, read_DFT, "read_DFT");	//顯示頻譜圖

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			read_lowpass[i][j].real = read_DFT[i][j].real * lowpass_butterD15N2[i][j];
			read_lowpass[i][j].image = read_DFT[i][j].image * lowpass_butterD15N2[i][j];
		}
	}
	//fftshow(256, read_lowpass, "read_lowpass");	//顯示頻譜圖
	
	IDFT2(read_lowpass, noise_image);	//二維傅立葉「反」轉換
	imshow(noise_image, "noise_image");	//顯示影像


	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((float *)(noise_image64F->imageData + (i)*noise_image64F->widthStep))[j] = ((uchar *)(noise_image->imageData + (i)*noise_image->widthStep))[j];
		}
	}
	DFT2(noise_image64F, noise_DFT);	//二維傅立葉轉換
	fftshow(256, noise_DFT, "noise_DFT");	//顯示頻譜圖


	/********************方法2 除式設限********************************/
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			if(lowpass_butterD15N2[i][j] < 0.001)
				lowpass_butterD15N2[i][j] = 1;
		}
	}
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			noise_lowpass[i][j].real = noise_DFT[i][j].real / lowpass_butterD15N2[i][j];
			noise_lowpass[i][j].image = noise_DFT[i][j].image / lowpass_butterD15N2[i][j];
		}
	}
	//fftshow(256, noise_lowpass, "noise_lowpass");	//顯示頻譜圖
	

	IDFT2(noise_lowpass, out_image);	//二維傅立葉「反」轉換


	for (i=0; i < read_image->width; i++) free(lowpass_butterD15N2[i]);
	free(lowpass_butterD15N2);

	for (i=0; i < read_image->width; i++) free(lowpass_butterD40N10[i]);
	free(lowpass_butterD40N10);

	for (i=0; i < read_image->width; i++) free(read_DFT[i]);
	free(read_DFT);

	for (i=0; i < read_image->width; i++) free(read_lowpass[i]);
	free(read_lowpass);

	cvReleaseImage(&noise_image);		//釋放圖像記憶體

	for (i=0; i < read_image->width; i++) free(noise_DFT[i]);
	free(noise_DFT);

	for (i=0; i < read_image->width; i++) free(noise_lowpass[i]);
	free(noise_lowpass);

	fclose(lowpass_butterD15N2_ptr);
	fclose(lowpass_butterD40N10_ptr);


	return;
}

void remove_motion(IplImage *read_image, IplImage *out_image)
{
	int i,j,tmpi,k;
	float a,b,c,d;
	//float motionMask[7] = {0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429};
	float tmp[7] = {0};
	float sum;
	float max = -999;

	IplImage *noise_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//分配記憶體給tmp_image	單通道
	IplImage *noise_image64F = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_64F,	1);	//分配記憶體給tmp_image	單通道
	IplImage *motion_filter64F = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_64F,	1);	//分配記憶體給tmp_image	單通道

	complex **motion_filter_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		motion_filter_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	complex **noise_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		noise_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	complex **recover_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		recover_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	/*complex **noise_lowpass = (complex **)malloc(sizeof(complex *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		noise_lowpass[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}*/

	
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			
			for(k=-3; k < 4; k++){
				if( (j+k) < 0 )
					tmp[k+3] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[0];
				else if( (j+k) >= read_image->width )
					tmp[k+3] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[read_image->width - 1];
				else
					tmp[k+3] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j+k];
			}

			sum = 0;
			for(k=0; k < 7; k++){
				sum = sum + tmp[k];
			}

			((uchar *)(noise_image->imageData + i*noise_image->widthStep))[j] = sum / 7;
		}
	}

	//imshow(noise_image, "noise_image");	//顯示影像

	for(i=0; i < motion_filter64F->height; i++){
		for(j=0; j < motion_filter64F->width; j++){
			if(i==0 & j < 7)
				((float *)(motion_filter64F->imageData + i*motion_filter64F->widthStep))[j] = 0.1429;
			else
				((float *)(motion_filter64F->imageData + i*motion_filter64F->widthStep))[j] = 0.0;
		}
	}

	DFT2(motion_filter64F, motion_filter_DFT);	//二維傅立葉轉換
	//fftshow(256, motion_filter_DFT, "motion_filter_DFT");	//顯示頻譜圖
	
	U8toF64(noise_image, noise_image64F);	//8U影像轉64F
	DFT2(noise_image64F, noise_DFT);	//二維傅立葉轉換
	//fftshow(256, noise_DFT, "noise_DFT");	//顯示頻譜圖


	
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			a = noise_DFT[i][j].real;
			b = noise_DFT[i][j].image;
			c = motion_filter_DFT[i][j].real;
			d = motion_filter_DFT[i][j].image;

			if( sqrt(c*c + d*d) < 0.02 ){
				recover_DFT[i][j].real = noise_DFT[i][j].real;
				recover_DFT[i][j].image = noise_DFT[i][j].image;
			}else{
				recover_DFT[i][j].real = ((a*c) + (b*d)) / (c*c + d*d);
				recover_DFT[i][j].image = ((-1)*(a*d) + (b*c)) / (c*c + d*d);
			}
			
			
		}
	}
	//fftshow(256, recover_DFT, "recover_DFT");	//顯示頻譜圖

	IDFT2(recover_DFT, out_image);	//二維傅立葉「反」轉換

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = ((uchar *)(out_image->imageData + i*out_image->widthStep))[j] * 3;
		}

	}


	cvReleaseImage(&noise_image);		//釋放圖像記憶體
	cvReleaseImage(&noise_image64F);		//釋放圖像記憶體
	cvReleaseImage(&motion_filter64F);		//釋放圖像記憶體

	for (i=0; i < read_image->width; i++) free(motion_filter_DFT[i]);
	free(motion_filter_DFT);

	for (i=0; i < read_image->width; i++) free(noise_DFT[i]);
	free(noise_DFT);

	for (i=0; i < read_image->width; i++) free(recover_DFT[i]);
	free(recover_DFT);


	return;

}

void wiener_filter(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	float K = 0.00001;
	float A, sumNoise, avgNoise, varNoise;
	float a,b,c,d;

	IplImage *read_image64F = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_64F,	1);	//分配記憶體給tmp_image	單通道
	complex **read_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		read_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	float **lowpass_butterD15N2 = (float **)malloc(sizeof(float *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		lowpass_butterD15N2[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	complex **noise_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		noise_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	IplImage *noise_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//分配記憶體給tmp_image	單通道

	IplImage *noise_image64F = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_64F,	1);	//分配記憶體給tmp_image	單通道
	complex **recover_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->width; i++){
		recover_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	FILE *lowpass_butterD15N2_ptr = fopen("butter_lowpass_D15N2_256.txt","r");

	


	U8toF64(read_image, read_image64F);	//8U影像轉64F
	DFT2(read_image64F, read_DFT);	//二維傅立葉轉換

	load_filter(lowpass_butterD15N2_ptr, lowpass_butterD15N2);	//讀取濾波器

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			noise_DFT[i][j].real = read_DFT[i][j].real * lowpass_butterD15N2[i][j];
			noise_DFT[i][j].image = read_DFT[i][j].image * lowpass_butterD15N2[i][j];
		}	
	}

	IDFT2(noise_DFT, noise_image);	//二維傅立葉「反」轉換

	/*****************************************************************/

	U8toF64(noise_image, noise_image64F);	//8U影像轉64F
	DFT2(noise_image64F, noise_DFT);	//二維傅立葉轉換

	sumNoise = 0;
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			sumNoise = sumNoise + ((uchar *)(noise_image->imageData + (i)*noise_image->widthStep))[j];
		}

	}
	avgNoise = sumNoise / (read_image->height * read_image->width);

	sumNoise = 0;
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			A = ( ((uchar *)(noise_image->imageData + (i)*noise_image->widthStep))[j] - avgNoise ) * ( ((uchar *)(noise_image->imageData + (i)*noise_image->widthStep))[j] - avgNoise );
			sumNoise = sumNoise + A;
		}

	}
	varNoise = sumNoise / (read_image->height * read_image->width - 1);

	
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			A = (lowpass_butterD15N2[i][j] * lowpass_butterD15N2[i][j]) / (lowpass_butterD15N2[i][j] * lowpass_butterD15N2[i][j] + K);
			recover_DFT[i][j].real = (A / lowpass_butterD15N2[i][j]) * noise_DFT[i][j].real;
			recover_DFT[i][j].image = (A / lowpass_butterD15N2[i][j]) * noise_DFT[i][j].image;
		}

	}

	IDFT2(recover_DFT, out_image);	//二維傅立葉「反」轉換

	return;


}

int Otsu(IplImage *read_image)
{
	int i,j,k;
	int Number;	//總pixel數
	int sum_of_number;	//用來統計累積次數
	int start,end;
	float min;
	float sum;
	int outValue;	//輸出值，也就是訓練好的閥值
	

	int *number = (int *)malloc(sizeof(int) * 256);

	float *probability = (float *)malloc(sizeof(float) * 256);
	float *accu_probabilityFore = (float *)malloc(sizeof(float) * 256);
	float *accu_probabilityBack = (float *)malloc(sizeof(float) * 256);

	float mean;
	float *meanFore = (float *)malloc(sizeof(float) * 256);
	float *meanBack = (float *)malloc(sizeof(float) * 256);

	float *varFore = (float *)malloc(sizeof(float) * 256);
	float *varBack = (float *)malloc(sizeof(float) * 256);

	float *sum_of_var = (float *)malloc(sizeof(float) * 256);

	//float *grayTH = (float *)malloc(sizeof(float) * 256);



	Number = read_image->height * read_image->width;

	//初始化
	for(i=0; i < 256; i++){
		number[i] = 0;
		probability[i] = 0;
		accu_probabilityFore[i] = 0;
		accu_probabilityBack[i] = 0;
		meanFore[i] = 0;
		meanBack[i] = 0;
		varFore[i] = 0;
		varBack[i] = 0;
	}

	//灰階值統計
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			number[ ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] ] ++;
		}
	}

	//找出開始點
	for(i=0; i < 256; i++){
		if(number[i]){
			start = i;
			break;
		}
	}
	//找出結束點
	for(i=255; i > 0; i--){
		if(number[i]){
			end = i;
			break;
		}
	}

	
	//計算灰階值機率
	for(i = start; i < end; i++){
		probability[i] = (float)number[i] / Number;
	}
	
	
	//計算各群集累積機率
	mean = 0;
	for(i = start; i <= end; i++){
		if(i == 0){
			accu_probabilityFore[i] = probability[i];
			accu_probabilityBack[i] = 1 - probability[i];
		}else{
			accu_probabilityFore[i] = accu_probabilityFore[i-1] + probability[i];	//前景累積機率
			accu_probabilityBack[i] = 1 - accu_probabilityFore[i];					//背景累積機率
		}
		mean = mean + (i * probability[i]);	//總平均值(期望值)
	}
	

	//計算各群集平均值
	for(i = start; i < end; i++){
		sum = 0;
		sum_of_number = 0;
		for(j = start; j <= i; j++){
			sum = sum + (j * number[j]);					//灰階值 * 該灰階值統計次數
			sum_of_number = sum_of_number + number[j];		//前景總累積次數
		}
		meanFore[i] = sum / sum_of_number;	//前景平均值

		sum = 0;
		sum_of_number = 0;
		for(k = j; k < end; k++){
			sum = sum + (k * number[k]);					//灰階值 * 該接值統計次數
			sum_of_number = sum_of_number + number[k];		//背景總累積次數
		}
		meanBack[i] = sum / sum_of_number;	//背景平均值
		
	}

	//計算各群集變異數
	for(i = start; i <= end; i++){
		sum = 0;
		sum_of_number = 0;
		for(j = start; j <= i; j++){
			sum = sum + (j - meanFore[i]) * (j - meanFore[i]) * number[j];	//(灰階值 - 前景平均值).^2 * 該灰階值統計次數
			sum_of_number = sum_of_number + number[j];						//前景總累積次數
		}
		varFore[i] = sum / sum_of_number;	//前景變異數

		sum = 0;
		sum_of_number = 0;
		for(k = j; k <= end; k++){
			sum = sum + (k - meanBack[i]) * (k - meanBack[i]) * number[k];	//(灰階值 - 背景平均值).^2 * 該灰階值統計次數
			sum_of_number = sum_of_number + number[k];						//背景總累積次數
		}
		varBack[i] = sum / sum_of_number;	//背景變異數
		
	}
	

	//計算群集變異的加權總合
	for(i = start; i < end; i++){
		sum_of_var[i] = accu_probabilityFore[i] * varFore[i] +  accu_probabilityBack[i] * varBack[i];	//前景累積機率 * 前景變異數 + 背景累積機率 * 背景變異數
	}
	
	//找出變異加權總合的最小值
	min = 99999999;
	for(i = start; i < end; i++){
		if(sum_of_var[i] < min){
			min = sum_of_var[i];
			outValue = i;	//i 就是我們要的閥值
		}
	}



	free(number);

	free(probability);
	free(accu_probabilityFore);
	free(accu_probabilityBack);

	free(meanFore);
	free(meanBack);

	free(varFore);
	free(varBack);
	free(sum_of_var);

	return outValue;

}


void fprintf_hist(IplImage *read_image, char *fileName)
{
	FILE *fp = fopen(fileName,"w");

	int i,j;
	int a1,a2,a3,a4,a5;
	int A[5];
	int sum;
	int time;	//欲平滑次數

	int *hist = (int *)malloc(sizeof(int) * 256);		//直條圖統計原始資料
	int *hist_in = (int *)malloc(sizeof(int) * 256);	//直條圖平滑前
	int *hist_out = (int *)malloc(sizeof(int) * 256);	//直條圖平滑後


	//初始化
	for(i=0; i < 256; i++){
		hist[i] = 0;
		hist_in[i] = 0;
		hist_out[i] = 0;
	}

	//灰階值統計
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			hist[ ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] ] ++;
		}
	}

	for(i = 0; i < 256; i++){
		hist_in[i] = hist[i];
		hist_out[i] = hist[i];
	}


	/*********請先設定平滑次數*************/
	time = 5;
	for(j = 0; j < time; j++){

		for(i = 0; i < 256; i++){
			if( (i - 2) < 0)	a1 = 0;
			else	a1 = i-2;

			if( (i - 1) < 0 )	a2 = 0;
			else a2 = i-1;

			a3 = i;

			if( (i + 1) >= 256 )	a4 = 255;
			else a4 = i+1;

			if( (i + 2) >= 256 )	a5 = 255;
			else a5 = i+1;

			sum = hist_in[a1] + hist_in[a2] + hist_in[a3] + hist_in[a4] + hist_in[a5];
			hist_out[i] = sum / 5;

		}

		for(i = 0; i < 256; i++){
			hist_in[i] = hist_out[i];
		}
	}


	for(i=0; i < 256; i++){
		fprintf(fp, "%d	", hist_out[i]);

	}

	fclose(fp);
	free(hist);
	free(hist_in);
	free(hist_out);

	return;
}


void difference_filter(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	int up, down, left, right;
	int X[9],Y[9];
	int sumx,sumy;
	float z;

	int cx_difference[3][3] = {0,	0,	0,
							   0,	1,	-1,
							   0,	0,	0};

	int cy_difference[3][3] = {0,	0,	0,
							   0,	1,	0,
							   0,	-1,	0};

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
			

			X[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * cx_difference[0][0];
			X[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * cx_difference[0][1];
			X[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * cx_difference[0][2];

			X[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * cx_difference[1][0];
			X[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * cx_difference[1][1];
			X[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * cx_difference[1][2];

			X[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * cx_difference[2][0];
			X[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * cx_difference[2][1];
			X[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * cx_difference[2][2];



			Y[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * cy_difference[0][0];
			Y[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * cy_difference[0][1];
			Y[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * cy_difference[0][2];

			Y[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * cy_difference[1][0];
			Y[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * cy_difference[1][1];
			Y[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * cy_difference[1][2];

			Y[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * cy_difference[2][0];
			Y[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * cy_difference[2][1];
			Y[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * cy_difference[2][2];

			sumx = X[0] + X[1] + X[2] + X[3] + X[4] + X[5] + X[6] + X[7] + X[8];
			sumy = Y[0] + Y[1] + Y[2] + Y[3] + Y[4] + Y[5] + Y[6] + Y[7] + Y[8];

			z = sqrt((float)sumx*sumx + (float)sumy*sumy);

			
			
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = is_overflow((int)z);
			
		}
	}

	return;


}

void roberts_filter(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	int up, down, left, right;
	int X[9],Y[9];
	int sumx,sumy;
	float z;

	int cx_sobel[3][3] = {0,	0,	0,
						  0,	-1,	0,
						  0,	0,	1};

	int cy_sobel[3][3] = {0,	0,	0,
						  0,	0,	-1,
						  0,	1,	0};

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){

			up = i-1;
			down = i+1;
			left = j-1;
			right = j+1;

			if( up < 0 )	up = 0;	
			if( down >= read_image->height )	down = read_image->height - 1;
			if( left < 0 ) left = 0;
			if( right >= read_image->width ) right = read_image->width - 1;
			

			X[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * cx_sobel[0][0];
			X[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * cx_sobel[0][1];
			X[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * cx_sobel[0][2];

			X[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * cx_sobel[1][0];
			X[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * cx_sobel[1][1];
			X[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * cx_sobel[1][2];

			X[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * cx_sobel[2][0];
			X[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * cx_sobel[2][1];
			X[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * cx_sobel[2][2];



			Y[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * cy_sobel[0][0];
			Y[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * cy_sobel[0][1];
			Y[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * cy_sobel[0][2];

			Y[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * cy_sobel[1][0];
			Y[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * cy_sobel[1][1];
			Y[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * cy_sobel[1][2];

			Y[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * cy_sobel[2][0];
			Y[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * cy_sobel[2][1];
			Y[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * cy_sobel[2][2];

			sumx = X[0] + X[1] + X[2] + X[3] + X[4] + X[5] + X[6] + X[7] + X[8];
			sumy = Y[0] + Y[1] + Y[2] + Y[3] + Y[4] + Y[5] + Y[6] + Y[7] + Y[8];

			z = sqrt((float)sumx*sumx + (float)sumy*sumy);


			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = is_overflow((int)z);
			
		}
	}

	return;


}

void sobel_filter(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	int up, down, left, right;
	int X[9],Y[9];
	int sumx,sumy;
	float z;

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

			X[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * cx_sobel[0][0];
			X[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * cx_sobel[0][1];
			X[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * cx_sobel[0][2];

			X[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * cx_sobel[1][0];
			X[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * cx_sobel[1][1];
			X[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * cx_sobel[1][2];

			X[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * cx_sobel[2][0];
			X[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * cx_sobel[2][1];
			X[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * cx_sobel[2][2];



			Y[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * cy_sobel[0][0];
			Y[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * cy_sobel[0][1];
			Y[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * cy_sobel[0][2];

			Y[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * cy_sobel[1][0];
			Y[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * cy_sobel[1][1];
			Y[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * cy_sobel[1][2];

			Y[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * cy_sobel[2][0];
			Y[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * cy_sobel[2][1];
			Y[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * cy_sobel[2][2];

			sumx = X[0] + X[1] + X[2] + X[3] + X[4] + X[5] + X[6] + X[7] + X[8];
			sumy = Y[0] + Y[1] + Y[2] + Y[3] + Y[4] + Y[5] + Y[6] + Y[7] + Y[8];

			z = sqrt((float)sumx*sumx + (float)sumy*sumy);


			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = is_overflow((int)z);
			
		}
	}

	return;


}

void prewitt(IplImage *read_image, IplImage *out_image)
{
	int i,j,k;
	int up, down, left, right;
	int d0, d1, d2, d3, d4, d5, d6, d7, d8;
	int A[8];
	int max;

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

			d0 = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left];
			d1 = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j];
			d2 = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right];

			d3 = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left];
			d4 = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
			d5 = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right];

			d6 = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left];
			d7 = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j];
			d8 = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right];

			A[0] = d0 + d1 + d2 + d3 -2*d4 + d5 - d6 - d7 - d8;
			A[1] = d0 + d1 + d2 + d3 -2*d4 - d5 + d6 - d7 - d8;
			A[2] = d0 + d1 - d2 + d3 -2*d4 - d5 + d6 + d7 - d8;

			A[3] = d0 - d1 - d2 + d3 -2*d4 - d5 + d6 + d7 + d8;
			A[4] = -d0 - d1 - d2 + d3 -2*d4 + d5 + d6 + d7 + d8;
			A[5] = -d0 - d1 + d2 - d3 -2*d4 + d5 + d6 + d7 + d8;

			A[6] = -d0 + d1 + d2 - d3 -2*d4 + d5 - d6 + d7 + d8;
			A[7] = d0 + d1 + d2 - d3 -2*d4 + d5 - d6 - d7 + d8;

			max = -99999;
			for(k = 0; k < 8; k++){
				if( max < A[k] )
					max = A[k];
			}

			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = is_overflow(max);
		}

	}
	return;
}

int conj_count(int p[9])
{
	int i, i1, i2;
	int q[9];
	int n = 0;

	for(i = 0; i < 9; i++){
		if( p[i] == 1 || p[i] == -1 )
			q[i] = 0;
		else
			q[i] = 1;
	}

	for(i = 1; i < 9; i+=2){
		i1 = i + 1;
		i2 = i + 2;

		if(i2 == 9)
			i2 = 1;

		n = n + q[i] - q[i] * q[i1] * q[i2];
	}
	return n;

}

void thinning(IplImage *read_image, IplImage *out_image)	//二值影像細線化
{
	int flag = 1;
	int i, j, k, n;
	int up, down, left, right;
	int p[9];	//圖形:1 背景:2 背景後補:-1;
	int TMP = 128;	//背景候補的濃度暫定值

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
		}
	}

	while( flag != 0){
		flag = 0;
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

				p[0] = ((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j];
				p[1] = ((uchar *)(out_image->imageData + (i)*out_image->widthStep))[right];
				p[2] = ((uchar *)(out_image->imageData + (up)*out_image->widthStep))[right];

				p[3] = ((uchar *)(out_image->imageData + (up)*out_image->widthStep))[j];
				p[4] = ((uchar *)(out_image->imageData + (up)*out_image->widthStep))[left];
				p[5] = ((uchar *)(out_image->imageData + (i)*out_image->widthStep))[left];

				p[6] = ((uchar *)(out_image->imageData + (down)*out_image->widthStep))[left];
				p[7] = ((uchar *)(out_image->imageData + (down)*out_image->widthStep))[j];
				p[8] = ((uchar *)(out_image->imageData + (down)*out_image->widthStep))[right];

				for(k = 0; k < 9; k++){
					if(p[k] == 255)	p[k] = 1;
					else if(p[k] == 0)	p[k] = 0;
					else	p[k] = -1;
				}

			

				/* 條件1:  */
				if(p[0] != 1)	continue;

				if(p[1] * p[3] * p[5] * p[7] != 0)	continue;

				n = 0;
				for(k = 1; k < 9; k++){
					if(p[k] != 0)
						n++;
				}
				if(n < 2)	continue;

				n = 0;
				for(k = 1; k < 9; k++){
					if(p[k] == 1)
						n++;
				}
				if(n < 1)	continue;

				if( conj_count(p) != 1 )	continue;

				n=0;
				for(k = 1; k < 9; k++){
					if(p[k] != -1)
						n++;
					else if(p[k] == -1){
						p[k] = 0;
						
						if( conj_count(p) == 1 )
							n++;

						p[k] = -1;
					}
				}
				if(n < 8)	continue;

				((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = TMP;
				flag++;
				
			}
		}

		for(i = 0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				if(((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] == TMP)
					((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 0;
			}
		}
	}

	
	return;

}


void laplacian_filater(IplImage *read_image, IplImage *out_image)
{
	
	int i,j;
	int up, down, left, right;
	int A[9];
	int sum;
	

	int lap[3][3] = {0, -1,	0,
					 -1, 4,	-1,
					 0,	-1,	0};

	int **tmp_matrix = (int **)malloc(sizeof(int *) * read_image->height);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->height; i++){
		tmp_matrix[i] = (int *)malloc(sizeof(int) * read_image->width);
	}


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

			A[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * lap[0][0];
			A[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * lap[0][1];
			A[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * lap[0][2];

			A[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * lap[1][0];
			A[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * lap[1][1];
			A[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * lap[1][2];

			A[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * lap[2][0];
			A[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * lap[2][1];
			A[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * lap[2][2];

			sum = A[0] + A[1] + A[2] + A[3] + A[4] + A[5] + A[6] + A[7] + A[8];

			tmp_matrix[i][j] = sum;

		}

	}

	mat2grat(tmp_matrix, out_image);	//矩陣資料 轉 灰階影像

	return;


}

void zero_cross(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	int up, down, left, right;
	int A[9];
	int sum;
	int TH;
	

	int lap[3][3] = {0, -1,	0,
					 -1, 4,	-1,
					 0,	-1,	0};

	int **tmp_matrix = (int **)malloc(sizeof(int *) * read_image->height);	//取得二維陣列動態記憶體;
	for(i=0; i < read_image->height; i++){
		tmp_matrix[i] = (int *)malloc(sizeof(int) * read_image->width);
	}


	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 0;
		}
	}



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

			A[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left] * lap[0][0];
			A[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j] * lap[0][1];
			A[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right] * lap[0][2];

			A[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left] * lap[1][0];
			A[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] * lap[1][1];
			A[5] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right] * lap[1][2];

			A[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left] * lap[2][0];
			A[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j] * lap[2][1];
			A[8] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right] * lap[2][2];

			sum = A[0] + A[1] + A[2] + A[3] + A[4] + A[5] + A[6] + A[7] + A[8];

			tmp_matrix[i][j] = sum;

		}
	}

	//print_int_mat_txt(tmp_matrix, read_image->height, read_image->width, "lap.txt");	//印出整數matrix

	TH = 150;
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

			if(tmp_matrix[i][j] == 0){
				if( (tmp_matrix[up][j] * tmp_matrix[down][j]) < 0 && abs((tmp_matrix[up][j] - tmp_matrix[down][j])) > TH )	((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
				if( (tmp_matrix[i][left] * tmp_matrix[i][right]) < 0 && abs(tmp_matrix[i][left] - tmp_matrix[i][right]) > TH)	((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
				if( (tmp_matrix[up][left] * tmp_matrix[down][right]) < 0 && abs(tmp_matrix[up][left] * tmp_matrix[down][right]) > TH)	((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
				if( (tmp_matrix[down][left] * tmp_matrix[up][right]) < 0 && abs(tmp_matrix[down][left] * tmp_matrix[up][right]) > TH)	((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
			}else{
				if( (tmp_matrix[up][left] * tmp_matrix[i][j]) < 0 && abs(tmp_matrix[up][left] * tmp_matrix[i][j]) > TH)	((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
				if( (tmp_matrix[up][j] * tmp_matrix[i][j]) < 0 && abs(tmp_matrix[up][j] * tmp_matrix[i][j]) > TH)	((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
				if( (tmp_matrix[up][right] * tmp_matrix[i][j]) < 0 && abs(tmp_matrix[up][right] * tmp_matrix[i][j]) > TH)	((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
				if( (tmp_matrix[i][left] * tmp_matrix[i][j]) < 0 && abs(tmp_matrix[i][left] * tmp_matrix[i][j]) > TH)	((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
			}

		}
	}

	

	return;



}

void dilation(IplImage *read_image, IplImage *out_image)	//膨脹
{
	int i,j;
	int up, down, left, right;

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 0;
		}
	}

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

			if(((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] == 255){

				((uchar *)(out_image->imageData + (up)*out_image->widthStep))[left] = 255;
				((uchar *)(out_image->imageData + (up)*out_image->widthStep))[j] = 255;
				((uchar *)(out_image->imageData + (up)*out_image->widthStep))[right] = 255;

				((uchar *)(out_image->imageData + (i)*out_image->widthStep))[left] = 255;
				((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
				((uchar *)(out_image->imageData + (i)*out_image->widthStep))[right] = 255;

				((uchar *)(out_image->imageData + (down)*out_image->widthStep))[left] = 255;
				((uchar *)(out_image->imageData + (down)*out_image->widthStep))[j] = 255;
				((uchar *)(out_image->imageData + (down)*out_image->widthStep))[right] = 255;


				
	
			}
		}
	}

	return;

}

void erosion(IplImage *read_image, IplImage *out_image)	//侵蝕
{
	int i, j, k;
	int up, down, left, right;
	int A[8];
	int min;

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 0;
		}
	}

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

			if(((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] == 255){
				
				A[0] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[left];
				A[1] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[j];
				A[2] = ((uchar *)(read_image->imageData + (up)*read_image->widthStep))[right];

				A[3] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[left];
				A[4] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[right];

				A[5] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[left];
				A[6] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[j];
				A[7] = ((uchar *)(read_image->imageData + (down)*read_image->widthStep))[right];

				if(A[0] == 255 && A[1] == 255 && A[2] == 255 && A[3] == 255 && A[4] == 255 && A[5] == 255 && A[6] == 255 && A[7] == 255)
					((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
				else
					((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 0;

				//imshow(out_image, "out");	//顯示影像
	
			}
		}
	}

	return;

}

void imopen(IplImage *read_image, IplImage *out_image)	//開運算
{
	int i, j, k;
	int up, down, left, right;
	int A[8];
	int min;

	IplImage *erosion_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 0;
		}
	}

	erosion(read_image, erosion_image);	//侵蝕
	dilation(erosion_image, out_image);	//膨脹

	cvReleaseImage(&erosion_image );	//釋放圖像記憶體

	return;

}

void regfill(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	int flag;

	IplImage *current_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//當前圖
	IplImage *last_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);		//過去圖
	IplImage *edge_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);		//邊界圖
	IplImage *invedge_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//邊界反向圖
	IplImage *erosion_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//侵蝕圖
	IplImage *dilation_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//膨脹圖

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(current_image->imageData + (i)*current_image->widthStep))[j] = 0;
			((uchar *)(last_image->imageData + (i)*last_image->widthStep))[j] = 0;
			((uchar *)(edge_image->imageData + (i)*edge_image->widthStep))[j] = 0;
			((uchar *)(invedge_image->imageData + (i)*invedge_image->widthStep))[j] = 0;
			((uchar *)(erosion_image->imageData + (i)*erosion_image->widthStep))[j] = 0;
			((uchar *)(dilation_image->imageData + (i)*dilation_image->widthStep))[j] = 0;
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 0;
		}
	}


	erosion(read_image, erosion_image);	//侵蝕
	
	//求邊界圖
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(edge_image->imageData + (i)*edge_image->widthStep))[j] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] - ((uchar *)(erosion_image->imageData + (i)*erosion_image->widthStep))[j];
		}
	}
	
	
	//求邊界反向圖
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){

			if(((uchar *)(edge_image->imageData + (i)*edge_image->widthStep))[j] == 255)
				((uchar *)(invedge_image->imageData + (i)*invedge_image->widthStep))[j] = 0;
			else
				((uchar *)(invedge_image->imageData + (i)*invedge_image->widthStep))[j] = 255;
		}
	}
	

	flag = 1;
	//決定要填充的區域
	((uchar *)(last_image->imageData + (26)*last_image->widthStep))[20] = 255;
	((uchar *)(last_image->imageData + (39)*last_image->widthStep))[93] = 255;
	((uchar *)(last_image->imageData + (166)*last_image->widthStep))[165] = 255;
	
	while(flag != 0){
		flag = 0;
		dilation(last_image, dilation_image);
		
		//膨脹圖 和 邊界反向圖 取交集，結果給「當前圖」
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				
				if(((uchar *)(dilation_image->imageData + (i)*dilation_image->widthStep))[j] == 255 && ((uchar *)(invedge_image->imageData + (i)*invedge_image->widthStep))[j] == 255)
					((uchar *)(current_image->imageData + (i)*current_image->widthStep))[j] = 255;
				else
					((uchar *)(current_image->imageData + (i)*current_image->widthStep))[j] = 0;
			}
		}

		

		//比對當前圖 和 過去圖 是否一樣
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				if(((uchar *)(current_image->imageData + (i)*current_image->widthStep))[j] == ((uchar *)(last_image->imageData + (i)*last_image->widthStep))[j])
					flag++;
			}
		}

		//把當前圖copy給過去圖
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				((uchar *)(last_image->imageData + (i)*last_image->widthStep))[j] = ((uchar *)(current_image->imageData + (i)*current_image->widthStep))[j];
			}
		}


		//如果flag 等於 所有pixel數，代表已經填充完成，將flag設為0，跳出迴圈
		if(flag == (read_image->height * read_image->width))
			flag = 0;
	}

	//將當前圖輸出
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			if(((uchar *)(current_image->imageData + (i)*current_image->widthStep))[j] == 255 || ((uchar *)(edge_image->imageData + (i)*edge_image->widthStep))[j] == 255)
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
		}
	}

	cvReleaseImage(&current_image );	//釋放圖像記憶體
	cvReleaseImage(&last_image );	//釋放圖像記憶體
	cvReleaseImage(&edge_image );	//釋放圖像記憶體
	cvReleaseImage(&invedge_image );	//釋放圖像記憶體
	cvReleaseImage(&erosion_image );	//釋放圖像記憶體
	cvReleaseImage(&dilation_image );	//釋放圖像記憶體

	return;
}


void skeleton(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	int flag;


	IplImage *wait_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//待處理圖
	IplImage *erosion_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//侵蝕圖
	IplImage *open_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//開運算圖
	IplImage *diffset_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//集合差圖

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(wait_image->imageData + (i)*wait_image->widthStep))[j] = 0;
			((uchar *)(erosion_image->imageData + (i)*erosion_image->widthStep))[j] = 0;
			((uchar *)(open_image->imageData + (i)*open_image->widthStep))[j] = 0;
			((uchar *)(diffset_image->imageData + (i)*diffset_image->widthStep))[j] = 0;
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 0;
		}
	}

	//原圖copy給待處理圖
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(wait_image->imageData + (i)*wait_image->widthStep))[j] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
		}
	}

	flag = 1;
	while(flag != 0){
		flag = 0;

		erosion(wait_image, erosion_image);	//侵蝕
		imopen(erosion_image, open_image);	//開運算

		//集合差圖 = 侵蝕圖 - 開運算圖
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				((uchar *)(diffset_image->imageData + (i)*diffset_image->widthStep))[j] = ((uchar *)(erosion_image->imageData + (i)*erosion_image->widthStep))[j] - ((uchar *)(open_image->imageData + (i)*open_image->widthStep))[j];
			}
		}

		//輸出圖 = 輸出圖 聯集 集合差圖
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				if(((uchar *)(diffset_image->imageData + (i)*diffset_image->widthStep))[j] == 255 || ((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] == 255)
					((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
			}
		}

		//檢驗開運算是否為空集合。flag != 0 : 不是空集合 ， flag == 0 : 為空集合
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				if(((uchar *)(open_image->imageData + (i)*open_image->widthStep))[j] == 255)
					flag++;
			}
		}

		//侵蝕圖copy給待處理圖
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				((uchar *)(wait_image->imageData + (i)*wait_image->widthStep))[j] = ((uchar *)(erosion_image->imageData + (i)*erosion_image->widthStep))[j];
			}
		}
	}


	cvReleaseImage(&wait_image );	//釋放圖像記憶體
	cvReleaseImage(&erosion_image );	//釋放圖像記憶體
	cvReleaseImage(&open_image );	//釋放圖像記憶體
	cvReleaseImage(&diffset_image );	//釋放圖像記憶體

	return;

}

void label(IplImage *read_image, IplImage *out_image, int part)
{
	int i,j;
	int label[99999] = {0};
	int labelCnt = 50;
	int labelMin;
	int up, left, rightUp, leftUp;

	int objCnt[99999] = {0};
	int area[9999] = {0};
	int sum = 0;

	IplImage *label_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//標記圖
	
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = 0;

			/*if(i == 0 || j == 0 || i == read_image->height - 1 || j == read_image->width - 1){
				((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] = 0;
				continue;
			}
			if(((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] == 0)
				((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] = 255;
			else
				((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] = 0;*/
		}
	}

	
	if(part == 4){
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				if(((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] == 255){

					if(i == 0 && j == 0){	//第1個pixel，直接給新標記
						((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
						label[labelCnt] = labelCnt;
						labelCnt ++;
					}else if(i == 0 && j != 0){	//第0列，只檢查左邊狀態
						left = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j-1];
						if(left == 255){	//左邊為前景，當格和左邊一樣
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1];
						}else{	//左邊為背景，給新標記
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
							label[labelCnt] = labelCnt;
							labelCnt ++;
						}
							
					}else if(i != 0 && j == 0){		//第0欄，檢查上面
						up = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j];
						if(up == 255){	//上面是前景，和上面一樣
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j];
						}else{			//上面是背景，給新標記
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
							label[labelCnt] = labelCnt;
							labelCnt ++;
						}
					}else{		//非邊邊點
						up = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j];
						left = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j-1];
						if(up == 0 && left == 0){	//上面、左邊都是背景，給新標記
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
							label[labelCnt] = labelCnt;
							labelCnt ++;
						}else{		//上面、左邊 至少一個是前景
							labelMin = 99999;

							//找出最小標籤
							if(up == 255)
								if(((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] <= labelMin)			
									labelMin = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j];			//最小標籤 = 上面
							if(left == 255)
								if(((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1] <= labelMin)		
									labelMin = ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1];			//最小標籤 = 左邊

							//註記所有標記相同
							if(label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] ] = label[ labelMin ];		//上面 = 最小標籤
							if(label[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1] ] = label[ labelMin ];		//左邊 = 最小標籤

							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelMin;
						}
						
						
					}
				}
			}
		}
	}else{
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				if(((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] == 255){
					
					if(i == 0 && j == 0){	//第1個pixel，直接給新標記
						((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
						label[labelCnt] = labelCnt;
						labelCnt ++;
					}else if(i == 0 && j != 0){	//第0列，只檢查左邊狀態
						left = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j-1];
						if(left == 255){	//左邊為前景，當格和左邊一樣
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1];
						}else{	//左邊為背景，給新標記
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
							label[labelCnt] = labelCnt;
							labelCnt ++;
						}	
					}else if(i != 0 && j == 0){		//第0欄，檢查上面、右上
						up = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j];
						rightUp = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j + 1];
						if(up == 0 && rightUp){		//上面、右上是背景，給新標記
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
							label[labelCnt] = labelCnt;
							labelCnt ++;
						}else{	//上面、右上 至少有一個是前景
							labelMin = 99999;
							//找出最小標籤
							if(up == 255)
								if(((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] <= labelMin)			
									labelMin = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j];			//最小標籤 = 上面
							if(rightUp == 255)
								if(((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1] <= labelMin)			
									labelMin = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1];			//最小標籤 = 右上

							//註記所有標記相同
							if(label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] ] = label[ labelMin ];		//上面 = 最小標記
							if(label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1] ] = label[ labelMin ];		//右上 = 最小標記

							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelMin;
						}
					}else{		//非邊邊點
						rightUp = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j + 1];
						up = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j];
						leftUp = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j - 1];
						left = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j-1];
						if(rightUp == 0 && up == 0 && leftUp == 0 && left == 0){	//右上、上面、左上、左邊都是背景，給新標記
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
							label[labelCnt] = labelCnt;
							labelCnt ++;
						}else{		//右上、上面、左上、左邊 至少一個是前景
							labelMin = 99999;

							//找出最小標籤
							if(rightUp == 255)
								if(((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1] <= labelMin)			
									labelMin = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1];		//最小標籤 = 右上
							if(up == 255)
								if(((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] <= labelMin)			
									labelMin = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j];			//最小標籤 = 上面
							if(leftUp == 255)
								if(((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j - 1] <= labelMin)			
									labelMin = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j - 1];		//最小標籤 = 左上
							if(left == 255)
								if(((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1] <= labelMin)		
									labelMin = ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1];			//最小標籤 = 左邊

							//註記所有標記相同
							if(label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1] ] = label[ labelMin ];	//右上 = 最小標記
							if(label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] ] = label[ labelMin ];		//上面 = 最小標記
							if(label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j - 1] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j - 1] ] = label[ labelMin ];	//左上 = 最小標記
							if(label[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1] ] = label[ labelMin ];		//左邊 = 最小標記
							

							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelMin;
						}
					}
				}
			}
		}
	}

	

	//將相同類別標記，改成同樣標記
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			labelCnt = ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j];
			while(1){
				if(labelCnt == label[ labelCnt ])
					break;
				else
					labelCnt = label[ labelCnt ];
			}
			((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
		}
	}

	
	

	//計算物件個數
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			if(((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j]){
				objCnt[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] ] = 1;
			}
		}
	}
	for(i = 0; i < 99999; i++){
		sum += objCnt[i];
	}

	//計算物件面積(pixel數)
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			area[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] ] ++;
		}
	}

	//去掉雜訊
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			if(area[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] ] < 500)
				((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 0;
			else
				((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
		}
	}

	if(!cvSaveImage("out_image.bmp",out_image)) printf("Could not save: %s\n", "out_image.bmp");

	cvReleaseImage(&label_image );	//釋放圖像記憶體
	return;

						
}



void MakeLut(Histogram *pLUT,int Min, int Max, int NumBins)
{
	int i;
	int BinSize;

	BinSize = 1 + (Max - Min) / NumBins;
	for(i = Min; i <= Max; i++)
		pLUT->arr[i] = (i-Min) / BinSize;

	return;
}

void MakeHistogram(IplImage *read_image, int Xstart, int Ystart, int Xsize, int Ysize, Histogram *pHist, int NumBins, Histogram *pLUT)
{
	int i,j;

	for(i=0; i < NumBins; i++)
		pHist->arr[i] = 0;

	for(i = Ystart; i < Ystart + Ysize; i++){
		for(j = Xstart; j < Xstart + Xsize; j++){
			pHist->arr[ pLUT->arr[ ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] ] ]++;
		}
	}

	FILE *fp = fopen("hist.txt","w");

	for(i=0; i < 256; i++){
		fprintf(fp, "%d ", pHist->arr[i]);
		
	}

	fclose(fp);

	return;
	

}

void ClipHistogram(Histogram *pHist, int NumBins, int ClipLimit)
{
	int i,j;
	Histogram *Histo;
	int NumExcess, OldNumExcess, BinExcess, Upper, BinIncr, StepSize;
	int Start, End;
	int sum;

	

	NumExcess = 0;
	for(i=0; i < NumBins; i++){
		BinExcess = pHist->arr[i] - ClipLimit;
		if(BinExcess > 0){
			NumExcess += BinExcess;
		}
	}

	BinIncr = NumExcess / NumBins;
	Upper = ClipLimit - BinIncr;

	for(i=0; i < NumBins; i++){
		if(pHist->arr[i] > ClipLimit)
			pHist->arr[i] = ClipLimit;
		else{
			if(pHist->arr[i] > Upper){
				NumExcess -= (ClipLimit - pHist->arr[i]);
				pHist->arr[i] = ClipLimit;
			}else{
				NumExcess -= BinIncr;
				pHist->arr[i] += BinIncr;
			}

		}
	}

	
	do{
		Start = 0;
		OldNumExcess = NumExcess;

		while(NumExcess && (Start < NumBins)){
			StepSize = NumBins / NumExcess;
			if(StepSize < 1)
				StepSize = 1;

			for(i = Start; i < NumBins && NumExcess; i += StepSize){
				if(pHist->arr[i] < ClipLimit){
					pHist->arr[i] ++;
					NumExcess --;
				}
			}
			Start ++;

		}

	}while(NumExcess && (NumExcess < OldNumExcess));


	/*

	Start = 0;
	while( NumExcess && (Start < NumBins)){
		StepSize = NumBins / NumExcess;
		if(StepSize < 1)	StepSize = 1;

		for(i = Start; (i < NumBins) && NumExcess; i += StepSize){
			if(pHist->arr[i] < ClipLimit){
				pHist->arr[i] ++;
				NumExcess --;
			}
		}
		Start ++;
	}

	sum = 0;
	for(i=0; i < NumBins; i++){
		sum += pHist->arr[i];
	}
	
	FILE *fp = fopen("hist.txt","w");

	for(i=0; i < 256; i++){
		fprintf(fp, "%d ", pHist->arr[i]);
		
	}

	fclose(fp);*/

	return;
}

void MapHistogram(Histogram *pHist, int Min, int Max, int NumBins, int NumPixels)
{
	int i;
	int sum = 0;

	for(i=0; i < NumBins; i++){
		sum += pHist->arr[i];
		pHist->arr[i] = Min + (sum * (Max - Min)) / NumPixels;
		if(pHist->arr[i] > Max)
			pHist->arr[i] = Max;
	}

	FILE *fp = fopen("hist.txt","w");

	for(i=0; i < 256; i++){
		fprintf(fp, "%d ", pHist->arr[i]);
		
	}

	fclose(fp);

	return;
}

void Interpolate(IplImage *read_image, IplImage *out_image, int Xstart, int Ystart, Histogram *pLU, Histogram *pLD, Histogram *pRU, Histogram *pRD, int Xsize, int Ysize, int NumPixels, Histogram *pLUT)
{
	int i,j;
	int Xcoef, Xbar, Ycoef, Ybar;
	int GreyValue;

	for(Ycoef = 0, Ybar = Ysize, i = Ystart; Ycoef < Ysize; Ycoef++, Ybar--, i++){
		for(Xcoef = 0, Xbar = Xsize, j = Xstart; Xcoef < Xsize; Xcoef++, Xbar--, j++){

			GreyValue = pLUT->arr[ ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] ];
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = (Xbar * Ybar * pLU->arr[GreyValue] + 
																			  Xcoef * Ybar * pRU->arr[GreyValue] +
																			  Xbar * Ycoef * pLD->arr[GreyValue] + 
																			  Xcoef * Ycoef * pRD->arr[GreyValue]) / NumPixels;	
		}
	}

	return;
}

void CLAHE(IplImage *read_image, int NumX, int NumY, int NumBins, float fClipLimit , IplImage *out_image)
{
	int i,j;

	int XRes = read_image->width;
	int YRes = read_image->height;
	
	int Min, Max;
	int Xsize, Ysize, subX, subY;
	int Xstart, Ystart;
	int ClipLimit, NumPixels;
	int up, down, left, right;
	Histogram *pLU, *pLD, *pRU, *pRD; 
	//unsigned int XL, XR, YU, YD;
	
	//IplImage *ImPointer;
	
	Histogram *pLUT = (Histogram *)malloc(sizeof(Histogram));
	Histogram **MapArray = (Histogram **)malloc(sizeof(Histogram *) * NumY);	
	for(i=0; i < NumY; i++){
		MapArray[i] = (Histogram *)malloc(sizeof(Histogram) * NumX);
	}
	Histogram *pHist;
	//Histogram *pLU, *pLD, *pRU, *pRD;



	Min = 99999;
	Max = -9999;
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			if(((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] < Min)	Min = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
			if(((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] > Max) Max = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
		}
	}
	Xsize = XRes / NumX;
	Ysize = YRes / NumY;
	NumPixels = Xsize * Ysize;
	ClipLimit = (fClipLimit * Xsize * Ysize) / NumBins;


	/*for(i=0; i < 256; i++)
		pLUT->arr[i] = 0;*/
	MakeLut(pLUT, Min, Max, NumBins);
	for(i=0; i < NumY; i++){
		for(j=0; j < NumX; j++){
			pHist = &MapArray[i][j];
			Xstart = j * Xsize;
			Ystart = i * Ysize;
			MakeHistogram(read_image, Xstart, Ystart, Xsize, Ysize, pHist, NumBins, pLUT);
			ClipHistogram(pHist, NumBins, ClipLimit);
			MapHistogram(pHist, Min, Max, NumBins, NumPixels);

		}
	}

	for(i=0; i < NumY; i++){

		if(i == 0){
				subY = Ysize >> 1;
				up = 0;
				down = 0;
		}else{
				
				if(i == NumY - 1){
					subY = Ysize >> 1;
					up = NumY - 1;
					down = up;
				}else{
					subY = Ysize;
					up = i - 1;
					down = i + 1;
				}
		}

		for(j = 0; j < NumX; j++){
			
			if(j == 0){
				subX = Xsize >> 1;
				left = 0;
				right = 0;
			}else{
				
				if(j == NumX - 1){
					subX = Xsize >> 1;
					left = NumX - 1;
					right = left;
				}else{
					subX = Xsize;
					left = j - 1;
					right = j + 1;
				}
			}

			pLU = &MapArray[up][left];
			pRU = &MapArray[up][right];
			pLD = &MapArray[down][left];
			pRD = &MapArray[down][right];

			Xstart = j * Xsize;
			Ystart = i * Ysize;

			Interpolate(read_image, out_image, Xstart, Ystart, pLU, pLD, pRU, pRD, Xsize, Ysize, NumPixels, pLUT);
		}
	}



		/*for(j=0; j <= NumX; j++){

			

			Xstart = j * Xsize;
			Ystart = i * Ysize;

			up = i-1;
			down = i+1;
			left = j-1;
			right = j+1;

			if( up < 0 )	up = 0;	
			if( down >= NumY ) down = NumY - 1;
			if( left < 0 ) left = 0;
			if( right >= NumX ) right = NumX - 1;

			pLU = &MapArray[up][left];
			pRU = &MapArray[up][right];
			pLD = &MapArray[down][left];
			pRD = &MapArray[down][right];

			Interpolate(read_image, out_image, Xstart, Ystart, pLU, pLD, pRU, pRD, Xsize, Ysize, NumPixels, pLUT);
		}
	}*/


	for (i=0; i < NumY; i++) free(MapArray[i]);
	free(MapArray);

	free(pLUT);


	return ;





}

void ALRHS(IplImage *read_image, int L, int H, int Z, IplImage *out_image)
{
	int i,j;
	int Le, He, Lnew, Hnew;
	int Hist[256] = {0};
	int LUT[256] = {0};


	Le = 99999;
	He = -9999;
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			if(((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] < Le) Le = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
			if(((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] > He) He = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
			Hist[ ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] ]++;
		}
	}

	Lnew = (L - Le + 1) / Z;
	Hnew = 255 - (He - H + 1) / Z;

	for(i=Le; i <= He; i++){
		if(i <= L){
			LUT[i] = (i - Le) / Z;
		}else if(i > L && i < H){
			LUT[i] = ((i - L) * (Hnew - Lnew)) / (H - L) + Lnew;
		}else{
			LUT[i] = 255 - (He - i) / Z;
		}
	}

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = LUT[((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j]];
		}
	}

	return;
}

void project_statistics(IplImage *read_image)	//取ROI後的血管中心定位
{
	int i,j;
	int grayProj[480] = {0};
	int sum;
	int MaxIndex;
	int Max;
	int cunt;
	int var;


	for(j=0; j < read_image->width; j++){
		for(i=0; i < read_image->height; i++){
			grayProj[i] += ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
		}
	}

	for(i=0; i < read_image->height; i++){
		grayProj[i] /= read_image->width;
	}

	Max = -999;
	for(i=0; i < read_image->height; i++){
		if(grayProj[i] > Max) Max = grayProj[i];
	}

	cunt = 0;
	MaxIndex = 0;
	for(i=0; i < read_image->height; i++){
		if(grayProj[i] == Max){
			cunt ++;
			MaxIndex += i;
		}
	}
	MaxIndex /= cunt;


	sum = 0;
	var = 0;
	for(i=0; i < read_image->height; i++){
		sum =  sum + (grayProj[i] - Max)*(grayProj[i] - Max);
	}
	sum /= (read_image->height - 1);
	var = sqrt(float(sum));







	FILE *fp = fopen("grayProj.txt","w");

	for(i=0; i < read_image->height; i++){
		fprintf(fp, "%d	", grayProj[i]);

	}

	fclose(fp);

	return;

}


/*
int hist_mode(IplImage *read_image)
{
	int i,j;
	int min,max;

	int *hist = (int *)malloc(sizeof(int) * 256);		//直條圖統計原始資料

	//初始化
	for(i=0; i < 256; i++){
		hist[i] = 0;
	}

	//灰階值統計
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			hist[ ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] ] ++;
		}
	}

	max = 0;
	for(i = 0; i < 256; i++){
		if(hist[i] >= max)
			max = hist[i];
		else
			break;

	}
}*/