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
�@�@�@�@�ΡG���o�@���}�C�ʺA�O����
�@�@�@�ѼơGsize = �@���}�C�j�p
			item_size = �}�C�����j�p
�@�@�@�Ǧ^�G�O����Ŷ����СA�Y��NULL��O���餣��
	  ����Gsize*item_size�����p��65536
�@�@�@------------------------------------------------------------ 
	int i, j;
	complex * ptr;

	ptr = (complex *) malloc(item_size * size);
	if (ptr == NULL){
		printf("�O���餣��, �������o���O����!\n");
		free(ptr);
		return NULL;
	}
	
	return ptr;
}


complex **Alloc2DArray(int xsize,int ysize,int item_size)
{
	/* ------------------------------------------------------------
�@�@�@�@�ΡG���o�G���}�C�ʺA�O����
�@�@�@�ѼơGxsize, ysize = �G���}�C�j�p
			item_size = �}�C�����j�p
�@�@�@�Ǧ^�G�O����Ŷ����СA�Y��NULL��O���餣��
	  ����Gxsize*item_size�����p��65536
�@�@�@------------------------------------------------------------ 
	
	int i, j;

	complex ** ptr;

	ptr = (complex**) malloc(sizeof(complex*) * ysize);

	if (ptr == NULL){
		printf("���o�O���饢��!\n");	
		return NULL;
	}
	
	for (i=0; i<ysize; i++){
		ptr[i] = (void*) malloc(item_size * xsize);
		if (ptr[i] == NULL){
			printf("�O���餣��, �������o���O����\n"); 
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
�@�@�@�ΡG����G���}�C�ʺA�O����
�@�@�ѼơGxsize, ysize = �G���}�C�j�p
�@�@�@�@�@ptr = �O����Ŷ�����
�@�@�Ƶ��Gxsize��ڤW�èS�Ψ�A�����e�����ѨϥΡA�G�[�W���@�Ѽ�
   ------------------------------------------------------------ 
	
	int i;

	if(ptr == NULL){
		printf("�ū���!!\n");
		return;
	}
	for (i=0; i<ysize; i++) free(ptr[i]);
	free(ptr);
}*/

//������ഫ�Ƕ���
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

//��ܼv��
void imshow(IplImage *read_image, char *windowName)
{
	cvNamedWindow(windowName,0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture1"
	cvResizeWindow(windowName, read_image->width, read_image->height);	//�]�w�����j�p
	cvShowImage(windowName, read_image);	//��ܹϹ�
	cvWaitKey(0);

	cvDestroyWindow(windowName);	//��������

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


//�ন��q�D�Ƕ�
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
			y=0.114*b + 0.587*g + 0.299*r;  										//�Ƕ��ഫ����
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j]=y;
		}

	}

	return ;
}

//�G�Ȥ�
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

//�줸����
void bit_plane(IplImage *read_image, IplImage *out_image, int bitNumber)
{
	/************************************************************************************************************************

	��J: bitNumber ���Q�n���X���줸�����C���]bitNumber = 3�A�N����X�ĤT�줸����(c3)�C�̧C�줸������bitNumber = 0(c0)�C

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

//�Ŷ��ѪR��
void im_resize(IplImage *read_image, IplImage *out_image, int resolution)
{
	/********************************************************************************************

	��J: resolution ���Q�n�վ㪺�ѪR�סC���]resolution = 128�A�N��վ㬰128 X 128 �ѪR��

	*********************************************************************************************/
	int i,j;
	int tmpi,tmpj;
	int maskN;

	maskN = 256 / resolution;	//maskN�N��O��N��N���B�n�U�h���ѪR�סA�]�N�ON��N�̭����ȳ��n�ۦP
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

//���öq��
void quantizatioin(IplImage *read_image, IplImage *out_image, int grayNumber)
{
	
	/*****************************************************************************************************************************

	��J: grayNumber ���Q�n�q�ƪ��Ƕ��ơC�p�GgrayNumber = 128�A�N��ϥ�128�Ƕ��V��A�j�����v���O�ϥ�256�Ƕ��V��(grayNumber = 256)

	******************************************************************************************************************************/
	int value;	//�v��piexl����
	int unit;	//�@�������h��piexl��
	int range;	//�q�ƪ��϶�
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

//�V��(�b���)
void dither(IplImage *read_image, IplImage *out_image)
{
	IplImage *dither_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);

	/**************************
		�ШM�w�V��x�}������
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

//�q�Ʀh�h�Ƕ��V��
void multi_dither(IplImage *read_image, IplImage *out_image, int grayNumber)
{
	/*****************************************************************************************************************************

	��J: grayNumber ���Q�n�q�ƪ��Ƕ��ơC�p�GgrayNumber = 128�A�N��ϥ�128�Ƕ��V��A�j�����v���O�ϥ�256�Ƕ��V��(grayNumber = 256)

	******************************************************************************************************************************/
	IplImage *normal_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	IplImage *dist_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);

	int i,j,k;
	int tmpi, tmpj;
	int normalMatrix[2][2];
	
	int dist;
	int unit;	//�@�������h��piexl��
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

//�~�t�X���k
void error_diff(IplImage *read_image, IplImage *out_image)
{
	IplImage *error_image = cvCreateImage(cvSize(read_image->width + 2, read_image->height + 2), IPL_DEPTH_8U,	1);

	int i,j;
	int e;	//�~�t�x�}�Ƕ��ȩM�q�ƭȪ��~�t

	//�N�~�t�x�}��l�Ƭ�0
	for(i=0; i < error_image->height; i++){
		for(j=0; j < error_image->width; j++){

			((uchar *)(error_image->imageData + i*error_image->widthStep))[j] = 0;

		}

	}

	//Copy��v��(n x n)�ܻ~�t�x�}(n+2 x n+2)
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j+1] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];

		}

	}

	
	for(i=1; i < error_image->height - 1; i++){
		for(j=1; j < error_image->width - 1; j++){
			if(((uchar *)(error_image->imageData + i*error_image->widthStep))[j] > 128){	//����~�t�x�}�Ƕ��ȩM128�A�p�G�j��128�A�N��q�ƭȬ�255
				((uchar *)(out_image->imageData + (i-1)*out_image->widthStep))[j-1] = 255;	//�N�~�t�x�}�q�ƫ�A��X�v��
				e = ((uchar *)(error_image->imageData + i*error_image->widthStep))[j] - 255;	//�p��~�t�x�}�Ƕ��ȩM�q�ƭȪ��~�t
			}else{
				((uchar *)(out_image->imageData + (i-1)*out_image->widthStep))[j-1] = 0;	
				e = ((uchar *)(error_image->imageData + i*error_image->widthStep))[j];
			}

			//�վ�~�t�x�}
			((uchar *)(error_image->imageData + i*error_image->widthStep))[j+1] = is_overflow(((uchar *)(error_image->imageData + i*error_image->widthStep))[j+1] + (e*7)/16);
			((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j-1] = is_overflow(((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j-1] + (e*3)/16);
			((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j] = is_overflow(((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j] + (e*5)/16);
			((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j+1] = is_overflow(((uchar *)(error_image->imageData + (i+1)*error_image->widthStep))[j+1] + e/16);
			
		}
	}
	cvReleaseImage(&error_image);
	return;
}

//�v���[�k
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

//�v����k
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

//�v�����k
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

//�v�����k
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

//�v���ɦ�
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

//�έp�����
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

//������X�i�k
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

//������X�iHS
void hist_stretch(IplImage *read_image, IplImage *out_image)
{
	

	int i,j;

	int Hist[256] = {0};	//�έp�Ƕ��X�{���Ưx�}
	int LUT[256] = {0};
	int greyMin, greyMax;
	int NumPixels;
	int sum;

	NumPixels = read_image->width * read_image->height;

	//�έp�Ƕ��X�{����
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

	

	/*********�ϥβֿn��ҩ��i***********/
	;
	sum = 0;
	for(i=0; i < 256; i++){
		if(Hist[i]){
			sum += Hist[i];
			LUT[i] = (255 * sum) / NumPixels;
		}
	}

	/*********�ϥνu�q��ҩ��i***********/
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
	int delt;	//�ھکP�򹳯������h��ܹ�����
	int low,high;	//�B�z���h���d��
	int BinIncr;	//���ƫ�A1�ӦǶ����h��������
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
		sort(buf_image, &buf[0], low);	//���өP�򹳯����Ƕ��Ѱ��V�C�ƦC�ܤ�

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

//LUT����ϵ���
void LUT_hist_eq(IplImage *read_image, IplImage *out_image)
{
	FILE *out_fp;
	

	int grayNumberArray[256] = {0};	//�έp�Ƕ��X�{���ư}�C
	int cdfArray[256] = {0};
	int equArray[256] = {0};	//����ϵ��ư}�C
	int LUT[256] = {0};
	int i,j;
	int flagMin, flagMax ;	//flagMin���Ƕ��̤p�ȵo�ͪ��a��AflagMax���Ƕ��̤j�ȵo�ͪ��a��
	int cdfNumber = 0, cdfMin, cdfMax;	//cdfNumber���ΨӲֿn�p��cdf���Ʀr�AcdfMin���ֿn������Ƴ̤p�ȡAcdfMax���ֿn������Ƴ̤j��(���]�N�Opiexl�`��)
	int changedGray;

	out_fp = fopen("out.txt", "w");
	

	//�ҥ�enginner.tif�m��
	;
	/***********************************
	FILE *out_divi_fp;
	out_divi_fp = fopen("out_divi.txt", "w");
	IplImage *tmp2_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	cvNamedWindow("showpicture3",0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture1"
	cvResizeWindow("showpicture3", read_image->width, read_image->height);	//�]�w�����j�p
	cvNamedWindow("showpicture4",0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture1"
	cvResizeWindow("showpicture4", read_image->width, read_image->height);	//�]�w�����j�p
	cvNamedWindow("showpicture5",0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture1"
	cvResizeWindow("showpicture5", read_image->width, read_image->height);	//�]�w�����j�p
	int tmpCalcArray[256] = {0};

	im_divide(read_image, tmp2_image, 4);
	cvShowImage("showpicture3", read_image);	//��ܹϹ�
	cvShowImage("showpicture4", tmp2_image);	//��ܹϹ�
	cvWaitKey(0);
	//�έp�Ƕ��X�{����
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

	//�έp�Ƕ��X�{����
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			grayNumberArray[ ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] ] ++;
		}

	}

	//�p����v�������pdf
	for(i=0; i < 256; i++){
		cdfNumber = cdfNumber + grayNumberArray[i];
		cdfArray[i] = cdfNumber;
		
	}

	//���Ƕ��̤p�ȵo�ͪ��a��B�ֿn������Ƴ̤p��
	for(i=0; i < 256; i++){
		if( grayNumberArray[i] ){
			flagMin = i;
			cdfMin = cdfArray[i];
			break;
		}
	}

	//���Ƕ��̤j�ȵo�ͪ��a��B�ֿn������Ƴ̤j��
	for(i=255; i >= 0; i--){
		if( grayNumberArray[i] ){
			flagMax = i;
			cdfMax = cdfArray[i];
			break;
		}
	}

	//�u��P���쪺�a�谵����ϵ���
	for(i = flagMin; i <= flagMax; i++){
		if(grayNumberArray[i]){
			changedGray = (int)((((float)cdfArray[i] - cdfMin) * 255) / ((float)cdfMax - cdfMin) + 0.5);
			LUT[i] = changedGray;
		}
	}

	//�d��LUT���A���۹������Ƕ���
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){

			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = LUT[ ((uchar *)(read_image->imageData + i*read_image->widthStep))[j] ];
		}

	}

	
	//�έp�d��LUT���ƫ᪺�X�{����(�]�N�O����ϵ���)
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


//�����o�i��
void avg_filter(IplImage *read_image, IplImage *out_image, int dim)
{
	
	/*****************************************************************************************************************************

	��J: shape �����w�v����t�������B�z��k�C-1:������t 0:�ɹs 1:��g
		  dim �����w�B�n�j�p�C�Ydim = 5�A�h����5x5�B�n

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
	

	
	
	//��l�ƾB�n
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


	cvReleaseImage(&copy_image);		//����Ϲ��O����

	for (i=0; i < dim; i++) free(avgFilter[i]);
	free(avgFilter);

	return;


}

//Laplacian�o�i��
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
			if(i == 0 && j == 0){	/*****���W��****/
				for(tmpi = i; tmpi < (dim / 2); tmpi++){
					for(tmpj = j; tmpj < (dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*tmpA_image->widthStep))[j];
					}
				}
			}else if(i == 0 && j == read_image->width - 1){	/*****�k�W��****/
				for(tmpi = 0; tmpi < (dim / 2); tmpi++){
					for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					}
				}
			}else if(i == read_image->height - 1 && j == 0){	/*****���U��****/
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					for(tmpj = 0; tmpj < (dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*tmpA_image->widthStep))[j];
					}
				}
			}else if(i == read_image->height - 1 && j == read_image->width - 1){	/*****�k�U��****/
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					}
				}
			}else if(i == 0){	/*****�W��****/
				for(tmpi = i; tmpi < (dim / 2); tmpi++){
					((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[j] = ((uchar *)(read_image->imageData + (0)*read_image->widthStep))[j];
				}
			}else if(j == read_image->width - 1){	/*****�k��****/
				for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
					((uchar *)(tmpA_image->imageData + i*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[read_image->width - 1];
				}
			}else if(i == read_image->height - 1){	/*****�U��****/
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[j] = ((uchar *)(read_image->imageData + (read_image->height - 1)*read_image->widthStep))[j];
				}
			}else if(j == 0){	/*****����****/
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

	
	cvReleaseImage(&tmpA_image);		//����Ϲ��O����
	
	return;
}

//LoG�o�i
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
			if(i == 0 && j == 0){	/*****���W��****/
				for(tmpi = i; tmpi < (dim / 2); tmpi++){
					for(tmpj = j; tmpj < (dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*tmpA_image->widthStep))[j];
					}
				}
			}else if(i == 0 && j == read_image->width - 1){	/*****�k�W��****/
				for(tmpi = 0; tmpi < (dim / 2); tmpi++){
					for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					}
				}
			}else if(i == read_image->height - 1 && j == 0){	/*****���U��****/
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					for(tmpj = 0; tmpj < (dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*tmpA_image->widthStep))[j];
					}
				}
			}else if(i == read_image->height - 1 && j == read_image->width - 1){	/*****�k�U��****/
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					}
				}
			}else if(i == 0){	/*****�W��****/
				for(tmpi = i; tmpi < (dim / 2); tmpi++){
					((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[j] = ((uchar *)(read_image->imageData + (0)*read_image->widthStep))[j];
				}
			}else if(j == read_image->width - 1){	/*****�k��****/
				for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
					((uchar *)(tmpA_image->imageData + i*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[read_image->width - 1];
				}
			}else if(i == read_image->height - 1){	/*****�U��****/
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[j] = ((uchar *)(read_image->imageData + (read_image->height - 1)*read_image->widthStep))[j];
				}
			}else if(j == 0){	/*****����****/
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

	
	cvReleaseImage(&tmpA_image);		//����Ϲ��O����
	return;
}

//���q�o�i��
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

	//�t�m�O����Ŷ�������ഫ�᪺image
	;
	/*
	tmpB_image = (int **)malloc(sizeof(int *) * read_image->height);	
	for(i=0; i < read_image->height; i++){
		tmpB_image[i] = (int *)malloc(sizeof(int ) * read_image->width);
	}*/

	;
	//copy�@�i����j����
	;
	/*
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			if(i == 0 && j == 0){	//*****���W��****
				for(tmpi = i; tmpi < (dim / 2); tmpi++){
					for(tmpj = j; tmpj < (dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*tmpA_image->widthStep))[j];
					}
				}
			}else if(i == 0 && j == read_image->width - 1){	//****�k�W��***
				for(tmpi = 0; tmpi < (dim / 2); tmpi++){
					for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					}
				}
			}else if(i == read_image->height - 1 && j == 0){	//*****���U��****
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					for(tmpj = 0; tmpj < (dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*tmpA_image->widthStep))[j];
					}
				}
			}else if(i == read_image->height - 1 && j == read_image->width - 1){	//*****�k�U��****
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
						((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
					}
				}
			}else if(i == 0){	//*****�W��****
				for(tmpi = i; tmpi < (dim / 2); tmpi++){
					((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[j] = ((uchar *)(read_image->imageData + (0)*read_image->widthStep))[j];
				}
			}else if(j == read_image->width - 1){	//*****�k��****
				for(tmpj = read_image->width + (dim / 2); tmpj < read_image->width + 2*(dim / 2); tmpj++){
					((uchar *)(tmpA_image->imageData + i*tmpA_image->widthStep))[tmpj] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[read_image->width - 1];
				}
			}else if(i == read_image->height - 1){	//*****�U��****
				for(tmpi = read_image->height + (dim / 2); tmpi < read_image->height + 2*(dim / 2); tmpi++){
					((uchar *)(tmpA_image->imageData + tmpi*tmpA_image->widthStep))[j] = ((uchar *)(read_image->imageData + (read_image->height - 1)*read_image->widthStep))[j];
				}
			}else if(j == 0){	*****����****
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

	
	
	
	//scal_trans(tmpB_image, out_image);	//������ഫ�Ƕ���

	return;
}

//�h�U�Q��
void unsharp(IplImage *read_image, IplImage *out_image, int gain, float k)
{
	/*****************************************************************************************************************************

	��J: gain ���վ��v������j��v
		  k ���վ��� k < 1

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

	
	

	//�h�U�Q�Ƽv�� = ������X�h�U�Q�e�B�n�A�ξB�n�U�h�o�i
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
	
	//�h�U�Q�Ƽv�� = ��l�v����jgain�� - �����o�i�v����j1/k��()
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

	cvNamedWindow("read",0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture1"
	cvResizeWindow("read", read_image->width, read_image->height);	//�]�w�����j�p
	cvShowImage("read", read_image);	//��ܹϹ�

	cvNamedWindow("avg",0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture1"
	cvResizeWindow("avg", read_image->width, read_image->height);	//�]�w�����j�p
	cvShowImage("avg", avg_image);	//��ܹϹ�

	cvNamedWindow("unsharp",0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture1"
	cvResizeWindow("unsharp", read_image->width, read_image->height);	//�]�w�����j�p
	cvShowImage("unsharp", out_image);	//��ܹϹ�

	if(!cvSaveImage("unsharp1.bmp",out_image)) printf("Could not save: %s\n", "out.bmp");

	cvWaitKey(0);

	cvDestroyWindow("read");	//��������
	cvDestroyWindow("unsharp");	//��������

	cvReleaseImage(&avg_image);		//����Ϲ��O����

	****************************************************/
	

	return;


}

//���W�T�o�i
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

	cvReleaseImage(&avg_image);		//����Ϲ��O����
	return;
}

//�̤j�̤p�o�i��
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

//�����o�i��
void mid_filter(IplImage *read_image, IplImage *out_image, int dim)
{
	/*****************************************************************************************************************************

	��J: shape �����w�v����t�������B�z��k�C-1:������t 0:�ɹs 1:��g
		  dim �����w�B�n�j�p�C�Ydim = 5�A�h����5x5�B�n

	******************************************************************************************************************************/

	int tmp;
	int A,B;
	int i,j,tmpi,tmpj;
	int row, col;
	int m,n;
	int filDex;

	IplImage *copy_image = cvCreateImage(cvSize(read_image->width + (2*(dim / 2)), read_image->height + (2*(dim / 2))), IPL_DEPTH_8U,	1);

	int *midFilter = (int *)malloc(sizeof(int ) * (dim*dim));
	
	

	
	
	//��l�ƾB�n
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

	cvReleaseImage(&copy_image);		//����Ϲ��O����

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

	cvReleaseImage(&copy_image);		//����Ϲ��O����

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

	cvReleaseImage(&copy_image);		//����Ϲ��O����

	return;
}

void band_reject_folter(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	float distance;
	float **filter = (float **)malloc(sizeof(float *) * read_image->height);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < in_no; i++){
		filter[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	complex **read_DFT = (complex **)malloc(sizeof(complex *) * read_image->height);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < in_no; i++){
		read_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	complex **read_pass = (complex **)malloc(sizeof(complex *) * read_image->height);	//���o�G���}�C�ʺA�O����;
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

	DFT2(read_image, read_DFT);	//�G���ť߸��ഫ
	

	for(i=0;i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			read_pass[i][j].real = read_DFT[i][j].real * filter[i][j];
			read_pass[i][j].image = read_DFT[i][j].image * filter[i][j];
		}
	}
	fftshow(256, read_DFT,"read_DFT");	//����W�й�
	IDFT2(read_pass, out_image);	//�G���ť߸��u�ϡv�ഫ

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
	float **filter = (float **)malloc(sizeof(float *) * read_image->height);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < in_no; i++){
		filter[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	complex **read_DFT = (complex **)malloc(sizeof(complex *) * read_image->height);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < in_no; i++){
		read_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	complex **read_pass = (complex **)malloc(sizeof(complex *) * read_image->height);	//���o�G���}�C�ʺA�O����;
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

	DFT2(read_image, read_DFT);	//�G���ť߸��ഫ
	

	for(i=0;i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			read_pass[i][j].real = read_DFT[i][j].real * filter[i][j];
			read_pass[i][j].image = read_DFT[i][j].image * filter[i][j];
		}
	}
	fftshow(256, read_pass,"read_pass");	//����W�й�
	IDFT2(read_pass, out_image);	//�G���ť߸��u�ϡv�ഫ

	for (i=0; i < read_image->width; i++) free(filter[i]);
	free(filter);

	for (i=0; i < read_image->width; i++) free(read_DFT[i]);
	free(read_DFT);

	for (i=0; i < read_image->width; i++) free(read_pass[i]);
	free(read_pass);

	return;

}

//ROI�o�i
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
		case 1:	//�����o�i��

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


//�L�X�Ƽư}�C
void print_complex(int data_no, complex *in_data)
{
	int i;
	for(i=0; i < data_no; i++){
		printf("%f + %fi\n", in_data[i].real, in_data[i].image);
	}
	printf("\n");
	return;
}

//�@���ƦC����
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

//�@���ť߸��ഫ
void DFT1(int data_no, complex *in_data, complex *out_data)
{
	
	float angle_step;	//�p��-2*PI / T
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



//�@���u�ϡv�ť߸��ഫ
void iDFT1(int data_no, complex *in_data, complex *out_data)
{
	
	float angle_step;	//�p��-2*PI / T
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

//�ۿn
void conv(int data1_M, complex *data1, int data2_N, complex *data2, complex *out_data)
{
	//complex tmpA[ data1_M + 2*(data2_N - 1) ];
	//complex tmpB[ data1_M + 2*(data2_N - 1) ];
	complex *tmpA = (complex *)malloc(sizeof(complex) * (data1_M + 2*(data2_N - 1)));
	complex *tmpB = (complex *)malloc(sizeof(complex) * (data1_M + 2*(data2_N - 1)));
	complex sum;
	complex tmp;
	

	int i,j,k;


	//tmpB����½��A�åB�ɹs
	for(j = M + 2*(N - 1) - 1; j >= 0 ; j--){
		if( j < N ){
			tmpB[j].real = data2[N - j - 1].real;
			tmpB[j].image = data2[N - j - 1].image;
		}else{
			tmpB[j].real = 0;
			tmpB[j].image = 0;
		}
	}

	//tmpA�W�U���ɹs�A������Ʀr
	for(j=0; j < (M + 2*(N - 1)); j++){
		if( N-1 <= j && j <= M+N-2 ){
			tmpA[j].real = data1[j-M+1].real;
			tmpA[j].image = data1[j-M+1].image;
		}else{
			tmpA[j].real = 0;
			tmpA[j].image = 0;
		}
	}

	//�x�}�ۭ�
	for(i=0; i < M+N-1; i++){
		sum.real = 0;
		sum.image = 0;
		for(j=0; j < (M + 2*(N - 1)); j++){
			sum.real = sum.real + ( tmpB[j].real * tmpA[j].real - tmpB[j].image * tmpA[j].image );
			sum.image = sum.image + ( tmpB[j].real * tmpA[j].image + tmpB[j].image * tmpA[j].real );
		}
		out_data[i].real = sum.real;
		out_data[i].image = sum.image;

		//�����@�����k����A�V�q�V�k��1��
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

//���αۿn
void cconv(int data_no, complex *data1, complex *data2, complex *out_data)
{
	complex *tmpB = (complex *)malloc(sizeof(complex) * data_no);
	complex sum;
	complex tmp;

	int i,j,k;

	//����½��
	for(j = data_no-1; j >= 0 ; j--){
		tmpB[j].real = data2[data_no - j - 1].real;
		tmpB[j].image = data2[data_no - j - 1].image;
	}

	//�M��tmpB�V�k�`��
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

		//tmpB�V�q����
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

//�G���ť߸��ഫ
void DFT2(IplImage *read_image, complex **out_data)
{
	/**************************************************************
	�\��: �G���ť߸��ഫ
	��J: data_no = ��J��ƪ��C/���ƭӼ�
		  read_image = ��J�v���A���A��IplImage *
	��X: out_data = �ť߸��ഫ���G�A���A��complex **

	**************************************************************/
	int i,j,k;
	int sign;
	int data_no = read_image->width;

	

	/****************************************************************
	�ʺA�O������t
	read_data = �Ƽƪ��G���x�}�A�ΨӦs��Ū�i�Ӫ��v��
	row		  = �Ƽƪ��@���x�}�A�ΨӼȦs���C�V�q
	col		  = �Ƽƪ��@���x�}�A�ΨӼȦs����V�q
	tmp		  = �Ƽƪ��@���x�}�A�ΨӼȦs�C���p��᪺���G
	****************************************************************/
	complex **read_data = (complex **)malloc(sizeof(complex *) * data_no);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < data_no; i++){
		read_data[i] = (complex *)malloc(sizeof(complex) * data_no);
	}

	complex *row = (complex *)malloc(sizeof(complex) * data_no);	//���o�@���}�C�ʺA�O����
	complex *col = (complex *)malloc(sizeof(complex) * data_no);	//���o�@���}�C�ʺA�O����
	complex *tmp = (complex *)malloc(sizeof(complex) * data_no);	//���o�@���}�C�ʺA�O����

	
	/****************************************************************
								������
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
	�ƻsread_data���C�V�q �� row
	��row���ť߸��ഫ �� tmp
	�ƻstmp �� out_data���C�V�q
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
	�ƻsout_data����V�q �� col
	��col���ť߸��ഫ �� tmp
	�ƻstmp �� out_data����V�q
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
				����O����
	******************************************/
	for (i=0; i < data_no; i++) free(read_data[i]);
	free(read_data);

	free(row);
	free(col);
	free(tmp);
	return;
}

//Ū���o�i��
void load_filter(FILE *datap, float **out_data)
{
	/*****************************************************************************************************************************
	�\��: Ū���o�i��
	��J: datap = �nŪ�����o�i���ɦW
	��X: out_data = �o�i���A���A��float **
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

//����W�й�
void fftshow(int data_no, complex **in_data, char *windowName)
{
	/*****************************************************************************************************************************
	�\��: ����W�й�
	��J: data_no = ��J��ƪ��C/���ƭӼ�
		  in_data = �n��ܪ��W�йϡA���A��complex **
	��X: �L
	******************************************************************************************************************************/
	IplImage *FreqImage = cvCreateImage(cvSize(data_no, data_no), IPL_DEPTH_8U,	1);

	

	int i,j;
	float max = -999;

	/****************************************************************
	�ʺA�O������t
	tmp = �B�I�ƪ��G���x�}�A�Ψ��x�s�ƼƱj�� |F(u,v)|
	****************************************************************/
	float **tmp = (float **)malloc(sizeof(float *) * data_no);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < data_no; i++){
		tmp[i] = (float *)malloc(sizeof(float) * data_no);
	}

	/****************************************************************
	�N�Ƽ��ഫ���j�� |F(u,v)| �A�P�ɧ��̤j��
	****************************************************************/
	for(i = 0; i < data_no; i++){
		for(j=0; j < data_no; j++){
			tmp[i][j] = log(1 + sqrt( in_data[i][j].real * in_data[i][j].real + in_data[i][j].image * in_data[i][j].image ));
			if(tmp[i][j] > max)
				max = tmp[i][j];
		}
	}

	/****************************************************************
	�N�C�� |F(u,v)| ���W�Ʀ�0~255�A��ܼv��
	****************************************************************/
	for(i=0; i < data_no; i++){
		for(j=0; j < data_no; j++){
			((uchar *)(FreqImage->imageData + i*FreqImage->widthStep))[j] = (tmp[i][j] * 255) / max;
		}
	}

	cvNamedWindow(windowName,0);	//�Ыؤ@�ӷs�����A�R�W��"FreqImage"
	cvResizeWindow(windowName, data_no, data_no);	//�]�w�����j�p
	cvShowImage(windowName, FreqImage);	//��ܹϹ�
	cvWaitKey(0);

	
	cvDestroyWindow(windowName);	//��������
	cvReleaseImage(&FreqImage);		//����Ϲ��O����
	for (i=0; i < data_no; i++) free(tmp[i]);	//����G���O����
	free(tmp);
    
	return;
}

//�G���ť߸��u�ϡv�ഫ
void IDFT2(complex **read_data, IplImage *out_image)
{
	/*****************************************************************************************************************************
	�\��: �G���ť߸��u�ϡv�ഫ
	��J: data_no = ��J��ƪ��C/���ƭӼ�
		  read_data = �o�i�᪺�W�йϡA���A��complex **
	��X: out_image = ��X�v���A���A��IplImage *

	******************************************************************************************************************************/
	int i,j,k;
	int sign;
	float max = -999;
	int data_no = out_image->width;

	/****************************************************************
	�ʺA�O������t
	out_data  = �Ƽƪ��G���x�}�A�Ψ��x�s���ഫ�᪺�ƭ�
	tmpInv	  = �B�I�ƪ��G���x�}�A�Ψ��x�s���W�Ƽƭ�
	row		  = �Ƽƪ��@���x�}�A�ΨӼȦs���C�V�q
	col		  = �Ƽƪ��@���x�}�A�ΨӼȦs����V�q
	tmp		  = �Ƽƪ��@���x�}�A�ΨӼȦs�C���p��᪺���G
	****************************************************************/

	complex **out_data = (complex **)malloc(sizeof(complex *) * data_no);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < data_no; i++){
		out_data[i] = (complex *)malloc(sizeof(complex) * data_no);
	}

	float **tmpInv = (float **)malloc(sizeof(float *) * data_no);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < data_no; i++){
		tmpInv[i] = (float *)malloc(sizeof(float) * data_no);
	}
	
	complex *row = (complex *)malloc(sizeof(complex) * data_no);
	complex *col = (complex *)malloc(sizeof(complex) * data_no);
	complex *tmp = (complex *)malloc(sizeof(complex) * data_no);



	/***************************************
	�ƻsread_data���C�V�q �� row
	��row���ť߸��ഫ �� tmp
	�ƻstmp �� out_data���C�V�q
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
	�ƻsout_data����V�q �� col
	��col���ť߸��ഫ �� tmp
	�ƻstmp �� out_data����V�q
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
			�N�ƭ��ഫ���쥻�b�Ŷ��쪺�ƭ� �A�P�ɧ��̤j��
	****************************************************************/
	for(i = 0; i < data_no; i++){
		for(j=0; j < data_no; j++){
			
			tmpInv[i][j] = sqrt( out_data[i][j].real * out_data[i][j].real + out_data[i][j].image * out_data[i][j].image );
			if(tmpInv[i][j] > max)
				max = tmpInv[i][j];
		}
	}

	/****************************************************************
				�N�C�Ӽƭȥ��W�Ʀ�0~255
	****************************************************************/
	for(i=0; i < data_no; i++){
		for(j=0; j < data_no; j++){
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = (tmpInv[i][j] * 255) / max;
		}
	}

	/****************************************************************
						����O����
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

	IplImage *read_image64F = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_64F,	1);	//���t�O���鵹tmp_image	��q�D

	float **lowpass_butterD15N2 = (float **)malloc(sizeof(float *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		lowpass_butterD15N2[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	float **lowpass_butterD40N10 = (float **)malloc(sizeof(float *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		lowpass_butterD40N10[i] = (float *)malloc(sizeof(float) * read_image->width);
	}
	

	complex **read_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		read_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	complex **read_lowpass = (complex **)malloc(sizeof(complex *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		read_lowpass[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	IplImage *noise_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//���t�O���鵹tmp_image	��q�D
	IplImage *noise_image64F = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_64F,	1);	//���t�O���鵹tmp_image	��q�D

	complex **noise_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		noise_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	complex **noise_lowpass = (complex **)malloc(sizeof(complex *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		noise_lowpass[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	FILE *lowpass_butterD15N2_ptr = fopen("butter_lowpass_D15N2_256.txt","r");
	FILE *lowpass_butterD40N10_ptr = fopen("butter_lowpass_D100N10_256.txt","r");


	load_filter(lowpass_butterD15N2_ptr, lowpass_butterD15N2);	//Ū���o�i��
	load_filter(lowpass_butterD40N10_ptr, lowpass_butterD40N10);	//Ū���o�i��


	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((float *)(read_image64F->imageData + (i)*read_image64F->widthStep))[j] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
		}
	}
	DFT2(read_image64F, read_DFT);	//�G���ť߸��ഫ
	fftshow(256, read_DFT, "read_DFT");	//����W�й�

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			read_lowpass[i][j].real = read_DFT[i][j].real * lowpass_butterD15N2[i][j];
			read_lowpass[i][j].image = read_DFT[i][j].image * lowpass_butterD15N2[i][j];
		}
	}
	//fftshow(256, read_lowpass, "read_lowpass");	//����W�й�
	
	IDFT2(read_lowpass, noise_image);	//�G���ť߸��u�ϡv�ഫ
	imshow(noise_image, "noise_image");	//��ܼv��


	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((float *)(noise_image64F->imageData + (i)*noise_image64F->widthStep))[j] = ((uchar *)(noise_image->imageData + (i)*noise_image->widthStep))[j];
		}
	}
	DFT2(noise_image64F, noise_DFT);	//�G���ť߸��ഫ
	fftshow(256, noise_DFT, "noise_DFT");	//����W�й�


	/********************��k2 �����]��********************************/
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
	//fftshow(256, noise_lowpass, "noise_lowpass");	//����W�й�
	

	IDFT2(noise_lowpass, out_image);	//�G���ť߸��u�ϡv�ഫ


	for (i=0; i < read_image->width; i++) free(lowpass_butterD15N2[i]);
	free(lowpass_butterD15N2);

	for (i=0; i < read_image->width; i++) free(lowpass_butterD40N10[i]);
	free(lowpass_butterD40N10);

	for (i=0; i < read_image->width; i++) free(read_DFT[i]);
	free(read_DFT);

	for (i=0; i < read_image->width; i++) free(read_lowpass[i]);
	free(read_lowpass);

	cvReleaseImage(&noise_image);		//����Ϲ��O����

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

	IplImage *noise_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//���t�O���鵹tmp_image	��q�D
	IplImage *noise_image64F = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_64F,	1);	//���t�O���鵹tmp_image	��q�D
	IplImage *motion_filter64F = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_64F,	1);	//���t�O���鵹tmp_image	��q�D

	complex **motion_filter_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		motion_filter_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	complex **noise_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		noise_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	complex **recover_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		recover_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	/*complex **noise_lowpass = (complex **)malloc(sizeof(complex *) * read_image->width);	//���o�G���}�C�ʺA�O����;
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

	//imshow(noise_image, "noise_image");	//��ܼv��

	for(i=0; i < motion_filter64F->height; i++){
		for(j=0; j < motion_filter64F->width; j++){
			if(i==0 & j < 7)
				((float *)(motion_filter64F->imageData + i*motion_filter64F->widthStep))[j] = 0.1429;
			else
				((float *)(motion_filter64F->imageData + i*motion_filter64F->widthStep))[j] = 0.0;
		}
	}

	DFT2(motion_filter64F, motion_filter_DFT);	//�G���ť߸��ഫ
	//fftshow(256, motion_filter_DFT, "motion_filter_DFT");	//����W�й�
	
	U8toF64(noise_image, noise_image64F);	//8U�v����64F
	DFT2(noise_image64F, noise_DFT);	//�G���ť߸��ഫ
	//fftshow(256, noise_DFT, "noise_DFT");	//����W�й�


	
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
	//fftshow(256, recover_DFT, "recover_DFT");	//����W�й�

	IDFT2(recover_DFT, out_image);	//�G���ť߸��u�ϡv�ഫ

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = ((uchar *)(out_image->imageData + i*out_image->widthStep))[j] * 3;
		}

	}


	cvReleaseImage(&noise_image);		//����Ϲ��O����
	cvReleaseImage(&noise_image64F);		//����Ϲ��O����
	cvReleaseImage(&motion_filter64F);		//����Ϲ��O����

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

	IplImage *read_image64F = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_64F,	1);	//���t�O���鵹tmp_image	��q�D
	complex **read_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		read_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	float **lowpass_butterD15N2 = (float **)malloc(sizeof(float *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		lowpass_butterD15N2[i] = (float *)malloc(sizeof(float) * read_image->width);
	}

	complex **noise_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		noise_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	IplImage *noise_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//���t�O���鵹tmp_image	��q�D

	IplImage *noise_image64F = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_64F,	1);	//���t�O���鵹tmp_image	��q�D
	complex **recover_DFT = (complex **)malloc(sizeof(complex *) * read_image->width);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < read_image->width; i++){
		recover_DFT[i] = (complex *)malloc(sizeof(complex) * read_image->width);
	}

	FILE *lowpass_butterD15N2_ptr = fopen("butter_lowpass_D15N2_256.txt","r");

	


	U8toF64(read_image, read_image64F);	//8U�v����64F
	DFT2(read_image64F, read_DFT);	//�G���ť߸��ഫ

	load_filter(lowpass_butterD15N2_ptr, lowpass_butterD15N2);	//Ū���o�i��

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			noise_DFT[i][j].real = read_DFT[i][j].real * lowpass_butterD15N2[i][j];
			noise_DFT[i][j].image = read_DFT[i][j].image * lowpass_butterD15N2[i][j];
		}	
	}

	IDFT2(noise_DFT, noise_image);	//�G���ť߸��u�ϡv�ഫ

	/*****************************************************************/

	U8toF64(noise_image, noise_image64F);	//8U�v����64F
	DFT2(noise_image64F, noise_DFT);	//�G���ť߸��ഫ

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

	IDFT2(recover_DFT, out_image);	//�G���ť߸��u�ϡv�ഫ

	return;


}

int Otsu(IplImage *read_image)
{
	int i,j,k;
	int Number;	//�`pixel��
	int sum_of_number;	//�ΨӲέp�ֿn����
	int start,end;
	float min;
	float sum;
	int outValue;	//��X�ȡA�]�N�O�V�m�n���֭�
	

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

	//��l��
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

	//�Ƕ��Ȳέp
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			number[ ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] ] ++;
		}
	}

	//��X�}�l�I
	for(i=0; i < 256; i++){
		if(number[i]){
			start = i;
			break;
		}
	}
	//��X�����I
	for(i=255; i > 0; i--){
		if(number[i]){
			end = i;
			break;
		}
	}

	
	//�p��Ƕ��Ⱦ��v
	for(i = start; i < end; i++){
		probability[i] = (float)number[i] / Number;
	}
	
	
	//�p��U�s���ֿn���v
	mean = 0;
	for(i = start; i <= end; i++){
		if(i == 0){
			accu_probabilityFore[i] = probability[i];
			accu_probabilityBack[i] = 1 - probability[i];
		}else{
			accu_probabilityFore[i] = accu_probabilityFore[i-1] + probability[i];	//�e���ֿn���v
			accu_probabilityBack[i] = 1 - accu_probabilityFore[i];					//�I���ֿn���v
		}
		mean = mean + (i * probability[i]);	//�`������(�����)
	}
	

	//�p��U�s��������
	for(i = start; i < end; i++){
		sum = 0;
		sum_of_number = 0;
		for(j = start; j <= i; j++){
			sum = sum + (j * number[j]);					//�Ƕ��� * �ӦǶ��Ȳέp����
			sum_of_number = sum_of_number + number[j];		//�e���`�ֿn����
		}
		meanFore[i] = sum / sum_of_number;	//�e��������

		sum = 0;
		sum_of_number = 0;
		for(k = j; k < end; k++){
			sum = sum + (k * number[k]);					//�Ƕ��� * �ӱ��Ȳέp����
			sum_of_number = sum_of_number + number[k];		//�I���`�ֿn����
		}
		meanBack[i] = sum / sum_of_number;	//�I��������
		
	}

	//�p��U�s���ܲ���
	for(i = start; i <= end; i++){
		sum = 0;
		sum_of_number = 0;
		for(j = start; j <= i; j++){
			sum = sum + (j - meanFore[i]) * (j - meanFore[i]) * number[j];	//(�Ƕ��� - �e��������).^2 * �ӦǶ��Ȳέp����
			sum_of_number = sum_of_number + number[j];						//�e���`�ֿn����
		}
		varFore[i] = sum / sum_of_number;	//�e���ܲ���

		sum = 0;
		sum_of_number = 0;
		for(k = j; k <= end; k++){
			sum = sum + (k - meanBack[i]) * (k - meanBack[i]) * number[k];	//(�Ƕ��� - �I��������).^2 * �ӦǶ��Ȳέp����
			sum_of_number = sum_of_number + number[k];						//�I���`�ֿn����
		}
		varBack[i] = sum / sum_of_number;	//�I���ܲ���
		
	}
	

	//�p��s���ܲ����[�v�`�X
	for(i = start; i < end; i++){
		sum_of_var[i] = accu_probabilityFore[i] * varFore[i] +  accu_probabilityBack[i] * varBack[i];	//�e���ֿn���v * �e���ܲ��� + �I���ֿn���v * �I���ܲ���
	}
	
	//��X�ܲ��[�v�`�X���̤p��
	min = 99999999;
	for(i = start; i < end; i++){
		if(sum_of_var[i] < min){
			min = sum_of_var[i];
			outValue = i;	//i �N�O�ڭ̭n���֭�
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
	int time;	//�����Ʀ���

	int *hist = (int *)malloc(sizeof(int) * 256);		//�����ϲέp��l���
	int *hist_in = (int *)malloc(sizeof(int) * 256);	//�����ϥ��ƫe
	int *hist_out = (int *)malloc(sizeof(int) * 256);	//�����ϥ��ƫ�


	//��l��
	for(i=0; i < 256; i++){
		hist[i] = 0;
		hist_in[i] = 0;
		hist_out[i] = 0;
	}

	//�Ƕ��Ȳέp
	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			hist[ ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] ] ++;
		}
	}

	for(i = 0; i < 256; i++){
		hist_in[i] = hist[i];
		hist_out[i] = hist[i];
	}


	/*********�Х��]�w���Ʀ���*************/
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

void thinning(IplImage *read_image, IplImage *out_image)	//�G�ȼv���ӽu��
{
	int flag = 1;
	int i, j, k, n;
	int up, down, left, right;
	int p[9];	//�ϧ�:1 �I��:2 �I�����:-1;
	int TMP = 128;	//�I���Ըɪ��@�׼ȩw��

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

			

				/* ����1:  */
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

	int **tmp_matrix = (int **)malloc(sizeof(int *) * read_image->height);	//���o�G���}�C�ʺA�O����;
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

	mat2grat(tmp_matrix, out_image);	//�x�}��� �� �Ƕ��v��

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

	int **tmp_matrix = (int **)malloc(sizeof(int *) * read_image->height);	//���o�G���}�C�ʺA�O����;
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

	//print_int_mat_txt(tmp_matrix, read_image->height, read_image->width, "lap.txt");	//�L�X���matrix

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

void dilation(IplImage *read_image, IplImage *out_image)	//����
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

void erosion(IplImage *read_image, IplImage *out_image)	//�I�k
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

				//imshow(out_image, "out");	//��ܼv��
	
			}
		}
	}

	return;

}

void imopen(IplImage *read_image, IplImage *out_image)	//�}�B��
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

	erosion(read_image, erosion_image);	//�I�k
	dilation(erosion_image, out_image);	//����

	cvReleaseImage(&erosion_image );	//����Ϲ��O����

	return;

}

void regfill(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	int flag;

	IplImage *current_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//��e��
	IplImage *last_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);		//�L�h��
	IplImage *edge_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);		//��ɹ�
	IplImage *invedge_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//��ɤϦV��
	IplImage *erosion_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//�I�k��
	IplImage *dilation_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//���ȹ�

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


	erosion(read_image, erosion_image);	//�I�k
	
	//�D��ɹ�
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(edge_image->imageData + (i)*edge_image->widthStep))[j] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j] - ((uchar *)(erosion_image->imageData + (i)*erosion_image->widthStep))[j];
		}
	}
	
	
	//�D��ɤϦV��
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){

			if(((uchar *)(edge_image->imageData + (i)*edge_image->widthStep))[j] == 255)
				((uchar *)(invedge_image->imageData + (i)*invedge_image->widthStep))[j] = 0;
			else
				((uchar *)(invedge_image->imageData + (i)*invedge_image->widthStep))[j] = 255;
		}
	}
	

	flag = 1;
	//�M�w�n��R���ϰ�
	((uchar *)(last_image->imageData + (26)*last_image->widthStep))[20] = 255;
	((uchar *)(last_image->imageData + (39)*last_image->widthStep))[93] = 255;
	((uchar *)(last_image->imageData + (166)*last_image->widthStep))[165] = 255;
	
	while(flag != 0){
		flag = 0;
		dilation(last_image, dilation_image);
		
		//���ȹ� �M ��ɤϦV�� ���涰�A���G���u��e�ϡv
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				
				if(((uchar *)(dilation_image->imageData + (i)*dilation_image->widthStep))[j] == 255 && ((uchar *)(invedge_image->imageData + (i)*invedge_image->widthStep))[j] == 255)
					((uchar *)(current_image->imageData + (i)*current_image->widthStep))[j] = 255;
				else
					((uchar *)(current_image->imageData + (i)*current_image->widthStep))[j] = 0;
			}
		}

		

		//����e�� �M �L�h�� �O�_�@��
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				if(((uchar *)(current_image->imageData + (i)*current_image->widthStep))[j] == ((uchar *)(last_image->imageData + (i)*last_image->widthStep))[j])
					flag++;
			}
		}

		//���e��copy���L�h��
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				((uchar *)(last_image->imageData + (i)*last_image->widthStep))[j] = ((uchar *)(current_image->imageData + (i)*current_image->widthStep))[j];
			}
		}


		//�p�Gflag ���� �Ҧ�pixel�ơA�N��w�g��R�����A�Nflag�]��0�A���X�j��
		if(flag == (read_image->height * read_image->width))
			flag = 0;
	}

	//�N��e�Ͽ�X
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			if(((uchar *)(current_image->imageData + (i)*current_image->widthStep))[j] == 255 || ((uchar *)(edge_image->imageData + (i)*edge_image->widthStep))[j] == 255)
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
		}
	}

	cvReleaseImage(&current_image );	//����Ϲ��O����
	cvReleaseImage(&last_image );	//����Ϲ��O����
	cvReleaseImage(&edge_image );	//����Ϲ��O����
	cvReleaseImage(&invedge_image );	//����Ϲ��O����
	cvReleaseImage(&erosion_image );	//����Ϲ��O����
	cvReleaseImage(&dilation_image );	//����Ϲ��O����

	return;
}


void skeleton(IplImage *read_image, IplImage *out_image)
{
	int i,j;
	int flag;


	IplImage *wait_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//�ݳB�z��
	IplImage *erosion_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//�I�k��
	IplImage *open_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//�}�B���
	IplImage *diffset_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//���X�t��

	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(wait_image->imageData + (i)*wait_image->widthStep))[j] = 0;
			((uchar *)(erosion_image->imageData + (i)*erosion_image->widthStep))[j] = 0;
			((uchar *)(open_image->imageData + (i)*open_image->widthStep))[j] = 0;
			((uchar *)(diffset_image->imageData + (i)*diffset_image->widthStep))[j] = 0;
			((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 0;
		}
	}

	//���copy���ݳB�z��
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			((uchar *)(wait_image->imageData + (i)*wait_image->widthStep))[j] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
		}
	}

	flag = 1;
	while(flag != 0){
		flag = 0;

		erosion(wait_image, erosion_image);	//�I�k
		imopen(erosion_image, open_image);	//�}�B��

		//���X�t�� = �I�k�� - �}�B���
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				((uchar *)(diffset_image->imageData + (i)*diffset_image->widthStep))[j] = ((uchar *)(erosion_image->imageData + (i)*erosion_image->widthStep))[j] - ((uchar *)(open_image->imageData + (i)*open_image->widthStep))[j];
			}
		}

		//��X�� = ��X�� �p�� ���X�t��
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				if(((uchar *)(diffset_image->imageData + (i)*diffset_image->widthStep))[j] == 255 || ((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] == 255)
					((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 255;
			}
		}

		//����}�B��O�_���Ŷ��X�Cflag != 0 : ���O�Ŷ��X �A flag == 0 : ���Ŷ��X
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				if(((uchar *)(open_image->imageData + (i)*open_image->widthStep))[j] == 255)
					flag++;
			}
		}

		//�I�k��copy���ݳB�z��
		for(i=0; i < read_image->height; i++){
			for(j = 0; j < read_image->width; j++){
				((uchar *)(wait_image->imageData + (i)*wait_image->widthStep))[j] = ((uchar *)(erosion_image->imageData + (i)*erosion_image->widthStep))[j];
			}
		}
	}


	cvReleaseImage(&wait_image );	//����Ϲ��O����
	cvReleaseImage(&erosion_image );	//����Ϲ��O����
	cvReleaseImage(&open_image );	//����Ϲ��O����
	cvReleaseImage(&diffset_image );	//����Ϲ��O����

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

	IplImage *label_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//�аO��
	
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

					if(i == 0 && j == 0){	//��1��pixel�A�������s�аO
						((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
						label[labelCnt] = labelCnt;
						labelCnt ++;
					}else if(i == 0 && j != 0){	//��0�C�A�u�ˬd���䪬�A
						left = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j-1];
						if(left == 255){	//���䬰�e���A���M����@��
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1];
						}else{	//���䬰�I���A���s�аO
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
							label[labelCnt] = labelCnt;
							labelCnt ++;
						}
							
					}else if(i != 0 && j == 0){		//��0��A�ˬd�W��
						up = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j];
						if(up == 255){	//�W���O�e���A�M�W���@��
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j];
						}else{			//�W���O�I���A���s�аO
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
							label[labelCnt] = labelCnt;
							labelCnt ++;
						}
					}else{		//�D�����I
						up = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j];
						left = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j-1];
						if(up == 0 && left == 0){	//�W���B���䳣�O�I���A���s�аO
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
							label[labelCnt] = labelCnt;
							labelCnt ++;
						}else{		//�W���B���� �ܤ֤@�ӬO�e��
							labelMin = 99999;

							//��X�̤p����
							if(up == 255)
								if(((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] <= labelMin)			
									labelMin = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j];			//�̤p���� = �W��
							if(left == 255)
								if(((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1] <= labelMin)		
									labelMin = ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1];			//�̤p���� = ����

							//���O�Ҧ��аO�ۦP
							if(label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] ] = label[ labelMin ];		//�W�� = �̤p����
							if(label[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1] ] = label[ labelMin ];		//���� = �̤p����

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
					
					if(i == 0 && j == 0){	//��1��pixel�A�������s�аO
						((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
						label[labelCnt] = labelCnt;
						labelCnt ++;
					}else if(i == 0 && j != 0){	//��0�C�A�u�ˬd���䪬�A
						left = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j-1];
						if(left == 255){	//���䬰�e���A���M����@��
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1];
						}else{	//���䬰�I���A���s�аO
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
							label[labelCnt] = labelCnt;
							labelCnt ++;
						}	
					}else if(i != 0 && j == 0){		//��0��A�ˬd�W���B�k�W
						up = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j];
						rightUp = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j + 1];
						if(up == 0 && rightUp){		//�W���B�k�W�O�I���A���s�аO
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
							label[labelCnt] = labelCnt;
							labelCnt ++;
						}else{	//�W���B�k�W �ܤ֦��@�ӬO�e��
							labelMin = 99999;
							//��X�̤p����
							if(up == 255)
								if(((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] <= labelMin)			
									labelMin = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j];			//�̤p���� = �W��
							if(rightUp == 255)
								if(((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1] <= labelMin)			
									labelMin = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1];			//�̤p���� = �k�W

							//���O�Ҧ��аO�ۦP
							if(label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] ] = label[ labelMin ];		//�W�� = �̤p�аO
							if(label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1] ] = label[ labelMin ];		//�k�W = �̤p�аO

							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelMin;
						}
					}else{		//�D�����I
						rightUp = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j + 1];
						up = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j];
						leftUp = ((uchar *)(read_image->imageData + (i - 1)*read_image->widthStep))[j - 1];
						left = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j-1];
						if(rightUp == 0 && up == 0 && leftUp == 0 && left == 0){	//�k�W�B�W���B���W�B���䳣�O�I���A���s�аO
							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelCnt;
							label[labelCnt] = labelCnt;
							labelCnt ++;
						}else{		//�k�W�B�W���B���W�B���� �ܤ֤@�ӬO�e��
							labelMin = 99999;

							//��X�̤p����
							if(rightUp == 255)
								if(((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1] <= labelMin)			
									labelMin = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1];		//�̤p���� = �k�W
							if(up == 255)
								if(((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] <= labelMin)			
									labelMin = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j];			//�̤p���� = �W��
							if(leftUp == 255)
								if(((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j - 1] <= labelMin)			
									labelMin = ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j - 1];		//�̤p���� = ���W
							if(left == 255)
								if(((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1] <= labelMin)		
									labelMin = ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1];			//�̤p���� = ����

							//���O�Ҧ��аO�ۦP
							if(label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j + 1] ] = label[ labelMin ];	//�k�W = �̤p�аO
							if(label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j] ] = label[ labelMin ];		//�W�� = �̤p�аO
							if(label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j - 1] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i - 1)*label_image->widthStep))[j - 1] ] = label[ labelMin ];	//���W = �̤p�аO
							if(label[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1] ] != 0)		
								label[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j - 1] ] = label[ labelMin ];		//���� = �̤p�аO
							

							((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] = labelMin;
						}
					}
				}
			}
		}
	}

	

	//�N�ۦP���O�аO�A�令�P�˼аO
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

	
	

	//�p�⪫��Ӽ�
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

	//�p�⪫�󭱿n(pixel��)
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			area[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] ] ++;
		}
	}

	//�h�����T
	for(i=0; i < read_image->height; i++){
		for(j = 0; j < read_image->width; j++){
			if(area[ ((uchar *)(label_image->imageData + (i)*label_image->widthStep))[j] ] < 500)
				((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = 0;
			else
				((uchar *)(out_image->imageData + (i)*out_image->widthStep))[j] = ((uchar *)(read_image->imageData + (i)*read_image->widthStep))[j];
		}
	}

	if(!cvSaveImage("out_image.bmp",out_image)) printf("Could not save: %s\n", "out_image.bmp");

	cvReleaseImage(&label_image );	//����Ϲ��O����
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

void project_statistics(IplImage *read_image)	//��ROI�᪺��ޤ��ߩw��
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

	int *hist = (int *)malloc(sizeof(int) * 256);		//�����ϲέp��l���

	//��l��
	for(i=0; i < 256; i++){
		hist[i] = 0;
	}

	//�Ƕ��Ȳέp
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