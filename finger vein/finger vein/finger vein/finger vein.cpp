// openpicture.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <cv.h>
#include <highgui.h>

#include <math.h>
#include "constant.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "fuzzy_function.h"
#include "fuzzy_retinex.h"
#include "retinex.h"
//#include "struct.h"

#include "function.h"

int _tmain(int argc, _TCHAR* argv[])
{
	int i,j;
	int Tn;
	IplImage *read_image = cvLoadImage("q03Hroip2.jpg" , -1);	//從文件中讀入圖像

	cvNamedWindow("showpicture1",0);	//創建一個新視窗，命名為"showpicture1"
	cvResizeWindow("showpicture1", read_image->width, read_image->height);	//設定視窗大小
	cvNamedWindow("showpicture2",0);	//創建一個新視窗，命名為"showpicture2"
	cvResizeWindow("showpicture2", read_image->width, read_image->height);	//設定視窗大小
	cvNamedWindow("showpicture3",0);	//創建一個新視窗，命名為"showpicture3"
	cvResizeWindow("showpicture3", read_image->width, read_image->height);	//設定視窗大小
	cvNamedWindow("showpicture4",0);	//創建一個新視窗，命名為"showpicture4"
	cvResizeWindow("showpicture4", read_image->width, read_image->height);	//設定視窗大小
	/*cvNamedWindow("showpicture5",0);	//創建一個新視窗，命名為"showpicture5"
	cvResizeWindow("showpicture5", read_image->width, read_image->height);	//設定視窗大小
	cvNamedWindow("showpicture6",0);	//創建一個新視窗，命名為"showpicture6"
	cvResizeWindow("showpicture6", read_image->width, read_image->height);	//設定視窗大小*/

	IplImage *gray_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//分配記憶體給tmp_image	單通道																							//單通道位元組
	IplImage *enhance_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	IplImage *mid_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	IplImage *avg_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	IplImage *segment_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	IplImage *erosion_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	IplImage *dilation_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	IplImage *compensate_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	IplImage *label_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	IplImage *bound_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	IplImage *thin_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	IplImage *out_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);
	
	
	

	/*complex **read_DFT = (complex **)malloc(sizeof(complex *) * in_no);	//取得二維陣列動態記憶體;
	for(i=0; i < in_no; i++){
		read_DFT[i] = (complex *)malloc(sizeof(complex) * in_no);
	}

	FILE *filterPtr = fopen("butter_lowpass_D15N1_256.txt","r");
	float **filter = (float **)malloc(sizeof(float *) * in_no);
	for(i=0; i < in_no; i++){
		filter[i] = (float *)malloc(sizeof(float) * in_no);
	}

	complex **read_low_pass = (complex **)malloc(sizeof(complex *) * in_no);	//取得二維陣列動態記憶體;
	for(i=0; i < in_no; i++){
		read_low_pass[i] = (complex *)malloc(sizeof(complex) * in_no);
	}*/
	



	/************************************************************************************************/

	
	
	gray(read_image, gray_image);	//轉成單通道灰階
	
	
	

	;
	;

	/*
	//binarization(tmp_image, out_image);	//二值化

	//bit_plane(tmp_image, out_image, 4);	//位元平面

	//im_resize(tmp_image, out_image, 16);	//空間解析度

	//quantizatioin(tmp_image, out_image, 4);	//均勻量化

	//dither(tmp_image, out_image);	//混色

	//multi_dither(tmp_image, out_image, 9);	//量化多層混色

	//error_diff(tmp_image, out_image);	//誤差擴散法

	//im_add(tmp_image, out_image, 128);	//影像加法
	//im_sub(tmp_image, out_image, 128);	//影像減法
	//im_multiply(tmp_image, out_image, 2);	//影像乘法
	//im_divide(tmp_image, out_image, 2);		//影像除法

	//im_complement(tmp_image, out_image);	//影像補色
	//calc_Hist(tmp_image);	//直方圖統計
	//hist_expand(tmp_image);	//直方圖擴展法
	//hist_eq(tmp_image);	//直方圖等化
	//LUT_hist_eq(tmp_image, out_image);	////LUT直方圖等化
	//avg_filter(tmp_image, out_image, 1, 25);	//平均濾波器
	//lapla_filter(tmp_image, out_image, 3);	//Laplacian濾波器
	//LoG_filter(tmp_image, out_image, 5);	//LoG濾波器
	//HPF_filter(tmp_image, out_image, 3);	//高通濾波器
	//unsharp(tmp_image, out_image, 4, (float)1/3);	//去銳利化
	//hb_filter(tmp_image, out_image, (float) 5/6);	//高增幅濾波器
	//max_min_filter(tmp_image, out_image);	//最大濾波器
	//mid_filter(tmp_image, out_image);	//中位濾波器
	//ROI_foiter(tmp_image, out_image, 0, 255, 0, 255, 4);	//ROI濾波

	//shift_array(in_no, &data1[0], &sh_data1[0]);	//一維數列平移
	//DFT1(in_no, &sh_data1[0], &sh_data1_DFT[0]);	//一維傅立葉轉換
	//print_complex(in_no, &sh_data1_DFT[0]);	//印出複數陣列

	//conv(M, &data1[0], N, &data2[0], &conv_out[0]);	//旋積
	//print_complex(M+N-1, &conv_out[0]);	//印出複數陣列

	//cconv(in_no, &data1[0], &data2[0], &cconv_out[0]);	//環形旋積
	//print_complex(in_no, &cconv_out[0]);	//印出複數陣列

	

	
	


	
	DFT2(in_no, tmp_image, read_DFT);	//二維傅立葉轉換
	fftshow(in_no, read_DFT);	//顯示頻譜圖
	load_filter(filterPtr, filter);	//讀取濾波器
	

	for(i=0;i < in_no; i++){
		for(j=0; j < in_no; j++){
			read_low_pass[i][j].real = read_DFT[i][j].real * filter[i][j];
			read_low_pass[i][j].image = read_DFT[i][j].image * filter[i][j];
		}
	}
	fftshow(in_no, read_low_pass);	//顯示頻譜圖

	IDFT2(in_no, read_low_pass, out_image);	//二維傅立葉「反」轉換
	
	


	int data_no = 8;
	complex *in_data[data_no];
	complex *out_data[data_no];
	in_data = (complex *)malloc(sizeof(complex) * data_no);
	out_data = (complex *)malloc(sizeof(complex) * data_no);

	in_data[0]->real = 2;	in_data[0]->image = 0;
	in_data[1]->real = 3;	in_data[1]->image = 0;
	in_data[2]->real = 4;	in_data[2]->image = 0;
	in_data[3]->real = 5;	in_data[3]->image = 0;
	in_data[4]->real = 6;	in_data[4]->image = 0;
	in_data[5]->real = 7;	in_data[5]->image = 0;
	in_data[6]->real = 8;	in_data[6]->image = 0;
	in_data[7]->real = 1;	in_data[7]->image = 0;

	//DFT1();	//一維離散傅立葉轉換
	//inv_DFT1();	//一維離散反傅立葉轉換
	//conv();	//旋積
	//cconv();	//環形旋積

	


	//avg_filter(tmp_image, out_image, 7);	//平均濾波器
	//mid_filter(tmp_image, out_image, 5);	//中位濾波器
	//outliers_filter(tmp_image, out_image, 5);	//歧異點濾波器
	//adaptive_filter(tmp_image, out_image, 7);	//可適性濾波
	//band_reject_folter(tmp_image, out_image);	//帶拒濾波
	//notch_filter(tmp_image, out_image);	//陷波濾波
	//remove_butter(tmp_image, out_image);	//去除butter雜訊
	//remove_motion(tmp_image, out_image);	//去除動態模糊
	//wiener_filter(tmp_image, out_image);	//wiener濾波去除butter雜訊
	//Tn = Otsu(tmp_image);	//Otsu演算法，找出最佳閥值方法
 	//binarization(tmp_image, out_image, 96);	//二值化

	//histprint(tmp_image);	//以數值與字母把直條圖列印出來
	//fprintf_hist(tmp_image, "hist.txt");	//統計影像灰階數，輸出至TXT檔

	//difference_filter(tmp_image, out_image);	//梯度(通常的差分)
	//roberts_filter(tmp_image, out_image);	//梯度(Rober)
	//sobel_filter(tmp_image, out_image);	//梯度(Sobel)
	//prewitt(tmp_image, grad_image);	//用prewitt擷取輪廓
	//Tn = Otsu(grad_image);	//Otsu演算法，找出最佳閥值方法
	//binarization(grad_image, binary_image, Tn);	//二值化
	//thinning(binary_image, out_image);	//二值影像細線化

	//laplacian_filater(gray_image, lap_image);	//laplacian二階微分
	//thinning(binary_image, out_image);	//二值影像細線化
	//zero_cross(gray_image, zero_image);	//利用零交錯點，擷取輪廓
	//Tn = Otsu(gray_image);	//Otsu演算法，找出最佳閥值方法
	//binarization(gray_image, binary_image,  110);	//二值化

	//binarization(gray_image, binary_image,  110);	//二值化
	
	//skeleton(gray_image, out_image);	//骨架化

	LUT_hist_eq(gray_image, out_image);	////LUT直方圖等化*/
	
	;
	;
	
	retinex(gray_image, enhance_image);
	mid_filter(enhance_image, mid_image, 3);	//中位濾波器

	;
	;
	//curve_line(mid_image, 0);	//橫切面曲線
	//smooth_line(gray_image, 0);	//原圖橫切面曲線(平滑化)

	;
	;

	segmentation(gray_image, mid_image, segment_image);	//指靜脈分割
	erosion(segment_image, erosion_image);	//侵蝕
	dilation(erosion_image, dilation_image);	//膨脹
	dilation(dilation_image, erosion_image);	//膨脹
	erosion(erosion_image, compensate_image);	//侵蝕
	
	;
	;
	//mid_filter(compensate_image, out_image, 3);	//中位濾波器
	
	//var_line(read_image, 99);	//變異量曲線
	
	//find_minpoint(avg_image, out_image);	//找血管中心點

	;
	;

	//project_statistics(gray_image);	//取ROI後的血管中心定位
	//Tn = Otsu(avg_image);	//Otsu演算法，找出最佳閥值方法
	//binarization(gray_image, enhance_image, 250);	//二值化
	//binarization2(gray_image, enhance_image, 250, 255);	//雙重閥值二值化
	
	//label(compensate_image, label_image,8);
	//thinning(enhance_image, thin_image);	//二值影像細線化
	//bound_mix_thin(thin_image, bound_image, out_image);	//把細線化和邊界圖 疊合

	


	
	
	
	

	
	
	;
	if(!cvSaveImage("q03Hroip2_seg.bmp",segment_image)) printf("Could not save: %s\n", "q03Hroip2_seg.bmp");

	
	cvShowImage("showpicture1", gray_image);	//顯示圖像
	cvShowImage("showpicture2", mid_image);	//顯示圖像
	cvShowImage("showpicture3", segment_image);	//顯示圖像
	cvShowImage("showpicture4", compensate_image);	//顯示圖像
	//cvShowImage("showpicture6", second_binary_image);	//顯示圖像

	cvWaitKey(0);

	
	cvDestroyWindow("showpicture1");	//關閉視窗
	cvDestroyWindow("showpicture2");	//關閉視窗
	cvDestroyWindow("showpicture3");	//關閉視窗
	cvDestroyWindow("showpicture4");	//關閉視窗
	//cvDestroyWindow("showpicture5");	//關閉視窗
	//cvDestroyWindow("showpicture6");	//關閉視窗
	cvReleaseImage(&read_image );	//釋放圖像記憶體
	cvReleaseImage(&gray_image);		//釋放圖像記憶體
	cvReleaseImage(&enhance_image);		//釋放圖像記憶體
	cvReleaseImage(&mid_image);		//釋放圖像記憶體
	cvReleaseImage(&avg_image);	//釋放圖像記憶體
	cvReleaseImage(&segment_image);	//釋放圖像記憶體
	cvReleaseImage(&erosion_image);		//釋放圖像記憶體
	cvReleaseImage(&dilation_image);		//釋放圖像記憶體
	cvReleaseImage(&compensate_image);		//釋放圖像記憶體
	cvReleaseImage(&bound_image);		//釋放圖像記憶體
	cvReleaseImage(&thin_image);		//釋放圖像記憶體
	cvReleaseImage(&out_image);		//釋放圖像記憶體

	/*for (i=0; i < in_no; i++) free(read_DFT[i]);
	free(read_DFT);

	for (i=0; i < in_no; i++) free(filter[i]);
	free(filter);

	for (i=0; i < in_no; i++) free(read_low_pass[i]);
	free(read_low_pass);*/


	
    
	

    return 0;
}

/*
	FILE *f_outR = fopen("outR.txt","w");
	for(i=0; i < in_no; i++){
		for(j=0; j < in_no; j++){
			fprintf(f_outR,"%f ",read_DFT[i][j].real);
		}
		fprintf(f_outR,"\n",read_DFT[i][j-1].real);
	}
	fclose(f_outR);

	
	FILE *f_outI = fopen("outI.txt","w");
	for(i=0; i < in_no; i++){
		for(j=0; j < in_no; j++){
			fprintf(f_outI,"%f ",read_DFT[i][j].image);
		}
		fprintf(f_outI,"\n",read_DFT[i][j-1].image);
	}
	fclose(f_outI);
*/


/*
	FILE *f_outR = fopen("out.txt","w");
	for(i=0; i < in_no; i++){
		for(j=0; j < in_no; j++){
			fprintf(f_outR,"%f ",idealFilter[i][j]);
		}
		fprintf(f_outR,"\n",idealFilter[i][j-1]);
	}
	fclose(f_outR);


*/

