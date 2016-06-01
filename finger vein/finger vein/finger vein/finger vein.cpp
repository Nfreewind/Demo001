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
	IplImage *read_image = cvLoadImage("q03Hroip2.jpg" , -1);	//�q���Ū�J�Ϲ�

	cvNamedWindow("showpicture1",0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture1"
	cvResizeWindow("showpicture1", read_image->width, read_image->height);	//�]�w�����j�p
	cvNamedWindow("showpicture2",0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture2"
	cvResizeWindow("showpicture2", read_image->width, read_image->height);	//�]�w�����j�p
	cvNamedWindow("showpicture3",0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture3"
	cvResizeWindow("showpicture3", read_image->width, read_image->height);	//�]�w�����j�p
	cvNamedWindow("showpicture4",0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture4"
	cvResizeWindow("showpicture4", read_image->width, read_image->height);	//�]�w�����j�p
	/*cvNamedWindow("showpicture5",0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture5"
	cvResizeWindow("showpicture5", read_image->width, read_image->height);	//�]�w�����j�p
	cvNamedWindow("showpicture6",0);	//�Ыؤ@�ӷs�����A�R�W��"showpicture6"
	cvResizeWindow("showpicture6", read_image->width, read_image->height);	//�]�w�����j�p*/

	IplImage *gray_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//���t�O���鵹tmp_image	��q�D																							//��q�D�줸��
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
	
	
	

	/*complex **read_DFT = (complex **)malloc(sizeof(complex *) * in_no);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < in_no; i++){
		read_DFT[i] = (complex *)malloc(sizeof(complex) * in_no);
	}

	FILE *filterPtr = fopen("butter_lowpass_D15N1_256.txt","r");
	float **filter = (float **)malloc(sizeof(float *) * in_no);
	for(i=0; i < in_no; i++){
		filter[i] = (float *)malloc(sizeof(float) * in_no);
	}

	complex **read_low_pass = (complex **)malloc(sizeof(complex *) * in_no);	//���o�G���}�C�ʺA�O����;
	for(i=0; i < in_no; i++){
		read_low_pass[i] = (complex *)malloc(sizeof(complex) * in_no);
	}*/
	



	/************************************************************************************************/

	
	
	gray(read_image, gray_image);	//�ন��q�D�Ƕ�
	
	
	

	;
	;

	/*
	//binarization(tmp_image, out_image);	//�G�Ȥ�

	//bit_plane(tmp_image, out_image, 4);	//�줸����

	//im_resize(tmp_image, out_image, 16);	//�Ŷ��ѪR��

	//quantizatioin(tmp_image, out_image, 4);	//���öq��

	//dither(tmp_image, out_image);	//�V��

	//multi_dither(tmp_image, out_image, 9);	//�q�Ʀh�h�V��

	//error_diff(tmp_image, out_image);	//�~�t�X���k

	//im_add(tmp_image, out_image, 128);	//�v���[�k
	//im_sub(tmp_image, out_image, 128);	//�v����k
	//im_multiply(tmp_image, out_image, 2);	//�v�����k
	//im_divide(tmp_image, out_image, 2);		//�v�����k

	//im_complement(tmp_image, out_image);	//�v���ɦ�
	//calc_Hist(tmp_image);	//����ϲέp
	//hist_expand(tmp_image);	//������X�i�k
	//hist_eq(tmp_image);	//����ϵ���
	//LUT_hist_eq(tmp_image, out_image);	////LUT����ϵ���
	//avg_filter(tmp_image, out_image, 1, 25);	//�����o�i��
	//lapla_filter(tmp_image, out_image, 3);	//Laplacian�o�i��
	//LoG_filter(tmp_image, out_image, 5);	//LoG�o�i��
	//HPF_filter(tmp_image, out_image, 3);	//���q�o�i��
	//unsharp(tmp_image, out_image, 4, (float)1/3);	//�h�U�Q��
	//hb_filter(tmp_image, out_image, (float) 5/6);	//���W�T�o�i��
	//max_min_filter(tmp_image, out_image);	//�̤j�o�i��
	//mid_filter(tmp_image, out_image);	//�����o�i��
	//ROI_foiter(tmp_image, out_image, 0, 255, 0, 255, 4);	//ROI�o�i

	//shift_array(in_no, &data1[0], &sh_data1[0]);	//�@���ƦC����
	//DFT1(in_no, &sh_data1[0], &sh_data1_DFT[0]);	//�@���ť߸��ഫ
	//print_complex(in_no, &sh_data1_DFT[0]);	//�L�X�Ƽư}�C

	//conv(M, &data1[0], N, &data2[0], &conv_out[0]);	//�ۿn
	//print_complex(M+N-1, &conv_out[0]);	//�L�X�Ƽư}�C

	//cconv(in_no, &data1[0], &data2[0], &cconv_out[0]);	//���αۿn
	//print_complex(in_no, &cconv_out[0]);	//�L�X�Ƽư}�C

	

	
	


	
	DFT2(in_no, tmp_image, read_DFT);	//�G���ť߸��ഫ
	fftshow(in_no, read_DFT);	//����W�й�
	load_filter(filterPtr, filter);	//Ū���o�i��
	

	for(i=0;i < in_no; i++){
		for(j=0; j < in_no; j++){
			read_low_pass[i][j].real = read_DFT[i][j].real * filter[i][j];
			read_low_pass[i][j].image = read_DFT[i][j].image * filter[i][j];
		}
	}
	fftshow(in_no, read_low_pass);	//����W�й�

	IDFT2(in_no, read_low_pass, out_image);	//�G���ť߸��u�ϡv�ഫ
	
	


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

	//DFT1();	//�@�������ť߸��ഫ
	//inv_DFT1();	//�@�������ϳť߸��ഫ
	//conv();	//�ۿn
	//cconv();	//���αۿn

	


	//avg_filter(tmp_image, out_image, 7);	//�����o�i��
	//mid_filter(tmp_image, out_image, 5);	//�����o�i��
	//outliers_filter(tmp_image, out_image, 5);	//�[���I�o�i��
	//adaptive_filter(tmp_image, out_image, 7);	//�i�A���o�i
	//band_reject_folter(tmp_image, out_image);	//�a���o�i
	//notch_filter(tmp_image, out_image);	//���i�o�i
	//remove_butter(tmp_image, out_image);	//�h��butter���T
	//remove_motion(tmp_image, out_image);	//�h���ʺA�ҽk
	//wiener_filter(tmp_image, out_image);	//wiener�o�i�h��butter���T
	//Tn = Otsu(tmp_image);	//Otsu�t��k�A��X�̨λ֭Ȥ�k
 	//binarization(tmp_image, out_image, 96);	//�G�Ȥ�

	//histprint(tmp_image);	//�H�ƭȻP�r���⪽���ϦC�L�X��
	//fprintf_hist(tmp_image, "hist.txt");	//�έp�v���Ƕ��ơA��X��TXT��

	//difference_filter(tmp_image, out_image);	//���(�q�`���t��)
	//roberts_filter(tmp_image, out_image);	//���(Rober)
	//sobel_filter(tmp_image, out_image);	//���(Sobel)
	//prewitt(tmp_image, grad_image);	//��prewitt�^������
	//Tn = Otsu(grad_image);	//Otsu�t��k�A��X�̨λ֭Ȥ�k
	//binarization(grad_image, binary_image, Tn);	//�G�Ȥ�
	//thinning(binary_image, out_image);	//�G�ȼv���ӽu��

	//laplacian_filater(gray_image, lap_image);	//laplacian�G���L��
	//thinning(binary_image, out_image);	//�G�ȼv���ӽu��
	//zero_cross(gray_image, zero_image);	//�Q�ιs����I�A�^������
	//Tn = Otsu(gray_image);	//Otsu�t��k�A��X�̨λ֭Ȥ�k
	//binarization(gray_image, binary_image,  110);	//�G�Ȥ�

	//binarization(gray_image, binary_image,  110);	//�G�Ȥ�
	
	//skeleton(gray_image, out_image);	//���[��

	LUT_hist_eq(gray_image, out_image);	////LUT����ϵ���*/
	
	;
	;
	
	retinex(gray_image, enhance_image);
	mid_filter(enhance_image, mid_image, 3);	//�����o�i��

	;
	;
	//curve_line(mid_image, 0);	//��������u
	//smooth_line(gray_image, 0);	//��Ͼ�������u(���Ƥ�)

	;
	;

	segmentation(gray_image, mid_image, segment_image);	//���R�ߤ���
	erosion(segment_image, erosion_image);	//�I�k
	dilation(erosion_image, dilation_image);	//����
	dilation(dilation_image, erosion_image);	//����
	erosion(erosion_image, compensate_image);	//�I�k
	
	;
	;
	//mid_filter(compensate_image, out_image, 3);	//�����o�i��
	
	//var_line(read_image, 99);	//�ܲ��q���u
	
	//find_minpoint(avg_image, out_image);	//���ޤ����I

	;
	;

	//project_statistics(gray_image);	//��ROI�᪺��ޤ��ߩw��
	//Tn = Otsu(avg_image);	//Otsu�t��k�A��X�̨λ֭Ȥ�k
	//binarization(gray_image, enhance_image, 250);	//�G�Ȥ�
	//binarization2(gray_image, enhance_image, 250, 255);	//�����֭ȤG�Ȥ�
	
	//label(compensate_image, label_image,8);
	//thinning(enhance_image, thin_image);	//�G�ȼv���ӽu��
	//bound_mix_thin(thin_image, bound_image, out_image);	//��ӽu�ƩM��ɹ� �|�X

	


	
	
	
	

	
	
	;
	if(!cvSaveImage("q03Hroip2_seg.bmp",segment_image)) printf("Could not save: %s\n", "q03Hroip2_seg.bmp");

	
	cvShowImage("showpicture1", gray_image);	//��ܹϹ�
	cvShowImage("showpicture2", mid_image);	//��ܹϹ�
	cvShowImage("showpicture3", segment_image);	//��ܹϹ�
	cvShowImage("showpicture4", compensate_image);	//��ܹϹ�
	//cvShowImage("showpicture6", second_binary_image);	//��ܹϹ�

	cvWaitKey(0);

	
	cvDestroyWindow("showpicture1");	//��������
	cvDestroyWindow("showpicture2");	//��������
	cvDestroyWindow("showpicture3");	//��������
	cvDestroyWindow("showpicture4");	//��������
	//cvDestroyWindow("showpicture5");	//��������
	//cvDestroyWindow("showpicture6");	//��������
	cvReleaseImage(&read_image );	//����Ϲ��O����
	cvReleaseImage(&gray_image);		//����Ϲ��O����
	cvReleaseImage(&enhance_image);		//����Ϲ��O����
	cvReleaseImage(&mid_image);		//����Ϲ��O����
	cvReleaseImage(&avg_image);	//����Ϲ��O����
	cvReleaseImage(&segment_image);	//����Ϲ��O����
	cvReleaseImage(&erosion_image);		//����Ϲ��O����
	cvReleaseImage(&dilation_image);		//����Ϲ��O����
	cvReleaseImage(&compensate_image);		//����Ϲ��O����
	cvReleaseImage(&bound_image);		//����Ϲ��O����
	cvReleaseImage(&thin_image);		//����Ϲ��O����
	cvReleaseImage(&out_image);		//����Ϲ��O����

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

