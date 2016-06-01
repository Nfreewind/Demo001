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

void retinex(IplImage *read_image, IplImage *out_image)
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
			fsum += (float)((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
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
				relation = (float)((uchar *)(read_image->imageData + i*read_image->widthStep))[j + widthStep] - (float)((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
				membership_tmp[i][j] = membership_input[i][j] - relation;
				membership_tmp[i][j + widthStep] = membership_input[i][j + widthStep] + relation;
			}
		}
		print_file(membership_tmp, read_image->height, read_image->width, "membership_tmp.txt");

		//再做垂直
		for(j=0; j < read_image->width; j++){
			for(i=0; i < read_image->height - heightStep; i++){
				relation = (float)((uchar *)(read_image->imageData + (i + heightStep)*read_image->widthStep))[j] - (float)((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
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
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 0;
			else if(membership_output[i][j] > 255)
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 255;
			else
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = (int)membership_output[i][j];
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

void curve_line(IplImage *read_image, int row)
{
	int i,j;
	int arr[228];

	
	for(j=0; j < read_image->width; j++){
		arr[j] = ((uchar *)(read_image->imageData + row*read_image->widthStep))[j];
	}


	FILE *fp = fopen("curve_line.txt","w");

	
	for(j=0 ; j < read_image->width; j++){
		fprintf(fp, "%d ", arr[j]);
	}
		
	

	fclose(fp);

	return;
}

void find_minpoint(IplImage *read_image, IplImage *out_image)
{
	int i,j,k;
	int Lpoint, Rpoint;
	bool Lflag, Rflag;
	int Lgray, Rgray;

	IplImage *midpoint_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//找到中心點影像
	IplImage *erosion_image = cvCreateImage(cvSize(read_image->width, read_image->height), IPL_DEPTH_8U,	1);	//侵蝕後影像

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((uchar *)(midpoint_image->imageData + i*out_image->widthStep))[j] = 0;
			((uchar *)(erosion_image->imageData + i*out_image->widthStep))[j] = 0;
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 0;
		}
	}

	Lflag = Rflag = 0;
	for(i=0; i < read_image->height; i++){
		for(j=1; j < read_image->width; j++){
			Lgray = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j-1];
			Rgray = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];

			if((Lgray - Rgray) > 0){
				Lflag = 1;
				Lpoint = j;
			}

			if((Lgray - Rgray) < 0){
				Rflag = 1;
				Rpoint = j-1;
			}

			if(Lflag == 1 && Rflag == 1){
				if((Rpoint - Lpoint) > 0){
					((uchar *)(out_image->imageData + i*out_image->widthStep))[(Lpoint + Rpoint) / 2] = 255;
				}
				Lflag = 0;
				Rflag = 0;
			}
		}
	}

	//erosion(midpoint_image, erosion_image);	//侵蝕
	//dilation(erosion_image, out_image);	//膨脹

	cvReleaseImage(&midpoint_image);		//釋放圖像記憶體
	cvReleaseImage(&erosion_image);		//釋放圖像記憶體
	

	return;



}

void var_line(IplImage *read_image, int row)
{
	/*******************計算兩格間的差距，arr[0]代表0和1之間的差距**************************/
	int i,j;
	int arr[206-1];
	int Left,Right;


	for(j=0; j < read_image->width-1; j++){
		
		arr[j] = ((uchar *)(read_image->imageData + row*read_image->widthStep))[j+1] - ((uchar *)(read_image->imageData + row*read_image->widthStep))[j];
	}


	FILE *fp = fopen("var_line.txt","w");

	
	for(j=0 ; j < read_image->width-1; j++){
		fprintf(fp, "%d ", arr[j]);
	}
		
	

	fclose(fp);

	return;

}
void smooth_line(IplImage *read_image, int row)
{
	/*******************找出在原灰階圖的橫切面，並且「平滑化」**************************/
	int i,j;
	int arr[228];
	int Left, Mid, Right;
	int mean;

	
	for(j=0; j < read_image->width; j++){
		if(j == 0){
			Left = ((uchar *)(read_image->imageData + row*read_image->widthStep))[j];
			Mid = ((uchar *)(read_image->imageData + row*read_image->widthStep))[j];
			Right = ((uchar *)(read_image->imageData + row*read_image->widthStep))[j+1];
		}else if(j == read_image->width - 1){
			Left = ((uchar *)(read_image->imageData + row*read_image->widthStep))[j-1];
			Mid = ((uchar *)(read_image->imageData + row*read_image->widthStep))[j];
			Right = ((uchar *)(read_image->imageData + row*read_image->widthStep))[j];
		}else{
			Left = ((uchar *)(read_image->imageData + row*read_image->widthStep))[j-1];
			Mid = ((uchar *)(read_image->imageData + row*read_image->widthStep))[j];
			Right = ((uchar *)(read_image->imageData + row*read_image->widthStep))[j+1];
		}
		mean = (Left + Mid + Right) / 3;
		arr[j] = mean;
	}


	FILE *fp = fopen("smooth_line.txt","w");

	
	for(j=0 ; j < read_image->width; j++){
		fprintf(fp, "%d ", arr[j]);
	}
		
	

	fclose(fp);

	return;

}



void segmentation(IplImage *gray_image, IplImage *read_image, IplImage *out_image)
{
	int i,j,k;
	bool minFlag, maxFlag;
	int current_gray, next_gray, before_gray;
	int next_index, before_index;
	int min, min_pixel,
		max, max_pixel;
	bool is_right;
	
	int *enhance_localmin = (int *)malloc(sizeof(int) * read_image->width);	//取得二維陣列動態記憶體
	int *gray_localmin = (int *)malloc(sizeof(int) * read_image->width);	//取得二維陣列動態記憶體
	int *gray_line = (int *)malloc(sizeof(int) * read_image->width);	//取得二維陣列動態記憶體
	int *slope_line = (int *)malloc(sizeof(int) * read_image->width);	//取得二維陣列動態記憶體
	int *tmp_line = (int *)malloc(sizeof(int) * read_image->width);	//取得二維陣列動態記憶體
	int *segmen_index = (int *)malloc(sizeof(int) * read_image->width);	//取得二維陣列動態記憶體

	for(i = 0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 0;
		}
	}
	
	
	min = 999;
	max = -999;

	for(i=0; i < read_image->height; i++){
		for(j=0; j < read_image->width; j++){
			enhance_localmin[j] = 0;
			gray_localmin[j] = 0;
			slope_line[j] = 0;
		}

		//先找出enhance的Local min
		maxFlag = 0;
		minFlag = 1;

		for(j=0; j < read_image->width; j++){
			next_index = j+1;
			if(next_index >= read_image->width)
				next_index = read_image->width - 1;

			current_gray = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
			next_gray = ((uchar *)(read_image->imageData + i*read_image->widthStep))[next_index];

			//找到第1個波谷
			if(minFlag == 1){
				if( next_gray > current_gray ){
					enhance_localmin[j] = 1;
					maxFlag = 1;
					minFlag = 0;
				}
			}

			//找到波谷後的第1個波峰
			if(maxFlag == 1){
				if( next_gray < current_gray ){
					maxFlag = 0;
					minFlag = 1;
				}
			}

		}

		//找出gray的Local min
		//gray_line平滑化
		for(j=0; j < read_image->width; j++){
			before_index = j-1;
			next_index = j+1;
			if(before_index < 0)
				before_index = 0;
			if(next_index >= read_image->width)
				next_index = read_image->width - 1;

			before_gray = ((uchar *)(gray_image->imageData + i*gray_image->widthStep))[before_index];
			current_gray = ((uchar *)(gray_image->imageData + i*gray_image->widthStep))[j];
			next_gray = ((uchar *)(gray_image->imageData + i*gray_image->widthStep))[next_index];

			tmp_line[j] = (before_gray + current_gray + next_gray) / 3;
		}
		for(j=0; j < read_image->width; j++){
			before_index = j-1;
			next_index = j+1;
			if(before_index < 0)
				before_index = 0;
			if(next_index >= read_image->width)
				next_index = read_image->width - 1;

			before_gray = tmp_line[before_index];
			current_gray = tmp_line[j];
			next_gray = tmp_line[next_index];

			gray_line[j] = (before_gray + current_gray + next_gray) / 3;
		}
		


		maxFlag = 0;
		minFlag = 1;
		for(j=0; j < read_image->width; j++){
			next_index = j+1;
			if(next_index >= read_image->width)
				next_index = read_image->width - 1;

			current_gray = gray_line[j];
			next_gray = gray_line[next_index];

			//找到第1個波谷
			if(minFlag == 1){
				if( next_gray > current_gray ){
					gray_localmin[j] = 1;
					maxFlag = 1;
					minFlag = 0;
				}
			}

			//找到波谷後的第1個波峰
			if(maxFlag == 1){
				if( next_gray < current_gray ){
					maxFlag = 0;
					minFlag = 1;
				}
			}

		}

		//再找出最低點兩旁的最高點
		for(j=0; j < read_image->width; j++){
			if(enhance_localmin[j] == 1){
				max = -999;
				min = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
				for(k=j; ; k--){
					before_index = k-1;
					if(before_index < 0)
						before_index = 0;
					before_gray = ((uchar *)(read_image->imageData + i*read_image->widthStep))[before_index];

					if(before_gray > min){
						if(before_gray > max){
							max = before_gray;
						}else{
							enhance_localmin[k] = 2;
							break;
						}
					}
					if(k < 0){
						enhance_localmin[0] = 2;
						break;
					}
					
				}

				max = -999;
				for(k=j; ; k++){
					next_index = k+1;
					if(next_index >= read_image->width)
						next_index = read_image->width - 1;
					next_gray = ((uchar *)(read_image->imageData + i*read_image->widthStep))[next_index];
					
					if(next_gray > min){
						if(next_gray > max){
							max = next_gray;
						}else{
							enhance_localmin[k] = 2;
							break;
						}
					}
					if(k >= read_image->width){
						enhance_localmin[k] = 2;
						break;
					}
					
				}
			}
		}

		//篩選正確的enhance 的 Local min
		
		for(j=0; j < read_image->width; j++){
			if(enhance_localmin[j] == 1){
				is_right = 0;

				for(k=j; enhance_localmin[k] != 2 ; k--){
					if(gray_localmin[k] == 1){
						is_right = 1;
						break;
					}
				}
				for(k=j; enhance_localmin[k] != 2 ; k++){
					if(gray_localmin[k] == 1){
						is_right = 1;
						break;
					}
				}

				if(is_right == 0)
					enhance_localmin[j] = 0;
			}
		}

		//找slope_line
		for(j=0; j < read_image->width; j++){
			next_index = j+1;
			if(next_index >= read_image->width)
				next_index = read_image->width - 1;
			current_gray = ((uchar *)(read_image->imageData + i*read_image->widthStep))[j];
			next_gray = ((uchar *)(read_image->imageData + i*read_image->widthStep))[next_index];

			slope_line[j] = next_gray - current_gray;
		}

		//透過var_line的max、min，找到指靜脈左右邊界，並且分割
		for(j=0; j < read_image->width; j++){
			segmen_index[j] = 0;
		}
		//先把要分割出來的部分做標記
		/*for(j=0; j < read_image->width; j++){
			segmen_index[k] = 0;
		}*/
		for(j=0; j < read_image->width; j++){
			if(enhance_localmin[j] == 1){
				max = -999;
				min = 999;
				max_pixel = 0;
				min_pixel = 0;
				for(k=j; ; k--){
					if(slope_line[k] <= min){
						min = slope_line[k];
						min_pixel = k;
					}
					if(enhance_localmin[k] == 2)
						break;
				}
				for(k=j; ; k++){
					if(slope_line[k] >= max){
						max = slope_line[k];
						max_pixel = k;
					}
					if(enhance_localmin[k] == 2)
						break;
				}

				for(k=min_pixel; k <= max_pixel; k++){
					segmen_index[k] = 1;
				}
			}
		}
		
		//對標記過的地方做靜脈分割
		for(j=0; j < read_image->width; j++){
			if(segmen_index[j] == 1){
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 255;
			}else{
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 0;
			}
			//printf("%d ",((uchar *)(out_image->imageData + i*out_image->widthStep))[j]);
		}
		//printf("\n\n");
		
		/*for(j=0; j < read_image->width; j++){
			printf("j:%d %d\n", j, segmen_index[j]);
		}
		
		printf("###########\n\n");*/

	}
	return;
}

void bound_mix_thin(IplImage *thin_image, IplImage *bound_image, IplImage *out_image)	////把細線化和邊界圖 疊合
{
	int i,j;

	for(i=0; i < thin_image->height; i++){
		for(j=1; j < thin_image->width; j++){
			((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 0;
		}
	}

	for(i=0; i < thin_image->height; i++){
		for(j=1; j < thin_image->width; j++){
			if(((uchar *)(thin_image->imageData + i*thin_image->widthStep))[j] == 255 || ((uchar *)(bound_image->imageData + i*bound_image->widthStep))[j] == 255)
				((uchar *)(out_image->imageData + i*out_image->widthStep))[j] = 255;
		}
	}

	return;
}