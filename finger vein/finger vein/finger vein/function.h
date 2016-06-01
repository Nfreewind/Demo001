
typedef struct { 
	float real;
	float image; 
}complex;

typedef struct Histogram{
	int arr[256];
};

typedef struct xyw{
	int x,y,w;
};

int is_overflow(int a);	//檢驗是否超出0~255範圍
void scal_trans(int **read_image, IplImage *out_image);	//按比例轉換灰階值
void imshow(IplImage *read_image, char *windowName);	//顯示影像
void U8toF64(IplImage *read_image, IplImage *out_image);	//8U影像轉64F
void mat2grat(int **read_matrix, IplImage *out_image);	//矩陣資料 轉 灰階影像
void print_int_mat_txt(int **in, int height, int width, char *fileName);	//印出整數matrix
//void print_inttxt(int *in, char *fileName);	//印出整數array

void complex_div(int dim, complex **dividend, complex **divisor, complex **quot);	//複數除法


void gray(IplImage *read_image, IplImage *out_image);	//轉成單通道灰階
void binarization2(IplImage *read_image, IplImage *out_image, int TH1, int TH2);	//雙重閥值二值化
void binarization(IplImage *read_image, IplImage *out_image, int TH);	//二值化
void bit_plane(IplImage *read_image, IplImage *out_image, int bitNumber);	//位元平面
void im_resize(IplImage *read_image, IplImage *out_image, int resolution);	//空間解析度
void quantizatioin(IplImage *read_image, IplImage *out_image, int number);	//均勻量化
void dither(IplImage *read_image, IplImage *out_image);	//混色
void multi_dither(IplImage *read_image, IplImage *out_image, int layer);	//量化多層混色

void error_diff(IplImage *read_image, IplImage *out_image);	//誤差擴散法
void im_add(IplImage *read_image, IplImage *out_image, int number);	//影像加法
void im_sub(IplImage *read_image, IplImage *out_image, int number);	//影像減法
void im_multiply(IplImage *read_image, IplImage *out_image, int number);	//影像乘法
void im_divide(IplImage *read_image, IplImage *out_image, int number);	//影像除法
void im_complement(IplImage *read_image, IplImage *out_image);	//影像補色
void calc_Hist(IplImage *read_image);	//直方圖統計
void hist_expand(IplImage *read_image);	//直方圖擴展法
void hist_stretch(IplImage *read_image, IplImage *out_image);	//直方圖擴展(HS)
void sort(IplImage *read_image, xyw* buf, int low);	//按照周圍像素的灰階由高向低排列變化
void hist_equ(IplImage *read_image, IplImage *out_image);	//直方圖均衡(HE)
void LUT_hist_eq(IplImage *read_image, IplImage *out_image);	////LUT直方圖等化
void avg_filter(IplImage *read_image, IplImage *out_image, int dim);	//平均濾波器
void lapla_filter(IplImage *read_image, IplImage *out_image, int dim);	//Laplacian濾波器
void LoG_filter(IplImage *read_image, IplImage *out_image, int dim);	//LoG濾波器
void HPF_filter(IplImage *read_image, IplImage *out_image, int dim);	//高通濾波器
void unsharp(IplImage *read_image, IplImage *out_image, int gain, float k);	//去銳利化
void hb_filter(IplImage *read_image, IplImage *out_image, float A);	//高增幅濾波器
void max_min_filter(IplImage *read_image, IplImage *out_image);	//最大最小濾波器
void mid_filter(IplImage *read_image, IplImage *out_image, int dim);	//中位濾波器
void ROI_foiter(IplImage *read_image, IplImage *out_image, int startI, int endI, int startJ, int endJ, int filterMode);	//ROI濾波

void print_complex(int data_no, complex *in_data);	//印出複數陣列
void DFT1(int data_no, complex *in_data, complex *out_data);	//一維傅立葉轉換
void iDFT1(int data_no, complex *in_data, complex *out_data);	//一維「反」傅立葉轉換
void shift_array(int data_no, complex *in_data, complex *out_data);	//一維數列平移
void conv(int data1_M, complex *data1, int data2_N, complex *data2, complex *out_data);	//旋積
void cconv(int data_no, complex *data1, complex *data2, complex *out_data);	//環形旋積

void DFT2(IplImage *read_image, complex **out_data);	//二維傅立葉轉換
void load_filter(FILE *datap, float **out_data);	//讀取濾波器
void fftshow(int data_no, complex **in_data, char *windowName);	//顯示頻譜圖
void IDFT2(complex **read_data, IplImage *out_image);	//二維傅立葉「反」轉換

void outliers_filter(IplImage *read_image, IplImage *out_image, int dim);	//歧異點濾波器
void adaptive_filter(IplImage *read_image, IplImage *out_image, int dim);	//可適性濾波
void band_reject_folter(IplImage *read_image, IplImage *out_image);	//帶拒濾波
void notch_filter(IplImage *read_image, IplImage *out_image);	//陷波濾波
void remove_butter(IplImage *read_image, IplImage *out_image);	//去除butter雜訊
void remove_motion(IplImage *read_image, IplImage *out_image);	//去除動態模糊
void wiener_filter(IplImage *read_image, IplImage *out_image);	//wiener濾波去除butter雜訊
int Otsu(IplImage *read_image);	//Otsu演算法，找出最佳閥值方法

void fprintf_hist(IplImage *read_image, char *fileName);	//統計影像灰階數，輸出至TXT檔
void difference_filter(IplImage *read_image, IplImage *out_image);	//梯度(通常的差分)
void roberts_filter(IplImage *read_image, IplImage *out_image);	//梯度(Rober)
void sobel_filter(IplImage *read_image, IplImage *out_image);	//梯度(Sobel)
void prewitt(IplImage *read_image, IplImage *out_image);	//用prewitt擷取輪廓
int conj_count(int p[9]);	//計算連結數
void thinning(IplImage *read_image, IplImage *out_image);	//二值影像細線化
void laplacian_filater(IplImage *read_image, IplImage *out_image);	//laplacian二階微分
void zero_cross(IplImage *read_image, IplImage *out_image);	//利用零交錯點，擷取輪廓

void dilation(IplImage *read_image, IplImage *out_image);	//膨脹
void erosion(IplImage *read_image, IplImage *out_image);	//侵蝕
void imopen(IplImage *read_image, IplImage *out_image);	//開運算
void regfill(IplImage *read_image, IplImage *out_image);	//區域填充
void skeleton(IplImage *read_image, IplImage *out_image);	//骨架化
void label(IplImage *read_image, IplImage *out_image, int part);

void MakeLut(Histogram *pLUT,int Min, int Max, int NumBins);
void MakeHistogram(IplImage *read_image, int Xstart, int Ystart, int Xsize, int Ysize, Histogram *pHist, int NumBins, Histogram *pLUT);
void ClipHistogram(Histogram *pHist, int NumBins, int ClipLimit);
void MapHistogram(Histogram *pHist, int Min, int Max, int NumBins, int NumPixels);
void Interpolate(IplImage *read_image, IplImage *out_image, int Xstart, int Ystart, Histogram *pLU, Histogram *pLD, Histogram *pRU, Histogram *pRD, int Xsize, int Ysize, int NumPixels, Histogram *pLUT);
void CLAHE(IplImage *read_image, int NumX, int NumY, int NumBins, float fClipLimit , IplImage *out_image);
void ALRHS(IplImage *read_image, int L, int H, int Z, IplImage *out_image);

void project_statistics(IplImage *read_image);	//取ROI後的血管中心定位















//complex *Alloc1DArray(int size, int item_size);	//取得一維陣列動態記憶體
//complex **Alloc2DArray(int xsize,int ysize,int item_size);	//取得二維陣列動態記憶體
//void Free2DArray(int xsize,int ysize, complex **ptr);	//釋放二維陣列動態記憶體

/*
void move_DFT1(void);	//一維傅立葉平移
void DFT1(void);	//一維離散傅立葉轉換
void inv_DFT1(void);	//一維離散反傅立葉轉換
void conv(void);	//旋積
void cconv(void);	//環形旋積*/

