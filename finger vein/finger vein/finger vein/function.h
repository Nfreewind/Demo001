
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

int is_overflow(int a);	//����O�_�W�X0~255�d��
void scal_trans(int **read_image, IplImage *out_image);	//������ഫ�Ƕ���
void imshow(IplImage *read_image, char *windowName);	//��ܼv��
void U8toF64(IplImage *read_image, IplImage *out_image);	//8U�v����64F
void mat2grat(int **read_matrix, IplImage *out_image);	//�x�}��� �� �Ƕ��v��
void print_int_mat_txt(int **in, int height, int width, char *fileName);	//�L�X���matrix
//void print_inttxt(int *in, char *fileName);	//�L�X���array

void complex_div(int dim, complex **dividend, complex **divisor, complex **quot);	//�Ƽư��k


void gray(IplImage *read_image, IplImage *out_image);	//�ন��q�D�Ƕ�
void binarization2(IplImage *read_image, IplImage *out_image, int TH1, int TH2);	//�����֭ȤG�Ȥ�
void binarization(IplImage *read_image, IplImage *out_image, int TH);	//�G�Ȥ�
void bit_plane(IplImage *read_image, IplImage *out_image, int bitNumber);	//�줸����
void im_resize(IplImage *read_image, IplImage *out_image, int resolution);	//�Ŷ��ѪR��
void quantizatioin(IplImage *read_image, IplImage *out_image, int number);	//���öq��
void dither(IplImage *read_image, IplImage *out_image);	//�V��
void multi_dither(IplImage *read_image, IplImage *out_image, int layer);	//�q�Ʀh�h�V��

void error_diff(IplImage *read_image, IplImage *out_image);	//�~�t�X���k
void im_add(IplImage *read_image, IplImage *out_image, int number);	//�v���[�k
void im_sub(IplImage *read_image, IplImage *out_image, int number);	//�v����k
void im_multiply(IplImage *read_image, IplImage *out_image, int number);	//�v�����k
void im_divide(IplImage *read_image, IplImage *out_image, int number);	//�v�����k
void im_complement(IplImage *read_image, IplImage *out_image);	//�v���ɦ�
void calc_Hist(IplImage *read_image);	//����ϲέp
void hist_expand(IplImage *read_image);	//������X�i�k
void hist_stretch(IplImage *read_image, IplImage *out_image);	//������X�i(HS)
void sort(IplImage *read_image, xyw* buf, int low);	//���өP�򹳯����Ƕ��Ѱ��V�C�ƦC�ܤ�
void hist_equ(IplImage *read_image, IplImage *out_image);	//����ϧ���(HE)
void LUT_hist_eq(IplImage *read_image, IplImage *out_image);	////LUT����ϵ���
void avg_filter(IplImage *read_image, IplImage *out_image, int dim);	//�����o�i��
void lapla_filter(IplImage *read_image, IplImage *out_image, int dim);	//Laplacian�o�i��
void LoG_filter(IplImage *read_image, IplImage *out_image, int dim);	//LoG�o�i��
void HPF_filter(IplImage *read_image, IplImage *out_image, int dim);	//���q�o�i��
void unsharp(IplImage *read_image, IplImage *out_image, int gain, float k);	//�h�U�Q��
void hb_filter(IplImage *read_image, IplImage *out_image, float A);	//���W�T�o�i��
void max_min_filter(IplImage *read_image, IplImage *out_image);	//�̤j�̤p�o�i��
void mid_filter(IplImage *read_image, IplImage *out_image, int dim);	//�����o�i��
void ROI_foiter(IplImage *read_image, IplImage *out_image, int startI, int endI, int startJ, int endJ, int filterMode);	//ROI�o�i

void print_complex(int data_no, complex *in_data);	//�L�X�Ƽư}�C
void DFT1(int data_no, complex *in_data, complex *out_data);	//�@���ť߸��ഫ
void iDFT1(int data_no, complex *in_data, complex *out_data);	//�@���u�ϡv�ť߸��ഫ
void shift_array(int data_no, complex *in_data, complex *out_data);	//�@���ƦC����
void conv(int data1_M, complex *data1, int data2_N, complex *data2, complex *out_data);	//�ۿn
void cconv(int data_no, complex *data1, complex *data2, complex *out_data);	//���αۿn

void DFT2(IplImage *read_image, complex **out_data);	//�G���ť߸��ഫ
void load_filter(FILE *datap, float **out_data);	//Ū���o�i��
void fftshow(int data_no, complex **in_data, char *windowName);	//����W�й�
void IDFT2(complex **read_data, IplImage *out_image);	//�G���ť߸��u�ϡv�ഫ

void outliers_filter(IplImage *read_image, IplImage *out_image, int dim);	//�[���I�o�i��
void adaptive_filter(IplImage *read_image, IplImage *out_image, int dim);	//�i�A���o�i
void band_reject_folter(IplImage *read_image, IplImage *out_image);	//�a���o�i
void notch_filter(IplImage *read_image, IplImage *out_image);	//���i�o�i
void remove_butter(IplImage *read_image, IplImage *out_image);	//�h��butter���T
void remove_motion(IplImage *read_image, IplImage *out_image);	//�h���ʺA�ҽk
void wiener_filter(IplImage *read_image, IplImage *out_image);	//wiener�o�i�h��butter���T
int Otsu(IplImage *read_image);	//Otsu�t��k�A��X�̨λ֭Ȥ�k

void fprintf_hist(IplImage *read_image, char *fileName);	//�έp�v���Ƕ��ơA��X��TXT��
void difference_filter(IplImage *read_image, IplImage *out_image);	//���(�q�`���t��)
void roberts_filter(IplImage *read_image, IplImage *out_image);	//���(Rober)
void sobel_filter(IplImage *read_image, IplImage *out_image);	//���(Sobel)
void prewitt(IplImage *read_image, IplImage *out_image);	//��prewitt�^������
int conj_count(int p[9]);	//�p��s����
void thinning(IplImage *read_image, IplImage *out_image);	//�G�ȼv���ӽu��
void laplacian_filater(IplImage *read_image, IplImage *out_image);	//laplacian�G���L��
void zero_cross(IplImage *read_image, IplImage *out_image);	//�Q�ιs����I�A�^������

void dilation(IplImage *read_image, IplImage *out_image);	//����
void erosion(IplImage *read_image, IplImage *out_image);	//�I�k
void imopen(IplImage *read_image, IplImage *out_image);	//�}�B��
void regfill(IplImage *read_image, IplImage *out_image);	//�ϰ��R
void skeleton(IplImage *read_image, IplImage *out_image);	//���[��
void label(IplImage *read_image, IplImage *out_image, int part);

void MakeLut(Histogram *pLUT,int Min, int Max, int NumBins);
void MakeHistogram(IplImage *read_image, int Xstart, int Ystart, int Xsize, int Ysize, Histogram *pHist, int NumBins, Histogram *pLUT);
void ClipHistogram(Histogram *pHist, int NumBins, int ClipLimit);
void MapHistogram(Histogram *pHist, int Min, int Max, int NumBins, int NumPixels);
void Interpolate(IplImage *read_image, IplImage *out_image, int Xstart, int Ystart, Histogram *pLU, Histogram *pLD, Histogram *pRU, Histogram *pRD, int Xsize, int Ysize, int NumPixels, Histogram *pLUT);
void CLAHE(IplImage *read_image, int NumX, int NumY, int NumBins, float fClipLimit , IplImage *out_image);
void ALRHS(IplImage *read_image, int L, int H, int Z, IplImage *out_image);

void project_statistics(IplImage *read_image);	//��ROI�᪺��ޤ��ߩw��















//complex *Alloc1DArray(int size, int item_size);	//���o�@���}�C�ʺA�O����
//complex **Alloc2DArray(int xsize,int ysize,int item_size);	//���o�G���}�C�ʺA�O����
//void Free2DArray(int xsize,int ysize, complex **ptr);	//����G���}�C�ʺA�O����

/*
void move_DFT1(void);	//�@���ť߸�����
void DFT1(void);	//�@�������ť߸��ഫ
void inv_DFT1(void);	//�@�������ϳť߸��ഫ
void conv(void);	//�ۿn
void cconv(void);	//���αۿn*/

