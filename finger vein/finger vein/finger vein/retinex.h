void retinex(IplImage *read_image, IplImage *out_image);
void curve_line(IplImage *read_image, int row);	//��������u
void find_minpoint(IplImage *read_image, IplImage *out_image);	//���ޤ����I
void var_line(IplImage *read_imageint, int row);	//�ܲ��q���u
void smooth_line(IplImage *read_image, int row);	//��Ͼ�������u(���Ƥ�)
void segmentation(IplImage *gray_image, IplImage *read_image, IplImage *out_image);	//���R�ߤ���
void bound_mix_thin(IplImage *thin_image, IplImage *bound_image, IplImage *out_image);	//��ӽu�ƩM��ɹ� �|�X