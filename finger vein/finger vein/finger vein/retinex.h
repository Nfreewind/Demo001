void retinex(IplImage *read_image, IplImage *out_image);
void curve_line(IplImage *read_image, int row);	//橫切面曲線
void find_minpoint(IplImage *read_image, IplImage *out_image);	//找血管中心點
void var_line(IplImage *read_imageint, int row);	//變異量曲線
void smooth_line(IplImage *read_image, int row);	//原圖橫切面曲線(平滑化)
void segmentation(IplImage *gray_image, IplImage *read_image, IplImage *out_image);	//指靜脈分割
void bound_mix_thin(IplImage *thin_image, IplImage *bound_image, IplImage *out_image);	//把細線化和邊界圖 疊合