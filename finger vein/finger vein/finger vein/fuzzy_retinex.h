void fuzzy_retinex(IplImage *read_image, IplImage *out_image);
void membership_function2(IplImage * read_image, int a, int b, int c, float r, float **membership_value);
void find_power_value_retinex(IplImage *read_image, float *fhist, int a, int b, int c, float r, float f, float **membership_value, float **window_fuzzy_entropy, float **power_value);
void defuzzification_retinex(IplImage *read_image, int a, int b, int c, float r, int Lmax, int Lmin, float **modified_mem_value, IplImage *out_image);
void retinex_contrast(IplImage *read_image, float **membership_value, float **membership_enhanced);	//retinex的對比強化