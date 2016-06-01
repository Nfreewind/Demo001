

void print_file(float **file_point, int height, int width, char *filename);
void membership_function1(IplImage *read_image, int a, int b, int c, float **membership_value);
float entropy(IplImage *read_image, float **membership_value);
void edge_function(IplImage *read_image, float **membership_value, float **edge_value);
void window_entropy(IplImage *read_image, int windowSize, float **membership_value, float **edge_value, float **window_fuzzy_entropy);
void find_power_value(IplImage *read_image, float *fhist, int a, int b, int c, float f, float **membership_value, float **window_fuzzy_entropy, float **power_value);
void find_avg_edge_value(IplImage *read_image, int windowSize, float **membership_value, float **edge_value, float **avg_edge_value);
void find_contrast_membership_value(IplImage *read_image, float **membership_value, float **avg_edge_value, float **contrast_mem_value);
void find_transform_contrast_membership_value(IplImage *read_image, float **contrast_mem_value, float **power_value, float **trans_contrast_mem_value);
void find_modified_membership_value(IplImage *read_image, float **membership_value, float **avg_edge_value, float **trans_contrast_mem_value, float **modified_mem_value);
void defuzzification(IplImage *read_image, int a, int b, int c, int Lmax, int Lmin, float **modified_mem_value, IplImage *out_image);

void fuzzy_contrast(IplImage *read_image, IplImage *out_image);