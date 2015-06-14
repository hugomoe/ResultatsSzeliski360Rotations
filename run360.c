#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "ftr.h"
#include "iio.h"
#include "ftr.c"
#include "iio.c"

/** description de run360.c un peu plus bas */

///choose an interpolation
//#define _IMG_TO_DOUBLE //cast all images as double
#define INTERPOLATION_TYPE 0 //reste à faire toutes les b splines et les tests en type double
    /*
     * interpolation type :
     * 0 -> raised cosine-weighted sinc
     * 1 -> cubic spline
     * 2 -> b-spline
     * 3 -> naive use of b-spline
     */
#define BETA 0.36 //beta for raised cosine (0.36+)
#define ORDER 3 //order of the B-Spline (3,5 float and 7,9,11 double)
#define NB_ROTATION 360 //must be even
///end of the choice

/**
  * Permet de tester 360 (ou NB_ROTATION) rotations et les comparer à l'identité.
  * 
  * Imprime une image .png contenant l'image tournée
  * Imprime un fichier .txt associé contenant les résultats numériques de cette comparaison
  * Écrit à la suite de all_result.txt ce qui a été imprimé dans le .txt précédemment mentionné
  * 
  * Avant de commencer ses propres tests sur une nouvelle machine, renommer le fichier all_results.txt afin de ne pas mélanger les résultats
  * Le choix de l'interpolation se fait via les #define ci-dessus
  * 
  * Les images imprimées sont au format .png
  */



//list of possible index for interpolation type
    #define RAISED_COSINE_WEIGHTED_SINC 0
    #define SPLINE_CUBIC 1
    #define B_SPLINE 2
    #define NAIVE_B_SPLINE 3



#ifdef _IMG_TO_DOUBLE
///process images with type double
#include "affineDouble360.h"
#include "affineDouble360.c"

int main(int argc, char *argv[]){
	if (argc != 2) {
		fprintf(stderr, "usage:\n\t%s [image.extension]\n", *argv);
		return 1;
	}
	char *filename_in = argv[1];



    ///get the input
	int w,h,pd;
    double *img_temp = iio_read_image_double_vec(filename_in,&w,&h,&pd);
    double *img_in = malloc(3*w*h*sizeof(double));
    if(pd!=3){
        for(int l=0;l<3;l++){
            for(int i=0;i<w;i++){
                for(int j=0;j<h;j++){
                    img_in[(i+w*j)*3+l] = img_temp[(i+w*j)*pd];
                }
            }
        }
    }else{
        for(int k=0;k<3*w*h;k++){
            img_in[k]=img_temp[k];
        }
    }
    pd = 3;
    free(img_temp);



    ///create intermediate images and matrices
    double theta = 2*3.14159265358979323/NB_ROTATION;
    double A[6] = {cos(theta),sin(theta),0.,-sin(theta),cos(theta),0.};

    int w2 = 3*fmax(w,h)/2; int h2 = w2;
    double *img_a = malloc(3*w2*h2*sizeof(double));
    double *img_b = malloc(3*w2*h2*sizeof(double));
    double *img_out = malloc(3*w*h*sizeof(double));

    double rms=0., moy=0., sup=0., duration = 0., temp;
    clock_t start, end;



    ///run the tests
    int k = 0;
    printf("---------- starting the experiments ----------\n");
    start = clock();

    apply_affinite(img_in,img_a,w,h,w2,h2,A);
    k++; printf("%d-",k); fflush(stdout);
    //for(int i=0;i<179;i++){
    for(int i=0;i<NB_ROTATION/2-1;i++){
        apply_affinite(img_a,img_b,w2,h2,w2,h2,A);
        k++; printf("%d-",k); fflush(stdout);
        apply_affinite(img_b,img_a,w2,h2,w2,h2,A);
        k++; printf("%d-",k); fflush(stdout);
    }
    apply_affinite(img_a,img_out,w2,h2,w,h,A);
    k++; printf("%d : done !\n",k); fflush(stdout);

    end = clock();
    printf("----------  end of the experiments  ----------\n");
    duration = (double)(end-start)/CLOCKS_PER_SEC;
    for(int l=0;l<3;l++){
        for(int i=0;i<w;i++){
            for(int j=0;j<h;j++){
                temp = fabs(img_in[(i+w*j)*3+l]-img_out[(i+w*j)*3+l]);
                rms += pow(temp,2);
                moy += temp;
                sup = (sup<temp) ? temp : sup;
            }
        }
    }
    rms = sqrt(rms/w/h/pd);
    moy = moy/w/h/pd;
    printf("erreur rms : %f\nerreur L1 : %f\nerreur max : %f\nduree : %fs\n",rms,moy,sup,duration);



    ///save output
    //output image
    char filename_out[FILENAME_MAX]; //name of the output image
    if(filename_out==NULL){
        printf("filename_out prend trop de place\n");
        exit(1);
    }
    char str_nb_rot[FILENAME_MAX];
    if(str_nb_rot==NULL){
        printf("str_nb_rot prend trop de place\n");
        exit(1);
    }
    switch(NB_ROTATION){
    case 360:
        {
        sprintf(str_nb_rot,"");
        }
        break;
    default:
        {
        sprintf(str_nb_rot,"_%drot",NB_ROTATION);
        }
        break;
    }
    switch(INTERPOLATION_TYPE){
    case RAISED_COSINE_WEIGHTED_SINC:
        sprintf(filename_out,"raised-cosine_beta%3.2f_double%s_%s",BETA,str_nb_rot,filename_in);
        break;
    case SPLINE_CUBIC:
        sprintf(filename_out,"cubic-spline_double%s_%s",str_nb_rot,filename_in);
        break;
    case B_SPLINE:
        sprintf(filename_out,"b-spline_order%d_double%s_%s",ORDER,str_nb_rot,filename_in);
        break;
    case NAIVE_B_SPLINE:
        sprintf(filename_out,"naive-b-spline_order%d_double%s_%s",ORDER,str_nb_rot,filename_in);
        break;
    default:
        sprintf(filename_out,"modified_double%s_%s",str_nb_rot,filename_in);
        break;
    }
    float *img_out_float = malloc(3*w*h*sizeof(float));
    if(img_out_float==NULL){
        printf("img_out_float prend trop de place\n");
        exit(1);
    }
    for(int i=0;i<3*w*h;i++){img_out_float[i] = (float)(img_out[i]);}
	iio_save_image_float_vec(filename_out, img_out_float, w, h, pd);

    //output txt and erase the latest test
    FILE *latest_result_file;
    char latest_result_filename[FILENAME_MAX]; //name of the output txt
    sprintf(latest_result_filename,"%s_results.txt",filename_out);
    latest_result_file = fopen(latest_result_filename,"w"); //w : erase the current file to create another
    if (latest_result_file==NULL) {
        printf("Can't create output file %s\n",latest_result_filename);
        exit(1);
    }
    fprintf(latest_result_file,"%s latest test :\n\terreur rms : %f\n\terreur L1 : %f\n\terreur max : %f\n\tduree : %fs\n",filename_out,rms,moy,sup,duration);
    fclose(latest_result_file);

    //output txt with all the others results already acquired
    FILE *results_file;
    char *results_filename = "all_results.txt";
    results_file = fopen(results_filename,"a"); //a : append the text (write it after the end of the file)
    if (results_file==NULL) {
        printf("Can't create output file %s\n",results_filename);
        exit(1);
    }
    fprintf(results_file,"%s :\n\terreur rms : %f\n\terreur L1 : %f\n\terreur max : %f\n\tduree : %fs\n",filename_out,rms,moy,sup,duration);
    fclose(results_file);



	///end the program
    free(img_in);
    free(img_a);
    free(img_b);
    free(img_out);
	free(img_out_float);
	return 0;
}



#else // _IMG_TO_DOUBLE
///process images with type float
#include "affine360.h"
#include "affine360.c"

int main(int argc, char *argv[]){
	if (argc != 2) {
		fprintf(stderr, "usage:\n\t%s [image.extension]\n", *argv);
		return 1;
	}
	char *filename_in = argv[1];



    ///get the input
	int w,h,pd;
    float *img_temp = iio_read_image_float_vec(filename_in,&w,&h,&pd);
    float *img_in = malloc(3*w*h*sizeof(float));
    if(pd!=3){
        for(int l=0;l<3;l++){
            for(int i=0;i<w;i++){
                for(int j=0;j<h;j++){
                    img_in[(i+w*j)*3+l] = img_temp[(i+w*j)*pd];
                }
            }
        }
    }else{
        for(int k=0;k<3*w*h;k++){
            img_in[k]=img_temp[k];
        }
    }
    pd = 3;
    free(img_temp);



    ///create intermediate images and matrices
    double theta = 2*3.14159265358979323/NB_ROTATION;
    double A[6] = {cos(theta),sin(theta),0.,-sin(theta),cos(theta),0.};

    int w2 = 3*fmax(w,h)/2; int h2 = w2;
    float *img_a = malloc(3*w2*h2*sizeof(float));
    float *img_b = malloc(3*w2*h2*sizeof(float));
    float *img_out = malloc(3*w*h*sizeof(float));

    float rms=0., moy=0., sup=0., duration = 0., temp;
    clock_t start, end;



    ///run the tests
    int k = 0;
    printf("---------- starting the experiments ----------\n");
    start = clock();

    apply_affinite(img_in,img_a,w,h,w2,h2,A);
    k++; printf("%d-",k); fflush(stdout);
    //for(int i=0;i<179;i++){
    for(int i=0;i<NB_ROTATION/2-1;i++){
        apply_affinite(img_a,img_b,w2,h2,w2,h2,A);
        k++; printf("%d-",k); fflush(stdout);
        apply_affinite(img_b,img_a,w2,h2,w2,h2,A);
        k++; printf("%d-",k); fflush(stdout);
    }
    apply_affinite(img_a,img_out,w2,h2,w,h,A);
    k++; printf("%d : done !\n",k); fflush(stdout);

    end = clock();
    printf("----------  end of the experiments  ----------\n");
    duration = (double)(end-start)/CLOCKS_PER_SEC;
    for(int l=0;l<3;l++){
        for(int i=0;i<w;i++){
            for(int j=0;j<h;j++){
                temp = fabs(img_in[(i+w*j)*3+l]-img_out[(i+w*j)*3+l]);
                rms += pow(temp,2);
                moy += temp;
                sup = (sup<temp) ? temp : sup;
            }
        }
    }
    rms = sqrt(rms/w/h/pd);
    moy = moy/w/h/pd;
    printf("erreur rms : %f\nerreur L1 : %f\nerreur max : %f\nduree : %fs\n",rms,moy,sup,duration);



    ///save output
    //output image
    char filename_out[FILENAME_MAX]; //name of the output image
    switch(NB_ROTATION){
    case 360:
        {
        char *str_nb_rot = "";
        }
        break;
    default:
        {
        char str_nb_rot[FILENAME_MAX];
        if(str_nb_rot==NULL){
            printf("str_nb_rot prend trop de place\n");
            exit(1);
        }
        sprintf(str_nb_rot,"_%drot",NB_ROTATION);
        }
        break;
    }
    switch(INTERPOLATION_TYPE){
    case RAISED_COSINE_WEIGHTED_SINC:
        sprintf(filename_out,"raised-cosine_beta%3.2f%s_%s",BETA,str_nb_rot,filename_in);
        break;
    case SPLINE_CUBIC:
        sprintf(filename_out,"cubic-spline%s_%s",str_nb_rot,filename_in);
        break;
    case B_SPLINE:
        sprintf(filename_out,"b-spline_order%d%s_%s",ORDER,str_nb_rot,filename_in);
        break;
    case NAIVE_B_SPLINE:
        sprintf(filename_out,"naive-b-spline_order%d%s_%s",ORDER,str_nb_rot,filename_in);
        break;
    default:
        sprintf(filename_out,"modified%s_%s",str_nb_rot,filename_in);
        break;
    }
	iio_save_image_float_vec(filename_out, img_out, w, h, pd);

    //output txt and erase the latest test
    FILE *latest_result_file;
    char latest_result_filename[FILENAME_MAX]; //name of the output txt
    sprintf(latest_result_filename,"%s_results.txt",filename_out);
    latest_result_file = fopen(latest_result_filename,"w"); //w : erase the current file to create another
    if (latest_result_file==NULL) {
        printf("Can't create output file %s\n",latest_result_filename);
        exit(1);
    }
    fprintf(latest_result_file,"%s latest test :\n\terreur rms : %f\n\terreur L1 : %f\n\terreur max : %f\n\tduree : %fs\n",filename_out,rms,moy,sup,duration);
    fclose(latest_result_file);

    //output txt with all the others results already acquired
    FILE *results_file;
    char *results_filename = "all_results.txt";
    results_file = fopen(results_filename,"a"); //a : append the text (write it after the end of the file)
    if (results_file==NULL) {
        printf("Can't create output file %s\n",results_filename);
        exit(1);
    }
    fprintf(results_file,"%s :\n\terreur rms : %f\n\terreur L1 : %f\n\terreur max : %f\n\tduree : %fs\n",filename_out,rms,moy,sup,duration);
    fclose(results_file);



	///end the program
    free(img_in);
    free(img_a);
    free(img_b);
    free(img_out);
	return 0;
}



#endif // _IMG_TO_DOUBLE
