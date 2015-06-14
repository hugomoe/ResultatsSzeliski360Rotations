#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "umax_vmax_double.h"
#include "spline_double.c"

#include "affineDouble360.h"

bool eq(double a,double b){if(a<b+pow(2,-PREC2)&& a>b-pow(2,-PREC2)){return true;}{return false;}}

void transpo_opt(double *img,double *a,int wh[2]){
	double c1,c2;
	c1 = sqrt(pow(a[0],2)+pow(a[1],2));
	c2 = sqrt(pow(a[3],2)+pow(a[4],2));
	double a0,a1,a3,a4;
	a1=fabs(a[1])/c1; a0=fabs(a[0])/c1; a3=fabs(a[3])/c2; a4=fabs(a[4])/c2;
	if(a0+a4<a1+a3){
		int w = wh[0];
		int h = wh[1];
	//inverse affinite
		double aa[6];
		aa[0]=a[3]; aa[1]=a[4]; aa[2]=a[5]; aa[3]=a[0]; aa[4]=a[1]; aa[5]=a[2];
		a[0]=aa[0]; a[1]=aa[1]; a[2]=aa[2]; a[3]=aa[3]; a[4]=aa[4]; a[5]=aa[5];

	//inverse img
		double *img_t = malloc(3*w*h*sizeof(double));
        if(img_t==NULL){
            printf("transpo_opt : img_t prend trop de place\n");
            exit(1);
        }
		int i,j,l;
		for(l=0;l<3;l++){
			for(i=0;i<w;i++){
				for(j=0;j<h;j++){
					img_t[(j+h*i)*3+l]=img[(i+j*w)*3+l];
				}
			}
		}
		for(i=0;i<w*h*3;i++){img[i]=img_t[i];}

	//inverse w / h
		wh[0]=h; wh[1]=w;

		free(img_t);
	}
}



double filter_fun(double x){
    switch(INTERPOLATION_TYPE){
    case RAISED_COSINE_WEIGHTED_SINC:
    default:
        if(eq(x,0)){return 1.;}
        else if(eq(fabs(x),PERIOD/2./BETA)){return BETA/2.*sin(PI/2./BETA);}
        else{return sin(PI*x/PERIOD)/(PI*x/PERIOD)*cos(PI*BETA*x/PERIOD)/(1.-pow(2.*BETA*x/PERIOD,2));}
        break;
    }
}



double filter_h(double *img,int w,int h,double *H,int N,int prec,double a0,double a1,int i,int j,int l){
	double x = i;
	double y = j;
	double xc = ((double)w-1.)/2.; //l'abscisse du centre, inchangé par l'opération
	double yc = ((double)h-1.)/2.; //l'ordonnée du centre, inchangé par l'opération
	double M = a0*(x-xc) + a1*(y-yc) + xc;
	int k = floor(M);
	int p = floor(prec*(M-k));
	double a = 0;

	int u = k-N;
	if(u<0){u=0;}

	for(;u<=k+N && u<w;u++){
        a += H[(k-u)*prec + p + N*prec]*img[(u+j*w)*3+l];
	}

    return a;
}

double filter_v(double *img,int w,int h,double *H,int N,int prec,double a0,double a1,int i,int j,int l){
	double x = i;
	double y = j;
	double xc = ((double)w-1.)/2.; //l'abscisse du centre, inchangé par l'opération
	double yc = ((double)h-1.)/2.; //l'ordonnée du centre, inchangé par l'opération
	double M = a1*(x-xc) + a0*(y-yc) + yc;
	int k = floor(M);
	int p = floor(prec*(M-k));
	double a = 0;

	int u = k-N;
	if(u<0){u=0;}

	for(;u<=k+N && u<h;u++){
        a += H[(k-u)*prec + p + N*prec]*img[(i+u*w)*3+l];
	}

    return a;
}



int apply_rh(double *img1,double *img2,int w,int h,double s,double a0,double a1){

	if(eq(a0,1.) && eq(a1,0.)){
		for(int i=0;i<3*w*h;i++){img2[i]=img1[i];}
		return 0;
	}

    switch(INTERPOLATION_TYPE){
    case SPLINE_CUBIC:
        {
        double xc = ((double)w-1.)/2., yc = ((double)h-1.)/2.;
        for(int i=0;i<w;i++){
            for(int j=0;j<h;j++){
                double x = i, y = j;

                double M = a0*(x-xc) + a1*(y-yc) + xc; //l'index en lequel on veut interpoler img1
                int k = floor(M); double phi = M - (double)(k);

                int index[4] = {k-1,k,k+1,k+2}; //index sont les points autour de M
                double v[4];
                for(int l=0;l<3;l++){
                    for (int u = 0;u<4;u++){
                        if(index[u]<0 || index[u]>w-1){v[u] = 0;}
                        else{v[u] = img1[(index[u]+w*j)*3+l];}
                    }

                    //cubic spline
                    img2[(i+j*w)*3 + l] = v[1] + 0.5 * phi*(v[2] - v[0]
                                        + phi*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
                                        + phi*(3.0*(v[1] - v[2]) + v[3] - v[0])));
                }
            }
        }
        }
        break; //fin SPLINE_CUBIC
    case B_SPLINE:
        {
        double *img_temp = malloc(3*w*h*sizeof(double));
        if(img_temp==NULL){
            printf("apply_rh : img_temp prend trop de place\n");
            exit(1);
        }
        double *img2ij = malloc(3*sizeof(double));
        if(img2ij==NULL){
            printf("apply_rh : img2ij prend trop de place\n");
            exit(1);
        }


        double x,y,M,xc,yc;
        if(img_temp==NULL){printf("apply_rh, b-spline : Impossible de recopier l'image\n");exit(1);}
        for(int i=0;i<3*w*h;i++){img_temp[i]=img1[i];}
        bool r = prepare_spline(img_temp,w,h,3,ORDER);
        if(!r){printf("apply_rh, b-spline : Impossible de preparer la spline\n");exit(1);}

        xc = ((double)w-1.)/2.; yc = ((double)h-1.)/2.;
        for(int i=0;i<w;i++){
            for(int j=0;j<h;j++){
                x = i; y = j;

                M = a0*(x-xc) + a1*(y-yc) + xc; //l'index en lequel on veut interpoler img1

                evaluate_spline_at(img2ij, img_temp, w, h, 3, ORDER, M, y);
                for(int l=0;l<3;l++){img2[(i+w*j)*3+l] = img2ij[l];}
            }
        }

        free(img_temp);
        free(img2ij);
        }
        break;
    case RAISED_COSINE_WEIGHTED_SINC:
    default:
        {
        int N = ceil(abs(s*TAPS));
        int prec = pow(2,PREC);
        double precf = (double) prec;
        //double sf = (double) s;
        double *H = malloc((2*N+1)*prec*sizeof(double));
        if(H==NULL){printf("apply_rh : H n'a pas ete cree");exit(1);}

        int k,p,i,j;

        for(p=0;p<prec;p++){
            double Htot = 0;
            double pf = (double) p;
            for(k=-N;k<=N;k++){
                    double kf = (double) k;
                Htot += H[k*prec + p + N*prec] = filter_fun((kf+pf/precf)/s);
            }
            for(k=-N;k<=N;k++){H[k*prec + p + N*prec] = H[k*prec + p + N*prec]/Htot;}
        }

        for(int l=0;l<3;l++){
            for(j=0;j<h;j++){
                for(i=0;i<w;i++){
                    img2[(i+j*w)*3 + l] = filter_h(img1,w,h,H,N,prec,a0,a1,i,j,l);
                }
            }
        }

        free(H);
        }
        break; //fin du cas filtre de convolution
    }
	return 0;
}

int apply_rv(double *img1,double *img2,int w,int h,double s,double a0,double a1){

	if(eq(a0,1.) && eq(a1,0.)){
		for(int i=0;i<3*w*h;i++){
            img2[i]=img1[i];
        }
		return 0;
	}

    switch(INTERPOLATION_TYPE){
    case SPLINE_CUBIC:
        {
        double xc = ((double)w-1.)/2., yc = ((double)h-1.)/2.;
            for(int j=0;j<h;j++){
                for(int i=0;i<w;i++){
                    double x = i, y = j;

                    double M = a1*(x-xc) + a0*(y-yc) + yc;
                    int k = floor(M); double phi = M - (double)(k);

                    int index[4] = {k-1,k,k+1,k+2};
                    double v[4];
                    for(int l=0;l<3;l++){
                        for (int u = 0;u<4;u++){
                            if(index[u]<0 || index[u]>h-1){v[u] = 0;}
                            else{v[u] = img1[(i+w*index[u])*3+l];}
                        }

                        //cubic spline
                        img2[(i+j*w)*3 + l] =  v[1] + 0.5 * phi*(v[2] - v[0]
                                            + phi*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
                                            + phi*(3.0*(v[1] - v[2]) + v[3] - v[0])));
                    }
                }
            }
        }
        break; //fin CUBIC_SPLINE
    case B_SPLINE:
        {
        double *img_temp = malloc(3*w*h*sizeof(double));
        if(img_temp==NULL){
            printf("apply_rv : img_temp prend trop de place\n");
            exit(1);
        }
        double *img2ij = malloc(3*sizeof(double));
        if(img2ij==NULL){
            printf("apply_rv : img2ij prend trop de place\n");
            exit(1);
        }
        double x,y,M,xc,yc;
        if(img_temp==NULL){printf("apply_rv, b-spline : Impossible de recopier l'image\n");exit(1);}
        for(int i=0;i<3*w*h;i++){img_temp[i]=img1[i];}
        bool r = prepare_spline(img_temp,w,h,3,ORDER);
        if(!r){printf("apply_rv, b-spline : Impossible de preparer la spline\n");exit(1);}

        xc = ((double)w-1.)/2.; yc = ((double)h-1.)/2.;
        for(int i=0;i<w;i++){
            for(int j=0;j<h;j++){
                x = i; y = j;

                M = a1*(x-xc) + a0*(y-yc) + yc;

                evaluate_spline_at(img2ij, img_temp, w, h, 3, ORDER, x, M);
                for(int l=0;l<3;l++){
                    img2[(i+w*j)*3+l] = img2ij[l];
                }
            }
        }

        free(img_temp);
        free(img2ij);
        }
        break;
    case RAISED_COSINE_WEIGHTED_SINC:
    default:
        {
        int N = ceil(abs(s*TAPS));
        int prec = pow(2,PREC);
        double sf = (double) s;
        double precf = (double) prec;
        double *H = malloc((2*N+1)*prec*sizeof(double));
        if(H==NULL){printf("apply_rv : H n'a pas ete cree");exit(1);}

        int k,p,i,j;

        for(p=0;p<prec;p++){
            double Htot = 0;
            double pf = (double) p;
            for(k=-N;k<=N;k++){
                double kf = (double) k;
                Htot += H[k*prec + p + N*prec] = filter_fun((kf+pf/precf)/s);
            }
            for(k=-N;k<=N;k++){H[k*prec + p + N*prec] = H[k*prec + p + N*prec]/Htot;}
        }

        for(int l=0;l<3;l++){
            for(i=0;i<w;i++){
                for(j=0;j<h;j++){
                    img2[(i+j*w)*3+l] = filter_v(img1,w,h,H,N,prec,a0,a1,i,j,l);
                }
            }
        }

        free(H);
        }
        break; //fin du cas filtre de convolution
    }
	return 0;
}



int apply_affinite(double *img,double *img_f,int w,int h,int w_f,int h_f,double *affinity){
    switch(INTERPOLATION_TYPE){
    case RAISED_COSINE_WEIGHTED_SINC:
    case SPLINE_CUBIC:
    case B_SPLINE:
        {
        int i,j,l;
        //copie des arguments
        //(pour ne pas modifier les arguments et pouvoir les rappeler lors des tests)
        double *img_i = malloc(3*w*h*sizeof(double));
        if(img_i==NULL){
            printf("apply_affinite : img_i prend trop de place\n");
            exit(1);
        }
        for(i=0;i<3*w*h;i++){
            img_i[i]=img[i];
        }
        double a[6];
        for(i=0;i<6;i++){a[i]=affinity[i];}

        int wh[2] = {w,h};
        transpo_opt(img_i,a,wh);
        w=wh[0];
        h=wh[1];

        ///umax,vmax
        double umax,vmax;
        double A[2][2] = {a[0],a[1],a[3],a[4]};
        int test = umax_vmax(&umax,&vmax,A);
        if(test==0){
            printf("@apply_affinite : erreur dans umax_vmax\n");
            exit(1);
        }

        ///coefficients des shear élémentaires
        double b0 = a[0] - a[1]*a[3]/a[4];
        double t2 = a[2] - a[1]*a[5]/a[4];
        //rv
        double aux_rv = fabs(a[4])*vmax;
        if(aux_rv > 1){aux_rv = 1;}
        double rv = fabs(a[1])*umax + aux_rv;
        if(rv>3){rv = 3;}else{if(rv<1){rv = 1;}}
        //rh
        double aux_rh = fabs(b0)*umax;
        if(aux_rh > 1){aux_rh = 1;}
        double rh = fabs(a[3]/a[4])*rv*vmax + aux_rh;
        if(rh>3){rh = 3;}else{if(rh<1){rh = 1;}}

        ///pour le calcul des tailles des images intermédiaires
        //ww et hh : plus grand que w et h pour conserver toute l'image
        int ww,hh;
        if(w<w_f){ww = w_f;} else {ww = w;}
        if(h<h_f){hh = h_f;} else {hh = h;}
        //dw et dh : différence de taille entre img_i et img_f
        int dw = (w_f-w)/2, dh = (h_f-h)/2; //x_f et x ont même parité
        //dwp et dhp leur partie positive (xp=max(0,x))
        int dwp = (dw<0) ? 0 : dw;
        int dhp = (dh<0) ? 0 : dh;
        //dwn et dhn leur partie négative (xn=max(0,-x))
        int dwn = (-dw<0) ? 0 : -dw;
        int dhn = (-dh<0) ? 0 : -dh;



        double *img1 = malloc(3*9*ww*hh*sizeof(double));
        if(img1==NULL){
            printf("apply_affinite : img1 prend trop de place\n");
            exit(1);
        }
        double *img2 = malloc(3*9*ww*hh*sizeof(double));
        if(img2==NULL){
            printf("apply_affinite : img2 prend trop de place\n");
            exit(1);
        }

        ///condition aux bords : symétrisation
        int i_sym,j_sym;
        for(l=0;l<3;l++){
            for(i=0;i<3*ww;i++){
                for(j=0;j<3*hh;j++){
                    i_sym = i-ww-dwp;
                    while(i_sym<0 || i_sym>w-1){i_sym = (i_sym<0) ? -1-i_sym : 2*w-1-i_sym;}
                    j_sym = j-hh-dhp;
                    while(j_sym<0 || j_sym>h-1){j_sym = (j_sym<0) ? -1-j_sym : 2*h-1-j_sym;}
                    img1[(i+3*ww*j)*3+l]=img_i[(i_sym+j_sym*w)*3+l];
                }
            }
        }

        ///application des opérations élémentaires
        apply_rv(img1,img2,3*ww,3*hh,1/vmax,a[4]/rv,0.);
        apply_rh(img2,img1,3*ww,3*hh,1/umax,b0/rh,a[1]/rv);
        apply_rv(img1,img2,3*ww,3*hh,rv,rv,a[3]*rv/a[4]/rh);
        apply_rh(img2,img1,3*ww,3*hh,rh,rh,0.);

        ///output
        for(l=0;l<3;l++){
            for(i=0;i<w_f;i++){
                for(j=0;j<h_f;j++){
                    img_f[(i+j*w_f)*3+l] = img1[(i+ww+dwn+(j+hh+dhn)*ww*3)*3+l];
                }
            }
        }

        free(img_i);
        free(img1);
        free(img2);
        }
	break;
	case NAIVE_B_SPLINE:
        {
        ///pour le calcul des tailles des images intermédiaires
        //ww et hh : plus grand que w et h pour conserver toute l'image
        int ww,hh;
        if(w<w_f){ww = w_f;} else {ww = w;}
        if(h<h_f){hh = h_f;} else {hh = h;}
        //dw et dh : différence de taille entre img_i et img_f
        int dw = (w_f-w)/2, dh = (h_f-h)/2; //x_f et x ont même parité
        //dwp et dhp leur partie positive (xp=max(0,x))
        int dwp = (dw<0) ? 0 : dw;
        int dhp = (dh<0) ? 0 : dh;
        //dwn et dhn leur partie négative (xn=max(0,-x))
        int dwn = (-dw<0) ? 0 : -dw;
        int dhn = (-dh<0) ? 0 : -dh;

        ///on place l'image originale au centre :
        double *img1 = malloc(3*ww*hh*sizeof(double));
        if(img1==NULL){
            printf("apply_affinite : img1 prend trop de place\n");
            exit(1);
        }
        int i_sym,j_sym;
        for(int l=0;l<3;l++){
            for(int i=0;i<ww;i++){
                for(int j=0;j<hh;j++){
                    i_sym = i-dwp;
                    while(i_sym<0 || i_sym>w-1){i_sym = (i_sym<0) ? -1-i_sym : 2*w-1-i_sym;}
                    j_sym = j-dhp;
                    while(j_sym<0 || j_sym>h-1){j_sym = (j_sym<0) ? -1-j_sym : 2*h-1-j_sym;}
                    img1[(i+ww*j)*3+l]=img[(i_sym+j_sym*w)*3+l];
                }
            }
        }

        ///prepare spline
        double *img2ij = malloc(3*sizeof(double));
        if(img2ij==NULL){
            printf("apply_affinite : img2ij prend trop de place\n");
            exit(1);
        }
        bool r = prepare_spline(img1,ww,hh,3,ORDER);
        if(!r){printf("apply_affinite, naive b-spline : Impossible de preparer la spline\n");exit(1);}

        ///evaluate spline
        double x,y,xx,yy,xc,yc;
        double *img2 = malloc(3*ww*hh*sizeof(double));
        if(img2==NULL){
            printf("apply_affinite : img2 prend trop de place\n");
            exit(1);
        }
        xc = ((double)ww-1.)/2.; yc = ((double)hh-1.)/2.;
        for(int i=0;i<ww;i++){
            for(int j=0;j<hh;j++){
                x = i; y = j;

                xx = affinity[0]*(x-xc)+affinity[1]*(y-yc)+xc;
                yy = affinity[3]*(x-xc)+affinity[4]*(y-yc)+yc;

                evaluate_spline_at(img2ij, img1, ww, hh, 3, ORDER, xx, yy);
                for(int l=0;l<3;l++){
                    img2[(i+ww*j)*3+l] = img2ij[l];
                }
            }
        }

        ///output
        for(int l=0;l<3;l++){
            for(int i=0;i<w_f;i++){
                for(int j=0;j<h_f;j++){
                    img_f[(i+j*w_f)*3+l] = img2[(i+dwn+(j+dhn)*ww)*3+l];
                }
            }
        }

        free(img2ij);
        free(img1);
        free(img2);
        }
        break;
    }
	return 0;
}
