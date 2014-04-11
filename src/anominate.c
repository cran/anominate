// R CMD SHLIB -I'/usr/include/R' -L'/usr/lib64/R/lib'  -lRlapack -lRblas nominate.c
// R CMD SHLIB nominate.c  -lRlapack -lRBlas

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#define MAXVOTES 100000
#define SLICE_W 8
#define SLICE_P 3

typedef struct{
    int row, col, vote;
//    int icpsr,period;
}cell;

/*  blockData is a generic structure used to pass information back and forth, mainly into the slicer
   'Data' is the raw data, 'ideal/yealoc/nayloc' are locations, 'length' is Data length
   'dim' is number of dimensions, 'paramDim' is dimension of parameter to be updated
   'beta/weight' are NOMINATE parameters.  'scaleParam' is the identifying parameter.
*/
typedef struct{
    cell **Data;
    double *ideal, *yealoc, *nayloc, *scaleParam;
    int length, dim, paramDim,nrow,ncol;
    double beta, alpha, weight;
} blockData;

typedef struct{
    cell ***rowData;
    double *ideal, *yealoc, *nayloc;
    int *rowLengths,nrow,ncol,dims;
    double weight,alpha,beta;
} betaBlock;

double slice(double (*fp)(double *, void *), double *init, void *ptr, double w, int p);
double slice_alpha(double (*fp)(double *, void *), double *init, void *ptr, double w, int p);
double nomLogLike(int vote, double *x, double *yea, double *nay, double beta, double weight, double alpha, int dims);
double legisLogLike(double *par, void *legis_ptr);
double nayLogLike(double *par, void *nay_ptr);
double yeaLogLike(double *par, void *yea_ptr);
void readDataFromVector(int *inputVector, cell **data, int *nvotes, int *nrow, int *ncol);
void formatData(cell *data, int nvotes, int nrow, int ncol, int **rowLengths, int **colLengths, cell ****rowData, cell ****colData);
void readInitsFromVector(double *initIdeal, double *initBill, double **idealpts, double **yealocs, double **naylocs, double **beta, double **alpha, int nrow, int ncol, int dims);
void writeDataOutput(double *output, int *op, int nrow, int ncol, double **idealpts, double **yealocs,
		     double **naylocs, double **beta, double **alpha, double *Sideal, double *Syea, double *Snay, int dims);
void updateYea(cell ****colData,int **colLengths,double **idealpts,double **yealocs,double **naylocs,double **beta,double **alpha,int nrow,int ncol,int dims,double *Syea);
void updateNay(cell ****colData,int **colLengths, double **idealpts,double **yealocs,double **naylocs,double **beta,double **alpha,int nrow,int ncol,int dims, double *Snay);
void updateLegis(cell ****rowData,int **rowLengths,double **idealpts,double **yealocs,double **naylocs,double **beta,double **alpha,int nrow,int ncol,int dims,double *Sideal);
void genScaleParams(int nrow, int ncol,int dims,double **idealpts, double **yealocs, double **naylocs,double *Syea,double *Snay,double *Sideal);
void sampleData(int Nsamples,int nrow,int ncol,int **rowLengths,int **colLengths,cell ****rowData,cell ****colData,double **idealpts,double **yealocs,double **naylocs, double **beta, double **alpha, int verbose, int thin, int dims, double *output, int constrain);
void freeData(int nrow, int ncol, cell **data, double **idealpts, double **yealocs, double **naylocs, int **rowLengths, int **colLengths, cell ****rowData, cell ****colData, double **beta, double **alpha);
void rwish(int v, double *S, int dims, double *output);
void riwish(int v, double *S, int dims, double *output);
double calcPrior(double *sample, int dims, double *scaleParam);

double betaLegisLL(double *par, void *legis_ptr);
double betaLogLike(double *par, void *beta_ptr);
void updateBeta(cell ****rowData,int **rowLengths,double **idealpts,double **yealocs,double **naylocs,double **beta,double **alpha,int nrow,int ncol,int dims);
void updateAlpha(cell ****rowData,int **rowLengths,double **idealpts,double **yealocs,double **naylocs,double **beta,double **alpha,int nrow,int ncol,int dims);
double alphaLogLike(double *par, void *alpha_ptr);
double alphaLegisLL(double *par, void *legis_ptr);
void Canominate(int *dataVector, double *initIdeal, double *initBill, double *output, int *thin, int *ncol, int *nrow, int *Nsamples, int *dims, int *verbose, int *constrain);



double slice(double (*fp)(double *, void *), double *init, void *ptr, double w, int p){
    double x,y,L,R;
    int K;
    short int flag;
    blockData *sliceData = (blockData *) ptr;

    y = fp(init, sliceData) - rexp(1);
    L = *init - w*runif(0,1);
    R = L + w;
    K = p;

    flag=0;
    if(y < fp(&L,sliceData)) flag=1;
    else if(y < fp(&R,sliceData)) flag=1;

    while(K>0 && flag==1) {
        if(runif(0,1)<0.5) {
            L = 2*L - R;
            if(y >=fp(&L,sliceData)) flag=0;
        }
        else {
            R = 2*R -L;
            if(y >=fp(&R,sliceData)) flag=0;
        }
        K--;
    }

    while(1){
        x = L + runif(0,1)*(R-L);
        if(fp(&x,sliceData)>y) break;
        if(*init > x) L = x;
        else R = x;
    }

    return(x);
}


double slice_alpha(double (*fp)(double *, void *), double *init, void *ptr, double w, int p){
    double x,y,L,R;
    int K;
    short int flag;
    blockData *sliceData = (blockData *) ptr;

    y = fp(init, sliceData) - rexp(1);
    L = *init - w*runif(0,1);
    R = L + w;
    K = p;
    if(L<0) L=0;
    if(R>1) R=1;

    flag=0;
    if(y < fp(&L,sliceData)) flag=1;
    else if(y < fp(&R,sliceData)) flag=1;

    while(K>0 && flag==1) {
        if(runif(0,1)<0.5) {
            L = 2*L - R;
	       if(L<0) L=0;
            if(y >=fp(&L,sliceData)) flag=0;
        }
        else {
            R = 2*R -L;
	       if(R>1) R=1;
            if(y >=fp(&R,sliceData)) flag=0;
        }
        K--;
    }

    while(1){
        x = L + runif(0,1)*(R-L);
        if(fp(&x,sliceData)>y) break;
        if(*init > x) L = x;
        else R = x;
    }

    return(x);
}



double nomLogLike(int vote, double *x, double *yea, double *nay, double beta, double weight, double alpha, int dims){

    int i;
    double nomyeaUtil, nomnayUtil, quadyeaUtil, quadnayUtil, nomdiffUtil, quaddiffUtil;
    double sqDistYea = 0, sqDistNay=0;

    for(i=0;i<dims;i++){
    sqDistYea += pow(x[i]-yea[i],2);
    sqDistNay += pow(x[i]-nay[i],2);
    }

    nomyeaUtil = beta*exp(-0.5*weight*weight*sqDistYea);
    nomnayUtil = beta*exp(-0.5*weight*weight*sqDistNay);
    nomdiffUtil = nomyeaUtil - nomnayUtil;

    quadyeaUtil = -0.5*beta*weight*weight*sqDistYea;
    quadnayUtil = -0.5*beta*weight*weight*sqDistNay;
    quaddiffUtil = quadyeaUtil - quadnayUtil;

    if(vote==1) return(pnorm(quaddiffUtil + alpha*(nomdiffUtil-quaddiffUtil), 0, 1, 1, 1));
    if(vote==0) return(pnorm( (quaddiffUtil + alpha*(nomdiffUtil-quaddiffUtil))*-1, 0, 1, 1, 1));

//    if(vote==0) return( pnorm(beta*(yeaUtil-nayUtil), 0, 1, 1, 1) );
//    if(vote==1) return( pnorm(beta*(nayUtil-yeaUtil), 0, 1, 1, 1) );

    return(0);
}

//Takes one legislator (pair of ideal pts), all roll call locations, gives log likelihood
double legisLogLike(double *par, void *legis_ptr){

    int i,j;
    double loglike=0.0, *yea, *nay, *templegis, prior;
    blockData *legisData = (blockData *)legis_ptr;

    yea = (double *) malloc((*legisData).dim * sizeof(double));
    nay = (double *) malloc((*legisData).dim * sizeof(double));
    templegis = (double *) malloc((*legisData).dim * sizeof(double));

    for(i=0; i<(*legisData).dim; i++) {
        templegis[i] = (*legisData).ideal[i];
    }
    templegis[(*legisData).paramDim -1] = *par;

    for(i=0; i<(*legisData).length; i++) {

        for(j=0;j<(*legisData).dim;j++){
            yea[j] = (*legisData).yealoc[(*(*legisData).Data[i]).col -1 + j*(*legisData).ncol];
            nay[j] = (*legisData).nayloc[(*(*legisData).Data[i]).col -1 + j*(*legisData).ncol];
        }

        loglike += nomLogLike((*(*legisData).Data[i]).vote, templegis, yea, nay, (*legisData).beta, (*legisData).weight, (*legisData).alpha, (*legisData).dim);
    }

    prior = calcPrior(templegis,(*legisData).dim,(*legisData).scaleParam);

    free(yea);
    free(nay);
    free(templegis);
    return(loglike - prior/2.0);

}

//Takes one pair yea/nay locations, all ideal pts, gives log likelihood
double nayLogLike(double *par, void *nay_ptr){

    int i,j;
    double loglike=0.0, *tempnay, *ideal, prior;
    blockData *nayData = (blockData *)nay_ptr;

    tempnay = (double *) malloc((*nayData).dim * sizeof(double));
    ideal = (double *) malloc((*nayData).dim * sizeof(double));

    for(i=0; i<(*nayData).dim; i++){
    tempnay[i] = (*nayData).nayloc[i];
    }
    tempnay[(*nayData).paramDim -1] = *par;

    for(i=0; i<(*nayData).length; i++){

        for(j=0;j<(*nayData).dim;j++){
            ideal[j] = (*nayData).ideal[(*(*nayData).Data[i]).row -1 + j*(*nayData).nrow];
        }
        loglike += nomLogLike((*(*nayData).Data[i]).vote, ideal,(*nayData).yealoc,tempnay,(*nayData).beta,(*nayData).weight,(*nayData).alpha,(*nayData).dim);
    }

    prior = calcPrior(tempnay,(*nayData).dim,(*nayData).scaleParam);

    free(tempnay);
    free(ideal);
    return(loglike - prior/2.0);
}


//Takes one pair yea/nay locations, all ideal pts, gives log likelihood
double yeaLogLike(double *par, void *yea_ptr){

    int i,j;
    double loglike=0.0, *tempyea, *ideal, prior;
    blockData *yeaData = (blockData *)yea_ptr;

    tempyea = (double *) malloc((*yeaData).dim * sizeof(double));
    ideal = (double *) malloc((*yeaData).dim * sizeof(double));

    for(i=0; i<(*yeaData).dim; i++){
    tempyea[i] = (*yeaData).yealoc[i];
    }
    tempyea[(*yeaData).paramDim -1] = *par;

    for(i=0; i<(*yeaData).length; i++){

        for(j=0; j<(*yeaData).dim; j++){
            ideal[j] = (*yeaData).ideal[(*(*yeaData).Data[i]).row -1 + j*(*yeaData).nrow];
        }
        loglike += nomLogLike((*(*yeaData).Data[i]).vote, ideal,tempyea,(*yeaData).nayloc,(*yeaData).beta,(*yeaData).weight,(*yeaData).alpha,(*yeaData).dim);
    }

    prior = calcPrior(tempyea,(*yeaData).dim,(*yeaData).scaleParam);

    free(tempyea);
    free(ideal);
    return(loglike - prior/2.0);

}

void readDataFromVector(int *inputVector, cell **data, int *nvotes, int *nrow, int *ncol) {
    int idx,i,j,nRow,nCol,nVotes;
    cell *p = (cell *) malloc((*nrow)*(*ncol)*sizeof(cell));

    nRow=*nrow;
    nCol=*ncol;
    Rprintf("nCol=%i\n",nCol);
    Rprintf("nRow=%i\n",nRow);

    Rprintf("ANOM::Reading roll call data ...\n");

    idx=0;
    nVotes = 0;
    for (i=0;i<nCol;i++) {
      for (j=0;j<nRow;j++) {
	if (inputVector[idx] != -1) {
	  p[nVotes].row = j+1;
	  p[nVotes].col = i+1;
	  p[nVotes].vote = inputVector[idx];
	  nVotes++;
	}
        idx++;
      }
    }
    *data = (cell *) realloc(p, nVotes * sizeof(cell));
    *nvotes = nVotes;
    Rprintf("\nAllocation OK, %i votes allocated.\n", nVotes);
    Rprintf("ANOM::Done reading rc data...\n");
}


void formatData(cell *data, int nvotes, int nrow, int ncol, int **rowLengths, int **colLengths, cell ****rowData, cell ****colData) {

    int i;
    int *rLengths, *cLengths, *rIndex, *cIndex;
    cell ***rData, ***cData;

    rLengths = calloc(nrow, sizeof(int));
    cLengths = calloc(ncol, sizeof(int));
    rIndex = calloc(nrow, sizeof(int));
    cIndex = calloc(ncol, sizeof(int));
    rData = (cell ***) malloc(nrow * sizeof(cell **));
    cData = (cell ***) malloc(ncol * sizeof(cell **));

    for(i=0; i<nvotes; i++) {
        rLengths[data[i].row - 1]++;
        cLengths[data[i].col - 1]++;
    }

    for(i=0; i<nrow; i++) rData[i] = (cell **) malloc(rLengths[i] * sizeof(cell *));
    for(i=0; i<ncol; i++) cData[i] = (cell **) malloc(cLengths[i] * sizeof(cell *));
    for(i=0; i<nvotes; i++) {
        rData[data[i].row-1][rIndex[data[i].row-1]] = &data[i];
        rIndex[data[i].row-1]++;
        cData[data[i].col-1][cIndex[data[i].col-1]] = &data[i];
        cIndex[data[i].col-1]++;

    }

    *rowLengths=rLengths;
    *colLengths=cLengths;
    *rowData=rData;
    *colData=cData;
    free(rIndex);
    free(cIndex);
}

void readInitsFromVector(double *inputIdeal,
			 double *inputBill,
			 double **idealpts, double **yealocs, double **naylocs,
			 double **beta, double **alpha, int nrow, int ncol, int dims){

    int i;
    double *idealptr,*yeaptr,*nayptr,*betaptr,*alphaptr;

    idealptr = (double *) malloc(dims * nrow * sizeof(double));
    yeaptr = (double *) malloc(dims * ncol * sizeof(double));
    nayptr = (double *) malloc(dims * ncol * sizeof(double));
    betaptr = (double *) malloc(sizeof(double));
    alphaptr = (double *) malloc(sizeof(double));

    for(i=0;i<(nrow*dims);i++) idealptr[i]=inputIdeal[i];
    Rprintf("\n%i legislator start coordinates read.\n", i);

    for(i=0;i<(ncol*dims);i++) {
      yeaptr[i] = inputBill[2*i];
      nayptr[i] = inputBill[2*i+1];
      //Rprintf("yea[%i] = %7.4f\n",i,inputBill[2*i]);
    }
    Rprintf("\n%i bill start coordinate pairs read.\n\n", i);

    betaptr[0] = 10.0;
    alphaptr[0] = 0.7;

    *idealpts = idealptr;
    *yealocs = yeaptr;
    *naylocs = nayptr;
    *beta = betaptr;
    *alpha = alphaptr;
}



void writeDataOutput(double *output, int *op, int nrow, int ncol, double **idealpts, double **yealocs,
		     double **naylocs, double **beta, double **alpha, double *Sideal, double *Syea, double *Snay, int dims) {
  int i;
  int outputpos = *op;
  output[outputpos++] = (*beta)[0];
  output[outputpos++] = (*alpha)[0];
  for(i=0;i<(nrow*dims);i++) output[outputpos++] = (*idealpts)[i];
  for(i=0;i<(ncol*dims);i++) output[outputpos++] = (*yealocs)[i];
  for(i=0;i<(ncol*dims);i++) output[outputpos++] = (*naylocs)[i];
  for(i=0;i<(dims*dims);i++) output[outputpos++] = Syea[i];
  for(i=0;i<(dims*dims);i++) output[outputpos++] = Snay[i];
  *op = outputpos;
}


void updateYea(cell ****colData,int **colLengths,double **idealpts,double **yealocs,double **naylocs,double **beta,double **alpha,int nrow,int ncol,int dims,double *Syea){

    int i,j;
    double *yea,*nay,par;
    blockData block;

    yea = (double *) malloc(dims * sizeof(double));
    nay = (double *) malloc(dims * sizeof(double));

    block.beta = (*beta)[0];
    block.alpha = (*alpha)[0];
    block.weight = 0.5;
    block.ideal = *idealpts;
    block.dim = dims;
    block.scaleParam = Syea;
    block.nrow = nrow;
    block.ncol = ncol;

    for(i=0;i<ncol;i++){

        for(j=0;j<dims;j++){
            yea[j] = (*yealocs)[i + j*ncol];
            nay[j] = (*naylocs)[i + j*ncol];
        }

        block.yealoc = yea;
        block.nayloc = nay;
        block.Data = (*colData)[i];
        block.length = (*colLengths)[i];

        for(j=0;j<dims;j++){

            par = (*yealocs)[i + j*ncol];
            block.paramDim = j+1;
            (*yealocs)[i + j*ncol] = slice(yeaLogLike,&par,&block,SLICE_W,SLICE_P);
        }
    }

    free(yea);
    free(nay);

}


void updateNay(cell ****colData,int **colLengths, double **idealpts,double **yealocs,double **naylocs,double **beta,double **alpha,int nrow,int ncol,int dims, double *Snay){

    int i,j;
    double *yea,*nay,par;
    blockData block;

    yea = (double *) malloc(dims * sizeof(double));
    nay = (double *) malloc(dims * sizeof(double));

    block.beta = (*beta)[0];
    block.alpha = (*alpha)[0];
    block.weight = 0.5;
    block.ideal = *idealpts;
    block.dim = dims;
    block.scaleParam = Snay;
    block.nrow = nrow;
    block.ncol = ncol;

    for(i=0;i<ncol;i++){

        for(j=0;j<dims;j++){
            yea[j] = (*yealocs)[i + j*ncol];
            nay[j] = (*naylocs)[i + j*ncol];
        }

        block.yealoc = yea;
        block.nayloc = nay;
        block.Data = (*colData)[i];
        block.length = (*colLengths)[i];

        for(j=0;j<dims;j++){

            par = (*naylocs)[i + j*ncol];
            block.paramDim = j+1;
            (*naylocs)[i + j*ncol] = slice(nayLogLike,&par,&block,SLICE_W,SLICE_P);
        }
    }

    free(yea);
    free(nay);
}


void updateLegis(cell ****rowData,int **rowLengths,double **idealpts,double **yealocs,double **naylocs,double **beta,double **alpha,int nrow,int ncol,int dims,double *Sideal){

    int i,j;
    double *ideal,par;
    blockData block;

    ideal = (double *) malloc(dims * sizeof(double));

    block.beta = (*beta)[0];
    block.alpha = (*alpha)[0];
    block.weight = 0.5;
    block.yealoc = *yealocs;
    block.nayloc = *naylocs;
    block.dim = dims;
    block.scaleParam = Sideal;
    block.nrow = nrow;
    block.ncol = ncol;

    for(i=0;i<nrow;i++){

        for(j=0;j<dims;j++){
            ideal[j] = (*idealpts)[i + j*nrow];
        }

        block.ideal = ideal;
        block.Data = (*rowData)[i];
        block.length = (*rowLengths)[i];

        for(j=0;j<dims;j++){
            par = (*idealpts)[i + j*nrow];
            block.paramDim = j+1;
            (*idealpts)[i + j*nrow] = slice(legisLogLike,&par,&block,SLICE_W,SLICE_P);
        }
    }

    free(ideal);

}

void genScaleParams(int nrow, int ncol,int dims,double **idealpts, double **yealocs, double **naylocs,double *Syea,double *Snay,double *Sideal){

    char trans = 't', notrans = 'n';
    double alpha=1.0, beta=0.0;
    double *S;

    S = (double *) malloc(dims * dims * sizeof(double));

    dgemm_(&trans,&notrans,&dims,&dims,&nrow,&alpha,*idealpts,&nrow,*idealpts,&nrow,&beta,S,&dims);         // Calculate outer product
    riwish(nrow-1,S,dims,Sideal);                                                                          //  Generate wishart draw
    dgemm_(&trans,&notrans,&dims,&dims,&ncol,&alpha,*yealocs,&ncol,*yealocs,&ncol,&beta,S,&dims);
    riwish(ncol-1,S,dims,Syea);
    dgemm_(&trans,&notrans,&dims,&dims,&ncol,&alpha,*naylocs,&ncol,*naylocs,&ncol,&beta,S,&dims);
    riwish(ncol-1,S,dims,Snay);

    free(S);

}


void sampleData(int Nsamples,int nrow,int ncol,int **rowLengths,int **colLengths,cell ****rowData,cell ****colData,double **idealpts,double **yealocs,double **naylocs, double **beta, double **alpha, int verbose, int thin, int dims, double *output, int constrain){

    int i;
    double *Syea,*Snay,*Sideal;

    int outputpos = 0; // keep track of current position in the output container

    Syea = (double *) malloc(dims * dims * sizeof(double));
    Snay = (double *) malloc(dims * dims * sizeof(double));
    Sideal = (double *) malloc(dims * dims * sizeof(double));

    writeDataOutput(output,&outputpos,nrow,ncol,idealpts,yealocs,naylocs,beta,alpha,Sideal,Syea,Snay,dims);
    for(i=1;i<Nsamples;i++){
        genScaleParams(nrow,ncol,dims,idealpts,yealocs,naylocs,Syea,Snay,Sideal);
	//Rprintf("here 1\n"); R_FlushConsole();
        updateYea(colData,colLengths,idealpts,yealocs,naylocs,beta,alpha,nrow,ncol,dims,Syea);
	//Rprintf("here 2\n"); R_FlushConsole();
        updateNay(colData,colLengths,idealpts,yealocs,naylocs,beta,alpha,nrow,ncol,dims,Snay);
	//Rprintf("here 3\n"); R_FlushConsole();
        updateLegis(rowData,rowLengths,idealpts,yealocs,naylocs,beta,alpha,nrow,ncol,dims,Sideal);
	//Rprintf("here 4\n"); R_FlushConsole();
        updateBeta(rowData,rowLengths,idealpts,yealocs,naylocs,beta,alpha,nrow,ncol,dims);
	//Rprintf("here 5\n"); R_FlushConsole();
        if(constrain !=1) updateAlpha(rowData,rowLengths,idealpts,yealocs,naylocs,beta,alpha,nrow,ncol,dims);
	//Rprintf("here 6\n"); R_FlushConsole();
        if(i % thin==0) writeDataOutput(output,&outputpos,nrow,ncol,idealpts,yealocs,naylocs,beta,alpha,Sideal,Syea,Snay,dims);
//        if(i % thin==0) writeData(nrow,ncol,idealpts,yealocs,naylocs,beta,alpha,Sideal,Syea,Snay,dims);
	//Rprintf("here 7\n"); R_FlushConsole();
        if(i % verbose==0) Rprintf("\tSample %i complete...\n", i);
	//Rprintf("here 8\n"); R_FlushConsole();
        R_FlushConsole();
    }

    free(Syea);
    free(Snay);
    free(Sideal);

}

void freeData(int nrow, int ncol, cell **data, double **idealpts, double **yealocs, double **naylocs, int **rowLengths, int **colLengths, cell ****rowData, cell ****colData, double **beta, double **alpha){

    int i;
    free(*data);
    free(*beta);
    free(*alpha);
    free(*idealpts);
    free(*yealocs);
    free(*naylocs);
    for(i=0;i<nrow;i++) free((*rowData)[i]);
    for(i=0;i<ncol;i++) free((*colData)[i]);
    free(*rowData);
    free(*colData);
    free(*rowLengths);
    free(*colLengths);

}


// warning: this function will overwrite S!!
void rwish(int v, double *S, int dims, double *output) {

    int i,j,info,*ipiv;
    double *temp, *Scopy;
    char uplo = 'U',trans = 't', notrans = 'n';
    double alpha=1.0, beta=0.0;

    ipiv = (int *) malloc(dims * sizeof(ipiv));
    temp = (double *) malloc(dims * dims * sizeof(double));
    Scopy = (double *) malloc(dims * dims * sizeof(double));

    for(i=0;i<(dims*dims);i++) Scopy[i] = S[i];                         //copy S, so it leaves unchanged
    dpotrf_(&uplo,&dims,Scopy,&dims,&info);                              // take Cholesky of S

    for(j=0;j<dims;j++){
        for(i=0;i<=j;i++){
            if(i==j) output[j*dims + i] = sqrt(rchisq(v--));   // sqrt(~chisq) on Z's diagnonals
            else {                                              // case of i<j, in lower triangle
                Scopy[i*dims + j] = 0;  // zero out S lower triangle to complete Cholesky
                output[i*dims + j] = 0;  // zero out Z's lower triangle
                output[j*dims + i] = rnorm(0.0,1.0); // normals on Z's upper triangle
            }
        }
    }

/*
    // Note: this gets replaced at the end with the line commented out above
    // This inefficiency is only included so the seeds are EXACTLY comparable to R results given identical seeds
    for(j=0;j<dims;j++){
        for(i=0;i<j;i++){
            if(i<j) output[j*dims + i] = rnorm(0.0,1.0); // normals on Z's upper triangle
        }
    }
*/

    dgemm_(&notrans,&notrans,&dims,&dims,&dims,&alpha,output,&dims,Scopy,&dims,&beta,temp,&dims);       // generates Z %*% S intermediate matrix
    dgemm_(&trans,&notrans,&dims,&dims,&dims,&alpha,temp,&dims,temp,&dims,&beta,output,&dims);   // generates Wishart of S via crossproduct

    free(ipiv);
    free(temp);

}


void riwish(int v, double *S, int dims, double *output){

    int info, *ipiv, i, j;
    double *Sinv, *temp, *Scopy;

    ipiv = (int *) malloc(dims * sizeof(ipiv));
    Sinv = (double *) calloc(dims * dims, sizeof(double));
    temp = (double *) malloc(dims * dims * sizeof(double));
    Scopy = (double *) malloc(dims * dims * sizeof(double));

    //make output and Sinv identity matrices
    for(j=0;j<dims;j++){
        for(i=0;i<=j;i++){
            if(i==j) {
                output[j*dims + i] = 1.0;
                Sinv[j*dims + i] = 1.0;
            }
            else {
                output[j*dims + i] = 0.0;
                output[i*dims + j] = 0.0;
            }
        }
    }

    for(i=0;i<(dims*dims);i++) Scopy[i] = S[i];                         //copy S
    dgesv_(&dims,&dims,Scopy,&dims,ipiv,Sinv,&dims,&info);              //invert S
    rwish(v,Sinv,dims,temp);                                            //generate Wishart
    dgesv_(&dims,&dims,temp,&dims,ipiv,output,&dims,&info);             //invert the Wishart

    free(ipiv);
    free(Sinv);
    free(temp);

}

double calcPrior(double *sample, int dims, double *scaleParam){

    // sample is 1 x dims. Is either idealpts, yealocs, or naylocs for a single row/column
    // scaleParam  is dims * dims
    // temp is 1 x dims, results from sample %*% scaleParam

    double *temp, alpha=1.0, beta=0.0, prior;
    char notrans='n';
    int one=1;

    temp = (double *) malloc(dims * sizeof(double));

    dgemm_(&notrans,&notrans,&one,&dims,&dims,&alpha,sample,&one,scaleParam,&dims,&beta,temp,&one);  //temp = sample %*% scaleParam
    dgemm_(&notrans,&notrans,&one,&one,&dims,&alpha,temp,&one,sample,&dims,&beta,&prior,&one);  //prior = temp %*% sample

    free(temp);
    return(prior);

}


//Takes one legislator (pair of ideal pts), all roll call locations, gives log likelihood
//Essentially legisLogLike, but without a prior.  Used in betaLogLike
double betaLegisLL(double *par, void *legis_ptr){

    int i,j;
    double loglike=0.0, *yea, *nay, *templegis;
    blockData *legisData = (blockData *)legis_ptr;

    yea = (double *) malloc((*legisData).dim * sizeof(double));
    nay = (double *) malloc((*legisData).dim * sizeof(double));
    templegis = (double *) malloc((*legisData).dim * sizeof(double));

    for(i=0; i<(*legisData).dim; i++) {
        templegis[i] = (*legisData).ideal[i];
    }

    for(i=0; i<(*legisData).length; i++) {

        for(j=0;j<(*legisData).dim;j++){
            yea[j] = (*legisData).yealoc[(*(*legisData).Data[i]).col -1 + j*(*legisData).ncol];
            nay[j] = (*legisData).nayloc[(*(*legisData).Data[i]).col -1 + j*(*legisData).ncol];
        }

        loglike += nomLogLike((*(*legisData).Data[i]).vote, templegis, yea, nay, *par, (*legisData).weight, (*legisData).alpha, (*legisData).dim);
    }

    free(yea);
    free(nay);
    free(templegis);
    return(loglike);
}


double betaLogLike(double *par, void *beta_ptr){

    int i, j;
    double loglike=0.0, *ideal,temppar;
    blockData block;

    temppar = *par;
    if(temppar < 0 ) return(-1e307);

    betaBlock *betaData = (betaBlock *)beta_ptr;
    block.beta = temppar;
    block.weight = (*betaData).weight;
    block.yealoc = (*betaData).yealoc;
    block.nayloc = (*betaData).nayloc;
    block.nrow = (*betaData).nrow;
    block.ncol = (*betaData).ncol;
    block.dim = (*betaData).dims;
    block.alpha = (*betaData).alpha;

    ideal = (double *) malloc((*betaData).dims * sizeof(double));

    for(i=0;i<block.nrow;i++){

            for(j=0;j<(*betaData).dims;j++){
                ideal[j] = (*betaData).ideal[i + j*block.nrow];
            }
            block.ideal = ideal;
            block.Data = (*betaData).rowData[i];
            block.length = (*betaData).rowLengths[i];

            loglike += betaLegisLL(&temppar, &block);
    }

    free(ideal);
    return(loglike);

}

void updateBeta(cell ****rowData,int **rowLengths,double **idealpts,double **yealocs,double **naylocs,double **beta,double **alpha,int nrow,int ncol,int dims){

    double temp;
    betaBlock block;

    temp = (*beta)[0];

    block.alpha = (*alpha)[0];
    block.weight = 0.5;
    block.yealoc = *yealocs;
    block.nayloc = *naylocs;
    block.dims = dims;
    block.nrow = nrow;
    block.ncol = ncol;
    block.ideal = *idealpts;
    block.rowData = *rowData;
    block.rowLengths = *rowLengths;

    (*beta)[0] = slice(betaLogLike,&temp,&block,SLICE_W,SLICE_P);

}


void updateAlpha(cell ****rowData,int **rowLengths,double **idealpts,double **yealocs,double **naylocs,double **beta,double **alpha,int nrow,int ncol,int dims){

    double temp;
    betaBlock block;

    temp = (*alpha)[0];
    block.beta = (*beta)[0];
    block.weight = 0.5;
    block.yealoc = *yealocs;
    block.nayloc = *naylocs;
    block.dims = dims;
    block.nrow = nrow;
    block.ncol = ncol;
    block.ideal = *idealpts;
    block.rowData = *rowData;
    block.rowLengths = *rowLengths;

    (*alpha)[0]= slice_alpha(alphaLogLike,&temp,&block,SLICE_W,SLICE_P);

}

double alphaLogLike(double *par, void *alpha_ptr){

    int i, j;
    double loglike=0.0, *ideal,temppar;
    blockData block;

    temppar = *par;
    // if(temppar < 0 ) return(-1e307);
    //if(temppar > 1 ) return(-1e307);

    betaBlock *alphaData = (betaBlock *)alpha_ptr;
    block.beta = temppar;
    block.weight = (*alphaData).weight;
    block.yealoc = (*alphaData).yealoc;
    block.nayloc = (*alphaData).nayloc;
    block.nrow = (*alphaData).nrow;
    block.ncol = (*alphaData).ncol;
    block.dim = (*alphaData).dims;
    block.beta = (*alphaData).beta;

    ideal = (double *) malloc((*alphaData).dims * sizeof(double));

    for(i=0;i<block.nrow;i++){

            for(j=0;j<(*alphaData).dims;j++){
                ideal[j] = (*alphaData).ideal[i + j*block.nrow];
            }
            block.ideal = ideal;
            block.Data = (*alphaData).rowData[i];
            block.length = (*alphaData).rowLengths[i];

            loglike += alphaLegisLL(&temppar, &block);
    }

    free(ideal);
    return(loglike);

}


//Takes one legislator (pair of ideal pts), all roll call locations, gives log likelihood
//Essentially legisLogLike, but without a prior.  Used in betaLogLike
double alphaLegisLL(double *par, void *legis_ptr){

    int i,j;
    double loglike=0.0, *yea, *nay, *templegis;
    blockData *legisData = (blockData *)legis_ptr;

    yea = (double *) malloc((*legisData).dim * sizeof(double));
    nay = (double *) malloc((*legisData).dim * sizeof(double));
    templegis = (double *) malloc((*legisData).dim * sizeof(double));

    for(i=0; i<(*legisData).dim; i++) {
        templegis[i] = (*legisData).ideal[i];
    }

    for(i=0; i<(*legisData).length; i++) {

        for(j=0;j<(*legisData).dim;j++){
            yea[j] = (*legisData).yealoc[(*(*legisData).Data[i]).col -1 + j*(*legisData).ncol];
            nay[j] = (*legisData).nayloc[(*(*legisData).Data[i]).col -1 + j*(*legisData).ncol];
        }

        loglike += nomLogLike((*(*legisData).Data[i]).vote, templegis, yea, nay, (*legisData).beta, (*legisData).weight, *par, (*legisData).dim);
    }

    free(yea);
    free(nay);
    free(templegis);
    return(loglike);
}

 void Canominate(int *dataVector, double *initIdeal, double *initBill,
		 double *output, int *thin, int *ncol, int *nrow,
	 	 int *Nsamples, int *dims, int *verbose, int *constrain) {
    cell *data;
    cell ***rowData, ***colData;
    int nvotes;
    int *rowLengths, *colLengths;
    double *idealpts, *yealocs, *naylocs, *beta, *alpha;

    readDataFromVector(dataVector,&data,&nvotes,nrow,ncol);
    formatData(data,nvotes,*nrow,*ncol,&rowLengths,&colLengths,&rowData,&colData);
    readInitsFromVector(initIdeal,initBill,&idealpts,&yealocs,&naylocs,&beta,&alpha,*nrow,*ncol,*dims);
    if(*constrain==1) *alpha=1;
    sampleData(*Nsamples,*nrow,*ncol,&rowLengths,&colLengths,&rowData,&colData,&idealpts,
	       &yealocs,&naylocs,&beta,&alpha,*verbose,*thin,*dims,output,*constrain);
    freeData(*nrow,*ncol,&data,&idealpts,&yealocs,&naylocs,&rowLengths,&colLengths,&rowData,&colData,&beta,&alpha);
}
