package org.moap.lbrs.mf;

import org.jblas.DoubleMatrix;

import java.util.Date;

/**
 * Created with IntelliJ IDEA.
 * User: igobrilhante
 * Date: 29/06/13
 * Time: 01:33
 * To change this template use File | Settings | File Templates.
 */
public class MatrixFactorization {

    public static DoubleMatrix P;
    public static DoubleMatrix Q;

    public static void main(String[] args){

        int k = 100;

        double[][] R = {
                {4,0,0,1,1,0,0},
                {5,5,4,0,0,0,0},
                {0,0,0,2,4,5,0},
                {0,3,0,0,0,0,3}
                    };

        DoubleMatrix matrixR = DoubleMatrix.rand(1000,3000);

        DoubleMatrix vectorP = DoubleMatrix.rand(matrixR.getRows(),k);
        DoubleMatrix vectorQ = DoubleMatrix.rand(matrixR.getColumns(),k);

        Date d1 = new Date();
        simpleFactorization(matrixR, vectorP, vectorQ, k, 5000, 0.0002, 0.02);
        Date d2 = new Date();
        System.out.println((d2.getTime()-d1.getTime()));

    }

    public static void simpleFactorization(DoubleMatrix matrixR, DoubleMatrix vectorP, DoubleMatrix vectorQ, int nFactors, int steps, double alpha, double beta){
        DoubleMatrix vectorQT = vectorQ.transpose();

        int nRows = matrixR.getRows();
        int nCols = matrixR.getColumns();

        for(int s=0;s<steps;s++){

            for(int i=0;i<nRows;i++){

               for(int j=0;j<nCols;j++){

                     if(matrixR.get(i,j) > 0){
                         double errorIJ = matrixR.get(i,j) - (vectorP.getRow(i).mmul(vectorQT.getColumn(j))).get(0);

                         // Update P and Q
                         for(int k=0;k<nFactors;k++ ){

                             double newPValue = vectorP.get(i,k) + alpha * (2*errorIJ * vectorQT.get(k,j) - beta * vectorP.get(i,k));
                             vectorP.put(i,k,newPValue);
                             
                             double newQValue = vectorQT.get(k,j) + alpha * (2*errorIJ * vectorP.get(i,k) - beta * vectorQT.get(k,j));
                             vectorQT.put(k,j,newQValue);
                         }
                     }
               }
            }
            double error = 0.0;

            // Compute the error
            for(int i=0;i<nRows;i++){
                for(int j=0;j<nCols;j++){
                    if(matrixR.get(i,j) > 0){
                        double aux = (vectorP.getRow(i).mmul(vectorQT.getColumn(j))).get(0);
                        error = error  + Math.pow(matrixR.get(i,j) - aux,2);
                        for(int k=0;k<nFactors;k++ ){
                            error = error + (beta/2) * (Math.pow(vectorP.get(i,k),2) + Math.pow(vectorQT.get(k,j),2) );
                        }
                    }
                }
            }

            // Tolerance ratio
            if(error < 0.001){
               break;
            }
        }
        P = vectorP;
        Q = vectorQT;
    }

}
