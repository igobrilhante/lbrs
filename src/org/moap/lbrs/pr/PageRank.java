package org.moap.lbrs.pr;

import gnu.trove.map.hash.TIntDoubleHashMap;
import org.jblas.DoubleMatrix;

import java.util.Arrays;
import java.util.Date;

/**
 * Created with IntelliJ IDEA.
 * User: igobrilhante
 * Date: 29/06/13
 * Time: 02:16
 * To change this template use File | Settings | File Templates.
 */
public class PageRank {

    private static double beta = 0.85;

    public static void main(String[] args){

        double[][] d = {
                {0,0.5,0,0},
                {1/3d,0,0,0.5},
                {1/3d,0,1,0.5},
                {1/3d,0.5,0,0}};
        int n = 10000;
        DoubleMatrix m = DoubleMatrix.rand(n,n);

        Date d1 = new Date();
        compute(m);
        Date d2 = new Date();
        System.out.println((d2.getTime()-d1.getTime()));

    }

    public static void compute(DoubleMatrix transitionMatrix){
        int nodeCount =  transitionMatrix.rows;

        DoubleMatrix pageRank = DoubleMatrix.ones(nodeCount,1).div(nodeCount);
        DoubleMatrix vector   = DoubleMatrix.ones(nodeCount,1).div(nodeCount).mul(1 - beta);

//        System.out.println("0 " + pageRank);

        /*
            PR = PR +
         */
        int count = 1;
        for(;;){
            pageRank = (transitionMatrix.mmul(pageRank)).mmul(beta).add(vector);
            count++;
            if(count == 100){
                break;
            }
        }

//        System.out.println(pageRank);






    }

}
