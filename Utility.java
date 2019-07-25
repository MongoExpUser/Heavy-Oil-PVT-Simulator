/* @License Starts
 *
 * Copyright © 2002 - present. MongoExpUser
 *
 * License: MIT - See: https://github.com/MongoExpUser/Heavy-Oil-PVT-Simulator/blob/master/LICENSE
 *
 * @License Ends
 
 *
 *
 * ...Ecotert's Utility.java (released as open-source under MIT License) implements:
 *
 *
 * Relevant general utility methods required as input into PREos.java (PVT simulator) e.g.
 *
 *
 */


package eosPVT;


public class Utility
{
  
 //contructor
 public void Utility(){}
 
 
  //a utility method for calculating complementary error function (erfc)
  double erfc(double x){
    double  coefficients = (1 + 0.0705230784*x + 0.0422820123*x*x + 0.0092705272*x*x*x + 0.0001520143*x*x*x*x + 0.0002765672*x*x*x*x*x + 0.0000430638*x*x*x*x*x*x);
    double value = Math.pow(coefficients, -16);
    return value; //Source: Butler, R.M. (1997). pgs. 35 & 36. error in computation  <= Math.abs(3E-7)
  }

  //a utility method for calculating error function (erf)
  double erf(double x){
    double  coefficients = (1 + 0.0705230784*x + 0.0422820123*x*x + 0.0092705272*x*x*x + 0.0001520143*x*x*x*x + 0.0002765672*x*x*x*x*x + 0.0000430638*x*x*x*x*x*x);
    double value = ( 1 - (Math.pow(coefficients, -16)) );
    return value; //Source: Butler, R.M. (1997). pgs 35 & 36. error in computation  <= Math.abs(3E-7)
  }
  
  
  //a utility method for finding the maximum (largest) of the 3 roots of an EOS = vapour z-factor
  double vapourCompressibilty(double [] A){
     double Value = maximumOfArrayValues(A);
     return Value;
   }

  //a utility method for finding the minimum (smallest) positive of the 3 roots of an EOS = liquid z-factor
  double liquidCompressibilty(double [] A){
     double [] B = new double[A.length];

     for(int i = 0; i < A.length; i++){
       if(A[i] > 0){
         B[i] = A[i];
       }
       else if(A[i] <= 0){
         B[i] = 1000;  //just any big value
       }
     }
     return minimumOfArrayValues(B);
   }

  //a utility method  for calculating natural log
  double nl(double value){
    return Math.log(value)/Math.log(2.718281828459045);
  }

  //a method for Newton-Jacobi solution method
  double [] newtonJacobi(double [] A, double [] B){
   int j = A.length;
   double [] E = new double[j];
   for(int i = 0; i < j; i++){
      E[i] = A[i]/B[i];
   }
   return E;
 }

  //implementation of thomas algorithm method for solving tri-diagonal matrix
  double [] thomasAlgorithmSolution (double upperElements [], double diagonalElements [], double lowerElements [], double rhsElements []){

    int matrixSize = diagonalElements.length; //For an m x n tri-diagonal matrix, number of rows (m) = number of cols n = (matrix size)

    //Arrays to hold (primary) main computed coefficients, secondary computed coefficients and solutions, x
    double A [] = new double [matrixSize]; 		//Primary
    double B [] = new double [matrixSize]; 		//Primary
    double C [] = new double [matrixSize]; 		//Primary
    double D [] = new double [matrixSize]; 		//Primary
    double W [] = new double [matrixSize]; 		//Secondary
    double G [] = new double [matrixSize]; 		//Secondary
    double E [] = new double [matrixSize]; 		//Solution
    double pivot = 0;                                   //pivoting value

    //Check for zero pivot element and Give a message box to signify exception
    if(diagonalElements[0] == 0){
      System.out.println("Zero Pivot Error 1 ! Change Time Step or Grid Number to Rectify.");
    }

    for(int i = 0; i < matrixSize; i++){
       //Primary
       if(i == 0)
         A[i] = 0;
       else
         A[i] = lowerElements[i - 1];
         B[i] = diagonalElements[i];
       if(i == (matrixSize - 1))
         C[i] = 0;
       else
         C[i] = upperElements[i];
         D[i] = rhsElements[i];

       //Secondary
       if(i == 0){
         W[i] = (C[i] / B[i]);
         G[i] = (D[i] / B[i]);
       }
       else if(i > 0){

       //Pivot condition testing
       pivot = (B[i] - (A[i] * W[i-1]));
         if(pivot == 0){
           //Give a message box to signify exception
           JOptionPane.showMessageDialog(null, "Zero Pivot Error 2 ! Change Time Step or Grid Number to Rectify.");
         }

         W[i] = (C[i] / (B[i] - (A[i] * W[i-1])));
         G[i] = ((D[i] - (A[i]*G[i-1])) / (B[i] - (A[i] * W[i-1])));
       }
     }

     //Solution
     for(int i = 0; i < matrixSize; i++){
       int j = (matrixSize - 1) - i;
       if(i == 0){
         E[j] = G[j];
       }
       else if(i > 0){
         E[j] = (G[j] - (W[j] * E[j+1]));
       }
     }
   return E;
 }

 //a utility method for copying array to array
 double [] copy(double [] A){
    int j = A.length;
    double [] E = new double [j];
    for(int i = 0; i < j; i++){
      E[i] = A[i];
    }
    return E;
  }

  //a utility method for copying array of integer values to array of integer values
 int [] copyInt(int [] A){
    int j = A.length;
    int [] E = new int [j];
    for(int i = 0; i < j; i++){
      E[i] = A[i];
    }
    return E;
  }

  //a utility method for summing element of an array
 double sumArray(double [] A){
    int j = A.length;
    double sum = 0;

    for(int i = 0; i < j; i++){
      sum = sum + A[i];
    }

    return sum;
  }

  //a utility for printing array elements to system prompt
 void printArrayToSystem(double [] A){
    for(int i = 0; i < A.length; i++){
          System.out.println( (i + 1)  +  ": "+ A[i] );
      }
  }

 
  //a utility method for finding the maximum value from an array of values
  double maximumOfArrayValues(double [] A){
     double staticMaximumValue = 0;
     double dynamicMaximumValue = 0;
     for(int i = 0; i < (A.length - 1); i++){
       if(i == 0){
         if(A[i] >= A[i+1]){
           dynamicMaximumValue = A[i];
           staticMaximumValue = dynamicMaximumValue;
         }
         else{
           dynamicMaximumValue = A[i+1];
           staticMaximumValue = dynamicMaximumValue;
         }
       }
       else if(i > 0){
         if(A[i] >= A[i+1]){
           dynamicMaximumValue = A[i];
           staticMaximumValue = Math.max(staticMaximumValue, dynamicMaximumValue);
         }
         else{
           dynamicMaximumValue = A[i+1];
           staticMaximumValue = Math.max(staticMaximumValue, dynamicMaximumValue);
         }
       }
     }
     double LastMaximumValue = staticMaximumValue;
     return LastMaximumValue;
   }

  //a utility method for finding the minimum value from an array of values
  double minimumOfArrayValues(double [] A){
     double staticMinimumValue = 0;
     double dynamicMinimumValue = 0;
     for(int i = 0; i < (A.length - 1); i++){
       if(i == 0){
         if(A[i] <= A[i+1]){
           dynamicMinimumValue = A[i];
           staticMinimumValue = dynamicMinimumValue;
         }
         else{
           dynamicMinimumValue = A[i+1];
           staticMinimumValue = dynamicMinimumValue;
         }
       }
       else if(i > 0){
         if(A[i] <= A[i+1]){
           dynamicMinimumValue = A[i];
           staticMinimumValue = Math.min(staticMinimumValue, dynamicMinimumValue);
         }
         else{
           dynamicMinimumValue = A[i+1];
           staticMinimumValue = Math.min(staticMinimumValue, dynamicMinimumValue);
         }
       }
     }
     double LastMinimumValue = staticMinimumValue;
     return LastMinimumValue;
   }
}
