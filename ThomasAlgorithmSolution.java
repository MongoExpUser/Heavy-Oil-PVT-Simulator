



public class ThomasAlgorithmSolution {

   //Thomas algorithm method
   static double [] ThomasAlgorithmSolution (double upperElements [], double diagonalElements [],
                                      double lowerElements [], double rhsElements []){

     //For an m x n tri-diagonal matrix, number of rows (m) = number of cols n = (matrix size)
     int matrixSize = diagonalElements.length;

     //Arrays to hold (primary) main computed coefficients, secondary computed coefficients and solutions, x
     double A [] = new double [matrixSize]; 		//Primary
     double B [] = new double [matrixSize]; 		//Primary
     double C [] = new double [matrixSize]; 		//Primary
     double D [] = new double [matrixSize]; 		//Primary
     double W [] = new double [matrixSize]; 		//Secondary
     double G [] = new double [matrixSize]; 		//Secondary
     double E [] = new double [matrixSize]; 		//Solution

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

}
