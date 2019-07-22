/* @License Starts
 *
 * Copyright © 2015 - present. MongoExpUser
 *
 * License: MIT - See: https://github.com/MongoExpUser/Heavy-Oil-PVT-Simulator/blob/master/LICENSE
 *
 * @License Ends
 *
 *
 * ...Ecotert's PREos.java (released as open-source under MIT License) implements: Heavy oil (non-isothermal) PVT simulator.
 *
 *
 *  The PVT simulator uses 2 other implemented classes (ThermodynamicProperties.java and Utility.java)
 *
 *  The implementation is based on the solution of the "CLASSICAL" Peng-Robison EOS for a 3-phase system.
 *
 *  Reference:
 *  Peng, D.-Y and Robinson, D.B. (1976). Two and three phase equilibrium calculations for systems containing water.
 *  Canadian Journal of Chem. Eng. Vol. 54, pg. 595-599.
 *
*
 * The PVT simulator calculates relevant non-isothermal PVT properties that can be used:
 * a) To describe or characterize a given heavy oil/bitumen sample.
 * b) As input into a thermal numerical reservoir simulator.
 *
 *
 */


package eosPVT;

import eosPVT.Utility;
import eosPVT.ThermodynamicProperties;


public class PREos
{
  //contructor
  public void PREos(){}
  

  public double componentZFactor()
  {
      return 0;
   }
  
  public double componentKValue()
  {
      return 0;
  }
  
  public double [] compute(int n, int [] typeOfHC, double [] MassC1_C10_CO2_H2S_N2, double T_in_K,
                            double P_in_MPa, boolean moleFraction,double [] moleFractionC1_C10_CO2_H2S_N2)
  {
                              
   ThermodynamicProperties tp = new ThermodynamicProperties();
   Utility ut = new Utility();
  
   double R = 8.31447; ///Gas universal constant in SI unit (N/m)
   double P = P_in_MPa;
   double T = T_in_K;
   double fugacityError = 0;
   double [] bi = new double[typeOfHC.length];
   double [] PpcArray = new double[typeOfHC.length];
   double [] TpcArray = new double[typeOfHC.length];
   double [] Zi = new double[typeOfHC.length];  //Mole fraction of component, i in the mixture
   double [] ki = new double[typeOfHC.length];  //ki-factor for alphai
   double [] alpha = new double[typeOfHC.length];
   double [] aci = new double[typeOfHC.length];
   double [] ai = new double[typeOfHC.length];
   double [] AcC1_C10_CO2_H2S_N2 = new double[typeOfHC.length];
   double [] ti = new double[typeOfHC.length];
   double [] w = new double[typeOfHC.length];
   double [] zFactors = new double [3];
   double sumbi = 0;
   double sumPpcArray = 0;
   double sumTpcArray = 0;
   double sumAlpha = 0;
   double sumaci = 0;
   double sumai = 0;
   double sumaij = 0;
   double x = 0;
   double y = 0;
   double b = 0;
   double a = 0;
   double fhiv = 0;   double fhiL = 0;
   double lnfhiv = 0; double lnfhiL = 0;
   double [] fugacityv = new double[typeOfHC.length];
   double [] fugacityL = new double[typeOfHC.length];
   double [] fugacityRatio = new double[typeOfHC.length];
   double [] vc = new double [typeOfHC.length];  //critical volume of component, i
   double [] Zc = new double [typeOfHC.length];  //critical z-factor of compenent, i
   double [] Kai = new double [typeOfHC.length]; //K-value for component, i in the first-liquid phase (oleic)
   double [] Kbi = new double [typeOfHC.length]; //K-value for component, i in the first-liquid phase (acqeous)
   double La = 0;
   double Lb = 0;
   double V = 0;
   double [] xia = new double[typeOfHC.length];  //mole fraction of component, i in the first-liquid phase
   double [] xib = new double[typeOfHC.length];  //mole fraction of component, i in the second- liquid phase
   double [] xi = new double[typeOfHC.length];   //mole fraction of component, i in the liquid phase
   double yi = 0;
   double Zii = La + Lb + V;
   double ip = 0; //interaction parameter
   double AA [] = new double[4]; //LHS
   double BB [] = new double[4]; //RHS
   double CC [] = new double[4]; //solution Matrix
   double fugacityCount = 0;
  
  
   if (moleFraction == false)
   {
     Zi = tp.globalMoleFraction(n, typeOfHC, MassC1_C10_CO2_H2S_N2);
   }
   else if (moleFraction == true)
   {
     Zi = ut.copy(moleFractionC1_C10_CO2_H2S_N2);
   }
  
    //copy acentric factors to w (a single character) for easy code writing
    w = ut.copy(tp.AcC1_C10_CO2_H2S_N2);
  
    //Pseudo-reduced properties
    for(int i = 0; i < typeOfHC.length; i++)
    {
      x = T / tp.TcC1_C10_CO2_H2S_N2[typeOfHC[i]]; //Pseudo-reduced temperature
      y = P / tp.PcC1_C10_CO2_H2S_N2[typeOfHC[i]]; //Pseudo-reduced pressure
    }
  
    do{
      
        fugacityCount = fugacityCount + 1;
  
        //Use initial mole fraction in phases here, and calculate initial K-factor and other contraints for input into EOS
        for(int i = 0; i < typeOfHC.length; i++){
          if(fugacityCount == 1)
          {
             Kai[i] = 1 / y * Math.exp(5.3727 * (1 + w[i]) * (1 - 1 / x)); //Intially guessed equilibrium values (K-values) - oleic
             Kbi[i] = 10E+6 * (y / x);                                     //Intially guessed equilibrium values (K-values) - aqueous
          }
          else
          {
             Kai[i] = Kai[i] + fugacityError;                               //Adjusted equilibrium values (K-values) - oleic
             Kbi[i] = Kai[i] + fugacityError;                               //Adjusted equilibrium values (K-values) - aqueous
          }
  
        //Solve constraints equations here
  
  
           xi[i] = xia[i] + xib[i];                                         //sum of liquid phase fractions (oleic+aqueous)
      }
  
  
      //average properties and coefficients computations
      for (int i = 0; i < typeOfHC.length; i++)
      {
         bi[i] = (0.07780 * R) / (tp.TcC1_C10_CO2_H2S_N2[typeOfHC[i]] / tp.PcC1_C10_CO2_H2S_N2[typeOfHC[i]]);
         sumbi = sumbi + bi[i];
         b = b + xia[i] * bi[i];
         PpcArray[i] = xia[i] * tp.PcC1_C10_CO2_H2S_N2[typeOfHC[i]];
         sumPpcArray = sumPpcArray + PpcArray[i];
         TpcArray[i] = xia[i] * tp.TcC1_C10_CO2_H2S_N2[typeOfHC[i]];
         sumTpcArray = sumTpcArray + TpcArray[i];
    
         //ki-factor: Source - CMG WinProp 2000 Manual
         if(typeOfHC[i] <= 13){
           //that is, C1_C10_CO2_H2S_N2_H2O
           ki[i] = 0.37464 + 1.54226 * w[i] - 0.26992 * w[i] * w[i];
         }
         else if(typeOfHC[i] >= 14)
         {
           //that is, heavy hydrocarbons
           ki[i] = 0.379642 + 1.48503 * w[i] - 0.164423 * w[i] * w[i] + 0.016666 * w[i] * w[i] * w[i];
         }
    
         alpha[i] = Math.pow( (1 + ki[i] * (1 - Math.pow(T / tp.TcC1_C10_CO2_H2S_N2[typeOfHC[i]], 0.5))),2);
         sumAlpha = sumAlpha + alpha[i];
         aci[i] = 0.45724 *  (Math.pow( (tp.TcC1_C10_CO2_H2S_N2[typeOfHC[i]]*R), 2) / tp.PcC1_C10_CO2_H2S_N2[typeOfHC[i]]);
         sumaci = sumaci + aci[i];
         ai[i] = aci[i] * alpha[i];
         Zc[i] = 0.2905 - 0.085 * tp.AcC1_C10_CO2_H2S_N2[typeOfHC[i]];
    
          if(i == 0){
            sumai = sumai + sumai;
            sumaij = sumaij + sumaij;
          }
          else if(i > 0){
           for(int j = 1; j <= i; j++){
              //determining interaction parameter here with these sets of loops
             if(typeOfHC[i] == 13){ //H20_HC = 0.5
               ip = 0.5;
             }
             else{
               ip = 1 - ( 2*Math.pow( (vc[i]*vc[j-1]), 1/6) / ( Math.pow(vc[i], 1/3) + Math.pow(vc[j-1], 1/3) ) );
             }
    
             sumaij = sumaij + Math.pow( (ai[i]*ai[j-1]), 0.5 ) * xia[i]*xia[j-1] * ip;
             sumai = sumai + Math.pow( (ai[i]*ai[j-1]), 0.5 ) * xia[j-1] * ip;
           }
    
       }
    }
  
     a = sumaij;
     b = b;
     double A = (a * P) / (R * R * T * T);
     double B = (b * P) / (R * T);
  
     //coefficents of cubic equations
     double a3 = 1;
     double a2 = (B - 1);
     double a1 = (A - 3 * B * B - 2 * B);
     double a0 = (A * B - B * B - B * B * B);
  
     //variables for intermediate coefficients and solution (roots) of equation
     //(Based on Cardano's formula for cubic equation)
     double b0 = 0;
     double y1 = 0;
     double solution1 = 0;
     double b1 = 0;
     double y2 = 0;
     double solution2 = 0;
     double d2 = 0;
     double y3 = 0;
     double solution3 = 0;
     double d = 0;
     double r = 0;
     double phi = 0;
  
     //intermediate coefficients computations
     b1 = (1 / 3) * (3 * a1 - a2 * a2);
     b0 = (1 / 27) * (2 * a2 * a2 * a2 - 9 * a1 * a2 + 27 * a0);
  
     //discriminant
     d2 = Math.pow(b1/3, 3)  + Math.pow(b0/2, 2);
     d = Math.pow(d2, 1 / 2);
  
     //roots computation (solution)
     if (d2 > 0) {
       //3 thesame real roots
       y1 = y2 = y3 = Math.pow(-(b0 / 2) + d, 1 / 3) + Math.pow( (-b0 / 2) - d, 1 / 3);
       solution1 = solution2 = solution3 = (y1 - a2 / 3);
     }
     else if (d2 == 0) {
       //2 thesame real roots
       y1 = 2 * Math.pow(-b0 / 2, 1/3);
       y2 = y3 = -y1 / 2;
       solution1 = (y1 - a2 / 3);
       solution2 = solution3 = (y2 - a2 / 3);
     }
     else if (d2 < 0)
     {
       //different 3 roots
       r = (b0 / Math.abs(b0)) * Math.pow(Math.abs(b1)/ 3, 1/2);
       phi = Math.acos(b0 / (2 * r * r * r));
       y1 = -2 * r * Math.cos( (phi / 3) + 0);
       y2 = -2 * r * Math.cos( (phi / 3) + ( (2 / 3) * Math.PI));
       y3 = -2 * r * Math.cos( (phi / 3) + ( (4 / 3) * Math.PI));
       solution1 = (y1 - a2 / 3);
       solution2 = (y2 - a2 / 3);
       solution3 = (y3 - a2 / 3); ;
     }
  
     //Put solution (roots) into an array.
     zFactors[0] = solution1; zFactors[1] = solution2; zFactors[2] = solution3;
  
     //z-factors for liquid and vapour phases
     double vZ = ut.vapourCompressibilty(zFactors);
     double lZ = ut.liquidCompressibilty(zFactors);
  
     //Test for fugacity condition here
     for (int i = 0; i < typeOfHC.length; i++)
     {
       fhiv = bi[i]/b * (vZ-1) - ut.nl(vZ-B) - A/(2.828*B)*(2*(sumai)/a-bi[i]/b)*(ut.nl((vZ+2.414*B)/(vZ-0.414*B)));
       fhiL = bi[i]/b * (lZ-1) - ut.nl(lZ-B) - A/(2.828*B)*(2*(sumai)/a-bi[i]/b)*(ut.nl((lZ+2.414*B)/(lZ-0.414*B)));
       lnfhiv = ut.nl(fhiv);
       lnfhiL = ut.nl(fhiL);
       fugacityv[i] = P * xia[i]*Math.exp(lnfhiv);
       fugacityL[i] = P * xia[i]*Math.exp(lnfhiL);
       fugacityRatio[i] = fugacityL[i]/fugacityv[i];
     }
  
     fugacityError = Math.abs(1 - Math.abs(ut.maximumOfArrayValues(fugacityRatio)));
  
    }
    while(fugacityError > 0.00001);
  
     return zFactors;
     //return (fugacityRatio);
   }
  
  
  //Computing phase components fractions
  double [] phaseComponentsFraction(double [] Zi, double [] Kai, double [] Kbi)
  {
     int j = Zi.length;
     double [] Xai = new double[j];
     double [] Xbi = new double[j];
     double [] Yi = new double[j];
     double Ao = 0;  double Bo = 0;  double Co = 0; double Do = 0;  double Eo = 0;
     double La = 0;  double Lb = 0;  double denomenator = 0;
     double [] results = new double[4];
  
     for(int i = 0; i < j; i++)
     {
       Ao = Ao + (Zi[i] - Kai[i]);
       Bo = Bo + (Kai[i]/Kbi[i] - Kai[i]);
       Co = Co + Kai[i];
       Do = Do + ((Kai[i]/Kbi[i])*Zi[i]);
       Eo = Eo + Kai[i]*Zi[i];
    }
  
     La = ((1 - Co)/2*Ao);
     Lb = ((Do - Co)/2*Bo);
  
     results[0] = La;
     results[1] = Lb;
     results[2] = ut.sumArray(Zi) - Co ;
     results[3] = Do - Co;
  
     for(int i = 0; i < j; i++)
     {
       denomenator = (Kai[i]/Kbi[i] - Kai[i]);
       Xai[i] = Ao + (1 - Kai[i]);
       Xbi[i] = Bo + (Kai[i]/Kbi[i] - Kai[i]);
       Yi[i] = Do + (Kai[i]/Kbi[i]);
    }

  return results;
  
 }

}

