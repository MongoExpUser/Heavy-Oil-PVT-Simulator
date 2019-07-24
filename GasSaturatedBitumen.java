/* @License Starts
 *
 * Copyright © 2002 - present. MongoExpUser
 *
 * License: MIT - See: https://github.com/MongoExpUser/Heavy-Oil-PVT-Simulator/blob/master/LICENSE
 *
 * @License Ends
 *
 *
 * ...Ecotert's GasSaturatedBitumen.java (released as open-source under MIT License) implements: 
 *
 * Simple viscosity calculation for gas-saturated bitumen.
 * 
 *
 */



package eosPVT;

import eosPVT.Utility;
import eosPVT.ThermodynamicProperties;

public class GasSaturatedBitumen {

  public static void main(String[] args) {

    ThermodynamicProperties tp = new ThermodynamicProperties();
    Utility ut = new Utility();
   
    System.out.println("Program for Calculating Viscosity of Gas-Saturated Bitumen");

    double T_in_K = 0;
    double P_in_Mpa = 0;
    double gasMoleFraction = 0;
    double viscosity_in_mPas = 0;
    double r1Vis_in_mPas = 0;
    double r2Vis_in_mPas = 0;
    double r1W = 0;
    double r2W = 0;
    boolean moleFraction = true;
    double Vc = 0;
    double VcSum = 0;
    double Tc = 0;
    double TcSum = 0;
    double M = 0;
    double MSum = 0;
    double Pc = 0;;
    double PcSum = 0;

    int solventType = 6;
    int numberOfInjectedComponents = 1;
    int totalComponents = numberOfInjectedComponents + 3;
    int [] typeOfHC = new int[totalComponents];
    String oilType = "cold lake"; oilType = "athabasca";
    double [] MassC1_C10_CO2_H2S_N2 = new double[totalComponents];
    double [] moleFractionC1_C10_CO2_H2S_N2 = new double[totalComponents];
    double oilMassFraction = 0.5;
    double oilMoleFraction = 0;
    double Zi [] = new double[totalComponents];
    double wi [] = new double[totalComponents];
    double w = 0;
    double wSum = 0;

    if(oilType == "athabasca"){
        typeOfHC[totalComponents-4] = solventType;
        typeOfHC[totalComponents-3] = 14;
        typeOfHC[totalComponents-2] = 15;
        typeOfHC[totalComponents-1] = 16;
        MassC1_C10_CO2_H2S_N2[totalComponents-4] = 1 - oilMoleFraction;
        MassC1_C10_CO2_H2S_N2[totalComponents-3] = tp.oilMassWeightFrac[0]*oilMassFraction;
        MassC1_C10_CO2_H2S_N2[totalComponents-2] = tp.oilMassWeightFrac[1]*oilMassFraction;
        MassC1_C10_CO2_H2S_N2[totalComponents-1] = tp.oilMassWeightFrac[2]*oilMassFraction;
    }
    else if(oilType == "cold lake"){
        typeOfHC[totalComponents-4] = solventType;
        typeOfHC[totalComponents-3] = 17;
        typeOfHC[totalComponents-2] = 18;
        typeOfHC[totalComponents-1] = 19;
        MassC1_C10_CO2_H2S_N2[totalComponents-4] = 1 - oilMoleFraction;
        MassC1_C10_CO2_H2S_N2[totalComponents-3] = tp.oilMassWeightFrac[3]*oilMassFraction;
        MassC1_C10_CO2_H2S_N2[totalComponents-2] = tp.oilMassWeightFrac[4]*oilMassFraction;
        MassC1_C10_CO2_H2S_N2[totalComponents-1] = tp.oilMassWeightFrac[5]*oilMassFraction;
    }


    if (moleFraction == false) {
      Zi = tp.globalMoleFraction(4, typeOfHC, MassC1_C10_CO2_H2S_N2);
    }
    else if (moleFraction == true) {
      Zi = ut.copy(moleFractionC1_C10_CO2_H2S_N2);
    }

    //copy acentric factors to w (a single character) for easy code writing
     wi = ut.copy(tp.AcC1_C10_CO2_H2S_N2);


    for (int i = 0; i < typeOfHC.length; i++) {
      Vc = Zi[i]*1;
      VcSum = VcSum + Vc;
      Tc = Zi[i]*tp.TcC1_C10_CO2_H2S_N2[typeOfHC[i]];
      TcSum = TcSum + Tc;
      M = Zi[i]*MassC1_C10_CO2_H2S_N2[i];
      MSum = MSum + M;
      Pc = Zi[i]*tp.PcC1_C10_CO2_H2S_N2[typeOfHC[i]];
      PcSum = PcSum + Pc;
      w = wi[typeOfHC[i]];
      wSum = wSum + w;
    }

    double netta = Math.pow(VcSum, 2/3) * Math.pow(TcSum, -1/2)  * Math.pow(MSum, -1/2) ;
    double a = ut.nl(netta*r1Vis_in_mPas);
    double b = (wSum - wi[0])/(1-wi[0]);
    double c = ut.nl(netta*r2Vis_in_mPas) - ut.nl(netta*r1Vis_in_mPas);;
    double d = a + b*c;
    viscosity_in_mPas = Math.exp(d);


    //testing enthalpy
    System.out.println("Viscosity of gas-satuared bitumen is" + " " + d + " mPas");

  }
}
