/* @License Starts
 *
 * Copyright © 2002 - present. MongoExpUser
 *
 * License: MIT - See: https://github.com/MongoExpUser/Heavy-Oil-PVT-Simulator/blob/master/LICENSE
 *
 * @License Ends
 *
 *
 * ...Ecotert's Test.java (released as open-source under MIT License) implements: 
 *
 *
 *  Testing of the following 2 classes: ThermodynamicProperties() and Utility()
 *
 *
 */


package eosPVT;

import eosPVT.Utility;
import eosPVT.ThermodynamicProperties;

public class Test
{
  
  //contructor
  public void Test(){}
  

  public static void main(String[] args)
  {
     //test 1: Results: MW-19.8676, vis-17.2 microPas, z-0.829, 2nd vis-130 micro Pas. density - 119.1 Kg/m3

      //Input data
     //Initial Parameters
      ThermodynamicProperties tp = new ThermodynamicProperties();
      Utility ut = new Utility();
      int numberOfInjectedComponents = 4;
      int totalComponents = numberOfInjectedComponents + 3;
      int [] typeOfHC = new int[totalComponents];
      String oilType = "cold lake"; oilType = "athabasca";
      double [] MassC1_C10_CO2_H2S_N2 = new double[totalComponents];
      double [] moleFractionC1_C10_CO2_H2S_N2 = new double[totalComponents];
      double oilMassFraction = 0;
      double oilMoleFraction = 0;

      if(oilType == "athabasca"){
        typeOfHC[totalComponents-3] = 14;
        typeOfHC[totalComponents-2] = 15;
        typeOfHC[totalComponents-1] = 16;
        MassC1_C10_CO2_H2S_N2[totalComponents-3] = tp.oilMassWeightFrac[0]*oilMassFraction;
        MassC1_C10_CO2_H2S_N2[totalComponents-2] = tp.oilMassWeightFrac[1]*oilMassFraction;
        MassC1_C10_CO2_H2S_N2[totalComponents-1] = tp.oilMassWeightFrac[2]*oilMassFraction;
      }
      else if(oilType == "cold lake"){
        typeOfHC[totalComponents-3] = 17;
        typeOfHC[totalComponents-2] = 18;
        typeOfHC[totalComponents-1] = 19;
        MassC1_C10_CO2_H2S_N2[totalComponents-3] = tp.oilMassWeightFrac[3]*oilMassFraction;
        MassC1_C10_CO2_H2S_N2[totalComponents-2] = tp.oilMassWeightFrac[4]*oilMassFraction;
        MassC1_C10_CO2_H2S_N2[totalComponents-1] = tp.oilMassWeightFrac[5]*oilMassFraction;
      }

      boolean moleFraction2 = true;
      int n2 = 7;
      int [] typeOfHC2 = {0, 1, 2, 3, 3, 4, 14};
      double [] MassC1_C10_CO2_H2S_N22 = {0.8080, 0.1255, 0.0380, 0.0110, 0.0065, 0.0030, 0.0080}; //Unit: grams
      double [] moleFractionC1_C10_CO2_H2S_N22 = {0.8080, 0.1255, 0.0380, 0.0110, 0.0065, 0.0030, 0.0080};

      double P = 14.5; double T = 351.15; char ch = 'v'; double tt = 330;
      double initialReservoirTemperature = 12; //deg C - also equals initial oil temperture = ref. temp
      double initiaReservoirPressure = 5000;   //Pa - also equals intial oil pressure = ref. pressure
      double initialOilDensity  = 990;         //Kg-cum
      double initialOilViscosity = 0;          //Pas
      double oilBoilingPoint = 513.85;         //deg C - normal boiling point at  atmospheric/standard conditions
      double oilMolecularWeight = 0.4465;      //Kg
      double oilCoefficientOfExpansion = 0.0001;   // 1/K or 1/degC
      double Swi = 0;
      double Soi = 0;
      double Kh = 0;
      double Kv = 0;
      double porosity = 0;                        //Also equals reference porosity
      double rockThermalConductivity = 0;         //W/m degC
      double overBurdenThermalConductivity = 0;   //W/m degC
      double underBurdenThermalConductivity = 0;  //W/m degC
      double rockHeatCapacity = 0;                //J/cu-m degC
      double overBurdenHeatCapacity = 0;          //J/cu-m degC
      double underBurdenHeatCapacity = 0;         //J/cu-m degC
      double interactionFactor = 0.5;             //remove later, just here as dummy
      double rockCompressibility = 0.001;         //1/pa

      /**test 2
      thermodynamicProperties tp = new thermodynamicProperties();
      boolean moleFraction2 = true;
      int n2 = 6;
      int [] typeOfHC2 = {0, 1, 2, 3, 3, 11};
      double [] MassC1_C10_CO2_H2S_N22 = {0.9432, 0.0390, 0.0117, 0.0080, 0.00130, 0.004};
      double [] moleFractionC1_C10_CO2_H2S_N22 = {0.9432, 0.0390, 0.0117, 0.0080, 0.00130, 0.0040};
     */

      //Mole-fraction
      System.out.println("mole fracs: ");
      for(int i = 0; i < n2; i++){
        System.out.println(tp.globalMoleFraction(n2, typeOfHC2, MassC1_C10_CO2_H2S_N22)[i]*100);
      }

      System.out.print(" ");

      //results
      //double z = tp.zFactor(n2, typeOfHC2, MassC1_C10_CO2_H2S_N22, T, P, moleFraction2, moleFractionC1_C10_CO2_H2S_N22);
      double MW = tp.gcMW(n2, typeOfHC2, tp.MWC1_C10_CO2_H2S_N2, moleFractionC1_C10_CO2_H2S_N22);
      double density = tp.gcDensity(MW*1000, 0.829, T, P*1000000);
      double viscosity = tp.gcViscosity(MW*1000, density, T);
      double aa = tp.steamViscosity(tt, ch);
      double bb = tp.steamThermalConductivity(tt, ch);
      double cc = tp.steamDensity(tt, ch);
      double dd = tp.steamSpecificEnthalpy(tt, ch);
      double ff = tp.steamSaturatedPressure(tt);
      double gg = tp.steamSaturatedTemperature(ff);
      double hh = tp.oilViscosity(26.67, oilMolecularWeight);
      double ii = tp.oilDensity(260,5200000);
      double jj = tp.oilSpecificEnthalpy(26.67);
      double ww = tp.hcEnthalpy('o', 736.5, 588.706, 544.44,  91.94, 0.306);
      double ww2 = tp.hcEnthalpy('v', 848.3, 533.15, 0,  281.28, 0);

      //z-factor
      //System.out.println("z-factor: " + z);

      //Sum of molecular weight
      System.out.println("avg. MW.: " + MW);

      //Density
      System.out.println("density: " + density);

      //gas viscosity
      System.out.println("visc: " + viscosity);

      //Steam viscosity computation (Correlation from Tortike  and Farouq-Ali (1989))
       System.out.println("steam visc: " + aa + " Pa.s");

      //Steam thermal conductivity computation (Correlation from Tortike  and Farouq-Ali (1989))
      System.out.println("steam thermal conduc: " + bb + " W/m.K" );

      //Steam density computation (Correlation from Tortike  and Farouq-Ali (1989))
      System.out.println("steam density: " + cc + " Kg/cu-m" );

      //Steam specific enthalpy computation (Correlation from Tortike  and Farouq-Ali (1989))
      System.out.println("steam Specific enthalpy: " + dd + " J/kg" );

      //Steam saturated pressure computation (Correlation from Tortike and Farouq-Ali (1989))
      System.out.println("steam saturated pressure : " + ff + " Pa" );

      //Steam saturated temperature computation (Correlation from Tortike and Farouq-Ali (1989))
      System.out.println("steam saturated temperature: " + gg + " K" );

      //Oil viscosity (Edmunds (2000))
      System.out.println("Oil viscosity: " + hh + " Pa.s" );

      //Oil density (Edmunds (2000))
      System.out.println("Oil density: " + ii + " kg/cu-m" );

      //Oil enthalpy (Edmunds (2000))
      System.out.println("Oil enthalpy: " + jj + " J/kg" );

      //Testing PREos - acentric factor
      PREos p = new PREos();

      double [] pp =  p.compute(n2, typeOfHC2, MassC1_C10_CO2_H2S_N22, T, P, moleFraction2, moleFractionC1_C10_CO2_H2S_N22, tp);

      for(int i = 0; i < pp.length; i++){
        System.out.println("Z-factor: " + pp[i]);
      }

      double [] kk = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

      for(int i = 0; i < kk.length; i++){
        System.out.println( "erf of " + kk[i] + " = " + (ut.erf(kk[i])) );
      }

      System.out.println("C ");
      System.out.println(ut.minimumOfArrayValues(kk));

      for(int i = 0; i < typeOfHC.length; i++){
        System.out.println( "C " + " " + i + " " + typeOfHC[i] );
      }

      //testing enthalpy
      System.out.println( "enthalphy - btu.ib -  " + ww );
      System.out.println( "enthalphy - btu.ib -  " + ww2 );

      //testing sum of array
      double [] aaa = {10,20,30,1,1};
      System.out.println( "sum of array -  " + ut.sumArray(aaa) );

      //testing phaseCcomponentFraction
      double [] a = {0.2, 0.3, 0.1, 0.4};
      double [] b = {0.12, 0.01, 0.45, 0.30};
      double [] c = {0.88, 0.99, 0.55, 0.70};
      PREos pr = new PREos();
      double  [] d = pr.phaseComponentsFraction(a, b, c);
      ut.printArrayToSystem(d);

  }

}
