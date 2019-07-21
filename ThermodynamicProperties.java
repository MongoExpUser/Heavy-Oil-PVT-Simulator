/* @License Starts
 *
 * Copyright © 2015 - present. MongoExpUser
 *
 * License: MIT - See: https://github.com/MongoExpUser/Heavy-Oil-PVT-Simulator/blob/master/LICENSE
 *
 * @License Ends
 *
 *
 * ...Ecotert's ThermodynamicProperties.java (released as open-source under MIT License) implements:
 *
 *
 * Relevant thermodynamic properties required as input into heavy oil (non-isothermal) PVT simulator e.g.
 *
 *
 * (1) steam thermal properties
 * (2) oil thermal properties
 * (3) gas thermal properties
 * (4) mode fractions
 * (4) etc.
 *
 */



package eosPVT;


/* not used here: comment out
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Random;
*/

public class ThermodynamicProperties {
  
  //contructor
  public void ThermodynamicProperties(){}
  

  //Create a new utility class for general use
  Utility ut = new Utility();

  //Global mole fraction computation
  double [] globalMoleFraction(int n, int [] typeOfHC, double [] MassC1_C10_CO2_H2S_N2){
    // n = number of hydrocarbon components; typeOfHC = 1,2...13 , e.g. (n-Hexane, n-Propane) = n = [6-1,3-1]
    // MassC1_C10_CO2_H2S_N2 = masses of n. unit: Kg

    double [] globalMoleFraction = new double[typeOfHC.length];
    double [] moleArray = new double [typeOfHC.length];
    double sumMoleArray = 0;

     for (int i = 0; i < typeOfHC.length; i++) {
       moleArray[i] = MassC1_C10_CO2_H2S_N2[i]/(MWC1_C10_CO2_H2S_N2[typeOfHC[i]])/1000;
       sumMoleArray = sumMoleArray + moleArray[i];
     }

     for (int i = 0; i < typeOfHC.length; i++) {
      globalMoleFraction[i] = moleArray[i]/sumMoleArray;
     }

    return globalMoleFraction;
  }

  /** Defining, criticical values and other thermodynamics properties of gas condensates mixture,
   *  impurities, water and oil-in-place
   *   C1:Methane-0,  C2:Ethane-1,    C3:Propane-2,  C4:n-Butane-3, C5:n-Pentane-4
   *   C6:n-Hexane-5, C7:n-Heptane-6, C8:n-Octane-7, C9:n-Nonane-8, C10:n-Decane-9
   *   CO2:carbon dioxide-10
   *   H2S:hydrogen sulphide-11
   *   N2:nitrogen-12
   *   H20: water-13
   *   Oil: heavy oil - 14 - 19
   */

  //Molecular weight
  double [] MWC1_C10_CO2_H2S_N2 = {16.04f, 30.07f, 44.10f, 58.12f, 72.15f, 86.18f, 100.20f, 114.23f,
                                         128.26f, 142.29f, 44.01f, 34.08f, 28.01f, 18.015f, 485f, 595f, 822f, 300f, 675f, 793f};
                                         //molecular weight - unit: g.  Source: 0-12: Chierici (1994), page 25,
                                         //13:From WinProp Version 2000 user Guide of CMG (Page B-175)
                                         //14-16 are for Athabasca Bitumen based on 3PC (3-pseudo components) scheme
                                         //of  Mehrotra, A.K. and Svrcek, W.Y. (1988). Characterization of Athabasca bitumen
                                         //for gas solubility calculations. JCPT. Nov-Dec 1988, Vol.27, No.6, pgs. 107-110.
                                         //17-19 are for Cold lake Bitumen based on 3PC (3-pseudo components) scheme
                                         //of  Mehrotra, A.K. and Svrcek, W.Y. (1988). Correlation and predicting of gas solubility
                                         //in cold lake bitumen. The Can J. of Chem. Eng. Aug. 1988, Vol.66, pgs. 666-669.

  //Critical pressure
  double [] PcC1_C10_CO2_H2S_N2 = {4.606f, 4.881f, 4.247f, 3.799f, 3.372f, 3.013f, 2.737f, 2.489f, 2.289f,
                                         2.096f, 7.384f, 9.004f, 3.399f, 22.09f, 1.024f, 0.797f, 0.359f, 1.7003f, 0.5873f, 0.403f};
                                        //critical pressure - unit: Mpa. Source: 0-12: Chierici (1994), page 25
                                        //13: From WinProp Version 2000 user Guide of CMG (Page B-175).
                                        //14-16 and 17-19: Refs. as above

  //Critical temperature
  double [] TcC1_C10_CO2_H2S_N2 = {190.6f, 305.6f, 370.0f, 425.0f, 469.4f, 507.2f, 540.0f, 568.9f, 594.4f,
                                         617.8f, 304.4f, 373.3f, 126.1f, 647.15f, 991.9f, 1101.6f, 1304.3, 843.25f, 1105.85f, 1284.35f};
                                         //critical temperature - unit: K. Source: 0-12: Chierici (1994), page 25
                                         //13:(From WinProp Version 2000 user Guide of CMG (Page B-175)
                                         //14-16 and 17-19: Refs. as above
  //Boiling point
  double [] Boil_C1_C10_CO2_H2S_N2_degC = {-161.5f, -88.602f, -42.10f, -0.48f, 36.042f, 68.744f, 98.429f, 125.675f,
                                            150.821f, 174.152f, -78.5f, -60.2f, -196f, 100f, 550.85f,671.85f, 931.15f, 375.6f, 705.9f, 902.3f};
                                            //boiling point at atmospheric condition  - unit: degC
                                            //Source: 0-9: Zwolinski and Wilhot (1971); 10,12: http://www.engineeringtoolbox.com/24_155.html
                                            //11: http://www.airliquide.com/en/business/products/gases/gasdata/index.asp?GasID=59
                                            //13: http://chemed.chem.purdue.edu/genchem/topicreview/bp/ch14/melting.html
                                            //14-16 and 17-19: Refs. as above

  //Accentric Factor
  double [] AcC1_C10_CO2_H2S_N2 = {0.008f, 0.098f, 0.152f, 0.193f, 0.251f, 0.296f, 0.351f, 0.394f, 0.444f, 0.490f,
                                   0.225f, 0.0970f, 0.0400f, 0.344, 1.232f, 1.409f, 1.792f, 0.78813f, 1.51613f, 1.73820f};
                                           //Accentric Factor - unit: none. Sources:
                                           //nC1-nC10 - CMG WinProp User's Guide Version 2000
                                           //H20 - From WinProp Version 2000 user Guide of CMG (Page B-175)
                                           //N2  - From WinProp Version 2000 user Guide of CMG (Page B-175)
                                           //H2S - From Keskinen and Aalto (2001) - ATPC2001 Paper
                                           //CO2 - From Valderrama, J.O. and Alfaro, M. (2000). Liquid volumes from Generalized Cubic Equation
                                           // of State: Take it with Care. Oil & Gas Science Tech.- Rev. IFP Vol.55, No5. Page 524-531. page. 526
                                           //14-16 and 17-19: Refs. as above

  //Oil specific gravity based on 3PC (Athabasca and Cold lake Bitumen) at standard conditions
  double [] oilSG = {1.0230f, 1.0890f, 1.1630f, 0.958f, 1.050f, 1.160f};
                    // OilSG - unit: none. Sources:
                    //1-2 are for Athabasca Bitumen based on 3PC (3-pseudo components) scheme of  Mehrotra, A.K. and Svrcek, W.Y. (1988). Characterization of Athabasca bitumen
                    //for gas solubility calculations. JCPT. Nov-Dec 1988, Vol.27, No.6, pgs. 107-110.
                    //3-5 are for Cold lake Bitumen based on 3PC (3-pseudo components) scheme
                    //of  Mehrotra, A.K. and Svrcek, W.Y. (1988). Correlation and predicting of gas solubility in cold lake bitumen. The Can J. of Chem. Eng. Aug. 1988, Vol.66, pgs. 666-669.

  //Oil weight/mass percenatge based on 3PC (Athabasca and Cold lake Bitumen)
  double [] oilMassWeightFrac = {0.396f, 0.412f, 0.192f, 0.4320f, 0.3985f, 0.1695f};
                                //oilMassWeightFrac - unit: none. Source: As in oilSG above.

  //Interaction Parameter
  double InteractionParameterC1_C10_CO2_H2S_N2(int [] typeOfHC, double T_in_K, double P_in_MPa){
    double [] ip = {0.091, 0.135, 0.07, 0.5469};  //Propane interaction with N2, CO2, H2S and H20
    double R = 8.31447; ///Gas universal constant in SI unit (N/m)
    double P = P_in_MPa;
    double T = T_in_K;
    double [] dij = new double[12];

    double [] vc = new double [typeOfHC.length];
    double [] Zc = new double [typeOfHC.length];

    for (int i = 0; i < typeOfHC.length; i++){
      Zc[i] = 0.2905 - 0.085 * AcC1_C10_CO2_H2S_N2[typeOfHC[i]];
      vc[i] = Zc[i] * R * (TcC1_C10_CO2_H2S_N2[typeOfHC[i]] / PcC1_C10_CO2_H2S_N2[typeOfHC[i]]);
    }

    //H20_HC = 0.5
    //H20_N2 = 0.275
    //H20_CO2 = 0.2
    //H20_H2S = 0.12
    //HC_N2-CO2 = 0.1
    //HC_H2S = 0.055
    //HC-HC = code above

    //dij[i] = 1 - ( 2*Math.pow( (vci*vcj), 1/6) / ( Math.pow(vci, 1/3) + Math.pow(vcj, 1/3) ) );

     double vci = 0;
     double vcj = 0;
     return 0;
  } //Interaction Parameter - unit: none. Source: Lee, B.I. and Kesler, M.G. (1975). A generalized thermodynamics correlation
                                                  //based on three-parameter corresponding states. AICHE Journal Vol. 21, No. 3
                                                  //pge.510-527.

  //Oil density
  double oilDensity(double Temperature_in_DegC, double refTemperature_in_DegC){
    double T = Temperature_in_DegC;
    double Tf = refTemperature_in_DegC;
    double value = 0;
    return value;  //Unit : kg/cu-m Correlation from Edmunds (2000) in JCPT, page 32.
                   //This correlation is for clean Athabasca oil sands, but also assumed for cold lake in the absence
                   //of correlation specific to Cold lake. Cold lake specifice correlation may be added later
  }

  //Oil viscosity
  double oilViscosity(double Temperature_in_DegC, double MW){
    double T = Temperature_in_DegC + 273.15;
    double b = -27.23 + 13.56*Math.log(MW) - 1.886*(Math.log(MW)*Math.log(MW));
    double value = Math.pow( 10, (100*Math.pow((0.01*T), b)) )  -  0.8;
    return value;  //Unit : Pa.s. Source: Mehrotra, A.K. (1992). A model for the viscosity of bitumen/bitumen fractions-diluent blends
                   //JCPT. November 1992. Vol. 31, No.9. Pgs.28-32.
  }

  //Oil specific enthalpy
  double oilSpecificEnthalpy(double Temperature_in_DegC){
    double T = Temperature_in_DegC;
    double value = (((-2.895E-6*T) + 2.610E-3)*T + 1.557)*T*1E+3;
    return value;  //Unit = J/kg. Correlation from Edmunds (2000) in JCPT, page 32
                   //This correlation is for clean Athabasca oil sands, but also assumed for cold lake in the absence
                   //of correlation specific to Cold lake. Cold lake specifice correlation may be added later
  }

  //Gas condensates mixture density computation
  public double gcDensity(double MW_in_g, double zfactor, double T_in_K, double P_in_Pa){
    double value = (P_in_Pa * MW_in_g/1000) / (zfactor * 8.31447 * T_in_K);
    return value; //Unit = kg/cubic-meter
  }

  //Gas condensates mixture viscosity computation
  double gcViscosity(double MW_in_g, double density_in_Kg_cum, double T_in_K){
    double T = T_in_K;
    double x = 3.5 + 547.8/T + 0.01 * MW_in_g;
    double y = 2.4 - 0.2 * x;
    double K = (12.61 + 0.027*MW_in_g)*Math.pow(T, 1.5)/(116.11 + 10.56*MW_in_g + T);
    double value = K * Math.exp(x * Math.pow( (density_in_Kg_cum/1000), y )) * 0.0001;
    return value/1000; //Unit = Pa.s. Correlation from Lee et al. reported in Chierici (1994), pages 25 & 26
  }

  //Gas condensates mixture specific enthalpy computation (Correlation from Zwolinski and Wilhoit (1971))
   double gcSpecificEnthalpy(double Temperature_in_DegC, String phase, int n){
    //n = component type (0 to 13)
    double value = 0;
    return 0; //Unit = J/kg.  Correlation is valid in the range of 0 to 371.85 degC
   }

  //Gas condensates mixture molecular weight computation
  double gcMW(int n, int [] typeOfHC, double [] MWC1_C10_CO2_H2S_N2,  double [] MoleFraction){
    double [] MWArray = new double [typeOfHC.length];
    double MWi = 0;

    for (int i = 0; i < typeOfHC.length; i++) {
       MWArray[i] = MoleFraction[i] * MWC1_C10_CO2_H2S_N2[typeOfHC[i]];
       MWi = MWi + MWArray[i];
     }

    double value =  MWi/1000;
    return value; //Unit = kg
  }

  //Steam density computation (Correlation from Tortike  and Farouq-Ali (1989))
  double steamDensity(double Temperature_in_DegC, char phase){
      //Note a = phase 3 - aqueous, v = phase 2 - vapour, o = phase 1 - oleic
      double value = 0;
      double T = Temperature_in_DegC + 273.15; //Temperature converted to Kelvin (K)

      if(phase == 'a'){
        value = 3786.31 - 37.2487*T + 0.196246*(T*T) - 5.047E-4*(T*T*T) + 6.29368E-7*(T*T*T*T) - 3.08480E-10*(T*T*T*T*T);
      }
      else if(phase == 'v'){
        value  = value = Math.exp(-93.7072 + 0.833941*T - 0.00320809*(T*T) + 6.57652E-6*(T*T*T) - 6.93747E-9*(T*T*T*T) + 2.97203E-12*(T*T*T*T*T));
      }
      return value; //Unit = kg/cubic-meter. Correlation is valid in the range of 0 to 371.85 degC
  }

  //Steam viscosity computation (Correlation from Tortike  and Farouq-Ali (1989))
  double steamViscosity(double Temperature_in_DegC, char phase){
   //Note a = phase 3 - aqueous, v = phase 2 - vapour, o = phase 1 - oleic
   double value = 0;
   double T = Temperature_in_DegC + 273.15; //Temperature converted to Kelvin (K)

   if(phase == 'a'){
     value = -0.0123274 + 27.1038/T - 23527.5/(T*T) + 1.01425E+7/(T*T*T) - 2.17342E+9/(T*T*T*T) + 1.86935E+11/(T*T*T*T*T);
   }
   else if(phase == 'v'){
     value  = -5.46807E-4 + 6.89490E-6*T - 3.39999E-8*(T*T) + 8.29842E-11*(T*T*T) - 9.97060E-14*(T*T*T*T) + 4.71914E-17*(T*T*T*T*T);
   }
   return value; //Unit = Pa.s. Correlation is valid in the range of 0 to 371.85 degC
  }

  //Steam specific enthalpy computation (Correlation from Tortike  and Farouq-Ali (1989))
  double steamSpecificEnthalpy(double Temperature_in_DegC, char phase){
   //Note a = phase 3 - aqueous, v = phase 2 - vapour, o = phase 1 - oleic
   double value = 0;
   double T = Temperature_in_DegC + 273.15; //Temperature converted to Kelvin (K)

   if(phase == 'a'){
      value = ( (23665.2 - 366.232*T + 2.26952*(T*T) - 0.00730365*(T*T*T) + 1.30241E-5*(T*T*T*T)
                 - 1.22103E-8*(T*T*T*T*T)  + 4.70878E-12*(T*T*T*T*T*T) ) * 1000 );
    }
    else if(phase == 'v'){
      value  = ( (-22026.9 + 365.317*T - 2.25837*(T*T) + 0.00737420*(T*T*T) - 1.33437E-5*(T*T*T*T)
                 + 1.26913E-8*(T*T*T*T*T)  - 4.96880E-12*(T*T*T*T*T*T) ) * 1000 );
    }
   return value; //Unit = J/kg.  Correlation is valid in the range of 0 to 371.85 degC
  }

  //Steam thermal conductivity computation (Correlation from Tortike  and Farouq-Ali (1989))
  double steamThermalConductivity(double Temperature_in_DegC, char phase){
    //Note a = phase 3 - aqueous, v = phase 2 - vapour, o = phase 1 - oleic
    double value = 0;
    double T = Temperature_in_DegC + 273.15; //Temperature converted to Kelvin (K)

    if(phase == 'a'){
      value = 3.51153 - 0.0443602*T + 2.41233E-4*(T*T) - 6.05099E-7*(T*T*T) + 7.22766E-10*(T*T*T*T) - 3.37136E-13*(T*T*T*T*T);
    }
    else if(phase == 'v'){
      value  = -2.35787 + 0.0297429*T - 1.46888E-4*(T*T) + 3.57767E-7*(T*T*T) - 4.29764E-10*(T*T*T*T) + 2.04511E-13*(T*T*T*T*T);
    }
    return value; //Unit = W/m.K or W/m.degC. Correlation is valid in the range of 0 to 371.85 degC
  }

 //Steam saturated pressure computation (Correlation from Tortike and Farouq-Ali (1989))
 double steamSaturatedPressure(double Temperature_in_DegC){
   double value = 0; double v = 0;
   double T = Temperature_in_DegC + 273.15; //Temperature converted to Kelvin (K)
   v = (-175.776 + 2.29272*T - 0.0113953*T*T + 2.62780E-5*T*T*T - 2.73726E-8*T*T*T*T + 1.13816E-11*T*T*T*T*T);
   value = v * v * 1000;
   return value; //Unit = Pascal.  Correlation is valid in the range of 6.85 to 374.15 degC (
 }

 //Steam saturated temperature computation (Correlation from Tortike and Farouq-Ali (1989))
 double steamSaturatedTemperature(double Pressure_in_Pa){
    double P = ut.nl(Pressure_in_Pa/1000); //Natural log of saturated pressure (KPa)
    double value = (280.034 + 14.0856*P + 1.38075*(P*P) - 0.101806*(P*P*P) + 0.019017*(P*P*P*P));
    //double value = (561.435 + 33.8866*P + 2.18893*P*P + 0.0808998*P*P*P + 0.0342030*P*P*P*P);
    return value; //Unit = K.  Correlation is valid in the range of 0.611kPa to 22.12 Mpa
  }

  //Steam quality specification
  double steamQuality(float a){
    return a;
  }

  //Phase density
  double phaseDensity(char c, double odensity_in_Kgcum [], double omoleFrac [], double Temperature_in_DegC,double Pressure_in_P, double zfactor, double MW_in_g){
    double value = 0;
    double T = Temperature_in_DegC;
    double P = Pressure_in_P;

    if(c == 'o'){
      for(int i = 0; i < omoleFrac.length; i++){
        value = value + omoleFrac[i] * odensity_in_Kgcum [i];
      }
    }
    else if(c == 'v'){
         value = gcDensity(MW_in_g, zfactor, T, P);
    }
    else if(c == 'a'){
          value =  steamDensity(T, 'a');
    }
    return value; //unit = Kg/cum
  }

  //Phase viscosity
  double phaseViscosity(char c, double oviscosity_in_Pas [], double odensity_in_Kgcum [], double omoleFrac [], double Temperature_in_DegC,
                      double Pressure_in_P, double zfactor, double MW_in_g, double density_in_Kg_p_Cum_at_25DegC){
    double value = 0;
    double T = Temperature_in_DegC + 273.15;
    double P = Pressure_in_P;
    double vdensity_in_Kg_cum =  0;
    double b = -40.95 + 58.74*density_in_Kg_p_Cum_at_25DegC - 21.67*(density_in_Kg_p_Cum_at_25DegC*density_in_Kg_p_Cum_at_25DegC);
    double sumMass = 0;

    if(c == 'o'){
      for(int i = 0; i < omoleFrac.length; i++){
        value = omoleFrac[i] * (value + Math.pow( 10, (100*Math.pow((0.01*T), b)) )  -  0.8);
                //Source: Adapted from Mehrotra, A.K. (1992). A model for the viscosity of bitumen/bitumen fractions-diluent blends. JCPT. November 1992. Vol. 31, No.9. Pgs.28-32.
      }
    }
    else if(c == 'v'){
       vdensity_in_Kg_cum =  phaseDensity(c, odensity_in_Kgcum, omoleFrac, T, P, zfactor, MW_in_g);
       value =  gcViscosity(zfactor, vdensity_in_Kg_cum, T);
    }
    else if(c == 'a'){
       value =  steamViscosity(T, 'a');
    }
    return value; //unit = Pa.s
 }

  //Phase specific enthalpy
  double phaseSpecificEnthalpy(char c, double enthalpy_in_JpKg [], double moleFrac [], double Temperature_in_DegC,
                      double Pressure_in_P, double zfactor, double MW_in_g){
    double value = 0;
    double T = Temperature_in_DegC;
    double P = Pressure_in_P;

    if(c == 'o'){
      for(int i = 0; i < moleFrac.length; i++){
        value = value + moleFrac[i] * enthalpy_in_JpKg[i];
      }
    }
    else if(c == 'v'){
      for(int i = 0; i < moleFrac.length; i++){
        value = value + moleFrac[i] * enthalpy_in_JpKg[i];
      }
    }
    else if(c == 'a'){
          value =  steamSpecificEnthalpy(T, 'a');
    }
    return value; //unit = J/kg
  }

  //Hydrocarbon enthalpy
  double hcEnthalpy(char c, double density_in_Kg_cum, double T_in_K, double Tc_in_K, double boilingPt_in_DegC, double AccFac){
    double SG = density_in_Kg_cum/1000;
    double Tbr = ((boilingPt_in_DegC + 273.15)*1.8) / (Tc_in_K*1.8); //Ratio in R
    double Tb_in_R = (boilingPt_in_DegC + 273.15)*1.8;  //boiling point in R
    double Tc_in_R = 1.8 * Tc_in_K;                    //Tc in R
    double T_in_F = (T_in_K - 273.15)*1.8 + 32;        //Temp. in F
    double T_in_R = (T_in_F + 459.67);                  //Temp. in R
    double K = ( Math.pow( (Tb_in_R), 1/3) / SG );     //Watson characterization factor
    double value = 0;
    double CF = Math.pow(( (12.8 - K) * (10 - K)/(10 * AccFac) ), 2);

    if(c == 'o'){
      value = ((0.35 + 0.055*K) * (0.6811 - 0.308*SG + (0.000815 + 0.000306*SG)*(T_in_F) )) * T_in_F;
    }
    else if(c == 'v'){
      if(Tbr <= 0.8){
          value = (( -0.32646 + 0.02678*K - (1.3892 - 1.2122*K + 0.03803*K*K) * 10E-4 * T_in_R  - 1.5393*10E-7 * T_in_R * T_in_R) )* T_in_R;
      }
      else if(Tbr > 0.8){
          value = (- 0.33886 + 0.02827 * K - (0.9291 - 1.1543*K + 0.0368*K*K)*10E-4 * T_in_R
                  - 1.6658 * 10E-7 * T_in_R * T_in_R - CF * (0.26105 - 0.59332 * AccFac - (4.56 - 9.48*AccFac)*10E-4 * T_in_R
                  - (0.536 - 0.6828 * AccFac * 10E-7) * T_in_R * T_in_R ) )*T_in_R;
     }
   }
    return (value);//*0.0004299226);     //converted to unit: J/kg; Correlations from Kesler, M.G. and Lee, B.I. (1976). Improve prediction of enthalpy of fractions. Hydrocarbon
                     //processing. Vol. 55, No.3, Pgs. 153-158, March 1976.
  }

  //Phase vapour pressure
  double phaseVapourPressure(){
    return 0;
  }

  //Phase liquid pressure
  double phaseLiquidPressure(){
    return 0;
  }

  //Oil and HC fraction specific heat capacity (liquid state)
  double OilHCSpecificHeatCapacity(double Temperature_in_degC, double SG){
    double T = Temperature_in_degC;
    return  ( (1.685 + 0.0034*T) / (Math.sqrt(SG)) )* 1000;  //Unit; J/Kg/K.
    //Source: Burger, J., Sourieau, P. and Combarnous, M. (1985). Thermal methods of oil recovery. Gulf Publishing Company, Houston, US. Pg.48.
  }

  //Rock specific heat capacity
  double rockSpecificHeatCapacity(double Temperature_in_degC){
    double T = Temperature_in_degC;
    return  (0.8 + 1.3E-3*T + 9E-7*T*T) * 1000;  //Unit; J/Kg/K.
    //Source: Burger, J., Sourieau, P. and Combarnous, M. (1985). Thermal methods of oil recovery. Gulf Publishing Company, Houston, US. Pg.48.
  }

  //Interfacial tension between two phases
  double interfacialTension(double pi, double xij){
    double  x = 0;
    double value = 0;
    return value; //Unit = Pa
  }

  boolean phaseAppearDisapper(double moleFraction){
    boolean state = true;

    //Use sum mole of fraction in a phase to estalish if a phase still exists
    if(moleFraction > 0 && moleFraction <=1){
      state = true;
    }
    else{
      state =false;
    }
    return state;
  }

}
