//System Parameters
int TSmax = ...; 
int Tmax = ...; 
int Smax = ...;
int Dmax = ...;
int Nmax = ...;
int Kmax = ...;

int ElectrolyzersAvailable=...;
int DeviceIndex = ...;
int nbOfScenarios=...;
range TS= 1..TSmax; //all stages
range S = 1..Smax; //stages                                                                                                                                                                                                                                                                                                                                                                                      es after stage 1
range D = 1..Nmax; //representative days (in this case 4)
range T = 1..Tmax; //time (24h)
range N = 1..ElectrolyzersAvailable; //Can buy up to 20 electrolysers in total (might not use all of them, but they are available)
range K = 1..Kmax; //2 device types with different parameters for every type of device used in the model (PV, Wind, Hydrogen tanks, Battery, Electrolysers)
range E= 1..nbOfScenarios; //2 different scenarios that have 2 different investment costs for everything as well as different hydrogen loads
int CapMaxPV = ...; //max installation capacities
int CapMaxWind = ...;
int CapMaxBA = ...;
int CapMaxSHST = ...;
float PenaltyElectricity = ...; //Penalty cost of unserved electricity
float PenaltyHydrogen = ...; //Penalty cost of unserved hydrogen
float Ir = ...; //Interest rate
int NbOfDataPoints = TSmax*Tmax*Nmax;
float ruPVWT[S][E]=...; //Investment factors for each scenario
float ruPTH[S][E]=...;
float ruBASHST[S][E]=...;
float ruELOAD[S][E]=...;
float ruHLOAD[S][E]=...;
float hprice[S][E]=...;
float eprice[S][E]=...;
float hprice1=...;
float eprice1=...;
float q[E]=...; //Probability of each scenario occuring
int budgetHN=...;
int budgetStages=...;
int CapPV0=...;
int CapWT0=...;

//Device parameters
tuple Electrolyzer {
  	 key string id;
     int Lifetime;
     float Capacity;
     float MaxEPI;
     float MinEPI;
     float ElectricityInputStandbyState;
     float a1;
     float a2;
     float InvestmentCost;
     float MaintenanceCost;
     float OperationCost;
     float EnergyRetention;
}

Electrolyzer ElectrolyzerParams = ...;

tuple PV {
  	 key string id;
     int Lifetime;
     float InvestmentCost;
     float MaintenanceCost;
     float OperationCost;
     float Efficiency;
}

PV PVParams = ...;

tuple Wind {
  	 key string id;
     int Lifetime;
     float InvestmentCost;
     float MaintenanceCost;
     float OperationCost;
     float Efficiency;
     float ConversionFactor;
     float BladeLength;
}

Wind WindParams= ...;

tuple Battery {
  	 key string id;
  	 int Lifetime;
     float Capacity;
     float ChargeEfficiency;
     float PowerCapacityRatio;
     float InvestmentCost;
     float MaintenanceCost;
     float OperationCost;
}

Battery BatteryParams = ...;

tuple HydrogenTank {
  	 key string id;
     int Lifetime;
     float Capacity;
     float ChargeEfficiency;
     float PowerCapacityRatio;
     float InvestmentCost;
     float MaintenanceCost;
     float OperationCost;
}

HydrogenTank HydrogenTankParams = ...;



//Defining ranges
int Stage[1..NbOfDataPoints] = ...; 
int Day[1..NbOfDataPoints] = ...;
int Time[1..NbOfDataPoints] = ...;
float PvOutput[1..NbOfDataPoints] = ...;
float WindOutput[1..NbOfDataPoints] = ...;
float ElectricityLoad[1..NbOfDataPoints] = ...;
float HydrogenLoad[1..NbOfDataPoints] = ...;
float ElectricityPrice[1..NbOfDataPoints] = ...;

//Extracting and defining system data in model
float yh[TS][D][T]; //hydrogen load data
float yep[TS][D][T]; //electricity price data
float ypv[TS][D][T]; //solar output data (per module in KW)
float ywt[TS][D][T]; //wind output data (per turbine in KW)
float ye[TS][D][T]; //electricity load data



execute{
  
    var dataPointCounter = 1;  // Initialize a counter for data points
    for (var ts in TS) {
        for (var d in D) {
            for (var t in T) {
                ye[ts][d][t] = ElectricityLoad[dataPointCounter];
                ywt[ts][d][t] = WindOutput[dataPointCounter];
                ypv[ts][d][t] = PvOutput[dataPointCounter];
                yep[ts][d][t] = ElectricityPrice[dataPointCounter];
                yh[ts][d][t] = HydrogenLoad[dataPointCounter];
                
                dataPointCounter++;

                // Check if dataPointCounter exceeds a certain value (e.g., NbOfDataPoints)
                if (dataPointCounter > NbOfDataPoints) {
                    break;
                }
            }

            // Check if dataPointCounter exceeds a certain value (e.g., NbOfDataPoints)
            if (dataPointCounter > NbOfDataPoints) {
                break;
            }
        }

        // Check if dataPointCounter exceeds a certain value (e.g., NbOfDataPoints)
        if (dataPointCounter > NbOfDataPoints) {
            break;
        }
    }

}


// Decision Variables
//Stage 1 decision variables (Here and now)
dvar float+ Ppth1[D][T][N]; // Represents the power produced by the nth electrolyzer unit at stage s during day D, time period T, unit K, and energy type E.
dvar float+ mpth1[D][T][N]; //Represents the hydrogen power output of the nth PtH
dvar float+ mpthshst1[D][T][N];
dvar float+ mpthload1[D][T][N];
dvar float+ SOCba1[D][T][K];//Represents the state of charge of BA at stage s on the dth representative day at time t, continuous variable
dvar float+ SOCshst1[D][T][K];//Represents the state of charge of hydrogen storage tank at stage s
dvar float+ Pbain1[D][T][K];//Represents the electricity power absorbed of BA at stage s on the dth representative day at time t, continuous variable
dvar float+ Pbaout1[D][T][K];//Represents the electricity power released of BA at stage s on the dth representative day at time t, continuous variable
dvar float+ Pbaoutgrid1[D][T][K];// Represents the electricity power released of BA to the grid at stage s on the dth representative day at time t, continuous variable
dvar float+ Pbaouth1[D][T][K];// Represents the electricity power released of BA to hydrogen production at stage s on the dth representative day at time t, continuous variable
dvar float+ mshstin1[D][T][K];//Represents the hydrogen power absorbed of SHST at stage s on the dth representative day at time t, continuous variable
dvar float+ mshstout1[D][T][K];//Represents the hydrogen power released of SHST at stage s on the dth representative day at time t, continuous variable

dvar float+ CapPV1;//Represents the installed capacity of PV at stage s, continuous variable
dvar boolean WT1[K];
dvar boolean BA1[K];
dvar boolean HT1[K];
dvar float+ Ppv1[D][T];//Represents the electricity power output of PV
dvar float+ Ppvgrid1[D][T];// Represents the electricity power output of PV to the grid
dvar float+ Ppvbattery1[D][T];//Represents the electricity power output of PV to the battery
dvar float+ Ppvh1[D][T];//Represents the electricity power output of PV to hydrogen production
dvar float+ Pwt1[D][T];//Represents the electricity power output of WT
dvar float+ Pwtgrid1[D][T];// Represents the electricity power output of WT to the grid
dvar float+ Pwtbattery1[D][T];//Represents the electricity power output of WT to the battery
dvar float+ Pwth1[D][T];//Represents the electricity power output of WT to hydrogen production
dvar float+ le1[D][T];//Represents the non-served electricity load
dvar float+ lh1[D][T];//Represents the non-served hydrogen load
dvar boolean zpth1[N];//Represents the installation status of the nth Pth at stage s, binary variable

//Stage 2+ decision variables (Wait and see)
dvar float+ Ppth[S][D][T][N][E]; // Represents the power produced by the nth electrolyzer unit at stage s during day D, time period T, unit K, and energy type E.
dvar float+ mpth[S][D][T][N][E]; //Represents the hydrogen power output of the nth PtH
dvar float+ mpthshst[S][D][T][N][E];
dvar float+ mpthload[S][D][T][N][E];
dvar float+ SOCba[S][D][T][K][E];//Represents the state of charge of BA at stage s on the dth representative day at time t, continuous variable
dvar float+ SOCshst[S][D][T][K][E];//Represents the state of charge of hydrogen storage tank at stage s
dvar float+ Pbain[S][D][T][K][E];//Represents the electricity power absorbed of BA at stage s on the dth representative day at time t, continuous variable
dvar float+ Pbaout[S][D][T][K][E];//Represents the electricity power released of BA at stage s on the dth representative day at time t, continuous variable
dvar float+ Pbaoutgrid[S][D][T][K][E];
dvar float+ Pbaouth[S][D][T][K][E];
dvar float+ mshstin[S][D][T][K][E];//Represents the hydrogen power absorbed of SHST at stage s on the dth representative day at time t, continuous variable
dvar float+ mshstout[S][D][T][K][E];//Represents the hydrogen power released of SHST at stage s on the dth representative day at time t, continuous variable
dvar float+ CapPV[S][E];//Represents the installed capacity of PV at stage s, continuous variable
dvar boolean WT[S][K][E];
dvar boolean BA[S][K][E];
dvar boolean HT[S][K][E];
dvar boolean zpth[S][N][E];//Represents the installation status of the nth Pth at stage s, binary variable
dvar float+ Ppv[S][D][T][E];//Represents the electricity power output of PV
dvar float+ Ppvgrid[S][D][T][E];
dvar float+ Ppvbattery[S][D][T][E];//Represents the electricity power output of PV
dvar float+ Ppvh[S][D][T][E];//Represents the electricity power output of PV
dvar float+ Pwt[S][D][T][E];//Represents the electricity power output of WT
dvar float+ Pwtgrid[S][D][T][E];
dvar float+ Pwtbattery[S][D][T][E];//Represents the electricity power output of WT
dvar float+ Pwth[S][D][T][E];//Represents the electricity power output of WT
dvar float+ le[S][D][T][E]; //Represents the non-served electricity load
dvar float+ lh[S][D][T][E]; //Represents the non-served hydrogen load
dvar float+ addCapPV[S][E];//Represents the installed capacity of PV at stage s, continuous variable
dvar int+ addCapWT[S][E];//Represents the installed capacity of Wind Turbines at stage s, continuous variable
dvar int+ addCapBA[S][E];//Represents the installed capacity of Batteries at stage s, continuous variable
dvar int+ addCapSHST[S][E];//Represents the installed capacity of Hydrogen tanks at stage s, continuous variable
dvar int+ addPth[S][E];
dvar float+ addCapPV1;
dvar int+ addCapWT1;
constraint hes[S][D][T][E];
//Cost function definition for objective function

// Adjusted Total Investment Cost Calculation for 3 Stages
//Here and now Investment cost 
     
dexpr float CinvPthHN =sum(n in N)ElectrolyzerParams.InvestmentCost * zpth1[n];
dexpr float CinvWindHN = WindParams.InvestmentCost * addCapWT1;
dexpr float CinvPVHN = PVParams.InvestmentCost * addCapPV1;  
dexpr float CinvBAHN = sum(k in K)BatteryParams.InvestmentCost * BA1[k];      
dexpr float CinvSHSTHN =sum(k in K)HydrogenTankParams.InvestmentCost * HT1[k];
dexpr float CinvHN =CinvPthHN+CinvWindHN+CinvPVHN+CinvBAHN+CinvSHSTHN;


// Wait and see Investment cost
dexpr float CinvPthWS[s in S][e in E] = 
    (1+ruPTH[s][e])*q[e]*((ElectrolyzerParams.InvestmentCost * addPth[s][e])/((1+Ir)^s));
dexpr float CinvWindWS[s in S][e in E] = 
    (1+ruPVWT[s][e])*q[e]*((WindParams.InvestmentCost * addCapWT[s][e])/((1+Ir)^s));
dexpr float CinvPVWS[s in S][e in E] = 
    (1+ruPVWT[s][e])*q[e]*((PVParams.InvestmentCost * addCapPV[s][e])/((1+Ir)^s));
dexpr float CinvBAWS[s in S][e in E] =  
    (1+ruBASHST[s][e])*q[e]*((BatteryParams.InvestmentCost * addCapBA[s][e])/((1+Ir)^s));
dexpr float CinvSHSTWS[s in S][e in E] =
    (1+ruBASHST[s][e])*q[e]*((HydrogenTankParams.InvestmentCost * addCapSHST[s][e])/((1+Ir)^s));
dexpr float CinvWS[s in S][e in E] = 
    CinvPthWS[s][e] + CinvWindWS[s][e] + CinvPVWS[s][e] + CinvBAWS[s][e] + CinvSHSTWS[s][e];

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Total Operations and Maintenance Cost
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

dexpr float ComPTHHN= 365/20*(sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T, n in N)
				  (ElectrolyzerParams.OperationCost*mpth1[d][t][n]
				 + ElectrolyzerParams.MaintenanceCost*mpth1[d][t][n]))
				 + 365/10* (sum(d in D: d==5 || d==7 || d==12 || d==14, t in T, n in N) 
				  (ElectrolyzerParams.OperationCost*mpth1[d][t][n]
				 + ElectrolyzerParams.MaintenanceCost*mpth1[d][t][n]));
			 
dexpr float ComBAHN= 365/20*(sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T, k in K) 
				  (BatteryParams.OperationCost*(Pbain1[d][t][k]+Pbaout1[d][t][k])
				 + BatteryParams.MaintenanceCost*(Pbain1[d][t][k]+Pbaout1[d][t][k])))
				 + 365/10*(sum(d in D: d==5 || d==7 || d==12 || d==14, t in T, k in K) 
				  (BatteryParams.OperationCost*(Pbain1[d][t][k]+Pbaout1[d][t][k])
				 + BatteryParams.MaintenanceCost*(Pbain1[d][t][k]+Pbaout1[d][t][k])));

dexpr float ComSHSTHN= 365/20*(sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T, k in K) 
				  (HydrogenTankParams.OperationCost*(mshstin1[d][t][k]+mshstout1[d][t][k])
				 + HydrogenTankParams.MaintenanceCost*(mshstin1[d][t][k]+mshstout1[d][t][k])))
				 + 365/10*(sum(d in D: d==5 || d==7 || d==12 || d==14, t in T, k in K) 
				  (HydrogenTankParams.OperationCost*(mshstin1[d][t][k]+mshstout1[d][t][k])
				 + HydrogenTankParams.MaintenanceCost*(mshstin1[d][t][k]+mshstout1[d][t][k]))) ;
				 
dexpr float ComPVHN= 365/20*(sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T)
				  (PVParams.OperationCost*Ppv1[d][t]
				 + PVParams.MaintenanceCost*Ppv1[d][t]))
				 + 365/10*(sum(d in D: d==5 || d==7 || d==12 || d==14, t in T)
				  (PVParams.OperationCost*Ppv1[d][t]
				 + PVParams.MaintenanceCost*Ppv1[d][t]));

dexpr float ComWindHN= 365/20*(sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T) 
				  (WindParams.OperationCost*Pwt1[d][t]
				 + WindParams.MaintenanceCost*Pwt1[d][t]))
				 + 365/10*(sum(d in D: d==5 || d==7 || d==12 || d==14, t in T) 
				  (WindParams.OperationCost*Pwt1[d][t]
				 + WindParams.MaintenanceCost*Pwt1[d][t]));
				 
dexpr float ComHN=ComWindHN+ComPVHN+ComSHSTHN+ComBAHN +ComPTHHN;
				 
// Electrolyzer Operation and Maintenance Costs for Wait and See
dexpr float ComPTHWS[s in S][e in E] = 
    q[e]*365/20 * (sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T, n in N) (
        ElectrolyzerParams.OperationCost * mpth[s][d][t][n][e] +
        ElectrolyzerParams.MaintenanceCost * mpth[s][d][t][n][e]))
        +  q[e]*365/10 * (sum(d in D: d==5 || d==7 || d==12 || d==14, t in T, n in N) (
        ElectrolyzerParams.OperationCost * mpth[s][d][t][n][e] +
        ElectrolyzerParams.MaintenanceCost * mpth[s][d][t][n][e]));

// Battery Operation and Maintenance Costs for Wait and See
dexpr float ComBAWS[s in S][e in E] = 
    q[e]*365/20 *(sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T, k in K) (
        BatteryParams.OperationCost * (Pbain[s][d][t][k][e]+Pbaout[s][d][t][k][e]) +
        BatteryParams.MaintenanceCost *  (Pbain[s][d][t][k][e]+Pbaout[s][d][t][k][e])))
        +q[e]*365/10 *(sum(d in D: d==5 || d==7 || d==12 || d==14, t in T, k in K) (
        BatteryParams.OperationCost * (Pbain[s][d][t][k][e]+Pbaout[s][d][t][k][e]) +
        BatteryParams.MaintenanceCost *  (Pbain[s][d][t][k][e]+Pbaout[s][d][t][k][e])));

// Hydrogen Storage Tank Operation and Maintenance Costs for Wait and See
dexpr float ComSHSTWS[s in S][e in E] = 
    q[e]*365/20 * (sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T, k in K) (
        HydrogenTankParams.OperationCost * (mshstin[s][d][t][k][e] + mshstout[s][d][t][k][e]) +
        HydrogenTankParams.MaintenanceCost * (mshstin[s][d][t][k][e] + mshstout[s][d][t][k][e])))
        + q[e]*365/10 * (sum(d in D: d==5 || d==7 || d==12 || d==14, t in T, k in K) (
        HydrogenTankParams.OperationCost * (mshstin[s][d][t][k][e] + mshstout[s][d][t][k][e]) +
        HydrogenTankParams.MaintenanceCost * (mshstin[s][d][t][k][e] + mshstout[s][d][t][k][e])));

// Photovoltaic (PV) Operation and Maintenance Costs for Wait and See
dexpr float ComPVWS[s in S][e in E] = 
    q[e]*365/20 * (sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T) (
        PVParams.OperationCost * Ppv[s][d][t][e] +
        PVParams.MaintenanceCost * Ppv[s][d][t][e]))
    + q[e]*365/10 * (sum(d in D: d==5 || d==7 || d==12 || d==14, t in T) (
        PVParams.OperationCost * Ppv[s][d][t][e] +
        PVParams.MaintenanceCost * Ppv[s][d][t][e]));

// Wind Turbine Operation and Maintenance Costs for Wait and See
dexpr float ComWindWS[s in S][e in E] = 
    q[e]*365/20 * (sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T) (
        WindParams.OperationCost * Pwt[s][d][t][e] +
        WindParams.MaintenanceCost * Pwt[s][d][t][e]))
+ q[e]*365/10 * (sum(d in D: d==5 || d==7 || d==12 || d==14, t in T) (
        WindParams.OperationCost * Pwt[s][d][t][e] +
        WindParams.MaintenanceCost * Pwt[s][d][t][e]));

dexpr float ComWS[s in S][e in E]=ComWindWS[s][e]+ComPVWS[s][e]+ComSHSTWS[s][e]+ComBAWS[s][e]+ComPTHWS[s][e];

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Total Penalty cost of unsatisfied loads
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

dexpr float CpenHN=365/20*(sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T)(PenaltyElectricity*le1[d][t]+PenaltyHydrogen*lh1[d][t]))
+365/10*(sum(d in D:d==5 || d==7 || d==12 || d==14, t in T)(PenaltyElectricity*le1[d][t]+PenaltyHydrogen*lh1[d][t]));

dexpr float CpenWS[s in S][e in E] =
    q[e]*365/20 * (sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T) (
        PenaltyElectricity * le[s][d][t][e] + PenaltyHydrogen * lh[s][d][t][e]))
        + q[e]*365/10 * (sum(d in D: d==5 || d==7 || d==12 || d==14, t in T) (
        PenaltyElectricity * le[s][d][t][e] + PenaltyHydrogen * lh[s][d][t][e]))
        ;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Revenue
dexpr float RevenueHN= 365/20*(sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T, n in N)(mpthload1[d][t][n]*hprice1)
+sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T, k in K)(mshstout1[d][t][k]*hprice1+ Pbaoutgrid1[d][t][k]*eprice1)
+sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T)(Ppvgrid1[d][t]+Pwtgrid1[d][t])*eprice1)
+ 365/10*(sum(d in D: d==5 || d==7 || d==12 || d==14, t in T, n in N) (mpthload1[d][t][n]*hprice1)
+sum(d in D: d==5 || d==7 || d==12 || d==14, t in T, k in K)(mshstout1[d][t][k]*hprice1+ Pbaoutgrid1[d][t][k]*eprice1)
+sum(d in D: d==5 || d==7 || d==12 || d==14, t in T)(Ppvgrid1[d][t]+Pwtgrid1[d][t])*eprice1);


dexpr float RevenueWS[s in S][e in E]= q[e]*365/20*(
sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T, n in N) (mpthload[s][d][t][n][e]*hprice[s][e])
+sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T, k in K)(mshstout[s][d][t][k][e]*hprice[s][e]+Pbaoutgrid[s][d][t][k][e]*eprice[s][e])
+sum(d in D: d>1 && d<18 && d!=5 && d!=7 && d!=12 && d!=14, t in T)(Ppvgrid[s][d][t][e]+Pwtgrid[s][d][t][e])*eprice[s][e])
+q[e]*365/10*(sum(d in D: d==5 || d==7 || d==12 || d==14, t in T, n in N)(mpthload[s][d][t][n][e]*hprice[s][e])
+sum(d in D: d==5 || d==7 || d==12 || d==14, t in T, k in K)(mshstout[s][d][t][k][e]*hprice[s][e]+Pbaoutgrid[s][d][t][k][e]*eprice[s][e])
+sum(d in D: d==5 || d==7 || d==12 || d==14, t in T)(Pwtgrid[s][d][t][e]+Ppvgrid[s][d][t][e])*eprice[s][e]);

//Objective Function
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dexpr float ProfitHN=RevenueHN-CpenHN-CinvHN-ComHN;
dexpr float ProfitInEveryCaseOfWS[s in S][e in E]=RevenueWS[s][e]-CpenWS[s][e]-ComWS[s][e]-CinvWS[s][e];
dexpr float TotalObj = ProfitHN +sum(s in S, e in E) ProfitInEveryCaseOfWS[s][e];
maximize TotalObj;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Constraints
 subject to {
forall(d in D, t in T, n in N) {
  /*
   Ppth1[d][t][n] <=ElectrolyzerParams.MaxEPI*zpth1[n];
   Ppth1[d][t][n]>= ElectrolyzerParams.ElectricityInputStandbyState*zpth1[n] + ElectrolyzerParams.MinEPI*zpth1[n];
   */
   Ppth1[d][t][n] <=ElectrolyzerParams.MaxEPI*zpth1[n];
   mpth1[d][t][n] <= ElectrolyzerParams.a1*Ppth1[d][t][n];
   mpth1[d][t][n] == mpthload1[d][t][n] + mpthshst1[d][t][n];
}


//forall(n in N) zpth1[n]==1;

forall(d in D, t in T, k in K){
   mshstin1[d][t][k] <= HydrogenTankParams.PowerCapacityRatio * HydrogenTankParams.Capacity*HT1[k];
   mshstout1[d][t][k] <= SOCshst1[d][t][k];
   mshstout1[d][t][k] <= HydrogenTankParams.PowerCapacityRatio * HydrogenTankParams.Capacity*HT1[k];
   Pbain1[d][t][k] <= BatteryParams.PowerCapacityRatio * BatteryParams.Capacity*BA1[k];
   Pbaout1[d][t][k]<=SOCba1[d][t][k]; 
   Pbaout1[d][t][k] <= BatteryParams.PowerCapacityRatio * BatteryParams.Capacity*BA1[k];
   Pbaout1[d][t][k] == Pbaoutgrid1[d][t][k] + Pbaouth1[d][t][k];
   SOCshst1[d][t][k] <= HydrogenTankParams.Capacity*HT1[k];
   SOCba1[d][t][k] <= BatteryParams.Capacity*BA1[k];
}
  

forall(d in D, t in T) {
   Pwt1[d][t] <= ywt[1][d][t]^3 * sum(k in K)WT1[k] * WindParams.BladeLength^2 * WindParams.ConversionFactor * WindParams.Efficiency;
   Pwt1[d][t] == Pwtgrid1[d][t] + Pwth1[d][t] + Pwtbattery1[d][t];
   Ppvgrid1[d][t] +  Pwtgrid1[d][t] +  sum(k in K)Pbaoutgrid1[d][t][k] == ye[1][d][t] - le1[d][t];
   sum(n in N) mpthload1[d][t][n] +  sum(k in K)mshstout1[d][t][k] == yh[1][d][t] - lh1[d][t];
   Ppv1[d][t] <= ypv[1][d][t] * CapPV1 * PVParams.Efficiency;
   Ppv1[d][t] == Ppvgrid1[d][t] + Ppvh1[d][t] + Ppvbattery1[d][t];
   sum(k in K)Pbain1[d][t][k] == Ppvbattery1[d][t] + Pwtbattery1[d][t];
   sum(k in K)mshstin1[d][t][k] <= sum(n in N) mpthshst1[d][t][n];
   sum(n in N) Ppth1[d][t][n] <= ElectrolyzerParams.EnergyRetention*(Ppvh1[d][t] + Pwth1[d][t] + sum(k in K)Pbaouth1[d][t][k]); 
}
  

// Hydrogen Storage Constraints
forall(k in K){
  SOCshst1[1][1][k] == mshstin1[1][1][k] * HydrogenTankParams.ChargeEfficiency - mshstout1[1][1][k] / HydrogenTankParams.ChargeEfficiency;
  SOCba1[1][1][k] == Pbain1[1][1][k] * BatteryParams.ChargeEfficiency - Pbaout1[1][1][k] / BatteryParams.ChargeEfficiency;
  } 
  
forall(d in D: d>1, k in K){
  SOCshst1[d][1][k] == SOCshst1[d-1][24][k] + (mshstin1[d][1][k] * HydrogenTankParams.ChargeEfficiency - mshstout1[d][1][k] / HydrogenTankParams.ChargeEfficiency);
  SOCba1[d][1][k] == SOCba1[d-1][24][k]+ (Pbain1[d][1][k] * BatteryParams.ChargeEfficiency - Pbaout1[d][1][k] / BatteryParams.ChargeEfficiency);
  } 
 
forall(d in D, t in T:t > 1, k in K){
  SOCshst1[d][t][k] == SOCshst1[d][t-1][k] + (mshstin1[d][t][k] * HydrogenTankParams.ChargeEfficiency - mshstout1[d][t][k] / HydrogenTankParams.ChargeEfficiency); 
  SOCba1[d][t][k] == SOCba1[d][t-1][k] + (Pbain1[d][t][k] * BatteryParams.ChargeEfficiency - Pbaout1[d][t][k] / BatteryParams.ChargeEfficiency);
  }

forall(k in K: k > 1){
  HT1[k-1] >= HT1[k];
  WT1[k-1] >= WT1[k];
  BA1[k-1] >= BA1[k];
  }

forall(n in N: n > 1) zpth1[n-1] >= zpth1[n];


// Capacity Constraints
CapPV1 <= CapMaxPV;
CapPV1==CapPV0+addCapPV1;
sum(k in K)WT1[k]==CapWT0+addCapWT1; 

sum(n in N) zpth1[n] <= ElectrolyzersAvailable;
sum(k in K) WT1[k]<=CapMaxWind;
sum(k in K) HT1[k]<=CapMaxSHST;
sum(k in K) BA1[k]<=CapMaxBA;

PVParams.InvestmentCost*addCapPV1+WindParams.InvestmentCost*addCapWT1+BatteryParams.InvestmentCost*sum(k in K)BA1[k]+HydrogenTankParams.InvestmentCost*sum(k in K)HT1[k] + ElectrolyzerParams.InvestmentCost*sum(n in N)zpth1[n]<=budgetHN;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Stages greater than 1  

//Pth constraints

   forall(s in S, d in D, t in T, n in N, e in E){
     /*
     Ppth[s][d][t][n][e]<=ElectrolyzerParams.MaxEPI*zpth[s][n][e];
     Ppth[s][d][t][n][e]>=(ElectrolyzerParams.ElectricityInputStandbyState+ElectrolyzerParams.MinEPI)*zpth[s][n][e];
     */
     Ppth[s][d][t][n][e]<=ElectrolyzerParams.MaxEPI*zpth[s][n][e];
     mpth[s][d][t][n][e]<=ElectrolyzerParams.a1*Ppth[s][d][t][n][e];
     mpth[s][d][t][n][e]==mpthload[s][d][t][n][e]+mpthshst[s][d][t][n][e];
     
   }
 
   forall(s in S, e in E){
     sum(n in N)zpth[s][n][e]<=ElectrolyzersAvailable;
     sum(k in K) BA[s][k][e]<=CapMaxBA;
     sum(k in K)HT[s][k][e]<= CapMaxSHST;
     sum (k in K) WT[s][k][e]<=CapMaxWind;
     }
     
   forall(s in S, n in N: n>1, e in E)zpth[s][n-1][e]>=zpth[s][n][e];
   
   forall(s in S, k in K: k>1, e in E){
     WT[s][k-1][e]>=WT[s][k][e];
     HT[s][k-1][e]>=HT[s][k][e];
     BA[s][k-1][e]>=BA[s][k][e];
     }
     
   forall(s in S:s>1, n in N, e in E) zpth[s][n][e]>=zpth[s-1][n][e];
   forall(s in S:s>1, k in K, e in E){
     WT[s][k][e]>=WT[s-1][k][e];
     HT[s][k][e]>=HT[s-1][k][e];
     BA[s][k][e]>=BA[s-1][k][e];
     }
  
   //Battery Storage Constraint
     
   forall(s in S, k in K, e in E){
     SOCba[s][1][1][k][e]==Pbain[s][1][1][k][e]*BatteryParams.ChargeEfficiency-Pbaout[s][1][1][k][e]/BatteryParams.ChargeEfficiency;
     SOCshst[s][1][1][k][e]==mshstin[s][1][1][k][e]*HydrogenTankParams.ChargeEfficiency-mshstout[s][1][1][k][e]/HydrogenTankParams.ChargeEfficiency;
     }
     
   forall(s in S, d in D: d>1, k in K, e in E){
     SOCba[s][d][1][k][e] == SOCba[s][d-1][24][k][e]+Pbain[s][d][1][k][e]*BatteryParams.ChargeEfficiency-Pbaout[s][d][1][k][e]/BatteryParams.ChargeEfficiency;
     SOCshst[s][d][1][k][e]==SOCshst[s][d-1][24][k][e]+(mshstin[s][d][1][k][e]*HydrogenTankParams.ChargeEfficiency-mshstout[s][d][1][k][e]/HydrogenTankParams.ChargeEfficiency);
   } 
   
   forall(s in S, d in D, t in T:t>1, k in K, e in E){
     SOCba[s][d][t][k][e]==SOCba[s][d][t-1][k][e]+Pbain[s][d][t][k][e]*BatteryParams.ChargeEfficiency-Pbaout[s][d][t][k][e]/BatteryParams.ChargeEfficiency;
     SOCshst[s][d][t][k][e]==SOCshst[s][d][t-1][k][e]+(mshstin[s][d][t][k][e]*HydrogenTankParams.ChargeEfficiency-mshstout[s][d][t][k][e]/HydrogenTankParams.ChargeEfficiency);
   }
   
   forall(s in S, d in D, t in T, k in K, e in E){
     Pbain[s][d][t][k][e]<=BatteryParams.PowerCapacityRatio*BA[s][k][e]*BatteryParams.Capacity;
     Pbaout[s][d][t][k][e]<=BatteryParams.PowerCapacityRatio*BA[s][k][e]*BatteryParams.Capacity;
     Pbaout[s][d][t][k][e]<=SOCba[s][d][t][k][e];
     SOCba[s][d][t][k][e]<=BatteryParams.Capacity;
     Pbaout[s][d][t][k][e]==Pbaoutgrid[s][d][t][k][e]+Pbaouth[s][d][t][k][e];
     
     mshstout[s][d][t][k][e]<=SOCshst[s][d][t][k][e];
     mshstin[s][d][t][k][e]<=HydrogenTankParams.PowerCapacityRatio*HT[s][k][e]*HydrogenTankParams.Capacity;
     mshstout[s][d][t][k][e]<=HydrogenTankParams.PowerCapacityRatio*HT[s][k][e]*HydrogenTankParams.Capacity;
     SOCshst[s][d][t][k][e]<=HydrogenTankParams.Capacity;
     }
     
     
     
   forall(s in S: s<Smax, e in E){
     sum(k in K) BA[s+1][k][e]==sum(k in K)BA[s][k][e]+addCapBA[s+1][e];
     sum(k in K) HT[s+1][k][e]==sum(k in K)HT[s][k][e]+addCapSHST[s+1][e];
     CapPV[s+1][e]==CapPV[s][e]+addCapPV[s+1][e];
     sum(k in K)WT[s+1][k][e]==sum(k in K)WT[s][k][e]+addCapWT[s+1][e]; 
     sum(n in N)zpth[s+1][n][e]==sum(n in N)zpth[s][n][e]+addPth[s+1][e];
   } 
 
//HES Constraint
    forall(s in S, d in D, t in T, e in E){
      Ppvgrid[s][d][t][e]+Pwtgrid[s][d][t][e] +sum(k in K)Pbaoutgrid[s][d][t][k][e]+le[s][d][t][e] == (1+ruELOAD[s][e])*ye[s+1][d][t];
      hes[s][d][t][e]: sum(n in N) mpthload[s][d][t][n][e]+sum(k in K)mshstout[s][d][t][k][e]+lh[s][d][t][e]== (1+ruHLOAD[s][e])*yh[s+1][d][t];
      Pwt[s][d][t][e]==Pwtgrid[s][d][t][e]+Pwth[s][d][t][e]+Pwtbattery[s][d][t][e];
      sum(k in K)Pbain[s][d][t][k][e]==Ppvbattery[s][d][t][e]+Pwtbattery[s][d][t][e];
      Ppv[s][d][t][e]==Ppvgrid[s][d][t][e]+Ppvh[s][d][t][e]+Ppvbattery[s][d][t][e];
      sum(k in K)mshstin[s][d][t][k][e]==sum(n in N)mpthshst[s][d][t][n][e];
      sum(n in N)Ppth[s][d][t][n][e]==ElectrolyzerParams.EnergyRetention*(Ppvh[s][d][t][e]+Pwth[s][d][t][e]+sum(k in K)Pbaouth[s][d][t][k][e]);
      }
        

    //PV and WT constraints
    forall(ts in TS:ts>1, s in S, d in D, t in T, e in E) {
    Ppv[s][d][t][e] <= ypv[ts][d][t] * CapPV[s][e]*PVParams.Efficiency;
    Pwt[s][d][t][e] <=ywt[ts][d][t]^3 * sum(k in K)WT[s][k][e] * WindParams.BladeLength^2 * WindParams.ConversionFactor * WindParams.Efficiency;
    }
  	
   
  
    //Constraints in relation to stage 1
    forall(e in E){
      sum(s in S)CapPV[s][e]<=CapMaxPV; 
      CapPV[1][e]==CapPV1+addCapPV[1][e];
      sum(k in K)WT[1][k][e]==sum(k in K)WT1[k]+addCapWT[1][e]; 
      sum(k in K)BA[1][k][e]==sum(k in K)BA1[k]+addCapBA[1][e];
      sum(k in K)HT[1][k][e]==sum(k in K)HT1[k]+addCapSHST[1][e];
      sum(n in N)zpth[1][n][e]==sum(n in N)zpth1[n]+addPth[1][e];
      sum(s in S)((1+ruPVWT[s][e])*PVParams.InvestmentCost*addCapPV[s][e]+(1+ruPVWT[s][e])*WindParams.InvestmentCost*addCapWT[s][e]
      +(1+ruBASHST[s][e])*BatteryParams.InvestmentCost*addCapBA[s][e]+(1+ruBASHST[s][e])*HydrogenTankParams.InvestmentCost*addCapSHST[s][e])
      +sum(s in S)((1+ruPTH[s][e])*ElectrolyzerParams.InvestmentCost*addPth[s][e])<=budgetStages;
      
   }

      
}

//Additional operational constraints
execute {
    var output = new IloOplOutputFile("ObjectiveValues.csv");
    var OptimalSolutionOfScenarioE;
    output.writeln("Scenario,ObjectiveValue");
    for (var e in E) { // Loop over days
    OptimalSolutionOfScenarioE = ProfitHN;
            for (var s in S) { // Loop over scenarios
               OptimalSolutionOfScenarioE = OptimalSolutionOfScenarioE + (ProfitInEveryCaseOfWS[s][e])/q[e]
            }
            output.writeln(e, ",", OptimalSolutionOfScenarioE);
            writeln("Total profit of scenario ", e, ", is ",OptimalSolutionOfScenarioE);
        }
    output.close();    
}


execute {
    // Function to write data for a variable to a CSV file
    function writeDataToCSV(fileName, variable, header) {
        var output = new IloOplOutputFile(fileName);
        output.writeln(header);
        for (var d in D)
        	for (var t in T)
            	for (var n in N)
            		output.writeln(d, ",", t, ",", n, "," ,variable[d][t][n]);
        output.close();
    }

    // Write data for the first set of variables
    writeDataToCSV("mpth1_file.csv", mpth1, "Day,Electrolyzer,Time,Supplier,Value");
    writeDataToCSV("mpthshst1_file.csv", mpthshst1, "Day,Electrolyzer,Time,Supplier,Value");
    writeDataToCSV("mpthload1_file.csv", mpthload1, "Day,Electrolyzer,Time,Supplier,Value");
    
}
execute {
    // Function to write data for a variable to a CSV file
    function writeDataToCSV(fileName, variable, header) {
        var output = new IloOplOutputFile(fileName);
        output.writeln(header);
        for (var d in D)
        	for (var t in T)
                for (var k in K)
                    output.writeln(d, ",", t, ",", k, ",", variable[d][t][k]);
        output.close();
    }

    // Write data for the second set of variables
    writeDataToCSV("SOCba1_file.csv", SOCba1, "Day,Time,Supplier,Value");
    writeDataToCSV("SOCshst1_file.csv", SOCshst1, "Day,Time,Supplier,Value");
    writeDataToCSV("Pbaout1_file.csv", Pbaout1, "Day,Time,Supplier,Value");
    writeDataToCSV("Pbaoutgrid1_file.csv", Pbaoutgrid1, "Day,Time,Supplier,Value");
    writeDataToCSV("Pbaouth1_file.csv", Pbaouth1, "Day,Time,Supplier,Value");
    writeDataToCSV("mshstout1_file.csv", mshstout1, "Day,Time,Supplier,Value");
}

execute {
    // Function to write data for a variable to a CSV file
    function writeDataToCSV(fileName, variable, header) {
        var output = new IloOplOutputFile(fileName);
        output.writeln(header);
        for (var d in D)
        	for (var t in T)
        		output.writeln(d, ",", t, ",", variable[d][t]);
        output.close();
    }

    // Write data for the second set of variables
   	writeDataToCSV("Ppv1_file.csv", Ppv1, "Day,Time,Value");
    writeDataToCSV("Ppvgrid1_file.csv", Ppvgrid1, "Day,Time,Value");
    writeDataToCSV("Ppvbattery1_file.csv", Ppvbattery1, "Day,Time,Value");
    writeDataToCSV("Ppvh1_file.csv", Ppvh1, "Day,Time,Value");
    writeDataToCSV("Pwt1_file.csv", Pwt1, "Day,Time,Value");
    writeDataToCSV("Pwtgrid1_file.csv", Pwtgrid1, "Day,Time,Value");
    writeDataToCSV("Pwtbattery1_file.csv", Pwtbattery1, "Day,Time,Value");
    writeDataToCSV("Pwth1_file.csv", Ppvh1, "Day,Time,Value");
    
    
}

execute {
    // Function to write data for a variable to a CSV file
    function writeDataToCSV(fileName, variable, header) {
        var output = new IloOplOutputFile(fileName);
        output.writeln(header);
        for (var s in S)
        	for (var d in D)
               	for (var t in T)
               		for (var n in N)
                    	for (var e in E)
                        	output.writeln(s, ",",d, ",", t, ",", n, ",", e, ",", variable[s][d][t][n][e]);
        output.close();
    }

    // Write data for the first set of variables
    writeDataToCSV("mpth_file.csv", mpth, "Stage,Day,Electrolyzer,Time,Supplier,Scenario,Value");
    writeDataToCSV("mpthshst_file.csv", mpthshst, "Stage,Day,Electrolyzer,Time,Supplier,Scenario,Value");
    writeDataToCSV("mpthload_file.csv", mpthload, "Stage,Day,Electrolyzer,Time,Supplier,Scenario,Value");
    
}

execute {
    // Function to write data for a variable to a CSV file
    function writeDataToCSV(fileName, variable, header) {
        var output = new IloOplOutputFile(fileName);
        output.writeln(header);
        for (var s in S)
        	for (var d in D)
            	for (var t in T)
                	for (var k in K)
                		for (var e in E)
                    		output.writeln(s, ",",d, ",", t, ",", k, ",", e, ",", variable[s][d][t][k][e]);
        output.close();
    }

    // Write data for the second set of variables
    writeDataToCSV("SOCba_file.csv", SOCba, "Stage,Day,Time,Supplier,Scenario,Value");
    writeDataToCSV("SOCshst_file.csv", SOCshst, "Stage,Day,Time,Supplier,Scenario,Value");
    writeDataToCSV("mshstout_file.csv", mshstout, "Stage,Day,Time,Supplier,Scenario,Value");
    writeDataToCSV("Pbaout_file.csv", Pbaout, "Stage,Day,Time,Supplier,Scenario,Value");
    writeDataToCSV("Pbaoutgrid_file.csv", Pbaoutgrid, "Stage,Day,Time,Supplier,Scenario,Value");
    writeDataToCSV("Pbaouth_file.csv", Pbaouth, "Stage,Day,Time,Supplier,Scenario,Value");
   
       
}

execute {
    // Function to write data for a variable to a CSV file
    function writeDataToCSV(fileName, variable, header) {
        var output = new IloOplOutputFile(fileName);
        output.writeln(header);
        for (var s in S)
        	for (var d in D)
            	for (var t in T)
                	for (var e in E)
                    	output.writeln(s, ",",d, ",", t, ",", e, ",", variable[s][d][t][e]);
        output.close();
    }

    // Write data for the second set of variables
    writeDataToCSV("le_file.csv", le, "Stage,Day,Time,Scenario,Value");
    writeDataToCSV("lh_file.csv", lh, "Stage,Day,Time,Scenario,Value");
    writeDataToCSV("Ppv_file.csv", Ppv, "Stage,Day,Time,Supplier,Scenario,Value");
    writeDataToCSV("Ppvgrid_file.csv", Ppvgrid, "Stage,Day,Time,Scenario,Value");
    writeDataToCSV("Ppvbattery_file.csv", Ppvbattery, "Stage,Day,Time,Scenario,Value");
    writeDataToCSV("Ppvh_file.csv", Ppvh, "Stage,Day,Time,Scenario,Value");
    writeDataToCSV("Pwt_file.csv", Pwt, "Stage,Day,Time,Scenario,Value");
    writeDataToCSV("Pwtgrid_file.csv", Pwtgrid, "Stage,Day,Time,Scenario,Value");
    writeDataToCSV("Pwtbattery_file.csv", Pwtbattery, "Stage,Day,Time,Scenario,Value");
    writeDataToCSV("Pwth_file.csv", Pwth, "Stage,Day,Time,Scenario,Value");
       
}

execute {
    // Function to write data for a variable to a CSV file
    function writeDataToCSV(fileName, variable, header) {
        var output = new IloOplOutputFile(fileName);
        output.writeln(header);
        for (var s in S)
        	for (var k in K)
            	for (var e in E)               	
                    output.writeln(s, ",",k, ",", e, ",", variable[s][k][e]);
        output.close();
    }

    // Write data for the second set of variables
    writeDataToCSV("CapWT.csv", WT, "Stage,Supplier,Scenario,Value");
    writeDataToCSV("CapBA.csv", BA, "Stage,Supplier,Scenario,Value");
    writeDataToCSV("CapSHST.csv", HT, "Stage,Supplier,Scenario,Value");
    
}

execute {
    // Function to write data for a variable to a CSV file
    function writeDataToCSV(fileName, variable, header) {
        var output = new IloOplOutputFile(fileName);
        output.writeln(header);
        for (var s in S)
        	for (var e in E)               	
        		output.writeln(s, ",", e, ",", variable[s][e]);
        output.close();
    }

    // Write data for the second set of variables
    writeDataToCSV("CapPV.csv", CapPV, "Stage, Scenario,Value")
    writeDataToCSV("addCapPV.csv", addCapPV, "Stage,Scenario,Value");
    writeDataToCSV("addCapWT.csv", addCapWT, "Stage,Scenario,Value");
    writeDataToCSV("addCapBA.csv", addCapBA, "Stage,Scenario,Value");
    writeDataToCSV("addCapSHST.csv", addCapSHST, "Stage,Scenario,Value");
    writeDataToCSV("addPth.csv", addPth, "Stage,Scenario,Value");
    
}



