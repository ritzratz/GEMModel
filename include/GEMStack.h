/* ===============================================================*/
/* GEM stack calculations v1.0                                    */
/*                                                                */
/* Written by Viktor Ratza, ratza@hiskp.uni-bonn.de, May 2019     */
/* ===============================================================*/


/* -----------------
 * General Remarks
 * -----------------
 * 
 * All distances in units of cm.
 * All electric potentials in units of Volts.
 * All electric fields in units of V/cm.
 * 
 * GEMId and VolumeId start counting from 0, e.g. first GEM or first volume in stack has ID 0.
 * 
 * Electron efficiencies are calculated from the model calculations using the tuning parameters s1, s2 and s3.
 * 
 * Ion efficiencies are calculated from polynomial fits to the simualted data.
 * 
 * Absolute gain are calculated from fitting exponential functions to the simulated gain curves.
 * 
 * Fluctuatios f are calculated from fitting Polya distribution to the single electron amplifications (from the gain simulations).
 * 
 */

#include "GEMEfficiencies.h"

#define CHRG_ION_AVAL 2
#define CHRG_ION_DRIFT 1

#define EFF_COLL 0
#define EFF_EXTR 1

#define GAS_NECO2 4
#define GAS_NECO2N2 2

#define GEM_ST 0
#define GEM_MP 1
#define GEM_LP 2
#define GEM_NO 3 //No GEM (used for avalanche ion efficiencies)

#define MODEL_ARC 0 //Single gain fluctuation arctan model
#define MODEL_EXP 1 //Single gain fluctuation exp model
#define MODEL_EXP_NOLIMIT 2 //Same as MODEL_EXP but can reach values > 1

#define MODEL_DRIFT_ION_SIMULATION 0 //Load drift ion efficiencies from simulations
#define MODEL_DRIFT_ION_CUT 1     //Cut after plataeau and scaling for values higher 

#define GEM_Efficiencies_N 200 //Default 200. Lower values will lead to faster calculations but loss of precision.

int COUT_MESSAGES = 0;
int COUT_WARNINGS = 0;

class GEMStack
{
    public:
	GEMStack();
	~GEMStack();

	//-------------------------	
	//Stack Definition
	//-------------------------
	
	    //Add a new field volume e.g. transfer field, induction field. Return VolumeId number.
	    int AddVolume(double VolumeDistance, double VolumeField);
	    void SetVolumeField(int VolumeId, double VolumeField) { mVolumeField[VolumeId] = VolumeField; };
	    double GetVolumeField(int VolumeId) { return mVolumeField[VolumeId]; };
	    
	    //Add a new GEM to the stack. Return GEMId number.
	    int AddGEM(int GEMType, double GEMVoltage);
	    void SetGEMVoltage(int GEMId, double GEMVoltage) { mGEMVoltage[GEMId] = GEMVoltage; };
	    double GetGEMVoltage(int GEMId) { return mGEMVoltage[GEMId]; };
	    
	    //Set the GEM type of a given GEM
	    void SetGEMType(int GEMId, int GEMType) { mGEMType[GEMId] = GEMType; };
	    
	    //Set the gas type of the detector
	    void SetGas(int GasType);
	    
	    //Set the photon energy which converges in drift region (in eV)
	    void SetPhotonEnergy(double PhotonEnery) { mPhotonEnergy = PhotonEnery; };
	    
	    //Enable / disable attachment: 1 on, 0 off (default)
	    void SetAttachment(int AttachmentStatus) { mAttachmentStatus = AttachmentStatus; };
	    
	    //Get the attachment coefficient (in 1/cm) for a given electric field
	    double GetAttachmentFactor(double VolumeField); //TODO: Currently this is just a fixed value! There is no field dep. implemented yet.
	    double GetErrorAttachmentFactor(double VolumeField);
	    
	    //Define a (yet constant) attachment factor (1/cm)
	    void SetAttachmentFactor(double AttachmentFactor) { mAttachmentFactor = AttachmentFactor; };
	
	//-------------------------
	//Electron efficiencies
	//-------------------------

	    //Get electron collection / extraction efficiency (results from fit model COLLTOP 6.0, EXTRBOT 6.1)
	    double GetGEMElectronCollectionCore(int GEMId, int Error);
	    double GetGEMElectronCollection(int GEMId);
	    double GetErrorGEMElectronCollection(int GEMId);
	    
	    double GetGEMElectronExtractionCore(int GEMId, int Error);
	    double GetGEMElectronExtraction(int GEMId);
	    double GetErrorGEMElectronExtraction(int GEMId);
	    
	    
	    //Fix Electron collection / extraction efficiency for a given GEM in detector
	    void SetGEMElectronCollection(int GEMId, double ElectronCollectionEfficiency) { mGEMFixElectronCollection[GEMId] = ElectronCollectionEfficiency; mGEMFixElectronCollectionStatus[GEMId] = 1; };
	    void SetGEMElectronExtraction(int GEMId, double ElectronExtractionEfficiency) { mGEMFixElectronExtraction[GEMId] = ElectronExtractionEfficiency; mGEMFixElectronExtractionStatus[GEMId] = 1; };

	//-------------------------
	//Avalanche ion efficiencies
	//-------------------------
	    
	    //Get avalanche ion collection efficiency (results from polynomial fits) for a specific GEM
	    double GetGEMAvalancheIonCollection(int GEMId);
	    double GetErrorGEMAvalancheIonCollection(int GEMId);
	    
	    //Get avalanche ion extraction efficiency (results from polynomial fits) for a specific GEM
	    double GetGEMAvalancheIonExtraction(int GEMId);
	    double GetErrorGEMAvalancheIonExtraction(int GEMId);
	    
	    //Get drift ion collection efficiency (results from polynomial fits) for a specific GEM
	    double GetGEMDriftIonCollection(int GEMId);
	    double GetErrorGEMDriftIonCollection(int GEMId);
	    
	    //Get drift ion extraction efficiency (results from polynomial fits) for a specific GEM
	    double GetGEMDriftIonExtraction(int GEMId);
	    double GetErrorGEMDriftIonExtraction(int GEMId);
	    
	    //Fix avalanche ion collection for a specific GEM within the stack
	    void SetGEMAvalancheIonCollection(int GEMId, double AvalancheIonCollection) { mGEMFixAvalancheIonCollectionStatus[GEMId] = 1; mGEMFixAvalancheIonCollection[GEMId] = AvalancheIonCollection; };

	    //Fix avalanche ion extraction for a specific GEM within the stack
	    void SetGEMAvalancheIonExtraction(int GEMId, double AvalancheIonExtraction) { mGEMFixAvalancheIonExtractionStatus[GEMId] = 1; mGEMFixAvalancheIonExtraction[GEMId] = AvalancheIonExtraction; };    

	    //Fix drift ion collection for a specific GEM within the stack
	    void SetGEMDriftIonCollection(int GEMId, double DriftIonCollection) { mGEMFixDriftIonCollectionStatus[GEMId] = 1; mGEMFixDriftIonCollection[GEMId] = DriftIonCollection; };

	    //Fix drift ion extraction for a specific GEM within the stack
	    void SetGEMDriftIonExtraction(int GEMId, double DriftIonExtraction) { mGEMFixDriftIonExtractionStatus[GEMId] = 1; mGEMFixDriftIonExtraction[GEMId] = DriftIonExtraction; };    
	    
	    //Set general avalanche ion collection efficiency for a specific gas and GEM type (p0+p1*Eta2+...)
	    void SetAvalancheIonCollection(int GasType, int GEMType1, double p0, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double FitEta2Min, double FitEta2Max);
	    void SetErrorAvalancheIonCollection(int GasType, int GEMType1, double dp0, double dp1, double dp2, double dp3, double dp4, double dp5, double dp6, double dp7);
	    
	    //Set general avalanche ion extraction efficiency for a specific gas and GEM type (p0+p1*Eta2+...)
	    void SetAvalancheIonExtraction(int GasType, int GEMType1, double p0, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double FitEta1Min, double FitEta1Max);
	    void SetErrorAvalancheIonExtraction(int GasType, int GEMType1, double dp0, double dp1, double dp2, double dp3, double dp4, double dp5, double dp6, double dp7);
	    	    
	    //Set general drift ion collection efficiency for a specific gas and GEM type (p0+p1*Eta2+...)
	    void SetDriftIonCollection(int GasType, int GEMType1, int GEMType2, double p0, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double FitEta2Min, double FitEta2Max);
	    void SetErrorDriftIonCollection(int GasType, int GEMType1, int GEMType2, double dp0, double dp1, double dp2, double dp3, double dp4, double dp5, double dp6, double dp7);
	    
	    //Set general drift ion extraction efficiency for a specific gas and GEM type (p0+p1*Eta2+...)
	    void SetDriftIonExtraction(int GasType, int GEMType1, int GEMType2, double p0, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double FitEta1Min, double FitEta1Max);
	    void SetErrorDriftIonExtraction(int GasType, int GEMType1, int GEMType2, double dp0, double dp1, double dp2, double dp3, double dp4, double dp5, double dp6, double dp7);
	    
	    void SetDriftIonCollectionScalingFactor(double DriftIonCollectionScalingFactor) { mDriftIonCollectionScalingFactor = DriftIonCollectionScalingFactor; };
	    void SetDriftIonExtractionScalingFactor(double DriftIonExtractionScalingFactor) { mDriftIonExtractionScalingFactor = DriftIonExtractionScalingFactor; };
	    
	    void SetAvalancheIonEfficienciesScalingFactor(double AvalancheIonEfficienciesScalingFactor) { mAvalancheIonEfficienciesScalingFactor = AvalancheIonEfficienciesScalingFactor; };
	    
	    
	    void SetDriftIonModel(int DriftIonModel) { mDriftIonModel = DriftIonModel; };
	    
	//-------------------------
	//Gain calculations
	//-------------------------	
	
	    //Get the absolute gain for a GEM with GEMId
	    double GetGEMAbsGainCore(int GEMId, int Error);
	    double GetGEMAbsGain(int GEMId);
	    double GetErrorGEMAbsGain(int GEMId);
	    
	    //TODO: Also implement errors here
	    //Fix the gain curve for an individual GEM within the stack: exp(Constant+Slope*UGEM)
	    void SetGEMAbsGain(int GEMId, double Slope, double Constant);

	    //Overwrite / Set the exponential gadin curve for a GEM type (all GEMs) and gas type: exp(Constant+Slope*UGEM)
	    void SetAbsGain(int GasType, int GEMType, double Slope, double Constant, double ErrorSlope, double ErrorConstant, double FitUGEMMin, double FitUGEMMax);
	    
	    //Get the effective gain for a specific GEM within the stack
	    double GetGEMEffGain(int GEMId);
	    double GetErrorGEMEffGain(int GEMId);
	
	    //Allows a scaling of the all calculated gain curves (by default this is 1.0, i.e. no scaling)
	    void SetAbsGainScalingFactor(double GainScaleFactor) { mGainScaleFactor = GainScaleFactor; };
	
	    //Get the single gain fluctuation f from Polya distribution (1 exponential, < 1 Polya distribution for higher GEM voltages)
	    double GetGEMSingleGainFluctuation(int GEMId);
	    
	    double GetErrorHiGEMSingleGainFluctuation(int GEMId);
	    double GetErrorLoGEMSingleGainFluctuation(int GEMId);
	    double GetErrorMaxGEMSingleGainFluctuation(int GEMId);
	    
	    //Set the model which will be used for the single gain calculation
	    //0 MODEL_ARC: 0.5(f0+1)+1/pi(1-f0)arctan(-Q(U-U0))
	    //1 MODEL_EXP (default): f0+exp(-(U-U0)*Q)
	    //2 MODEL_EXP_NOLIMIT: Same as MODEL_EXP but can reach values > 1
	    void SetSingleGainFluctuationModel(int ModelType) { mModelType = ModelType; };

	    //Same as GetGEMSingleGainFluctuation(int GEMId) but OutOfRange will be 1 if the value for the single gain fluctuation
	    //is outside of the available data region for the fit	    
	    double GetGEMSingleGainFluctuation(int GEMId, int *OutOfRange);	    
	    double GetGEMSingleGainFluctuationCore(int GEMId, int *OutOfRange, int Error);
	
	    //Set the single gain fluctuation for a specific GEM within the stack
	    void SetGEMSingleGainFluctuation(int GEMId, double FixSingleGainFluctuation) { mGEMFixSingleGainFluctuationStatus[GEMId] = 1; mGEMFixSingleGainFluctuation[GEMId] = FixSingleGainFluctuation; }; 
  
	    //Set the general single gain fluctuation coefficients for a specific GEM and GAS Type and Model type
	    //0 MODEL_ARC: 0.5(f0+1)+1/pi(1-f0)arctan(-Q(U-U0))
	    //1 MODEL_EXP (default): f0+exp(-(U-U0)*Q)
	    //2 MODEL_EXP_NOLIMIT: Same as MODEL_EXP but can reach values > 1
	    void SetSingleGainFluctuationCoefficients(int GasType, int GEMType, int ModelType, double f0, double df0, double U0, double dU0, double Q, double dQ);
	    
	    //Get the total effective gain for the full stack
	    double GetTotalEffGain();
	    double GetErrorTotalEffGain();
	//---------------------------
	//Stack properties & analysis
	//---------------------------	
	
	    //Get the energy resolution of the stack (sigma/mu)
	    double GetEnergyResolutionCore(int Error);
	    double GetEnergyResolution();
	    double GetErrorHiEnergyResolution();
	    double GetErrorLoEnergyResolution();
	    
	    //Get number of primary charges in drift region
	    double GetPrimaryCharges() { return mPhotonEnergy/mWi; };
	    double GetErrorPrimaryCharges() { return TMath::Sqrt(mPhotonEnergy/mWi); };
	    
	    //Get number of final charges after extraction of the last amplification stage
	    double GetFinalCharges();
	    double GetErrorFinalCharges();
	    
	    //Get the ion backflow of the stack
	    double GetIonBackflowCore(int Error);
	    double GetIonBackflow();
	    double GetErrorIonBackflow();
	    
	    //Get the epsilon contribution of a single GEM to the ion backflow of the stack
	    double GetGEMEpsilonContribution(int GEMId) { return mGEMEpsilonContribution[GEMId]; };
	    
	    //Get the number of electrons which are collected at a specific GEM stage within the stack
	    double GetGEMElectronContribution(int GEMId) { return mN[GEMId]; };
	    
	    //Returns a histogram of the electrons which are collected at all GEM stages + at the anode
	    TH1F *DrawElectronContribution();
	    
	    //Generates a plot of the single gain fluctuation of a GEM and indicates the working point
	    void DrawGEMSingleGainFluctuation(int GEMId, TCanvas *cPlot);
	    
	    //Generates an electron efficiency plot for a given GEM within the stack and indicates
	    //the working point for used stack settings (voltages and fields)
	    void DrawGEMEfficiencyCurves(int GEMId, TCanvas *cPlot, int PlotElectron, int PlotAvalancheIon, int PlotDriftIon);
	    
	    //Set how the tuning parameters s1, s2 and s3 should be calculated
	    //0 (default): Use polynomial fit and dependence tuningparameter(pitch)
	    //1: Use direct fit results for s1, s2, and s3
	    void SetTuneParameterMethode(int TuneParameterMethode) { mTuneParameterMethode = TuneParameterMethode; };
	    
	
    private:
	
	//Results from model calculations in order to calculate efficiencies
	GEM_Efficiencies_v3 *GEMEfficiency;
	
	//Tuning parameters for the electron collection calculations (TuneParameterID: 0 S1, 1 S2, 2 S3)
	double GetTuneParameter(int TuneParameterID, double Pitch);
	double GetErrorTuneParameter(int TuneParameterID, double Pitch);
	
	double ConvertGEMTypeToPitch(int GEMType);
	
	//Calculate ion collection / extraction (EffType) for avalanche or drift ions (ChargeType) and a specific GEM (GEMId)
	double GetGEMIonEfficiency(int GEMId, int ChargeType, int EffType, int CalcError);

	//Volumes in GEM detector (e.g. induction field, transfer field ...)
	std::vector<double> mVolumeDistance; //mm
	std::vector<double> mVolumeField; //V/cm
	int mVolumeEnties;
	
	//GEMs which are used in the detector
	std::vector<int> mGEMType;
	std::vector<double> mGEMVoltage;
	int mGEMEntries;
	
	//Fix single GEM electron collection / extraction efficiency
	std::vector<int> mGEMFixElectronCollectionStatus;
	std::vector<double> mGEMFixElectronCollection;
	std::vector<int> mGEMFixElectronExtractionStatus; 
	std::vector<double> mGEMFixElectronExtraction;

	//Fix single GEM avalanche ion collection / extraction efficiency
	std::vector<int> mGEMFixAvalancheIonCollectionStatus; 
	std::vector<double> mGEMFixAvalancheIonCollection;
	std::vector<int> mGEMFixAvalancheIonExtractionStatus; 
	std::vector<double> mGEMFixAvalancheIonExtraction;	
	
	//Fix single GEM drift ion collection / extraction efficiency
	std::vector<int> mGEMFixDriftIonCollectionStatus; 
	std::vector<double> mGEMFixDriftIonCollection;
	std::vector<int> mGEMFixDriftIonExtractionStatus; 
	std::vector<double> mGEMFixDriftIonExtraction;	

	//Fix absolute gain for a single GEM
	std::vector<int> mGEMFixAbsGainStatus; 
	std::vector<double> mGEMFixAbsGainSlope;
	std::vector<double> mGEMFixAbsGainConstant;

	//Fix single gain fluctuation for a single GEM 
	std::vector<int> mGEMFixSingleGainFluctuationStatus; 
	std::vector<double> mGEMFixSingleGainFluctuation;


	//Gas which is used in detector
	int mGasType;
	
	//Fano factor for used gas
	double mFano;
	double mErrorFano;
	
	//Mean ionization energy for used gas (values from ALICE-TDR-016)
	double mWi; //eV
	
	//Attachment enables or disabled
	int mAttachmentStatus;
	
	double mAttachmentFactor; //Fixed attachment factor (1/cm)
	
	//Energy of incident photon
	double mPhotonEnergy; //eV
	
	//Number of electrons which are collected at the individual GEM stages
	std::vector<double> mN;
	std::vector<double> mErrorN;
	
	//Gain Coefficients for absolute gain calculations (for exponential fit)   
	//Index 1: GasType 
	//Index 2: GEM Type 
	//Index 3: Coefficients of exponential function (0 slope, 1 constant)
	//Index 4: 0 value, 1 error
	double mGainCoefficients[5][3][2][2];
	
	//Fit range for gain coefficients
	//Index 1: GasType
	//Index 2: GEM Type
	//Index 3: 0 Fit range min, 1 Fit range max
	double mGainCoefficientsRange[5][3][2];
	
	//Allows a scaling of the calculated gain curve (by default this is 1.0, i.e. no scaling)
	double mGainScaleFactor;
	
	//Polynomial coefficients for avalanche / drift ion collection and extraction efficiency
	//Index 1: GasType (0 not used, 1 not used, 2 NeCO2N2 90-10-5, 3 not used, 4 NeCO2 90-10)
	//Index 2: GEM Type 1 (0 Standard, 1 Medium, 2 Large)
	//Index 3: GEM Type 2 (for drift ions)
	//Index 4: 0 Collection, 1 Extraction
	//Index 5: 0 Avalanche, 1 Drift
	//Index 6: Polynomial coefficients
	//Index 7: 0 value, 1 error
	double mGEMIonEfficiencyCoefficients[5][3][3][2][2][8][2];
	
	//Limits for the polynomial fits
	//Index 1 - 5 same as for mGEMIonEfficiencyCoefficientsRange
	//Index 6: 0 Fit range min, fit range max
	double mGEMIonEfficiencyCoefficientsRange[5][3][3][2][2][2];
	
	double mDriftIonCollectionScalingFactor;
	double mDriftIonExtractionScalingFactor;
	
	double mAvalancheIonEfficienciesScalingFactor;
	
	//Save the ion backflow epsilon contribution for each GEM 
	std::vector<double> mGEMEpsilonContribution;
	
	//Tuning parameter polynomial coefficients for S1 S2 and S3
	//Index 1: GasType
	//Index 2: Tune Parameter (0 S1, 1 S2, 2 S3)
	//Index 3: Polynomial coefficients (0 p0, 1 p1, 2 p2)
	//Index 4: 0 Value, 1 Error
	double TuneParameterCoefficients[5][3][3][2];
	
	//Alternative approach to obtain tuning parameters:
	//Use the fit results for s1, s2 and s3 directly fom the fit, i.e. not from the 
	//relation tuningparameter(pitch). So this is only possible for St, MP and LP.
	//By default this way of calculation is deactivated. You need to set this by
	//SetTuneParameterMethode(Methode) where: 
	//Methode: 0 (default): use (polynomial fit) relation tuningparameter(pitch)
	//Methode: 1: use direct fit results for s1, s2 and s3
	//Index 1: GasType
	//Index 2: Tune Parameter (0 S1, 1 S2, 2 S3)
	//Index 3: GEMType
	//Index 4: 0 Value, 1 Error
	double TuneParameter[5][3][3][2];
	
	//Methode how the tuning parameters should be obtained
	//0 (default): Use polynomial fit and use relation tuningparameter(pitch)
	//1: Use direct fit results for s1, s2 and s3
	int mTuneParameterMethode;
	
	//Coefficients for the single gain fluctuation factor f
	//Index 1: Gas Type
	//Index 2: GEM Type
	//Index 3: Model Type (0 MODEL_ARC, 1 MODEL_EXP, 2 MODEL_EXP_NOLIMIT)
	//Index 4: Coefficients (0 f0, 1 U0, 2 Q)
	//Index 5: 0 Value, 1 Error
	double mSingleGainFluctuationCoefficients[5][3][3][3][2];
	
	//Assumed distribution function for the single gain fluctuation
	double SingleGainFluctuationDistribution(double f0, double U0, double Q, double GEMVoltage);
	
	//Used model to calculate single gain fluctuation
	//0 MODEL_ARC : 0.5(f0+1)+1/pi(1-f0)arctan(-Q(U-U0)) 
	//1 MODEL_EXP (default): f0+exp(-(U-U0)*Q)
	//2 MODEL_EXP_NOLIMIT: Same as MODEL_EXP but can reach values > 1
	int mModelType;
	
	//Which model do we use to calculate drift ion efficiencies?
	//MODEL_DRIFT_ION_SIMULATION: Load the values from simulations (default)
	//MODEL_DRIFT_ION_CUT: Cut after plataeau, scaling after
	int mDriftIonModel;

			
};

void GEMStack::DrawGEMSingleGainFluctuation(int GEMId, TCanvas *cPlot)
{
    std::vector<double> SingleGainFluctuationOutOfRange, SingleGainFluctuation;
    std::vector<double> ErrorHiSingleGainFluctuationOutOfRange, ErrorHiSingleGainFluctuation;
    std::vector<double> ErrorLoSingleGainFluctuationOutOfRange, ErrorLoSingleGainFluctuation;
    std::vector<double> GEMVoltageOutOfRange, GEMVoltage;
    int NOutOfRange = 0;
    int N = 0;
    
    double GEMVoltage_Old = GetGEMVoltage(GEMId);
    
    int COUT_MESSAGES_OLD = COUT_MESSAGES;
    int COUT_WARNINGS_OLD = COUT_WARNINGS;
    
    COUT_MESSAGES = 0;
    COUT_WARNINGS = 0;
    
    double UGEM = 100.0;
    int OutOfRange;
    
    while ( UGEM <= 400.0 )
    {
	OutOfRange = 0;
	
	SetGEMVoltage(GEMId, UGEM);
	double SingleGainFluctuationBuffer = GetGEMSingleGainFluctuation(GEMId, &OutOfRange);
	double ErrorHiSingleGainFluctuationBuffer = GetErrorHiGEMSingleGainFluctuation(GEMId);
	double ErrorLoSingleGainFluctuationBuffer = GetErrorLoGEMSingleGainFluctuation(GEMId);

	if ( OutOfRange )
	{
	    SingleGainFluctuationOutOfRange.push_back(SingleGainFluctuationBuffer);
	    ErrorHiSingleGainFluctuationOutOfRange.push_back(ErrorHiSingleGainFluctuationBuffer);
	    ErrorLoSingleGainFluctuationOutOfRange.push_back(ErrorLoSingleGainFluctuationBuffer);
	    GEMVoltageOutOfRange.push_back(UGEM);
	    ++NOutOfRange;
	}else{
	    SingleGainFluctuation.push_back(SingleGainFluctuationBuffer);
	    ErrorHiSingleGainFluctuation.push_back(ErrorHiSingleGainFluctuationBuffer);
	    ErrorLoSingleGainFluctuation.push_back(ErrorLoSingleGainFluctuationBuffer);
	    GEMVoltage.push_back(UGEM);
	    ++N;
	}
	
	UGEM += 1.0;
    }
    
    SetGEMVoltage(GEMId, GEMVoltage_Old);
    
    COUT_MESSAGES = COUT_MESSAGES_OLD;
    COUT_WARNINGS = COUT_WARNINGS_OLD;
    
    TH1F *hr = cPlot->DrawFrame(100.0,0.0,400.0,1.0);    
    hr->SetXTitle("U_{GEM}");
    hr->SetYTitle("Single Gain Fluctuatio f"); 

    TGraphAsymmErrors *grSingleGainFluctuationOutOfRange = new TGraphAsymmErrors(NOutOfRange, &GEMVoltageOutOfRange[0], &SingleGainFluctuationOutOfRange[0], NULL, NULL, &ErrorLoSingleGainFluctuationOutOfRange[0], &ErrorHiSingleGainFluctuationOutOfRange[0]);
    grSingleGainFluctuationOutOfRange->SetLineWidth(2.0);
    grSingleGainFluctuationOutOfRange->SetLineColor(kBlue);
    //grSingleGainFluctuationOutOfRange->SetMarkerStyle(20);
    grSingleGainFluctuationOutOfRange->SetFillColorAlpha(kBlue-10,1.0);
    grSingleGainFluctuationOutOfRange->Draw("3L");
    grSingleGainFluctuationOutOfRange->Draw("Lsamex");

    
    TGraphAsymmErrors *grSingleGainFluctuation = new TGraphAsymmErrors(N, &GEMVoltage[0], &SingleGainFluctuation[0], NULL, NULL, &ErrorLoSingleGainFluctuation[0], &ErrorHiSingleGainFluctuation[0]);
    grSingleGainFluctuation->SetLineWidth(2.0);
    grSingleGainFluctuation->SetLineColor(kBlack);
    grSingleGainFluctuation->SetFillColorAlpha(kGray,1.0);
    grSingleGainFluctuation->Draw("3L");
    grSingleGainFluctuation->Draw("Lsamex");
    
    TLine *lineSingleGainFluctuation = new TLine(GetGEMVoltage(GEMId),GetGEMSingleGainFluctuation(GEMId)-0.02,GetGEMVoltage(GEMId),GetGEMSingleGainFluctuation(GEMId)+0.02);
    lineSingleGainFluctuation->SetLineColor(kBlue);
    lineSingleGainFluctuation->SetLineWidth(7);
    lineSingleGainFluctuation->Draw("el");
}

TH1F *GEMStack::DrawElectronContribution()
{
    TH1F *th = new TH1F("th","Electron Contribution",mGEMEntries+1,0,mGEMEntries+1);
    
    for ( int n = 0; n <= mGEMEntries; ++n )
    {
	th->SetBinContent(n+1,mN[n]);

    }
    
    th->SetFillColor(kBlue);
    th->SetXTitle("GEMId");
    th->SetYTitle("Number of collected electrons");
    
    return th;
}

void GEMStack::DrawGEMEfficiencyCurves(int GEMId, TCanvas *cPlot, int PlotElectron, int PlotAvalancheIon, int PlotDriftIon)
{
    std::vector<double> ElectronCollection, ElectronExtraction, ErrorElectronCollection, ErrorElectronExtraction;
    std::vector<double> AvalancheIonCollection, AvalancheIonExtraction, ErrorAvalancheIonCollection, ErrorAvalancheIonExtraction;
    std::vector<double> DriftIonCollection, DriftIonExtraction, ErrorDriftIonCollection, ErrorDriftIonExtraction;
    std::vector<double> Eta;
    int NElectron = 0;
    int NAvalancheIon = 0;
    int NDriftIon = 0;
    double EExtern = 0.0;

    //Is this the last GEM in the stack and do we want to plot drift ion efficiencies? If so then there are no drift ions.
    if ( GEMId == mGEMEntries-1 && PlotDriftIon )
    {
	if ( COUT_WARNINGS ) std::cerr << "Warning: Trying to plot drift ion efficiencies for last GEM with GEMId=" << GEMId << " (there are no drift ions)!" << std::endl;
    }
    
    //Save old values since we change them to create the curves
    double VolumeFieldAbove = GetVolumeField(GEMId);
    double VolumeFieldBelow = GetVolumeField(GEMId+1);
    
    int COUT_WARNINGS_OLD = COUT_WARNINGS;
    int COUT_MESSAGES_OLD = COUT_MESSAGES;
    
    //Disable messages and errors for the plot curve generation
    COUT_WARNINGS = 0;
    COUT_MESSAGES = 0;

    //Generate the curves
    while ( EExtern <= 15000.0 )
    {
	//Collection and extraction are not functions of the opposite fields, that is why we can change 
	//EBelow and EBelow at the same time.
	SetVolumeField(GEMId, EExtern);
	SetVolumeField(GEMId+1, EExtern);

	if ( PlotElectron )
	{
	    ElectronCollection.push_back(GetGEMElectronCollection(GEMId));
	    ElectronExtraction.push_back(GetGEMElectronExtraction(GEMId));
	    ErrorElectronCollection.push_back(GetErrorGEMElectronCollection(GEMId));
	    ErrorElectronExtraction.push_back(GetErrorGEMElectronExtraction(GEMId));
	    ++NElectron;
	}
	
	if ( PlotAvalancheIon )
	{
	    AvalancheIonCollection.push_back(GetGEMAvalancheIonCollection(GEMId));
	    AvalancheIonExtraction.push_back(GetGEMAvalancheIonExtraction(GEMId));
	    ErrorAvalancheIonCollection.push_back(GetErrorGEMAvalancheIonCollection(GEMId));
	    ErrorAvalancheIonExtraction.push_back(GetErrorGEMAvalancheIonExtraction(GEMId));
	    ++NAvalancheIon;
	}
	
	if ( GEMId < mGEMEntries-1 && PlotDriftIon )
	{
	    DriftIonCollection.push_back(GetGEMDriftIonCollection(GEMId));
	    DriftIonExtraction.push_back(GetGEMDriftIonExtraction(GEMId));
	    ErrorDriftIonCollection.push_back(GetErrorGEMDriftIonCollection(GEMId));
	    ErrorDriftIonExtraction.push_back(GetErrorGEMDriftIonExtraction(GEMId));
	    ++NDriftIon;
	}
	
	Eta.push_back(EExtern/(GetGEMVoltage(GEMId)/0.005));

	EExtern += 250.0;
    }  
    
    TH1F *hr = cPlot->DrawFrame(0.0,0.0,0.25,1.0);    
    hr->SetXTitle("E_{Extern}/E_{GEM}");
    hr->SetYTitle("Efficiencies"); 
    
    if ( PlotElectron )
    {
	
	//TGraphErrors *grElectronCollection = new TGraphErrors(NElectron, &Eta[0], &ElectronCollection[0], NULL, &ErrorElectronCollection[0]);
	TGraphErrors *grElectronCollection = new TGraphErrors(NElectron, &Eta[0], &ElectronCollection[0], NULL, NULL);
	grElectronCollection->SetLineWidth(3);
	grElectronCollection->SetLineColor(kBlue);
	grElectronCollection->SetFillColorAlpha(kBlue-10,0.5);
	
	//TGraphErrors *grElectronExtraction = new TGraphErrors(NElectron, &Eta[0], &ElectronExtraction[0], NULL, &ErrorElectronExtraction[0]);
	TGraphErrors *grElectronExtraction = new TGraphErrors(NElectron, &Eta[0], &ElectronExtraction[0], NULL, NULL);
	grElectronExtraction->SetLineWidth(3);
	grElectronExtraction->SetLineColor(kBlue);
	grElectronExtraction->SetLineStyle(2);
	grElectronExtraction->SetFillColorAlpha(kBlue-10,0.5);
	
	//grElectronCollection->Draw("3L");   
	grElectronCollection->Draw("l");
	//grElectronExtraction->Draw("3L"); 
	grElectronExtraction->Draw("l");
	
    }
    
    if ( PlotAvalancheIon )
    {
	//TGraphErrors *grAvalancheIonCollection = new TGraphErrors(NAvalancheIon, &Eta[0], &AvalancheIonCollection[0], NULL, &ErrorAvalancheIonCollection[0]);
	TGraphErrors *grAvalancheIonCollection = new TGraphErrors(NAvalancheIon, &Eta[0], &AvalancheIonCollection[0], NULL, NULL);
	grAvalancheIonCollection->SetLineWidth(3);
	grAvalancheIonCollection->SetLineColor(kBlack);  
	grAvalancheIonCollection->SetFillColorAlpha(kBlack-10,0.5);
	
	//TGraphErrors *grAvalancheIonExtraction = new TGraphErrors(NAvalancheIon, &Eta[0], &AvalancheIonExtraction[0], NULL, &ErrorAvalancheIonExtraction[0]);
	TGraphErrors *grAvalancheIonExtraction = new TGraphErrors(NAvalancheIon, &Eta[0], &AvalancheIonExtraction[0], NULL, NULL);
	grAvalancheIonExtraction->SetLineWidth(3);
	grAvalancheIonExtraction->SetLineColor(kBlack);
	grAvalancheIonExtraction->SetLineStyle(2);
	grAvalancheIonExtraction->SetFillColorAlpha(kBlack-10,0.5);
    
	//grAvalancheIonCollection->Draw("3L");
	//grAvalancheIonExtraction->Draw("3L");
	grAvalancheIonCollection->Draw("l");
	grAvalancheIonExtraction->Draw("l");
    }
    
    if ( GEMId < mGEMEntries-1 && PlotDriftIon )
    {
	//TGraphErrors *grDriftIonCollection = new TGraphErrors(NDriftIon, &Eta[0], &DriftIonCollection[0], NULL, &ErrorDriftIonCollection[0]);
	TGraphErrors *grDriftIonCollection = new TGraphErrors(NDriftIon, &Eta[0], &DriftIonCollection[0], NULL, NULL);
	grDriftIonCollection->SetLineWidth(3);
	grDriftIonCollection->SetLineColor(kRed);  
	grDriftIonCollection->SetFillColorAlpha(kRed-10,0.5);
	
	//TGraphErrors *grDriftIonExtraction = new TGraphErrors(NDriftIon, &Eta[0], &DriftIonExtraction[0], NULL, &ErrorDriftIonExtraction[0]);
	TGraphErrors *grDriftIonExtraction = new TGraphErrors(NDriftIon, &Eta[0], &DriftIonExtraction[0], NULL, NULL);
	grDriftIonExtraction->SetLineWidth(3);
	grDriftIonExtraction->SetLineColor(kRed);
	grDriftIonExtraction->SetLineStyle(2);
	grDriftIonExtraction->SetFillColorAlpha(kRed-10,0.5);
    
	//grDriftIonCollection->Draw("3L");
	//grDriftIonExtraction->Draw("3L");
	grDriftIonCollection->Draw("l");
	grDriftIonExtraction->Draw("l");
    }
    
    //Restore old values
    SetVolumeField(GEMId, VolumeFieldAbove);
    SetVolumeField(GEMId+1, VolumeFieldBelow);
    
    COUT_WARNINGS = COUT_WARNINGS_OLD;
    COUT_MESSAGES = COUT_MESSAGES_OLD;
    
    double DeltaEfficiency = 0.02;
    
    int DrawWorkingPoint = 1;
    
    //Indicate where the GEM is operating
    if ( PlotElectron && DrawWorkingPoint )
    {
	TLine *lineElectronCollection = new TLine(VolumeFieldAbove/(GetGEMVoltage(GEMId)/0.005),GetGEMElectronCollection(GEMId)-DeltaEfficiency,VolumeFieldAbove/(GetGEMVoltage(GEMId)/0.005),GetGEMElectronCollection(GEMId)+DeltaEfficiency);
	lineElectronCollection->SetLineColor(kBlue);
	lineElectronCollection->SetLineWidth(7);
	lineElectronCollection->Draw();
	
	TLine *lineElectronExtraction = new TLine(VolumeFieldBelow/(GetGEMVoltage(GEMId)/0.005),GetGEMElectronExtraction(GEMId)-DeltaEfficiency,VolumeFieldBelow/(GetGEMVoltage(GEMId)/0.005),GetGEMElectronExtraction(GEMId)+DeltaEfficiency);
	lineElectronExtraction->SetLineColor(kBlue);
	lineElectronExtraction->SetLineWidth(7);
	lineElectronExtraction->Draw();
    }
    
    if ( PlotAvalancheIon && DrawWorkingPoint)
    {
	TLine *lineAvalancheIonCollection = new TLine(VolumeFieldBelow/(GetGEMVoltage(GEMId)/0.005),GetGEMAvalancheIonCollection(GEMId)-DeltaEfficiency+0.005,VolumeFieldBelow/(GetGEMVoltage(GEMId)/0.005),GetGEMAvalancheIonCollection(GEMId)+DeltaEfficiency-0.005);
	lineAvalancheIonCollection->SetLineColor(kBlack);
	lineAvalancheIonCollection->SetLineWidth(7);
	lineAvalancheIonCollection->Draw();
	
	TLine *lineAvalancheIonExtraction = new TLine(VolumeFieldAbove/(GetGEMVoltage(GEMId)/0.005),GetGEMAvalancheIonExtraction(GEMId)-DeltaEfficiency+0.005,VolumeFieldAbove/(GetGEMVoltage(GEMId)/0.005),GetGEMAvalancheIonExtraction(GEMId)+DeltaEfficiency-0.005);
	lineAvalancheIonExtraction->SetLineColor(kBlack);
	lineAvalancheIonExtraction->SetLineWidth(7);
	lineAvalancheIonExtraction->Draw();    
    }
    
    if ( GEMId < mGEMEntries-1 && PlotDriftIon && DrawWorkingPoint )
    {
	TLine *lineDriftIonCollection = new TLine(VolumeFieldBelow/(GetGEMVoltage(GEMId)/0.005),GetGEMDriftIonCollection(GEMId)-DeltaEfficiency+0.010,VolumeFieldBelow/(GetGEMVoltage(GEMId)/0.005),GetGEMDriftIonCollection(GEMId)+DeltaEfficiency-0.010);
	lineDriftIonCollection->SetLineColor(kRed);
	lineDriftIonCollection->SetLineWidth(7);
	lineDriftIonCollection->Draw();
	
	TLine *lineDriftIonExtraction = new TLine(VolumeFieldAbove/(GetGEMVoltage(GEMId)/0.005),GetGEMDriftIonExtraction(GEMId)-DeltaEfficiency+0.010,VolumeFieldAbove/(GetGEMVoltage(GEMId)/0.005),GetGEMDriftIonExtraction(GEMId)+DeltaEfficiency-0.010);
	lineDriftIonExtraction->SetLineColor(kRed);
	lineDriftIonExtraction->SetLineWidth(7);
	lineDriftIonExtraction->Draw();      
    }
}

//Default contructor
GEMStack::GEMStack()
{
    GEMEfficiency = new GEM_Efficiencies_v3(GEM_Efficiencies_N);
    
    SetAttachment(0); //default no attachment
    SetAttachmentFactor(0.13);
    
    mVolumeEnties = 0;
    mGEMEntries = 0;

    //Tuning parameter polynomial coefficients for S1 S2 and S3
    mTuneParameterMethode = 0; //By default we use the polynomial fits which give the relation tuningparameter(pitch)
    
    TuneParameterCoefficients[GAS_NECO2N2][0][0][0] = -1.80809; 	//NeCO2N2 S1 p0
    TuneParameterCoefficients[GAS_NECO2N2][0][0][1] = 0.0850362; 	//NeCO2N2 S1 dp0
    TuneParameterCoefficients[GAS_NECO2N2][0][1][0] = 0.047546; 	//NeCO2N2 S1 p1
    TuneParameterCoefficients[GAS_NECO2N2][0][1][1] = 0.00113236; 	//NeCO2N2 S1 dp1
    TuneParameterCoefficients[GAS_NECO2N2][0][2][0] = 9.82832e-05; 	//NeCO2N2 S1 p2
    TuneParameterCoefficients[GAS_NECO2N2][0][2][1] = 3.39821e-06; 	//NeCO2N2 S1 dp2
    
    TuneParameterCoefficients[GAS_NECO2N2][1][0][0] = -3.8041; 		//NeCO2N2 S2 p0
    TuneParameterCoefficients[GAS_NECO2N2][1][0][1] = 0.44283; 		//NeCO2N2 S2 dp0
    TuneParameterCoefficients[GAS_NECO2N2][1][1][0] = 0.0941137; 	//NeCO2N2 S2 p1
    TuneParameterCoefficients[GAS_NECO2N2][1][1][1] = 0.00534736; 	//NeCO2N2 S2 dp1
    TuneParameterCoefficients[GAS_NECO2N2][1][2][0] = -0.00015573; 	//NeCO2N2 S2 p2
    TuneParameterCoefficients[GAS_NECO2N2][1][2][1] = 1.29112e-05; 	//NeCO2N2 S2 dp2
    
    TuneParameterCoefficients[GAS_NECO2N2][2][0][0] = 1.29655; 		//NeCO2N2 S3 p0
    TuneParameterCoefficients[GAS_NECO2N2][2][0][1] = 0.00212428; 	//NeCO2N2 S3 dp0
    TuneParameterCoefficients[GAS_NECO2N2][2][1][0] = 0.0;	 	//NeCO2N2 S3 p1
    TuneParameterCoefficients[GAS_NECO2N2][2][1][1] = 0.0;	 	//NeCO2N2 S3 dp1
    TuneParameterCoefficients[GAS_NECO2N2][2][2][0] = 0.0;	 	//NeCO2N2 S3 p2
    TuneParameterCoefficients[GAS_NECO2N2][2][2][1] = 0.0;	 	//NeCO2N2 S3 dp2
 
    TuneParameterCoefficients[GAS_NECO2][0][0][0] = -1.21223; 		//NeCO2 S1 p0
    TuneParameterCoefficients[GAS_NECO2][0][0][1] = 0.659692; 		//NeCO2 S1 dp0
    TuneParameterCoefficients[GAS_NECO2][0][1][0] = 0.0523091; 		//NeCO2 S1 p1
    TuneParameterCoefficients[GAS_NECO2][0][1][1] = 0.0072931; 		//NeCO2 S1 dp1
    TuneParameterCoefficients[GAS_NECO2][0][2][0] = 7.80905e-05; 	//NeCO2 S1 p2
    TuneParameterCoefficients[GAS_NECO2][0][2][1] = 1.90625e-05; 	//NeCO2 S1 dp2
    
    TuneParameterCoefficients[GAS_NECO2][1][0][0] = 4.27833; 		//NeCO2 S2 p0
    TuneParameterCoefficients[GAS_NECO2][1][0][1] = 4.14004; 		//NeCO2 S2 dp0
    TuneParameterCoefficients[GAS_NECO2][1][1][0] = 0.0479822; 		//NeCO2 S2 p1
    TuneParameterCoefficients[GAS_NECO2][1][1][1] = 0.0434204; 		//NeCO2 S2 dp1
    TuneParameterCoefficients[GAS_NECO2][1][2][0] = -0.000103869; 	//NeCO2 S2 p2
    TuneParameterCoefficients[GAS_NECO2][1][2][1] = 0.00010573; 	//NeCO2 S2 dp2
 
    TuneParameterCoefficients[GAS_NECO2][2][0][0] = 1.34333; 		//NeCO2 S3 p0
    TuneParameterCoefficients[GAS_NECO2][2][0][1] = 0.0057735; 		//NeCO2 S3 dp0
    TuneParameterCoefficients[GAS_NECO2][2][1][0] = 0.0; 		//NeCO2 S3 p1
    TuneParameterCoefficients[GAS_NECO2][2][1][1] = 0.0; 		//NeCO2 S3 dp1
    TuneParameterCoefficients[GAS_NECO2][2][2][0] = 0.0; 		//NeCO2 S3 p2
    TuneParameterCoefficients[GAS_NECO2][2][2][1] = 0.0; 		//NeCO2 S3 dp2
    
    TuneParameter[GAS_NECO2N2][0][GEM_ST][0] = 6.78972; 	//NeCO2N2 S1 St value
    TuneParameter[GAS_NECO2N2][0][GEM_ST][1] = 0.0234872; 	//NeCO2N2 S1 St error
    TuneParameter[GAS_NECO2N2][1][GEM_ST][0] = 6.88737; 	//NeCO2N2 S2 St value
    TuneParameter[GAS_NECO2N2][1][GEM_ST][1] = 0.196986; 	//NeCO2N2 S2 St error
    TuneParameter[GAS_NECO2N2][2][GEM_ST][0] = 1.30061; 	//NeCO2N2 S3 St value
    TuneParameter[GAS_NECO2N2][2][GEM_ST][1] = 0.00516626; 	//NeCO2N2 S3 St error
    
    TuneParameter[GAS_NECO2N2][0][GEM_MP][0] = 11.5214; 	//NeCO2N2 S1 MP value
    TuneParameter[GAS_NECO2N2][0][GEM_MP][1] = 0.0554358; 	//NeCO2N2 S1 MP error
    TuneParameter[GAS_NECO2N2][1][GEM_MP][0] = 8.98337; 	//NeCO2N2 S2 MP value
    TuneParameter[GAS_NECO2N2][1][GEM_MP][1] = 0.327606; 	//NeCO2N2 S2 MP error
    TuneParameter[GAS_NECO2N2][2][GEM_MP][0] = 1.30285; 	//NeCO2N2 S3 MP value
    TuneParameter[GAS_NECO2N2][2][GEM_MP][1] = 0.00639796; 	//NeCO2N2 S3 MP error
    
    TuneParameter[GAS_NECO2N2][0][GEM_LP][0] = 18.788; 		//NeCO2N2 S1 LP value
    TuneParameter[GAS_NECO2N2][0][GEM_LP][1] = 0.143816; 	//NeCO2N2 S1 LP error
    TuneParameter[GAS_NECO2N2][1][GEM_LP][0] = 9.90459; 	//NeCO2N2 S2 LP value
    TuneParameter[GAS_NECO2N2][1][GEM_LP][1] = 0.417301; 	//NeCO2N2 S2 LP error
    TuneParameter[GAS_NECO2N2][2][GEM_LP][0] = 1.30125; 	//NeCO2N2 S3 LP value
    TuneParameter[GAS_NECO2N2][2][GEM_LP][1] = 0.00738531; 	//NeCO2N2 S3 LP error
    
    TuneParameter[GAS_NECO2][0][GEM_ST][0] = 7.64161; 	//NeCO2 S1 St value
    TuneParameter[GAS_NECO2][0][GEM_ST][1] = 0.03; 	//NeCO2 S1 St error
    TuneParameter[GAS_NECO2][1][GEM_ST][0] = 8.96; 	//NeCO2 S2 St value
    TuneParameter[GAS_NECO2][1][GEM_ST][1] = 0.28; 	//NeCO2 S2 St error
    TuneParameter[GAS_NECO2][2][GEM_ST][0] = 1.37; 	//NeCO2 S3 St value
    TuneParameter[GAS_NECO2][2][GEM_ST][1] = 0.01; 	//NeCO2 S3 St error
    
    TuneParameter[GAS_NECO2][0][GEM_MP][0] = 12.3732; 	//NeCO2 S1 MP value
    TuneParameter[GAS_NECO2][0][GEM_MP][1] = 0.0606857; 	//NeCO2 S1 MP error
    TuneParameter[GAS_NECO2][1][GEM_MP][0] = 9.72; 	//NeCO2 S2 MP value
    TuneParameter[GAS_NECO2][1][GEM_MP][1] = 0.42; 	//NeCO2 S2 MP error
    TuneParameter[GAS_NECO2][2][GEM_MP][0] = 1.34; 	//NeCO2 S3 MP value
    TuneParameter[GAS_NECO2][2][GEM_MP][1] = 0.01; 	//NeCO2 S3 MP error
    
    TuneParameter[GAS_NECO2][0][GEM_LP][0] = 19.5566; 		//NeCO2 S1 LP value
    TuneParameter[GAS_NECO2][0][GEM_LP][1] = 0.154699; 	//NeCO2 S1 LP error
    TuneParameter[GAS_NECO2][1][GEM_LP][0] = 9.57; 	//NeCO2 S2 LP value
    TuneParameter[GAS_NECO2][1][GEM_LP][1] = 0.55; 	//NeCO2 S2 LP error
    TuneParameter[GAS_NECO2][2][GEM_LP][0] = 1.32; 	//NeCO2 S3 LP value
    TuneParameter[GAS_NECO2][2][GEM_LP][1] = 0.01; 	//NeCO2 S3 LP error
 
    //Gain curves (Eff. from model) //Slop Const
    SetAbsGain(GAS_NECO2N2, GEM_ST, 0.0183423, -1.91668, 1.99349e-05, 0.00614254, 200.0, 400.0); 
    SetAbsGain(GAS_NECO2N2, GEM_MP, 0.0185194, -1.95479, 1.98767e-05, 0.00610268, 200.0, 400.0); 
    SetAbsGain(GAS_NECO2N2, GEM_LP, 0.0186415, -1.98842, 2.50665e-05, 0.00788618, 200.0, 400.0);

    //TODO: Implement new gain simulations for NeCO2 (these values are outdated) FOR MP AND LP! (St is updated)
    //SetAbsGain(GAS_NECO2, GEM_ST, 0.0204875, -2.1974, 0.0, 0.0, 160.0, 400.0); //Old Gain curve
    SetAbsGain(GAS_NECO2, GEM_ST, 0.0198924, -1.75538, 1.99934e-05, 0.00616886, 200.0, 400.0);
    SetAbsGain(GAS_NECO2, GEM_MP, 0.0211416, -2.2612, 0.0, 0.0, 160.0, 400.0);
    SetAbsGain(GAS_NECO2, GEM_LP, 0.0196198, -1.41057, 0.0, 0.0, 160.0, 400.0);
    
    //Single gain fluctuation factor f (from fits to absolute gain) arctan model
    //SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_ST, MODEL_ARC, 0.450676, 0.00168362, 210.0, 25.0, 0.10, 0.0);
    //SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_MP, MODEL_ARC, 0.457851, 0.00171114, 210.0, 25.0, 0.10, 0.0);
    //SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_LP, MODEL_ARC, 0.465322, 0.00212891, 210.0, 25.0, 0.10, 0.0);
    //SetSingleGainFluctuationCoefficients(GAS_NECO2, GEM_ST, MODEL_ARC, 0.445784, 0.00145876, 210.0, 25.0, 0.10, 0.0);
    //SetSingleGainFluctuationCoefficients(GAS_NECO2, GEM_MP, MODEL_ARC, 0.0, 0.0, 210.0, 25.0, 0.10, 0.0);
    //SetSingleGainFluctuationCoefficients(GAS_NECO2, GEM_LP, MODEL_ARC, 0.0, 0.0, 210.0, 25.0, 0.10, 0.0);

    //Single gain fluctuation factor f (from fits to effective gain, VIK_gaincurve_eff.cpp) arctan model
    SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_ST, MODEL_ARC, 0.541429, 0.00172418, 230.0, 25.0, 0.10, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_MP, MODEL_ARC, 0.549837, 0.00175462, 230.0, 25.0, 0.10, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_LP, MODEL_ARC, 0.554301, 0.00214115, 230.0, 25.0, 0.10, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2, GEM_ST, MODEL_ARC, 0.528631, 0.00166785, 210.0, 25.0, 0.10, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2, GEM_MP, MODEL_ARC, 0.531467, 0.00168922, 210.0, 25.0, 0.10, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2, GEM_LP, MODEL_ARC, 0.541278, 0.00207736, 210.0, 25.0, 0.10, 0.0);
    
    SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_ST, MODEL_EXP, 0.541429, 0.00172418, 190.0, 25.0, 0.1, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_MP, MODEL_EXP, 0.549837, 0.00175462, 190.0, 25.0, 0.1, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_LP, MODEL_EXP, 0.554301, 0.00214115, 190.0, 25.0, 0.1, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2, GEM_ST, MODEL_EXP, 0.528631, 0.00166785, 140.0, 25.0, 0.027, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2, GEM_MP, MODEL_EXP, 0.531467, 0.00168922, 140.0, 25.0, 0.027, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2, GEM_LP, MODEL_EXP, 0.541278, 0.00207736, 140.0, 25.0, 0.027, 0.0);
    
    //Config 1
    //SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_ST, MODEL_EXP_NOLIMIT, 0.541429, 0.00172418, 220.0, 25.0, 0.02, 0.0);
    //SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_MP, MODEL_EXP_NOLIMIT, 0.549837, 0.00175462, 220.0, 25.0, 0.02, 0.0);
    //SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_LP, MODEL_EXP_NOLIMIT, 0.554301, 0.00214115, 220.0, 25.0, 0.02, 0.0);
    
    //Config 2
    //SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_ST, MODEL_EXP_NOLIMIT, 0.541429, 0.00172418, 200.0, 25.0, 0.1, 0.0);
    //SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_MP, MODEL_EXP_NOLIMIT, 0.549837, 0.00175462, 200.0, 25.0, 0.1, 0.0);
    //SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_LP, MODEL_EXP_NOLIMIT, 0.554301, 0.00214115, 200.0, 25.0, 0.1, 0.0);
    
    //Config 3
    SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_ST, MODEL_EXP_NOLIMIT, 0.541429, 0.00172418, 230.0, 25.0, 0.03, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_MP, MODEL_EXP_NOLIMIT, 0.549837, 0.00175462, 230.0, 25.0, 0.03, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2N2, GEM_LP, MODEL_EXP_NOLIMIT, 0.554301, 0.00214115, 230.0, 25.0, 0.03, 0.0);
    
    SetSingleGainFluctuationCoefficients(GAS_NECO2, GEM_ST, MODEL_EXP_NOLIMIT, 0.528631, 0.00166785, 160.0, 25.0, 0.027, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2, GEM_MP, MODEL_EXP_NOLIMIT, 0.531467, 0.00168922, 160.0, 25.0, 0.027, 0.0);
    SetSingleGainFluctuationCoefficients(GAS_NECO2, GEM_LP, MODEL_EXP_NOLIMIT, 0.541278, 0.00207736, 160.0, 25.0, 0.027, 0.0);
   
    //Exp model to calculate single gain fluctuation as default
    mModelType = MODEL_EXP;
    
    //Avalanche ion efficiencies
    SetAvalancheIonCollection(GAS_NECO2,GEM_ST,0.999762,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.064);
    SetErrorAvalancheIonCollection(GAS_NECO2,GEM_ST,0.000107637,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    //SetAvalancheIonExtraction(GAS_NECO2,GEM_ST,0.0199939,22.7327,-182.224,-22.6817,0.0,0.0,0.0,0.0,0.0,0.064);
    SetAvalancheIonExtraction(GAS_NECO2,GEM_ST,0.0253249,21.6261,-132.257,-640.726,0.0,0.0,0.0,0.0,0.0,0.060);
    SetErrorAvalancheIonExtraction(GAS_NECO2,GEM_ST,0.0112062,1.56692,59.4524,646.926,0.0,0.0,0.0,0.0);
    
    SetAvalancheIonCollection(GAS_NECO2,GEM_MP,0.99954,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.064);
    SetErrorAvalancheIonCollection(GAS_NECO2,GEM_MP,0.000198875,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    SetAvalancheIonExtraction(GAS_NECO2,GEM_MP,0.00767916,55.1862,-2036.81,39477.7,-370641,1.32788e+06,0.0,0.0,0.0,0.07);
    SetErrorAvalancheIonExtraction(GAS_NECO2,GEM_MP,0.0111851,2.85861,229.495,7734.31,114896,620562,0.0,0.0);
    
    SetAvalancheIonCollection(GAS_NECO2N2,GEM_ST,0.999621,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.08);
    SetErrorAvalancheIonCollection(GAS_NECO2N2,GEM_ST,0.000146674,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    SetAvalancheIonExtraction(GAS_NECO2N2,GEM_ST,0.0621603,15.0545,395.82,-18531.1,242227,-1.0615e+06,0.0,0.0,0.0,0.08);
    SetErrorAvalancheIonExtraction(GAS_NECO2N2,GEM_ST,0.0117168,2.69559,194.113,5851,77649.8,374772,0.0,0.0);
    
    SetAvalancheIonCollection(GAS_NECO2N2,GEM_MP,0.999442,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.08);
    SetErrorAvalancheIonCollection(GAS_NECO2N2,GEM_MP,0.000208842,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    SetAvalancheIonExtraction(GAS_NECO2N2,GEM_MP,0.0290267,51.6565,-1860.61,35805.5,-340681,1.26835e+06,0.0,0.0,0.0,0.08);
    SetErrorAvalancheIonExtraction(GAS_NECO2N2,GEM_MP,0.00775146,1.72954,120.695,3522.94,45255.9,211374,0.0,0.0);
    
    SetAvalancheIonCollection(GAS_NECO2N2,GEM_LP,0.999376,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.08);
    SetErrorAvalancheIonCollection(GAS_NECO2N2,GEM_LP,0.000235986,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    //pol6 0 - 0.09
    //SetAvalancheIonExtraction(GAS_NECO2N2,GEM_LP,0.13123,52.8226,-2789.52,88155.2,-1.54013e+06,1.3692e+07,-4.82221e+07,0.0,0.0,0.08);
    //SetErrorAvalancheIonExtraction(GAS_NECO2N2,GEM_LP,0.0119363,3.54255,340.252,14514,304857,3.08717e+06,1.20333e+07,0.0);
    //pol5 0 - 0.08
    SetAvalancheIonExtraction(GAS_NECO2N2,GEM_LP,0.158632,42.4768,-1666.55,36498.1,-397872,1.68003e+06,0.0,0.0,0.0,0.08);
    SetErrorAvalancheIonExtraction(GAS_NECO2N2,GEM_LP,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);    
    
    //Drift ion efficiencies
    SetDriftIonCollection(GAS_NECO2,GEM_ST,GEM_ST,1.00801,-2.53942,178.708,-4505.92,44285.9,-212679,502728,-469370,0.0,0.25);
    SetErrorDriftIonCollection(GAS_NECO2,GEM_ST,GEM_ST,0.00309738,0.501365,22.4847,433.669,4231.6,21884.2,57177.5,59372.5);
    SetDriftIonExtraction(GAS_NECO2,GEM_ST,GEM_ST,0.969835,3.09599,-106.075,1701.71,-14453,66915.1,-159509,153177,0.0,0.25);
    SetErrorDriftIonExtraction(GAS_NECO2,GEM_ST,GEM_ST,0.00314911,0.39795,16.2232,297.56,2814.31,14250.6,36673.7,37652.1);
    
    SetDriftIonCollection(GAS_NECO2N2,GEM_ST,GEM_ST,0.994521,3.5591,-231.682,2075.77,-7475.91,9701.22,0.0,0.0,0.0,0.25);
    SetErrorDriftIonCollection(GAS_NECO2N2,GEM_ST,GEM_ST,0.00503649,0.454507,11.0506,105.409,430.437,629.849,0.0,0.0);
    SetDriftIonExtraction(GAS_NECO2N2,GEM_ST,GEM_ST,0.964093,3.68216,-126.101,2022.42,-17173.8,79501.8,-189494,181958,0.0,0.25);
    SetErrorDriftIonExtraction(GAS_NECO2N2,GEM_ST,GEM_ST,0.00372995,0.47135,19.2155,352.444,3333.4,16879.1,43438,44596.9);
    
    SetDriftIonCollection(GAS_NECO2N2,GEM_ST,GEM_LP,0.996456,3.26179,-225.997,2032.41,-7322.73,9493.02,0.0,0.0,0.0,0.25);
    SetErrorDriftIonCollection(GAS_NECO2N2,GEM_ST,GEM_LP,0.00536217,0.483897,11.7652,112.225,458.271,670.578,0.0,0.0);
    SetDriftIonExtraction(GAS_NECO2N2,GEM_ST,GEM_LP,0.961471,3.96089,-135.827,2180.18,-18523.6,85783.5,-204527,196438,0.0,0.25);
    SetErrorDriftIonExtraction(GAS_NECO2N2,GEM_ST,GEM_LP,0.00410894,0.519243,21.1679,388.255,3672.1,18594.1,47851.7,49128.3);
    
    SetDriftIonCollection(GAS_NECO2N2,GEM_LP,GEM_ST,1.191,-80.6387,2133.22,-28806.2,216744,-918111,2.04557e+06,-1.86331e+06,0.0,0.25);
    SetErrorDriftIonCollection(GAS_NECO2N2,GEM_LP,GEM_ST,0.00878856,1.11634,45.3783,831.617,7869.08,39883.6,102746,105591);
    SetDriftIonExtraction(GAS_NECO2N2,GEM_LP,GEM_ST,0.999999,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.25);
    SetErrorDriftIonExtraction(GAS_NECO2N2,GEM_LP,GEM_ST,1.11093e-06,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    
    SetDriftIonCollection(GAS_NECO2N2,GEM_LP,GEM_LP,1.18765,-80.7014,2143.03,-29028,218927,-929025,2.07277e+06,-1.89015e+06,0.0,0.25);
    SetErrorDriftIonCollection(GAS_NECO2N2,GEM_LP,GEM_LP,0.00885388,1.12382,45.5672,833.693,7885.1,39972.7,103023,105930);
    SetDriftIonExtraction(GAS_NECO2N2,GEM_LP,GEM_LP,0.999998,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.25);
    SetErrorDriftIonExtraction(GAS_NECO2N2,GEM_LP,GEM_LP,2.19527e-06,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    
    mDriftIonCollectionScalingFactor = 1.0;
    mDriftIonExtractionScalingFactor = 1.0;
    
    mDriftIonModel = MODEL_DRIFT_ION_SIMULATION;
    
    mAvalancheIonEfficienciesScalingFactor = 1.0;
    
    SetAbsGainScalingFactor(1.0);
    
}

//Default destructor
GEMStack::~GEMStack()
{
    delete GEMEfficiency;
}

int GEMStack::AddVolume(double VolumeDistance, double VolumeField)
{
    mVolumeDistance.push_back(VolumeDistance);
    mVolumeField.push_back(VolumeField);
    ++mVolumeEnties;
    
    return mVolumeEnties-1;
}


int GEMStack::AddGEM(int GEMType, double GEMVoltage)
{
    mGEMType.push_back(GEMType);
    mGEMVoltage.push_back(GEMVoltage);
    
    mGEMFixElectronCollectionStatus.push_back(0); //0 -> no fixed values by default
    mGEMFixElectronCollection.push_back(0.0);
    mGEMFixElectronExtractionStatus.push_back(0); //0 -> no fixed values by default
    mGEMFixElectronExtraction.push_back(0.0);
    
    mGEMFixAvalancheIonCollectionStatus.push_back(0); //0 -> no fixed values by default
    mGEMFixAvalancheIonCollection.push_back(0.0);
    mGEMFixAvalancheIonExtractionStatus.push_back(0); //0 -> no fixed values by default
    mGEMFixAvalancheIonExtraction.push_back(0.0);
    
    mGEMFixDriftIonCollectionStatus.push_back(0); //0 -> no fixed values by default
    mGEMFixDriftIonCollection.push_back(0.0);
    mGEMFixDriftIonExtractionStatus.push_back(0); //0 -> no fixed values by default
    mGEMFixDriftIonExtraction.push_back(0.0);
    
    mGEMFixAbsGainStatus.push_back(0);
    mGEMFixAbsGainSlope.push_back(0.0);
    mGEMFixAbsGainConstant.push_back(0.0);
    
    mGEMFixSingleGainFluctuationStatus.push_back(0);
    mGEMFixSingleGainFluctuation.push_back(0);
    
    //Save the ion backflow epsilon contribution for each GEM 
    mGEMEpsilonContribution.push_back(0.0);
    
    //Number of electrons which are collected at the individual GEM stage
    //Note: -1.0 as initial value -> Used in GetFinalCharges() to crosscheck if the energy reolusion (ans thus mN[..]) has already been calculated.
    mN.push_back(-1.0);
    mErrorN.push_back(-1.0);
    
    ++mGEMEntries;
    
    return mGEMEntries-1;
}

void GEMStack::SetGas(int GasType)
{ 
    int ValidGas = 0;
    
    if ( GasType == GAS_NECO2N2 )
    {
	ValidGas = 1;

	mWi = 37.3; //Mean ionization energy / eV [TDR]
	mFano = 0.13; //TODO: Check this!
	mErrorFano = 0.0;
    }
    
    if ( GasType == GAS_NECO2 )
    {
	ValidGas = 1;
	
	mWi = 38.1; //Mean ionization energy / eV [TDR]
	mFano = 0.2; //TODO: Check this!
	mErrorFano = 0.0;
    }

    if ( !ValidGas )
    {
	if ( COUT_WARNINGS ) std::cerr << "Warning: Can not set gas parameters for GasType=" << GasType;
	if ( COUT_WARNINGS ) std::cerr << " -> switch to NeCO2N2" << std::endl;
	
	GasType = GAS_NECO2N2;
    }

    mGasType = GasType;
}

double GEMStack::ConvertGEMTypeToPitch(int GEMType)
{
    double Pitch = 0;
    
    switch ( GEMType )
    {
	case GEM_ST: //standard
	    Pitch = 140.0;
	break;
	case GEM_MP: //medium
	    Pitch = 200.0;
	break;
	case GEM_LP: //large
	    Pitch = 280.0;
	break;
    }
    
    return Pitch;
}

double GEMStack::GetGEMElectronCollectionCore(int GEMId, int Error)
{
    double result;
    
    //Is there a fixed GEM collection for this GEM?
    if ( mGEMFixElectronCollectionStatus[GEMId] == 0 )
    {
	//no, let's calculate it
	double EAbove = mVolumeField[GEMId];
	double EGEM = mGEMVoltage[GEMId]/0.005;
	
	double Eta1Raw = EAbove/EGEM;
	
	//Are we within the fitting rang?
	if ( Eta1Raw < 0.0 || Eta1Raw > 0.16 )
	{
	    if ( COUT_WARNINGS ) std::cerr << "Warning: Electron collection efficiency out of fit range for Eta1=" << Eta1Raw << std::endl;
	}
	
	double Pitch = ConvertGEMTypeToPitch(mGEMType[GEMId]);
	
	GEMEfficiency->SetGeometry(50.0, Pitch, 60.0, 2110.0, 2110.0);
	
	double S1 = GetTuneParameter(0, Pitch);
	
	double Eta1 = S1*Eta1Raw; 
	double Eta2 = 0.0/60000.0;
	
	double c7b = GEMEfficiency->GetParameter_C7Bar(Eta1,Eta2);
	double c8b = GEMEfficiency->GetParameter_C8Bar(Eta1,Eta2);
	double c9b = GEMEfficiency->GetParameter_C9Bar(Eta1,Eta2);
	
	double c1 = GEMEfficiency->GetParameter_C1();
	double c2 = GEMEfficiency->GetParameter_C2();
	double c3 = GEMEfficiency->GetParameter_C3();  

	//Do we want efficiency or error?
	if ( Error == 0 )
	{
	    //Fit model: CollTop 6.0
	    result = 2.0*3.1415927*(c7b+c9b*Eta1+c8b*Eta2)/(c1+c3*Eta1+c2*Eta2); 
	}else{
	  
	    double h = 0.00001;
	  
	    double dS1 = GetErrorTuneParameter(0, Pitch);
	    
	    double dEta1 = Eta1Raw*dS1;
	    
	    double dc7b = (GEMEfficiency->GetParameter_C7Bar(Eta1+h,Eta2)-GEMEfficiency->GetParameter_C7Bar(Eta1,Eta2))/h;
	    double dc8b = (GEMEfficiency->GetParameter_C8Bar(Eta1+h,Eta2)-GEMEfficiency->GetParameter_C8Bar(Eta1,Eta2))/h;
	    double dc9b = (GEMEfficiency->GetParameter_C9Bar(Eta1+h,Eta2)-GEMEfficiency->GetParameter_C9Bar(Eta1,Eta2))/h;
	    
	    double dCollTop1 = 2.0*3.1415927*(dc7b+dc9b*Eta1+c9b+dc8b*Eta2)/(c2*Eta2+c3*Eta1+c1);
	    double dCollTop2 = -2.0*3.1415927*c3*(c7b+c9b*Eta1+c8b*Eta2)/TMath::Power(c2*Eta2+c3*Eta1+c1,2.0);

	    result = TMath::Abs((dCollTop1+dCollTop2)*dEta1);
	}

    }else{
	//Yes. Do we want fixed efficiency or error?
	if ( Error == 0 )
	{
	    //Output fixed efficiency
	    result = mGEMFixElectronCollection[GEMId];
	}else{
	    //Output error of fixed efficiency
	    if ( COUT_MESSAGES ) std::cout << "Trying to compute error for fixed electron collection efficiency -> Set to 0" << std::endl;
	    
	    result = 0.0;
	}
    }
    
    return result;
}

double GEMStack::GetGEMElectronCollection(int GEMId)
{
    return GetGEMElectronCollectionCore(GEMId, 0);
}

double GEMStack::GetErrorGEMElectronCollection(int GEMId)
{
    return GetGEMElectronCollectionCore(GEMId, 1);
}




double GEMStack::GetGEMElectronExtractionCore(int GEMId, int Error)
{
    double result;
    
    //Is there a fixed GEM extraction for this GEM?
    if ( mGEMFixElectronExtractionStatus[GEMId] == 0 )
    {
	//no, let's calculate it
	double EBelow = mVolumeField[GEMId+1];
	double EGEM = mGEMVoltage[GEMId]/0.005;
	
	double Eta2Raw = EBelow/EGEM;
	
	//Are we within the fitting range?
	if ( Eta2Raw < 0.0 || Eta2Raw > 0.16 )
	{
	    if ( COUT_WARNINGS ) std::cerr << "Warning: Electron extraction efficiency out of fit range for Eta2=" << Eta2Raw << std::endl;
	}

	double Pitch = ConvertGEMTypeToPitch(mGEMType[GEMId]);
	
	GEMEfficiency->SetGeometry(50.0, Pitch, 60.0, 2110.0, 2110.0);
	
	double S2 = GetTuneParameter(1, Pitch);
	double S3 = GetTuneParameter(2, Pitch); //Actually S3 is only a constant  
	
	double Eta1 = 2000.0/60000.0; 
	double Eta2 = S2*Eta2Raw;
	
	double c7 = GEMEfficiency->GetParameter_C7(Eta1,Eta2);
	double c8 = GEMEfficiency->GetParameter_C8(Eta1,Eta2);
	double c9 = GEMEfficiency->GetParameter_C9(Eta1,Eta2);
	
	double c4 = GEMEfficiency->GetParameter_C4();
	double c5 = GEMEfficiency->GetParameter_C5();
	double c6 = GEMEfficiency->GetParameter_C6();    
	
	//Do we want efficiency or error?
	if ( Error == 0 )
	{
	    //Fit model: ExtrBot 6.1
	    result = 2.0*3.1415927*(c7+c8*Eta1+c9*Eta2)/(TMath::Power(S3,3.5)*c4+c5*Eta1+1.0/S3*c6*Eta2);
	}else{
	  
	    double h = 0.00001;	    
	    
	    double dS2 = GetErrorTuneParameter(1, Pitch);
	    double dS3 = GetErrorTuneParameter(2, Pitch);
	    
	    double dEta2 = Eta2Raw*dS2;
	    
	    double dc7 = (GEMEfficiency->GetParameter_C7(Eta1,Eta2+h)-GEMEfficiency->GetParameter_C7(Eta1,Eta2))/h;
	    double dc8 = (GEMEfficiency->GetParameter_C8(Eta1,Eta2+h)-GEMEfficiency->GetParameter_C8(Eta1,Eta2))/h;
	    double dc9 = (GEMEfficiency->GetParameter_C9(Eta1,Eta2+h)-GEMEfficiency->GetParameter_C9(Eta1,Eta2))/h;
	    
	    double dExtrBot1 = 2.0*3.1415927*(dc7+dc8*Eta1+dc9*Eta2+c9)/(TMath::Power(S3,3.5)*c4+c5*Eta1+(1.0/S3)*Eta2*c6);
	    double dExtrBot2 = -2.0*3.1415927*c6*(c7+c8*Eta1+c9*Eta2)/(TMath::Power(TMath::Power(S3,3.5)*c4+c5*Eta1+(1.0/S3)*Eta2*c6,2.0)*S3);
	    double dExtrBot3 = -2.0*3.1415927*(c7+c8*Eta1+c9*Eta2)/TMath::Power(TMath::Power(S3,3.5)*c4+c5*Eta1+(1.0/S3)*Eta2*c6,2.0);
	    dExtrBot3 *= ((TMath::Power(S3,3.5)*3.5*c4)/S3-(Eta2*c6)/(S3*S3));

	    result = TMath::Power((dExtrBot1+dExtrBot2)*(dExtrBot1+dExtrBot2)*dEta2*dEta2+dExtrBot3*dExtrBot3*dS3*dS3,0.5);
	}
	  
    }else{
	//Yes. Do we want fixed efficiency or error?
	if ( Error == 0 )
	{
	    //Output fixed efficiency
	    result = mGEMFixElectronExtraction[GEMId];
	}else{
	    //Output error of fixed efficiency
	    if ( COUT_MESSAGES ) std::cout << "Trying to compute error for fixed electron extracition efficiency -> Set to 0" << std::endl;
	    
	    result = 0.0;
	}
    }
    
    return result;
}

double GEMStack::GetGEMElectronExtraction(int GEMId)
{
    return GetGEMElectronExtractionCore(GEMId, 0);
}

double GEMStack::GetErrorGEMElectronExtraction(int GEMId)
{
    return GetGEMElectronExtractionCore(GEMId, 1);
}

//TODO: Unify both following functions
double GEMStack::GetTuneParameter(int TuneParameterID, double Pitch)
{    
    if ( mGasType == GAS_NECO2N2 || mGasType == GAS_NECO2 )
    {
	if ( mTuneParameterMethode == 0 )
	{
	    //Use polynomial fit and relation tuningparameter(pitch)
	    double p0 = TuneParameterCoefficients[mGasType][TuneParameterID][0][0];
	    double p1 = TuneParameterCoefficients[mGasType][TuneParameterID][1][0];
	    double p2 = TuneParameterCoefficients[mGasType][TuneParameterID][2][0];
	    
	    return p0 + p1 * Pitch + p2 * Pitch * Pitch;
	    
	}else{
	  
	    //Use direct fit results for tuning parameters
	    int GEMType;
	    
	    switch ( (int) Pitch )
	    {
		  case 140:
		      GEMType = GEM_ST;
		      break;
		  case 200:
		      GEMType = GEM_MP;
		      break;
		  case 280:
		      GEMType = GEM_LP;
		      break;
	    }
	    
	    return TuneParameter[mGasType][TuneParameterID][GEMType][0];
	}
	
	
	
    }else{
      
	if ( COUT_WARNINGS ) std::cerr << "Warning: Can not calculate tuning parameter s" << TuneParameterID+1 << " for GasType=" << mGasType;
	if ( COUT_WARNINGS ) std::cerr << " -> Set to 1" << std::endl;
	
	return 1.0;
    }
}

double GEMStack::GetErrorTuneParameter(int TuneParameterID, double Pitch)
{
    if ( mGasType == GAS_NECO2N2 || mGasType == GAS_NECO2 )
    {
	if ( mTuneParameterMethode == 0 )
	{
	    //Use polynomial fit and relation tuningparameter(pitch)
	    double p0err = TuneParameterCoefficients[mGasType][TuneParameterID][0][1];
	    double p1err = TuneParameterCoefficients[mGasType][TuneParameterID][1][1];
	    double p2err = TuneParameterCoefficients[mGasType][TuneParameterID][2][1];
	
	    return TMath::Sqrt(p0err*p0err+(p1err*Pitch)*(p1err*Pitch)+(Pitch*Pitch*p2err)*(Pitch*Pitch*p2err));
	}else{
	    //Use direct fit results for tuning parameters
	    int GEMType;
	    
	    switch ( (int) Pitch )
	    {
		  case 140:
		      GEMType = GEM_ST;
		      break;
		  case 200:
		      GEMType = GEM_MP;
		      break;
		  case 280:
		      GEMType = GEM_LP;
		      break;
	    }
	    
	    return TuneParameter[mGasType][TuneParameterID][GEMType][1];	  
	}
	
    }else{
      
	if ( COUT_WARNINGS ) std::cerr << "Warning: Can not calculate error for tuning parameter s" << TuneParameterID+1 << " for GasType=" << mGasType;
	if ( COUT_WARNINGS ) std::cerr << " -> Set to 0" << std::endl;
	
	return 0.0;
    }  
}

double GEMStack::GetGEMAbsGainCore(int GEMId, int Error)
{
    double Slope;
    double Constant;
    
    double ErrorSlope;
    double ErrorConstant;
    
    int GEMType = mGEMType[GEMId];
    double GEMVoltage = mGEMVoltage[GEMId];
    
    //Is the gain curve user-defined for this specific GEM?
    if ( mGEMFixAbsGainStatus[GEMId] == 0 )
    {
	int ValidGas = 0;
	if ( mGasType == GAS_NECO2N2 || mGasType == GAS_NECO2 )
	{
	    Slope = mGainCoefficients[mGasType][GEMType][0][0];
	    Constant = mGainCoefficients[mGasType][GEMType][1][0];
	    
	    ErrorSlope = mGainCoefficients[mGasType][GEMType][0][1];
	    ErrorConstant = mGainCoefficients[mGasType][GEMType][1][1];	    
	    
	    ValidGas = 1;
	    
	    if ( GEMVoltage < mGainCoefficientsRange[mGasType][GEMType][0] || GEMVoltage > mGainCoefficientsRange[mGasType][GEMType][1] )
	    {
		if ( COUT_WARNINGS ) std::cerr << "Warning: Absolute gain out of fit range for GEMVoltage=" << GEMVoltage << "V (GasType=" << mGasType << ", GEMType=" << GEMType << ") -> still using fit" << std::endl;
	    }	
	}
	    
	if ( !ValidGas )
	{
	    if ( COUT_WARNINGS ) std::cerr << "Warning: Can not calculate gain for GasType=" << mGasType;   
	    if ( COUT_WARNINGS ) std::cerr << " -> Gain is fixed to 1 +/- 0" << std::endl;
	    
	    Slope = 0.0;
	    Constant = 0.0;
	    
	    ErrorSlope = 0.0;
	    ErrorConstant = 0.0;
	}
	
    }else{
      
      Slope = mGEMFixAbsGainSlope[GEMId];
      Constant = mGEMFixAbsGainConstant[GEMId];

      ErrorSlope = 0.0;
      ErrorConstant = 0.0;
    }
    
    //Do we want value or error?
    if ( Error == 0 )
    {
	return mGainScaleFactor*TMath::Exp(Constant+Slope*GEMVoltage);
    }else{
      
	if ( mGEMFixAbsGainStatus[GEMId] == 1 && COUT_MESSAGES ) std::cout << "Info: Trying to compute error for a fixed GEM gain curve. Error not supported yet and assumed to be zero!" << std::endl;
      
	double Error1 = mGainScaleFactor*mGainScaleFactor*TMath::Power(TMath::Exp(Constant+Slope*GEMVoltage),2.0)*ErrorConstant*ErrorConstant;
	double Error2 = mGainScaleFactor*mGainScaleFactor*GEMVoltage*GEMVoltage*TMath::Power(TMath::Exp(Constant+Slope*GEMVoltage),2.0)*ErrorSlope*ErrorSlope;
	return TMath::Power(Error1+Error2,0.5);
    }
}

double GEMStack::GetGEMAbsGain(int GEMId)
{
    return GetGEMAbsGainCore(GEMId, 0);
}


double GEMStack::GetErrorGEMAbsGain(int GEMId)
{
    return GetGEMAbsGainCore(GEMId, 1);
}

//Set gain curve for gas and GEM type
void GEMStack::SetAbsGain(int GasType, int GEMType, double Slope, double Constant, double ErrorSlope, double ErrorConstant, double FitUGEMMin, double FitUGEMMax)
{
    mGainCoefficients[GasType][GEMType][0][0] = Slope;
    mGainCoefficients[GasType][GEMType][1][0] = Constant;
    
    mGainCoefficients[GasType][GEMType][0][1] = ErrorSlope;
    mGainCoefficients[GasType][GEMType][1][1] = ErrorConstant;
    
    mGainCoefficientsRange[GasType][GEMType][0] = FitUGEMMin;
    mGainCoefficientsRange[GasType][GEMType][1] = FitUGEMMax;
}

//Set gain curve for a specific GEM within the stack, indep. of the general definition for GEM and gastype
void GEMStack::SetGEMAbsGain(int GEMId, double Slope, double Constant)
{
    mGEMFixAbsGainStatus[GEMId] = 1;
    mGEMFixAbsGainSlope[GEMId] = Slope;
    mGEMFixAbsGainConstant[GEMId] = Constant;
}


void GEMStack::SetSingleGainFluctuationCoefficients(int GasType, int GEMType, int ModelType, double f0, double df0, double U0, double dU0, double Q, double dQ)    
{
    mSingleGainFluctuationCoefficients[GasType][GEMType][ModelType][0][0] = f0;
    mSingleGainFluctuationCoefficients[GasType][GEMType][ModelType][0][1] = df0;
    mSingleGainFluctuationCoefficients[GasType][GEMType][ModelType][1][0] = U0;
    mSingleGainFluctuationCoefficients[GasType][GEMType][ModelType][1][1] = dU0;
    mSingleGainFluctuationCoefficients[GasType][GEMType][ModelType][2][0] = Q;
    mSingleGainFluctuationCoefficients[GasType][GEMType][ModelType][2][1] = dQ;  
}

double GEMStack::GetGEMSingleGainFluctuation(int GEMId)
{
    int Dummy;
    double result = GetGEMSingleGainFluctuationCore(GEMId, &Dummy, 0);
    return result;
}

double GEMStack::GetGEMSingleGainFluctuation(int GEMId, int *OutOfRange)
{
    return GetGEMSingleGainFluctuationCore(GEMId, OutOfRange, 0);
}

double GEMStack::GetErrorHiGEMSingleGainFluctuation(int GEMId)
{
    int Dummy;
    return GetGEMSingleGainFluctuationCore(GEMId, &Dummy, 1);  
}

double GEMStack::GetErrorLoGEMSingleGainFluctuation(int GEMId)
{
    int Dummy;
    return GetGEMSingleGainFluctuationCore(GEMId, &Dummy, 2);  
}

double GEMStack::GetErrorMaxGEMSingleGainFluctuation(int GEMId)
{
    int Dummy;
    return GetGEMSingleGainFluctuationCore(GEMId, &Dummy, 3);    
}

double GEMStack::SingleGainFluctuationDistribution(double f0, double U0, double Q, double GEMVoltage)
{
    if ( mModelType == MODEL_ARC )
    {
	return 0.5*(f0+1.0)+(1.0-f0)/3.1415927*TMath::ATan(-Q*(GEMVoltage-U0));
    }
    
    if ( mModelType == MODEL_EXP )
    {
	double result = f0+TMath::Exp(-(GEMVoltage-U0)*Q);
	
	if ( result > 1.0 ) result = 1.0;
	
	return result; 
    }
    
    if ( mModelType == MODEL_EXP_NOLIMIT )
    {
	return f0+TMath::Exp(-(GEMVoltage-U0)*Q);
    }
    
}


double GEMStack::GetGEMSingleGainFluctuationCore(int GEMId, int *OutOfRange, int Error)
{ 
  
    //Error:
    //0 Value
    //1 Error High
    //2 Error Low
  
    int GEMType = mGEMType[GEMId];
    double GEMVoltage = mGEMVoltage[GEMId];
    
    double result;
    
    double f0, U0, Q;
    double df0, dU0, dQ;
    
    *OutOfRange = 0;
    
    //Fixed value by user?
    if ( mGEMFixSingleGainFluctuationStatus[GEMId] == 0 )
    {
	int ValidGas = 0;
	
	if ( mGasType == GAS_NECO2N2 ||  mGasType == GAS_NECO2 )
	{
	      ValidGas = 1;
	      
	      f0 = mSingleGainFluctuationCoefficients[mGasType][GEMType][mModelType][0][0];
	      U0 = mSingleGainFluctuationCoefficients[mGasType][GEMType][mModelType][1][0];
	      Q  = mSingleGainFluctuationCoefficients[mGasType][GEMType][mModelType][2][0];
	      	      
	      //df0 = mSingleGainFluctuationCoefficients[mGasType][GEMType][mModelType][0][1];
	      dU0 = mSingleGainFluctuationCoefficients[mGasType][GEMType][mModelType][1][1];
	      //dQ  = mSingleGainFluctuationCoefficients[mGasType][GEMType][mModelType][2][1];
	}
	 
	
	if ( ValidGas )
	{
	    if ( Error == 0 )
	    {
		result = SingleGainFluctuationDistribution(f0, U0, Q, GEMVoltage);
	    }
	    
	    //TODO: Implement error propagation for f0 and Q
	    if ( Error == 1 )
	    {
		result = SingleGainFluctuationDistribution(f0, U0+dU0, Q, GEMVoltage)-SingleGainFluctuationDistribution(f0, U0, Q, GEMVoltage);
	    }
	   
	    if ( Error == 2 )
	    {
		result = SingleGainFluctuationDistribution(f0, U0, Q, GEMVoltage)-SingleGainFluctuationDistribution(f0, U0-dU0, Q, GEMVoltage);
	    }
	   
	}else{
	    if ( COUT_WARNINGS ) std::cerr << "Warning: No single gain fluctuation known for GasType=" << mGasType << "! -> Set to 1" << std::endl;
	    result = 1.0;
	}
	
	if ( GEMVoltage <= 200 ) 
	{
	    if ( COUT_WARNINGS ) std::cerr << "Warning: Single gain fluctuation out of fit range for GEMId=" << GEMId << " and voltage " << GEMVoltage << "! -> Still using fit model" << std::endl;
	    
	    *OutOfRange = 1;
	}
	
    }else{
	//Fixed value by user
      
	//Value or error?
	if ( Error == 0 )
	{
	    result = mGEMFixSingleGainFluctuation[GEMId];
	}else{
	  
	    if ( COUT_MESSAGES ) std::cout << "Info: Trying to compute error for fixed single gain fluctuation -> Set to 0" << std::endl;
	    result = 0.0;
	}
    }
    
    return result;
}

//Old model calculation (from wrong gain simulations)
/*
double GEMStack::GetGEMSingleGainFluctuation(int GEMId, int *OutOfRange)
{   
    int GEMType = mGEMType[GEMId];
    double GEMVoltage = mGEMVoltage[GEMId];
    
    double result;
    
    *OutOfRange = 0;

    //Fixed value by user?
    if ( mGEMFixSingleGainFluctuationStatus[GEMId] == 0 )
    {
	//No.. calculate it 
      
	double p0, p1, p2; //Results from GEM_Gain_Analysis
	
	int error = 0;
	int ValidGas = 0;
	
	if ( mGasType == GAS_NECO2N2 )
	{
	    ValidGas = 1;
	    int ValidGEMType = 0;

	    if ( GEMType == GEM_ST )
	    {
		//Data range: 280V - 500V
		p0 = 1.56162;
		p1 = -0.00468792;
		p2 = 4.91418e-06;
		ValidGEMType = 1;

		if ( GEMVoltage >= 180.0 && GEMVoltage < 280.0  )
		{
		    if ( COUT_WARNINGS ) std::cerr << "Warning: Gain fluctuation f out of fit range for GEMVoltage=" << GEMVoltage << "V, standard pitch, NeCO2N2 90-10-5";
		    if ( COUT_WARNINGS ) std::cerr << " -> still using fit" << std::endl; 
		    *OutOfRange = 1;
		}
		
		if ( GEMVoltage < 180.0  )
		{
		    if ( COUT_WARNINGS ) std::cerr << "Warning: Gain fluctuation f out of fit range for GEMVoltage=" << GEMVoltage << "V, standard pitch, NeCO2N2 90-10-5";
		    if ( COUT_WARNINGS ) std::cerr << " -> set to 1.0" << std::endl; 
		    
		    *OutOfRange = 1;
		    error = 1;
		    result = 1.0;
		}
		
		if ( GEMVoltage > 500  )
		{
		    if ( COUT_WARNINGS ) std::cerr << "Warning: Gain fluctuation f out of fit range for GEMVoltage=" << GEMVoltage << "V, standard pitch, NeCO2N2 90-10-5";
		    if ( COUT_WARNINGS ) std::cerr << " -> set to 0.45" << std::endl; 
		    
		    *OutOfRange = 1;
		    error = 1;
		    result = 0.45;
		}

	    }
	    
	    if ( GEMType == GEM_LP )
	    {
		//Data range: 280V - 460V
		p0 = 2.16906;
		p1 = -0.00816123;
		p2 = 9.74361e-06;
		ValidGEMType = 1;
		
		if ( GEMVoltage >= 180.0 && GEMVoltage < 280.0  )
		{
		    if ( COUT_WARNINGS ) std::cerr << "Warning: Gain fluctuation f out of fit range for GEMVoltage=" << GEMVoltage << "V, large pitch, NeCO2N2 90-10-5";
		    if ( COUT_WARNINGS ) std::cerr << " -> still using fit" << std::endl; 
		    *OutOfRange = 1;
		}
		
		if ( GEMVoltage < 180.0  )
		{
		    if ( COUT_WARNINGS )  std::cerr << "Warning: Gain fluctuation f out of fit range for GEMVoltage=" << GEMVoltage << "V, large pitch, NeCO2N2 90-10-5";
		    if ( COUT_WARNINGS )  std::cerr << " -> set to 1.0" << std::endl; 
		    
		    *OutOfRange = 1;
		    error = 1;
		    result = 1.0;
		}
		
		if ( GEMVoltage > 460  )
		{
		    if ( COUT_WARNINGS )  std::cerr << "Warning: Gain fluctuation f out of fit range for GEMVoltage=" << GEMVoltage << "V, large pitch, NeCO2N2 90-10-5";
		    if ( COUT_WARNINGS )  std::cerr << " -> set to 0.45" << std::endl; 
		    
		    *OutOfRange = 1;
		    error = 1;
		    result = 0.45;
		}
		
		
	    }
	    
	    if ( !ValidGEMType )
	    {
		if ( COUT_WARNINGS ) std::cerr << "Warning: Can not calculate single gain fluctuation for NeCO2N2 90-10-5 and GEMType=" << GEMType;  
		if ( COUT_WARNINGS ) std::cerr << " -> set to f=1.0" << std::endl; 	
		
		*OutOfRange = 1;
		error = 1;
		result = 1.0;
	    }
	}
	
	if ( mGasType == GAS_NECO2 )
	{
	    ValidGas = 1;
	    int ValidGEMType = 0;
	    
	    if ( GEMType == GEM_ST )
	    {
		//Data range: 240V - 360V
		p0 = 2.80162;
		p1 = -0.0129278;
		p2 = 1.78423e-05;
		ValidGEMType = 1;

	
		if ( GEMVoltage >= 180 && GEMVoltage < 240.0  )
		{
		    if ( COUT_WARNINGS ) std::cerr << "Warning: Gain fluctuation f out of fit range for GEMVoltage=" << GEMVoltage << "V, standard pitch, NeCO2 90-10";
		    if ( COUT_WARNINGS ) std::cerr << " -> still using fit" << std::endl; 
		    *OutOfRange = 1;
		}
		
		if ( GEMVoltage < 180.0 )
		{
		    if ( COUT_WARNINGS ) std::cerr << "Warning: Gain fluctuation f out of fit range for GEMVoltage=" << GEMVoltage << "V, standard pitch, NeCO2 90-10";
		    if ( COUT_WARNINGS ) std::cerr << " -> set to f=1.0" << std::endl; 	
		    
		    *OutOfRange = 1;
		    error = 1;
		    result = 1.0;
		}
		
		if ( GEMVoltage > 360.0 )
		{
		    if ( COUT_WARNINGS ) std::cerr << "Warning: Gain fluctuation f out of fit range for GEMVoltage=" << GEMVoltage << "V, standard pitch, NeCO2 90-10";
		    if ( COUT_WARNINGS ) std::cerr << " -> set to f=0.46" << std::endl; 	
		    
		    *OutOfRange = 1;
		    error = 1;
		    result = 0.46;
		}
	    }
	    
	    
	    if ( !ValidGEMType )
	    {
		if ( COUT_WARNINGS ) std::cerr << "Warning: Can not calculate single gain fluctuation for NeCO2 90-10 and GEMType=" << GEMType;  
		if ( COUT_WARNINGS ) std::cerr << " -> set to f=1.0" << std::endl; 	
		
		*OutOfRange = 1;
		error = 1;
		result = 1.0;
	    }
	}
	
	if ( !ValidGas  )
	{
	    if ( COUT_WARNINGS ) std::cerr << "WARNING: Can not calculate single gain fluctuation for GasType=" << mGasType;  
	    if ( COUT_WARNINGS ) std::cerr << " -> set to f=1.0" << std::endl; 
	    
	    *OutOfRange = 1;
	    error = 1;
	    result = 1.0;
	}
	
	if ( !error )
	{
	    result = p0 + p1 * GEMVoltage + p2 * GEMVoltage * GEMVoltage;
	}
	
    }else{
	//Fixed value by user
	result = mGEMFixSingleGainFluctuation[GEMId];
    }
    
    if ( result > 1.0 )
      result = 1.0;
    
    return result;
  
}*/

double GEMStack::GetAttachmentFactor(double VolumeField)
{
    //TODO: This here!
    double result = 0.0; 
    
    if ( mAttachmentStatus )
    {
	result = mAttachmentFactor;
    }
    
    return result;
}

//TODO: This here as well
double GEMStack::GetErrorAttachmentFactor(double VolumeField)
{
    return 0.0;
}



double GEMStack::GetEnergyResolutionCore(int Error)
{
    //Calculate number of charges in front of each amplification stage (start counting from 0!!!)
    
    double PrimaryCharges = GetPrimaryCharges();
    //double ErrorPrimaryCharges = GetErrorPrimaryCharges();
    
    if ( COUT_MESSAGES ) std::cout << "====================================" << std::endl;
    if ( COUT_MESSAGES ) std::cout << "Calculating stack energy resolution:" << std::endl;

    for ( int n = 0; n <= mGEMEntries-1; ++n )
    {
	if ( COUT_MESSAGES ) std::cout << "\nGEM ID " << n << ":" << std::endl;
	
	double Attachment = 1.0 - GetAttachmentFactor(mVolumeField[n]) * mVolumeDistance[n]; //1.0 - 1/cm * cm => unitless
	//double ErrorAttachment = GetErrorAttachmentFactor(mVolumeField[n]) * mVolumeDistance[n];
	double ErrorAttachment = 0.0;
	
	if ( COUT_MESSAGES ) std::cout << "Attachment = 1.0 - GetAttachmentFactor(mVolumeField["<<n<<"]) * mVolumeDistance["<<n<<"]" << std::endl;
	if ( COUT_MESSAGES ) std::cout << "           = 1.0 - " << GetAttachmentFactor(mVolumeField[n]) << " * " << mVolumeDistance[n] << " = " << Attachment << " +/- " << ErrorAttachment << std::endl;

	//First GEM stage
	if ( n == 0 )
	{
	    mN[n] = PrimaryCharges*GetGEMElectronCollection(n)*Attachment;  
	    
	    //mErrorN[n] =  TMath::Power(GetGEMElectronCollection(n)*Attachment*ErrorPrimaryCharges,2.0);
	    //mErrorN[n] += TMath::Power(PrimaryCharges*Attachment*GetErrorGEMElectronCollection(n),2.0);
	    //mErrorN[n] += TMath::Power(PrimaryCharges*GetGEMElectronCollection(n)*ErrorAttachment,2.0);
	    //mErrorN[n] = TMath::Power(mErrorN[n],0.5);
	    
	    mErrorN[n] = TMath::Sqrt(mN[n]*(1.0-GetGEMElectronCollection(n)*Attachment));
	    
	    if ( COUT_MESSAGES ) std::cout << "mN["<<n<<"] = PrimaryCharges * GetGEMElectronCollection("<<n<<") * Attachment" << std::endl;
	    if ( COUT_MESSAGES ) std::cout << "mN["<<n<<"] = " << PrimaryCharges << " * " << GetGEMElectronCollection(n) << " * " << Attachment << " = " << mN[n] << " +/- " << mErrorN[n] << std::endl;
	    
	}else{	    
	    mN[n] = mN[n-1]*GetGEMAbsGain(n-1)*GetGEMElectronExtraction(n-1)*GetGEMElectronCollection(n)*Attachment;
	    
	    //mErrorN[n] =  TMath::Power(mErrorN[n-1]*GetGEMAbsGain(n-1)*GetGEMElectronExtraction(n-1)*GetGEMElectronCollection(n)*Attachment,2.0);
	    //mErrorN[n] += TMath::Power(mN[n-1]*GetErrorGEMAbsGain(n-1)*GetGEMElectronExtraction(n-1)*GetGEMElectronCollection(n)*Attachment,2.0);
	    //mErrorN[n] += TMath::Power(mN[n-1]*GetGEMAbsGain(n-1)*GetErrorGEMElectronExtraction(n-1)*GetGEMElectronCollection(n)*Attachment,2.0);
	    //mErrorN[n] += TMath::Power(mN[n-1]*GetGEMAbsGain(n-1)*GetGEMElectronExtraction(n-1)*GetErrorGEMElectronCollection(n)*Attachment,2.0);
	    //mErrorN[n] += TMath::Power(mN[n-1]*GetGEMAbsGain(n-1)*GetGEMElectronExtraction(n-1)*GetGEMElectronCollection(n)*ErrorAttachment,2.0);
	    //mErrorN[n] = TMath::Power(mErrorN[n],0.5);
	    
	    mErrorN[n] = TMath::Sqrt(mN[n]*(1.0-GetGEMElectronExtraction(n-1)*GetGEMElectronCollection(n)*Attachment));
	    
	    if ( COUT_MESSAGES ) std::cout << "mN["<<n<<"] = mN["<<n-1<<"] * GetGEMAbsGain("<<n-1<<") * GetGEMElectronExtraction("<<n-1<<") * GetGEMElectronCollection("<<n<<") * Attachment" << std::endl;
	    if ( COUT_MESSAGES ) std::cout << "mN["<<n<<"] = "<<mN[n-1]<<" * "<<GetGEMAbsGain(n-1)<<" * "<<GetGEMElectronExtraction(n-1)<<" * "<<GetGEMElectronCollection(n)<<" * "<<Attachment<<" = "<<mN[n]<<" +/- "<<mErrorN[n]<<std::endl;
	  
	}

    }
  
    if ( COUT_MESSAGES ) std::cout << "\nAfter GEM ID " << mGEMEntries-1 << ":" <<  std::endl;
  
    //Number of charges after extraction from last amplification stage
    double Attachment = 1.0 - GetAttachmentFactor(mVolumeField[mGEMEntries]) * mVolumeDistance[mGEMEntries];
    //double ErrorAttachment = GetErrorAttachmentFactor(mVolumeField[mGEMEntries]) * mVolumeDistance[mGEMEntries];
    double ErrorAttachment = 0.0;
    
    mN[mGEMEntries] = mN[mGEMEntries-1] * GetGEMAbsGain(mGEMEntries-1) * GetGEMElectronExtraction(mGEMEntries-1) * Attachment;
    
    //mErrorN[mGEMEntries] =  TMath::Power(mErrorN[mGEMEntries-1] * GetGEMAbsGain(mGEMEntries-1) * GetGEMElectronExtraction(mGEMEntries-1) * Attachment,2.0);
    //mErrorN[mGEMEntries] += TMath::Power(mN[mGEMEntries-1] * GetErrorGEMAbsGain(mGEMEntries-1) * GetGEMElectronExtraction(mGEMEntries-1) * Attachment,2.0);
    //mErrorN[mGEMEntries] += TMath::Power(mN[mGEMEntries-1] * GetGEMAbsGain(mGEMEntries-1) * GetErrorGEMElectronExtraction(mGEMEntries-1) * Attachment,2.0);
    //mErrorN[mGEMEntries] += TMath::Power(mN[mGEMEntries-1] * GetGEMAbsGain(mGEMEntries-1) * GetGEMElectronExtraction(mGEMEntries-1) * ErrorAttachment,2.0);
    //mErrorN[mGEMEntries] = TMath::Power(mErrorN[mGEMEntries],0.5);
    
    mErrorN[mGEMEntries] = TMath::Sqrt(mN[mGEMEntries]*(1.0-GetGEMElectronExtraction(mGEMEntries-1) * Attachment));
    
    if ( COUT_MESSAGES ) std::cout << "Attachment = 1.0 - GetAttachmentFactor(mVolumeField["<<mGEMEntries<<"]) * mVolumeDistance["<<mGEMEntries<<"]" << std::endl;
    if ( COUT_MESSAGES ) std::cout << "           = 1.0 - " << GetAttachmentFactor(mVolumeField[mGEMEntries]) << " * " << mVolumeDistance[mGEMEntries] << " = " << Attachment << " +/- " << ErrorAttachment << std::endl;
    
    if ( COUT_MESSAGES ) std::cout << "mN["<<mGEMEntries<<"] = mN["<<mGEMEntries-1<<"] * GetGEMAbsGain("<<mGEMEntries-1<<") * GetGEMElectronExtraction("<<mGEMEntries-1<<") * Attachment" << std::endl;
    if ( COUT_MESSAGES ) std::cout << "                     = " << mN[mGEMEntries-1] << " * " << GetGEMAbsGain(mGEMEntries-1) << " * " << GetGEMElectronExtraction(mGEMEntries-1) << " * " << Attachment << " = " << mN[mGEMEntries] << " +/- " << mErrorN[mGEMEntries] << std::endl;

    //Calculate Energy Resolution
    double SigmaOverMuSquare = 0.0;
    double HiSigmaOverMuSquare = 0.0;
    double LoSigmaOverMuSquare = 0.0;
    
    double SingleGainFluctuation;
    double ErrorHiSingleGainFluctuation;
    double ErrorLoSingleGainFluctuation;
    
    if ( COUT_MESSAGES ) std::cout << "\nSigmaOverMuSquare = " << SigmaOverMuSquare << std::endl;
    
    for ( int n = 0; n <= mGEMEntries-1; ++n )
    {
	SingleGainFluctuation = GetGEMSingleGainFluctuation(n);
	ErrorHiSingleGainFluctuation = GetErrorHiGEMSingleGainFluctuation(n);
	ErrorLoSingleGainFluctuation = GetErrorLoGEMSingleGainFluctuation(n);
	
	//First stage needs Fano Factor
	if ( n == 0 )
	{
	    if ( COUT_MESSAGES ) std::cout << "SigmaOverMuSquare = SigmaOverMuSquare + (mFano + GetGEMSingleGainFluctuation("<<n<<"))/mN["<<n<<"]" << std::endl;
	    if ( COUT_MESSAGES ) std::cout << "                  = "<<SigmaOverMuSquare<<" + ("<<mFano<<" + "<<SingleGainFluctuation<<")/"<<mN[n]<<std::endl;
	    
	    SigmaOverMuSquare += (mFano + SingleGainFluctuation)/mN[n];
	    HiSigmaOverMuSquare += (mFano+mErrorFano + SingleGainFluctuation+ErrorHiSingleGainFluctuation)/(mN[n]-mErrorN[n]);
	    LoSigmaOverMuSquare += (mFano-mErrorFano + SingleGainFluctuation-ErrorLoSingleGainFluctuation)/(mN[n]+mErrorN[n]);
	    
	    if ( COUT_MESSAGES ) std::cout << "                  = "<<SigmaOverMuSquare<< std::endl;
	    
	}else{

	    if ( COUT_MESSAGES ) std::cout << "SigmaOverMuSquare = SigmaOverMuSquare + GetGEMSingleGainFluctuation("<<n<<")/mN["<<n<<"]" << std::endl;
	    if ( COUT_MESSAGES ) std::cout << "                  = "<<SigmaOverMuSquare<<" + "<<SingleGainFluctuation<<"/"<<mN[n]<<std::endl;
	    
	    SigmaOverMuSquare += SingleGainFluctuation/mN[n];
	    HiSigmaOverMuSquare += (SingleGainFluctuation+ErrorHiSingleGainFluctuation)/(mN[n]-mErrorN[n]);
	    LoSigmaOverMuSquare += (SingleGainFluctuation-ErrorLoSingleGainFluctuation)/(mN[n]+mErrorN[n]);
	    
	    if ( COUT_MESSAGES ) std::cout << "                  = "<<SigmaOverMuSquare << std::endl;
	}
    }
    
    if ( COUT_MESSAGES ) std::cout << "\nEnergy resolution = sqrt(SigmaOverMuSquare) = sqrt("<<SigmaOverMuSquare<<") -> "<<100.0*TMath::Sqrt(SigmaOverMuSquare)<<"%"<<std::endl;
    
    if ( COUT_MESSAGES ) std::cout << "====================================" << std::endl;
    
    double Eres = TMath::Sqrt(SigmaOverMuSquare);
    double HiEres = TMath::Sqrt(HiSigmaOverMuSquare);
    double LoEres = TMath::Sqrt(LoSigmaOverMuSquare);
    
    double ErrorHiEres = HiEres-Eres;
    double ErrorLoEres = Eres-LoEres;
    
    if ( Error == 0 ) return Eres;
    
    if ( Error == 1 ) return ErrorHiEres;
    
    if ( Error == 2 ) return ErrorLoEres;
}

double GEMStack::GetEnergyResolution()
{
    return GetEnergyResolutionCore(0);
}

double GEMStack::GetErrorHiEnergyResolution()
{
    return GetEnergyResolutionCore(1);
}

double GEMStack::GetErrorLoEnergyResolution()
{
    return GetEnergyResolutionCore(2);
}


double GEMStack::GetGEMEffGain(int GEMId)
{
    //Remark: There is no difference between attachment below or above the GEM.
    //Indeed this should be taken into account but this means that the attachment couples
    //in quadratically. The total product of all eff. gains for a GEM stack differs than
    //from GetTotalEffGain(). In GetTotalEffGain() the transfers are described correctly
    //what I think ... so where's the rub?
    double Attachment = 1.0 - GetAttachmentFactor(mVolumeField[GEMId]) * mVolumeDistance[GEMId];
    return GetGEMElectronCollection(GEMId)*GetGEMAbsGain(GEMId)*GetGEMElectronExtraction(GEMId)*Attachment;
}

double GEMStack::GetErrorGEMEffGain(int GEMId)
{
    double Collection = GetGEMElectronCollection(GEMId);
    double ErrorCollection = GetErrorGEMElectronCollection(GEMId);
    
    double AbsGain = GetGEMAbsGain(GEMId);
    double ErrorAbsGain = GetErrorGEMAbsGain(GEMId);
    
    double Extraction = GetGEMElectronExtraction(GEMId);
    double ErrorExtraction = GetErrorGEMElectronExtraction(GEMId);

    double Error1 = TMath::Power(AbsGain*Extraction*ErrorCollection,2.0);
    double Error2 = TMath::Power(Collection*Extraction*ErrorAbsGain,2.0);
    double Error3 = TMath::Power(Collection*AbsGain*ErrorExtraction,2.0);
    
    return TMath::Power(Error1+Error2+Error3,0.5);
}

double GEMStack::GetIonBackflowCore(int Error)
{
 
    if ( COUT_MESSAGES ) std::cout << "====================================" << std::endl;
    if ( COUT_MESSAGES ) std::cout << "Calculating stack ion backflow:" << std::endl;
    
    //Epsilon contribution from each stage
    double epsilon[mGEMEntries-1];
    double Errorepsilon[mGEMEntries-1];
    
    //Stack epsilon
    double EpsilonStack = 0.0;
    double ErrorEpsilonStack = 0.0;
    
    double EpsilonPrevStages, ErrorEpsilonPrevStages;
    
    //Stack effective gain
    double EffGainStack = 1.0;
    double ErrorEffGainStack = 0.0;
  
    for ( int n = 0; n <= mGEMEntries-1; ++n )
    {
	if ( COUT_MESSAGES ) std::cout << "\nGEM ID " << n << ":" << std::endl;
	
	double EpsilonThisStage = GetGEMElectronCollection(n)*(GetGEMAbsGain(n)-1.0)*GetGEMAvalancheIonCollection(n)*GetGEMAvalancheIonExtraction(n);
	double ErrorEpsilonThisStage;
	ErrorEpsilonThisStage  = TMath::Power((GetGEMAbsGain(n)-1.0)*GetGEMAvalancheIonCollection(n)*GetGEMAvalancheIonExtraction(n)*GetErrorGEMElectronCollection(n),2.0);
	ErrorEpsilonThisStage += TMath::Power(GetGEMElectronCollection(n)*GetGEMAvalancheIonCollection(n)*GetGEMAvalancheIonExtraction(n)*GetErrorGEMAbsGain(n),2.0);
	ErrorEpsilonThisStage += TMath::Power(GetGEMElectronCollection(n)*(GetGEMAbsGain(n)-1.0)*GetGEMAvalancheIonExtraction(n)*GetErrorGEMAvalancheIonCollection(n),2.0);
	ErrorEpsilonThisStage += TMath::Power(GetGEMElectronCollection(n)*(GetGEMAbsGain(n)-1.0)*GetGEMAvalancheIonCollection(n)*GetErrorGEMAvalancheIonExtraction(n),2.0);
	ErrorEpsilonThisStage = TMath::Power(ErrorEpsilonThisStage,0.5);
	
	if ( COUT_MESSAGES ) std::cout << "EpsilonThisStage = GetGEMElectronCollection("<<n<<")*(GetGEMAbsGain("<<n<<")-1.0)*GetGEMAvalancheIonCollection("<<n<<")*GetGEMAvalancheIonExtraction("<<n<<")" << std::endl;
	if ( COUT_MESSAGES ) std::cout << "EpsilonThisStage = " << GetGEMElectronCollection(n) << "*("<<GetGEMAbsGain(n)<<"-1.0)*"<<GetGEMAvalancheIonCollection(n)<<"*"<<GetGEMAvalancheIonExtraction(n)<<"="<<EpsilonThisStage<<std::endl;
	
	//stages before
	EpsilonPrevStages = 1.0;
	ErrorEpsilonPrevStages = 0.0;
	
	if ( COUT_MESSAGES ) std::cout << "\nEpsilonPrevStages = " << EpsilonPrevStages << std::endl;
	
	for ( int nPrev = 0; nPrev <= n-1; ++nPrev )
	{
	    if ( COUT_MESSAGES ) std::cout << "EpsilonPrevStages = EpsilonPrevStages*GetGEMEffGain("<<nPrev<<")*GetGEMDriftIonCollection("<<nPrev<<")*GetGEMDriftIonExtraction("<<nPrev<<")" << std::endl;
	    if ( COUT_MESSAGES ) std::cout << "                  = "<<EpsilonPrevStages<<"*"<<GetGEMEffGain(nPrev)<<"*"<<GetGEMDriftIonCollection(nPrev)<<"*"<<GetGEMDriftIonExtraction(nPrev)<<std::endl;

	    EpsilonPrevStages *= GetGEMEffGain(nPrev)*GetGEMDriftIonCollection(nPrev)*GetGEMDriftIonExtraction(nPrev);
	    
	    if ( COUT_MESSAGES ) std::cout << "                  = "<<EpsilonPrevStages<<std::endl;
	  
	}
	
	for ( int nPrev = 0; nPrev <= n-1; ++nPrev )
	{
	      ErrorEpsilonPrevStages += TMath::Power(EpsilonPrevStages/GetGEMEffGain(nPrev) * GetErrorGEMEffGain(nPrev),2.0);
	      ErrorEpsilonPrevStages += TMath::Power(EpsilonPrevStages/GetGEMDriftIonCollection(nPrev) * GetErrorGEMDriftIonCollection(nPrev),2.0);
	      ErrorEpsilonPrevStages += TMath::Power(EpsilonPrevStages/GetGEMDriftIonExtraction(nPrev) * GetErrorGEMDriftIonExtraction(nPrev),2.0);
	}
	
	ErrorEpsilonPrevStages = TMath::Power(ErrorEpsilonPrevStages,0.5);
	
	epsilon[n] = EpsilonThisStage*EpsilonPrevStages;
	Errorepsilon[n] = TMath::Power(TMath::Power(EpsilonPrevStages*ErrorEpsilonThisStage,2.0)+TMath::Power(EpsilonThisStage*ErrorEpsilonPrevStages,2.0),0.5);
	
	if ( COUT_MESSAGES ) std::cout << "\nepsilon["<<n<<"] = EpsilonThisStage*EpsilonPrevStages = " << epsilon[n] << std::endl;

	if ( COUT_MESSAGES ) std::cout << "EffGainStack = EffGainStack*GetGEMEffGain("<<n<<")" << std::endl;
	if ( COUT_MESSAGES ) std::cout << "             = "<<EffGainStack<<"*"<<GetGEMEffGain(n)<< std::endl;
	
	EffGainStack *= GetGEMEffGain(n); //TODO: error here
	
	if ( COUT_MESSAGES ) std::cout << "             = " << EffGainStack << std::endl;
	
	
	
	if ( COUT_MESSAGES ) std::cout << "EpsilonStack = EpsilonStack + epsilon["<<n<<"]" << std::endl;
	if ( COUT_MESSAGES ) std::cout << "             = " << EpsilonStack << " + " << epsilon[n] << std::endl;
	
	mGEMEpsilonContribution[n] = epsilon[n]; //Save this for upcoming analysis
	EpsilonStack += epsilon[n];
	
	ErrorEpsilonStack += Errorepsilon[n]*Errorepsilon[n];
	
	if ( COUT_MESSAGES ) std::cout << "             = " << EpsilonStack << std::endl;
    }
    
    for ( int n = 0; n <= mGEMEntries-1; ++n )
    {
	ErrorEffGainStack += TMath::Power(EffGainStack/GetGEMEffGain(n)*GetErrorGEMEffGain(n),2.0);
    }
    
    if ( COUT_MESSAGES ) std::cout << "\nIon backflow = (1.0 + EpsilonStack)/EffGainStack = (1.0 + "<<EpsilonStack<<")/"<<EffGainStack<<" -> "<<(1.0 + EpsilonStack)/EffGainStack*100.0<<"%"<<std::endl;
    if ( COUT_MESSAGES ) std::cout << "====================================" << std::endl;
    
    ErrorEpsilonStack = TMath::Power(ErrorEpsilonStack,0.5);
    ErrorEffGainStack = TMath::Power(ErrorEffGainStack,0.5);
    
    double result;
    
    if ( Error )
    {
	//ErrorEpsilonStack = 0.0;
	//ErrorEffGainStack = 0.0;
	result = TMath::Power(ErrorEpsilonStack/EffGainStack,2.0)+TMath::Power(((1.0+EpsilonStack)*ErrorEffGainStack)/(EffGainStack*EffGainStack),2.0);
	result = TMath::Power(result,0.5);
    }else{
	result = (1.0 + EpsilonStack)/EffGainStack;
    }
    return result;
    
}

double GEMStack::GetIonBackflow()
{
    return GetIonBackflowCore(0);
}

double GEMStack::GetErrorIonBackflow()
{
    return GetIonBackflowCore(1);
}


double GEMStack::GetFinalCharges()
{
    if ( mN[mGEMEntries] == -1.0 )
	GetEnergyResolution();
    
    return mN[mGEMEntries];
}

double GEMStack::GetErrorFinalCharges()
{
    if ( mErrorN[mGEMEntries] == -1.0 )
	GetEnergyResolution();
    
    return mErrorN[mGEMEntries];
}

double GEMStack::GetTotalEffGain()
{
    return GetFinalCharges()/GetPrimaryCharges();
}

double GEMStack::GetErrorTotalEffGain()
{
    double N0 = GetPrimaryCharges();
    double ErrorN0 = GetErrorPrimaryCharges();
    
    double N = GetFinalCharges();
    double ErrorN = GetErrorFinalCharges();
    
    //G=N/N0
    return TMath::Sqrt(TMath::Power((1.0/N0)*ErrorN,2.0)+TMath::Power((N*ErrorN0)/(N0*N0),2.0));
}

double GEMStack::GetGEMIonEfficiency(int GEMId, int ChargeType, int EffType, int CalcError)
{  
    double result = 0;

    //Is there a fixed value for this GEMid?
    int FixedEff = 0;
    double FixedEffValue;
    
    if ( EffType == EFF_COLL )
    {
	if ( ChargeType == CHRG_ION_AVAL )
	{
	    if ( mGEMFixAvalancheIonCollectionStatus[GEMId] )
	    {
		FixedEff = 1;
		if ( !CalcError ) FixedEffValue = mGEMFixAvalancheIonCollection[GEMId];
		if ( CalcError ) FixedEffValue = 0.0; //Fixed values do not have an error
	    }
	}
	
	if ( ChargeType == CHRG_ION_DRIFT )
	{
	    if ( mGEMFixDriftIonCollectionStatus[GEMId] )
	    {
		FixedEff = 1;
		if ( !CalcError ) FixedEffValue = mGEMFixDriftIonCollection[GEMId];
		if ( CalcError ) FixedEffValue = 0.0; //Fixed values do not have an error
	    }
	}	
	
    }else{
	if ( ChargeType == CHRG_ION_AVAL )
	{
	    if ( mGEMFixAvalancheIonExtractionStatus[GEMId] )
	    {
		FixedEff = 1;
		if ( !CalcError ) FixedEffValue = mGEMFixAvalancheIonExtraction[GEMId];
		if ( CalcError ) FixedEffValue = 0.0; //Fixed values do not have an error
	    }
	}    
	
	if ( ChargeType == CHRG_ION_DRIFT )
	{
	    if ( mGEMFixDriftIonExtractionStatus[GEMId] )
	    {
		FixedEff = 1;
		if ( !CalcError ) FixedEffValue = mGEMFixDriftIonExtraction[GEMId];
		if ( CalcError ) FixedEffValue = 0.0; //Fixed values do not have an error
	    }
	} 
    }
 
    if ( !FixedEff )
    {
	int ValidGas = 0;
	int ValidGEMType1 = 0;
	int ValidGEMType2 = 0;
	
	double p[8], dp[8];
	
	double EGEM = mGEMVoltage[GEMId]/0.005; 
	double EExtern;
      
	if ( EffType == EFF_COLL )
	{
	    EExtern = mVolumeField[GEMId+1]; //Below
	}else{
	    EExtern = mVolumeField[GEMId]; //Above
	}
	
	double Eta = EExtern/EGEM;
	
	int GEMType1 = mGEMType[GEMId];
	int GEMType2;

	//In case of avalanche ions we are not interested in GEMType2, so we use GEMType "No GEM"
	if ( ChargeType == CHRG_ION_AVAL )
	{
	    GEMType2 = GEM_NO;
	}
	
	//In case oft drift ions ...
	if ( ChargeType == CHRG_ION_DRIFT )
	{
	    //.. what is the next GEM?
	    if ( GEMId <= mGEMEntries-2 )
	    {
		GEMType2 = mGEMType[GEMId+1]; 
	    }else{
		if ( COUT_WARNINGS ) std::cerr << "WARNING: Trying to get drift ion efficiency for GEMId=" << GEMId << " but there is no following GEM in stack";
		if ( COUT_WARNINGS ) std::cerr << " -> Next \"GEM\" set to standard" << std::endl;  
		
		GEMType2 = 0;
	    }
	}
	
	//Do we have the data?
	if ( mGasType == GAS_NECO2  )
	{
	    ValidGas = 1;

	    if ( ChargeType == CHRG_ION_AVAL )
	    {
		if ( GEMType1 == GEM_ST )
		{
		    ValidGEMType1 = 1;
		}
		
		if ( GEMType1 == GEM_MP )
		{
		    ValidGEMType1 = 1;
		}
	    }
	    
	    if ( ChargeType == CHRG_ION_DRIFT )
	    {
		if ( GEMType1 == GEM_ST && GEMType2 == GEM_ST )
		{
		    ValidGEMType1 = 1;
		    ValidGEMType2 = 1;
		}
	    }
	}
	
	if ( mGasType == GAS_NECO2N2  )
	{
	    ValidGas = 1;

	    if ( ChargeType == CHRG_ION_AVAL )
	    {
		if ( GEMType1 == GEM_ST || GEMType1 == GEM_MP || GEMType1 == GEM_LP )
		{
		    ValidGEMType1 = 1;
		}
	    }
	    
	    if ( ChargeType == CHRG_ION_DRIFT )
	    {
		if ( GEMType1 == GEM_ST && GEMType2 == GEM_ST )
		{
		    ValidGEMType1 = 1;
		    ValidGEMType2 = 1;
		}
		
		if ( GEMType1 == GEM_ST && GEMType2 == GEM_LP )
		{
		    ValidGEMType1 = 1;
		    ValidGEMType2 = 1;
		}
		
		if ( GEMType1 == GEM_LP && GEMType2 == GEM_ST )
		{
		    ValidGEMType1 = 1;
		    ValidGEMType2 = 1;
		}
		
		if ( GEMType1 == GEM_LP && GEMType2 == GEM_LP )
		{
		    ValidGEMType1 = 1;
		    ValidGEMType2 = 1;
		}
	    }
	}
	
	int Error = 0;
	
	
	//Error output
	if ( !ValidGas  )
	{
	    if ( COUT_WARNINGS ) std::cerr << "WARNING: No ion efficiency data for GasType=" << mGasType;
	    if ( COUT_WARNINGS ) std::cerr << " -> set to 0!" << std::endl;
	    
	    result = 0.0;
	    Error = 1;
	}
	
	if ( ValidGas && ChargeType == CHRG_ION_AVAL && !ValidGEMType1 )
	{
	    if ( COUT_WARNINGS ) std::cerr << "WARNING: No avalanche ion efficiency data for GasType=" << mGasType <<  " and GEMType1=" << GEMType1;
	    if ( COUT_WARNINGS ) std::cerr << " -> set to 0!" << std::endl;
	    
	    
	    result = 0.0;
	    Error = 1;
	}
	
	if ( ValidGas && ChargeType == CHRG_ION_DRIFT && (!ValidGEMType1 || !ValidGEMType2) )
	{
	    if ( COUT_WARNINGS ) std::cerr << "WARNING: No drift ion efficiency data for GasType=" << mGasType <<  ", GEMType1=" << GEMType1 << " and GEMType2=" << GEMType2;
	    if ( COUT_WARNINGS ) std::cerr << " -> set to 0!" << std::endl;
	 	    
	    result = 0.0;
	    Error = 1;
	}
	
	//Calculate if everything is fine
	if ( !Error )
	{
	    result = 0.0;
	    
	    if ( CalcError )
	    {
		for ( int i = 0; i <= 7; ++i )
		{
		    result += TMath::Power(mGEMIonEfficiencyCoefficients[mGasType][GEMType1][GEMType2][EffType][ChargeType][i][1]*TMath::Power(Eta, i),2.0);
		}
		result = TMath::Power(result,0.5);
		
	    }else{
		for ( int i = 0; i <= 7; ++i )
		{
		    result += mGEMIonEfficiencyCoefficients[mGasType][GEMType1][GEMType2][EffType][ChargeType][i][0]*TMath::Power(Eta, i);	
		}	 
		
	    }

	
	    //Are we within the fitting range for the polynomial fit?
	    double FitEtaMin = mGEMIonEfficiencyCoefficientsRange[mGasType][GEMType1][GEMType2][EffType][ChargeType][0];
	    double FitEtaMax = mGEMIonEfficiencyCoefficientsRange[mGasType][GEMType1][GEMType2][EffType][ChargeType][1];
		
	    if ( Eta < FitEtaMin || FitEtaMax < Eta )
	    {
		//We are out of the fit range!
		if ( ChargeType == CHRG_ION_DRIFT || (ChargeType == CHRG_ION_AVAL && EffType == EFF_COLL) )
		{
		    //For drift ions or avalanche ion collection: Force the eff. to be 1
		    if ( COUT_WARNINGS ) std::cerr << "Warning: Out of ion fit range for Eta=" << Eta << " EffType=" << EffType << " ChargeType=" << ChargeType << " GEMId=" << GEMId << " GasType=" << mGasType;
		    if ( COUT_WARNINGS ) std::cerr << " -> Set to 1" << std::endl;
		    
		    if ( !CalcError ) result = 1.0;
		    
		    if ( CalcError ) result = 0.0; //No error
		}else{
		    //For avalanche ions extraction: Take the maximum value for the largest valid eta and do a linear regression to eta=0.25 where the eff. is forced to be 1
		    double MaxIonEfficiency = 0;
		    
		    for ( int i = 0; i <= 7; ++i )
		    {
			MaxIonEfficiency += mGEMIonEfficiencyCoefficients[mGasType][GEMType1][GEMType2][EffType][ChargeType][i][0]*TMath::Power(FitEtaMax, i);	
		    }		 

		    double LinearSlope = (1.0 - MaxIonEfficiency)/(0.25 - FitEtaMax);
		    
		    if ( COUT_WARNINGS ) std::cerr << "Warning: Out of ion fit range for Eta=" << Eta << " EffType=" << EffType << " ChargeType=" << ChargeType << " GEMId=" << GEMId << " GasType=" << mGasType;
		    if ( COUT_WARNINGS ) std::cerr << " -> Using linear regression to calculate efficiency" << std::endl;
		    
		    if ( !CalcError ) result = MaxIonEfficiency+LinearSlope*(Eta-FitEtaMax);
		    
		    if ( CalcError ) result = 0.0; //No error
		}
	    }
	    
	    //Crosscheck: Forbid values below 0 and above 1
	    if ( result < 0.0 && !CalcError )
	    {
		result = 0.0;
	    }
		    
	    if ( result > 1.0 && !CalcError )
	    {
		result = 1.0;
	    }
	}
    }else{
	result = FixedEffValue;
    }
    
    if ( ChargeType == CHRG_ION_AVAL  )
    {
	return mAvalancheIonEfficienciesScalingFactor * result; //also correct for error calculation
    }else{
	if ( EffType == EFF_COLL )
	{
	    return mDriftIonCollectionScalingFactor * result; //also correct for error calculation
	}else{
	    return mDriftIonExtractionScalingFactor * result;
	}
    }
}

double GEMStack::GetGEMAvalancheIonCollection(int GEMId)
{
    return GetGEMIonEfficiency(GEMId, CHRG_ION_AVAL, EFF_COLL, 0);
}

double GEMStack::GetGEMAvalancheIonExtraction(int GEMId)
{
    return GetGEMIonEfficiency(GEMId, CHRG_ION_AVAL, EFF_EXTR, 0);
}

double GEMStack::GetGEMDriftIonCollection(int GEMId)
{
    if ( mDriftIonModel == MODEL_DRIFT_ION_SIMULATION )
    {
	return GetGEMIonEfficiency(GEMId, CHRG_ION_DRIFT, EFF_COLL, 0);
    }else{
	return GetGEMIonEfficiency(GEMId, CHRG_ION_DRIFT, EFF_COLL, 0);
    }
}

double GEMStack::GetGEMDriftIonExtraction(int GEMId)
{
    if ( mDriftIonModel == MODEL_DRIFT_ION_SIMULATION )
    {
	return GetGEMIonEfficiency(GEMId, CHRG_ION_DRIFT, EFF_EXTR, 0);
    }else{
	return 1.0; //Drift ions are always extracted
    }
}

double GEMStack::GetErrorGEMAvalancheIonCollection(int GEMId)
{
    return GetGEMIonEfficiency(GEMId, CHRG_ION_AVAL, EFF_COLL, 1);
}

double GEMStack::GetErrorGEMAvalancheIonExtraction(int GEMId)
{
    return GetGEMIonEfficiency(GEMId, CHRG_ION_AVAL, EFF_EXTR, 1);
}

double GEMStack::GetErrorGEMDriftIonCollection(int GEMId)
{
    return GetGEMIonEfficiency(GEMId, CHRG_ION_DRIFT, EFF_COLL, 1);
}

double GEMStack::GetErrorGEMDriftIonExtraction(int GEMId)
{
    return GetGEMIonEfficiency(GEMId, CHRG_ION_DRIFT, EFF_EXTR, 1);
}


void GEMStack::SetAvalancheIonCollection(int GasType, int GEMType1, double p0, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double FitEta2Min, double FitEta2Max)
{
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][0][0] = p0;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][1][0] = p1;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][2][0] = p2;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][3][0] = p3;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][4][0] = p4;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][5][0] = p5;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][6][0] = p6;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][7][0] = p7;
    
    mGEMIonEfficiencyCoefficientsRange[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][0] = FitEta2Min;
    mGEMIonEfficiencyCoefficientsRange[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][1] = FitEta2Max;
}

void GEMStack::SetErrorAvalancheIonCollection(int GasType, int GEMType1, double dp0, double dp1, double dp2, double dp3, double dp4, double dp5, double dp6, double dp7)
{
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][0][1] = dp0;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][1][1] = dp1;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][2][1] = dp2;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][3][1] = dp3;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][4][1] = dp4;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][5][1] = dp5;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][6][1] = dp6;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_COLL][CHRG_ION_AVAL][7][1] = dp7;
}

void GEMStack::SetAvalancheIonExtraction(int GasType, int GEMType1, double p0, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double FitEta1Min, double FitEta1Max)
{
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][0][0] = p0;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][1][0] = p1;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][2][0] = p2;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][3][0] = p3;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][4][0] = p4;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][5][0] = p5;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][6][0] = p6;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][7][0] = p7;
   
    mGEMIonEfficiencyCoefficientsRange[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][0] = FitEta1Min;
    mGEMIonEfficiencyCoefficientsRange[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][1] = FitEta1Max;
}

void GEMStack::SetErrorAvalancheIonExtraction(int GasType, int GEMType1, double dp0, double dp1, double dp2, double dp3, double dp4, double dp5, double dp6, double dp7)
{
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][0][1] = dp0;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][1][1] = dp1;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][2][1] = dp2;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][3][1] = dp3;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][4][1] = dp4;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][5][1] = dp5;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][6][1] = dp6;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEM_NO][EFF_EXTR][CHRG_ION_AVAL][7][1] = dp7;
}


void GEMStack::SetDriftIonCollection(int GasType, int GEMType1, int GEMType2, double p0, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double FitEta2Min, double FitEta2Max)
{
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][0][0] = p0;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][1][0] = p1;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][2][0] = p2;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][3][0] = p3;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][4][0] = p4;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][5][0] = p5;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][6][0] = p6;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][7][0] = p7;
    
    mGEMIonEfficiencyCoefficientsRange[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][0] = FitEta2Min;
    mGEMIonEfficiencyCoefficientsRange[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][1] = FitEta2Max;
}

void GEMStack::SetErrorDriftIonCollection(int GasType, int GEMType1, int GEMType2, double dp0, double dp1, double dp2, double dp3, double dp4, double dp5, double dp6, double dp7)
{
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][0][1] = dp0;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][1][1] = dp1;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][2][1] = dp2;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][3][1] = dp3;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][4][1] = dp4;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][5][1] = dp5;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][6][1] = dp6;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_COLL][CHRG_ION_DRIFT][7][1] = dp7;
}

void GEMStack::SetDriftIonExtraction(int GasType, int GEMType1, int GEMType2, double p0, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double FitEta1Min, double FitEta1Max)
{
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][0][0] = p0;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][1][0] = p1;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][2][0] = p2;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][3][0] = p3;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][4][0] = p4;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][5][0] = p5;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][6][0] = p6;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][7][0] = p7;
   
    mGEMIonEfficiencyCoefficientsRange[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][0] = FitEta1Min;
    mGEMIonEfficiencyCoefficientsRange[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][1] = FitEta1Max;
}

void GEMStack::SetErrorDriftIonExtraction(int GasType, int GEMType1, int GEMType2, double dp0, double dp1, double dp2, double dp3, double dp4, double dp5, double dp6, double dp7)
{
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][0][1] = dp0;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][1][1] = dp1;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][2][1] = dp2;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][3][1] = dp3;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][4][1] = dp4;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][5][1] = dp5;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][6][1] = dp6;
    mGEMIonEfficiencyCoefficients[GasType][GEMType1][GEMType2][EFF_EXTR][CHRG_ION_DRIFT][7][1] = dp7;
   
}

