/* ===============================================================*/
/* GEM stack calculations                                         */
/* Example Root macro for some basic GEM stack operations.        */
/* Simply run by: root Example.c                                  */
/*                                                                */
/* Written by Viktor Ratza, ratza@hiskp.uni-bonn.de, May 2019     */
/* ===============================================================*/

#include "include/GEMStack.h"


void DrawEfficiencies()
{
    GEMStack *Detector = new GEMStack();
    
    Detector->SetGas(GAS_NECO2N2);
    Detector->SetPhotonEnergy(5900); 
    Detector->SetAttachment(0); 
    Detector->SetAttachmentFactor(1.0);
    Detector->SetDriftIonCollectionScalingFactor(1.0);
    Detector->SetTuneParameterMethode(1); //0: use relation parameter(pitch), 1: direct fit result

    int Vol1 = Detector->AddVolume(4.0, 400.0);
    int GEM1 = Detector->AddGEM(GEM_ST, 300.0);
    int Vol2 = Detector->AddVolume(4.0, 2000.0);
    int GEM2 = Detector->AddGEM(GEM_ST, 350.0);
    int Vol3 = Detector->AddVolume(4.0, 2000.0);
   
    TCanvas *c = new TCanvas("c","c",0,0,1024,768);
    c->SetGrid(); 

    Detector->DrawGEMEfficiencyCurves(GEM1, c, 1, 1, 1);
}

void DrawGaincurves()
{
    GEMStack *Detector = new GEMStack();
    
    Detector->SetGas(GAS_NECO2N2);
    Detector->SetPhotonEnergy(5900); 
    Detector->SetAttachment(0); 

    //How should eff. be calculated?
    //0: Use polynomial fit tuningparameter(pitch), this is bad for NeCO2 as we need to simulate more pitches in order to get a better polynomial description
    //1: Use direct fit results for s1, s2, and s3
    Detector->SetTuneParameterMethode(0); 
    
    int EAbove = Detector->AddVolume(4.0, 400.0);
    int GEM = Detector->AddGEM(GEM_LP, 0.0);
    int EBelow = Detector->AddVolume(4.0, 2000.0);    
    
    std::vector<double> AbsGain[2], EffGain[2], vecUGEM; //0 NeCO2N2 1 NeCO2
    std::vector<double> ErrorAbsGain[2], ErrorEffGain[2];
    int GainCount = 0;

    
    double UGEM = 200.0; //Volts
    
    while ( UGEM <= 400 )
    {
	Detector->SetGEMVoltage(GEM, UGEM);
	
	Detector->SetGas(GAS_NECO2N2);
	AbsGain[0].push_back(Detector->GetGEMAbsGain(GEM));
	EffGain[0].push_back(Detector->GetGEMEffGain(GEM));
	ErrorAbsGain[0].push_back(Detector->GetErrorGEMAbsGain(GEM));
	ErrorEffGain[0].push_back(Detector->GetErrorGEMEffGain(GEM));
	
	Detector->SetGas(GAS_NECO2);
	AbsGain[1].push_back(Detector->GetGEMAbsGain(GEM));
	EffGain[1].push_back(Detector->GetGEMEffGain(GEM));
	ErrorAbsGain[1].push_back(Detector->GetErrorGEMAbsGain(GEM));
	ErrorEffGain[1].push_back(Detector->GetErrorGEMEffGain(GEM));
	
	vecUGEM.push_back(UGEM);
	++GainCount;
	UGEM += 20.0;
    }
    
    TCanvas *cGain = new TCanvas("cGain","cGain",0,0,1024,768);
    cGain->SetGrid(); 
    
    gPad->SetLogy();
    
    TH1F *hrGain = cGain->DrawFrame(140.0,0.1,500.0,10000.0);    
    hrGain->SetXTitle("U_{GEM}");
    hrGain->SetYTitle("Gain"); 
    
    TGraph *grAbsGain[2];
    TGraph *grEffGain[2];
    
    for ( int i = 0; i <= 1; ++i )
    {
	grAbsGain[i] = new TGraphErrors(GainCount, &vecUGEM[0], &AbsGain[i][0], NULL, &ErrorAbsGain[i][0]);
	grAbsGain[i]->SetLineWidth(3);
	grAbsGain[i]->SetLineColor(i+1);
	grAbsGain[i]->Draw("le");
	
	grEffGain[i] = new TGraphErrors(GainCount, &vecUGEM[0], &EffGain[i][0], NULL, &ErrorEffGain[i][0]);
	grEffGain[i]->SetLineStyle(2);
	grEffGain[i]->SetLineWidth(3);
	grEffGain[i]->SetLineColor(i+1);
	grEffGain[i]->Draw("le");
    }
    
    TLegend *legGain = new TLegend(0.40,0.15,0.40+0.35,0.15+0.20);
    legGain->SetTextSize(0.02);
    legGain->AddEntry(grAbsGain[0],"Abs. gain NeCO2N2", "l");
    legGain->AddEntry(grAbsGain[1],"Abs. gain NeCO2", "l");
    legGain->AddEntry(grEffGain[0],"Eff. gain NeCO2N2", "l");
    legGain->AddEntry(grEffGain[1],"Eff. gain NeCO2", "l");
    legGain->Draw();

    
}

void SimpleStack()
{
    GEMStack *Detector = new GEMStack();
    
    Detector->SetGas(GAS_NECO2N2);
    Detector->SetPhotonEnergy(5900); 
    
    Detector->SetAttachment(1); 
    Detector->SetAttachmentFactor(0.06);
    Detector->SetDriftIonCollectionScalingFactor(0.30);
    Detector->SetAbsGainScalingFactor(1.51);
    Detector->SetSingleGainFluctuationModel(MODEL_EXP_NOLIMIT); 
   
    Detector->AddVolume(4.0, 400.0);
    Detector->AddGEM(GEM_ST, 265.0);
    Detector->AddVolume(0.2, 4000.0);
    Detector->AddGEM(GEM_LP, 235.0);
    Detector->AddVolume(0.2, 4000.0);
    Detector->AddGEM(GEM_LP, 287.0);
    Detector->AddVolume(0.2, 100.0);
    Detector->AddGEM(GEM_ST, 360.0);
    Detector->AddVolume(0.2, 4000.0);
    
    std::cout << "Resolution: " << Detector->GetEnergyResolution() << std::endl;
    std::cout << "Ion backflow: " << Detector->GetIonBackflow() << std::endl;
    std::cout << "Gain: " << Detector->GetFinalCharges()/Detector->GetPrimaryCharges() << std::endl;
   
}


void Example()
{
    DrawEfficiencies();
    
    DrawGaincurves();

    SimpleStack();

}
