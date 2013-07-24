/////////////////////////////////////////////////////////////////////////
//
// SPlot tutorial
// author: Kyle Cranmer
// date Dec. 2008
//
// This tutorial shows an example of using SPlot to unfold two distributions.
// The physics context for the example is that we want to know
// the isolation distribution for real electrons from Z events
// and fake electrons from QCD.  Isolation is our 'control' variable
// To unfold them, we need a model for an uncorrelated variable that
// discriminates between Z and QCD.  To do this, we use the invariant
// mass of two electrons.  We model the Z with a Gaussian and the QCD
// with a falling exponential.
//
// Note, since we don't have real data in this tutorial, we need to generate
// toy data.  To do that we need a model for the isolation variable for
// both Z and QCD.  This is only used to generate the toy data, and would
// not be needed if we had real data.
/////////////////////////////////////////////////////////////////////////


#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

//#define isMC 1


/*#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
*/


// use this order for safety on library loading
using namespace RooFit ;
using namespace RooStats ;



TFile *myFileMC = new TFile("histoFiles_mc.root"); // the Z resonnances in MC
//TFile *myFileDATA = new TFile("histoFiles_data.root"); // the Z resonnaces in DATA (for tests purpose)
TFile *saveHisto; // the root files to store the result
TChain *chain = new TChain("treeClusterShape"); //the TCHAIN with the tree

////////////
int nbSpectator = 33;

TString theVars[33] ={//name of the tree branches were are store the vars of interest
    "CL_R9","CL_fBrem","CL_hkfchi2","CL_hkfhits","CL_deltaPhiIn","CL_deltaEtaIn","CL_detacalo","CL_see","CL_spp","CL_etawidth","CL_phiwidth","CL_e1x5e5x5","CL_HoE","CL_EoP","CL_IoEmIoP","CL_eleEoPout","CL_d0","CL_ip3d","CL_sigmaIetaIeta","CL_passConversionVeto","CL_isEcalDriven","CL_hnHits","CL_dZ", "CL_isoECAL", "CL_isoHCAL", "CL_isoTracker", "CL_isoECALRelat","CL_isoHCALRelat", "CL_isoTrackerRelat", "CL_isoECALRelatModif", "CL_MVA","CL_CombIsoHWW","CL_nonTrigMVA"
};
float lowerVal[33] = {
    0,      0.,         0.,         -1.,        -0.05,           -0.02,          -0.02,          0.,     0.,     0.,             0.,         0.,             0.,     0.,     0.,            0.,            -0.02,  -0.03,      0.,                 0.,                    0.,              0.,         -0.05, -1.,          -.1,            -.1,            -0.1,             0,                 0,                       -0.2,                -1.1,    0, -1.1
};
float higherVal[33] = {
    1.1,    1.3,         4.,         21.,        0.05,            0.02,           0.02,           0.05,   0.1,    0.05,            0.2,        0.8,           0.2,    5,    0.01,           5,           0.02,     0.03,    0.05,               2,                      2.,             5.,         0.05,     10.,           2,              2,              0.5,              0.1,                  0.1,                      0.3,             1.1,    50, 1.1
};
float binning[33] = {
    15,     15,        20,         11,         20,             20,             20,             20,     20,     20,             15,          15,            15,     20,     15,             20,            20,      20,     20,                 2,                      2,              5,          20,     40,             40,         40,                  40,             40,                 40,                     40,                     15,         200, 15
    
};




RooWorkspace* fitTheMC(int i, int j, RooWorkspace* ws){
    TH1F *theHisto = (TH1F*) myFileMC->Get(Form("mass_ptbin_%i_etabin_%i",i,j));
    TString signalPDF ="signalPass";
    TString tailPDF ="expP";

    if (i==5) {signalPDF ="signalHighBin"; tailPDF="landp";} // for the highest bin, the signal model is different
    RooDataHist histoMC("histoMC","histoMC",*(ws->var("mass")),Import(*theHisto));
    ws->pdf(signalPDF)->fitTo(histoMC);
    RooPlot* xframe = ws->var("mass")->frame(Title("MC model"),Bins(40)) ;
    histoMC.plotOn(xframe);
    ws->pdf(signalPDF)->plotOn(xframe,LineColor(kBlue));
    ws->pdf(signalPDF)->plotOn(xframe,Components(*ws->pdf(tailPDF)), LineStyle(kDashed),LineColor(kRed));

    TCanvas *c0 = new TCanvas("c0","coucou",600,600);
    xframe->Draw();
    c0->Print(Form("MCplots/mass_ptbin_%i_etabin_%i.png",i,j));
    delete c0;

    // now create the WS for fitting the data
    RooWorkspace* theWS = new RooWorkspace(Form("dataWS_%i_etabin_%i",i,j));
    theWS->factory(Form("alpha[%f]",ws->var("alpha")->getVal()));
    theWS->factory("fSigAll[0.7,0,1]");
    theWS->factory("largerResPass[1.,0.5,2.]");
    theWS->factory("numTot[61655,0,1e10]");
    theWS->factory("scaleTp[1,0.9,1.1]");
    theWS->factory(Form("lep[%f]",ws->var("lep")->getVal()));
    theWS->factory(Form("mean[%f]",ws->var("mean")->getVal()));
    theWS->factory(Form("n[%f]",ws->var("n")->getVal()));
    theWS->factory(Form("scale[%f]",ws->var("scale")->getVal()));
    theWS->factory(Form("sigma[%f]",ws->var("sigma")->getVal()));
    theWS->factory(Form("vFrac[%f]",ws->var("vFrac")->getVal()));
    theWS->factory(Form("Lmp[%f]",ws->var("Lmp")->getVal()));
    theWS->factory(Form("wp[%f]",ws->var("wp")->getVal()));
    theWS->factory("nSignalPass[10,0,1e10]");
    theWS->factory("nBkgPass[10,0,1e10]");
    theWS->factory("mass[60,120]");
    theWS->factory("expr::NewMean1p('mean*scaleTp',mean,scaleTp)");
    theWS->factory("expr::NewSigma1p('sigma*largerResPass',sigma,largerResPass)");
    theWS->factory("CBShape::cbs(mass, scale, NewSigma1p, alpha, n)");
    theWS->factory("RooBreitWigner::vs(mass, NewMean1p, NewSigma1p)");
    theWS->factory("FCONV::convPass(mass,vs,cbs)");
    theWS->factory("Exponential::expP(mass, lep)");
    theWS->factory("SUM::signalPass(vFrac*convPass, expP)");
    theWS->factory("RooBernstein::backgroundPass(mass,{a0[10,0,50],a1[1,0,50],a2[1,0,50],a3[1,0,50]})");
    theWS->factory("SUM::tot(nSignalPass*signalPass,nBkgPass*backgroundPass)");
    theWS->factory("RooLandau::landp(mass, Lmp,wp)");
    theWS->factory("SUM::signalPassHighBin(vFrac*convPass, landp)");
    theWS->factory("SUM::totHighBin(nSignalPass*signalPassHighBin,nBkgPass*backgroundPass)");

    
   /* RooPlot* xframe1 = theWS->var("mass")->frame(Title("MC model"),Bins(40)) ;
    theWS->pdf("totHighBin")->plotOn(xframe1,LineColor(kBlue));
                   
    TCanvas *c0 = new TCanvas("c0","coucou",600,600);
    xframe1->Draw();*/
    return theWS;

}


fitTheDATA(int i, int j, RooWorkspace* ws){ //try the binned fit to see if it works fine
    TH1F *theHisto = (TH1F*) myFileDATA->Get(Form("mass_ptbin_%i_etabin_%i",i,j));
    TString signalPDF ="tot";
    if (i==5) signalPDF ="totHighBin";
    RooDataHist histoDATA("histoDATA","histoDATA",*(ws->var("mass")),Import(*theHisto));
    ws->pdf(signalPDF)->fitTo(histoDATA);
    RooPlot* xframe = ws->var("mass")->frame(Title("DATA fit"),Bins(40)) ;
    histoDATA.plotOn(xframe);
    ws->pdf(signalPDF)->plotOn(xframe,LineColor(kBlue));
    ws->pdf(signalPDF)->plotOn(xframe,Components(*ws->pdf("backgroundPass")), LineStyle(kDashed),LineColor(kRed));
    
    TCanvas *c0 = new TCanvas("c0","coucou",600,600);
    xframe->Draw();
    c0->Print(Form("dataPlots/mass_ptbin_%i_etabin_%i.png",i,j));
    delete c0;
}



fitTheDATAunbinned(int i, int j, RooWorkspace* ws){
    TString signalPDF ="tot";
    if (i==5) signalPDF ="totHighBin";
    
    int nbEntree = chain->GetEntries();
    float ptBins[7] = {7,10,15,20,30,50,200};
    float etaBins[4] = {0, 1.4442, 1.556, 2.5};
    TString theCut = Form("pt>%f&&pt<%f&&absSCeta>%f&&absSCeta<%f",ptBins[i], ptBins[i+1], etaBins[j],etaBins[j+1]); //eventMatched
    if (i==0) theCut+="&&CL_PFchargedIso<0.4";
    TTree *theTree = chain->CopyTree(theCut,"",nbEntree,0);
    TH1F *theMass = new TH1F("theMass","",60,60,120);
    theTree->Draw("mass>>theMass");
    theMass->Draw();

    RooArgSet allTheVars;
    allTheVars.add(*(ws->var("mass")));
#if defined(isMC)
    allTheVars.add(*(ws->var("weight")));
#endif
   for (int k=0 ; k < nbSpectator ; k++){
        if (k>=23&&k<=29) continue;
        allTheVars.add(*(ws->var(theVars[k])));
    }
    
#if defined(isMC)
    RooDataSet treeDATA("treeDATA","treeDATA",allTheVars,Import(*theTree),WeightVar(*(ws->var("weight"))));
#else
    RooDataSet treeDATA("treeDATA","treeDATA",allTheVars,Import(*theTree));
#endif
    ws->pdf(signalPDF)->fitTo(treeDATA,Extended() );
    RooPlot* xframe2 = ws->var("mass")->frame(Title("DATA fit"),Bins(60)) ;
    treeDATA.plotOn(xframe2);
    ws->pdf(signalPDF)->plotOn(xframe2,LineColor(kBlue));
    ws->pdf(signalPDF)->plotOn(xframe2,Components(*ws->pdf("backgroundPass")), LineStyle(kDashed),LineColor(kRed));
   
    TCanvas *c0 = new TCanvas("c0","coucou",600,600);
    xframe2->Draw();
    c0->Print(Form("dataPlotsUnbined/mass_ptbin_%i_etabin_%i.png",i,j));

    // now preparing for sPlot
    //recup the signal and bg yields the parameter :
    RooRealVar *nSignalPass = ws->var("nSignalPass"); 
    RooRealVar *nBkgPass = ws->var("nBkgPass"); 

    //fix all the other parameter of the model
    RooRealVar *fSigAll = ws->var("fSigAll"); fSigAll->setConstant();
    RooRealVar *largerResPass = ws->var("largerResPass"); largerResPass->setConstant();
    RooRealVar *numTot = ws->var("numTot"); numTot->setConstant();
    RooRealVar *scaleTp = ws->var("scaleTp"); scaleTp->setConstant();
    RooRealVar *a0 = ws->var("a0"); a0->setConstant();
    RooRealVar *a1 = ws->var("a1"); a1->setConstant();
    RooRealVar *a2 = ws->var("a2"); a2->setConstant();
    RooRealVar *a3 = ws->var("a3"); a3->setConstant();
    
    // compute the sPlot weights ! 
    RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
                                                 treeDATA, ws->pdf(signalPDF), RooArgList(*nSignalPass,*nBkgPass) );
    
    //create the dataset with signal weights
    RooDataSet * dataw_z = new RooDataSet(treeDATA.GetName(),treeDATA.GetTitle(),&treeDATA,*treeDATA.get(),0,"nSignalPass_sw") ;
    //plot and recup the histos for signal shape
    for (int k=0 ; k < nbSpectator ; k++){
        if (k>=23&&k<=29) continue;
        RooPlot* frame2 = ws->var(theVars[k])->frame(Title("COUCOU fit"),Bins(binning[k])) ;
        dataw_z->plotOn(frame2, DataError(RooAbsData::SumW2), Name("hugues") ) ;
        RooHist* histo = (RooHist*) frame2->findObject("hugues") ;
        Double_t xdata,ydata ;
        saveHisto->cd();
        histo->Write(theVars[k]+Form("_ptbin_%i_etabin_%i",i,j)+"_signal");
        delete frame2;
    }
    
    //create the dataset with bg weights
    RooDataSet * dataw_bg = new RooDataSet(treeDATA.GetName(),treeDATA.GetTitle(),&treeDATA,*treeDATA.get(),0,"nBkgPass_sw") ;
    //plot and recup the histos for bg shape
    for (int k=0 ; k < nbSpectator ; k++){
        if (k>=23&&k<=29) continue;
        RooPlot* frame2 = ws->var(theVars[k])->frame(Title("COUCOU fit"),Bins(binning[k])) ;
        dataw_bg->plotOn(frame2, DataError(RooAbsData::SumW2), Name("hugues")  ) ;
        RooHist* histo = (RooHist*) frame2->findObject("hugues") ;
        saveHisto->cd();
        histo->Write(theVars[k]+Form("_ptbin_%i_etabin_%i",i,j)+"_bg");
        delete frame2;
    }
    
}

doTheFit(){
    // fill the tree with the data
    TString theWorkingDir = "/afs/cern.ch/work/h/hbrun/CLshape/";
    
    saveHisto = new TFile("outputHistos_runAll.root","RECREATE");
    
#if !defined(isMC)
        chain->Add(theWorkingDir+"TnP_CLtree_electrons_runA.root");
        chain->Add(theWorkingDir+"TnP_CLtree_electrons_runB_part0_div0.root");
        chain->Add(theWorkingDir+"TnP_CLtree_electrons_runB_part0_div1.root");
        chain->Add(theWorkingDir+"TnP_CLtree_electrons_runB_part1.root");

        chain->Add(theWorkingDir+"TnP_CLtree_electrons_runCv1_RERECO.root");
        chain->Add(theWorkingDir+"TnP_CLtree_electrons_runCv2_part0.root");
        chain->Add(theWorkingDir+"TnP_CLtree_electrons_runCv2_part1.root");
        chain->Add(theWorkingDir+"TnP_CLtree_electrons_runCv2_part2.root");
        chain->Add(theWorkingDir+"TnP_CLtree_electrons_runCv2_part3.root");

        chain->Add(theWorkingDir+"TnP_CLtree_electrons_runD_part1.root");
        chain->Add(theWorkingDir+"TnP_CLtree_electrons_runD_part2.root");
        chain->Add(theWorkingDir+"TnP_CLtree_electrons_runD_part3.root");
    
#else
    chain->Add(theWorkingDir+"TnP_CLtree_DYJet_ter_part1.root");
    chain->Add(theWorkingDir+"TnP_CLtree_DYJet_ter_part2.root");
    chain->Add(theWorkingDir+"TnP_CLtree_DYJet_ter_part3.root");
#endif
    
    //workspace with the signal model
    RooWorkspace* myWorkSpace = new RooWorkspace("myWorkSpace");
    myWorkSpace->factory("mass[60,120]");
    myWorkSpace->factory("mean[92,85,95]");
    myWorkSpace->factory("scale[1,0.2,1.8]");
    myWorkSpace->factory("sigma[3,1,10]");
    myWorkSpace->factory("CBShape::cbs(mass, scale, sigma,alpha[3., 0.5, 5.], n[1, 0., 100.])");
    myWorkSpace->factory("RooBreitWigner::vs(mass, mean, sigma)");
    myWorkSpace->factory("FCONV::convPass(mass,vs,cbs)");
    myWorkSpace->factory("Exponential::expP(mass, lep[0,-5,5])");
    myWorkSpace->factory("SUM::signalPass(vFrac[0.8,0.5,1]*convPass, expP)");
    myWorkSpace->factory("RooLandau::landp(mass, Lmp[100,90,105],wp[1,0,10])");
    myWorkSpace->factory("SUM::signalHighBin(vFrac*convPass, landp)");


    // fit the MC
   /*RooWorkspace* theWsForData = fitTheMC(0,0,myWorkSpace);
    for (int k=0 ; k<nbSpectator ; k++){
        if (k>=23&&k<=29) continue;
        float middle = 0.5*(higherVal[k]-lowerVal[k]);
        theWsForData->factory(theVars[k]+Form("[%f,%f,%f]",middle,lowerVal[k],higherVal[k]));
        
    }
    fitTheDATAunbinned(0,0,theWsForData);*/
    //theWsForData->Print();
 for (int i = 0 ; i < 6 ; i++){//6
        for (int j = 0 ; j < 3 ; j++){
	    cout << "poisson, i=" << i << ", j=" << j << endl;
            if (j==1) continue;
            RooWorkspace* theWsForData = fitTheMC(i,j,myWorkSpace);
            for (int k=0 ; k<nbSpectator ; k++){
            float middle = 0.5*(higherVal[k]-lowerVal[k]);
            theWsForData->factory(theVars[k]+Form("[%f,%f,%f]",middle,lowerVal[k],higherVal[k]));
            }
#if defined(isMC)
            theWsForData->factory("weight[1,0,1.5]");
#endif
            fitTheDATAunbinned(i,j,theWsForData);
        }
    }


    
    
    
    
    
    
    
    
    
}

