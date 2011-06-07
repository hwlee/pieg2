// @(#)ClasTool/TPiEG2.cc:$Name:  $:$Id: TPiEG2.cc,v 1.10 2005/07/28 19:14:17 holtrop Exp $
// Author: Maurik Holtrop <http://www.physics.unh.edu/~maurik>

/****************************************************************************
 * Copyright (C) CopyLeft  This code is freely available to all. *
 *                                                                          *
 * Documentation  : TPiEG2.html                                          *
 *                  (available also at:                                     *
 *                   http://www.physics.unh.edu/maurik/ClasTool)            *
 * Created on     :  30/1/2006 (UNH)                                        *
 * Initial Authors:  Maurik Holtrop (UNH)                                   *
 ***************************************************************************/
// Author:  Maurik Holtrop <http://www.physics.unh.edu/~maurik>
//____________________
//Begin_Html <!--
/* -->
</pre><H1> TPiEG2 </H1>
<pre>
<!-- */
// --> End_Html

#include "TPiEG2.h"
#define MaxNpi 20

;
ClassImp(TPiEG2) // A class for checking the PID.
;  

TPiEG2::TPiEG2()
{
  // Initialize code.

  // Add a signal handler so you can stop running with crtl-C
  Interrupt = new TInterrupt();
  gSystem->AddSignalHandler(Interrupt);  // Setup an intterupt handler.
  
  // Allow for timing the loop
  Time = new TStopwatch(); 

  H = new TObjArray();
  treeZCut = new TTree("ZCut","For all data"); 

//  rc = new TRadCor();
//  fC = new THapradConfig();
  InitHistos();

  Sieve = new TPartSieve(this);
  Sieve->SortOn();
  Sieve->ConfidenceOn();

  SHists = new TPartSieveHists(Sieve); // Turn off by commenting.

  // Initialize the PARTICLE ID Cuts.
  Sieve->fEC_U_Cut=40.;
  Sieve->fEC_U_Cut_width=5.;
  Sieve->fEC_V_Cut=360.;
  Sieve->fEC_V_Cut_width=5.;
  Sieve->fEC_W_Cut=390.;
  Sieve->fEC_W_Cut_width=5.;
  Sieve->fEC_Eio_slope=1.0;
  Sieve->fEC_Eio_off=0.156;
  Sieve->fEC_Ein_cut=0.06;
  Sieve->fEC_Eout_cut=0.06;
  Sieve->fCC_Phe_cut=25.0;
  Sieve->fCC_Phe_width=3.;
  Sieve->fConf_Pim_DT_min=-0.5;
  Sieve->fConf_Pim_DT_max=0.5;
  Sieve->fConf_Pip_DT_min=-0.5;
  Sieve->fConf_Pip_DT_max=0.5;
  Sieve->fConf_Km_DT_min=-1.0; //-0.25;
  Sieve->fConf_Km_DT_max=1.0; //0.25;
  Sieve->fConf_Kp_DT_min=-1.0; //-0.25;
  Sieve->fConf_Kp_DT_max=1.0; //0.25;
  Sieve->fConf_Prot_DT_min=-0.5;
  Sieve->fConf_Prot_DT_max=1.0;
  Sieve->fConf_Neut_DT_min=-0.8;
  Sieve->fConf_Neut_DT_max=1.0;

  Sieve->fConf_El_DC_chi2=4.;
  Sieve->fConf_Prot_DC_chi2=5.;
  Sieve->fConf_Pim_DC_chi2=5.;
  Sieve->fConf_Pip_DC_chi2=5.;
  Sieve->fConf_Km_DC_chi2=5.;
  Sieve->fConf_Kp_DC_chi2=5.;

  Sieve->fConf_El_DC_chi2_width=0.;
  Sieve->fConf_Prot_DC_chi2_width=0.;
  Sieve->fConf_Pim_DC_chi2_width=0.;
  Sieve->fConf_Pip_DC_chi2_width=0.;
  Sieve->fConf_Km_DC_chi2_width=0.;
  Sieve->fConf_Kp_DC_chi2_width=0.;

  Sieve->SetBadPaddle(1,25);
  Sieve->SetBadPaddle(3,11);
  Sieve->SetBadPaddle(3,37);
  Sieve->SetBadPaddle(5,20);
  Sieve->SetBadPaddle(5,24);

  for(Int_t i=1;i<7;i++)
  {
    Sieve->SetBadPaddle(i,40); // Cut all backward paddles.
    Sieve->SetBadPaddle(i,41);
    Sieve->SetBadPaddle(i,42);
    Sieve->SetBadPaddle(i,43);
    Sieve->SetBadPaddle(i,44);
    Sieve->SetBadPaddle(i,45);
    Sieve->SetBadPaddle(i,46);
    Sieve->SetBadPaddle(i,47);
    Sieve->SetBadPaddle(i,48);
  }

  gFiltered_data=1;

}

void TPiEG2::InitHistos(void)
{
  // Initialize the Histograms.
  H->Add(Elec_vert_z = new TH1D("Elec_vert_z","Z Vertex electron",500,-40.,0.));
  H->Add(Elec_vert_xy = new TH2D("Elec_vert_xy","X Y Vertex electron",300,-10.,10.,300,-10.,10.));
  H->Add(EMass = new TH1F("EMass", "Electron Mass from EVNT Bank", 200, -.001, 0.001));
  H->Add(SampFrac = new TH1F("SampFrac", "Sampling Fraction for the EC", 200, 0, 2.0));

  H->Add(PiPMom = new TH1D("PiPMom", "Pi+ Momentum", 200, 0, 5.0));
  H->Add(PiPCos = new TH1D("PiPCos", "Cos of Pi+ Angular Distribution", 50, -1.0, 1.0));
  H->Add(PiMMom = new TH1D("PiMMom", "Pi- Momentum", 200, 0, 5.0));
  H->Add(PiMCos = new TH1D("PiMCos", "Cos of Pi- Angular Distribution", 50, -1.0, 1.0));
}

void TPiEG2::DeleteHistos(void)
{
  // Delete all histograms by calling H->Delete()
  H->Delete();
  treeZCut->Delete(); 
}

void TPiEG2::ClearHistos(void)
{
  // Clear all histograms.
  TIter next(H);
  TH1 *hist;

  while((hist = (TH1 *)next()))
  {
    hist->Reset();
  }
}

void TPiEG2::Write(void)
{
  // Write all histograms in H 
  // to the currently open file

  treeZCut->Write();

  TIter next(H);
  TH1 *hist;
  
  gDirectory->mkdir("Histograms");
  gDirectory->cd("Histograms");

  while((hist = (TH1 *)next()))
  {
    hist->Write();
  }
  gDirectory->cd("../");

  if(SHists!=NULL)
  {
    SHists->Write();
  }
}


Int_t TPiEG2::Run(Int_t Run_num, Int_t Max_evt) 
{
  Int_t Eg2Mc=41000;  //run #  >Eg2Mc : data  <Eg2Mc : MC
  Int_t runnum=Run_num;
  Float_t BeamEnergy=5.01483;
  TVector3 zaxis(0., 0., 1.0);

  TEVNTClass *elec;
  TEVNTClass *pion;
  TSCPBClass *sc;
  TECPBClass *ec;

  Bool_t Good_Electron;
  Int_t iEvent, iEvent_ori;

  TVector3 EMom3v, EinMom3v, q3v, pim3v;
  TVector3 EV3, EV3_corr;
  TLorentzVector V4Elin,V4Elout,V4VirtualPhoton,V4MMIni,V4MMSub,V4Proton,V4Neutron;
  TLorentzVector pi4vec,V4piMM,V4qminusph;

  Float_t p_elec;
  Float_t u,v,w;
  Float_t sampfrac;
  Float_t e_mass;
  Float_t Q2,W;
  Double_t EloutCos, EloutMom;
  Float_t nproton;
  Float_t nneutron;
  Int_t target;  //0 for nothing, 1 for D2, 2 for Solid
  Double_t Tree_EV3_x,Tree_EV3_y,Tree_EV3_z;
  Double_t Tree_EV3_x_cor,Tree_EV3_y_cor,Tree_EV3_z_cor;
  Float_t npiplus;
  Float_t npiminus;
  Double_t magMM, MM2, xbj, nu;  //Missing mass
  Int_t Sector;
  Double_t pip[MaxNpi];
  Double_t picos[MaxNpi];
  Int_t piq[MaxNpi];
  Double_t pix[MaxNpi], piy[MaxNpi], piz[MaxNpi], zh[MaxNpi], th[MaxNpi], pt[MaxNpi], phih[MaxNpi];
  Int_t npip;
  Int_t npim;
  Int_t npi;
  Int_t index;
  Double_t mhp, mhm;
  Double_t phi1, phi2;

//  TRadCor *rc;
//  Double_t m;
  Double_t factor[MaxNpi];

  treeZCut->Branch("runnum",&runnum,"runnum/I");
  treeZCut->Branch("ievent",&iEvent_ori,"ievent/I");

  treeZCut->Branch("W",&W,"W/F"); 
  treeZCut->Branch("Q2",&Q2,"Q2/F");
  treeZCut->Branch("x",&xbj,"x/D");
  treeZCut->Branch("EloutMom",&EloutMom,"EloutMom/D");
  treeZCut->Branch("EloutCos",&EloutCos,"EloutCos/D");   

  treeZCut->Branch("Sector",&Sector,"Sector/I");
  treeZCut->Branch("target",&target,"target/I"); 
  treeZCut->Branch("EV3X",&Tree_EV3_x,"EV3X/D");
  treeZCut->Branch("EV3Y",&Tree_EV3_y,"EV3Y/D");
  treeZCut->Branch("EV3Z",&Tree_EV3_z,"EV3Z/D");
  treeZCut->Branch("EV3X_cor",&Tree_EV3_x_cor,"EV3X_cor/D");
  treeZCut->Branch("EV3Y_cor",&Tree_EV3_y_cor,"EV3Y_cor/D");
  treeZCut->Branch("EV3Z_cor",&Tree_EV3_z_cor,"EV3Z_cor/D");

  treeZCut->Branch("magMM",&magMM,"magMM/D");
  treeZCut->Branch("MM2",&MM2,"MM2/D");

  treeZCut->Branch("npip",&npip,"npip/I");
  treeZCut->Branch("npim",&npim,"npim/I");
  treeZCut->Branch("npi",&npi,"npi/I");
  treeZCut->Branch("piq",piq,"piq[npi]/I");
  treeZCut->Branch("pip",pip,"pip[npi]/D");
  treeZCut->Branch("pix",pix,"pix[npi]/D");
  treeZCut->Branch("piy",piy,"piy[npi]/D");
  treeZCut->Branch("piz",piz,"piz[npi]/D");
  treeZCut->Branch("picos",picos,"picos[npi]/D");
  treeZCut->Branch("z",zh,"z[npi]/D");
  treeZCut->Branch("pt",pt,"pt[npi]/D");
  treeZCut->Branch("phih",phih,"phih[npi]/D");

  treeZCut->Branch("rc",factor,"rc[npi]/D");

  // Start up function; initialize counter.

  cout << "Running ...." << endl;

  Time->Start(kTRUE);

  // Loop through events...
  for(iEvent=0;iEvent<Max_evt && Next()==0; iEvent++)
  {
    iEvent_ori=GetHEADER()->GetNEvent();

    // Print out number of events analyzed every 100K events.
    if(iEvent % 10000 == 0)
    {
      printf("Analyzed %9d events, now at %9d \n",iEvent,iEvent_ori);
    }
    
    // Make sure that there are particles associated with the event.
    if(GetNPart()<1)
    {
      printf("No events found.\n");
      continue;
    }

    if(!gFiltered_data)
    {
      cout << "Refilter?\n";
      cout << "Are you nuts!\n";
      return(0);
    }

    Good_Electron=true; // Need to change this to some checking routine

    if(Good_Electron)  // From here on, only interested in events with good electron.
    {

      // Start by looking for electrons.  Assume electron is first particle in 
      // EVNT bank.  "elec" is pointer to electron's info; store electron 
      // momentum.
      elec = GetEVNT(0);
      p_elec = elec->GetMomentum();

      EMom3v=GetPart4Vector(0).Vect();

      sc= (TSCPBClass *)GetBankRow("SCPB",elec->GetSCidx());   

      // Look at raw particle mass as a diagnostic.
      e_mass = elec->GetMass();
      EMass->Fill(e_mass);

      // Get u,v,w coordinates of calorimeter hit; calculate various 
      ec =(TECPBClass *)GetBankRow("ECPB",elec->GetECidx());

      // cout << "Electron ec entry:  " << elec->GetECidx() << endl;

      // Look at sampling fraction.
      sampfrac=ec->Etot/p_elec/0.29;
      SampFrac->Fill(sampfrac);
      
      ec->GetUVW(&u,&v,&w);

      // kinematic quantities (W, Q2, etc).
      SetBeamEnergy(BeamEnergy);
      GetInvariants(&Q2,&W);
      if(W>0) 
      {
        W=TMath::Sqrt(W); 
      }
      else 
      {
        W=-1.;
      }

      // Calculate incident, scattered electron p vectors and p vector of 
      // virtual photon.
      EloutCos = EMom3v.Z()/EMom3v.Mag();  
      EloutMom = EMom3v.Mag();  
      V4Elout = GetPart4Vector(0);

      // Record scattering vertex of electron.
      EV3 = GetPartVertex(0);  

      EV3_corr = vertex_corr(EV3, EMom3v);

      Tree_EV3_x = EV3.X();  
      Tree_EV3_y = EV3.Y();  
      Tree_EV3_z = EV3.Z();  
      
      if (runnum>Eg2Mc)
      {
        Tree_EV3_x_cor = EV3_corr.X();  
        Tree_EV3_y_cor = EV3_corr.Y();  
        Tree_EV3_z_cor = EV3_corr.Z();  
      }
      else //MC data : don't need correction
      {
        Tree_EV3_x_cor = EV3.X();  
        Tree_EV3_y_cor = EV3.Y();  
        Tree_EV3_z_cor = EV3.Z();  

        EV3_corr = EV3;
      }

      if(EV3_corr.Z()<-31.8) target=0; 
      if(EV3_corr.Z()>-31.8 && EV3_corr.Z()<-28.4) target=1; 
      if(EV3_corr.Z()>-28.4 && EV3_corr.Z()<-25.7) target=0; 
      if(EV3_corr.Z()>-25.7 && EV3_corr.Z()<-24.0) target=2; 
      if(EV3_corr.Z()>-24.0) target=0; 

      Elec_vert_z->Fill(EV3_corr.Z());
      Elec_vert_xy->Fill(EV3_corr.X(),EV3_corr.Y());

      if(target > 0) //Z cut from here
      {	       	
	
	Sieve->SieveEvent(); // Run the particle sorter. (cost = 12 micro sec/event taro CPU).
	if(SHists!=NULL) SHists->Fill();      // Fill Standard Histograms.

        nproton=0;
        nneutron=0;

        if (Sieve->fNPart[kProton]>0) nproton=nproton+Sieve->fNPart[kProton];
        if (Sieve->fNPart[kNeutron]>0) nneutron=nneutron+Sieve->fNPart[kNeutron];
        if (Sieve->fNPart[kDeuteron]>0)
        {
          nproton=nproton+Sieve->fNPart[kDeuteron];
          nneutron=nneutron+Sieve->fNPart[kDeuteron];
        }
        if (Sieve->fNPart[kHe3]>0)
        {
          nproton=nproton+2*Sieve->fNPart[kHe3];
          nneutron=nneutron+Sieve->fNPart[kHe3];
        }
        if (Sieve->fNPart[kHe4]>0)
        {
          nproton=nproton+2*Sieve->fNPart[kHe4];
          nneutron=nneutron+2*Sieve->fNPart[kHe4];
        }

	  
        if (nproton+nneutron==0) nproton=1;

        V4Elin.SetXYZM(0,0,GetBeamEnergy(),ClasTool::fgParticle_Mass[kElectron]);
        V4Proton.SetXYZM(0,0,0,nproton*(ClasTool::fgParticle_Mass[kProton]));  
        V4Neutron.SetXYZM(0,0,0,nneutron*(ClasTool::fgParticle_Mass[kNeutron]));  
        V4VirtualPhoton=V4Elin-V4Elout;
        V4MMIni = V4VirtualPhoton+V4Proton+V4Neutron;

        //Elastic peak Missing Mass  
        V4MMSub = V4MMIni;  
        for(Int_t imm=1;imm<GetNPart();imm++)
        {   // : mm=0 -> outgoing electron
          if (Sieve->GetPartConf(imm)>0)
          {
            V4MMSub=V4MMSub - GetPart4Vector(imm);  
          }
        }  
        MM2 = V4MMSub.M2();
        magMM = V4MMSub.M(); 

        xbj = Q2/2./ClasTool::fgParticle_Mass[kProton]/(V4Elin.E()-V4Elout.E());
        Sector = sc->GetSector();
        nu=V4Elin.E()-V4Elout.E();
        EinMom3v=V4Elin.Vect();
        q3v=V4VirtualPhoton.Vect();
        phi1=(V4Elout.Vect()).Phi();

        // Here's where I start looking for Pions.  Let's first record the
        // number of pions of each charge species in the event.
        npiplus=Sieve->GetNIdx(kPion_Plus);
        npiminus=Sieve->GetNIdx(kPion_Minus);
  
        npip=0;
        npim=0;
        mhp=ClasTool::fgParticle_Mass[kPion_Plus];
        mhm=ClasTool::fgParticle_Mass[kPion_Minus];

        if(npiplus>0)
        {  
          for(Int_t ipip=0; ipip<npiplus; ipip++)
          {  
            index = Sieve->GetIndexIdx(kPion_Plus, ipip);  
            pi4vec=GetPart4Vector(index);  
            pion=GetEVNT(Sieve->GetIndexIdx(kPion_Plus, ipip));  
            pip[ipip]=pion->GetMomentum(); 
            piq[ipip]=1;
            pim3v=pion->GetMomVec();
            pix[ipip]=pim3v.X();
            piy[ipip]=pim3v.Y();
            piz[ipip]=pim3v.Z();
            picos[ipip]=cos(zaxis.Angle(pim3v)); 
            phi2=pim3v.Phi();
            npip=ipip+1;
            zh[ipip]=pi4vec.E()/nu;
            V4qminusph=V4VirtualPhoton-pi4vec;
            th[ipip]=V4qminusph.M2();
            pt[ipip]=sqrt(pi4vec.E()*pi4vec.E()-mhp*mhp-pow((th[ipip]-mhp*mhp+Q2+2*pi4vec.E()*nu),2)/4/(nu*nu+Q2));
            phih[ipip]=TMath::RadToDeg()*(EinMom3v.Cross(EMom3v)).Angle(q3v.Cross(pim3v));
            if (phi1>phi2)
            {
              phih[ipip]=-phih[ipip];
            }

//            m = TMath::Power((kMassNeutron + kMassPion),2);
//            rc->SetICHARGE(1);
//            fC = rc->GetConfig();
//            cout << "pip  " << fC->ICHARGE() << endl;
//if (npiplus+npiminus==1){
//            factor[ipip] = rc->GetRCFactor(GetBeamEnergy(),xbj,Q2,zh[ipip],pt[ipip],phih[ipip],m);
//            std::cout << std::endl << "RC factor: " << factor << std::endl;
//}
//else  factor[ipip] = -1;
            factor[ipip] = -1;
            PiPMom->Fill(pip[ipip]);
            PiPCos->Fill(picos[ipip]);
	  }  
	}  
        if(npiminus>0)
        {  
          for(Int_t ipim=0; ipim<npiminus; ipim++)
          {  
            index = Sieve->GetIndexIdx(kPion_Minus, ipim);  
            pi4vec=GetPart4Vector(index);  
            pion=GetEVNT(Sieve->GetIndexIdx(kPion_Minus, ipim));  
            pip[ipim+npip]=pion->GetMomentum(); 
            pim3v=pion->GetMomVec();
            pix[ipim+npip]=pim3v.X();
            piy[ipim+npip]=pim3v.Y();
            piz[ipim+npip]=pim3v.Z();
            picos[ipim+npip]=cos(zaxis.Angle(pim3v));           
            piq[ipim+npip]=-1; 
            phi2=pim3v.Phi();
            npim=ipim+1;
            zh[ipim+npip]=pi4vec.E()/nu;
            V4qminusph=V4VirtualPhoton-pi4vec;
            th[ipim+npip]=V4qminusph.M2();
            pt[ipim+npip]=sqrt(pi4vec.E()*pi4vec.E()-mhm*mhm-pow((th[ipim+npip]-mhm*mhm+Q2+2*pi4vec.E()*nu),2)/4/(nu*nu+Q2));
            phih[ipim+npip]=TMath::RadToDeg()*(EinMom3v.Cross(EMom3v)).Angle(q3v.Cross(pim3v));
            if (phi1>phi2)
            {
              phih[ipim+npip]=-phih[ipim+npip];
            }

//            m = TMath::Power((kMassNeutron + kMassPion),2);
//            rc->SetICHARGE(2);
//            fC = rc->GetConfig();
//            cout << "pim  " << fC->ICHARGE() << endl;
//            factor[ipim+npip] = rc->GetRCFactor(GetBeamEnergy(),xbj,Q2,zh[ipim+npip],pt[ipim+npip],phih[ipim+npip],m);
//            std::cout << std::endl << "RC factor: " << factor << std::endl;
            factor[ipim+npip] = -1;
            PiMMom->Fill(pip[ipim+npip]);
            PiMCos->Fill(picos[ipim+npip]);
 	  }  
        }  
        npi=npip+npim;
        treeZCut->Fill(); // Fill tree

      } // Vertex Test.
    } // End good electron test.
    
    
    if(Interrupt->IsInterrupted()) break;
    
  } // END EVENT LOOP.

  Time->Stop();

  printf("Processed %d events in %6.2f sec, %6.2f CPU sec. \n",iEvent,Time->RealTime(),Time->CpuTime());
  printf("Rate is %7.1f events/sec, %7.1f events/CPUsec \n",iEvent/Time->RealTime(),iEvent/Time->CpuTime());
  
  return(iEvent);
}
  
TVector3 TPiEG2::vertex_corr(TVector3 v, TVector3 p)
{
  // corrected vertex
  TVector3 v_corr;
  // page 26 of:
  // http://www.jlab.org/Hall-B/secure/e1-6/mauri/web/pi0e16/main/anote.pdf

  Double_t s0,sp,sv;
  /*static*/ Double_t n[3][6];
  n[0][0]=1.;
  n[1][0]=0.;

  n[0][1]=0.5 ;
  n[1][1]=0.866025388;

  n[0][2]=-0.5 ;
  n[1][2]=0.866025388;

  n[0][3]=-1. ;
  n[1][3]=0.;

  n[0][4]=-0.5 ;
  n[1][4]=-0.866025388;

  n[0][5]=0.5 ;
  n[1][5]=-0.866025388;

  // BEAM POSITION
  /*static*/ Double_t x0=0.043;
  /*static*/ Double_t y0=-0.33;
//  /*static*/ Double_t z0=0.;

  Double_t phi;
  phi=TMath::RadToDeg()*TMath::ATan2(p.Y(),p.X());

  int s=sect(phi);

  Double_t A;

  s0=x0*n[0][s]+y0*n[1][s]; //+z0*n[2][s];
  sp=p.X()*n[0][s]+p.Y()*n[1][s]; //+p.Z()*n[2][s];
  sv=v.X()*n[0][s]+v.Y()*n[1][s]; //+v.Z()*n[2][s];

  if(sp)
  {
    A=(s0-sv)/sp;
    v_corr.SetX(v.X()+ A*p.X());
    v_corr.SetY(v.Y()+ A*p.Y());
    v_corr.SetZ(v.Z()+ A*p.Z());
  }
  return v_corr;
}

int TPiEG2::sect(Double_t phi)
{
  phi+=30.;
  while(phi<0) phi+=360.;
  return (int) TMath::Floor(phi/60.);
} 
