//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Dec 11 18:31:07 2022 by ROOT version 6.22/00
// from TTree CoalData/CoalData DST Tree
// found on file: CoalData_18.root
//////////////////////////////////////////////////////////

#ifndef CoalData_h
#define CoalData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class CoalData {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Event_nevent;
   Int_t           Event_multi;
   Int_t           BaryonNum[26];   //[multi]
   Float_t         HadronX[26];   //[multi]
   Float_t         HadronY[26];   //[multi]
   Float_t         HadronPx[26];   //[multi]
   Float_t         HadronPy[26];   //[multi]
   Float_t         Distance[26];   //[multi]
   Float_t         Parton0BaryonNum[26];   //[multi]
   Float_t         Parton0X[26];   //[multi]
   Float_t         Parton0Y[26];   //[multi]
   Float_t         Parton1BaryonNum[26];   //[multi]
   Float_t         Parton1X[26];   //[multi]
   Float_t         Parton1Y[26];   //[multi]
   Float_t         Parton2BaryonNum[26];   //[multi]
   Float_t         Parton2X[26];   //[multi]
   Float_t         Parton2Y[26];   //[multi]

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_BaryonNum;   //!
   TBranch        *b_HadronX;   //!
   TBranch        *b_HadronY;   //!
   TBranch        *b_HadronPx;   //!
   TBranch        *b_HadronPy;   //!
   TBranch        *b_Distance;   //!
   TBranch        *b_Parton0BaryonNum;   //!
   TBranch        *b_Parton0X;   //!
   TBranch        *b_Parton0Y;   //!
   TBranch        *b_Parton1BaryonNum;   //!
   TBranch        *b_Parton1X;   //!
   TBranch        *b_Parton1Y;   //!
   TBranch        *b_Parton2BaryonNum;   //!
   TBranch        *b_Parton2X;   //!
   TBranch        *b_Parton2Y;   //!

   CoalData(TTree *tree=0);
   virtual ~CoalData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef CoalData_cxx
CoalData::CoalData(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("CoalData_18.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("CoalData_18.root");
      }
      f->GetObject("CoalData",tree);

   }
   Init(tree);
}

CoalData::~CoalData()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CoalData::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CoalData::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void CoalData::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event_nevent, &b_Event);
   fChain->SetBranchAddress("BaryonNum", BaryonNum, &b_BaryonNum);
   fChain->SetBranchAddress("HadronX", HadronX, &b_HadronX);
   fChain->SetBranchAddress("HadronY", HadronY, &b_HadronY);
   fChain->SetBranchAddress("HadronPx", HadronPx, &b_HadronPx);
   fChain->SetBranchAddress("HadronPy", HadronPy, &b_HadronPy);
   fChain->SetBranchAddress("Distance", Distance, &b_Distance);
   fChain->SetBranchAddress("Parton0BaryonNum", Parton0BaryonNum, &b_Parton0BaryonNum);
   fChain->SetBranchAddress("Parton0X", Parton0X, &b_Parton0X);
   fChain->SetBranchAddress("Parton0Y", Parton0Y, &b_Parton0Y);
   fChain->SetBranchAddress("Parton1BaryonNum", Parton1BaryonNum, &b_Parton1BaryonNum);
   fChain->SetBranchAddress("Parton1X", Parton1X, &b_Parton1X);
   fChain->SetBranchAddress("Parton1Y", Parton1Y, &b_Parton1Y);
   fChain->SetBranchAddress("Parton2BaryonNum", Parton2BaryonNum, &b_Parton2BaryonNum);
   fChain->SetBranchAddress("Parton2X", Parton2X, &b_Parton2X);
   fChain->SetBranchAddress("Parton2Y", Parton2Y, &b_Parton2Y);
   Notify();
}

Bool_t CoalData::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CoalData::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CoalData::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CoalData_cxx
