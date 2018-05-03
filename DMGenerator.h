#include <math.h>
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TRandom.h"
#include "FairPrimaryGenerator.h"
#include "DMGenerator.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoEltu.h"
#include "TGeoCompositeShape.h"

using std::cout;
using std::endl;
//Read events from Ntuples produced with MadDump

// -------------------------------------------------------------------------
// -----   Default constructor   -------------------------------------------
DMGenerator::DMGenerator() {}


// -------------------------------------------------------------------------
// -----     Constructor with name  ----------------------------------------
Bool_t DMGenerator::Init(const char* fileName) {
	return Init(fileName, 0);
}

// -------------------------------------------------------------------------
// -----     Constructor with name and firstEvent   ------------------------

Bool_t DMGenerator::Init(const char* fileName, const int firstEvent) {
  
	fLogger = FairLogger::GetLogger();
	if (0 == strncmp("/eos",fileName,4)){
		TString tmp = gSystem->Getenv("EOSSHIP");
		tmp+=fileName;
		fInputFile  = TFile::Open(tmp); 
		fLogger->Info(MESSAGE_ORIGIN,"Opening input file on eos %s",tmp.Data());
	}else{
		fInputFile = new TFile(fileName);
		fLogger->Info(MESSAGE_ORIGIN,"Opening input file %s",fileName);}
	
	if (fInputFile->IsZombie() or !fInputFile){
     		fLogger->Fatal(MESSAGE_ORIGIN, "Error opening input file");
	return kFALSE;}

	fTree = (TTree *)fInputFile->Get("dmtree");
	fNevents = fTree->GetEntries();
	fn = firstEvent;
	fTree->SetBranchAddress("Edm",&Edm);		// incoming DM energy
	fTree->SetBranchAddress("pxdm",&pxdm);
	fTree->SetBranchAddress("pydm",&pydm);
	fTree->SetBranchAddress("pzdm",&pzdm);
	fTree->SetBranchAddress("dmpdg",&dmpdg);	// incoming DM PDG code
	fTree->SetBranchAddress("dis",&dis);		// Is it a DIS event?
	fTree->SetBranchAddress("el",&el);		// Is it a EL event?
	fTree->SetBranchAddress("E_2ry",&E_2ry);	// outgoing particles 4-Momenta
	fTree->SetBranchAddress("px_2ry",&px_2ry);    
	fTree->SetBranchAddress("py_2ry",&py_2ry);
	fTree->SetBranchAddress("pz_2ry",&pz_2ry);
	fTree->SetBranchAddress("pdg_2ry",&pdg_2ry);	// pdg code of hadron
	fTree->SetBranchAddress("n_2ry",&n_2ry);	// nr of outgoing hadrons
	fFirst=kTRUE;
	return kTRUE;
}

// -----   Destructor   ----------------------------------------------------
DMGenerator::~DMGenerator()
{
	fInputFile->Close();
	fInputFile->Delete();
	delete fInputFile;
}

// -------------------------------------------------------------------------

Int_t DMGenerator::GetNevents()
{
	return fNevents;
}

// -------------------------------------------------------------------------

Double_t DMGenerator::MeanMaterialBudget(const Double_t *start, const Double_t *end, Double_t *mparam)
{
  //
  // Calculate mean material budget and material properties between
  //    the points "start" and "end".
  //
  // "mparam" - parameters used for the energy and multiple scattering
  //  corrections:
  //
  // mparam[0] - mean density: sum(x_i*rho_i)/sum(x_i) [g/cm3]
  // mparam[1] - equivalent rad length fraction: sum(x_i/X0_i) [adimensional]
  // mparam[2] - mean A: sum(x_i*A_i)/sum(x_i) [adimensional]
  // mparam[3] - mean Z: sum(x_i*Z_i)/sum(x_i) [adimensional]
  // mparam[4] - length: sum(x_i) [cm]
  // mparam[5] - Z/A mean: sum(x_i*Z_i/A_i)/sum(x_i) [adimensional]
  // mparam[6] - number of boundary crosses
  // mparam[7] - maximum density encountered (g/cm^3)
  // mparam[8] - equivalent interaction length fraction: sum(x_i/I0_i) [adimensional]
  // mparam[9] - maximum cross section encountered (mbarn)
  //
  //  Origin:  Marian Ivanov, Marian.Ivanov@cern.ch
  //
  //  Corrections and improvements by
  //        Andrea Dainese, Andrea.Dainese@lnl.infn.it,
  //        Andrei Gheata,  Andrei.Gheata@cern.ch
  //        Thomas Ruf,  Thomas.Ruf@cern.ch
  //

	mparam[0]=0; mparam[1]=1; mparam[2]=0; mparam[3]=0; mparam[4]=0;
	mparam[5]=0; mparam[6]=0; mparam[7]=0; mparam[8]=0; mparam[9]=0;

	Double_t bparam[7]; // total parameters
	Double_t lparam[7]; // local parameters
	Double_t mbarn = 1E-3*1E-24*TMath::Na(); // cm^2 * Avogadro

	for (Int_t i=0;i<7;i++) bparam[i]=0;

		if (!gGeoManager) {
    			//AliFatalClass("No TGeo\n");
			return 0.;}

	Double_t length;
	Double_t dir[3];
	length = TMath::Sqrt((end[0]-start[0])*(end[0]-start[0])+(end[1]-start[1])*(end[1]-start[1])+(end[2]-start[2])*(end[2]-start[2]));
	mparam[4]=length;
	if (length<TGeoShape::Tolerance()) return 0.0;
	Double_t invlen = 1./length;
	dir[0] = (end[0]-start[0])*invlen;
	dir[1] = (end[1]-start[1])*invlen;
	dir[2] = (end[2]-start[2])*invlen;

	// Initialize start point and direction
	TGeoNode *currentnode = 0;
	TGeoNode *startnode = gGeoManager->InitTrack(start, dir);
	if (!startnode) {
		//AliErrorClass(Form("start point out of geometry: x %f, y %f, z %f",
    //		  start[0],start[1],start[2]));
    		return 0.0;}

	TGeoMaterial *material = startnode->GetVolume()->GetMedium()->GetMaterial();
	lparam[0]   = material->GetDensity();
	if (lparam[0] > mparam[7]) mparam[7]=lparam[0];
	lparam[1]   = material->GetRadLen();
	lparam[2]   = material->GetA();
	lparam[3]   = material->GetZ();
	lparam[4]   = length;
	lparam[5]   = lparam[3]/lparam[2];
	lparam[6]   = material->GetIntLen();
	Double_t  n = lparam[0]/lparam[2];
	Double_t  sigma = 1./(n*lparam[6])/mbarn;
	if (sigma > mparam[9]) mparam[9]=sigma;
	if (material->IsMixture()) {
		TGeoMixture * mixture = (TGeoMixture*)material;
    		lparam[5] =0;
    		Double_t sum =0;
    		for (Int_t iel=0;iel<mixture->GetNelements();iel++){
      			sum  += mixture->GetWmixt()[iel];
      			lparam[5]+= mixture->GetZmixt()[iel]*mixture->GetWmixt()[iel]/mixture->GetAmixt()[iel];
    		}
    	lparam[5]/=sum;}

	// Locate next boundary within length without computing safety.
	// Propagate either with length (if no boundary found) or just cross boundary
	gGeoManager->FindNextBoundaryAndStep(length, kFALSE);
	Double_t step = 0.0; // Step made
	Double_t snext = gGeoManager->GetStep();
	// If no boundary within proposed length, return current density
	if (!gGeoManager->IsOnBoundary()) {
		mparam[0] = lparam[0];
		mparam[1] = lparam[4]/lparam[1];
		mparam[2] = lparam[2];
		mparam[3] = lparam[3];
		mparam[4] = lparam[4];
		mparam[8] = lparam[4]/lparam[6];
		return lparam[0];}

	// Try to cross the boundary and see what is next
	Int_t nzero = 0;
	while (length>TGeoShape::Tolerance()) {
		currentnode = gGeoManager->GetCurrentNode();
		if (snext<2.*TGeoShape::Tolerance()) nzero++;
		else nzero = 0;
		if (nzero>3) {
      			// This means navigation has problems on one boundary
      			// Try to cross by making a small step
      			//AliErrorClass("Cannot cross boundary\n");
			mparam[0] = bparam[0]/step;
			mparam[1] = bparam[1];
			mparam[2] = bparam[2]/step;
			mparam[3] = bparam[3]/step;
			mparam[5] = bparam[5]/step;
			mparam[8] = bparam[6];
			mparam[4] = step;
			mparam[0] = 0.;             // if crash of navigation take mean density 0
			mparam[1] = 1000000;        // and infinite rad length
			return bparam[0]/step;
    			}
    		mparam[6]+=1.;
    		step += snext;
    		bparam[1]    += snext/lparam[1];
    		bparam[2]    += snext*lparam[2];
    		bparam[3]    += snext*lparam[3];
		bparam[5]    += snext*lparam[5];
    		bparam[6]    += snext/lparam[6];
    		bparam[0]    += snext*lparam[0];

    		if (snext>=length) break;
    		if (!currentnode) break;
    		length -= snext;
    		material = currentnode->GetVolume()->GetMedium()->GetMaterial();
    		lparam[0] = material->GetDensity();
    		if (lparam[0] > mparam[7]) mparam[7]=lparam[0];
    		lparam[1]  = material->GetRadLen();
    		lparam[2]  = material->GetA();
    		lparam[3]  = material->GetZ();
    		lparam[5]   = lparam[3]/lparam[2];
    		lparam[6]   = material->GetIntLen();
    		n = lparam[0]/lparam[2];
   		sigma = 1./(n*lparam[6])/mbarn;
   		if (sigma > mparam[9]) mparam[9]=sigma;
    		if (material->IsMixture()) {
      			TGeoMixture * mixture = (TGeoMixture*)material;
      			lparam[5]=0;
      			Double_t sum =0;
      			for (Int_t iel=0;iel<mixture->GetNelements();iel++){
        			sum+= mixture->GetWmixt()[iel];
        			lparam[5]+= mixture->GetZmixt()[iel]*mixture->GetWmixt()[iel]/mixture->GetAmixt()[iel];
      			}
      			lparam[5]/=sum;
    		}
    		gGeoManager->FindNextBoundaryAndStep(length, kFALSE);
    		snext = gGeoManager->GetStep();
 	}//end of while
  	mparam[0] = bparam[0]/step;
  	mparam[1] = bparam[1];
  	mparam[2] = bparam[2]/step;
 	mparam[3] = bparam[3]/step;
  	mparam[5] = bparam[5]/step;
  	mparam[8] = bparam[6];
  	return bparam[0]/step;
}


// -------------------------------------------------------------------------
// -----   Passing the event   ---------------------------------------------
Bool_t DMGenerator::ReadEvent(FairPrimaryGenerator* cpg)
{
	//some start/end positions in z (emulsion to Tracker 1)
	Double_t start[3]={0.,0.,startZ};
	Double_t end[3]={0.,0.,endZ};
	
	//cout << "Enter DMGenerator " << endl;
    	if (fFirst){
     		Double_t bparam=0.;
      		Double_t mparam[10];
      		bparam=MeanMaterialBudget(start, end, mparam);
      		cout << "Info DMGenerator: MaterialBudget " << start[2] << " - "<< end[2] <<  endl;
      		cout << "Info DMGenerator: MaterialBudget " << bparam <<  endl;
      		cout << "Info DMGenerator: MaterialBudget 0 " << mparam[0] <<  endl;
		cout << "Info DMGenerator: MaterialBudget 1 " << mparam[1] <<  endl;
		cout << "Info DMGenerator: MaterialBudget 2 " << mparam[2] <<  endl;
		cout << "Info DMGenerator: MaterialBudget 3 " << mparam[3] <<  endl;
		cout << "Info DMGenerator: MaterialBudget 4 " << mparam[4] <<  endl;
		cout << "Info DMGenerator: MaterialBudget 5 " << mparam[5] <<  endl;
		cout << "Info DMGenerator: MaterialBudget 6 " << mparam[6] <<  endl;
		cout << "Info DMGenerator: MaterialBudget " << mparam[0]*mparam[4] <<  endl;

      		fFirst = kFALSE;}

	if(fn==fNevents)	fLogger->Warning(MESSAGE_ORIGIN, "End of input file. Rewind.");
	fTree->GetEntry(fn);
    	fn++;
	if(fn%100==0)		 cout << "Info DMGenerator: Dark Matter event-nr "<< fn << endl;

	// Incoming DM particle, get a random px,py
    	//cout << "Info GenieGenerator: neutrino " << neu << "p-in "<< pzv << " nf "<< nf << endl;
    	//cout << "Info GenieGenerator: ztarget " << ztarget << endl;
    	Double_t bparam=0.;
    	Double_t mparam[10];
    	Double_t txdm=0;
    	Double_t tydm=0;
    	
	txdm=pxdm/pzdm;
	tydm=pydm/pzdm;
	start[0]=(txdm)*(start[2]-ztarget);
	start[1]=(tydm)*(start[2]-ztarget);
	//cout << "Info DMGenerator: dark matter xyz-start " << start[0] << "-" << start[1] << "-" << start[2] << endl;
	end[0]=txdm*(end[2]-ztarget);
        end[1]=tydm*(end[2]-ztarget);
        //cout << "Info DMGenerator: dark matter xyz-end " << end[0] << "-" << end[1] << "-" << end[2] << endl;
        //get material density between these two points
	bparam=MeanMaterialBudget(start, end, mparam);


	//loop over trajectory between start and end to pick an interaction point
	Double_t prob2int = -1.;
	Double_t x;
	Double_t y;
	Double_t z;
	Int_t count=0;
	Bool_t insideCone = false;
	

	while (prob2int<gRandom->Uniform(0.,1.)) {
      		//place x,y,z uniform along path
      		z=gRandom->Uniform(start[2],end[2]);
      		x=txdm*(z-ztarget);
      		y=tydm*(z-ztarget);

      		if (mparam[6]<0.5){
        		//mparam is number of boundaries along path. mparam[6]=0.: uniform material budget along path, use present x,y,z
        		prob2int=2.;
      		}else{
        		//get local material at this point, to calculate probability that interaction is at this point.
        		TGeoNode *node = gGeoManager->FindNode(x,y,z);
        		TGeoMaterial *mat = 0;
        		if (node && !gGeoManager->IsOutside()) {
          			mat = node->GetVolume()->GetMaterial();
         			//cout << "Info DMGenerator: mat " <<  count << ", " << mat->GetName() << ", " << mat->GetDensity() << endl;
        			//density relative to Prob largest density along this trajectory, i.e. use rho(Pt)
         			prob2int= mat->GetDensity()/mparam[7];
         			if (prob2int>1.) cout << "***WARNING*** DMGenerator: prob2int > Maximum density????" << prob2int << " maxrho:" << mparam[7] << " material: " <<  mat->GetName() << endl;
         				count+=1;
        			}else{
          				prob2int=0.;
			}
     		}//end if-else mparam

   	 }//end of while

	//cout << "Info DMGenerator: prob2int " << prob2int << ", " << count << endl;
	//Add the track corresponding to dark matter particle
	Double_t zrelative=z-ztarget;
	Double_t tof=TMath::Sqrt(x*x+y*y+zrelative*zrelative)/2.99792458e+6;
	cpg->AddTrack(dmpdg,pxdm,pydm,pzdm,x,y,z,-1,false,Edm,tof,mparam[0]*mparam[4]);
	
	//Add the track corresponding to outgoing particle
	for(int i=0; i<n_2ry; ++i){
		cpg->AddTrack(pdg_2ry[i],px_2ry[i],py_2ry[i],pz_2ry[i],x,y,z,0,true,E_2ry[i],tof,mparam[0]*mparam[4]);
	}
	
	//cout << "Info DMGenerator Return from DMGenerator" << endl;

  return kTRUE;
}



