#include <iostream>

#include "ElectronFunctions.h"

// globals 
TFile *f_2011, *f_2012; 
TH2D  *sf_2011_central, *sf_2012_central; 
TH2D  *sf_2011_upper, *sf_2012_upper; 
TH2D  *sf_2011_lower, *sf_2012_lower; 

//---------------------------------------------------------------------- 
void initializeEfficiencyWeights(const char * fname2011, 
				 const char * fname2012) 
//---------------------------------------------------------------------- 
{
  f_2011 = new TFile(fname2011);
  if( !(f_2011->IsOpen()) ) { 
    std::cerr << "error: bad 2011 filename, quitting ..." << std::endl;
    exit(1);
  }
  sf_2011_central = (TH2D*)(f_2011->Get("sfEtaPt"));
  sf_2011_upper   = (TH2D*)(f_2011->Get("sfEtaPtErrHi"));
  sf_2011_lower   = (TH2D*)(f_2011->Get("sfEtaPtErrLo"));

  f_2012 = new TFile(fname2012);
  if( !(f_2012->IsOpen()) ) { 
    std::cerr << "error: bad 2012 filename, quitting ..." << std::endl;
    exit(1);
  }
  sf_2012_central = (TH2D*)(f_2012->Get("sfEtaPt"));	 
  sf_2012_upper   = (TH2D*)(f_2012->Get("sfEtaPtErrHi"));
  sf_2012_lower   = (TH2D*)(f_2012->Get("sfEtaPtErrLo"));
}

//---------------------------------------------------------------------- 
float getEfficiencyWeight( TLorentzVector * v, int type, int era ) 
//---------------------------------------------------------------------- 
{
  TH2D *h_central, *h_error;
  
  if( era == 2011 ) { 
    h_central = sf_2011_central;
    if( type == CorrectionType::LOWER )      h_error   = sf_2011_lower;
    else if( type == CorrectionType::UPPER ) h_error   = sf_2011_upper;
    else if ( type != CorrectionType::CENTRAL)  { 
      std::cerr << "error: bad SF type, quitting ... " << std::endl; 
      exit(1);
    }
  } 
  else if( era == 2012 ) { 
    h_central = sf_2012_central;
    if( type == CorrectionType::LOWER )      h_error = sf_2012_lower;
    else if( type == CorrectionType::UPPER ) h_error = sf_2012_upper;
    else if ( type != CorrectionType::CENTRAL)  { 
      std::cerr << "error: bad SF type, quitting ... " << std::endl; 
      exit(1);
    }
  } 
  else {
    std::cerr << "error: unknown era, quitting ..."  << std::endl;
  }

  float pt  = v->Pt();
  float eta = fabs(v->Eta());
  float sf = h_central->GetBinContent( h_central->FindBin(eta,pt) );
  if( type == CorrectionType::LOWER ) sf -= h_error->GetBinContent( h_error->FindBin(eta,pt) );
  if( type == CorrectionType::UPPER ) sf += h_error->GetBinContent( h_error->FindBin(eta,pt) );
  return sf;
}



//---------------------------------------------------------------------- 
std::pair<float,float> getScaleSigmaPair(int type, int era, int isEE)
//---------------------------------------------------------------------- 
{
  if( era == 2011 ) { 
    if( type == CorrectionType::CENTRAL ) { 
      if( isEE ) { return std::pair<float,float>(SCALE_2011_EE,SIGMA_2011_EE); }   
      else       { return std::pair<float,float>(SCALE_2011_EB,SIGMA_2011_EB); }   
    }
    if( type == CorrectionType::UPPER ) { 
      if( isEE ) { return std::pair<float,float>(SCALE_2011_EE+SCALE_UNCERT_2011_EE,
						 SIGMA_2011_EE+SIGMA_UNCERT_2011_EE); }   
      else       { return std::pair<float,float>(SCALE_2011_EB+SCALE_UNCERT_2011_EB,
						 SIGMA_2011_EB+SIGMA_UNCERT_2011_EB); }   
    }
    if( type == CorrectionType::LOWER ) { 
      if( isEE ) { return std::pair<float,float>(SCALE_2011_EE-SCALE_UNCERT_2011_EE,
						 SIGMA_2011_EE-SIGMA_UNCERT_2011_EE); }   
      else       { return std::pair<float,float>(SCALE_2011_EB-SCALE_UNCERT_2011_EB,
						 SIGMA_2011_EB-SIGMA_UNCERT_2011_EB); }   
    } else { 
      std::cerr << "error: unknown correction type, quitting ... " << std::endl;
    }
  } // 2011 

  else if( era == 2012 ) { 
    if( type == CorrectionType::CENTRAL ) { 
      if( isEE ) { return std::pair<float,float>(SCALE_2012_EE,SIGMA_2012_EE); }   
      else       { return std::pair<float,float>(SCALE_2012_EB,SIGMA_2012_EB); }   
    }
    if( type == CorrectionType::UPPER ) { 
      if( isEE ) { return std::pair<float,float>(SCALE_2012_EE+SCALE_UNCERT_2012_EE,
						 SIGMA_2012_EE+SIGMA_UNCERT_2012_EE); }   
      else       { return std::pair<float,float>(SCALE_2012_EB+SCALE_UNCERT_2012_EB,
						 SIGMA_2012_EB+SIGMA_UNCERT_2012_EB); }   
    }
    if( type == CorrectionType::LOWER ) { 
      if( isEE ) { return std::pair<float,float>(SCALE_2012_EE-SCALE_UNCERT_2012_EE,
						 SIGMA_2012_EE-SIGMA_UNCERT_2012_EE); }   
      else       { return std::pair<float,float>(SCALE_2012_EB-SCALE_UNCERT_2012_EB,
						 SIGMA_2012_EB-SIGMA_UNCERT_2012_EB); }   
    } else { 
      std::cerr << "error: unknown correction type, quitting ... " << std::endl;
    }
  } // 2012 

  else { 
    std::cerr << "error: unknown era, quitting ..."  << std::endl;
  }
}


//---------------------------------------------------------------------- 
float getCorrectedPt( TLorentzVector * v, int type, int era, TRandom3 * r) 
//---------------------------------------------------------------------- 
{
  bool isEE = ( fabs(v->Eta()) >= 1.566 ); 
  std::pair<float,float> scale_sigma_pair = getScaleSigmaPair(type,era,isEE);
  float newpt = r->Gaus(scale_sigma_pair.first*v->Pt(),
		       scale_sigma_pair.second*v->Pt());
  return newpt;
}
