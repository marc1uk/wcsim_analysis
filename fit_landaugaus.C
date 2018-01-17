=================
.L $ROOTSYS/tutorials/fit/langaus.C

Double_t fr[2];                   // fit range
fr[0]=0                           // fit lower limit
fr[1]=5                           // fit upper limit

Double_t sv[4], pllo[4], plhi[4]; // start values, lower limits, upper limits on fit params

sv[0]=0.16;                       // start val 1: width of landau component (suitable?)
sv[1]=2.;                         // start val 2: peak centre
sv[2]=ahist.Integral();           // start val 3: integral
sv[3]=0.2;                        // start val 3: width of gaussian component (suitable?)

pllo[0]=0.01; pllo[1]=0.01; pllo[2]=1.0; pllo[3]=0.01;    // lower limits of parameter values
plhi[0]=5.0; plhi[1]=10.0; plhi[2]=100.*ahist.Integral(); plhi[3]=10.0;  // upper limits of parameter values

Double_t fp[4], fpe[4];           // fit result returns, fit result errors
Double_t chisqr;                  // fit chi^2
Int_t    ndf;                     // fit NDF

TF1 *fitsnr = langaufit(&ahist,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);

Double_t SNRPeak, SNRFWHM;
langaupro(fp,SNRPeak,SNRFWHM);
fitsnr->Draw("lsame");
=======================

.L $ROOTSYS/tutorials/fit/langaus.C before compiling...

declare in the caller source code:
        TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);

then use by:
         // fit range, lower limit, upper limit
        std::vector<double> fitrange{-25,5};
        // starting fit vals: landau width, landau centre, integral, width of gaussian component
        std::vector<double> start_vals{1.,0.,likelihood->Integral(),20.};
        // lower limits on fit parameter ranges
        std::vector<double> param_lowerlims{0.,-100.,1.,0.01};
        // upper limits on fit parameter ranges
        std::vector<double> param_upperlims{10.,100.,100.*likelihood->Integral(),100.};
        // fit results
        std::vector<double> fit_vals(4), fit_val_errors(4);
        Double_t chisqr;                  // fit chi^2
        Int_t    ndf;                     // fit NDF

        //gROOT->ProcessLine(".L $ROOTSYS/tutorials/fit/langaus.C");
        TF1 *fitlangaus = langaufit(likelihood,fitrange.data(),start_vals.data(),param_lowerlims.data(),
                              param_upperlims.data(),fit_vals.data(),fit_val_errors.data(),&chisqr,&ndf);
        fitlangaus->Draw("lsame");

print resulting expression formula with:
TString formstring = fitlangaus->GetExpFormula();

===========================

// the actual 'fit function' is a c++ function as follows:
// note this means it's not possible to get the 'fit expression' - because that would be all the source code here!
Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}
