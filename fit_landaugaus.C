
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
