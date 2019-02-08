#include "Framework/Conventions/Constants.h"
#include "Framework/Utils/cohrhosigxy.h"

using namespace std;
using namespace genie::constants;
double genie::utils::cohrhosigxy::cohrhosigxy_cal(double ENU, double UU, double ELEP, double QQ2, double TMIN, double TMAX){


int ICOHPIN=2;
double RFAC=0;
double RFLG1=1;
double RFLG2=0;
//YY=1.;
double ATM=12;
double ATM2=ATM*ATM;
double RA=R0*TMath::Power(ATM, 1.0/3.0);
double AMT=11.17793;
//cout<<"libo test in cohrhosig! "<<endl;

double BSLOPE=(1.0/3.0)*pow((RA*FMTOGEV), 2.0);
double BTNUM=exp(-BSLOPE*TMIN)-exp(-BSLOPE*TMAX);
//cout<<"R0="<<R0<<endl;
//cout<<"RA="<<RA<<endl;
//cout<<"BSLOPE= "<<BSLOPE<<" "<<"BTNUM= "<<BTNUM<<endl;
//------------------------------------------------------
double ENU2=ENU*ENU;
double UU2=UU*UU;
//cout<<"in cohrhosigxy: "<<"ENU2="<<ENU2<<" "<<"UU2="<<UU2<<endl;
double EPSILON=(4.0*ENU*ELEP-QQ2)/(4.0*ENU*ELEP+QQ2+2.0*UU2);


double PROPRHO=AMRHO2/(QQ2+AMRHO2);
//double SIGOUT;

double GRHO2=2.4*4.0*kPi;
double RHOSIG=0.0;
if(ICOHPIN==1){
    RHOSIG=24.0;
}
else if(ICOHPIN==2){
RHOSIG=24.0*(1.0+0.5/sqrt(UU));
}
//double SIGOUT;

//cout<<"EPSILON="<<EPSILON<<endl;
double R=RFAC*(RFLG1+RFLG2*QQ2/AMRHO2);
//cout<<"RFAC="<<RFAC<<" "<<RFLG1<<" "<<RFLG2<<endl;
double FABS=0;
if(R>1.0) {R=1.0;}
if(ATM<1.1){
FABS=1.0;
}
else{
FABS=exp(-9.0/16.0/kPi*RA/(pow(R0, 3.0))/FM2TOMB*19.0*(1.0+0.44891/sqrt(UU)));
}
double Q = TMath::Sqrt(UU*UU+QQ2);

double SIGOUT = GFERMI2/(2.0*kPi2)*
        1.0/GRHO2*
        PROPRHO*PROPRHO*
        Q/ENU2*
        QQ2/(1.0 - EPSILON)*
        (1.0 + EPSILON*R)*
        ATM2/(16.0*kPi)*
        (RHOSIG*RHOSIG)*
        FABS*
        BTNUM/BSLOPE*
        MBTOCM2/GEVTOMB/MBTOCM2/GEVTOMB;

			  
double YY=UU/ENU;
SIGOUT=2.0*AMT*YY*ENU2*SIGOUT;
//cout<<"MBTOCM2="<<MBTOCM2<<" "<<"GEVTOMB="<<GEVTOMB<<endl;
//cout<<"GFERMI2="<<GFERMI2<<" "<<"GRHO2="<<GRHO2<<" "<<PROPRHO<<" "<<Q<<" "<<ENU2<<" "<<RHOSIG<<" "<<FABS<<endl;
//cout<<"SIGOUT="<<SIGOUT<<""<<"YY="<<YY<<endl;


return SIGOUT;       
//------------------------------------------------------------
}
