
const double GARHO=0.1507;
const double GARHO2=GARHO*GARHO;

using namespace std;
double fgrhom(double RHOM){

//get input variables-----------------------
//double CMSMOM;
double FGRHOM;
double Q1, Q0, GARHOP;
const double AMPI=0.13957;
const double AMRHO=0.7700;

//include breit-wigner function for coherent rho
double breitw(double, double, double);


double cmsmom(double, double, double);


Q1=cmsmom(RHOM, AMPI, AMPI);
Q0=cmsmom(AMRHO, AMPI, AMPI);
GARHOP=GARHO*pow((Q1/Q0),3.0)*AMRHO/RHOM;


FGRHOM=breitw(RHOM, AMRHO, GARHOP);
FGRHOM=FGRHOM*RHOM/Q1;
//cout<<"libo test in fgrhom! "<<endl;
return FGRHOM;

}


//calculate breit wigner function
//XM = Mass at which to calculate function
//XM0 Central mass of resonance
//GA Width of resonance
double breitw(double XM, double XM0, double GA){
double BREITW;
BREITW=GA/(pow((XM0*XM0-XM*XM),2.0)+pow((XM0*GA),2.0));
return BREITW;
}

// calculate the cms momentum of the rho decay system
// XM0: the mass of parent
// XM1: mass of daughter for which momentum is required
// XM2: mass of the other daughter
double cmsmom(double XM0, double XM1, double XM2){
double CMSMOM;
CMSMOM=sqrt((XM0*XM0-(XM1+XM2)*(XM1+XM2))*(XM0*XM0-(XM1-XM2)*(XM1-XM2)));
return CMSMOM;
}
