#include <Stdlib.h>
#include <String.h>

#define PI() 3.141592

 int g_nInputNodes=0;
 int g_nOutputNodes=0;

 int g_nStepCount=0;

int limitFlag = 0;
int wCnt;
double fPwmFreq =5000.;
double fPwmT = 0.;
double fDTheta = 0.;
double fTheta = 0.;
double fTheta2 = 0.;
double shiftMM = 0.;
double shift2ph = 0.;

double fIas = 0.;
double fIbs = 0.;
double fIcs = 0.;
double fIds = 0.;
double fIqs = 0.;
double fIqe = 0.;
double fIde = 0.;
double fIs = 0.;
double fTeCalc = 0.;

double fLamdsFilt = 0.;
double fLamdsOld = 0.;
double fLamqsFilt = 0.;
double fLamqsOld = 0.;
double fLamHpfTc = 0.001;
double fHpfConst =  0.;


double fFreqRef = 0.;
double fVrmsRef = 0.;

double fSin = 0.;
double fCos = 0.;

double fVphpk = 0.;

double fVqsRef = 0.;
double fVdsRef = 0.;

double fVqeRef = 0.;
double fVdeRef = 0.;

double fVaRef = 0.;
double fVbRef = 0.;
double fVcRef = 0.;
double fVaRefMM = 0.;
double fVbRefMM = 0.;
double fVcRefMM = 0.;


double fPwmU = 0.;
double fPwmV = 0.;
double fPwmW = 0.;

double fPwmUReg = 0.;
double fPwmVReg = 0.;
double fPwmWReg = 0.;

double fPwmUmm = 0.;
double fPwmVmm = 0.;
double fPwmWmm = 0.;

double fRs  =0.294;
double fLls = 0.00139;
double fLm = 0.041;
double fLs = 0.04239;
double fsigLs = 0.002734;

double fLamds = 0.;
double fLamqs = 0.;

double fVdc = 0.;
double shift = 0.;
int sector;

double fVqeCalc;
double fVdeCalc;
double fDV;
double fVsCalc;
double fIqeEst;
double fIdeEst;

void dqe2dqs(double qe, double de, double *qs, double *ds, double sin, double cos)
{
	*qs = sin * de + cos * qe;
	*ds = cos * de - sin * qe;
}

void dqs2dqe(double qs, double ds, double *qe, double *de, double sin, double cos)
{
	*qe = - sin * ds + cos * qs;
	*de = cos * ds + sin * qs;
}


void dqs2abs(double qs, double ds, double *as, double *bs, double *cs)
{
	*as =  qs;
	*bs =( -0.5 * qs - sqrt(3.)/2. * ds );
	*cs =( -0.5 * qs + sqrt(3.)/2. * ds );
}

void abs2dqs(double as, double bs, double cs, double *qs, double *ds)
{
	*qs =  as;
	*ds = 1./sqrt(3.) * (cs-bs);
}

void real2mod(double fVa,double fVb,double fVc,double *PwmU,double *PwmV,double *PwmW,double Vdc)
{

	*PwmU = 0.5 + 0.5 * fVa  / (0.5*Vdc);
	*PwmV = 0.5 + 0.5 * fVb  / (0.5*Vdc);
	*PwmW = 0.5 + 0.5 * fVc  / (0.5*Vdc);
}

void minmax(double phU,double phV,double phW,double  *phUmm,double  *phVmm,double  *phWmm)
{
	double min,max;
	
	if(phU>phV)
	{
		max = phU;
		min = phV;
	}
	else
	{
		max = phV;
		min = phU;
	}
	if(phW>max)
	{
		max = phW;
	}
	if(phW<min)
	{
		min = phW;
	}

	shiftMM  =  (min+max)*0.5;
//	shiftMM  =  (min+max)*0.5 - 0.5;

	*phUmm = phU  -shiftMM;
	*phVmm = phV - shiftMM ; 
	*phWmm = phW - shiftMM ;
//	*phUmm = phU;
//	*phVmm = phV;
//	*phWmm = phW;
}

void getShft2ph(double phU,double phV,double phW)
{
	if(phU>1.0)
	{
		shift2ph =  phU - 1.;
	}
	if(phV>1.0)
	{
		shift2ph =  phV - 1.;
	}
	if(phW>1.0)
	{
		shift2ph =  phW - 1.;
	}

	if(phU<0.)
	{
		shift2ph =  phU;
	}
	if(phV<0.)
	{
		shift2ph =  phV;
	}
	if(phW<0.)
	{
		shift2ph = phW;
	}

}

void fHpf(double *out, double in, double *prev_in, double a)
{
	*out = a * (*out) + a * (in-*prev_in);
	*prev_in = in;
}

void Bound(double *in, double low, double hi)
{
	if(*in>hi)
		*in = hi;
	if(*in<low)
		*in = low;

}
double GetHpfConst(double Tc, double Tsamp)
{
	return Tc/(Tc+Tsamp);
}

void CalcTeEst(void)
{
	
	fLamds = fLamds + ( fVdsRef - fRs * fIds ) * fPwmT;
	fLamqs = fLamqs + ( fVqsRef - fRs * fIqs ) * fPwmT;

	fHpf( &fLamdsFilt , fLamds , &fLamdsOld, fHpfConst );
	fHpf( &fLamqsFilt , fLamqs , &fLamqsOld, fHpfConst );

	fTeCalc = 3./2. * 6. / 2. * (fLamdsFilt * fIqs - fLamqsFilt * fIds);

}

double getTheta(double ThetaCur, double freq,double Tsamp)
{
	double DTheta,Theta;

	Theta = ThetaCur;

	DTheta = freq * 2. * PI() * Tsamp;

	Theta = ThetaCur+ DTheta;

	if( Theta>=2.*PI()) Theta -= (2.*PI());
	if( Theta<= 0 ) Theta += (2.*PI());

	return Theta;
}