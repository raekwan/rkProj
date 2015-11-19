 g_nStepCount++;

// In case of error, uncomment next two lines. Set *pnError to 1 and copy Error message to szErrorMsg
 //*pnError=1;
 //strcpy(szErrorMsg, "Place Error description here.");

fFreqRef = in[0];
fVrmsRef = in[1]; 
fVdc  =  in[2];
fIas = in[3];
fIbs = in[4];
fIcs = in[5];

abs2dqs(fIas, fIbs, fIcs, &fIqs, &fIds);
dqs2dqe(fIqs,fIds,&fIqe,&fIde,fSin,fCos);

if(t>0.1)
{

	fIs = sqrt(fIqe*fIqe+fIde*fIde);

	fIdeEst = fIs;
	fIqeEst = sqrt(fIs*fIs-fIdeEst*fIdeEst);

	fDV = 100.;

	while ( fDV > .1 )
	{
		fVqeCalc = fRs * fIqeEst + 2*PI()*fFreqRef *fLm*fIdeEst;
		fVdeCalc = fRs * fIdeEst -  2*PI()*fFreqRef *fsigLs*fIqeEst;
		fVsCalc = sqrt(fVqeCalc*fVqeCalc+fVdeCalc*fVdeCalc) * sqrt(3) / sqrt(2);
		fDV = fVsCalc - fVrmsRef ;
		fIdeEst = fIdeEst - 1.0;
		fIqeEst = sqrt(fIs*fIs-fIdeEst*fIdeEst);
	}
}

CalcTeEst();

fTheta = getTheta( fTheta, fFreqRef, fPwmT);

sector = int(( fTheta2 / (PI()/3)));


fSin = sin(fTheta);
fCos = cos(fTheta);

fVphpk = fVrmsRef * sqrt(2) / sqrt(3);

fVqeRef = fVphpk;
fVdeRef = 0.;

dqe2dqs(fVqeRef,fVdeRef,&fVqsRef,&fVdsRef,fSin,fCos);

dqs2abs(fVqsRef,fVdsRef,&fVaRef,&fVbRef,&fVcRef);

real2mod(fVaRef,fVbRef,fVcRef,&fPwmU,&fPwmV,&fPwmW,fVdc);

minmax(fPwmU,fPwmV,fPwmW,&fPwmUmm,&fPwmVmm,&fPwmWmm);
getShft2ph(fPwmUmm,fPwmVmm,fPwmWmm);

switch(sector)
{
	case 0:
//		shift = fPwmWmm;
		shift =  ( fPwmUmm -1 );
	break;
	case 1:
		shift = fPwmWmm;
	break;
	case 2:
//		shift = fPwmUmm;
		shift = fPwmVmm -1 ;
	break;
	case 3:
		shift = fPwmUmm;
	break;
	case 4:
//		shift = fPwmVmm;
		shift = fPwmWmm -1 ;
	break;
	case 5:
		shift =fPwmVmm;
	break;
}
/*if(t<0.3)
{
	fPwmUmm -= shift2ph;
	fPwmVmm -= shift2ph;
	fPwmWmm -= shift2ph;
}
else
{
	fPwmUmm -= shiftMM;
	fPwmVmm -= shiftMM;
	fPwmWmm -= shiftMM;
}*/

	fPwmUmm -= shift2ph;
	fPwmVmm -= shift2ph;
	fPwmWmm -= shift2ph;


if(fPwmUmm>1.0||fPwmVmm >1.0||fPwmWmm >1.0||
	fPwmUmm < 0.0 || fPwmVmm < 0.00 || fPwmWmm < 0.0 )
{
	limitFlag = 1;
}
else
{
	limitFlag = 0;
}
	
Bound(&fPwmU,0,1);
Bound(&fPwmV,0,1);
Bound(&fPwmW,0,1);


if(0)
{
	out[0] = fPwmUmm;
	out[1] = fPwmVmm;
	out[2] = fPwmWmm;
}
else
{
	out[0] = fPwmU;
	out[1] = fPwmV;
	out[2] = fPwmW;
}
out[3] = fTeCalc;
out[4] = fIde;
out[5] = fIqeEst;
out[6] = fIdeEst;
