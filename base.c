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

fTheta = getTheta( fTheta, fFreqRef, fPwmT);

sector = int(( fTheta2 / (PI()/3)));

fSin = sin(fTheta);
fCos = cos(fTheta);

fVphpk = fVrmsRef * sqrt(2) / sqrt(3);

fVqeRef = fVphpk;
fVdeRef = 0.;

dqe2dqs(fVqeRef,fVdeRef,&fVqsRef,&fVdsRef,fSin,fCos);

dqs2abs(fVqsRef,fVdsRef,&fVaRef,&fVbRef,&fVcRef);

minmax(fVaRef,fVbRef,fVcRef,&fVaRefMM,&fVbRefMM,&fVcRefMM);

real2mod(fVaRefMM,fVbRefMM,fVcRefMM,&fPwmU,&fPwmV,&fPwmW,fVdc);
	
Bound(&fPwmU,0,1);
Bound(&fPwmV,0,1);
Bound(&fPwmW,0,1);

out[0] = fPwmU;
out[1] = fPwmV;
out[2] = fPwmW;
out[3] = shiftMM ;
out[4] = fVaRefMM;
out[5] = fVbRefMM;
out[6] = fVcRefMM;

