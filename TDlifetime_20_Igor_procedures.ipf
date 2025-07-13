#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Global Fit 2>
#include <SaveGraph>

Function Decay1_func (w, yw, xw) : FitFunc		//Calculate the decay using parameter wave w
	wave w, yw, xw
	wave Decay = root:Decay1 	//Declare all needed global variables
	wave fit_Decay=root:fit_decay1
	wave IRF = root:IRF1
	wave Timex=root:Time1
	wave Coef_DecayB=root:Coef_DecayB1
	wave Res_Decay=root:Res_Decay1
	wave pr = root:pr
	wave tda1 = root:tda1
	wave tda2 = root:tda2
	wave Cursors=root:Cursors
	
	Duplicate/FREE/O/R=[Cursors[2]/Timex[1],Cursors[3]/Timex[1]] IRF IRFx
	SetScale/P x 0,Timex[1],"s", IRFx 	//Prepare IRF data for convolution
	Duplicate/O Decay Res_Decay, fit_Decay
	SetScale/P x 0,deltax(IRFx),"s", fit_Decay
	Make/D/FREE/O/N=(numpnts(Decay)/10+1) fit_Decay_short
	SetScale/P x 0,deltax(IRFx)*10,"s", fit_Decay_short
	Make/D/FREE/O/N=(numpnts(IRFx)+w[10]/deltax(IRFx)) IRF_shifted
	SetScale/P x 0,deltax(IRFx),"s", IRF_shifted
	IRF_shifted=0

		
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = **
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 16
	//CurveFitDialog/ w[0] = frac_D
	//CurveFitDialog/ w[1] = tauD1
	//CurveFitDialog/ w[2] = alphaD1
	//CurveFitDialog/ w[3] = tauD2
	//CurveFitDialog/ w[4] = r0
	//CurveFitDialog/ w[5] = avg_r1
	//CurveFitDialog/ w[6] = sig_r1
	//CurveFitDialog/ w[7] = frac_r2
	//CurveFitDialog/ w[8] = avg_r2
	//CurveFitDialog/ w[9] = sig_r2
	//CurveFitDialog/ w[10] = shift_irf
	//CurveFitDialog/ w[11] = frac_buf
	//CurveFitDialog/ w[12] = a0
	//CurveFitDialog/ w[13] = bkgr_dec
	//CurveFitDialog/ w[14] = bkgr_IRF
	//CurveFitDialog/ w[15] = V_chisq

	
	pr=((1-w[7])*gauss(x, w[5], w[6])+w[7]*gauss(x, w[8], w[9]))	//Double Gaussian distance distribution
	tda1=1/(1/w[1]+1/w[1]*(w[4]/x)^6)
	tda2=1/(1/w[3]+1/w[3]*(w[4]/x)^6)
	
	fit_Decay_short=w[12]*(w[0]*(w[2]*exp(-x/w[1])+(1-w[2])*exp(-x/w[3]))+(1-w[0])*it_func (w[2],x))
	fit_Decay=fit_Decay_short(x)+w[11]*(Coef_DecayB[0]*exp(-x/Coef_DecayB[1])+Coef_DecayB[2]*exp(-x/Coef_DecayB[3])+Coef_DecayB[4]*exp(-x/Coef_DecayB[5])+Coef_DecayB[6]*exp(-x/Coef_DecayB[7])) //Add buffer-only fit to dacay
	
	IRF_shifted[w[10]/deltax(IRFx), ] = IRFx[max(p-w[10]/deltax(IRFx),0)]-w[14]
	K0=sum(IRF_shifted)
	IRF_shifted=IRF_shifted/K0
	convolve IRF_shifted, fit_Decay
	Redimension/N=(numpnts(Decay)) fit_Decay
	Fit_Decay=Fit_Decay+w[13]
	Res_Decay=Decay-Fit_Decay
	yw = fit_Decay(xw[p])
End

Function Decay2_func (w, yw, xw) : FitFunc		//Calculate the decay using parameter wave w
	wave w, yw, xw
	wave Decay = root:Decay2 	//Declare all needed global variables
	wave fit_Decay=root:fit_decay2
	wave IRF = root:IRF2
	wave Timex=root:Time2
	wave Coef_DecayB=root:Coef_DecayB2
	wave Res_Decay=root:Res_Decay2
	wave pr = root:pr
	wave tda1 = root:tda1
	wave tda2 = root:tda2
	wave Cursors=root:Cursors
	
	Duplicate/FREE/O/R=[Cursors[2]/Timex[1],Cursors[3]/Timex[1]] IRF IRFx
	SetScale/P x 0,Timex[1],"s", IRFx 	//Prepare IRF data for convolution
	Duplicate/O Decay Res_Decay, fit_Decay
	SetScale/P x 0,deltax(IRFx),"s", fit_Decay
	Make/D/FREE/O/N=(numpnts(Decay)/10+1) fit_Decay_short
	SetScale/P x 0,deltax(IRFx)*10,"s", fit_Decay_short
	Make/D/FREE/O/N=(numpnts(IRFx)+w[10]/deltax(IRFx)) IRF_shifted
	SetScale/P x 0,deltax(IRFx),"s", IRF_shifted
	IRF_shifted=0

		
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = **
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 16
	//CurveFitDialog/ w[0] = frac_D
	//CurveFitDialog/ w[1] = tauD1
	//CurveFitDialog/ w[2] = alphaD1
	//CurveFitDialog/ w[3] = tauD2
	//CurveFitDialog/ w[4] = r0
	//CurveFitDialog/ w[5] = avg_r1
	//CurveFitDialog/ w[6] = sig_r1
	//CurveFitDialog/ w[7] = frac_r2
	//CurveFitDialog/ w[8] = avg_r2
	//CurveFitDialog/ w[9] = sig_r2
	//CurveFitDialog/ w[10] = shift_irf
	//CurveFitDialog/ w[11] = frac_buf
	//CurveFitDialog/ w[12] = a0
	//CurveFitDialog/ w[13] = bkgr_dec
	//CurveFitDialog/ w[14] = bkgr_IRF
	//CurveFitDialog/ w[15] = V_chisq

	
	pr=((1-w[7])*gauss(x, w[5], w[6])+w[7]*gauss(x, w[8], w[9]))	//Double Gaussian distance distribution
	tda1=1/(1/w[1]+1/w[1]*(w[4]/x)^6)
	tda2=1/(1/w[3]+1/w[3]*(w[4]/x)^6)
	
	fit_Decay_short=w[12]*(w[0]*(w[2]*exp(-x/w[1])+(1-w[2])*exp(-x/w[3]))+(1-w[0])*it_func (w[2],x))
	fit_Decay=fit_Decay_short(x)+w[11]*(Coef_DecayB[0]*exp(-x/Coef_DecayB[1])+Coef_DecayB[2]*exp(-x/Coef_DecayB[3])+Coef_DecayB[4]*exp(-x/Coef_DecayB[5])+Coef_DecayB[6]*exp(-x/Coef_DecayB[7])) //Add buffer-only fit to dacay
	
	IRF_shifted[w[10]/deltax(IRFx), ] = IRFx[max(p-w[10]/deltax(IRFx),0)]-w[14]
	K0=sum(IRF_shifted)
	IRF_shifted=IRF_shifted/K0
	convolve IRF_shifted, fit_Decay
	Redimension/N=(numpnts(Decay)) fit_Decay
	Fit_Decay=Fit_Decay+w[13]
	Res_Decay=Decay-Fit_Decay
	yw = fit_Decay(xw[p])
End

Function Decay3_func (w, yw, xw) : FitFunc		//Calculate the decay using parameter wave w
	wave w, yw, xw
	wave Decay = root:Decay3 	//Declare all needed global variables
	wave fit_Decay=root:fit_decay3
	wave IRF = root:IRF3
	wave Timex=root:Time3
	wave Coef_DecayB=root:Coef_DecayB3
	wave Res_Decay=root:Res_Decay3
	wave pr = root:pr
	wave tda1 = root:tda1
	wave tda2 = root:tda2
	wave Cursors=root:Cursors
	
	Duplicate/FREE/O/R=[Cursors[2]/Timex[1],Cursors[3]/Timex[1]] IRF IRFx
	SetScale/P x 0,Timex[1],"s", IRFx 	//Prepare IRF data for convolution
	Duplicate/O Decay Res_Decay, fit_Decay
	SetScale/P x 0,deltax(IRFx),"s", fit_Decay
	Make/D/FREE/O/N=(numpnts(Decay)/10+1) fit_Decay_short
	SetScale/P x 0,deltax(IRFx)*10,"s", fit_Decay_short
	Make/D/FREE/O/N=(numpnts(IRFx)+w[10]/deltax(IRFx)) IRF_shifted
	SetScale/P x 0,deltax(IRFx),"s", IRF_shifted
	IRF_shifted=0

		
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = **
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 16
	//CurveFitDialog/ w[0] = frac_D
	//CurveFitDialog/ w[1] = tauD1
	//CurveFitDialog/ w[2] = alphaD1
	//CurveFitDialog/ w[3] = tauD2
	//CurveFitDialog/ w[4] = r0
	//CurveFitDialog/ w[5] = avg_r1
	//CurveFitDialog/ w[6] = sig_r1
	//CurveFitDialog/ w[7] = frac_r2
	//CurveFitDialog/ w[8] = avg_r2
	//CurveFitDialog/ w[9] = sig_r2
	//CurveFitDialog/ w[10] = shift_irf
	//CurveFitDialog/ w[11] = frac_buf
	//CurveFitDialog/ w[12] = a0
	//CurveFitDialog/ w[13] = bkgr_dec
	//CurveFitDialog/ w[14] = bkgr_IRF
	//CurveFitDialog/ w[15] = V_chisq

	
	pr=((1-w[7])*gauss(x, w[5], w[6])+w[7]*gauss(x, w[8], w[9]))	//Double Gaussian distance distribution
	tda1=1/(1/w[1]+1/w[1]*(w[4]/x)^6)
	tda2=1/(1/w[3]+1/w[3]*(w[4]/x)^6)
	
	fit_Decay_short=w[12]*(w[0]*(w[2]*exp(-x/w[1])+(1-w[2])*exp(-x/w[3]))+(1-w[0])*it_func (w[2],x))
	fit_Decay=fit_Decay_short(x)+w[11]*(Coef_DecayB[0]*exp(-x/Coef_DecayB[1])+Coef_DecayB[2]*exp(-x/Coef_DecayB[3])+Coef_DecayB[4]*exp(-x/Coef_DecayB[5])+Coef_DecayB[6]*exp(-x/Coef_DecayB[7])) //Add buffer-only fit to dacay
	
	IRF_shifted[w[10]/deltax(IRFx), ] = IRFx[max(p-w[10]/deltax(IRFx),0)]-w[14]
	K0=sum(IRF_shifted)
	IRF_shifted=IRF_shifted/K0
	convolve IRF_shifted, fit_Decay
	Redimension/N=(numpnts(Decay)) fit_Decay
	Fit_Decay=Fit_Decay+w[13]
	Res_Decay=Decay-Fit_Decay
	yw = fit_Decay(xw[p])
End

Function Decay4_func (w, yw, xw) : FitFunc		//Calculate the decay using parameter wave w
	wave w, yw, xw
	wave Decay = root:Decay4 	//Declare all needed global variables
	wave fit_Decay=root:fit_decay4
	wave IRF = root:IRF4
	wave Timex=root:Time4
	wave Coef_DecayB=root:Coef_DecayB4
	wave Res_Decay=root:Res_Decay4
	wave pr = root:pr
	wave tda1 = root:tda1
	wave tda2 = root:tda2
	wave Cursors=root:Cursors
	
	Duplicate/FREE/O/R=[Cursors[2]/Timex[1],Cursors[3]/Timex[1]] IRF IRFx
	SetScale/P x 0,Timex[1],"s", IRFx 	//Prepare IRF data for convolution
	Duplicate/O Decay Res_Decay, fit_Decay
	SetScale/P x 0,deltax(IRFx),"s", fit_Decay
	Make/D/FREE/O/N=(numpnts(Decay)/10+1) fit_Decay_short
	SetScale/P x 0,deltax(IRFx)*10,"s", fit_Decay_short
	Make/D/FREE/O/N=(numpnts(IRFx)+w[10]/deltax(IRFx)) IRF_shifted
	SetScale/P x 0,deltax(IRFx),"s", IRF_shifted
	IRF_shifted=0

		
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = **
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 16
	//CurveFitDialog/ w[0] = frac_D
	//CurveFitDialog/ w[1] = tauD1
	//CurveFitDialog/ w[2] = alphaD1
	//CurveFitDialog/ w[3] = tauD2
	//CurveFitDialog/ w[4] = r0
	//CurveFitDialog/ w[5] = avg_r1
	//CurveFitDialog/ w[6] = sig_r1
	//CurveFitDialog/ w[7] = frac_r2
	//CurveFitDialog/ w[8] = avg_r2
	//CurveFitDialog/ w[9] = sig_r2
	//CurveFitDialog/ w[10] = shift_irf
	//CurveFitDialog/ w[11] = frac_buf
	//CurveFitDialog/ w[12] = a0
	//CurveFitDialog/ w[13] = bkgr_dec
	//CurveFitDialog/ w[14] = bkgr_IRF
	//CurveFitDialog/ w[15] = V_chisq

	
	pr=((1-w[7])*gauss(x, w[5], w[6])+w[7]*gauss(x, w[8], w[9]))	//Double Gaussian distance distribution
	tda1=1/(1/w[1]+1/w[1]*(w[4]/x)^6)
	tda2=1/(1/w[3]+1/w[3]*(w[4]/x)^6)
	
	fit_Decay_short=w[12]*(w[0]*(w[2]*exp(-x/w[1])+(1-w[2])*exp(-x/w[3]))+(1-w[0])*it_func (w[2],x))
	fit_Decay=fit_Decay_short(x)+w[11]*(Coef_DecayB[0]*exp(-x/Coef_DecayB[1])+Coef_DecayB[2]*exp(-x/Coef_DecayB[3])+Coef_DecayB[4]*exp(-x/Coef_DecayB[5])+Coef_DecayB[6]*exp(-x/Coef_DecayB[7])) //Add buffer-only fit to dacay
	
	IRF_shifted[w[10]/deltax(IRFx), ] = IRFx[max(p-w[10]/deltax(IRFx),0)]-w[14]
	K0=sum(IRF_shifted)
	IRF_shifted=IRF_shifted/K0
	convolve IRF_shifted, fit_Decay
	Redimension/N=(numpnts(Decay)) fit_Decay
	Fit_Decay=Fit_Decay+w[13]
	Res_Decay=Decay-Fit_Decay
	yw = fit_Decay(xw[p])
End

Function it_func (alphaD1,t)		//Calculate the intensity for time t using parameter wave w
	variable alphaD1, t
	wave pr = root:pr 			//Declare all needed global variables
	wave tda1 = root:tda1
	wave tda2 = root:tda2
	wave itr = root:itr
	itr=pr(x)*(alphaD1*exp(-t/tda1(x))+(1-alphaD1)*exp(-t/tda2(x)))
	return area(itr,0,inf)
End

Function jDA (w) : FitFunc		//Calculate normalizaton factor J for the donor/acceptor sample using parameter wave w
	wave w
	wave pr = root:pr		//Declare all needed global variables
	wave tda1 = root:tda1
	wave tda2 = root:tda2
	wave dfr = root:dfr
	
	pr=((1-w[7])*gauss(x, w[5], w[6])+w[7]*gauss(x, w[8], w[9]))	//Double Gaussian distance distribution
	tda1=1/(1/w[1]+1/w[1]*(w[4]/x)^6)
	tda2=1/(1/w[3]+1/w[3]*(w[4]/x)^6)
	dfr=(1-w[0])*(pr(x)*w[2]*tda1(x)+pr(x)*(1-w[2])*tda2(x))
	return area(dfr,0,inf)+w[0]*(w[2]*w[1]+(1-w[2])*w[3])	//J is just the unnormalized D at freq=0
End

Function jD (w) : FitFunc		//Calculate normalizaton factor J for the donor-only sample using parameter wave w
	wave w
	
	return w[2]*w[1]+(1-w[2])*w[3]
End

Function eff_func (w,c) : FitFunc	//Calculate the FRET efficiency as measured from the intensity using parameter wave w
	wave w; variable c	//need vaiable c to use as a fitting function

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = **
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 15
	//CurveFitDialog/ w[0] = frac_D
	//CurveFitDialog/ w[1] = tauD1
	//CurveFitDialog/ w[2] = alphaD1
	//CurveFitDialog/ w[3] = tauD2
	//CurveFitDialog/ w[4] = r0
	//CurveFitDialog/ w[5] = avg_r1
	//CurveFitDialog/ w[6] = sig_r1
	//CurveFitDialog/ w[7] = frac_r2
	//CurveFitDialog/ w[8] = avg_r2
	//CurveFitDialog/ w[9] = sig_r2
	//CurveFitDialog/ w[10] = shift_irf
	//CurveFitDialog/ w[11] = frac_buf
	//CurveFitDialog/ w[12] = a0
	//CurveFitDialog/ w[13] = bkgr_dec
	//CurveFitDialog/ w[14] = bkgr_irf
	
	return 1-jDA(w)/jD(w)
End

Function DecayB1_func (w, yw, xw) : FitFunc		//Fit the decay in buffer only sample using parameter wave w
	wave w, yw, xw
	wave DecayB = root:DecayB1 	//Declare all needed global variables
	wave fit_DecayB=root:fit_decayB1
	wave IRF = root:IRF1
	wave Timex=root:Time1
	wave Cursors=root:Cursors
	
	Duplicate/FREE/O/R=[Cursors[2]/Timex[1],Cursors[3]/Timex[1]] IRF IRFx
	SetScale/P x 0,Timex[1],"s", IRFx 	//Prepare IRF data for convolution
	Duplicate/O DecayB fit_DecayB			// Double-precision free waves
	SetScale/P x 0,deltax(IRFx),"s", fit_DecayB
	Make/D/FREE/O/N=(numpnts(IRFx)+w[8]/deltax(IRFx)) IRF_shifted
	SetScale/P x 0,deltax(IRFx),"s", IRF_shifted
	IRF_shifted=0
		
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = **
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 12
	//CurveFitDialog/ w[0] = A_1
	//CurveFitDialog/ w[1] = tau_1
	//CurveFitDialog/ w[2] = A_2
	//CurveFitDialog/ w[3] = tau_2
	//CurveFitDialog/ w[4] = A_3
	//CurveFitDialog/ w[5] = tau_3
	//CurveFitDialog/ w[6] = A_4
	//CurveFitDialog/ w[7] = tau_4
	//CurveFitDialog/ w[8] = shift_irf
	//CurveFitDialog/ w[9] = bkgr_dec
	//CurveFitDialog/ w[10] = bkgr_irf
	//CurveFitDialog/ w[11] = V_chisq
	
	fit_DecayB=w[0]*exp(-x/w[1])+w[2]*exp(-x/w[3])+w[4]*exp(-x/w[5])+w[6]*exp(-x/w[7])
	
	IRF_shifted[w[8]/deltax(IRFx), ] = IRFx[max(p-w[8]/deltax(IRFx),0)]-w[10]
	K0=sum(IRF_shifted)
	IRF_shifted=IRF_shifted/K0
	convolve IRF_shifted, fit_DecayB
	fit_DecayB=fit_DecayB+w[9]
	yw = fit_DecayB(xw[p])
End

Function DecayB2_func (w, yw, xw) : FitFunc		//Fit the decay in buffer only sample using parameter wave w
	wave w, yw, xw
	wave DecayB = root:DecayB2 	//Declare all needed global variables
	wave fit_DecayB=root:fit_decayB2
	wave IRF = root:IRF2
	wave Timex=root:Time2
	wave Cursors=root:Cursors
	
	Duplicate/FREE/O/R=[Cursors[2]/Timex[1],Cursors[3]/Timex[1]] IRF IRFx
	SetScale/P x 0,Timex[1],"s", IRFx 	//Prepare IRF data for convolution
	Duplicate/O DecayB fit_DecayB			// Double-precision free waves
	SetScale/P x 0,deltax(IRFx),"s", fit_DecayB
	Make/D/FREE/O/N=(numpnts(IRFx)+w[8]/deltax(IRFx)) IRF_shifted
	SetScale/P x 0,deltax(IRFx),"s", IRF_shifted
	IRF_shifted=0
		
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = **
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 12
	//CurveFitDialog/ w[0] = A_1
	//CurveFitDialog/ w[1] = tau_1
	//CurveFitDialog/ w[2] = A_2
	//CurveFitDialog/ w[3] = tau_2
	//CurveFitDialog/ w[4] = A_3
	//CurveFitDialog/ w[5] = tau_3
	//CurveFitDialog/ w[6] = A_4
	//CurveFitDialog/ w[7] = tau_4
	//CurveFitDialog/ w[8] = shift_irf
	//CurveFitDialog/ w[9] = bkgr_dec
	//CurveFitDialog/ w[10] = bkgr_irf
	//CurveFitDialog/ w[11] = V_chisq
	
	fit_DecayB=w[0]*exp(-x/w[1])+w[2]*exp(-x/w[3])+w[4]*exp(-x/w[5])+w[6]*exp(-x/w[7])
	
	IRF_shifted[w[8]/deltax(IRFx), ] = IRFx[max(p-w[8]/deltax(IRFx),0)]-w[10]
	K0=sum(IRF_shifted)
	IRF_shifted=IRF_shifted/K0
	convolve IRF_shifted, fit_DecayB
	fit_DecayB=fit_DecayB+w[9]
	yw = fit_DecayB(xw[p])
End

Function DecayB3_func (w, yw, xw) : FitFunc		//Fit the decay in buffer only sample using parameter wave w
	wave w, yw, xw
	wave DecayB = root:DecayB3 	//Declare all needed global variables
	wave fit_DecayB=root:fit_decayB3
	wave IRF = root:IRF3
	wave Timex=root:Time3
	wave Cursors=root:Cursors
	
	Duplicate/FREE/O/R=[Cursors[2]/Timex[1],Cursors[3]/Timex[1]] IRF IRFx
	SetScale/P x 0,Timex[1],"s", IRFx 	//Prepare IRF data for convolution
	Duplicate/O DecayB fit_DecayB			// Double-precision free waves
	SetScale/P x 0,deltax(IRFx),"s", fit_DecayB
	Make/D/FREE/O/N=(numpnts(IRFx)+w[8]/deltax(IRFx)) IRF_shifted
	SetScale/P x 0,deltax(IRFx),"s", IRF_shifted
	IRF_shifted=0
		
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = **
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 12
	//CurveFitDialog/ w[0] = A_1
	//CurveFitDialog/ w[1] = tau_1
	//CurveFitDialog/ w[2] = A_2
	//CurveFitDialog/ w[3] = tau_2
	//CurveFitDialog/ w[4] = A_3
	//CurveFitDialog/ w[5] = tau_3
	//CurveFitDialog/ w[6] = A_4
	//CurveFitDialog/ w[7] = tau_4
	//CurveFitDialog/ w[8] = shift_irf
	//CurveFitDialog/ w[9] = bkgr_dec
	//CurveFitDialog/ w[10] = bkgr_irf
	//CurveFitDialog/ w[11] = V_chisq
	
	fit_DecayB=w[0]*exp(-x/w[1])+w[2]*exp(-x/w[3])+w[4]*exp(-x/w[5])+w[6]*exp(-x/w[7])
	
	IRF_shifted[w[8]/deltax(IRFx), ] = IRFx[max(p-w[8]/deltax(IRFx),0)]-w[10]
	K0=sum(IRF_shifted)
	IRF_shifted=IRF_shifted/K0
	convolve IRF_shifted, fit_DecayB
	fit_DecayB=fit_DecayB+w[9]
	yw = fit_DecayB(xw[p])
End

Function DecayB4_func (w, yw, xw) : FitFunc		//Fit the decay in buffer only sample using parameter wave w
	wave w, yw, xw
	wave DecayB = root:DecayB4 	//Declare all needed global variables
	wave fit_DecayB=root:fit_decayB4
	wave IRF = root:IRF4
	wave Timex=root:Time4
	wave Cursors=root:Cursors
	
	Duplicate/FREE/O/R=[Cursors[2]/Timex[1],Cursors[3]/Timex[1]] IRF IRFx
	SetScale/P x 0,Timex[1],"s", IRFx 	//Prepare IRF data for convolution
	Duplicate/O DecayB fit_DecayB			// Double-precision free waves
	SetScale/P x 0,deltax(IRFx),"s", fit_DecayB
	Make/D/FREE/O/N=(numpnts(IRFx)+w[8]/deltax(IRFx)) IRF_shifted
	SetScale/P x 0,deltax(IRFx),"s", IRF_shifted
	IRF_shifted=0
		
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = **
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 12
	//CurveFitDialog/ w[0] = A_1
	//CurveFitDialog/ w[1] = tau_1
	//CurveFitDialog/ w[2] = A_2
	//CurveFitDialog/ w[3] = tau_2
	//CurveFitDialog/ w[4] = A_3
	//CurveFitDialog/ w[5] = tau_3
	//CurveFitDialog/ w[6] = A_4
	//CurveFitDialog/ w[7] = tau_4
	//CurveFitDialog/ w[8] = shift_irf
	//CurveFitDialog/ w[9] = bkgr_dec
	//CurveFitDialog/ w[10] = bkgr_irf
	//CurveFitDialog/ w[11] = V_chisq
	
	fit_DecayB=w[0]*exp(-x/w[1])+w[2]*exp(-x/w[3])+w[4]*exp(-x/w[5])+w[6]*exp(-x/w[7])
	
	IRF_shifted[w[8]/deltax(IRFx), ] = IRFx[max(p-w[8]/deltax(IRFx),0)]-w[10]
	K0=sum(IRF_shifted)
	IRF_shifted=IRF_shifted/K0
	convolve IRF_shifted, fit_DecayB
	fit_DecayB=fit_DecayB+w[9]
	yw = fit_DecayB(xw[p])
End

Function Fit_dataset_1 ()
	wave Decay = root:Decay1 	//Declare all needed global variables
	wave fit_Decay = root:fit_Decay1
	wave Timex=root:Time1
	wave Coef_Decay=root:Coef_Decay1
	wave SD_Decay=root:SD_decay1
	wave eps_Coef_Decay=root:eps_Coef_Decay
	wave Cursors=root:Cursors
	wave/T Fix=root:Fix1
	SVAR fixstr

	Duplicate/O Decay SD_Decay
	If (numpnts(fit_Decay)!=numpnts(Decay))
		Duplicate/O Decay fit_Decay
	Endif	
	fixstr=Fix[0]+Fix[1]+Fix[2]+Fix[3]+Fix[4]+Fix[5]+Fix[6]+Fix[7]+Fix[8]+Fix[9]+Fix[10]+Fix[11]+Fix[12]+Fix[13]+Fix[14]+"1"
	SD_Decay= sqrt(Fit_Decay)
	FuncFit/L=2000 /H=fixstr Decay1_func Coef_Decay Decay[Cursors[0]/Timex[1],Cursors[1]/Timex[1]] /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_Decay
	SD_Decay= sqrt(fit_Decay)
	FuncFit/L=2000 /H=fixstr Decay1_func Coef_Decay Decay[Cursors[0]/Timex[1],Cursors[1]/Timex[1]] /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_Decay
	SD_Decay= sqrt(fit_Decay)
	FuncFit/L=2000 /H=fixstr Decay1_func Coef_Decay Decay[Cursors[0]/Timex[1],Cursors[1]/Timex[1]] /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_Decay
	Coef_Decay[15]=V_chisq
End

Function Fit_dataset_2 ()
	wave Decay = root:Decay2 	//Declare all needed global variables
	wave fit_Decay = root:fit_Decay2
	wave Timex=root:Time2
	wave Coef_Decay=root:Coef_Decay2
	wave SD_Decay=root:SD_decay2
	wave eps_Coef_Decay=root:eps_Coef_Decay
	wave Cursors=root:Cursors
	wave/T Fix=root:Fix2
	SVAR fixstr

	Duplicate/O Decay SD_Decay
	If (numpnts(fit_Decay)!=numpnts(Decay))
		Duplicate/O Decay fit_Decay
	Endif	
	fixstr=Fix[0]+Fix[1]+Fix[2]+Fix[3]+Fix[4]+Fix[5]+Fix[6]+Fix[7]+Fix[8]+Fix[9]+Fix[10]+Fix[11]+Fix[12]+Fix[13]+Fix[14]+"1"
	SD_Decay= sqrt(Fit_Decay)
	FuncFit/L=2000 /H=fixstr Decay2_func Coef_Decay Decay[Cursors[0]/Timex[1],Cursors[1]/Timex[1]] /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_Decay
	SD_Decay= sqrt(fit_Decay)
	FuncFit/L=2000 /H=fixstr Decay2_func Coef_Decay Decay[Cursors[0]/Timex[1],Cursors[1]/Timex[1]] /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_Decay
	SD_Decay= sqrt(fit_Decay)
	FuncFit/L=2000 /H=fixstr Decay2_func Coef_Decay Decay[Cursors[0]/Timex[1],Cursors[1]/Timex[1]] /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_Decay
	Coef_Decay[15]=V_chisq
End

Function Fit_dataset_3 ()
	wave Decay = root:Decay3 	//Declare all needed global variables
	wave fit_Decay = root:fit_Decay3
	wave Timex=root:Time3
	wave Coef_Decay=root:Coef_Decay3
	wave SD_Decay=root:SD_decay3
	wave eps_Coef_Decay=root:eps_Coef_Decay
	wave Cursors=root:Cursors
	wave/T Fix=root:Fix3
	SVAR fixstr

	Duplicate/O Decay SD_Decay
	If (numpnts(fit_Decay)!=numpnts(Decay))
		Duplicate/O Decay fit_Decay
	Endif	
	fixstr=Fix[0]+Fix[1]+Fix[2]+Fix[3]+Fix[4]+Fix[5]+Fix[6]+Fix[7]+Fix[8]+Fix[9]+Fix[10]+Fix[11]+Fix[12]+Fix[13]+Fix[14]+"1"
	SD_Decay= sqrt(Fit_Decay)
	FuncFit/L=2000 /H=fixstr Decay3_func Coef_Decay Decay[Cursors[0]/Timex[1],Cursors[1]/Timex[1]] /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_Decay
	SD_Decay= sqrt(fit_Decay)
	FuncFit/L=2000 /H=fixstr Decay3_func Coef_Decay Decay[Cursors[0]/Timex[1],Cursors[1]/Timex[1]] /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_Decay
	SD_Decay= sqrt(fit_Decay)
	FuncFit/L=2000 /H=fixstr Decay3_func Coef_Decay Decay[Cursors[0]/Timex[1],Cursors[1]/Timex[1]] /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_Decay
	Coef_Decay[15]=V_chisq
End

Function Fit_dataset_4 ()
	wave Decay = root:Decay4 	//Declare all needed global variables
	wave fit_Decay = root:fit_Decay4
	wave Timex=root:Time4
	wave Coef_Decay=root:Coef_Decay4
	wave SD_Decay=root:SD_decay4
	wave eps_Coef_Decay=root:eps_Coef_Decay
	wave Cursors=root:Cursors
	wave/T Fix=root:Fix4
	SVAR fixstr

	Duplicate/O Decay SD_Decay
	If (numpnts(fit_Decay)!=numpnts(Decay))
		Duplicate/O Decay fit_Decay
	Endif	
	fixstr=Fix[0]+Fix[1]+Fix[2]+Fix[3]+Fix[4]+Fix[5]+Fix[6]+Fix[7]+Fix[8]+Fix[9]+Fix[10]+Fix[11]+Fix[12]+Fix[13]+Fix[14]+"1"
	SD_Decay= sqrt(Fit_Decay)
	FuncFit/L=2000 /H=fixstr Decay4_func Coef_Decay Decay[Cursors[0]/Timex[1],Cursors[1]/Timex[1]] /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_Decay
	SD_Decay= sqrt(fit_Decay)
	FuncFit/L=2000 /H=fixstr Decay4_func Coef_Decay Decay[Cursors[0]/Timex[1],Cursors[1]/Timex[1]] /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_Decay
	SD_Decay= sqrt(fit_Decay)
	FuncFit/L=2000 /H=fixstr Decay4_func Coef_Decay Decay[Cursors[0]/Timex[1],Cursors[1]/Timex[1]] /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_Decay
	Coef_Decay[15]=V_chisq
End

Function Fit_buffer_1 ()
	wave DecayB = root:DecayB1 	//Declare all needed global variables
	wave fit_DecayB = root:fit_DecayB1
	wave Timex=root:Time1
	wave Coef_DecayB=root:Coef_DecayB1
	wave SD_Decay=root:SD_decay1
	wave eps_Coef_DecayB=root:eps_Coef_DecayB
	wave/T T_Constraints =root:T_Constraints
	wave/T Fix=root:FixB
	SVAR fixstr

	Duplicate/O DecayB SD_Decay
	If (numpnts(fit_DecayB)!=numpnts(DecayB))
		Duplicate/O DecayB fit_DecayB
	Endif
	fixstr=Fix[0]+Fix[1]+Fix[2]+Fix[3]+Fix[4]+Fix[5]+Fix[6]+Fix[7]+Fix[8]+Fix[9]+Fix[10]+"1"
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB1_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB1_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB1_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB1_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB1_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	Coef_DecayB[11]=V_chisq
End

Function Fit_buffer_2 ()
	wave DecayB = root:DecayB2 	//Declare all needed global variables
	wave fit_DecayB = root:fit_DecayB2
	wave Timex=root:Time2
	wave Coef_DecayB=root:Coef_DecayB2
	wave SD_Decay=root:SD_decay2
	wave eps_Coef_DecayB=root:eps_Coef_DecayB
	wave/T T_Constraints =root:T_Constraints
	wave/T Fix=root:FixB
	SVAR fixstr

	Duplicate/O DecayB SD_Decay
	If (numpnts(fit_DecayB)!=numpnts(DecayB))
		Duplicate/O DecayB fit_DecayB
	Endif
	fixstr=Fix[0]+Fix[1]+Fix[2]+Fix[3]+Fix[4]+Fix[5]+Fix[6]+Fix[7]+Fix[8]+Fix[9]+Fix[10]+"1"
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB2_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB2_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB2_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB2_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB2_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	Coef_DecayB[11]=V_chisq
End

Function Fit_buffer_3 ()
	wave DecayB = root:DecayB3 	//Declare all needed global variables
	wave fit_DecayB = root:fit_DecayB3
	wave Timex=root:Time3
	wave Coef_DecayB=root:Coef_DecayB3
	wave SD_Decay=root:SD_decay3
	wave eps_Coef_DecayB=root:eps_Coef_DecayB
	wave/T T_Constraints =root:T_Constraints
	wave/T Fix=root:FixB
	SVAR fixstr

	Duplicate/O DecayB SD_Decay
	If (numpnts(fit_DecayB)!=numpnts(DecayB))
		Duplicate/O DecayB fit_DecayB
	Endif
	fixstr=Fix[0]+Fix[1]+Fix[2]+Fix[3]+Fix[4]+Fix[5]+Fix[6]+Fix[7]+Fix[8]+Fix[9]+Fix[10]+"1"
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB3_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB3_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB3_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB3_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB3_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	Coef_DecayB[11]=V_chisq
End

Function Fit_buffer_4 ()
	wave DecayB = root:DecayB4 	//Declare all needed global variables
	wave fit_DecayB = root:fit_DecayB4
	wave Timex=root:Time4
	wave Coef_DecayB=root:Coef_DecayB4
	wave SD_Decay=root:SD_decay4
	wave eps_Coef_DecayB=root:eps_Coef_DecayB
	wave/T T_Constraints =root:T_Constraints
	wave/T Fix=root:FixB
	SVAR fixstr

	Duplicate/O DecayB SD_Decay
	If (numpnts(fit_DecayB)!=numpnts(DecayB))
		Duplicate/O DecayB fit_DecayB
	Endif
	fixstr=Fix[0]+Fix[1]+Fix[2]+Fix[3]+Fix[4]+Fix[5]+Fix[6]+Fix[7]+Fix[8]+Fix[9]+Fix[10]+"1"
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB4_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB4_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB4_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB4_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	SD_Decay= sqrt(fit_DecayB)
	FuncFit/L=2000 /H=fixstr DecayB4_func Coef_DecayB DecayB /X=Timex /W=SD_Decay /I=1 /E=eps_Coef_DecayB /C=T_Constraints
	Coef_DecayB[11]=V_chisq
End

Function ChiSquareSurface(a,b,c,d, name)
   Variable a,b,c,d
   String name
   wave FitFuncNames, DataSets, CoefDataSetLinkage,CoefWave, CoefNames
 	Variable/G GF_chisq
 	wave FitFuncNames=root:Packages:NewGlobalFit:NewGF_FitFuncNames
 	wave DataSets=root:Packages:NewGlobalFit:NewGF_DataSetsList
 	wave CoefDataSetLinkage=root:Packages:NewGlobalFit:NewGF_LinkageMatrix
 	wave CoefWave=root:Packages:NewGlobalFit:NewGF_CoefWave
 	wave CoefNames=root:Packages:NewGlobalFit:NewGF_CoefficientNames
   WAVE/T GF_ConstraintWave=root:GF_ConstraintWave
  
   Duplicate/T/O GF_ConstraintWave ConstraintWave
    
   String searchString   // Substring to search for
   searchString="K"+num2str(d)+" "

   Variable numPoints = numpnts(ConstraintWave)
   Variable i

   for (i = numPoints - 1; i >= 0; i -= 1)
       if (strsearch(ConstraintWave[i], searchString, 0) >= 0)  // Check if substring exists in the cell
          DeletePoints i,1, ConstraintWave
       endif
   endfor

 	 	 	 	
   Variable id= (c-b)/a
   Variable ic
   
   Duplicate/FREE/O CoefWave, CoefWave1 
   String endwave= "ChiSqWave_" + name 
   Make/FREE/O/N=(DimSize(CoefWave,0) + 1, a+1), ChiSqWave
   for (ic=0;ic<a+1; ic+=1)
 	  CoefWave1[d][1]= 1
     CoefWave1[d][0]= b + (ic * id)
     DoNewGlobalFit(FitFuncNames, DataSets, CoefDataSetLinkage, CoefWave1, CoefNames, ConstraintWave, NewGFOptionWTISSTD+NewGFOptionMAKE_FIT_WAVES+NewGFOptionCOR_MATRIX+NewGFOptionGLOBALFITVARS, 1000, 1, maxIters=100)
     Print ic, GF_chisq
     ChiSqWave[0,DimSize(CoefWave,0)-1][ic]= CoefWave1[p][0]
     ChiSqWave[DimSize(CoefWave,0)][ic]= GF_chisq
     CoefWave1=CoefWave
   endfor
   Duplicate/O ChiSqWave, $(endwave)
end

Function p2ChiSquareSurface(a,b,c,d,e,f,g,h, name)
   Variable a,b,c,d,e,f,g,h 
   String name
   wave FitFuncNames, DataSets, CoefDataSetLinkage,CoefWave, CoefNames
 	Variable/G GF_chisq
 	wave FitFuncNames=root:Packages:NewGlobalFit:NewGF_FitFuncNames
 	wave DataSets=root:Packages:NewGlobalFit:NewGF_DataSetsList
 	wave CoefDataSetLinkage=root:Packages:NewGlobalFit:NewGF_LinkageMatrix
 	wave CoefWave=root:Packages:NewGlobalFit:NewGF_CoefWave
 	wave CoefNames=root:Packages:NewGlobalFit:NewGF_CoefficientNames
   WAVE/T GF_ConstraintWave=root:GF_ConstraintWave
  
   Duplicate/T/O GF_ConstraintWave ConstraintWave
    
   String searchString1, searchString2   // Substring to search for
   searchString1="K"+num2str(d)+" "
   searchString2="K"+num2str(h)+" "
  
   Variable numPoints = numpnts(ConstraintWave)
   Variable i

   for (i = numPoints - 1; i >= 0; i -= 1)
       if (strsearch(ConstraintWave[i], searchString1, 0) >= 0)  // Check if substring exists in the cell
          DeletePoints i,1, ConstraintWave
       endif
   endfor
   
	numPoints = numpnts(ConstraintWave)
   for (i = numPoints - 1; i >= 0; i -= 1)
       if (strsearch(ConstraintWave[i], searchString2, 0) >= 0)  // Check if substring exists in the cell
          DeletePoints i,1, ConstraintWave
       endif
   endfor
   
   Variable id= (c-b)/a
   Variable ic
   
   Variable ih= (g-f)/e
   Variable ig
   
   Variable qq = 0
    
   Duplicate/O CoefWave, CoefWave1
   
   String endwave= "p2ChiSqWave_" + name 
   String endwave1="p2ChiSqSurf_" + name 
   Make/O/N=(DimSize(CoefWave,0) + 1, a+1, e+1), p2ChiSqWave
   Make/O/N=((a+1)*(e+1),3), p2ChiSqSurf
  
   for (ic=0;ic<a+1; ic+=1)
 	  CoefWave1[d][1]= 1
     CoefWave1[d][0]= b + (ic * id)

     for (ig=0;ig<e+1; ig+=1)
  		CoefWave1[h][1]= 1
      	CoefWave1[h][0]= f + (ig * ih)
      DoNewGlobalFit(FitFuncNames, DataSets, CoefDataSetLinkage, CoefWave1, CoefNames, ConstraintWave, NewGFOptionWTISSTD+NewGFOptionMAKE_FIT_WAVES+NewGFOptionCOR_MATRIX+NewGFOptionGLOBALFITVARS, 1000, 1, maxIters=100)

    	p2ChiSqWave[0,DimSize(CoefWave,0)-1][ic][ig]= CoefWave1[p][0]
    	p2ChiSqWave[DimSize(CoefWave,0)][ic][ig]= GF_chisq
    	p2ChiSqSurf[qq][0]= CoefWave1[d][0]
    	p2ChiSqSurf[qq][1]= CoefWave1[h][0]
    	p2ChiSqSurf[qq][2]= GF_chisq 
      Print ic, ig, p2ChiSqSurf[qq][0], p2ChiSqSurf[qq][1], GF_chisq
    	qq+=1
    endfor 
  	endfor
   Duplicate/O p2ChiSqWave, $(endwave)
   Duplicate/O p2ChiSqSurf, $(endwave1)
end

Menu "Macros"		//Put macros in the Macro menu
	"Fit dataset 1", Fit_dataset_1 ()
	"Fit dataset 2", Fit_dataset_2 ()
	"Fit dataset 3", Fit_dataset_3 ()
	"Fit dataset 4", Fit_dataset_4 ()
	"Fit buffer 1", Fit_buffer_1 ()
	"Fit buffer 2", Fit_buffer_2 ()
	"Fit buffer 3", Fit_buffer_3 ()
	"Fit buffer 4", Fit_buffer_4 ()
	"FRET_efficiency 1", print eff_func(Coef_Decay1, 0)
	"FRET_efficiency 2", print eff_func(Coef_Decay2, 0)
	"FRET_efficiency 3", print eff_func(Coef_Decay3, 0)
	"FRET_efficiency 4", print eff_func(Coef_Decay4, 0)
End


Window Graph13() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(479.25,40.25,1020,564.5)
	AppendXYZContour 'p2ChiSqSurf_avg_r1,sig_r1'
	ModifyContour 'p2ChiSqSurf_avg_r1,sig_r1' moreLevels={23700,24000,25000,26000,27000}
	ModifyContour 'p2ChiSqSurf_avg_r1,sig_r1' moreLevels={28000,29000,30000,40000,50000}
	ModifyContour 'p2ChiSqSurf_avg_r1,sig_r1' moreLevels={60000,70000,80000}
	ModifyGraph width={Aspect,1}
	ModifyGraph hideTrace('p2ChiSqSurf_avg_r=boundary')=1
	ModifyGraph mirror=0
	Label left "standard deviation (\\f02σ\\B1\\M\\f00), Å"
	Label bottom "average distance (\\f02 ̅r\\B1\\M\\f00), Å"
EndMacro
