#ifndef __TRIALS_H__
#define __TRIALS_H__

#include<iostream>
#include"data_type.h"
#include"SysUtil.h"

using namespace std;



// ------- updated 7/13/2016 -------
txProfileType	GetTx_G1_SOF_LDV12	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_G1_SOF_DCV12_24	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_G3_SOF_DCV12_24	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_G4_SOF_LDV12	(const modelParamType & argModelParam,const patientType & argPat);


// ------- updated 12/12/2015 -------
txProfileType	GetTx_TN_TOL_G1b_SOF_LDV_ASV3	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TN_TOL_G1b_SOF_DCV_SMV3	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TN_TOL_G1b_SOF_DCV_ASV3	(const modelParamType & argModelParam,const patientType & argPat);


// ------- updated 8/28/2015 --------

txProfileType	GetTx_TN_TOL_G1a_DCV_SOF12_24	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TN_TOL_G1a_PrOD_RBV12_24	(const modelParamType & argModelParam,const patientType & argPat);

txProfileType	GetTx_TN_TOL_G1b_DCV_SOF12_24	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TN_TOL_G1b_PrOD12	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TN_TOL_G1b_SOF_SMV12	(const modelParamType & argModelParam,const patientType & argPat);

txProfileType	GetTx_TN_TOL_G2_SOF_RBV12_16	(const modelParamType & argModelParam,const patientType & argPat);


txProfileType	GetTx_TN_TOL_G3_DCV_SOF12_24	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TN_TOL_G3_SOF_PEG_RBV12	(const modelParamType & argModelParam,const patientType & argPat);


txProfileType	GetTx_TN_TOL_G4_LDV_SOF12	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TN_TOL_G4_PrOD_RBV12	(const modelParamType & argModelParam,const patientType & argPat);
//txProfileType	GetTx_TN_TOL_G4_SOF_PEG_RBV12	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TN_TOL_G4_SOF_RBV24	(const modelParamType & argModelParam,const patientType & argPat);


txProfileType	GetTx_TE_G1a_DCV_SOF12_24	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TE_G1a_PrOD_RBV12_24	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TE_G1a_SMV_SOF12	(const modelParamType & argModelParam,const patientType & argPat);


txProfileType	GetTx_TE_G1b_DCV_SOF12_24	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TE_G1b_PrOD12	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TE_G1b_SMV_SOF12	(const modelParamType & argModelParam,const patientType & argPat);

txProfileType	GetTx_TE_G2_LDV_RBV12_16	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TE_G2_SOF_PEG_RBV12	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TE_G2_SMV_SOF12	(const modelParamType & argModelParam,const patientType & argPat);

txProfileType	GetTx_TE_G3_DCV_SOF12_24	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TE_G3_SOF_PEG_RBV12	(const modelParamType & argModelParam,const patientType & argPat);


txProfileType	GetTx_TE_G4_LDV_SOF12	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TE_G4_PrOD_RBV12	(const modelParamType & argModelParam,const patientType & argPat);
txProfileType	GetTx_TE_G4_SOF_RBV24	(const modelParamType & argModelParam,const patientType & argPat);
//txProfileType	GetTx_TE_G4_SOF_PEG_RBV12	(const modelParamType & argModelParam,const patientType & argPat);


// ------- new treatment ------------
// TN TOL
txProfileType GetTx_TN_G1_LDP_SOF8(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TN_G1_LDP_SOF12(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TN_TOL_G2_SOF_RBV12(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TN_TOL_G3_SOF_RBV24(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TN_TOL_G4_SOF_PEG_RBV12(const modelParamType & argModelParam,const patientType & argPat);

// TN INTOL
txProfileType GetTx_TN_NOT_G1_SOF_SMV12(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TN_NOT_G2_SOF_RBV12(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TN_NOT_G3_SOF_RBV24(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TN_NOT_G4_SOF_RBV24(const modelParamType & argModelParam,const patientType & argPat);

// TE
txProfileType GetTx_TE_G1_LDP_SOF12_24(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TE_G2_SOF_RBV12(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TE_G3_SOF_RBV24(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TE_G4_SOF_PEG_RBV12(const modelParamType & argModelParam,const patientType & argPat);

// ------- old SOC ------------
// TN TOL
txProfileType GetTx_TN_TOL_G1_BOC(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TN_TOL_G1_TEL(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TN_TOL_G2_PEG_RBV24( const modelParamType & argModelParam,const patientType & argPat );
txProfileType GetTx_TN_TOL_G3_PEG_RBV24( const modelParamType & argModelParam,const patientType & argPat );
txProfileType GetTx_TN_TOL_G4_PEG_RBV48( const modelParamType & argModelParam,const patientType & argPat );

// TN INTOL
inline txProfileType GetNoTx(){txProfileType theTx; theTx._txDuration=0; return theTx;}	// return an "empty" object

// TE
txProfileType GetTx_TE_G1_BOC(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TE_G1_TEL(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TE_G2_PEG_RBV24( const modelParamType & argModelParam,const patientType & argPat );
txProfileType GetTx_TE_G3_PEG_RBV24( const modelParamType & argModelParam,const patientType & argPat );
txProfileType GetTx_TE_G4_PEG_RBV48( const modelParamType & argModelParam,const patientType & argPat );





// ----- others -----
txProfileType GetTx_TN_TOL_G1_SOF_PEG_RBV12(const modelParamType & argModelParam,const patientType & argPat);
txProfileType GetTx_TE_G1_SOF_SMV12(const modelParamType & argModelParam,const patientType & argPat);


txProfileType GetTx_FBOP(const modelParamType & argModelParam,const patientType & argPat);
					





//txProfileType GetTx_TN_G1_F0F3_LDP_SOF8(const modelParamType & argModelParam,const patientType & argPat);
//txProfileType GetTx_TN_G1_F4_LDP_SOF12(const modelParamType & argModelParam,const patientType & argPat);
//txProfileType GetTx_TE_G1_F0F3_LDP_SOF12(const modelParamType & argModelParam,const patientType & argPat);
//txProfileType GetTx_TE_G1_F4_LDP_SOF24(const modelParamType & argModelParam,const patientType & argPat);
//
//txProfileType GetTx_TN_G1_LDP_SOF(const modelParamType & argModelParam,const patientType & argPat);
//txProfileType GetTx_TE_G1_LDP_SOF(const modelParamType & argModelParam,const patientType & argPat);


txProfileType GetTx(const modelParamType & argModelParam, const patientType & argPat);


// ----------------- BELOW THIS POINT: 2016/12/2 added for retreatment -------------------------------
txProfileType GetTx_Retreatment(const modelParamType & argModelParam, const patientType & argPat);
txProfileType GetReTx_TreatmentNameX(const modelParamType & argModelParam, const patientType & argPat);

// ----------------- BELOW THIS POINT: 2016/11/22 added for acute HCV phase --------------------------

txProfileType GetTx_Acute(const modelParamType & argModelParam, const patientType & argPat);
txProfileType GetTx_Acute_G1_TreatmentX(const modelParamType & argModelParam, const patientType & argPat);
txProfileType GetTx_Acute_G1_TreatmentX_ChangableDuration(const modelParamType & argModelParam, const patientType & argPat);

txProfileType GetTx_G1_SOF_LDV8_NoDiscontNoAE(const modelParamType & argModelParam, const patientType & argPat);
txProfileType GetTx_G1_SOF_LDV12_NoDiscontNoAE_ChangableDuration(const modelParamType & argModelParam, const patientType & argPat);

#endif