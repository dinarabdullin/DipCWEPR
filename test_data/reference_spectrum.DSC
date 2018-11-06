#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	MAN
BSEQ	BIG
IKKF	CPLX
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
IIFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	2500
XMIN	3114.000000
XWID	499.800000
*
*	Documentational Text:
*
TITL	'2018_07_26_NFL17_d8-Toluene_200uM_HCl_30dB_MA2G_100K_bckg'
IRNAM	'Intensity'
IINAM	''
XNAM	'Field'
IRUNI	''
IIUNI	''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    07/26/18
TIME    16:24:35
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    0.3364
A1SW    0.05
MWFQ    9.43978e+09
MWPW    0.0001915
AVGS    4
SPTP    0.048
RCAG    40
RCHM    1
B0MA    0.0002
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    2500
RCTC    0.01024
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*
.DVC     acqStart, 1.0
.DVC     ezcw, 1.0
MetalIonsDef       False
OrgRadicalsDef     False
.DVC     fieldCtrl, 1.0
AllegroMode        True
CenterField        3364.00 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
MeasuringHall      False
SetToSampleG       False
StaticFieldMon     3114.000 G
SweepDirection     Up
SweepWidth         500.0 G
WidthTM            200.0 G
.DVC     fieldSweep, 1.0
.DVC     freqCounter, 1.0
FrequencyMon       9.439780 GHz
QMonitBridge       On
.DVC     mwBridge, 1.0
AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
BridgeCalib        100.0
Power              0.1915 mW
PowerAtten         30.0 dB
.DVC     recorder, 1.0
BaselineCorr       Off
NbScansAcc         4
NbScansDone        4
NbScansToDo        4
ReplaceMode        Off
SmoothMode         Manual
SmoothPoints       0
.DVC     scanEnd, 1.0
.DVC     signalChannel, 1.0
AFCTrap            True
AllegroMode        True
AllowShortCt       False
Calibrated         True
ConvFact           9.271e-09
ConvTime           48.00 ms
DModDetectSCT      First
DualDetect         0/90 Degrees
EliDelay           1.0 us
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               40 dB
Harmonic           1
HighPass           True
Integrator         False
IsCalibExp         False
ModAmp             2.000 G
ModFreq            100.00 kHz
ModPhase           0.0
ModResol           10
Offset             0.0 %
PolyCof            {2;3,9;0} 0,0.99652,0,0,0.00737177,0,0,-0.00559614,0,0,-2.88221e-05,0,0,1.00404e-05,0,0,3.43695e-08,0,0,-5.0404e-09,0,0,-1.4783e-11,0,0,-1.29132e-12,0
QuadMode           True
ResCenter          62.5 mm
ResLength          40 mm
Resolution         2500
Resonator          1
SctNorm            True
SctRevision        Allegro
SpuExtension       False
SweepTime          120.00 s
TimeConst          10.24 ms
TimeExp            False
TuneCaps           1
*
************************************************************
*
#MHL	1.0 * MANIPULATION HISTORY LAYER by BRUKER
*
*
*       Input Dataset:
*
SOURCE_PRIM
'2018_07_26_NFL17_d8-Toluene_200uM_HCl_30dB_MA2G_100K'
END_SOURCE_PRIM
*
*       Input Dataset:
*
SOURCE_SCND
  *
  SOURCE
  '2018_07_26_NFL17_d8-Toluene_200uM_HCl_30dB_MA2G_100K'
  END_SOURCE
  SELECT  'qualiRegions'
  'qualiRegions'
  END_SELECT
  PROCESS 'prAverage'
  REPORT
                  a = -0.120765
  reduced chi-square = 3.173e-05
  END_REPORT
END_SOURCE_SCND
*
*       Selection:
*
SELECT  'qualiRegions'
'qualiRegions'
END_SELECT
*
*       Process:
*
PROCESS 'prDiff'
*
*       Parameter:
*
PAR_VAL Gain(Sec.) =  1.000e+00
PAR_VAL x-Shift(Sec.) =  0.000e+00
PAR_VAL x-Stretch(Sec.) =  1.000e+00
MDATE   07/27/18 15:22:57
*
************************************************************
