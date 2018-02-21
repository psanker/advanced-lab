@echo off

REM *********************************************************
REM 	CPMG.bat 
REM 	This file is intended as an example of using the CPMG executable with a batch file.
REM 	Details about the specific CPMG parameters are included in the README-CPMG.html file. 

REM 	SpinCore Technologies, Inc.
REM 	2009/10/12 12:56:00
REM *********************************************************

REM BOARD_NUMBER is the number of the board in your system to be used by CPMG. If you have multiple boards attached to your system, please make sure this value is correct.
SET BOARD_NUMBER=0

REM SPECTROMETER_FREQUENCY (MHz) must be between 0 and 100.
set SPECTROMETER_FREQUENCY=21.12

REM SPECTRAL_WIDTH (kHz) must be between 0.150 and 10000.  Equals the number of points acquired per millisecond.
set SPECTRAL_WIDTH=30

REM P2_TIME (microseconds): 90 degree pulse time must be atleast 0.065. 
set P2_TIME=4.4

REM RINGDOWN_TIME (microseconds) must be atleast 0.065.
set RINGDOWN_TIME=150

REM P2_PHASE (degrees): Output (Tx) phase of 90 degree pulse.
set P2_PHASE=285.0

REM P1_PHASE (degrees): Output (Tx) phase of 180 degree pulse.  Should be P2_PHASE + 90 degrees in most cases
set P1_PHASE=15.0

REM TAU (microseconds) must be at least 0.065.  Equals half the time between 180 degree echos and should be set approximately to the FID decay time of the 90 degree pulse.
set TAU=20000

REM NUMBER_OF_POINTS_PER_ECHO must be greater than or equal to 0.  If set to zero the program will perform a continuous scan, starting immediately after the ringdown time and stopping after the last echo loop.
set NUMBER_OF_POINTS_PER_ECHO=0

REM NUMBER_OF_ECHOES must be greater than or equal to 0.
set NUMBER_OF_ECHOES=4

REM NUMBER_OF_SCANS must be greater than or equal to 1.  Multiple scans are autonomously averaged.
set NUMBER_OF_SCANS=4

REM FILENAME is the name of the output file the data will be acquired data will be stored in. File extensions will be appended automatically.
set FILENAME=cpmg

REM BYPASS_FIR will disable the FIR filter if set to 1. Setting BYPASS_FIR to 0 will enable the FIR filter.
set BYPASS_FIR=1

REM ADC_FREQUENCY (MHz) is the analog to digital converter frequency of the RadioProcessor board selected.
SET ADC_FREQUENCY=75.0

REM REPETITION_DELAY (s) is the time between each consecutive scan. It must be greater than 0.065.
SET REPETITION_DELAY=1.0

REM AMPLITUDE is the (non-linear) output scaling factor for the transmitter output.  Must be between 0.0 and 1.0.
SET AMPLITUDE=0.2

REM BLANKING_BIT is the number of the TTL bit to use for the power amplifier blanking signal.
SET BLANKING_BIT=2

REM BLANKING_DELAY (ms) is the time before the excitation pulse necessary to blank the power amplifier. 3.0ms is a good number for the 10W PA.
SET BLANKING_DELAY=3.0

REM DEBUG Enables the debug output log.
SET DEBUG=0

CPMG %BOARD_NUMBER% %SPECTROMETER_FREQUENCY% %SPECTRAL_WIDTH% %P2_TIME% %RINGDOWN_TIME% %P2_PHASE% %P1_PHASE% %TAU% %NUMBER_OF_POINTS_PER_ECHO% %NUMBER_OF_ECHOES% %NUMBER_OF_SCANS% %FILENAME% %BYPASS_FIR% %ADC_FREQUENCY% %REPETITION_DELAY% %AMPLITUDE% %BLANKING_BIT% %BLANKING_DELAY% %DEBUG%
pause
