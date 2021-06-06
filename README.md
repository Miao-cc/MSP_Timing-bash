# MSP_Timing-bash
Timing bash for MSP Timing

## update RA and DEC (for FAST)
python update_pos.py filenam ra dec

## combine fits and update the combined file header
combine_fits 
update

## 


## fold bash script
 - fold-example-2021.sh
 - fold-example-Cal.sh
fold pulsar and polarization calibrataion 


## check the polarization calibrataion
 - check_Noise_Time.py
python check_Noise_Time.py cal.fits


## check the polarization calibrataion fold result

```bash
cd calCheck
mkdir freqPlot gainPlot
ls ../*.ucf > calFilename.list
sh plotCal.sh
```
or 

```bash
cd calCheck
plotCal-single.sh ../cal.ucf
```

## check the pulsar fold result

```bash
cd foldCheck
mkdir freqPlot timePlot
ls ../*.calibP > foldFilename.list
sh plotFold.sh
```

or 

```bash
cd foldCheck
sh plotFold-single.sh pulsar.calibP
```

## Rotation Measurement

```bash
# sum in the time
ls *.calibP | xargs -i getFT.sh {}
# get the I, Q, U, I_err, Q_err, U_err
ls *.T | xargs -i polarDegree.py {}
# get the fit result
ls *.calibP.T.txt | xargs -i dofit.sh {}
```

### for timing-mcc.sh
1. Example:
- Cal File 
 > Cal File: J0631+4142_20190724Cal1.fits
 > cal TXT : J0631+4142_20190724Cal1.txt

- PSR File
 > PSR File: J0631+4142_20190724Timing.fits
 > par File: J0631+4142.par

> psrzap: J0631+4142_20190724Timing.zap.pazi
> calzap: J0631+4142_20190724Cal1.zap.pazi


 #### Fold PSR file

1. fold with par file and fits file
- with cal file and par file
 > ./timing-mcc.sh -psr J0631+4142_20190724Timing.fits -par J0631+4142.par -Calfile J0631+4142_20190724Cal1.txt -plot 1

  or

- with cal file and period, DM
 > ./timing-mcc.sh -psr J0631+4142_20190724Timing.fits -psrp 0.1 -psrDM 10 -Calfile J0631+4142_20190724Cal1.txt -plot 1

2. fold with pazi file
- with pazi file and cal file
 > ./timing-mcc.sh -psrzap J0631+4142_20190724Timing.fits -Calfile J0631+4142_20190724Cal1.txt -plot 1


#### Fold CAL file
1. fold with 1 cal fits file
 - with 1 cal fits
> ./timing-mcc.sh -cal1 J0631+4142_20190724Cal1.fits -calp 1.00663296 -plot 1

2. fold wit zap file
 - with pazi file
> ./timing-mcc.sh -calzap J0631+4142_20190724Cal1.zap.pazi -plot 1
