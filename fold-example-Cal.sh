###########################example
#./timing-mcc.sh -cpus 30 -cal1 J0203-0150_20190517Cal.fits -calp 1.00663296 -plot 1
#./timing-mcc.sh -nbins 128 -foldl 20 -cpus 30 -par C8_working.par -psr J0203-0150_20190517Timing.fits -plot 1 -Calfile J0203-0150_20190517Cal.txt
#./timing-mcc.sh -cpus 30 -par C8_working.par -psrzap J0203-0150_20190517Timing.zap.pazi -plot 1 -Calfile J0203-0150_20190517Cal.txt
#cp J0203-0150_20190614Timing.zap.pTF pTF1.std 
#pat -f "tempo2" -s pTF1.std J0203-0150_20190517Timing.zap.pF >  J0203-0150_20190517Timing-C8_working.tim
###########################example

#dspsr -c 0.201326592 -L 10.0663296 -t 32 -A -O polCal/Cal /Cal1/*.fits
#dspsr -c 0.100663296 -L 10.0663296 -t 32 -A -O polCal/Cal /Cal1/*.fits
