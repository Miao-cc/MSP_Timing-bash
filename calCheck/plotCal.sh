cat calFilename.list | while read filename
do
  echo "Plot" ${filename} 
  outName=${filename:3:-4}
  pacv -D gainPlot/${outName}-gain.eps/CPS  ${filename}
  pav -G -g freqPlot/${outName}-freq.eps/CPS  ${filename}
  convert gainPlot/${outName}-gain.eps -background white -flatten -rotate 90 gainPlot/${outName}-gain.eps.png
  convert freqPlot/${outName}-freq.eps -background white -flatten -rotate 90 freqPlot/${outName}-freq.eps.png
done
