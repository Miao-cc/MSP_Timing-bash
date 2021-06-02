filename=$1
outName=${filename:3:-7}
echo "Plot" ${filename} ${outName}

pav -pdTG -g freqPlot/${outName}.eps/CPS  ${filename}
pav -pdFY -g timePlot/${outName}.eps/CPS  ${filename}
convert freqPlot/${outName}.eps -background white -flatten -rotate 90 freqPlot/${outName}.eps.png
convert timePlot/${outName}.eps -background white -flatten -rotate 90 timePlot/${outName}.eps.png
