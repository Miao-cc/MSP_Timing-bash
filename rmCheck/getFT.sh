filename=$1
pam -T -e calibP.T ${filename}
#pam -R -73 -e T.derot ${filename}.T
#pam -R -73 -F -e derot.TF ${filename}.T
#pam -pF -e pTF ${filename}.derot.TF
#pam -F -e F.derot ${filename}.T.derot
#pam -pF -e Fp.derot ${filename}.T.derot
#pdv -t -Z ${filename}.FT.derot > ${filename}.FT.derot.txt
#pdv -Z -t -A -K ${filename}.derot.TF > ${filename}.derot.TF.txt
