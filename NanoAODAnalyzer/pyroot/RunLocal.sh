#!/bin/sh

SAMPLES=(
# "WprimeWZM600"
"WprimeWZM800"
# "WprimeWZM1000"
# "WprimeWZM1200"
# "WprimeWZM1400"
# "WprimeWZM1600"
# "WprimeWZM2000"
# "WprimeWZM2500"
# "WprimeWZM3000"
# "WprimeWZM3500"
# "WprimeWZM4000"
# "WprimeWZM4500"
# "WprimeWZM5000"
# "WprimeWZM5500"
# "WprimeWZM6000"
# "WprimeWZM6500"
# "WprimeWZM7000"
# "WprimeWZM7500"
# "WprimeWZM8000"
# "QCDPt15To7000"
)

for s in "${SAMPLES[@]}"
do
  FULLPATH="../../SampleListNanoAOD/${s}.txt"
  python ReadNanoAOD.py --input ${FULLPATH} --outputDir "root://eosuser.cern.ch//eos/user/s/ssyedoma/AnaJetTagging/Ntuples/"
done
