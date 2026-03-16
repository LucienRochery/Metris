set -e

#wget -P /progs/OCC/ https://acdl.mit.edu/ESP/OCC781-linux-aarch64.tgz
#cd /progs/OCC/
## Get rid of top level directory so we don't need to handle version names
#tar --strip-components=1 -xvf OCC781-linux-aarch64.tgz 
#OCC_DIR=/progs/OCC



#wget -P /progs/ESP/ https://acdl.mit.edu/ESP/ESP.tgz
#cd /progs/ESP/
#tar --strip-components=1 -xvf ESP.tgz

echo "-- START configuring EGADS"
cd /progs/ESP/config/
./makeEnv
echo "-- DONE configuring EGADS"

cd ../src/
. ../ESPenv.sh
make -j12 
