#!/bin/bash

POSTGRESQL_LIB_DIR=/usr/lib/postgresql/9.1/lib
OB_INSTALL_DIR=$POSTGRESQL_LIB_DIR/openbabel

# delete old libraries
sudo rm $POSTGRESQL_LIB_DIR/libinchi*
sudo rm $POSTGRESQL_LIB_DIR/libopenbabel.*
sudo rm $POSTGRESQL_LIB_DIR/libbarsoi.so
sudo rm $POSTGRESQL_LIB_DIR/inchiformat.so
sudo rm $POSTGRESQL_LIB_DIR/openbabel/ -rf
rm postgres* -rf
apt-get source postgresql-9.1
#had to be changed accordingly
mv postgresql-9.1-9.1.6 postgresql
cd postgresql
./configure


cd contrib
git clone https://github.com/bgruening/pgchem_tigress.git
mv pgchem_tigress pgchem
cd pgchem/src

# compile openbabel
mkdir openbabel-2.3.2/build
cd openbabel-2.3.2/build
cmake .. -DBUILD_SHARED=ON -DBUILD_GUI=OFF -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=$OB_INSTALL_DIR
make
sudo make install

# compile barsoi
cd ../../
cd barsoi
make -f Makefile.linux

# compile pgchem
cd ..
mv openbabel-2.3.2/include/openbabel/locale.h openbabel-2.3.2/include/openbabel/_locale.h
cd openbabel-2.3.2/build/lib/
ln -s ./inchiformat.so ./libinchiformat.so
cd ../../../
make -f Makefile.linux.x64

# copy all libraries
sudo cp libpgchem.so $POSTGRESQL_LIB_DIR/
sudo cp barsoi/libbarsoi.so $POSTGRESQL_LIB_DIR/
sudo cp openbabel-2.3.2/build/lib/libinchi.so.0.4.1 $POSTGRESQL_LIB_DIR/libinchi.so.0.4.1
sudo cp openbabel-2.3.2/build/lib/libopenbabel.so.4.0.2 $POSTGRESQL_LIB_DIR/libopenbabel.so.4.0.2
sudo cp openbabel-2.3.2/build/lib/inchiformat.so $POSTGRESQL_LIB_DIR/inchiformat.so
sudo cp setup/tigress/obdata/dictionary* $POSTGRESQL_LIB_DIR/openbabel/share/openbabel/2.3.2/
cd $POSTGRESQL_LIB_DIR
sudo ln -s libinchi.so.0.4.1 libinchi.so.0
sudo ln -s libinchi.so.0 libinchi.so 
sudo ln -s libopenbabel.so.4.0.2 libopenbabel.so.4
sudo ln -s libopenbabel.so.4 libopenbabel.so
sudo ln -s inchiformat.so libinchiformat.so 

# if you havn't configured ldconf to use your postgres libdir run the following commented lines
#echo '$POSTGRESQL_LIB_DIR/' | sudo tee -a /etc/ld.so.conf.d/libc.conf
#sudo ldconfig

# or preferably, set the shared libary path
# on server start (example below is CentOS)
# $ echo "LD_LIBRARY_PATH=/usr/lib64/pgsql" >> /etc/sysconfig/pgsql
# the edit postgresql.service and add the following under '[service]'
# EnvironmentFile=-/etc/sysconfig/pgsql

# remove unused files, only the lib dir is needed
sudo rm $OB_INSTALL_DIR/bin/ -rf
sudo rm $OB_INSTALL_DIR/include/ -rf

sudo /etc/init.d/postgresql restart


