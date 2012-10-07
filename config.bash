export PYTHONPATH=/home/huskeypm/localTemp/srcs/homogenization/:$PYTHONPATH
export PYTHONPATH=/home/huskeypm/localTemp/srcs/smolfin/:$PYTHONPATH
export PYTHONPATH=/home/huskeypm/bin/grids/:$PYTHONPATH

#Is this for my virtual machine?
export BOOST_DIR=/home/huskeypm/Work/FEniCS/
source $BOOST_DIR/share/dolfin/dolfin.conf
export PREFIX=~/localTemp/srcs/gamer-src/
export FETK_INCLUDE=$PREFIX/include
export FETK_LIBRARY=$PREFIX/lib
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PREFIX/lib/python2.6/site-packages:$PYTHONPATH
export PYTHONPATH=$PYTHONPATH:/home/huskeypm/localTemp/srcs/smol/
export PYTHONPATH=$PYTHONPATH:/home/huskeypm/bin/dolfin/
export PYTHONPATH=$PYTHONPATH:/home/huskeypm/bin/grids/
export BOOST_DIR=/home/huskeypm/Work/FEniCS/

# added for rocce
export PYTHONPATH=$PYTHONPATH:/home/huskeypm/sources/dolfin_smol/
#source ~/localTemp/srcs/smol/config.bash


