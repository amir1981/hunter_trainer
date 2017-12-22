# file: asc.sh
#
# a simple driver script for the sound comparison application:
#

# demonstrate training
#
echo "... demonstrating training ..."

./asc \
train \
./models/0000_model.dat \
"./demo/f1.raw" \
"./demo/f2.raw" \
"./demo/f3.raw" \
"./demo/f4.raw"

# demonstrate evaluation
#
echo "... demonstrating evaluation ..."
./asc \
evaluate \
./models/0000_model.dat \
"./demo/f1.raw" \

# 
# end of file
