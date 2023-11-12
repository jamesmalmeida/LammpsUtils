#Generate a charge file for ReaxFF Simulations. All charges are set to zero.

input="initial.xyz"

echo "xyz2charge" > charge.top
echo " " >> charge.top
natoms=$(head -n1 ${input})
echo "${natoms} atoms" >> charge.top
echo "3 atom types" >> charge.top
echo " " >> charge.top
echo "0.0 25.0 xlo xhi" >> charge.top
echo "0.0 25.0 ylo yhi" >> charge.top
echo "0.0 25.0 zlo zhi" >> charge.top
echo " " >> charge.top
echo "Masses" >> charge.top
echo "      1    12.01100" >> charge.top
echo "      2     1.00800" >> charge.top
echo "      3    15.99940" >> charge.top
echo " " >> charge.top
echo "Atoms" >> charge.top
echo " " >> charge.top
cat ${input} | grep "C\|H\|O" | awk '{print NR,"1",$1,"0.00",$2,$3,$4}' | sed 's/ C / 1 /g' | sed 's/ H / 2 /g' | sed 's/ O / 3 /g' >> charge.top

