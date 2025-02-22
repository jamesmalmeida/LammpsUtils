#Generate a charge file for ReaxFF Simulations. All charges are set to zero.

input="initial.xyz"

echo "xyz2charge" > charge.top
echo " " >> charge.top
natoms=$(head -n1 ${input})
echo "${natoms} atoms" >> charge.top
echo "14 atom types" >> charge.top
echo " " >> charge.top
echo "0.0 40.0 xlo xhi" >> charge.top
echo "0.0 40.0 ylo yhi" >> charge.top
echo "0.0 65.0 zlo zhi" >> charge.top
echo " " >> charge.top
echo "Masses" >> charge.top
echo "" >> charge.top
echo "      1    12.01100 #C " >> charge.top
echo "      2     1.00800 #H " >> charge.top
echo "      3    15.99940 #O " >> charge.top
echo "      4    14.00700 #N " >> charge.top
echo "      5    32.06000 #S " >> charge.top
echo "      6    24.30500 #Mg" >> charge.top
echo "      7    30.97400 #P " >> charge.top
echo "      8    22.99000 #Na" >> charge.top
echo "      9    47.86700 #Ti" >> charge.top
echo "     10    35.45000 #Cl" >> charge.top
echo "     11    18.99800 #F " >> charge.top
echo "     12    39.09800 #K " >> charge.top
echo "     13     6.94000 #Li" >> charge.top
echo "     14    99.99999 #X " >> charge.top
echo " " >> charge.top
echo "Atoms" >> charge.top
echo " " >> charge.top
cat ${input} | grep "C \|H \|O \|N \|S \|Mg \|P \|Na \|Ti \|Cl \|F \|K \|Li " | awk '{print NR,"1",$1,"0.00",$2,$3,$4}' | sed 's/ C / 1 /g' | sed 's/ H / 2 /g' | sed 's/ O / 3 /g' | sed 's/ N / 4 /g' | sed 's/ S / 5 /g' | sed 's/ Mg / 6 /g' | sed 's/ P / 7 /g' | sed 's/ Na / 8 /g' | sed 's/ Ti / 9 /g' | sed 's/ Cl / 10 /g' | sed 's/ F / 11 /g' | sed 's/ K / 12 /g' | sed 's/ Li / 13 /g' >> charge.top
