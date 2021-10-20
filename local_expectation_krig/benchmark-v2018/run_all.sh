mkdir DATA RESULTS
cd input
for name in inho_temp_*.dat
do
  pkl=${name/.dat/.pkl}
  ~/anaconda3/bin/python3 ../makepkl.py $name
  cp $pkl ../DATA/df_temp.pkl
  ~/anaconda3/bin/python3 ../../local_expectation_krig/calc_homogenization_pelt.py -years=1900,1999 -fourier=1
  ~/anaconda3/bin/python3 ../makedat.py
  mv df_temp_homog.dat ../RESULTS/$name
done

for name in inho_temp_*.dat
do
  armsd=${name/.dat/.rmsd}
  inv=${name/.dat/.inv}
  tru=${name/inho/orig}
  ~/anaconda3/bin/python3 ../rmsd.py $inv $tru ../RESULTS/$name > ../RESULTS/$armsd
done
