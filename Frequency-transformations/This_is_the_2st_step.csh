#!/bin/csh
# For automatically execute step1.macro for each SACs.
# 2010.04.16 By suntian

if ( -e ./put_sac_txt) then
 echo "./put_sac_txt is exsited."
 rm ./put_sac_txt/*
else
 mkdir ./put_sac_txt
endif

set runtime = 0

while ($runtime < 5)
echo $runtime
sac << EOF
macro ./divide_fft.macro day $runtime file ${runtime}day.SAC
q
EOF

rm temp.sac

./sacr.csh

mv -f *part.sac ./put_sac_txt

./pick_amplitude_from_band_spectrum.csh

mv pick_amp_for_each_windows_0.004-0.01Hz.txt ${runtime}day_pick_amp_for_each_windows_0.004-0.01Hz.txt
mv pick_amp_for_each_windows_0.05-0.5Hz.txt ${runtime}day_pick_amp_for_each_windows_0.05-0.5Hz.txt
mv pick_amp_for_each_windows_1.3-1.5Hz.txt ${runtime}day_pick_amp_for_each_windows_1.3-1.5Hz.txt

mv -f *part.sac.txt ./put_sac_txt

@ runtime ++

end

cat [0-9]day_pick_amp_for_each_windows_0.004-0.01Hz.txt > allday_pick_amp_0.004-0.01Hz.txt
cat [0-9]day_pick_amp_for_each_windows_0.05-0.5Hz.txt > allday_pick_amp_0.05-0.5Hz.txt
cat [0-9]day_pick_amp_for_each_windows_1.3-1.5Hz.txt > allday_pick_amp_1.3-1.5Hz.txt

rm *day_pick_amp_for_each_windows_0.004-0.01Hz.txt
rm *day_pick_amp_for_each_windows_0.05-0.5Hz.txt
rm *day_pick_amp_for_each_windows_1.3-1.5Hz.txt

mkdir data_we_obtain
mv -f allday_pick_amp_0.004-0.01Hz.txt ./data_we_obtain
mv -f allday_pick_amp_0.05-0.5Hz.txt ./data_we_obtain
mv -f allday_pick_amp_1.3-1.5Hz.txt ./data_we_obtain

