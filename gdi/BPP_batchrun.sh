#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=10:00:00
#SBATCH --ntasks=40

#A10 algorithm0
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_0_1/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_0_1/kentucki.bpp.A10.0.1.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_0_2/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_0_2/kentucki.bpp.A10.0.2.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_0_3/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_0_3/kentucki.bpp.A10.0.3.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_0_4/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_0_4/kentucki.bpp.A10.0.4.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_0_5/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_0_5/kentucki.bpp.A10.0.5.ctl

#A10 algorithm1
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_1_1/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10__1_1/kentucki.bpp.A10.1.1.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_1_2/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_1_2/kentucki.bpp.A10.1.2.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_1_3/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A10_0_3/kentucki.bpp.A10.1.3.ctl

#A11 algorithm0
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_0_1/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_0_1/kentucki.bpp.A11.0.1.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_0_2/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_0_2/kentucki.bpp.A11.0.2.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_0_3/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_0_3/kentucki.bpp.A11.0.3.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_0_4/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_0_4/kentucki.bpp.A11.0.4.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_0_5/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_0_5/kentucki.bpp.A11.0.5.ctl

##A11 algorithm1
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_1_1/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11__1_1/kentucki.bpp.A11.1.1.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_1_2/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_1_2/kentucki.bpp.A11.1.2.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_1_3/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_0_3/kentucki.bpp.A11.1.3.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_1_4/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11__1_4/kentucki.bpp.A11.1.4.ctl
cd /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_1_5/
/users/PHS0328/efwatts/bpp-4.4.1-linux-x86_64/bin/bpp --cfile /fs/scratch/PHS0328/kentucki_EFW/bpp/Final_round/kentucki_A11_1_5/kentucki.bpp.A11.1.5.ctl
