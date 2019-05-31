!--------------------------------------------------------------------------------
!
!  Copyright (C) 2017  L. J. Allen, H. G. Brown, A. J. DAlfonso, S.D. Findlay, B. D. Forbes
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!--------------------------------------------------------------------------------

module m_electron

    use m_precision
    implicit none


    interface elsa_ext
        module procedure single_elsa_ext
        module procedure double_elsa_ext
    end interface

    interface Peng_ionic_FF
        module procedure Peng_ionic_FF_fractional_ionicity
        module procedure Peng_ionic_FF_integer_ionicity
    end interface


!X-RAY FORM FACTORS
character(2), parameter :: element(92) = &
(/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ', &
  'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd', &
  'Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm', &
  'Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U '  /)
real(fp_kind), parameter :: xrayFF(14,92) = reshape( &
[1.00_fp_kind,0.73235400_fp_kind,0.7538960_fp_kind,0.2838190_fp_kind,0.1900030_fp_kind,0.0391390_fp_kind,&
               11.553918_fp_kind,4.59583100_fp_kind,1.54629900_fp_kind,26.4639640_fp_kind,0.37752300_fp_kind,&
              0.00048700_fp_kind,2.0_fp_kind,8.769017454422158E-3_fp_kind,&
  2.00_fp_kind,0.73235400_fp_kind,0.7538960_fp_kind,0.2838190_fp_kind,0.1900030_fp_kind,0.0391390_fp_kind,&
              11.553918_fp_kind,4.59583100_fp_kind,1.54629900_fp_kind,26.4639640_fp_kind,0.37752300_fp_kind,&
              0.00048700_fp_kind,2.0_fp_kind,8.769017454422158E-3_fp_kind,&
  3.00_fp_kind,0.97463700_fp_kind,0.1584720_fp_kind,0.8118550_fp_kind,0.2624160_fp_kind,0.7901080_fp_kind,&
              4.3349460_fp_kind,0.34245100_fp_kind,97.1029690_fp_kind,201.363824_fp_kind,1.40923400_fp_kind,&
              0.00254200_fp_kind,3.0_fp_kind,3.053829815542358E-2_fp_kind,&
  4.00_fp_kind,1.53371200_fp_kind,0.6382830_fp_kind,0.6010520_fp_kind,0.1061390_fp_kind,1.1184140_fp_kind,&
              42.662078_fp_kind,0.59542000_fp_kind,99.1065010_fp_kind,0.15134000_fp_kind,1.84309300_fp_kind,&
              0.00251100_fp_kind,4.0_fp_kind,2.673011412693727E-2_fp_kind,&
  5.00_fp_kind,2.08518500_fp_kind,1.0645800_fp_kind,1.0627880_fp_kind,0.1405150_fp_kind,0.6417840_fp_kind,&
              23.494069_fp_kind,1.13789400_fp_kind,61.2389750_fp_kind,0.11488600_fp_kind,0.39903600_fp_kind,&
               0.00382300_fp_kind,5.0_fp_kind,4.375691177284862E-2_fp_kind,&
  6.00_fp_kind,2.65750600_fp_kind,1.0780790_fp_kind,1.4909090_fp_kind,-4.241070_fp_kind,0.7137910_fp_kind,&
              14.780758_fp_kind,0.77677500_fp_kind,42.0868430_fp_kind,-0.0002940_fp_kind,0.23953500_fp_kind,&
              4.29798300_fp_kind, 6.0_fp_kind, 7.163539380752063E-2_fp_kind,&
  7.0_fp_kind,1.18937800e+01_fp_kind,3.27747900e+00_fp_kind,1.85809200e+00_fp_kind,8.58927000e-01_fp_kind,9.12985000e-01_fp_kind,&
               1.58000000e-04_fp_kind,1.02327230e+01_fp_kind,3.03446900e+01_fp_kind,6.56065000e-01_fp_kind,2.17287000e-01_fp_kind,&
               -1.18049020e+01_fp_kind,7.00000000e+00_fp_kind,1.12381925e-01_fp_kind,&
   8.0_fp_kind,2.96042700e+00_fp_kind,2.50881800e+00_fp_kind,6.37853000e-01_fp_kind,7.22838000e-01_fp_kind,1.14275600e+00_fp_kind,&
               1.41822590e+01_fp_kind,5.93685800e+00_fp_kind,1.12726000e-01_fp_kind,3.49584810e+01_fp_kind,3.90240000e-01_fp_kind,&
               2.70140000e-02_fp_kind,8.00000000e+00_fp_kind,1.71987724e-01_fp_kind,&
   9.0_fp_kind,3.51195400e+00_fp_kind,2.77224400e+00_fp_kind,6.78385000e-01_fp_kind,9.15159000e-01_fp_kind,1.08926100e+00_fp_kind,&
               1.06878590e+01_fp_kind,4.38046600e+00_fp_kind,9.39820000e-02_fp_kind,2.72552030e+01_fp_kind,3.13066000e-01_fp_kind,&
               3.25570000e-02_fp_kind,9.00000000e+00_fp_kind,2.23745514e-01_fp_kind,&
  10.00_fp_kind,4.18374900e+00_fp_kind,2.90572600e+00_fp_kind,5.20513000e-01_fp_kind,1.13564100e+00_fp_kind,1.22806500e+00_fp_kind,&
               8.17545700e+00_fp_kind,3.25253600e+00_fp_kind,6.32950000e-02_fp_kind,2.18139090e+01_fp_kind,2.24952000e-01_fp_kind,&
               2.55760000e-02_fp_kind,1.00000000e+01_fp_kind,2.87627317e-01_fp_kind,&
  11.00_fp_kind,4.91012700e+00_fp_kind,3.08178300e+00_fp_kind,1.26206700e+00_fp_kind,1.09893800e+00_fp_kind,5.60991000e-01_fp_kind,&
               3.28143400e+00_fp_kind,9.11917800e+00_fp_kind,1.02763000e-01_fp_kind,1.32013942e+02_fp_kind,4.05878000e-01_fp_kind,&
               7.97120000e-02_fp_kind,1.10000000e+01_fp_kind,3.66747421e-01_fp_kind,&
  12.00_fp_kind,4.70897100e+00_fp_kind,1.19481400e+00_fp_kind,1.55815700e+00_fp_kind,1.17041300e+00_fp_kind,3.23940300e+00_fp_kind,&
               4.87520700e+00_fp_kind,1.08506079e+02_fp_kind,1.11516000e-01_fp_kind,4.82924070e+01_fp_kind,1.92817100e+00_fp_kind,&
               1.26842000e-01_fp_kind,1.20000000e+01_fp_kind,4.70985328e-01_fp_kind,&
  13.00_fp_kind,4.73079600e+00_fp_kind,2.31395100e+00_fp_kind,1.54198000e+00_fp_kind,1.11756400e+00_fp_kind,3.15475400e+00_fp_kind,&
               3.62893100e+00_fp_kind,4.30511660e+01_fp_kind,9.59600000e-02_fp_kind,1.08932389e+02_fp_kind,1.55591800e+00_fp_kind,&
               1.39509000e-01_fp_kind,1.30000000e+01_fp_kind,5.28931702e-01_fp_kind,&
  14.00_fp_kind,5.27532900e+00_fp_kind,3.19103800e+00_fp_kind,1.51151400e+00_fp_kind,1.35684900e+00_fp_kind,2.51911400e+00_fp_kind,&
               2.63133800e+00_fp_kind,3.37307280e+01_fp_kind,8.11190000e-02_fp_kind,8.62886400e+01_fp_kind,1.17008700e+00_fp_kind,&
               1.45073000e-01_fp_kind,1.40000000e+01_fp_kind,5.92196011e-01_fp_kind,&
  15.00_fp_kind,1.95054100e+00_fp_kind,4.14693000e+00_fp_kind,1.49456000e+00_fp_kind,1.52204200e+00_fp_kind,5.72971100e+00_fp_kind,&
               9.08139000e-01_fp_kind,2.70449530e+01_fp_kind,7.12800000e-02_fp_kind,6.75201900e+01_fp_kind,1.98117300e+00_fp_kind,&
               1.55233000e-01_fp_kind,1.50000000e+01_fp_kind,6.60050129e-01_fp_kind,&
  16.00_fp_kind,6.37215700e+00_fp_kind,5.15456800e+00_fp_kind,1.47373200e+00_fp_kind,1.63507300e+00_fp_kind,1.20937200e+00_fp_kind,&
               1.51434700e+00_fp_kind,2.20925280e+01_fp_kind,6.13730000e-02_fp_kind,5.54451760e+01_fp_kind,6.46925000e-01_fp_kind,&
               1.54722000e-01_fp_kind,1.60000000e+01_fp_kind,7.26458198e-01_fp_kind,&
  17.00_fp_kind,1.44607100e+00_fp_kind,6.87060900e+00_fp_kind,6.15180100e+00_fp_kind,1.75034700e+00_fp_kind,6.34168000e-01_fp_kind,&
               5.23570000e-02_fp_kind,1.19316500e+00_fp_kind,1.83434160e+01_fp_kind,4.63983940e+01_fp_kind,4.01005000e-01_fp_kind,&
               1.46773000e-01_fp_kind,1.70000000e+01_fp_kind,7.92912007e-01_fp_kind,&
  18.00_fp_kind,7.18800400e+00_fp_kind,6.63845400e+00_fp_kind,4.54180000e-01_fp_kind,1.92959300e+00_fp_kind,1.52365400e+00_fp_kind,&
               9.56221000e-01_fp_kind,1.53398770e+01_fp_kind,1.53398620e+01_fp_kind,3.90438240e+01_fp_kind,6.24090000e-02_fp_kind,&
               2.65954000e-01_fp_kind,1.80000000e+01_fp_kind,8.74904036e-01_fp_kind,&
  19.00_fp_kind,8.16399100e+00_fp_kind,7.14694500e+00_fp_kind,1.07014000e+00_fp_kind,8.77316000e-01_fp_kind,1.48643400e+00_fp_kind,&
               1.28163230e+01_fp_kind,8.08945000e-01_fp_kind,2.10327009e+02_fp_kind,3.95976510e+01_fp_kind,5.28210000e-02_fp_kind,&
               2.53614000e-01_fp_kind,1.90000000e+01_fp_kind,9.24257331e-01_fp_kind,&
  20.00_fp_kind,8.59365500e+00_fp_kind,1.47732400e+00_fp_kind,1.43625400e+00_fp_kind,1.18283900e+00_fp_kind,7.11325800e+00_fp_kind,&
               1.04606440e+01_fp_kind,4.18910000e-02_fp_kind,8.13903820e+01_fp_kind,1.69847839e+02_fp_kind,6.88098000e-01_fp_kind,&
               1.96255000e-01_fp_kind,2.00000000e+01_fp_kind,9.67132503e-01_fp_kind,&
  21.00_fp_kind,1.47656600e+00_fp_kind,1.48727800e+00_fp_kind,1.60018700e+00_fp_kind,9.17746300e+00_fp_kind,7.09975000e+00_fp_kind,&
               5.31310220e+01_fp_kind,3.53250000e-02_fp_kind,1.37319495e+02_fp_kind,9.09803100e+00_fp_kind,6.02102000e-01_fp_kind,&
               1.57765000e-01_fp_kind,2.10000000e+01_fp_kind,1.01297507e+00_fp_kind,&
  22.00_fp_kind,9.81852400e+00_fp_kind,1.52264600e+00_fp_kind,1.70310100e+00_fp_kind,1.76877400e+00_fp_kind,7.08255500e+00_fp_kind,&
               8.00187900e+00_fp_kind,2.97630000e-02_fp_kind,3.98854230e+01_fp_kind,1.20158000e+02_fp_kind,5.32405000e-01_fp_kind,&
               1.02473000e-01_fp_kind,2.20000000e+01_fp_kind,1.05087643e+00_fp_kind,&
  23.00_fp_kind,1.04735750e+01_fp_kind,1.54788100e+00_fp_kind,1.98638100e+00_fp_kind,1.86561600e+00_fp_kind,7.05625000e+00_fp_kind,&
               7.08194000e+00_fp_kind,2.60400000e-02_fp_kind,3.19096720e+01_fp_kind,1.08022844e+02_fp_kind,4.74882000e-01_fp_kind,&
               6.77440000e-02_fp_kind,2.30000000e+01_fp_kind,1.08670420e+00_fp_kind,&
  24.00_fp_kind,1.10070690e+01_fp_kind,1.55547700e+00_fp_kind,2.98529300e+00_fp_kind,1.34785500e+00_fp_kind,7.03477900e+00_fp_kind,&
               6.36628100e+00_fp_kind,2.39870000e-02_fp_kind,2.32448380e+01_fp_kind,1.05774500e+02_fp_kind,4.29369000e-01_fp_kind,&
               6.55100000e-02_fp_kind,2.40000000e+01_fp_kind,1.11564928e+00_fp_kind,&
  25.00_fp_kind,1.17095420e+01_fp_kind,1.73341400e+00_fp_kind,2.67314100e+00_fp_kind,2.02336800e+00_fp_kind,7.00318000e+00_fp_kind,&
               5.59712000e+00_fp_kind,1.78000000e-02_fp_kind,2.17884190e+01_fp_kind,8.95179150e+01_fp_kind,3.83054000e-01_fp_kind,&
               -1.47293000e-01_fp_kind,2.50000000e+01_fp_kind,1.13790440e+00_fp_kind,&
  26.00_fp_kind,1.23110980e+01_fp_kind,1.87662300e+00_fp_kind,3.06617700e+00_fp_kind,2.07045100e+00_fp_kind,6.97518500e+00_fp_kind,&
               5.00941500e+00_fp_kind,1.44610000e-02_fp_kind,1.87430410e+01_fp_kind,8.27678740e+01_fp_kind,3.46506000e-01_fp_kind,&
               -3.04931000e-01_fp_kind,2.60000000e+01_fp_kind,1.15778197e+00_fp_kind,&
  27.00_fp_kind,1.29145100e+01_fp_kind,2.48190800e+00_fp_kind,3.46689400e+00_fp_kind,2.10635100e+00_fp_kind,6.96089200e+00_fp_kind,&
               4.50713800e+00_fp_kind,9.12600000e-03_fp_kind,1.64381300e+01_fp_kind,7.69873170e+01_fp_kind,3.14418000e-01_fp_kind,&
               -9.36572000e-01_fp_kind,2.70000000e+01_fp_kind,1.17078217e+00_fp_kind,&
  28.00_fp_kind,1.35218650e+01_fp_kind,6.94728500e+00_fp_kind,3.86602800e+00_fp_kind,2.13590000e+00_fp_kind,4.28473100e+00_fp_kind,&
               4.07727700e+00_fp_kind,2.86763000e-01_fp_kind,1.46226340e+01_fp_kind,7.19660780e+01_fp_kind,4.43700000e-03_fp_kind,&
               -2.76269700e+00_fp_kind,2.80000000e+01_fp_kind,1.18146331e+00_fp_kind,&
  29.00_fp_kind,1.40141920e+01_fp_kind,4.78457700e+00_fp_kind,5.05680600e+00_fp_kind,1.45797100e+00_fp_kind,6.93299600e+00_fp_kind,&
               3.73828000e+00_fp_kind,3.74400000e-03_fp_kind,1.30349820e+01_fp_kind,7.25547930e+01_fp_kind,2.65666000e-01_fp_kind,&
               -3.25447700e+00_fp_kind,2.90000000e+01_fp_kind,1.18912072e+00_fp_kind,&
  30.00_fp_kind,1.47410020e+01_fp_kind,6.90774800e+00_fp_kind,4.64233700e+00_fp_kind,2.19176600e+00_fp_kind,3.84240420e+01_fp_kind,&
               3.38823200e+00_fp_kind,2.43315000e-01_fp_kind,1.19036890e+01_fp_kind,6.33121300e+01_fp_kind,3.97000000e-04_fp_kind,&
               -3.69158280e+01_fp_kind,3.00000000e+01_fp_kind,1.19526702e+00_fp_kind,&
  31.00_fp_kind,1.57589460e+01_fp_kind,6.84112300e+00_fp_kind,4.12101600e+00_fp_kind,2.71468100e+00_fp_kind,2.39524600e+00_fp_kind,&
               3.12175400e+00_fp_kind,2.26057000e-01_fp_kind,1.24821960e+01_fp_kind,6.62036220e+01_fp_kind,7.23800000e-03_fp_kind,&
               -8.47395000e-01_fp_kind,3.10000000e+01_fp_kind,1.20051278e+00_fp_kind,&
  32.00_fp_kind,1.65406140e+01_fp_kind,1.56790000e+00_fp_kind,3.72782900e+00_fp_kind,3.34509800e+00_fp_kind,6.78507900e+00_fp_kind,&
               2.86661800e+00_fp_kind,1.21980000e-02_fp_kind,1.34321630e+01_fp_kind,5.88660460e+01_fp_kind,2.10974000e-01_fp_kind,&
               1.87260000e-02_fp_kind,3.20000000e+01_fp_kind,1.20065585e+00_fp_kind,&
  33.00_fp_kind,1.70256430e+01_fp_kind,4.50344100e+00_fp_kind,3.71590400e+00_fp_kind,3.93720000e+00_fp_kind,6.79017500e+00_fp_kind,&
               2.59773900e+00_fp_kind,3.01200000e-03_fp_kind,1.42721190e+01_fp_kind,5.04379970e+01_fp_kind,1.93015000e-01_fp_kind,&
               -2.98411700e+00_fp_kind,3.30000000e+01_fp_kind,1.19831414e+00_fp_kind,&
  34.00_fp_kind,1.73540710e+01_fp_kind,4.65324800e+00_fp_kind,4.25948900e+00_fp_kind,4.13645500e+00_fp_kind,6.74916300e+00_fp_kind,&
               2.34978700e+00_fp_kind,2.55000000e-03_fp_kind,1.55794600e+01_fp_kind,4.51812010e+01_fp_kind,1.77432000e-01_fp_kind,&
               -3.16098200e+00_fp_kind,3.40000000e+01_fp_kind,1.19852893e+00_fp_kind,&
  35.00_fp_kind,1.75505700e+01_fp_kind,5.41188200e+00_fp_kind,3.93718000e+00_fp_kind,3.88064500e+00_fp_kind,6.70779300e+00_fp_kind,&
               2.11922600e+00_fp_kind,1.65571850e+01_fp_kind,2.48100000e-03_fp_kind,4.21640090e+01_fp_kind,1.62121000e-01_fp_kind,&
               -2.49208800e+00_fp_kind,3.50000000e+01_fp_kind,1.19916154e+00_fp_kind,&
  36.00_fp_kind,1.76552790e+01_fp_kind,6.84810500e+00_fp_kind,4.17100400e+00_fp_kind,3.44676000e+00_fp_kind,6.68520000e+00_fp_kind,&
               1.90823100e+00_fp_kind,1.66062350e+01_fp_kind,1.59800000e-03_fp_kind,3.99174710e+01_fp_kind,1.46896000e-01_fp_kind,&
               -2.81059200e+00_fp_kind,3.60000000e+01_fp_kind,1.19968779e+00_fp_kind,&
  37.00_fp_kind,8.12313400e+00_fp_kind,2.13804200e+00_fp_kind,6.76170200e+00_fp_kind,1.15605100e+00_fp_kind,1.76795470e+01_fp_kind,&
               1.51423850e+01_fp_kind,3.35426660e+01_fp_kind,1.29372000e-01_fp_kind,2.24132506e+02_fp_kind,1.71336800e+00_fp_kind,&
               1.13954800e+00_fp_kind,3.70000000e+01_fp_kind,1.21057522e+00_fp_kind,&
  38.00_fp_kind,1.77302190e+01_fp_kind,9.79586700e+00_fp_kind,6.09976300e+00_fp_kind,2.62002500e+00_fp_kind,6.00053000e-01_fp_kind,&
               1.56306000e+00_fp_kind,1.43108680e+01_fp_kind,1.20574000e-01_fp_kind,1.35771318e+02_fp_kind,1.20574000e-01_fp_kind,&
               1.14025100e+00_fp_kind,3.80000000e+01_fp_kind,1.20174901e+00_fp_kind,&
  39.00_fp_kind,1.77920400e+01_fp_kind,1.02532520e+01_fp_kind,5.71494900e+00_fp_kind,3.17051600e+00_fp_kind,9.18251000e-01_fp_kind,&
               1.42969100e+00_fp_kind,1.31328160e+01_fp_kind,1.12173000e-01_fp_kind,1.08197029e+02_fp_kind,1.12173000e-01_fp_kind,&
               1.31787000e+00_fp_kind,3.90000000e+01_fp_kind,1.37502186e+00_fp_kind,&
  40.00_fp_kind,1.78597710e+01_fp_kind,1.09110380e+01_fp_kind,5.82111500e+00_fp_kind,3.51251300e+00_fp_kind,7.46965000e-01_fp_kind,&
               1.31069200e+00_fp_kind,1.23192850e+01_fp_kind,1.04353000e-01_fp_kind,9.17775440e+01_fp_kind,1.04353000e-01_fp_kind,&
               1.12485900e+00_fp_kind,4.00000000e+01_fp_kind,1.18844292e+00_fp_kind,&
  41.00_fp_kind,1.79583980e+01_fp_kind,1.20630540e+01_fp_kind,5.00701500e+00_fp_kind,3.28766700e+00_fp_kind,1.53101900e+00_fp_kind,&
               1.21159000e+00_fp_kind,1.22466870e+01_fp_kind,9.86150000e-02_fp_kind,7.50119440e+01_fp_kind,9.86150000e-02_fp_kind,&
               1.12345200e+00_fp_kind,4.10000000e+01_fp_kind,1.18935960e+00_fp_kind,&
  42.00_fp_kind,6.23621800e+00_fp_kind,1.79877110e+01_fp_kind,1.29731270e+01_fp_kind,3.45142600e+00_fp_kind,2.10899000e-01_fp_kind,&
               9.07800000e-02_fp_kind,1.10831000e+00_fp_kind,1.14687200e+01_fp_kind,6.66841530e+01_fp_kind,9.07800000e-02_fp_kind,&
               1.10877000e+00_fp_kind,4.20000000e+01_fp_kind,1.19948206e+00_fp_kind,&
  43.00_fp_kind,1.78409640e+01_fp_kind,3.42823600e+00_fp_kind,1.37301200e+00_fp_kind,1.29473640e+01_fp_kind,6.33546900e+00_fp_kind,&
               1.00572900e+00_fp_kind,4.19013830e+01_fp_kind,1.19320541e+02_fp_kind,9.78154200e+00_fp_kind,8.33910000e-02_fp_kind,&
               1.07478400e+00_fp_kind,4.30000000e+01_fp_kind,1.20219703e+00_fp_kind,&
  44.00_fp_kind,6.27162400e+00_fp_kind,1.79067390e+01_fp_kind,1.41232690e+01_fp_kind,3.74600800e+00_fp_kind,9.08235000e-01_fp_kind,&
               7.70400000e-02_fp_kind,9.28222000e-01_fp_kind,9.55534500e+00_fp_kind,3.58606780e+01_fp_kind,1.23552247e+02_fp_kind,&
               1.04399200e+00_fp_kind,4.40000000e+01_fp_kind,1.21422835e+00_fp_kind,&
  45.00_fp_kind,6.21664800e+00_fp_kind,1.79197380e+01_fp_kind,3.85425200e+00_fp_kind,8.40326000e-01_fp_kind,1.51734980e+01_fp_kind,&
               7.07890000e-02_fp_kind,8.56121000e-01_fp_kind,3.38894840e+01_fp_kind,1.21686688e+02_fp_kind,9.02951700e+00_fp_kind,&
               9.95452000e-01_fp_kind,4.50000000e+01_fp_kind,1.22566200e+00_fp_kind,&
  46.00_fp_kind,6.12151100e+00_fp_kind,4.78406300e+00_fp_kind,1.66316830e+01_fp_kind,4.31825800e+00_fp_kind,1.32467730e+01_fp_kind,&
               6.25490000e-02_fp_kind,7.84031000e-01_fp_kind,8.75139100e+00_fp_kind,3.44899830e+01_fp_kind,7.84031000e-01_fp_kind,&
               8.83099000e-01_fp_kind,4.60000000e+01_fp_kind,1.23621345e+00_fp_kind,&
  47.00_fp_kind,6.07387400e+00_fp_kind,1.71554370e+01_fp_kind,4.17334400e+00_fp_kind,8.52238000e-01_fp_kind,1.79886850e+01_fp_kind,&
               5.53330000e-02_fp_kind,7.89651200e+00_fp_kind,2.84437390e+01_fp_kind,1.10376108e+02_fp_kind,7.16809000e-01_fp_kind,&
               7.56603000e-01_fp_kind,4.70000000e+01_fp_kind,1.25659754e+00_fp_kind,&
  48.00_fp_kind,6.08098600e+00_fp_kind,1.80194680e+01_fp_kind,4.01819700e+00_fp_kind,1.30351000e+00_fp_kind,1.79746690e+01_fp_kind,&
               4.89900000e-02_fp_kind,7.27364600e+00_fp_kind,2.91192830e+01_fp_kind,9.58312080e+01_fp_kind,6.61231000e-01_fp_kind,&
               6.03504000e-01_fp_kind,4.80000000e+01_fp_kind,1.27826001e+00_fp_kind,&
  49.00_fp_kind,6.19647700e+00_fp_kind,1.88161830e+01_fp_kind,4.05047900e+00_fp_kind,1.63892900e+00_fp_kind,1.79629120e+01_fp_kind,&
               4.20720000e-02_fp_kind,6.69566500e+00_fp_kind,3.10097910e+01_fp_kind,1.03284350e+02_fp_kind,6.10714000e-01_fp_kind,&
               3.33097000e-01_fp_kind,4.90000000e+01_fp_kind,1.29047036e+00_fp_kind,&
  50.00_fp_kind,1.93251710e+01_fp_kind,6.28157100e+00_fp_kind,4.49886600e+00_fp_kind,1.85693400e+00_fp_kind,1.79173180e+01_fp_kind,&
               6.11810400e+00_fp_kind,3.69150000e-02_fp_kind,3.25290470e+01_fp_kind,9.50371820e+01_fp_kind,5.65651000e-01_fp_kind,&
               1.19024000e-01_fp_kind,5.00000000e+01_fp_kind,1.33054938e+00_fp_kind,&
  51.00_fp_kind,5.39495600e+00_fp_kind,6.54957000e+00_fp_kind,1.96506810e+01_fp_kind,1.82782000e+00_fp_kind,1.78678330e+01_fp_kind,&
               3.33265230e+01_fp_kind,3.09740000e-02_fp_kind,5.56492900e+00_fp_kind,8.71309650e+01_fp_kind,5.23992000e-01_fp_kind,&
               -2.90506000e-01_fp_kind,5.10000000e+01_fp_kind,1.36040269e+00_fp_kind,&
  52.00_fp_kind,6.66030200e+00_fp_kind,6.94075600e+00_fp_kind,1.98470150e+01_fp_kind,1.55717500e+00_fp_kind,1.78024270e+01_fp_kind,&
               3.30316560e+01_fp_kind,2.57500000e-02_fp_kind,5.06554700e+00_fp_kind,8.41016130e+01_fp_kind,4.87660000e-01_fp_kind,&
               -8.06668000e-01_fp_kind,5.20000000e+01_fp_kind,1.39517164e+00_fp_kind,&
  53.00_fp_kind,1.98845020e+01_fp_kind,6.73659300e+00_fp_kind,8.11051600e+00_fp_kind,1.17095300e+00_fp_kind,1.75487150e+01_fp_kind,&
               4.62859100e+00_fp_kind,2.77540000e-02_fp_kind,3.18490960e+01_fp_kind,8.44063910e+01_fp_kind,4.63550000e-01_fp_kind,&
               -4.48811000e-01_fp_kind,5.30000000e+01_fp_kind,1.43493499e+00_fp_kind,&
  54.00_fp_kind,1.99789200e+01_fp_kind,1.17749450e+01_fp_kind,9.33218200e+00_fp_kind,1.24474900e+00_fp_kind,1.77375010e+01_fp_kind,&
               4.14335600e+00_fp_kind,1.01420000e-02_fp_kind,2.87961990e+01_fp_kind,7.52806880e+01_fp_kind,4.13616000e-01_fp_kind,&
               -6.06590200e+00_fp_kind,5.40000000e+01_fp_kind,1.46192815e+00_fp_kind,&
  55.00_fp_kind,1.74186750e+01_fp_kind,8.31444400e+00_fp_kind,1.03231930e+01_fp_kind,1.38383400e+00_fp_kind,1.98762520e+01_fp_kind,&
               3.99828000e-01_fp_kind,1.68720000e-02_fp_kind,2.56058280e+01_fp_kind,2.33339674e+02_fp_kind,3.82691500e+00_fp_kind,&
               -2.32280200e+00_fp_kind,5.50000000e+01_fp_kind,1.50473467e+00_fp_kind,&
  56.00_fp_kind,1.97473440e+01_fp_kind,1.73684760e+01_fp_kind,1.04657180e+01_fp_kind,2.59260200e+00_fp_kind,1.10036530e+01_fp_kind,&
               3.48182300e+00_fp_kind,3.71224000e-01_fp_kind,2.12266410e+01_fp_kind,1.73834271e+02_fp_kind,1.07190000e-02_fp_kind,&
               -5.18349700e+00_fp_kind,5.60000000e+01_fp_kind,1.54005179e+00_fp_kind,&
  57.00_fp_kind,1.99660180e+01_fp_kind,2.73296540e+01_fp_kind,1.10184250e+01_fp_kind,3.08669600e+00_fp_kind,1.73354540e+01_fp_kind,&
               3.19740800e+00_fp_kind,3.44600000e-03_fp_kind,1.99554920e+01_fp_kind,1.41381979e+02_fp_kind,3.41817000e-01_fp_kind,&
               -2.17454890e+01_fp_kind,5.70000000e+01_fp_kind,1.57946005e+00_fp_kind,&
  58.00_fp_kind,1.73551210e+01_fp_kind,4.39884980e+01_fp_kind,2.05466500e+01_fp_kind,3.13067000e+00_fp_kind,1.13536650e+01_fp_kind,&
               3.28369000e-01_fp_kind,2.04700000e-03_fp_kind,3.08819600e+00_fp_kind,1.34907661e+02_fp_kind,1.88329610e+01_fp_kind,&
               -3.83860170e+01_fp_kind,5.80000000e+01_fp_kind,1.60642148e+00_fp_kind,&
  59.00_fp_kind,2.15513110e+01_fp_kind,1.71617290e+01_fp_kind,1.19038590e+01_fp_kind,2.67910300e+00_fp_kind,9.56419700e+00_fp_kind,&
               2.99567500e+00_fp_kind,3.12491000e-01_fp_kind,1.77167050e+01_fp_kind,1.52192827e+02_fp_kind,1.04680000e-02_fp_kind,&
               -3.87106800e+00_fp_kind,5.90000000e+01_fp_kind,1.72002615e+00_fp_kind,&
  60.00_fp_kind,1.73312440e+01_fp_kind,6.27839230e+01_fp_kind,1.21600970e+01_fp_kind,2.66348300e+00_fp_kind,2.22399510e+01_fp_kind,&
               3.00269000e-01_fp_kind,1.32000000e-03_fp_kind,1.70260010e+01_fp_kind,1.48748986e+02_fp_kind,2.91026800e+00_fp_kind,&
               -5.71898440e+01_fp_kind,6.00000000e+01_fp_kind,1.68365235e+00_fp_kind,&
  61.00_fp_kind,1.72863880e+01_fp_kind,5.15601610e+01_fp_kind,1.24785570e+01_fp_kind,2.67551500e+00_fp_kind,2.29609470e+01_fp_kind,&
               2.86620000e-01_fp_kind,1.55000000e-03_fp_kind,1.62237550e+01_fp_kind,1.43984513e+02_fp_kind,2.79648000e+00_fp_kind,&
               -4.59736810e+01_fp_kind,6.10000000e+01_fp_kind,1.72469340e+00_fp_kind,&
  62.00_fp_kind,2.37003640e+01_fp_kind,2.30722150e+01_fp_kind,1.27777820e+01_fp_kind,2.68421700e+00_fp_kind,1.72043660e+01_fp_kind,&
               2.68953900e+00_fp_kind,3.49100000e-03_fp_kind,1.54954370e+01_fp_kind,1.39862475e+02_fp_kind,2.74536000e-01_fp_kind,&
               -1.74521660e+01_fp_kind,6.20000000e+01_fp_kind,1.76401216e+00_fp_kind,&
  63.00_fp_kind,1.71861950e+01_fp_kind,3.71568390e+01_fp_kind,1.31033870e+01_fp_kind,2.70724600e+00_fp_kind,2.44192710e+01_fp_kind,&
               2.61678000e-01_fp_kind,1.99500000e-03_fp_kind,1.47873600e+01_fp_kind,1.34816293e+02_fp_kind,2.58188300e+00_fp_kind,&
               -3.15866870e+01_fp_kind,6.30000000e+01_fp_kind,1.79780574e+00_fp_kind,&
  64.00_fp_kind,2.48981180e+01_fp_kind,1.71049510e+01_fp_kind,1.32225810e+01_fp_kind,3.26615200e+00_fp_kind,4.89952140e+01_fp_kind,&
               2.43502800e+00_fp_kind,2.46961000e-01_fp_kind,1.39963250e+01_fp_kind,1.10863093e+02_fp_kind,1.38300000e-03_fp_kind,&
               -4.35056840e+01_fp_kind,6.40000000e+01_fp_kind,1.84011864e+00_fp_kind,&
  65.00_fp_kind,2.59100130e+01_fp_kind,3.23441390e+01_fp_kind,1.37651170e+01_fp_kind,2.75140400e+00_fp_kind,1.70644050e+01_fp_kind,&
               2.37391200e+00_fp_kind,2.03400000e-03_fp_kind,1.34819690e+01_fp_kind,1.25836511e+02_fp_kind,2.36916000e-01_fp_kind,&
               -2.68519700e+01_fp_kind,6.50000000e+01_fp_kind,1.87131310e+00_fp_kind,&
  66.00_fp_kind,2.66717850e+01_fp_kind,8.86875770e+01_fp_kind,1.40654450e+01_fp_kind,2.76849700e+00_fp_kind,1.70677820e+01_fp_kind,&
               2.28259300e+00_fp_kind,6.65000000e-04_fp_kind,1.29202300e+01_fp_kind,1.21937188e+02_fp_kind,2.25531000e-01_fp_kind,&
               -8.32798310e+01_fp_kind,6.60000000e+01_fp_kind,1.90372064e+00_fp_kind,&
  67.00_fp_kind,2.71501900e+01_fp_kind,1.69998190e+01_fp_kind,1.40593340e+01_fp_kind,3.38697900e+00_fp_kind,4.65464710e+01_fp_kind,&
               2.16966000e+00_fp_kind,2.15414000e-01_fp_kind,1.22131480e+01_fp_kind,1.00506781e+02_fp_kind,1.21100000e-03_fp_kind,&
               -4.11652530e+01_fp_kind,6.70000000e+01_fp_kind,1.92623516e+00_fp_kind,&
  68.00_fp_kind,2.81748860e+01_fp_kind,8.24932690e+01_fp_kind,1.46240020e+01_fp_kind,2.80275600e+00_fp_kind,1.70185150e+01_fp_kind,&
               2.12099500e+00_fp_kind,6.40000000e-04_fp_kind,1.19152560e+01_fp_kind,1.14529936e+02_fp_kind,2.07519000e-01_fp_kind,&
               -7.71352210e+01_fp_kind,6.80000000e+01_fp_kind,1.94691327e+00_fp_kind,&
  69.00_fp_kind,2.89258940e+01_fp_kind,7.61737960e+01_fp_kind,1.49047040e+01_fp_kind,2.81481200e+00_fp_kind,1.69981170e+01_fp_kind,&
               2.04620300e+00_fp_kind,6.56000000e-04_fp_kind,1.14653750e+01_fp_kind,1.11411979e+02_fp_kind,1.99376000e-01_fp_kind,&
               -7.08398130e+01_fp_kind,6.90000000e+01_fp_kind,1.96372421e+00_fp_kind,&
  70.00_fp_kind,2.96767600e+01_fp_kind,6.56240680e+01_fp_kind,1.51608540e+01_fp_kind,2.83028800e+00_fp_kind,1.69978500e+01_fp_kind,&
               1.97763000e+00_fp_kind,7.20000000e-04_fp_kind,1.10446220e+01_fp_kind,1.08139150e+02_fp_kind,1.92110000e-01_fp_kind,&
               -6.03138120e+01_fp_kind,7.00000000e+01_fp_kind,1.97925898e+00_fp_kind,&
  71.00_fp_kind,3.01228650e+01_fp_kind,1.50993460e+01_fp_kind,5.63148990e+01_fp_kind,3.54098000e+00_fp_kind,1.69437300e+01_fp_kind,&
               1.88309000e+00_fp_kind,1.03427640e+01_fp_kind,7.80000000e-04_fp_kind,8.95592480e+01_fp_kind,1.83849000e-01_fp_kind,&
               -5.10494170e+01_fp_kind,7.10000000e+01_fp_kind,1.99544635e+00_fp_kind,&
  72.00_fp_kind,3.06170330e+01_fp_kind,1.51453510e+01_fp_kind,5.49335480e+01_fp_kind,4.09625300e+00_fp_kind,1.68961570e+01_fp_kind,&
               1.79561300e+00_fp_kind,9.93446900e+00_fp_kind,7.39000000e-04_fp_kind,7.61897070e+01_fp_kind,1.75914000e-01_fp_kind,&
               -4.97198380e+01_fp_kind,7.20000000e+01_fp_kind,2.00672838e+00_fp_kind,&
  73.00_fp_kind,3.10663580e+01_fp_kind,1.53418230e+01_fp_kind,4.92782960e+01_fp_kind,4.57766500e+00_fp_kind,1.68283210e+01_fp_kind,&
               1.70873200e+00_fp_kind,9.61845500e+00_fp_kind,7.60000000e-04_fp_kind,6.63462020e+01_fp_kind,1.68002000e-01_fp_kind,&
               -4.41190250e+01_fp_kind,7.30000000e+01_fp_kind,2.01481064e+00_fp_kind,&
  74.00_fp_kind,3.15079010e+01_fp_kind,1.56824980e+01_fp_kind,3.79601270e+01_fp_kind,4.88550900e+00_fp_kind,1.67921130e+01_fp_kind,&
               1.62948500e+00_fp_kind,9.44644800e+00_fp_kind,8.98000000e-04_fp_kind,5.99806750e+01_fp_kind,1.60798000e-01_fp_kind,&
               -3.28645760e+01_fp_kind,7.40000000e+01_fp_kind,2.02423309e+00_fp_kind,&
  75.00_fp_kind,3.18884560e+01_fp_kind,1.61171030e+01_fp_kind,4.23902960e+01_fp_kind,5.21166900e+00_fp_kind,1.67675910e+01_fp_kind,&
               1.54923800e+00_fp_kind,9.23347400e+00_fp_kind,6.89000000e-04_fp_kind,5.45163710e+01_fp_kind,1.52815000e-01_fp_kind,&
               -3.74126810e+01_fp_kind,7.50000000e+01_fp_kind,2.03220395e+00_fp_kind,&
  76.00_fp_kind,3.22102980e+01_fp_kind,1.66784400e+01_fp_kind,4.85599070e+01_fp_kind,5.45583900e+00_fp_kind,1.67355320e+01_fp_kind,&
               1.47353100e+00_fp_kind,9.04969500e+00_fp_kind,5.19000000e-04_fp_kind,5.02102010e+01_fp_kind,1.45771000e-01_fp_kind,&
               -4.36779540e+01_fp_kind,7.60000000e+01_fp_kind,2.03756230e+00_fp_kind,&
  77.00_fp_kind,3.20044370e+01_fp_kind,1.97545400e+00_fp_kind,1.70701040e+01_fp_kind,1.59394540e+01_fp_kind,5.99000300e+00_fp_kind,&
               1.35376700e+00_fp_kind,8.10141720e+01_fp_kind,1.28093000e-01_fp_kind,7.66119600e+00_fp_kind,2.66594030e+01_fp_kind,&
               4.01889300e+00_fp_kind,7.70000000e+01_fp_kind,2.07093147e+00_fp_kind,&
  78.00_fp_kind,3.12738910e+01_fp_kind,1.84454410e+01_fp_kind,1.70637450e+01_fp_kind,5.55593300e+00_fp_kind,1.57527000e+00_fp_kind,&
               1.31699200e+00_fp_kind,8.79715400e+00_fp_kind,1.24741000e-01_fp_kind,4.01779940e+01_fp_kind,1.31699700e+00_fp_kind,&
               4.05039400e+00_fp_kind,7.80000000e+01_fp_kind,2.07030766e+00_fp_kind,&
  79.00_fp_kind,1.67773890e+01_fp_kind,1.93171560e+01_fp_kind,3.29796820e+01_fp_kind,5.59545300e+00_fp_kind,1.05768540e+01_fp_kind,&
               1.22737000e-01_fp_kind,8.62157000e+00_fp_kind,1.25690200e+00_fp_kind,3.80088210e+01_fp_kind,6.01000000e-04_fp_kind,&
               -6.27907800e+00_fp_kind,7.90000000e+01_fp_kind,2.05883489e+00_fp_kind,&
  80.00_fp_kind,1.68398890e+01_fp_kind,2.00238230e+01_fp_kind,2.84285650e+01_fp_kind,5.88156400e+00_fp_kind,4.71470600e+00_fp_kind,&
               1.15905000e-01_fp_kind,8.25692700e+00_fp_kind,1.19525000e+00_fp_kind,3.92472260e+01_fp_kind,1.19525000e+00_fp_kind,&
               4.07647800e+00_fp_kind,8.00000000e+01_fp_kind,2.06302643e+00_fp_kind,&
  81.00_fp_kind,1.66307950e+01_fp_kind,1.93866150e+01_fp_kind,3.28085700e+01_fp_kind,1.74719100e+00_fp_kind,6.35686200e+00_fp_kind,&
               1.10704000e-01_fp_kind,7.18140100e+00_fp_kind,1.11973000e+00_fp_kind,9.06602620e+01_fp_kind,2.60149780e+01_fp_kind,&
               4.06693900e+00_fp_kind,8.10000000e+01_fp_kind,2.05598134e+00_fp_kind,&
  82.00_fp_kind,1.64195670e+01_fp_kind,3.27385920e+01_fp_kind,6.53024700e+00_fp_kind,2.34274200e+00_fp_kind,1.99164750e+01_fp_kind,&
               1.05499000e-01_fp_kind,1.05504900e+00_fp_kind,2.50258900e+01_fp_kind,8.09065960e+01_fp_kind,6.66444900e+00_fp_kind,&
               4.04982400e+00_fp_kind,8.20000000e+01_fp_kind,2.05000970e+00_fp_kind,&
  83.00_fp_kind,1.62822740e+01_fp_kind,3.27251370e+01_fp_kind,6.67830200e+00_fp_kind,2.69475000e+00_fp_kind,2.05765590e+01_fp_kind,&
               1.01180000e-01_fp_kind,1.00228700e+00_fp_kind,2.57141450e+01_fp_kind,7.70575500e+01_fp_kind,6.29188200e+00_fp_kind,&
               4.04091400e+00_fp_kind,8.30000000e+01_fp_kind,2.04784654e+00_fp_kind,&
  84.00_fp_kind,1.62891640e+01_fp_kind,3.28071700e+01_fp_kind,2.10951640e+01_fp_kind,2.50590100e+00_fp_kind,7.25458900e+00_fp_kind,&
               9.81210000e-02_fp_kind,9.66265000e-01_fp_kind,6.04662200e+00_fp_kind,7.65980710e+01_fp_kind,2.80961280e+01_fp_kind,&
               4.04655600e+00_fp_kind,8.40000000e+01_fp_kind,2.04864086e+00_fp_kind,&
  85.00_fp_kind,1.60114610e+01_fp_kind,3.26155490e+01_fp_kind,8.11389900e+00_fp_kind,2.88408200e+00_fp_kind,2.13778670e+01_fp_kind,&
               9.26390000e-02_fp_kind,9.04416000e-01_fp_kind,2.65432560e+01_fp_kind,6.83729610e+01_fp_kind,5.49951200e+00_fp_kind,&
               3.99568400e+00_fp_kind,8.50000000e+01_fp_kind,2.04357711e+00_fp_kind,&
  86.00_fp_kind,1.60702280e+01_fp_kind,3.26411050e+01_fp_kind,2.14896590e+01_fp_kind,2.29921800e+00_fp_kind,9.48018400e+00_fp_kind,&
               9.04370000e-02_fp_kind,8.76409000e-01_fp_kind,5.23968700e+00_fp_kind,6.91884770e+01_fp_kind,2.76326400e+01_fp_kind,&
               4.02097700e+00_fp_kind,8.60000000e+01_fp_kind,2.05334332e+00_fp_kind,&
  87.00_fp_kind,1.60073860e+01_fp_kind,3.26638300e+01_fp_kind,2.15943510e+01_fp_kind,1.59849700e+00_fp_kind,1.11211920e+01_fp_kind,&
               8.70310000e-02_fp_kind,8.40187000e-01_fp_kind,4.95446700e+00_fp_kind,1.99805805e+02_fp_kind,2.69051060e+01_fp_kind,&
               4.00347200e+00_fp_kind,8.70000000e+01_fp_kind,2.05640286e+00_fp_kind,&
  88.00_fp_kind,3.25636910e+01_fp_kind,2.13966710e+01_fp_kind,1.12980930e+01_fp_kind,2.83468800e+00_fp_kind,1.59149650e+01_fp_kind,&
               8.01980000e-01_fp_kind,4.59066600e+00_fp_kind,2.27589730e+01_fp_kind,1.60404392e+02_fp_kind,8.35440000e-02_fp_kind,&
               3.98177300e+00_fp_kind,8.80000000e+01_fp_kind,2.06235047e+00_fp_kind,&
  89.00_fp_kind,1.59140530e+01_fp_kind,3.25350420e+01_fp_kind,2.15539760e+01_fp_kind,1.14333940e+01_fp_kind,3.61240900e+00_fp_kind,&
               8.05110000e-02_fp_kind,7.70669000e-01_fp_kind,4.35220600e+00_fp_kind,2.13816220e+01_fp_kind,1.30500748e+02_fp_kind,&
               3.93921200e+00_fp_kind,8.90000000e+01_fp_kind,2.05960633e+00_fp_kind,&
  90.00_fp_kind,1.57840240e+01_fp_kind,3.24548980e+01_fp_kind,2.18492220e+01_fp_kind,4.23907700e+00_fp_kind,1.17361910e+01_fp_kind,&
               7.70670000e-02_fp_kind,7.35137000e-01_fp_kind,4.09797600e+00_fp_kind,1.09464113e+02_fp_kind,2.05121380e+01_fp_kind,&
               3.92253300e+00_fp_kind,9.00000000e+01_fp_kind,2.07609371e+00_fp_kind,&
  91.00_fp_kind,3.27402080e+01_fp_kind,2.19736740e+01_fp_kind,1.29573980e+01_fp_kind,3.68383200e+00_fp_kind,1.57440580e+01_fp_kind,&
               7.09545000e-01_fp_kind,4.05088100e+00_fp_kind,1.92315420e+01_fp_kind,1.17255006e+02_fp_kind,7.40400000e-02_fp_kind,&
               3.88606600e+00_fp_kind,9.10000000e+01_fp_kind,2.08476870e+00_fp_kind,&
  92.00_fp_kind,1.56792750e+01_fp_kind,3.28243050e+01_fp_kind,1.36604590e+01_fp_kind,3.68726100e+00_fp_kind,2.22794350e+01_fp_kind,&
               7.12060000e-02_fp_kind,6.81177000e-01_fp_kind,1.82361570e+01_fp_kind,1.12500040e+02_fp_kind,3.93032500e+00_fp_kind,&
               3.85444400e+00_fp_kind,9.20000000e+01_fp_kind,2.09628663e+00_fp_kind],shape(xrayFF))
  !Ionic form factors
!! Ionic form factor objects
real(fp_kind), parameter :: ionicKnot(27) = &
  real([0.00e0_fp_kind, 0.05e0_fp_kind, 0.10e0_fp_kind, 0.15e0_fp_kind, 0.20e0_fp_kind, 0.25e0_fp_kind, 0.30e0_fp_kind,&
        0.35e0_fp_kind, 0.40e0_fp_kind, 0.45e0_fp_kind, 0.50e0_fp_kind, 0.60e0_fp_kind, 0.70e0_fp_kind, 0.80e0_fp_kind,&
        0.90e0_fp_kind, 1.00e0_fp_kind, 1.20e0_fp_kind, 1.40e0_fp_kind, 1.60e0_fp_kind, 1.80e0_fp_kind, 2.00e0_fp_kind,&
        2.50e0_fp_kind, 3.00e0_fp_kind, 3.50e0_fp_kind, 4.00e0_fp_kind, 5.00e0_fp_kind, 6.00e0_fp_kind],kind=fp_kind)

type :: t_ionicFF
  character(2) :: element
  integer      :: Z,dZ
  real(fp_kind)    :: fx(27)
  real(fp_kind)    :: fx__    !extrapolation parameter
  real(fp_kind)    :: fe0
end type t_ionicFF

 !! Ionic elastic form factors
type(t_ionicFF) :: ionicXrayFF(2,92)
!no data for H
data ionicXrayFF( 1, 1 ) / t_ionicFF('H ',  1 ,  0, &
     [1.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind], 0.0000000000000000e+00_fp_kind,0.000000e+00_fp_kind) /
data ionicXrayFF( 2, 1 ) / t_ionicFF('H ',  1 ,  0, &
     [1.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind], 0.0000000000000000e+00_fp_kind,0.000000e+00_fp_kind) /
data ionicXrayFF( 1, 2 ) / t_ionicFF('He',  2 ,  0, &
     [2.000000e+00_fp_kind,1.957200e+00_fp_kind,1.837300e+00_fp_kind,1.662800e+00_fp_kind,1.460700e+00_fp_kind,&
     1.254500e+00_fp_kind,1.060600e+00_fp_kind,8.877000e-01_fp_kind,7.387000e-01_fp_kind,6.134000e-01_fp_kind,&
     5.093000e-01_fp_kind,3.532000e-01_fp_kind,2.483000e-01_fp_kind,1.774000e-01_fp_kind,1.290000e-01_fp_kind,&
     9.540000e-02_fp_kind,5.460000e-02_fp_kind,3.310000e-02_fp_kind,2.100000e-02_fp_kind,1.390000e-02_fp_kind,&
     9.500000e-03_fp_kind,4.200000e-03_fp_kind,2.100000e-03_fp_kind,1.200000e-03_fp_kind,7.000000e-04_fp_kind,&
     3.000000e-04_fp_kind,1.000000e-04_fp_kind], 1.8000900045009871e-03_fp_kind,4.173000e-01_fp_kind) /
data ionicXrayFF( 2, 2 ) / t_ionicFF('He',  2 ,  0, &
     [2.000000e+00_fp_kind,1.957200e+00_fp_kind,1.837300e+00_fp_kind,1.662800e+00_fp_kind,1.460700e+00_fp_kind,&
     1.254500e+00_fp_kind,1.060600e+00_fp_kind,8.877000e-01_fp_kind,7.387000e-01_fp_kind,6.134000e-01_fp_kind,&
     5.093000e-01_fp_kind,3.532000e-01_fp_kind,2.483000e-01_fp_kind,1.774000e-01_fp_kind,1.290000e-01_fp_kind,&
     9.540000e-02_fp_kind,5.460000e-02_fp_kind,3.310000e-02_fp_kind,2.100000e-02_fp_kind,1.390000e-02_fp_kind,&
     9.500000e-03_fp_kind,4.200000e-03_fp_kind,2.100000e-03_fp_kind,1.200000e-03_fp_kind,7.000000e-04_fp_kind,&
     3.000000e-04_fp_kind,1.000000e-04_fp_kind], 1.8000900045009871e-03_fp_kind,4.173000e-01_fp_kind) /
data ionicXrayFF( 1, 3 ) / t_ionicFF('Li',  3 ,  0, &
     [3.000000e+00_fp_kind,2.707700e+00_fp_kind,2.215300e+00_fp_kind,1.903800e+00_fp_kind,1.741700e+00_fp_kind,&
     1.626200e+00_fp_kind,1.512900e+00_fp_kind,1.393400e+00_fp_kind,1.270700e+00_fp_kind,1.149300e+00_fp_kind,&
     1.033000e+00_fp_kind,8.237000e-01_fp_kind,6.508000e-01_fp_kind,5.127000e-01_fp_kind,4.045000e-01_fp_kind,&
     3.205000e-01_fp_kind,2.048000e-01_fp_kind,1.345000e-01_fp_kind,9.080000e-02_fp_kind,6.300000e-02_fp_kind,&
     4.480000e-02_fp_kind,2.090000e-02_fp_kind,1.090000e-02_fp_kind,6.200000e-03_fp_kind,3.700000e-03_fp_kind,&
     1.600000e-03_fp_kind,8.000000e-04_fp_kind], 9.6025606828433752e-03_fp_kind,3.255600e+00_fp_kind) /
data ionicXrayFF( 2, 3 ) / t_ionicFF('Li',  3 ,  1, &
     [2.000000e+00_fp_kind,1.983700e+00_fp_kind,1.936100e+00_fp_kind,1.860700e+00_fp_kind,1.762700e+00_fp_kind,&
     1.648200e+00_fp_kind,1.523600e+00_fp_kind,1.394600e+00_fp_kind,1.266100e+00_fp_kind,1.141800e+00_fp_kind,&
     1.024200e+00_fp_kind,8.151000e-01_fp_kind,6.435000e-01_fp_kind,5.069000e-01_fp_kind,3.999000e-01_fp_kind,&
     3.169000e-01_fp_kind,2.025000e-01_fp_kind,1.330000e-01_fp_kind,8.980000e-02_fp_kind,6.230000e-02_fp_kind,&
     4.430000e-02_fp_kind,2.070000e-02_fp_kind,1.080000e-02_fp_kind,6.100000e-03_fp_kind,3.700000e-03_fp_kind,&
     1.600000e-03_fp_kind,8.000000e-04_fp_kind], 9.6025606828433752e-03_fp_kind,1.569000e-01_fp_kind) /
data ionicXrayFF( 1, 4 ) / t_ionicFF('Be',  4 ,  0, &
     [4.000000e+00_fp_kind,3.706900e+00_fp_kind,3.065700e+00_fp_kind,2.463100e+00_fp_kind,2.060200e+00_fp_kind,&
     1.827900e+00_fp_kind,1.692500e+00_fp_kind,1.599800e+00_fp_kind,1.520700e+00_fp_kind,1.442900e+00_fp_kind,&
     1.362600e+00_fp_kind,1.195700e+00_fp_kind,1.030300e+00_fp_kind,8.768000e-01_fp_kind,7.403000e-01_fp_kind,&
     6.225000e-01_fp_kind,4.389000e-01_fp_kind,3.112000e-01_fp_kind,2.233000e-01_fp_kind,1.625000e-01_fp_kind,&
     1.202000e-01_fp_kind,6.030000e-02_fp_kind,3.280000e-02_fp_kind,1.920000e-02_fp_kind,1.190000e-02_fp_kind,&
     5.200000e-03_fp_kind,2.600000e-03_fp_kind], 2.3415219892932271e-02_fp_kind,3.038300e+00_fp_kind) /
data ionicXrayFF( 2, 4 ) / t_ionicFF('Be',  4 ,  2, &
     [2.000000e+00_fp_kind,1.991500e+00_fp_kind,1.966300e+00_fp_kind,1.925400e+00_fp_kind,1.870400e+00_fp_kind,&
     1.803200e+00_fp_kind,1.726100e+00_fp_kind,1.641300e+00_fp_kind,1.551400e+00_fp_kind,1.458400e+00_fp_kind,&
     1.364400e+00_fp_kind,1.179800e+00_fp_kind,1.007500e+00_fp_kind,8.530000e-01_fp_kind,7.182000e-01_fp_kind,&
     6.031000e-01_fp_kind,4.248000e-01_fp_kind,3.012000e-01_fp_kind,2.162000e-01_fp_kind,1.574000e-01_fp_kind,&
     1.164000e-01_fp_kind,5.840000e-02_fp_kind,3.180000e-02_fp_kind,1.860000e-02_fp_kind,1.150000e-02_fp_kind,&
     5.000000e-03_fp_kind,2.500000e-03_fp_kind], 2.2514071294558850e-02_fp_kind,8.170000e-02_fp_kind) /
data ionicXrayFF( 1, 5 ) / t_ionicFF('B ',  5 ,  0, &
     [5.000000e+00_fp_kind,4.724600e+00_fp_kind,4.060300e+00_fp_kind,3.316600e+00_fp_kind,2.699400e+00_fp_kind,&
     2.263200e+00_fp_kind,1.979300e+00_fp_kind,1.799200e+00_fp_kind,1.681000e+00_fp_kind,1.596000e+00_fp_kind,&
     1.526700e+00_fp_kind,1.401900e+00_fp_kind,1.276100e+00_fp_kind,1.147500e+00_fp_kind,1.020800e+00_fp_kind,&
     9.007000e-01_fp_kind,6.906000e-01_fp_kind,5.247000e-01_fp_kind,3.984000e-01_fp_kind,3.039000e-01_fp_kind,&
     2.335000e-01_fp_kind,1.258000e-01_fp_kind,7.190000e-02_fp_kind,4.340000e-02_fp_kind,2.750000e-02_fp_kind,&
     1.240000e-02_fp_kind,6.300000e-03_fp_kind], 4.5417225704397879e-02_fp_kind,2.785000e+00_fp_kind) /
data ionicXrayFF( 2, 5 ) / t_ionicFF('B ',  5 ,  0, &
     [5.000000e+00_fp_kind,4.724600e+00_fp_kind,4.060300e+00_fp_kind,3.316600e+00_fp_kind,2.699400e+00_fp_kind,&
     2.263200e+00_fp_kind,1.979300e+00_fp_kind,1.799200e+00_fp_kind,1.681000e+00_fp_kind,1.596000e+00_fp_kind,&
     1.526700e+00_fp_kind,1.401900e+00_fp_kind,1.276100e+00_fp_kind,1.147500e+00_fp_kind,1.020800e+00_fp_kind,&
     9.007000e-01_fp_kind,6.906000e-01_fp_kind,5.247000e-01_fp_kind,3.984000e-01_fp_kind,3.039000e-01_fp_kind,&
     2.335000e-01_fp_kind,1.258000e-01_fp_kind,7.190000e-02_fp_kind,4.340000e-02_fp_kind,2.750000e-02_fp_kind,&
     1.240000e-02_fp_kind,6.300000e-03_fp_kind], 4.5417225704397879e-02_fp_kind,2.785000e+00_fp_kind) /
data ionicXrayFF( 1, 6 ) / t_ionicFF('C ',  6 ,  0, &
     [6.000000e+00_fp_kind,5.751800e+00_fp_kind,5.115300e+00_fp_kind,4.322200e+00_fp_kind,3.569800e+00_fp_kind,&
     2.955800e+00_fp_kind,2.497700e+00_fp_kind,2.173400e+00_fp_kind,1.949500e+00_fp_kind,1.794900e+00_fp_kind,&
     1.685100e+00_fp_kind,1.536800e+00_fp_kind,1.425800e+00_fp_kind,1.322400e+00_fp_kind,1.218700e+00_fp_kind,&
     1.114500e+00_fp_kind,9.141000e-01_fp_kind,7.367000e-01_fp_kind,5.884000e-01_fp_kind,4.687000e-01_fp_kind,&
     3.736000e-01_fp_kind,2.160000e-01_fp_kind,1.296000e-01_fp_kind,8.100000e-02_fp_kind,5.260000e-02_fp_kind,&
     2.450000e-02_fp_kind,1.280000e-02_fp_kind], 7.6964190272583721e-02_fp_kind,2.470900e+00_fp_kind) /
data ionicXrayFF( 2, 6 ) / t_ionicFF('C ',  6 ,  0, &
     [6.000000e+00_fp_kind,5.751800e+00_fp_kind,5.115300e+00_fp_kind,4.322200e+00_fp_kind,3.569800e+00_fp_kind,&
     2.955800e+00_fp_kind,2.497700e+00_fp_kind,2.173400e+00_fp_kind,1.949500e+00_fp_kind,1.794900e+00_fp_kind,&
     1.685100e+00_fp_kind,1.536800e+00_fp_kind,1.425800e+00_fp_kind,1.322400e+00_fp_kind,1.218700e+00_fp_kind,&
     1.114500e+00_fp_kind,9.141000e-01_fp_kind,7.367000e-01_fp_kind,5.884000e-01_fp_kind,4.687000e-01_fp_kind,&
     3.736000e-01_fp_kind,2.160000e-01_fp_kind,1.296000e-01_fp_kind,8.100000e-02_fp_kind,5.260000e-02_fp_kind,&
     2.450000e-02_fp_kind,1.280000e-02_fp_kind], 7.6964190272583721e-02_fp_kind,2.470900e+00_fp_kind) /
data ionicXrayFF( 1, 7 ) / t_ionicFF('N ',  7 ,  0, &
     [7.000000e+00_fp_kind,6.776500e+00_fp_kind,6.181000e+00_fp_kind,5.386600e+00_fp_kind,4.564800e+00_fp_kind,&
     3.826800e+00_fp_kind,3.219700e+00_fp_kind,2.747500e+00_fp_kind,2.392900e+00_fp_kind,2.132200e+00_fp_kind,&
     1.941800e+00_fp_kind,1.697400e+00_fp_kind,1.551000e+00_fp_kind,1.445100e+00_fp_kind,1.353300e+00_fp_kind,&
     1.265000e+00_fp_kind,1.090000e+00_fp_kind,9.218000e-01_fp_kind,7.691000e-01_fp_kind,6.367000e-01_fp_kind,&
     5.251000e-01_fp_kind,3.246000e-01_fp_kind,2.044000e-01_fp_kind,1.324000e-01_fp_kind,8.830000e-02_fp_kind,&
     4.270000e-02_fp_kind,2.270000e-02_fp_kind], 1.1712266922735640e-01_fp_kind,2.203400e+00_fp_kind) /
data ionicXrayFF( 2, 7 ) / t_ionicFF('N ',  7 ,  0, &
     [7.000000e+00_fp_kind,6.776500e+00_fp_kind,6.181000e+00_fp_kind,5.386600e+00_fp_kind,4.564800e+00_fp_kind,&
     3.826800e+00_fp_kind,3.219700e+00_fp_kind,2.747500e+00_fp_kind,2.392900e+00_fp_kind,2.132200e+00_fp_kind,&
     1.941800e+00_fp_kind,1.697400e+00_fp_kind,1.551000e+00_fp_kind,1.445100e+00_fp_kind,1.353300e+00_fp_kind,&
     1.265000e+00_fp_kind,1.090000e+00_fp_kind,9.218000e-01_fp_kind,7.691000e-01_fp_kind,6.367000e-01_fp_kind,&
     5.251000e-01_fp_kind,3.246000e-01_fp_kind,2.044000e-01_fp_kind,1.324000e-01_fp_kind,8.830000e-02_fp_kind,&
     4.270000e-02_fp_kind,2.270000e-02_fp_kind], 1.1712266922735640e-01_fp_kind,2.203400e+00_fp_kind) /
data ionicXrayFF( 1, 8 ) / t_ionicFF('O ',  8 ,  0, &
     [8.000000e+00_fp_kind,7.797400e+00_fp_kind,7.243900e+00_fp_kind,6.470600e+00_fp_kind,5.621500e+00_fp_kind,&
     4.806400e+00_fp_kind,4.087200e+00_fp_kind,3.487300e+00_fp_kind,3.005400e+00_fp_kind,2.627900e+00_fp_kind,&
     2.337100e+00_fp_kind,1.945500e+00_fp_kind,1.714400e+00_fp_kind,1.568200e+00_fp_kind,1.463500e+00_fp_kind,&
     1.377200e+00_fp_kind,1.221400e+00_fp_kind,1.070700e+00_fp_kind,9.261000e-01_fp_kind,7.927000e-01_fp_kind,&
     6.739000e-01_fp_kind,4.434000e-01_fp_kind,2.926000e-01_fp_kind,1.963000e-01_fp_kind,1.345000e-01_fp_kind,&
     6.750000e-02_fp_kind,3.680000e-02_fp_kind], 1.6636528028933381e-01_fp_kind,1.983900e+00_fp_kind) /
data ionicXrayFF( 2, 8 ) / t_ionicFF('O ',  8 , -2, &
     [1.000000e+01_fp_kind,9.588400e+00_fp_kind,8.534500e+00_fp_kind,7.218800e+00_fp_kind,5.955800e+00_fp_kind,&
     4.893400e+00_fp_kind,4.057400e+00_fp_kind,3.419000e+00_fp_kind,2.936400e+00_fp_kind,2.572000e+00_fp_kind,&
     2.296300e+00_fp_kind,1.927500e+00_fp_kind,1.708500e+00_fp_kind,1.567800e+00_fp_kind,1.464400e+00_fp_kind,&
     1.377400e+00_fp_kind,1.219100e+00_fp_kind,1.066800e+00_fp_kind,9.215000e-01_fp_kind,7.881000e-01_fp_kind,&
     6.696000e-01_fp_kind,4.404000e-01_fp_kind,2.907000e-01_fp_kind,1.951000e-01_fp_kind,1.332000e-01_fp_kind,&
     6.860000e-02_fp_kind,3.510000e-02_fp_kind], 1.5864605958643099e-01_fp_kind,4.099200e+00_fp_kind) /
data ionicXrayFF( 1, 9 ) / t_ionicFF('F ',  9 ,  0, &
     [9.000000e+00_fp_kind,8.815200e+00_fp_kind,8.301100e+00_fp_kind,7.558800e+00_fp_kind,6.708000e+00_fp_kind,&
     5.849900e+00_fp_kind,5.052700e+00_fp_kind,4.351500e+00_fp_kind,3.758100e+00_fp_kind,3.269100e+00_fp_kind,&
     2.873800e+00_fp_kind,2.308600e+00_fp_kind,1.956500e+00_fp_kind,1.734700e+00_fp_kind,1.588100e+00_fp_kind,&
     1.482500e+00_fp_kind,1.323900e+00_fp_kind,1.186800e+00_fp_kind,1.054900e+00_fp_kind,9.284000e-01_fp_kind,&
     8.103000e-01_fp_kind,5.643000e-01_fp_kind,3.894000e-01_fp_kind,2.703000e-01_fp_kind,1.902000e-01_fp_kind,&
     9.930000e-02_fp_kind,5.520000e-02_fp_kind], 2.2216259726320689e-01_fp_kind,1.801700e+00_fp_kind) /
data ionicXrayFF( 2, 9 ) / t_ionicFF('F ',  9 , -1, &
     [1.000000e+01_fp_kind,9.733800e+00_fp_kind,9.015700e+00_fp_kind,8.032400e+00_fp_kind,6.974200e+00_fp_kind,&
     5.971200e+00_fp_kind,5.087800e+00_fp_kind,4.343000e+00_fp_kind,3.731800e+00_fp_kind,3.238900e+00_fp_kind,&
     2.846100e+00_fp_kind,2.291000e+00_fp_kind,1.947500e+00_fp_kind,1.731200e+00_fp_kind,1.587200e+00_fp_kind,&
     1.482500e+00_fp_kind,1.323600e+00_fp_kind,1.185600e+00_fp_kind,1.053000e+00_fp_kind,9.262000e-01_fp_kind,&
     8.080000e-01_fp_kind,5.625000e-01_fp_kind,3.881000e-01_fp_kind,2.695000e-01_fp_kind,1.895000e-01_fp_kind,&
     9.940000e-02_fp_kind,5.470000e-02_fp_kind], 2.2013794953775090e-01_fp_kind,2.615900e+00_fp_kind) /
data ionicXrayFF( 1, 10) / t_ionicFF('Ne',  10,  0, &
     [1.000000e+01_fp_kind,9.830300e+00_fp_kind,9.351900e+00_fp_kind,8.644200e+00_fp_kind,7.806200e+00_fp_kind,&
     6.929100e+00_fp_kind,6.080900e+00_fp_kind,5.303900e+00_fp_kind,4.618700e+00_fp_kind,4.030800e+00_fp_kind,&
     3.536500e+00_fp_kind,2.791200e+00_fp_kind,2.296400e+00_fp_kind,1.972000e+00_fp_kind,1.757200e+00_fp_kind,&
     1.609800e+00_fp_kind,1.418000e+00_fp_kind,1.280800e+00_fp_kind,1.158600e+00_fp_kind,1.041600e+00_fp_kind,&
     9.292000e-01_fp_kind,6.808000e-01_fp_kind,4.896000e-01_fp_kind,3.514000e-01_fp_kind,2.540000e-01_fp_kind,&
     1.377000e-01_fp_kind,7.870000e-02_fp_kind], 2.8556741556045750e-01_fp_kind,1.649400e+00_fp_kind) /
data ionicXrayFF( 2, 10) / t_ionicFF('Ne',  10,  0, &
     [1.000000e+01_fp_kind,9.830300e+00_fp_kind,9.351900e+00_fp_kind,8.644200e+00_fp_kind,7.806200e+00_fp_kind,&
     6.929100e+00_fp_kind,6.080900e+00_fp_kind,5.303900e+00_fp_kind,4.618700e+00_fp_kind,4.030800e+00_fp_kind,&
     3.536500e+00_fp_kind,2.791200e+00_fp_kind,2.296400e+00_fp_kind,1.972000e+00_fp_kind,1.757200e+00_fp_kind,&
     1.609800e+00_fp_kind,1.418000e+00_fp_kind,1.280800e+00_fp_kind,1.158600e+00_fp_kind,1.041600e+00_fp_kind,&
     9.292000e-01_fp_kind,6.808000e-01_fp_kind,4.896000e-01_fp_kind,3.514000e-01_fp_kind,2.540000e-01_fp_kind,&
     1.377000e-01_fp_kind,7.870000e-02_fp_kind], 2.8556741556045750e-01_fp_kind,1.649400e+00_fp_kind) /
data ionicXrayFF( 1, 11) / t_ionicFF('Na',  11,  0, &
     [1.100000e+01_fp_kind,1.056800e+01_fp_kind,9.760700e+00_fp_kind,9.027400e+00_fp_kind,8.336400e+00_fp_kind,&
     7.619600e+00_fp_kind,6.882600e+00_fp_kind,6.157200e+00_fp_kind,5.473000e+00_fp_kind,4.849100e+00_fp_kind,&
     4.294900e+00_fp_kind,3.398900e+00_fp_kind,2.754700e+00_fp_kind,2.305900e+00_fp_kind,1.997500e+00_fp_kind,&
     1.784600e+00_fp_kind,1.524100e+00_fp_kind,1.367100e+00_fp_kind,1.246500e+00_fp_kind,1.137500e+00_fp_kind,&
     1.033000e+00_fp_kind,7.917000e-01_fp_kind,5.920000e-01_fp_kind,4.393000e-01_fp_kind,3.234000e-01_fp_kind,&
     1.844000e-01_fp_kind,1.067000e-01_fp_kind], 3.5262041805513888e-01_fp_kind,4.738700e+00_fp_kind) /
data ionicXrayFF( 2, 11) / t_ionicFF('Na',  11,  1, &
     [1.000000e+01_fp_kind,9.883200e+00_fp_kind,9.546300e+00_fp_kind,9.026200e+00_fp_kind,8.374600e+00_fp_kind,&
     7.647500e+00_fp_kind,6.896100e+00_fp_kind,6.161400e+00_fp_kind,5.472400e+00_fp_kind,4.846400e+00_fp_kind,&
     4.291700e+00_fp_kind,3.396500e+00_fp_kind,2.753500e+00_fp_kind,2.305800e+00_fp_kind,1.997800e+00_fp_kind,&
     1.785200e+00_fp_kind,1.524600e+00_fp_kind,1.367300e+00_fp_kind,1.246400e+00_fp_kind,1.137100e+00_fp_kind,&
     1.032500e+00_fp_kind,7.909000e-01_fp_kind,5.912000e-01_fp_kind,4.380000e-01_fp_kind,3.247000e-01_fp_kind,&
     1.827000e-01_fp_kind,1.074000e-01_fp_kind], 3.5495657602409381e-01_fp_kind,1.128500e+00_fp_kind) /
data ionicXrayFF( 1, 12) / t_ionicFF('Mg',  12,  0, &
     [1.200000e+01_fp_kind,1.150760e+01_fp_kind,1.047310e+01_fp_kind,9.502700e+00_fp_kind,8.736400e+00_fp_kind,&
     8.078700e+00_fp_kind,7.447900e+00_fp_kind,6.818300e+00_fp_kind,6.197200e+00_fp_kind,5.597000e+00_fp_kind,&
     5.035400e+00_fp_kind,4.067300e+00_fp_kind,3.298100e+00_fp_kind,2.729600e+00_fp_kind,2.317400e+00_fp_kind,&
     2.022800e+00_fp_kind,1.660100e+00_fp_kind,1.459100e+00_fp_kind,1.326100e+00_fp_kind,1.219000e+00_fp_kind,&
     1.121000e+00_fp_kind,8.921000e-01_fp_kind,6.912000e-01_fp_kind,5.284000e-01_fp_kind,3.990000e-01_fp_kind,&
     2.372000e-01_fp_kind,1.394000e-01_fp_kind], 4.2311518810178228e-01_fp_kind,5.178100e+00_fp_kind) /
data ionicXrayFF( 2, 12) / t_ionicFF('Mg',  12,  2, &
     [1.000000e+01_fp_kind,9.913900e+00_fp_kind,9.662500e+00_fp_kind,9.265500e+00_fp_kind,8.752300e+00_fp_kind,&
     8.157100e+00_fp_kind,7.514900e+00_fp_kind,6.857600e+00_fp_kind,6.211500e+00_fp_kind,5.596800e+00_fp_kind,&
     5.026800e+00_fp_kind,4.048000e+00_fp_kind,3.289000e+00_fp_kind,2.724800e+00_fp_kind,2.316000e+00_fp_kind,&
     2.023400e+00_fp_kind,1.662100e+00_fp_kind,1.460500e+00_fp_kind,1.326500e+00_fp_kind,1.218500e+00_fp_kind,&
     1.119800e+00_fp_kind,8.902000e-01_fp_kind,6.892000e-01_fp_kind,5.260000e-01_fp_kind,3.998000e-01_fp_kind,&
     2.336000e-01_fp_kind,1.409000e-01_fp_kind], 4.2772217115970790e-01_fp_kind,8.300000e-01_fp_kind) /
data ionicXrayFF( 1, 13) / t_ionicFF('Al',  13,  0, &
     [1.300000e+01_fp_kind,1.243860e+01_fp_kind,1.122900e+01_fp_kind,1.005920e+01_fp_kind,9.158700e+00_fp_kind,&
     8.465800e+00_fp_kind,7.874300e+00_fp_kind,7.317200e+00_fp_kind,6.768000e+00_fp_kind,6.224300e+00_fp_kind,&
     5.694300e+00_fp_kind,4.714200e+00_fp_kind,3.884900e+00_fp_kind,3.222000e+00_fp_kind,2.712700e+00_fp_kind,&
     2.331100e+00_fp_kind,1.841700e+00_fp_kind,1.571700e+00_fp_kind,1.408500e+00_fp_kind,1.292600e+00_fp_kind,&
     1.196100e+00_fp_kind,9.799000e-01_fp_kind,7.836000e-01_fp_kind,6.160000e-01_fp_kind,4.769000e-01_fp_kind,&
     2.936000e-01_fp_kind,1.772000e-01_fp_kind], 4.9748884798951559e-01_fp_kind,5.866500e+00_fp_kind) /
data ionicXrayFF( 2, 13) / t_ionicFF('Al',  13,  3, &
     [1.000000e+01_fp_kind,9.933600e+00_fp_kind,9.738300e+00_fp_kind,9.425700e+00_fp_kind,9.013700e+00_fp_kind,&
     8.523900e+00_fp_kind,7.980000e+00_fp_kind,7.405300e+00_fp_kind,6.820900e+00_fp_kind,6.245000e+00_fp_kind,&
     5.691600e+00_fp_kind,4.690500e+00_fp_kind,3.860400e+00_fp_kind,3.204000e+00_fp_kind,2.701800e+00_fp_kind,&
     2.326000e+00_fp_kind,1.842700e+00_fp_kind,1.574000e+00_fp_kind,1.410200e+00_fp_kind,1.293100e+00_fp_kind,&
     1.195500e+00_fp_kind,9.776000e-01_fp_kind,7.808000e-01_fp_kind,6.128000e-01_fp_kind,4.770000e-01_fp_kind,&
     2.890000e-01_fp_kind,1.789000e-01_fp_kind], 5.0232819336874712e-01_fp_kind,6.392000e-01_fp_kind) /
data ionicXrayFF( 1, 14) / t_ionicFF('Si',  14,  0, &
     [1.400000e+01_fp_kind,1.343810e+01_fp_kind,1.214150e+01_fp_kind,1.077410e+01_fp_kind,9.674100e+00_fp_kind,&
     8.858600e+00_fp_kind,8.231000e+00_fp_kind,7.698700e+00_fp_kind,7.204100e+00_fp_kind,6.720800e+00_fp_kind,&
     6.242200e+00_fp_kind,5.314600e+00_fp_kind,4.472000e+00_fp_kind,3.752000e+00_fp_kind,3.165300e+00_fp_kind,&
     2.703000e+00_fp_kind,2.076500e+00_fp_kind,1.717100e+00_fp_kind,1.505300e+00_fp_kind,1.367700e+00_fp_kind,&
     1.264800e+00_fp_kind,1.055900e+00_fp_kind,8.677000e-01_fp_kind,6.997000e-01_fp_kind,5.551000e-01_fp_kind,&
     3.529000e-01_fp_kind,2.195000e-01_fp_kind], 5.7341896157614514e-01_fp_kind,5.753500e+00_fp_kind) /
data ionicXrayFF( 2, 14) / t_ionicFF('Si',  14,  0, &
     [1.400000e+01_fp_kind,1.343810e+01_fp_kind,1.214150e+01_fp_kind,1.077410e+01_fp_kind,9.674100e+00_fp_kind,&
     8.858600e+00_fp_kind,8.231000e+00_fp_kind,7.698700e+00_fp_kind,7.204100e+00_fp_kind,6.720800e+00_fp_kind,&
     6.242200e+00_fp_kind,5.314600e+00_fp_kind,4.472000e+00_fp_kind,3.752000e+00_fp_kind,3.165300e+00_fp_kind,&
     2.703000e+00_fp_kind,2.076500e+00_fp_kind,1.717100e+00_fp_kind,1.505300e+00_fp_kind,1.367700e+00_fp_kind,&
     1.264800e+00_fp_kind,1.055900e+00_fp_kind,8.677000e-01_fp_kind,6.997000e-01_fp_kind,5.551000e-01_fp_kind,&
     3.529000e-01_fp_kind,2.195000e-01_fp_kind], 5.7341896157614514e-01_fp_kind,5.753500e+00_fp_kind) /
data ionicXrayFF( 1, 15) / t_ionicFF('P ',  15,  0, &
     [1.500000e+01_fp_kind,1.445750e+01_fp_kind,1.313690e+01_fp_kind,1.162780e+01_fp_kind,1.032540e+01_fp_kind,&
     9.334100e+00_fp_kind,8.599800e+00_fp_kind,8.030100e+00_fp_kind,7.548300e+00_fp_kind,7.105500e+00_fp_kind,&
     6.676500e+00_fp_kind,5.831700e+00_fp_kind,5.021500e+00_fp_kind,4.285400e+00_fp_kind,3.650100e+00_fp_kind,&
     3.123200e+00_fp_kind,2.364700e+00_fp_kind,1.903300e+00_fp_kind,1.626400e+00_fp_kind,1.453100e+00_fp_kind,&
     1.333900e+00_fp_kind,1.122500e+00_fp_kind,9.428000e-01_fp_kind,7.780000e-01_fp_kind,6.315000e-01_fp_kind,&
     4.139000e-01_fp_kind,2.655000e-01_fp_kind], 6.4868166547897488e-01_fp_kind,5.473800e+00_fp_kind) /
data ionicXrayFF( 2, 15) / t_ionicFF('P ',  15,  0, &
     [1.500000e+01_fp_kind,1.445750e+01_fp_kind,1.313690e+01_fp_kind,1.162780e+01_fp_kind,1.032540e+01_fp_kind,&
     9.334100e+00_fp_kind,8.599800e+00_fp_kind,8.030100e+00_fp_kind,7.548300e+00_fp_kind,7.105500e+00_fp_kind,&
     6.676500e+00_fp_kind,5.831700e+00_fp_kind,5.021500e+00_fp_kind,4.285400e+00_fp_kind,3.650100e+00_fp_kind,&
     3.123200e+00_fp_kind,2.364700e+00_fp_kind,1.903300e+00_fp_kind,1.626400e+00_fp_kind,1.453100e+00_fp_kind,&
     1.333900e+00_fp_kind,1.122500e+00_fp_kind,9.428000e-01_fp_kind,7.780000e-01_fp_kind,6.315000e-01_fp_kind,&
     4.139000e-01_fp_kind,2.655000e-01_fp_kind], 6.4868166547897488e-01_fp_kind,5.473800e+00_fp_kind) /
data ionicXrayFF( 1, 16) / t_ionicFF('S ',  16,  0, &
     [1.600000e+01_fp_kind,1.548280e+01_fp_kind,1.417350e+01_fp_kind,1.257810e+01_fp_kind,1.110410e+01_fp_kind,&
     9.923900e+00_fp_kind,9.038200e+00_fp_kind,8.376000e+00_fp_kind,7.857300e+00_fp_kind,7.419000e+00_fp_kind,&
     7.019400e+00_fp_kind,6.256000e+00_fp_kind,5.507100e+00_fp_kind,4.791400e+00_fp_kind,4.139500e+00_fp_kind,&
     3.571200e+00_fp_kind,2.700100e+00_fp_kind,2.133200e+00_fp_kind,1.779100e+00_fp_kind,1.557200e+00_fp_kind,&
     1.411000e+00_fp_kind,1.182700e+00_fp_kind,1.009400e+00_fp_kind,8.501000e-01_fp_kind,7.045000e-01_fp_kind,&
     4.759000e-01_fp_kind,3.146000e-01_fp_kind], 7.2204725413442361e-01_fp_kind,5.164000e+00_fp_kind) /
data ionicXrayFF( 2, 16) / t_ionicFF('S ',  16,  0, &
     [1.600000e+01_fp_kind,1.548280e+01_fp_kind,1.417350e+01_fp_kind,1.257810e+01_fp_kind,1.110410e+01_fp_kind,&
     9.923900e+00_fp_kind,9.038200e+00_fp_kind,8.376000e+00_fp_kind,7.857300e+00_fp_kind,7.419000e+00_fp_kind,&
     7.019400e+00_fp_kind,6.256000e+00_fp_kind,5.507100e+00_fp_kind,4.791400e+00_fp_kind,4.139500e+00_fp_kind,&
     3.571200e+00_fp_kind,2.700100e+00_fp_kind,2.133200e+00_fp_kind,1.779100e+00_fp_kind,1.557200e+00_fp_kind,&
     1.411000e+00_fp_kind,1.182700e+00_fp_kind,1.009400e+00_fp_kind,8.501000e-01_fp_kind,7.045000e-01_fp_kind,&
     4.759000e-01_fp_kind,3.146000e-01_fp_kind], 7.2204725413442361e-01_fp_kind,5.164000e+00_fp_kind) /
data ionicXrayFF( 1, 17) / t_ionicFF('Cl',  17,  0, &
     [1.700000e+01_fp_kind,1.650990e+01_fp_kind,1.523190e+01_fp_kind,1.359390e+01_fp_kind,1.198870e+01_fp_kind,&
     1.063180e+01_fp_kind,9.575800e+00_fp_kind,8.782400e+00_fp_kind,8.182100e+00_fp_kind,7.707600e+00_fp_kind,&
     7.306500e+00_fp_kind,6.596800e+00_fp_kind,5.916900e+00_fp_kind,5.247300e+00_fp_kind,4.608500e+00_fp_kind,&
     4.024300e+00_fp_kind,3.071200e+00_fp_kind,2.405400e+00_fp_kind,1.967800e+00_fp_kind,1.686400e+00_fp_kind,&
     1.502600e+00_fp_kind,1.240800e+00_fp_kind,1.068900e+00_fp_kind,9.156000e-01_fp_kind,7.730000e-01_fp_kind,&
     5.378000e-01_fp_kind,3.656000e-01_fp_kind], 7.9122781705382295e-01_fp_kind,4.856600e+00_fp_kind) /
data ionicXrayFF( 2, 17) / t_ionicFF('Cl',  17, -1, &
     [1.800000e+01_fp_kind,1.736130e+01_fp_kind,1.571400e+01_fp_kind,1.380270e+01_fp_kind,1.201250e+01_fp_kind,&
     1.058500e+01_fp_kind,9.521500e+00_fp_kind,8.742800e+00_fp_kind,8.159500e+00_fp_kind,7.697800e+00_fp_kind,&
     7.304900e+00_fp_kind,6.601900e+00_fp_kind,5.922200e+00_fp_kind,5.251100e+00_fp_kind,4.610900e+00_fp_kind,&
     4.025800e+00_fp_kind,3.072100e+00_fp_kind,2.406200e+00_fp_kind,1.968500e+00_fp_kind,1.687000e+00_fp_kind,&
     1.503000e+00_fp_kind,1.240800e+00_fp_kind,1.068700e+00_fp_kind,9.154000e-01_fp_kind,7.721000e-01_fp_kind,&
     5.387000e-01_fp_kind,3.645000e-01_fp_kind], 7.8879504673739831e-01_fp_kind,6.385800e+00_fp_kind) /
data ionicXrayFF( 1, 18) / t_ionicFF('Ar',  18,  0, &
     [1.800000e+01_fp_kind,1.753610e+01_fp_kind,1.629880e+01_fp_kind,1.464890e+01_fp_kind,1.295100e+01_fp_kind,&
     1.144260e+01_fp_kind,1.021790e+01_fp_kind,9.273600e+00_fp_kind,8.559200e+00_fp_kind,8.012400e+00_fp_kind,&
     7.576600e+00_fp_kind,6.876200e+00_fp_kind,6.253300e+00_fp_kind,5.640900e+00_fp_kind,5.037300e+00_fp_kind,&
     4.461500e+00_fp_kind,3.463100e+00_fp_kind,2.713800e+00_fp_kind,2.192900e+00_fp_kind,1.845000e+00_fp_kind,&
     1.614500e+00_fp_kind,1.301100e+00_fp_kind,1.123000e+00_fp_kind,9.747000e-01_fp_kind,8.364000e-01_fp_kind,&
     5.987000e-01_fp_kind,4.179000e-01_fp_kind], 8.5566570546180287e-01_fp_kind,4.571200e+00_fp_kind) /
data ionicXrayFF( 2, 18) / t_ionicFF('Ar',  18,  0, &
     [1.800000e+01_fp_kind,1.753610e+01_fp_kind,1.629880e+01_fp_kind,1.464890e+01_fp_kind,1.295100e+01_fp_kind,&
     1.144260e+01_fp_kind,1.021790e+01_fp_kind,9.273600e+00_fp_kind,8.559200e+00_fp_kind,8.012400e+00_fp_kind,&
     7.576600e+00_fp_kind,6.876200e+00_fp_kind,6.253300e+00_fp_kind,5.640900e+00_fp_kind,5.037300e+00_fp_kind,&
     4.461500e+00_fp_kind,3.463100e+00_fp_kind,2.713800e+00_fp_kind,2.192900e+00_fp_kind,1.845000e+00_fp_kind,&
     1.614500e+00_fp_kind,1.301100e+00_fp_kind,1.123000e+00_fp_kind,9.747000e-01_fp_kind,8.364000e-01_fp_kind,&
     5.987000e-01_fp_kind,4.179000e-01_fp_kind], 8.5566570546180287e-01_fp_kind,4.571200e+00_fp_kind) /
data ionicXrayFF( 1, 19) / t_ionicFF('K ',  19,  0, &
     [1.900000e+01_fp_kind,1.820460e+01_fp_kind,1.673410e+01_fp_kind,1.524440e+01_fp_kind,1.373010e+01_fp_kind,&
     1.227060e+01_fp_kind,1.097860e+01_fp_kind,9.909900e+00_fp_kind,9.063100e+00_fp_kind,8.404700e+00_fp_kind,&
     7.890000e+00_fp_kind,7.126400e+00_fp_kind,6.524300e+00_fp_kind,5.963200e+00_fp_kind,5.407300e+00_fp_kind,&
     4.860300e+00_fp_kind,3.856200e+00_fp_kind,3.046800e+00_fp_kind,2.450000e+00_fp_kind,2.033000e+00_fp_kind,&
     1.749400e+00_fp_kind,1.368200e+00_fp_kind,1.175100e+00_fp_kind,1.030100e+00_fp_kind,8.933000e-01_fp_kind,&
     6.590000e-01_fp_kind,4.709000e-01_fp_kind], 9.1490682224176823e-01_fp_kind,8.896600e+00_fp_kind) /
data ionicXrayFF( 2, 19) / t_ionicFF('K ',  19,  1, &
     [1.800000e+01_fp_kind,1.764870e+01_fp_kind,1.667800e+01_fp_kind,1.529870e+01_fp_kind,1.376180e+01_fp_kind,&
     1.227730e+01_fp_kind,1.097440e+01_fp_kind,9.903100e+00_fp_kind,9.057200e+00_fp_kind,8.400700e+00_fp_kind,&
     7.887900e+00_fp_kind,7.126600e+00_fp_kind,6.525100e+00_fp_kind,5.963900e+00_fp_kind,5.407700e+00_fp_kind,&
     4.860300e+00_fp_kind,3.856000e+00_fp_kind,3.046700e+00_fp_kind,2.450000e+00_fp_kind,2.033300e+00_fp_kind,&
     1.749700e+00_fp_kind,1.367700e+00_fp_kind,1.174600e+00_fp_kind,1.028600e+00_fp_kind,8.951000e-01_fp_kind,&
     6.578000e-01_fp_kind,4.714000e-01_fp_kind], 9.1590298241636958e-01_fp_kind,3.430000e+00_fp_kind) /
data ionicXrayFF( 1, 20) / t_ionicFF('Ca',  20,  0, &
     [2.000000e+01_fp_kind,1.909150e+01_fp_kind,1.733210e+01_fp_kind,1.572450e+01_fp_kind,1.430630e+01_fp_kind,&
     1.296340e+01_fp_kind,1.170740e+01_fp_kind,1.059240e+01_fp_kind,9.651800e+00_fp_kind,8.886900e+00_fp_kind,&
     8.276600e+00_fp_kind,7.393000e+00_fp_kind,6.763300e+00_fp_kind,6.229500e+00_fp_kind,5.718600e+00_fp_kind,&
     5.210800e+00_fp_kind,4.234800e+00_fp_kind,3.392600e+00_fp_kind,2.733400e+00_fp_kind,2.250300e+00_fp_kind,&
     1.909800e+00_fp_kind,1.444600e+00_fp_kind,1.226600e+00_fp_kind,1.080600e+00_fp_kind,9.458000e-01_fp_kind,&
     7.189000e-01_fp_kind,5.229000e-01_fp_kind], 9.6648885100965032e-01_fp_kind,9.837100e+00_fp_kind) /
data ionicXrayFF( 2, 20) / t_ionicFF('Ca',  20,  2, &
     [1.800000e+01_fp_kind,1.772140e+01_fp_kind,1.693550e+01_fp_kind,1.577550e+01_fp_kind,1.441440e+01_fp_kind,&
     1.301880e+01_fp_kind,1.171490e+01_fp_kind,1.057660e+01_fp_kind,9.630400e+00_fp_kind,8.868500e+00_fp_kind,&
     8.264000e+00_fp_kind,7.390600e+00_fp_kind,6.765400e+00_fp_kind,6.232300e+00_fp_kind,5.720300e+00_fp_kind,&
     5.211300e+00_fp_kind,4.233800e+00_fp_kind,3.391600e+00_fp_kind,2.732900e+00_fp_kind,2.250400e+00_fp_kind,&
     1.910200e+00_fp_kind,1.444600e+00_fp_kind,1.225900e+00_fp_kind,1.078200e+00_fp_kind,9.487000e-01_fp_kind,&
     7.150000e-01_fp_kind,5.244000e-01_fp_kind], 9.6933598964858447e-01_fp_kind,2.707500e+00_fp_kind) /
data ionicXrayFF( 1, 21) / t_ionicFF('Sc',  21,  0, &
     [2.100000e+01_fp_kind,2.013120e+01_fp_kind,1.835420e+01_fp_kind,1.664160e+01_fp_kind,1.513150e+01_fp_kind,&
     1.373020e+01_fp_kind,1.242280e+01_fp_kind,1.124490e+01_fp_kind,1.022740e+01_fp_kind,9.379300e+00_fp_kind,&
     8.689100e+00_fp_kind,7.683900e+00_fp_kind,6.998000e+00_fp_kind,6.461400e+00_fp_kind,5.976700e+00_fp_kind,&
     5.502900e+00_fp_kind,4.571900e+00_fp_kind,3.724400e+00_fp_kind,3.023900e+00_fp_kind,2.485700e+00_fp_kind,&
     2.091400e+00_fp_kind,1.534000e+00_fp_kind,1.280400e+00_fp_kind,1.127300e+00_fp_kind,9.954000e-01_fp_kind,&
     7.737000e-01_fp_kind,5.747000e-01_fp_kind], 1.0129202508653561e+00_fp_kind,9.250200e+00_fp_kind) /
data ionicXrayFF( 2, 21) / t_ionicFF('Sc',  21,  3, &
     [1.800000e+01_fp_kind,1.777210e+01_fp_kind,1.712070e+01_fp_kind,1.613470e+01_fp_kind,1.493630e+01_fp_kind,&
     1.365340e+01_fp_kind,1.239610e+01_fp_kind,1.124290e+01_fp_kind,1.023750e+01_fp_kind,9.393500e+00_fp_kind,&
     8.702800e+00_fp_kind,7.692100e+00_fp_kind,7.002700e+00_fp_kind,6.466500e+00_fp_kind,5.984000e+00_fp_kind,&
     5.512800e+00_fp_kind,4.584300e+00_fp_kind,3.735900e+00_fp_kind,3.033000e+00_fp_kind,2.492300e+00_fp_kind,&
     2.095900e+00_fp_kind,1.535000e+00_fp_kind,1.279800e+00_fp_kind,1.125000e+00_fp_kind,9.979000e-01_fp_kind,&
     7.694000e-01_fp_kind,5.765000e-01_fp_kind], 1.0161823389722640e+00_fp_kind,2.207800e+00_fp_kind) /
data ionicXrayFF( 1, 22) / t_ionicFF('Ti',  22,  0, &
     [2.200000e+01_fp_kind,2.117090e+01_fp_kind,1.940750e+01_fp_kind,1.762970e+01_fp_kind,1.603740e+01_fp_kind,&
     1.456670e+01_fp_kind,1.319530e+01_fp_kind,1.194770e+01_fp_kind,1.085200e+01_fp_kind,9.920800e+00_fp_kind,&
     9.149100e+00_fp_kind,8.008100e+00_fp_kind,7.241800e+00_fp_kind,6.677800e+00_fp_kind,6.201700e+00_fp_kind,&
     5.753800e+00_fp_kind,4.873900e+00_fp_kind,4.040500e+00_fp_kind,3.317100e+00_fp_kind,2.735300e+00_fp_kind,&
     2.292100e+00_fp_kind,1.637900e+00_fp_kind,1.339200e+00_fp_kind,1.172900e+00_fp_kind,1.041500e+00_fp_kind,&
     8.250000e-01_fp_kind,6.254000e-01_fp_kind], 1.0533249745024520e+00_fp_kind,8.727000e+00_fp_kind) /
data ionicXrayFF( 2, 22) / t_ionicFF('Ti',  22,  4, &
     [1.800000e+01_fp_kind,1.780930e+01_fp_kind,1.725950e+01_fp_kind,1.641230e+01_fp_kind,1.535610e+01_fp_kind,&
     1.418900e+01_fp_kind,1.300270e+01_fp_kind,1.187070e+01_fp_kind,1.084320e+01_fp_kind,9.946500e+00_fp_kind,&
     9.187000e+00_fp_kind,8.039000e+00_fp_kind,7.258400e+00_fp_kind,6.688100e+00_fp_kind,6.213600e+00_fp_kind,&
     5.770900e+00_fp_kind,4.900000e+00_fp_kind,4.068200e+00_fp_kind,3.340700e+00_fp_kind,2.753100e+00_fp_kind,&
     2.304300e+00_fp_kind,1.641200e+00_fp_kind,1.339000e+00_fp_kind,1.170600e+00_fp_kind,1.043500e+00_fp_kind,&
     8.206000e-01_fp_kind,6.274000e-01_fp_kind], 1.0567923415962570e+00_fp_kind,1.842900e+00_fp_kind) /
data ionicXrayFF( 1, 23) / t_ionicFF('V ',  23,  0, &
     [2.300000e+01_fp_kind,2.220800e+01_fp_kind,2.047160e+01_fp_kind,1.865380e+01_fp_kind,1.699420e+01_fp_kind,&
     1.545670e+01_fp_kind,1.401970e+01_fp_kind,1.270150e+01_fp_kind,1.152830e+01_fp_kind,1.051470e+01_fp_kind,&
     9.660400e+00_fp_kind,8.373600e+00_fp_kind,7.507400e+00_fp_kind,6.893700e+00_fp_kind,6.407500e+00_fp_kind,&
     5.973600e+00_fp_kind,5.141500e+00_fp_kind,4.335500e+00_fp_kind,3.605600e+00_fp_kind,2.992800e+00_fp_kind,&
     2.507600e+00_fp_kind,1.756900e+00_fp_kind,1.405100e+00_fp_kind,1.219200e+00_fp_kind,1.085100e+00_fp_kind,&
     8.729000e-01_fp_kind,6.747000e-01_fp_kind], 1.0879674629232310e+00_fp_kind,8.263300e+00_fp_kind) /
data ionicXrayFF( 2, 23) / t_ionicFF('V ',  23,  5, &
     [1.800000e+01_fp_kind,1.783770e+01_fp_kind,1.736670e+01_fp_kind,1.663140e+01_fp_kind,1.569720e+01_fp_kind,&
     1.463990e+01_fp_kind,1.353410e+01_fp_kind,1.244560e+01_fp_kind,1.142390e+01_fp_kind,1.050160e+01_fp_kind,&
     9.694800e+00_fp_kind,8.427300e+00_fp_kind,7.544000e+00_fp_kind,6.915300e+00_fp_kind,6.425300e+00_fp_kind,&
     5.996100e+00_fp_kind,5.179100e+00_fp_kind,4.379900e+00_fp_kind,3.646800e+00_fp_kind,3.025900e+00_fp_kind,&
     2.531700e+00_fp_kind,1.764700e+00_fp_kind,1.406000e+00_fp_kind,1.217100e+00_fp_kind,1.086500e+00_fp_kind,&
     8.686000e-01_fp_kind,6.768000e-01_fp_kind], 1.0914564220183460e+00_fp_kind,1.565800e+00_fp_kind) /
data ionicXrayFF( 1, 24) / t_ionicFF('Cr',  24,  0, &
     [2.400000e+01_fp_kind,2.332680e+01_fp_kind,2.178160e+01_fp_kind,2.000930e+01_fp_kind,1.824580e+01_fp_kind,&
     1.654800e+01_fp_kind,1.495590e+01_fp_kind,1.350620e+01_fp_kind,1.222330e+01_fp_kind,1.111600e+01_fp_kind,&
     1.017950e+01_fp_kind,8.756300e+00_fp_kind,7.791900e+00_fp_kind,7.119000e+00_fp_kind,6.607900e+00_fp_kind,&
     6.173800e+00_fp_kind,5.374700e+00_fp_kind,4.599800e+00_fp_kind,3.876200e+00_fp_kind,3.245900e+00_fp_kind,&
     2.728900e+00_fp_kind,1.889400e+00_fp_kind,1.479500e+00_fp_kind,1.267700e+00_fp_kind,1.127800e+00_fp_kind,&
     9.164000e-01_fp_kind,7.225000e-01_fp_kind], 1.1173880356567449e+00_fp_kind,6.955000e+00_fp_kind) /
data ionicXrayFF( 2, 24) / t_ionicFF('Cr',  24,  4, &
     [2.000000e+01_fp_kind,1.980150e+01_fp_kind,1.922780e+01_fp_kind,1.833890e+01_fp_kind,1.722070e+01_fp_kind,&
     1.596890e+01_fp_kind,1.467430e+01_fp_kind,1.341210e+01_fp_kind,1.223700e+01_fp_kind,1.118240e+01_fp_kind,&
     1.026350e+01_fp_kind,8.824400e+00_fp_kind,7.828100e+00_fp_kind,7.133200e+00_fp_kind,6.613600e+00_fp_kind,&
     6.180400e+00_fp_kind,5.392700e+00_fp_kind,4.626000e+00_fp_kind,3.903300e+00_fp_kind,3.269400e+00_fp_kind,&
     2.747000e+00_fp_kind,1.896200e+00_fp_kind,1.480800e+00_fp_kind,1.266600e+00_fp_kind,1.128900e+00_fp_kind,&
     9.136000e-01_fp_kind,7.240000e-01_fp_kind], 1.1197800309331500e+00_fp_kind,1.916700e+00_fp_kind) /
data ionicXrayFF( 1, 25) / t_ionicFF('Mn',  25,  0, &
     [2.500000e+01_fp_kind,2.427390e+01_fp_kind,2.260790e+01_fp_kind,2.075690e+01_fp_kind,1.900210e+01_fp_kind,&
     1.735260e+01_fp_kind,1.579550e+01_fp_kind,1.434500e+01_fp_kind,1.302510e+01_fp_kind,1.185350e+01_fp_kind,&
     1.083690e+01_fp_kind,9.244100e+00_fp_kind,8.137200e+00_fp_kind,7.369200e+00_fp_kind,6.808800e+00_fp_kind,&
     6.360400e+00_fp_kind,5.587600e+00_fp_kind,4.851500e+00_fp_kind,4.145600e+00_fp_kind,3.507500e+00_fp_kind,&
     2.964700e+00_fp_kind,2.038400e+00_fp_kind,1.564000e+00_fp_kind,1.320500e+00_fp_kind,1.169100e+00_fp_kind,&
     9.596000e-01_fp_kind,7.675000e-01_fp_kind], 1.1402042711234801e+00_fp_kind,7.477200e+00_fp_kind) /
data ionicXrayFF( 2, 25) / t_ionicFF('Mn',  25,  2, &
     [2.300000e+01_fp_kind,2.270580e+01_fp_kind,2.187130e+01_fp_kind,2.062220e+01_fp_kind,1.911740e+01_fp_kind,&
     1.750770e+01_fp_kind,1.591140e+01_fp_kind,1.440950e+01_fp_kind,1.304930e+01_fp_kind,1.185280e+01_fp_kind,&
     1.082370e+01_fp_kind,9.227500e+00_fp_kind,8.127500e+00_fp_kind,7.366100e+00_fp_kind,6.809300e+00_fp_kind,&
     6.361800e+00_fp_kind,5.587000e+00_fp_kind,4.848500e+00_fp_kind,4.141600e+00_fp_kind,3.503800e+00_fp_kind,&
     2.961900e+00_fp_kind,2.037700e+00_fp_kind,1.563900e+00_fp_kind,1.319600e+00_fp_kind,1.171000e+00_fp_kind,&
     9.560000e-01_fp_kind,7.691000e-01_fp_kind], 1.1426566904242179e+00_fp_kind,2.855400e+00_fp_kind) /
data ionicXrayFF( 1, 26) / t_ionicFF('Fe',  26,  0, &
     [2.600000e+01_fp_kind,2.530310e+01_fp_kind,2.367520e+01_fp_kind,2.182280e+01_fp_kind,2.003620e+01_fp_kind,&
     1.834230e+01_fp_kind,1.673350e+01_fp_kind,1.522390e+01_fp_kind,1.383700e+01_fp_kind,1.259190e+01_fp_kind,&
     1.149770e+01_fp_kind,9.751100e+00_fp_kind,8.511400e+00_fp_kind,7.645300e+00_fp_kind,7.023900e+00_fp_kind,&
     6.546600e+00_fp_kind,5.776900e+00_fp_kind,5.072700e+00_fp_kind,4.389900e+00_fp_kind,3.754800e+00_fp_kind,&
     3.196900e+00_fp_kind,2.198000e+00_fp_kind,1.658700e+00_fp_kind,1.378300e+00_fp_kind,1.211900e+00_fp_kind,&
     9.990000e-01_fp_kind,8.108000e-01_fp_kind], 1.1587823352865529e+00_fp_kind,7.140300e+00_fp_kind) /
data ionicXrayFF( 2, 26) / t_ionicFF('Fe',  26,  2, &
     [2.400000e+01_fp_kind,2.371010e+01_fp_kind,2.288500e+01_fp_kind,2.164180e+01_fp_kind,2.013100e+01_fp_kind,&
     1.849850e+01_fp_kind,1.686160e+01_fp_kind,1.530310e+01_fp_kind,1.387400e+01_fp_kind,1.260020e+01_fp_kind,&
     1.148990e+01_fp_kind,9.734700e+00_fp_kind,8.499900e+00_fp_kind,7.640300e+00_fp_kind,7.023300e+00_fp_kind,&
     6.547800e+00_fp_kind,5.777000e+00_fp_kind,5.070300e+00_fp_kind,4.386100e+00_fp_kind,3.750900e+00_fp_kind,&
     3.193600e+00_fp_kind,2.196900e+00_fp_kind,1.658400e+00_fp_kind,1.377600e+00_fp_kind,1.213800e+00_fp_kind,&
     9.957000e-01_fp_kind,8.123000e-01_fp_kind], 1.1609952476804111e+00_fp_kind,2.810900e+00_fp_kind) /
data ionicXrayFF( 1, 27) / t_ionicFF('Co',  27,  0, &
     [2.700000e+01_fp_kind,2.633040e+01_fp_kind,2.474180e+01_fp_kind,2.289540e+01_fp_kind,2.108620e+01_fp_kind,&
     1.935540e+01_fp_kind,1.770110e+01_fp_kind,1.613760e+01_fp_kind,1.468840e+01_fp_kind,1.337380e+01_fp_kind,&
     1.220520e+01_fp_kind,1.030670e+01_fp_kind,8.929500e+00_fp_kind,7.954900e+00_fp_kind,7.259200e+00_fp_kind,&
     6.739300e+00_fp_kind,5.952100e+00_fp_kind,5.272000e+00_fp_kind,4.615100e+00_fp_kind,3.990500e+00_fp_kind,&
     3.425900e+00_fp_kind,2.367400e+00_fp_kind,1.763800e+00_fp_kind,1.442400e+00_fp_kind,1.256800e+00_fp_kind,&
     1.036500e+00_fp_kind,8.518000e-01_fp_kind], 1.1727308189473860e+00_fp_kind,6.830500e+00_fp_kind) /
data ionicXrayFF( 2, 27) / t_ionicFF('Co',  27,  2, &
     [2.500000e+01_fp_kind,2.471520e+01_fp_kind,2.390180e+01_fp_kind,2.266850e+01_fp_kind,2.115750e+01_fp_kind,&
     1.950870e+01_fp_kind,1.783820e+01_fp_kind,1.622990e+01_fp_kind,1.473780e+01_fp_kind,1.339170e+01_fp_kind,&
     1.220370e+01_fp_kind,1.029140e+01_fp_kind,8.916500e+00_fp_kind,7.948100e+00_fp_kind,7.257300e+00_fp_kind,&
     6.739900e+00_fp_kind,5.952700e+00_fp_kind,5.270300e+00_fp_kind,4.611500e+00_fp_kind,3.986500e+00_fp_kind,&
     3.422300e+00_fp_kind,2.365800e+00_fp_kind,1.763400e+00_fp_kind,1.441800e+00_fp_kind,1.258500e+00_fp_kind,&
     1.033400e+00_fp_kind,8.533000e-01_fp_kind], 1.1748633670788351e+00_fp_kind,2.759500e+00_fp_kind) /
data ionicXrayFF( 1, 28) / t_ionicFF('Ni',  28,  0, &
     [2.800000e+01_fp_kind,2.735570e+01_fp_kind,2.580620e+01_fp_kind,2.397100e+01_fp_kind,2.214610e+01_fp_kind,&
     2.038500e+01_fp_kind,1.869110e+01_fp_kind,1.707910e+01_fp_kind,1.557280e+01_fp_kind,1.419340e+01_fp_kind,&
     1.295430e+01_fp_kind,1.090810e+01_fp_kind,9.391500e+00_fp_kind,8.301000e+00_fp_kind,7.520100e+00_fp_kind,&
     6.945300e+00_fp_kind,6.119600e+00_fp_kind,5.453000e+00_fp_kind,4.821200e+00_fp_kind,4.212400e+00_fp_kind,&
     3.648600e+00_fp_kind,2.544000e+00_fp_kind,1.878900e+00_fp_kind,1.513600e+00_fp_kind,1.304800e+00_fp_kind,&
     1.072600e+00_fp_kind,8.907000e-01_fp_kind], 1.1828118025917360e+00_fp_kind,6.546600e+00_fp_kind) /
data ionicXrayFF( 2, 28) / t_ionicFF('Ni',  28,  2, &
     [2.600000e+01_fp_kind,2.572050e+01_fp_kind,2.492000e+01_fp_kind,2.369920e+01_fp_kind,2.219190e+01_fp_kind,&
     2.053240e+01_fp_kind,1.883450e+01_fp_kind,1.718290e+01_fp_kind,1.563400e+01_fp_kind,1.422110e+01_fp_kind,&
     1.295970e+01_fp_kind,1.089460e+01_fp_kind,9.377700e+00_fp_kind,8.292600e+00_fp_kind,7.516800e+00_fp_kind,&
     6.945200e+00_fp_kind,6.120500e+00_fp_kind,5.451900e+00_fp_kind,4.818100e+00_fp_kind,4.208500e+00_fp_kind,&
     3.644700e+00_fp_kind,2.542100e+00_fp_kind,1.878300e+00_fp_kind,1.513000e+00_fp_kind,1.306400e+00_fp_kind,&
     1.069700e+00_fp_kind,8.921000e-01_fp_kind], 1.1847321260591950e+00_fp_kind,2.705700e+00_fp_kind) /
data ionicXrayFF( 1, 29) / t_ionicFF('Cu',  29,  0, &
     [2.900000e+01_fp_kind,2.837930e+01_fp_kind,2.686820e+01_fp_kind,2.504790e+01_fp_kind,2.321320e+01_fp_kind,&
     2.142780e+01_fp_kind,1.969970e+01_fp_kind,1.804470e+01_fp_kind,1.648630e+01_fp_kind,1.504680e+01_fp_kind,&
     1.374140e+01_fp_kind,1.155280e+01_fp_kind,9.897100e+00_fp_kind,8.685800e+00_fp_kind,7.810700e+00_fp_kind,&
     7.170200e+00_fp_kind,6.285400e+00_fp_kind,5.619800e+00_fp_kind,5.009500e+00_fp_kind,4.419500e+00_fp_kind,&
     3.862100e+00_fp_kind,2.725500e+00_fp_kind,2.003300e+00_fp_kind,1.592200e+00_fp_kind,1.356800e+00_fp_kind,&
     1.107900e+00_fp_kind,9.276000e-01_fp_kind], 1.1895527279463001e+00_fp_kind,6.285100e+00_fp_kind) /
data ionicXrayFF( 2, 29) / t_ionicFF('Cu',  29,  2, &
     [2.700000e+01_fp_kind,2.672600e+01_fp_kind,2.593900e+01_fp_kind,2.473250e+01_fp_kind,2.323230e+01_fp_kind,&
     2.156680e+01_fp_kind,1.984710e+01_fp_kind,1.815820e+01_fp_kind,1.655860e+01_fp_kind,1.508440e+01_fp_kind,&
     1.375420e+01_fp_kind,1.154200e+01_fp_kind,9.883000e+00_fp_kind,8.675900e+00_fp_kind,7.805900e+00_fp_kind,&
     7.169100e+00_fp_kind,6.286500e+00_fp_kind,5.619300e+00_fp_kind,5.007000e+00_fp_kind,4.415800e+00_fp_kind,&
     3.858200e+00_fp_kind,2.723100e+00_fp_kind,2.002500e+00_fp_kind,1.591600e+00_fp_kind,1.358300e+00_fp_kind,&
     1.105200e+00_fp_kind,9.290000e-01_fp_kind], 1.1914075024046080e+00_fp_kind,2.650700e+00_fp_kind) /
data ionicXrayFF( 1, 30) / t_ionicFF('Zn',  30,  0, &
     [3.000000e+01_fp_kind,2.960430e+01_fp_kind,2.792780e+01_fp_kind,2.612540e+01_fp_kind,2.428580e+01_fp_kind,&
     2.248120e+01_fp_kind,2.072400e+01_fp_kind,1.903080e+01_fp_kind,1.742530e+01_fp_kind,1.593050e+01_fp_kind,&
     1.456300e+01_fp_kind,1.223840e+01_fp_kind,1.044530e+01_fp_kind,9.110200e+00_fp_kind,8.134000e+00_fp_kind,&
     7.418600e+00_fp_kind,6.455300e+00_fp_kind,5.776600e+00_fp_kind,5.182000e+00_fp_kind,4.611300e+00_fp_kind,&
     4.064900e+00_fp_kind,2.909100e+00_fp_kind,2.135900e+00_fp_kind,1.678500e+00_fp_kind,1.413600e+00_fp_kind,&
     1.143000e+00_fp_kind,9.627000e-01_fp_kind], 1.1935407217613230e+00_fp_kind,6.043400e+00_fp_kind) /
data ionicXrayFF( 2, 30) / t_ionicFF('Zn',  30,  2, &
     [2.800000e+01_fp_kind,2.773150e+01_fp_kind,2.695840e+01_fp_kind,2.576750e+01_fp_kind,2.427700e+01_fp_kind,&
     2.260940e+01_fp_kind,2.087310e+01_fp_kind,1.915260e+01_fp_kind,1.750800e+01_fp_kind,1.597780e+01_fp_kind,&
     1.458360e+01_fp_kind,1.223090e+01_fp_kind,1.043150e+01_fp_kind,9.099200e+00_fp_kind,8.127800e+00_fp_kind,&
     7.416400e+00_fp_kind,6.456300e+00_fp_kind,5.776700e+00_fp_kind,5.180100e+00_fp_kind,4.608000e+00_fp_kind,&
     4.061000e+00_fp_kind,2.906400e+00_fp_kind,2.134700e+00_fp_kind,1.677800e+00_fp_kind,1.415000e+00_fp_kind,&
     1.140600e+00_fp_kind,9.641000e-01_fp_kind], 1.1953340519839131e+00_fp_kind,2.595700e+00_fp_kind) /
data ionicXrayFF( 1, 31) / t_ionicFF('Ga',  31,  0, &
     [3.100000e+01_fp_kind,3.030290e+01_fp_kind,2.866700e+01_fp_kind,2.677800e+01_fp_kind,2.493460e+01_fp_kind,&
     2.317710e+01_fp_kind,2.148600e+01_fp_kind,1.985190e+01_fp_kind,1.828360e+01_fp_kind,1.679870e+01_fp_kind,&
     1.541500e+01_fp_kind,1.299990e+01_fp_kind,1.107650e+01_fp_kind,9.606600e+00_fp_kind,8.512100e+00_fp_kind,&
     7.704100e+00_fp_kind,6.634600e+00_fp_kind,5.927900e+00_fp_kind,5.344000e+00_fp_kind,4.793700e+00_fp_kind,&
     4.262400e+00_fp_kind,3.098600e+00_fp_kind,2.278600e+00_fp_kind,1.774100e+00_fp_kind,1.475500e+00_fp_kind,&
     1.179300e+00_fp_kind,9.957000e-01_fp_kind], 1.1946687641438061e+00_fp_kind,7.141300e+00_fp_kind) /
data ionicXrayFF( 2, 31) / t_ionicFF('Ga',  31,  0, &
     [3.100000e+01_fp_kind,3.030290e+01_fp_kind,2.866700e+01_fp_kind,2.677800e+01_fp_kind,2.493460e+01_fp_kind,&
     2.317710e+01_fp_kind,2.148600e+01_fp_kind,1.985190e+01_fp_kind,1.828360e+01_fp_kind,1.679870e+01_fp_kind,&
     1.541500e+01_fp_kind,1.299990e+01_fp_kind,1.107650e+01_fp_kind,9.606600e+00_fp_kind,8.512100e+00_fp_kind,&
     7.704100e+00_fp_kind,6.634600e+00_fp_kind,5.927900e+00_fp_kind,5.344000e+00_fp_kind,4.793700e+00_fp_kind,&
     4.262400e+00_fp_kind,3.098600e+00_fp_kind,2.278600e+00_fp_kind,1.774100e+00_fp_kind,1.475500e+00_fp_kind,&
     1.179300e+00_fp_kind,9.957000e-01_fp_kind], 1.1946687641438061e+00_fp_kind,7.141300e+00_fp_kind) /
data ionicXrayFF( 1, 32) / t_ionicFF('Ge',  32,  0, &
     [3.200000e+01_fp_kind,3.127350e+01_fp_kind,2.952730e+01_fp_kind,2.749400e+01_fp_kind,2.555990e+01_fp_kind,&
     2.378940e+01_fp_kind,2.213890e+01_fp_kind,2.056600e+01_fp_kind,1.905370e+01_fp_kind,1.760540e+01_fp_kind,&
     1.623380e+01_fp_kind,1.377520e+01_fp_kind,1.174850e+01_fp_kind,1.015350e+01_fp_kind,8.938600e+00_fp_kind,&
     8.029400e+00_fp_kind,6.831200e+00_fp_kind,6.078100e+00_fp_kind,5.494900e+00_fp_kind,4.962700e+00_fp_kind,&
     4.448600e+00_fp_kind,3.288300e+00_fp_kind,2.428900e+00_fp_kind,1.877900e+00_fp_kind,1.544200e+00_fp_kind,&
     1.216500e+00_fp_kind,1.027500e+00_fp_kind], 1.1942852530470560e+00_fp_kind,7.374100e+00_fp_kind) /
data ionicXrayFF( 2, 32) / t_ionicFF('Ge',  32,  0, &
     [3.200000e+01_fp_kind,3.127350e+01_fp_kind,2.952730e+01_fp_kind,2.749400e+01_fp_kind,2.555990e+01_fp_kind,&
     2.378940e+01_fp_kind,2.213890e+01_fp_kind,2.056600e+01_fp_kind,1.905370e+01_fp_kind,1.760540e+01_fp_kind,&
     1.623380e+01_fp_kind,1.377520e+01_fp_kind,1.174850e+01_fp_kind,1.015350e+01_fp_kind,8.938600e+00_fp_kind,&
     8.029400e+00_fp_kind,6.831200e+00_fp_kind,6.078100e+00_fp_kind,5.494900e+00_fp_kind,4.962700e+00_fp_kind,&
     4.448600e+00_fp_kind,3.288300e+00_fp_kind,2.428900e+00_fp_kind,1.877900e+00_fp_kind,1.544200e+00_fp_kind,&
     1.216500e+00_fp_kind,1.027500e+00_fp_kind], 1.1942852530470560e+00_fp_kind,7.374100e+00_fp_kind) /
data ionicXrayFF( 1, 33) / t_ionicFF('As',  33,  0, &
     [3.300000e+01_fp_kind,3.226770e+01_fp_kind,3.045710e+01_fp_kind,2.828870e+01_fp_kind,2.622230e+01_fp_kind,&
     2.438060e+01_fp_kind,2.272520e+01_fp_kind,2.119070e+01_fp_kind,1.973290e+01_fp_kind,1.833450e+01_fp_kind,&
     1.699620e+01_fp_kind,1.454050e+01_fp_kind,1.244530e+01_fp_kind,1.074320e+01_fp_kind,9.412400e+00_fp_kind,&
     8.397800e+00_fp_kind,7.051200e+00_fp_kind,6.232500e+00_fp_kind,5.637900e+00_fp_kind,5.118800e+00_fp_kind,&
     4.622600e+00_fp_kind,3.476000e+00_fp_kind,2.585200e+00_fp_kind,1.989800e+00_fp_kind,1.619900e+00_fp_kind,&
     1.255000e+00_fp_kind,1.058700e+00_fp_kind], 1.1932263245390831e+00_fp_kind,7.368600e+00_fp_kind) /
data ionicXrayFF( 2, 33) / t_ionicFF('As',  33,  0, &
     [3.300000e+01_fp_kind,3.226770e+01_fp_kind,3.045710e+01_fp_kind,2.828870e+01_fp_kind,2.622230e+01_fp_kind,&
     2.438060e+01_fp_kind,2.272520e+01_fp_kind,2.119070e+01_fp_kind,1.973290e+01_fp_kind,1.833450e+01_fp_kind,&
     1.699620e+01_fp_kind,1.454050e+01_fp_kind,1.244530e+01_fp_kind,1.074320e+01_fp_kind,9.412400e+00_fp_kind,&
     8.397800e+00_fp_kind,7.051200e+00_fp_kind,6.232500e+00_fp_kind,5.637900e+00_fp_kind,5.118800e+00_fp_kind,&
     4.622600e+00_fp_kind,3.476000e+00_fp_kind,2.585200e+00_fp_kind,1.989800e+00_fp_kind,1.619900e+00_fp_kind,&
     1.255000e+00_fp_kind,1.058700e+00_fp_kind], 1.1932263245390831e+00_fp_kind,7.368600e+00_fp_kind) /
data ionicXrayFF( 1, 34) / t_ionicFF('Se',  34,  0, &
     [3.400000e+01_fp_kind,3.327270e+01_fp_kind,3.142970e+01_fp_kind,2.915120e+01_fp_kind,2.694250e+01_fp_kind,&
     2.499050e+01_fp_kind,2.328530e+01_fp_kind,2.175510e+01_fp_kind,2.033470e+01_fp_kind,1.898510e+01_fp_kind,&
     1.769070e+01_fp_kind,1.527550e+01_fp_kind,1.314900e+01_fp_kind,1.136460e+01_fp_kind,9.929100e+00_fp_kind,&
     8.809700e+00_fp_kind,7.299900e+00_fp_kind,6.396400e+00_fp_kind,5.776600e+00_fp_kind,5.263900e+00_fp_kind,&
     4.784200e+00_fp_kind,3.659800e+00_fp_kind,2.746200e+00_fp_kind,2.109500e+00_fp_kind,1.702700e+00_fp_kind,&
     1.295700e+00_fp_kind,1.098500e+00_fp_kind], 1.2019512788170770e+00_fp_kind,7.266100e+00_fp_kind) /
data ionicXrayFF( 2, 34) / t_ionicFF('Se',  34,  0, &
     [3.400000e+01_fp_kind,3.327270e+01_fp_kind,3.142970e+01_fp_kind,2.915120e+01_fp_kind,2.694250e+01_fp_kind,&
     2.499050e+01_fp_kind,2.328530e+01_fp_kind,2.175510e+01_fp_kind,2.033470e+01_fp_kind,1.898510e+01_fp_kind,&
     1.769070e+01_fp_kind,1.527550e+01_fp_kind,1.314900e+01_fp_kind,1.136460e+01_fp_kind,9.929100e+00_fp_kind,&
     8.809700e+00_fp_kind,7.299900e+00_fp_kind,6.396400e+00_fp_kind,5.776600e+00_fp_kind,5.263900e+00_fp_kind,&
     4.784200e+00_fp_kind,3.659800e+00_fp_kind,2.746200e+00_fp_kind,2.109500e+00_fp_kind,1.702700e+00_fp_kind,&
     1.295700e+00_fp_kind,1.098500e+00_fp_kind], 1.2019512788170770e+00_fp_kind,7.266100e+00_fp_kind) /
data ionicXrayFF( 1, 35) / t_ionicFF('Br',  35,  0, &
     [3.500000e+01_fp_kind,3.428700e+01_fp_kind,3.244030e+01_fp_kind,3.008280e+01_fp_kind,2.773790e+01_fp_kind,&
     2.565130e+01_fp_kind,2.385580e+01_fp_kind,2.229030e+01_fp_kind,2.087910e+01_fp_kind,1.956470e+01_fp_kind,&
     1.831390e+01_fp_kind,1.596400e+01_fp_kind,1.384190e+01_fp_kind,1.200470e+01_fp_kind,1.048190e+01_fp_kind,&
     9.263900e+00_fp_kind,7.582000e+00_fp_kind,6.575500e+00_fp_kind,5.915100e+00_fp_kind,5.400000e+00_fp_kind,&
     4.933900e+00_fp_kind,3.837700e+00_fp_kind,2.910200e+00_fp_kind,2.236500e+00_fp_kind,1.792900e+00_fp_kind,&
     1.339500e+00_fp_kind,1.120600e+00_fp_kind], 1.1907412764098591e+00_fp_kind,7.080500e+00_fp_kind) /
data ionicXrayFF( 2, 35) / t_ionicFF('Br',  35, -1, &
     [3.600000e+01_fp_kind,3.511460e+01_fp_kind,3.289970e+01_fp_kind,3.022740e+01_fp_kind,2.772240e+01_fp_kind,&
     2.559330e+01_fp_kind,2.380690e+01_fp_kind,2.226140e+01_fp_kind,2.086620e+01_fp_kind,1.956160e+01_fp_kind,&
     1.831600e+01_fp_kind,1.596880e+01_fp_kind,1.384540e+01_fp_kind,1.200660e+01_fp_kind,1.048270e+01_fp_kind,&
     9.264100e+00_fp_kind,7.581900e+00_fp_kind,6.575400e+00_fp_kind,5.915200e+00_fp_kind,5.400100e+00_fp_kind,&
     4.934100e+00_fp_kind,3.838000e+00_fp_kind,2.910600e+00_fp_kind,2.237100e+00_fp_kind,1.792500e+00_fp_kind,&
     1.341100e+00_fp_kind,1.119600e+00_fp_kind], 1.1896435697335330e+00_fp_kind,8.873600e+00_fp_kind) /
data ionicXrayFF( 1, 36) / t_ionicFF('Kr',  36,  0, &
     [3.600000e+01_fp_kind,3.530380e+01_fp_kind,3.346800e+01_fp_kind,3.105710e+01_fp_kind,2.859320e+01_fp_kind,&
     2.636700e+01_fp_kind,2.445660e+01_fp_kind,2.282330e+01_fp_kind,2.139130e+01_fp_kind,2.009120e+01_fp_kind,&
     1.887430e+01_fp_kind,1.659820e+01_fp_kind,1.450870e+01_fp_kind,1.264840e+01_fp_kind,1.106010e+01_fp_kind,&
     9.754900e+00_fp_kind,7.899800e+00_fp_kind,6.774800e+00_fp_kind,6.058200e+00_fp_kind,5.530100e+00_fp_kind,&
     5.072900e+00_fp_kind,4.008300e+00_fp_kind,3.075700e+00_fp_kind,2.369900e+00_fp_kind,1.890100e+00_fp_kind,&
     1.387000e+00_fp_kind,1.152100e+00_fp_kind], 1.1901893657867431e+00_fp_kind,6.881100e+00_fp_kind) /
data ionicXrayFF( 2, 36) / t_ionicFF('Kr',  36,  0, &
     [3.600000e+01_fp_kind,3.530380e+01_fp_kind,3.346800e+01_fp_kind,3.105710e+01_fp_kind,2.859320e+01_fp_kind,&
     2.636700e+01_fp_kind,2.445660e+01_fp_kind,2.282330e+01_fp_kind,2.139130e+01_fp_kind,2.009120e+01_fp_kind,&
     1.887430e+01_fp_kind,1.659820e+01_fp_kind,1.450870e+01_fp_kind,1.264840e+01_fp_kind,1.106010e+01_fp_kind,&
     9.754900e+00_fp_kind,7.899800e+00_fp_kind,6.774800e+00_fp_kind,6.058200e+00_fp_kind,5.530100e+00_fp_kind,&
     5.072900e+00_fp_kind,4.008300e+00_fp_kind,3.075700e+00_fp_kind,2.369900e+00_fp_kind,1.890100e+00_fp_kind,&
     1.387000e+00_fp_kind,1.152100e+00_fp_kind], 1.1901893657867431e+00_fp_kind,6.881100e+00_fp_kind) /
data ionicXrayFF( 1, 37) / t_ionicFF('Rb',  37,  0, &
     [3.700000e+01_fp_kind,3.594840e+01_fp_kind,3.390870e+01_fp_kind,3.168300e+01_fp_kind,2.937130e+01_fp_kind,&
     2.715110e+01_fp_kind,2.516160e+01_fp_kind,2.343590e+01_fp_kind,2.193790e+01_fp_kind,2.060900e+01_fp_kind,&
     1.939490e+01_fp_kind,1.717170e+01_fp_kind,1.513080e+01_fp_kind,1.327610e+01_fp_kind,1.164840e+01_fp_kind,&
     1.027340e+01_fp_kind,8.253900e+00_fp_kind,6.999000e+00_fp_kind,6.210700e+00_fp_kind,5.657400e+00_fp_kind,&
     5.202500e+00_fp_kind,4.170600e+00_fp_kind,3.240700e+00_fp_kind,2.509300e+00_fp_kind,1.991900e+00_fp_kind,&
     1.438700e+00_fp_kind,1.184600e+00_fp_kind], 1.1907056740955040e+00_fp_kind,1.166900e+01_fp_kind) /
data ionicXrayFF( 2, 37) / t_ionicFF('Rb',  37,  1, &
     [3.600000e+01_fp_kind,3.543500e+01_fp_kind,3.389210e+01_fp_kind,3.174280e+01_fp_kind,2.939590e+01_fp_kind,&
     2.715070e+01_fp_kind,2.515330e+01_fp_kind,2.342770e+01_fp_kind,2.193250e+01_fp_kind,2.060660e+01_fp_kind,&
     1.939480e+01_fp_kind,1.717330e+01_fp_kind,1.513220e+01_fp_kind,1.327680e+01_fp_kind,1.164840e+01_fp_kind,&
     1.027320e+01_fp_kind,8.253600e+00_fp_kind,6.998900e+00_fp_kind,6.210700e+00_fp_kind,5.657700e+00_fp_kind,&
     5.202600e+00_fp_kind,4.169900e+00_fp_kind,3.240000e+00_fp_kind,2.507800e+00_fp_kind,1.993800e+00_fp_kind,&
     1.437800e+00_fp_kind,1.185200e+00_fp_kind], 1.1913287244379449e+00_fp_kind,5.534400e+00_fp_kind) /
data ionicXrayFF( 1, 38) / t_ionicFF('Sr',  38,  0, &
     [3.800000e+01_fp_kind,3.680200e+01_fp_kind,3.445930e+01_fp_kind,3.217350e+01_fp_kind,2.999060e+01_fp_kind,&
     2.786660e+01_fp_kind,2.587820e+01_fp_kind,2.409390e+01_fp_kind,2.252580e+01_fp_kind,2.114480e+01_fp_kind,&
     1.990630e+01_fp_kind,1.770030e+01_fp_kind,1.570650e+01_fp_kind,1.387670e+01_fp_kind,1.223390e+01_fp_kind,&
     1.080920e+01_fp_kind,8.642100e+00_fp_kind,7.251400e+00_fp_kind,6.377300e+00_fp_kind,5.786400e+00_fp_kind,&
     5.325500e+00_fp_kind,4.322600e+00_fp_kind,3.403100e+00_fp_kind,2.652100e+00_fp_kind,2.099700e+00_fp_kind,&
     1.497100e+00_fp_kind,1.217500e+00_fp_kind], 1.1915992659552830e+00_fp_kind,1.300520e+01_fp_kind) /
data ionicXrayFF( 2, 38) / t_ionicFF('Sr',  38,  2, &
     [3.600000e+01_fp_kind,3.552440e+01_fp_kind,3.419870e+01_fp_kind,3.228270e+01_fp_kind,3.009160e+01_fp_kind,&
     2.789510e+01_fp_kind,2.586410e+01_fp_kind,2.406800e+01_fp_kind,2.250390e+01_fp_kind,2.113180e+01_fp_kind,&
     1.990170e+01_fp_kind,1.770480e+01_fp_kind,1.571190e+01_fp_kind,1.387980e+01_fp_kind,1.223460e+01_fp_kind,&
     1.080840e+01_fp_kind,8.640900e+00_fp_kind,7.250800e+00_fp_kind,6.377300e+00_fp_kind,5.786800e+00_fp_kind,&
     5.325800e+00_fp_kind,4.321800e+00_fp_kind,3.402000e+00_fp_kind,2.649500e+00_fp_kind,2.103400e+00_fp_kind,&
     1.493900e+00_fp_kind,1.219100e+00_fp_kind], 1.1932171317178160e+00_fp_kind,4.635300e+00_fp_kind) /
data ionicXrayFF( 1, 39) / t_ionicFF('Y ',  39,  0, &
     [3.900000e+01_fp_kind,3.781570e+01_fp_kind,3.536310e+01_fp_kind,3.290370e+01_fp_kind,3.063440e+01_fp_kind,&
     2.849110e+01_fp_kind,2.649230e+01_fp_kind,2.468360e+01_fp_kind,2.308210e+01_fp_kind,2.167110e+01_fp_kind,&
     2.041490e+01_fp_kind,1.821550e+01_fp_kind,1.625800e+01_fp_kind,1.445580e+01_fp_kind,1.281080e+01_fp_kind,&
     1.135200e+01_fp_kind,9.058200e+00_fp_kind,7.532600e+00_fp_kind,6.562300e+00_fp_kind,5.921200e+00_fp_kind,&
     5.444600e+00_fp_kind,4.463100e+00_fp_kind,3.559600e+00_fp_kind,2.795000e+00_fp_kind,2.213900e+00_fp_kind,&
     1.559300e+00_fp_kind,1.251900e+00_fp_kind], 1.1939249922512689e+00_fp_kind,1.259530e+01_fp_kind) /
data ionicXrayFF( 2, 39) / t_ionicFF('Y ',  39,  3, &
     [3.600000e+01_fp_kind,3.559020e+01_fp_kind,3.443250e+01_fp_kind,3.271660e+01_fp_kind,3.068770e+01_fp_kind,&
     2.857760e+01_fp_kind,2.655750e+01_fp_kind,2.472210e+01_fp_kind,2.309990e+01_fp_kind,2.167510e+01_fp_kind,&
     2.041100e+01_fp_kind,1.820660e+01_fp_kind,1.625060e+01_fp_kind,1.445100e+01_fp_kind,1.280780e+01_fp_kind,&
     1.135000e+01_fp_kind,9.057100e+00_fp_kind,7.531900e+00_fp_kind,6.562100e+00_fp_kind,5.921500e+00_fp_kind,&
     5.445200e+00_fp_kind,4.463900e+00_fp_kind,3.560000e+00_fp_kind,2.793700e+00_fp_kind,2.218500e+00_fp_kind,&
     1.555100e+00_fp_kind,1.254500e+00_fp_kind], 1.1964869984501481e+00_fp_kind,3.980200e+00_fp_kind) /
data ionicXrayFF( 1, 40) / t_ionicFF('Zr',  40,  0, &
     [4.000000e+01_fp_kind,3.884580e+01_fp_kind,3.635120e+01_fp_kind,3.375190e+01_fp_kind,3.135910e+01_fp_kind,&
     2.914250e+01_fp_kind,2.709870e+01_fp_kind,2.525310e+01_fp_kind,2.361760e+01_fp_kind,2.217830e+01_fp_kind,&
     2.090410e+01_fp_kind,1.870480e+01_fp_kind,1.677910e+01_fp_kind,1.500920e+01_fp_kind,1.337450e+01_fp_kind,&
     1.189660e+01_fp_kind,9.498700e+00_fp_kind,7.842900e+00_fp_kind,6.768600e+00_fp_kind,6.065800e+00_fp_kind,&
     5.563500e+00_fp_kind,4.593800e+00_fp_kind,3.710200e+00_fp_kind,2.938200e+00_fp_kind,2.332800e+00_fp_kind,&
     1.625600e+00_fp_kind,1.288800e+00_fp_kind], 1.1985368575502731e+00_fp_kind,1.210940e+01_fp_kind) /
data ionicXrayFF( 2, 40) / t_ionicFF('Zr',  40,  4, &
     [3.600000e+01_fp_kind,3.564120e+01_fp_kind,3.461750e+01_fp_kind,3.307230e+01_fp_kind,3.119880e+01_fp_kind,&
     2.919290e+01_fp_kind,2.721490e+01_fp_kind,2.537010e+01_fp_kind,2.370810e+01_fp_kind,2.223520e+01_fp_kind,&
     2.093030e+01_fp_kind,1.869220e+01_fp_kind,1.675520e+01_fp_kind,1.498820e+01_fp_kind,1.336020e+01_fp_kind,&
     1.188820e+01_fp_kind,9.496200e+00_fp_kind,7.841700e+00_fp_kind,6.767900e+00_fp_kind,6.065700e+00_fp_kind,&
     5.564300e+00_fp_kind,4.596400e+00_fp_kind,3.712900e+00_fp_kind,2.938900e+00_fp_kind,2.338300e+00_fp_kind,&
     1.621400e+00_fp_kind,1.291900e+00_fp_kind], 1.2015159617754421e+00_fp_kind,3.477200e+00_fp_kind) /
data ionicXrayFF( 1, 41) / t_ionicFF('Nb',  41,  0, &
     [4.100000e+01_fp_kind,3.996380e+01_fp_kind,3.758640e+01_fp_kind,3.489120e+01_fp_kind,3.228610e+01_fp_kind,&
     2.987310e+01_fp_kind,2.769410e+01_fp_kind,2.576820e+01_fp_kind,2.408930e+01_fp_kind,2.262980e+01_fp_kind,&
     2.135050e+01_fp_kind,1.916960e+01_fp_kind,1.728140e+01_fp_kind,1.554640e+01_fp_kind,1.392770e+01_fp_kind,&
     1.288020e+01_fp_kind,9.957400e+00_fp_kind,8.179700e+00_fp_kind,6.998000e+00_fp_kind,6.223800e+00_fp_kind,&
     5.685800e+00_fp_kind,4.714700e+00_fp_kind,3.852700e+00_fp_kind,3.078900e+00_fp_kind,2.454900e+00_fp_kind,&
     1.695200e+00_fp_kind,1.328800e+00_fp_kind], 1.2058319385347569e+00_fp_kind,1.069910e+01_fp_kind) /
data ionicXrayFF( 2, 41) / t_ionicFF('Nb',  41,  5, &
     [3.600000e+01_fp_kind,3.568190e+01_fp_kind,3.476780e+01_fp_kind,3.336870e+01_fp_kind,3.163910e+01_fp_kind,&
     2.974380e+01_fp_kind,2.782770e+01_fp_kind,2.599770e+01_fp_kind,2.431590e+01_fp_kind,2.280540e+01_fp_kind,&
     2.146060e+01_fp_kind,1.917170e+01_fp_kind,1.723450e+01_fp_kind,1.549350e+01_fp_kind,1.388730e+01_fp_kind,&
     1.241510e+01_fp_kind,9.951200e+00_fp_kind,8.178400e+00_fp_kind,6.996800e+00_fp_kind,6.223000e+00_fp_kind,&
     5.686400e+00_fp_kind,4.720100e+00_fp_kind,3.859700e+00_fp_kind,3.084000e+00_fp_kind,2.462000e+00_fp_kind,&
     1.692800e+00_fp_kind,1.331800e+00_fp_kind], 1.2086457162160120e+00_fp_kind,3.077200e+00_fp_kind) /
data ionicXrayFF( 1, 42) / t_ionicFF('Mo',  42,  0, &
     [4.200000e+01_fp_kind,4.099670e+01_fp_kind,3.863630e+01_fp_kind,3.587880e+01_fp_kind,3.316800e+01_fp_kind,&
     3.064620e+01_fp_kind,2.837120e+01_fp_kind,2.636510e+01_fp_kind,2.462170e+01_fp_kind,2.311400e+01_fp_kind,&
     2.180260e+01_fp_kind,1.960280e+01_fp_kind,1.773820e+01_fp_kind,1.604160e+01_fp_kind,1.445210e+01_fp_kind,&
     1.297150e+01_fp_kind,1.043150e+01_fp_kind,8.543400e+00_fp_kind,7.252000e+00_fp_kind,6.397800e+00_fp_kind,&
     5.814600e+00_fp_kind,4.829200e+00_fp_kind,3.989600e+00_fp_kind,3.219400e+00_fp_kind,2.579600e+00_fp_kind,&
     1.770300e+00_fp_kind,1.371300e+00_fp_kind], 1.2150721042022070e+00_fp_kind,1.028300e+01_fp_kind) /
data ionicXrayFF( 2, 42) / t_ionicFF('Mo',  42,  6, &
     [3.600000e+01_fp_kind,3.571520e+01_fp_kind,3.489240e+01_fp_kind,3.361920e+01_fp_kind,3.202070e+01_fp_kind,&
     3.023570e+01_fp_kind,2.839330e+01_fp_kind,2.659600e+01_fp_kind,2.491250e+01_fp_kind,2.337760e+01_fp_kind,&
     2.199860e+01_fp_kind,1.965170e+01_fp_kind,1.769700e+01_fp_kind,1.597140e+01_fp_kind,1.438770e+01_fp_kind,&
     1.292490e+01_fp_kind,1.041520e+01_fp_kind,8.538600e+00_fp_kind,7.249500e+00_fp_kind,6.396200e+00_fp_kind,&
     5.814800e+00_fp_kind,4.836100e+00_fp_kind,3.999600e+00_fp_kind,3.227700e+00_fp_kind,2.588800e+00_fp_kind,&
     1.769200e+00_fp_kind,1.374500e+00_fp_kind], 1.2180034707265150e+00_fp_kind,2.750800e+00_fp_kind) /
data ionicXrayFF( 1, 43) / t_ionicFF('Tc',  43,  0, &
     [4.300000e+01_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind], 0.0000000000000000e+00_fp_kind,0.000000e+00_fp_kind) /
data ionicXrayFF( 2, 43) / t_ionicFF('Tc',  43,  0, &
     [4.300000e+01_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind], 0.0000000000000000e+00_fp_kind,0.000000e+00_fp_kind) /
data ionicXrayFF( 1, 44) / t_ionicFF('Ru',  44,  0, &
     [4.400000e+01_fp_kind,4.305820e+01_fp_kind,4.076080e+01_fp_kind,3.794410e+01_fp_kind,3.507160e+01_fp_kind,&
     3.234440e+01_fp_kind,2.986030e+01_fp_kind,2.766250e+01_fp_kind,2.575510e+01_fp_kind,2.411620e+01_fp_kind,&
     2.270880e+01_fp_kind,2.041740e+01_fp_kind,1.856630e+01_fp_kind,1.693720e+01_fp_kind,1.542000e+01_fp_kind,&
     1.398250e+01_fp_kind,1.139900e+01_fp_kind,9.335800e+00_fp_kind,7.833900e+00_fp_kind,6.803400e+00_fp_kind,&
     6.104100e+00_fp_kind,5.042500e+00_fp_kind,4.425000e+00_fp_kind,3.491000e+00_fp_kind,2.832400e+00_fp_kind,&
     1.934100e+00_fp_kind,1.465900e+00_fp_kind], 1.2407080436637941e+00_fp_kind,9.552100e+00_fp_kind) /
data ionicXrayFF( 2, 44) / t_ionicFF('Ru',  44,  0, &
     [4.400000e+01_fp_kind,4.305820e+01_fp_kind,4.076080e+01_fp_kind,3.794410e+01_fp_kind,3.507160e+01_fp_kind,&
     3.234440e+01_fp_kind,2.986030e+01_fp_kind,2.766250e+01_fp_kind,2.575510e+01_fp_kind,2.411620e+01_fp_kind,&
     2.270880e+01_fp_kind,2.041740e+01_fp_kind,1.856630e+01_fp_kind,1.693720e+01_fp_kind,1.542000e+01_fp_kind,&
     1.398250e+01_fp_kind,1.139900e+01_fp_kind,9.335800e+00_fp_kind,7.833900e+00_fp_kind,6.803400e+00_fp_kind,&
     6.104100e+00_fp_kind,5.042500e+00_fp_kind,4.425000e+00_fp_kind,3.491000e+00_fp_kind,2.832400e+00_fp_kind,&
     1.934100e+00_fp_kind,1.465900e+00_fp_kind], 1.2407080436637941e+00_fp_kind,9.552100e+00_fp_kind) /
data ionicXrayFF( 1, 45) / t_ionicFF('Rh',  45,  0, &
     [4.500000e+01_fp_kind,4.408700e+01_fp_kind,4.182920e+01_fp_kind,3.900680e+01_fp_kind,3.607790e+01_fp_kind,&
     3.326140e+01_fp_kind,3.067410e+01_fp_kind,2.837260e+01_fp_kind,2.637030e+01_fp_kind,2.465040e+01_fp_kind,&
     2.317900e+01_fp_kind,2.081220e+01_fp_kind,1.894620e+01_fp_kind,1.734040e+01_fp_kind,1.586030e+01_fp_kind,&
     1.445460e+01_fp_kind,1.188130e+01_fp_kind,9.757000e+00_fp_kind,8.160100e+00_fp_kind,7.038000e+00_fp_kind,&
     6.270100e+00_fp_kind,5.145100e+00_fp_kind,4.356700e+00_fp_kind,3.620500e+00_fp_kind,2.958800e+00_fp_kind,&
     2.022100e+00_fp_kind,1.518300e+00_fp_kind], 1.2570529671102970e+00_fp_kind,9.224400e+00_fp_kind) /
data ionicXrayFF( 2, 45) / t_ionicFF('Rh',  45,  0, &
     [4.500000e+01_fp_kind,4.408700e+01_fp_kind,4.182920e+01_fp_kind,3.900680e+01_fp_kind,3.607790e+01_fp_kind,&
     3.326140e+01_fp_kind,3.067410e+01_fp_kind,2.837260e+01_fp_kind,2.637030e+01_fp_kind,2.465040e+01_fp_kind,&
     2.317900e+01_fp_kind,2.081220e+01_fp_kind,1.894620e+01_fp_kind,1.734040e+01_fp_kind,1.586030e+01_fp_kind,&
     1.445460e+01_fp_kind,1.188130e+01_fp_kind,9.757000e+00_fp_kind,8.160100e+00_fp_kind,7.038000e+00_fp_kind,&
     6.270100e+00_fp_kind,5.145100e+00_fp_kind,4.356700e+00_fp_kind,3.620500e+00_fp_kind,2.958800e+00_fp_kind,&
     2.022100e+00_fp_kind,1.518300e+00_fp_kind], 1.2570529671102970e+00_fp_kind,9.224400e+00_fp_kind) /
data ionicXrayFF( 1, 46) / t_ionicFF('Pd',  46,  0, &
     [4.600000e+01_fp_kind,4.523220e+01_fp_kind,4.317450e+01_fp_kind,4.036230e+01_fp_kind,3.729350e+01_fp_kind,&
     3.429360e+01_fp_kind,3.153510e+01_fp_kind,2.909090e+01_fp_kind,2.697490e+01_fp_kind,2.516740e+01_fp_kind,&
     2.363100e+01_fp_kind,2.119090e+01_fp_kind,1.930960e+01_fp_kind,1.772460e+01_fp_kind,1.628030e+01_fp_kind,&
     1.490830e+01_fp_kind,1.235810e+01_fp_kind,1.018850e+01_fp_kind,8.506300e+00_fp_kind,7.294000e+00_fp_kind,&
     6.453100e+00_fp_kind,5.248000e+00_fp_kind,4.464300e+00_fp_kind,3.743400e+00_fp_kind,3.083900e+00_fp_kind,&
     2.112500e+00_fp_kind,1.574400e+00_fp_kind], 1.2758049412951000e+00_fp_kind,7.564500e+00_fp_kind) /
data ionicXrayFF( 2, 46) / t_ionicFF('Pd',  46,  0, &
     [4.600000e+01_fp_kind,4.523220e+01_fp_kind,4.317450e+01_fp_kind,4.036230e+01_fp_kind,3.729350e+01_fp_kind,&
     3.429360e+01_fp_kind,3.153510e+01_fp_kind,2.909090e+01_fp_kind,2.697490e+01_fp_kind,2.516740e+01_fp_kind,&
     2.363100e+01_fp_kind,2.119090e+01_fp_kind,1.930960e+01_fp_kind,1.772460e+01_fp_kind,1.628030e+01_fp_kind,&
     1.490830e+01_fp_kind,1.235810e+01_fp_kind,1.018850e+01_fp_kind,8.506300e+00_fp_kind,7.294000e+00_fp_kind,&
     6.453100e+00_fp_kind,5.248000e+00_fp_kind,4.464300e+00_fp_kind,3.743400e+00_fp_kind,3.083900e+00_fp_kind,&
     2.112500e+00_fp_kind,1.574400e+00_fp_kind], 1.2758049412951000e+00_fp_kind,7.564500e+00_fp_kind) /
data ionicXrayFF( 1, 47) / t_ionicFF('Ag',  47,  0, &
     [4.700000e+01_fp_kind,4.613950e+01_fp_kind,4.396550e+01_fp_kind,4.116040e+01_fp_kind,3.815860e+01_fp_kind,&
     3.519710e+01_fp_kind,3.242130e+01_fp_kind,2.991470e+01_fp_kind,2.771160e+01_fp_kind,2.580990e+01_fp_kind,&
     2.418500e+01_fp_kind,2.161100e+01_fp_kind,1.966520e+01_fp_kind,1.807320e+01_fp_kind,1.665560e+01_fp_kind,&
     1.532100e+01_fp_kind,1.281780e+01_fp_kind,1.062630e+01_fp_kind,8.871600e+00_fp_kind,7.571600e+00_fp_kind,&
     6.653700e+00_fp_kind,5.353300e+00_fp_kind,4.568200e+00_fp_kind,3.864300e+00_fp_kind,3.207500e+00_fp_kind,&
     2.208800e+00_fp_kind,1.633500e+00_fp_kind], 1.2962428223468829e+00_fp_kind,8.641400e+00_fp_kind) /
data ionicXrayFF( 2, 47) / t_ionicFF('Ag',  47,  2, &
     [4.500000e+01_fp_kind,4.446250e+01_fp_kind,4.294470e+01_fp_kind,4.069290e+01_fp_kind,3.801850e+01_fp_kind,&
     3.521150e+01_fp_kind,3.249000e+01_fp_kind,2.998880e+01_fp_kind,2.777160e+01_fp_kind,2.585140e+01_fp_kind,&
     2.420980e+01_fp_kind,2.161450e+01_fp_kind,1.966040e+01_fp_kind,1.806700e+01_fp_kind,1.665060e+01_fp_kind,&
     1.531780e+01_fp_kind,1.281690e+01_fp_kind,1.062610e+01_fp_kind,8.871400e+00_fp_kind,7.571300e+00_fp_kind,&
     6.653500e+00_fp_kind,5.353300e+00_fp_kind,4.568300e+00_fp_kind,3.864000e+00_fp_kind,3.209100e+00_fp_kind,&
     2.207000e+00_fp_kind,1.635100e+00_fp_kind], 1.2975582443695399e+00_fp_kind,5.221100e+00_fp_kind) /
data ionicXrayFF( 1, 48) / t_ionicFF('Cd',  48,  0, &
     [4.800000e+01_fp_kind,4.708510e+01_fp_kind,4.479860e+01_fp_kind,4.192600e+01_fp_kind,3.893400e+01_fp_kind,&
     3.601220e+01_fp_kind,3.325620e+01_fp_kind,3.073060e+01_fp_kind,2.847350e+01_fp_kind,2.649640e+01_fp_kind,&
     2.478860e+01_fp_kind,2.206750e+01_fp_kind,2.003140e+01_fp_kind,1.840940e+01_fp_kind,1.700430e+01_fp_kind,&
     1.570250e+01_fp_kind,1.325700e+01_fp_kind,1.106370e+01_fp_kind,9.251400e+00_fp_kind,7.869500e+00_fp_kind,&
     6.873300e+00_fp_kind,5.462900e+00_fp_kind,4.667300e+00_fp_kind,3.979700e+00_fp_kind,3.329400e+00_fp_kind,&
     2.308400e+00_fp_kind,1.696100e+00_fp_kind], 1.3186707815108389e+00_fp_kind,9.202400e+00_fp_kind) /
data ionicXrayFF( 2, 48) / t_ionicFF('Cd',  48,  2, &
     [4.600000e+01_fp_kind,4.547120e+01_fp_kind,4.397180e+01_fp_kind,4.173150e+01_fp_kind,3.904550e+01_fp_kind,&
     3.619680e+01_fp_kind,3.340610e+01_fp_kind,3.081690e+01_fp_kind,2.850350e+01_fp_kind,2.648840e+01_fp_kind,&
     2.476040e+01_fp_kind,2.203440e+01_fp_kind,2.001270e+01_fp_kind,1.840520e+01_fp_kind,1.700890e+01_fp_kind,&
     1.571040e+01_fp_kind,1.326260e+01_fp_kind,1.106510e+01_fp_kind,9.250900e+00_fp_kind,7.869000e+00_fp_kind,&
     6.873200e+00_fp_kind,5.463100e+00_fp_kind,4.666400e+00_fp_kind,3.977700e+00_fp_kind,3.329800e+00_fp_kind,&
     2.304900e+00_fp_kind,1.698000e+00_fp_kind], 1.3202021510949820e+00_fp_kind,5.132400e+00_fp_kind) /
data ionicXrayFF( 1, 49) / t_ionicFF('In',  49,  0, &
     [4.900000e+01_fp_kind,4.796470e+01_fp_kind,4.550960e+01_fp_kind,4.259110e+01_fp_kind,3.964170e+01_fp_kind,&
     3.678350e+01_fp_kind,3.406950e+01_fp_kind,3.154670e+01_fp_kind,2.925420e+01_fp_kind,2.721400e+01_fp_kind,&
     2.542870e+01_fp_kind,2.255450e+01_fp_kind,2.041090e+01_fp_kind,1.874010e+01_fp_kind,1.733370e+01_fp_kind,&
     1.605770e+01_fp_kind,1.367470e+01_fp_kind,1.149590e+01_fp_kind,9.641300e+00_fp_kind,8.185700e+00_fp_kind,&
     7.112000e+00_fp_kind,5.579200e+00_fp_kind,4.762800e+00_fp_kind,4.090000e+00_fp_kind,3.448200e+00_fp_kind,&
     2.410700e+00_fp_kind,1.762300e+00_fp_kind], 1.3430543824106640e+00_fp_kind,1.059280e+01_fp_kind) /
data ionicXrayFF( 2, 49) / t_ionicFF('In',  49,  0, &
     [4.900000e+01_fp_kind,4.796470e+01_fp_kind,4.550960e+01_fp_kind,4.259110e+01_fp_kind,3.964170e+01_fp_kind,&
     3.678350e+01_fp_kind,3.406950e+01_fp_kind,3.154670e+01_fp_kind,2.925420e+01_fp_kind,2.721400e+01_fp_kind,&
     2.542870e+01_fp_kind,2.255450e+01_fp_kind,2.041090e+01_fp_kind,1.874010e+01_fp_kind,1.733370e+01_fp_kind,&
     1.605770e+01_fp_kind,1.367470e+01_fp_kind,1.149590e+01_fp_kind,9.641300e+00_fp_kind,8.185700e+00_fp_kind,&
     7.112000e+00_fp_kind,5.579200e+00_fp_kind,4.762800e+00_fp_kind,4.090000e+00_fp_kind,3.448200e+00_fp_kind,&
     2.410700e+00_fp_kind,1.762300e+00_fp_kind], 1.3430543824106640e+00_fp_kind,1.059280e+01_fp_kind) /
data ionicXrayFF( 1, 50) / t_ionicFF('Sn',  50,  0, &
     [5.000000e+01_fp_kind,4.891560e+01_fp_kind,4.632210e+01_fp_kind,4.327600e+01_fp_kind,4.029070e+01_fp_kind,&
     3.746840e+01_fp_kind,3.480820e+01_fp_kind,3.231860e+01_fp_kind,3.002360e+01_fp_kind,2.794720e+01_fp_kind,&
     2.610190e+01_fp_kind,2.308240e+01_fp_kind,2.081630e+01_fp_kind,1.907560e+01_fp_kind,1.764990e+01_fp_kind,&
     1.638880e+01_fp_kind,1.406730e+01_fp_kind,1.191800e+01_fp_kind,1.003720e+01_fp_kind,8.517800e+00_fp_kind,&
     7.369300e+00_fp_kind,5.703700e+00_fp_kind,4.855800e+00_fp_kind,4.194800e+00_fp_kind,3.564300e+00_fp_kind,&
     2.515100e+00_fp_kind,1.832100e+00_fp_kind], 1.3692853539390410e+00_fp_kind,1.103760e+01_fp_kind) /
data ionicXrayFF( 2, 50) / t_ionicFF('Sn',  50,  0, &
     [5.000000e+01_fp_kind,4.891560e+01_fp_kind,4.632210e+01_fp_kind,4.327600e+01_fp_kind,4.029070e+01_fp_kind,&
     3.746840e+01_fp_kind,3.480820e+01_fp_kind,3.231860e+01_fp_kind,3.002360e+01_fp_kind,2.794720e+01_fp_kind,&
     2.610190e+01_fp_kind,2.308240e+01_fp_kind,2.081630e+01_fp_kind,1.907560e+01_fp_kind,1.764990e+01_fp_kind,&
     1.638880e+01_fp_kind,1.406730e+01_fp_kind,1.191800e+01_fp_kind,1.003720e+01_fp_kind,8.517800e+00_fp_kind,&
     7.369300e+00_fp_kind,5.703700e+00_fp_kind,4.855800e+00_fp_kind,4.194800e+00_fp_kind,3.564300e+00_fp_kind,&
     2.515100e+00_fp_kind,1.832100e+00_fp_kind], 1.3692853539390410e+00_fp_kind,1.103760e+01_fp_kind) /
data ionicXrayFF( 1, 51) / t_ionicFF('Sb',  51,  0, &
     [5.100000e+01_fp_kind,4.989460e+01_fp_kind,4.720230e+01_fp_kind,4.401130e+01_fp_kind,4.093270e+01_fp_kind,&
     3.809950e+01_fp_kind,3.547770e+01_fp_kind,3.303240e+01_fp_kind,3.075990e+01_fp_kind,2.867440e+01_fp_kind,&
     2.679130e+01_fp_kind,2.364750e+01_fp_kind,2.125350e+01_fp_kind,1.942490e+01_fp_kind,1.796100e+01_fp_kind,&
     1.670130e+01_fp_kind,1.443450e+01_fp_kind,1.232590e+01_fp_kind,1.043460e+01_fp_kind,8.862800e+00_fp_kind,&
     7.644300e+00_fp_kind,5.838300e+00_fp_kind,4.947600e+00_fp_kind,4.294500e+00_fp_kind,3.677200e+00_fp_kind,&
     2.620700e+00_fp_kind,1.905600e+00_fp_kind], 1.3973406335549541e+00_fp_kind,1.117720e+01_fp_kind) /
data ionicXrayFF( 2, 51) / t_ionicFF('Sb',  51,  0, &
     [5.100000e+01_fp_kind,4.989460e+01_fp_kind,4.720230e+01_fp_kind,4.401130e+01_fp_kind,4.093270e+01_fp_kind,&
     3.809950e+01_fp_kind,3.547770e+01_fp_kind,3.303240e+01_fp_kind,3.075990e+01_fp_kind,2.867440e+01_fp_kind,&
     2.679130e+01_fp_kind,2.364750e+01_fp_kind,2.125350e+01_fp_kind,1.942490e+01_fp_kind,1.796100e+01_fp_kind,&
     1.670130e+01_fp_kind,1.443450e+01_fp_kind,1.232590e+01_fp_kind,1.043460e+01_fp_kind,8.862800e+00_fp_kind,&
     7.644300e+00_fp_kind,5.838300e+00_fp_kind,4.947600e+00_fp_kind,4.294500e+00_fp_kind,3.677200e+00_fp_kind,&
     2.620700e+00_fp_kind,1.905600e+00_fp_kind], 1.3973406335549541e+00_fp_kind,1.117720e+01_fp_kind) /
data ionicXrayFF( 1, 52) / t_ionicFF('Te',  52,  0, &
     [5.200000e+01_fp_kind,5.088830e+01_fp_kind,4.812890e+01_fp_kind,4.480010e+01_fp_kind,4.159560e+01_fp_kind,&
     3.870680e+01_fp_kind,3.609450e+01_fp_kind,3.368810e+01_fp_kind,3.145090e+01_fp_kind,2.937800e+01_fp_kind,&
     2.747980e+01_fp_kind,2.424180e+01_fp_kind,2.172430e+01_fp_kind,1.979530e+01_fp_kind,1.827540e+01_fp_kind,&
     1.700160e+01_fp_kind,1.477680e+01_fp_kind,1.271660e+01_fp_kind,1.082940e+01_fp_kind,9.217600e+00_fp_kind,&
     7.935400e+00_fp_kind,5.984500e+00_fp_kind,5.039800e+00_fp_kind,4.389600e+00_fp_kind,3.786200e+00_fp_kind,&
     2.727100e+00_fp_kind,1.982700e+00_fp_kind], 1.4270502406167509e+00_fp_kind,1.117470e+01_fp_kind) /
data ionicXrayFF( 2, 52) / t_ionicFF('Te',  52,  0, &
     [5.200000e+01_fp_kind,5.088830e+01_fp_kind,4.812890e+01_fp_kind,4.480010e+01_fp_kind,4.159560e+01_fp_kind,&
     3.870680e+01_fp_kind,3.609450e+01_fp_kind,3.368810e+01_fp_kind,3.145090e+01_fp_kind,2.937800e+01_fp_kind,&
     2.747980e+01_fp_kind,2.424180e+01_fp_kind,2.172430e+01_fp_kind,1.979530e+01_fp_kind,1.827540e+01_fp_kind,&
     1.700160e+01_fp_kind,1.477680e+01_fp_kind,1.271660e+01_fp_kind,1.082940e+01_fp_kind,9.217600e+00_fp_kind,&
     7.935400e+00_fp_kind,5.984500e+00_fp_kind,5.039800e+00_fp_kind,4.389600e+00_fp_kind,3.786200e+00_fp_kind,&
     2.727100e+00_fp_kind,1.982700e+00_fp_kind], 1.4270502406167509e+00_fp_kind,1.117470e+01_fp_kind) /
data ionicXrayFF( 1, 53) / t_ionicFF('I ',  53,  0, &
     [5.300000e+01_fp_kind,5.190100e+01_fp_kind,4.911720e+01_fp_kind,4.567300e+01_fp_kind,4.231980e+01_fp_kind,&
     3.932660e+01_fp_kind,3.667950e+01_fp_kind,3.428930e+01_fp_kind,3.208680e+01_fp_kind,3.004150e+01_fp_kind,&
     2.814970e+01_fp_kind,2.485590e+01_fp_kind,2.222980e+01_fp_kind,2.019440e+01_fp_kind,1.860210e+01_fp_kind,&
     1.729730e+01_fp_kind,1.509580e+01_fp_kind,1.308710e+01_fp_kind,1.121760e+01_fp_kind,9.579000e+00_fp_kind,&
     8.241200e+00_fp_kind,6.143800e+00_fp_kind,5.133900e+00_fp_kind,4.480500e+00_fp_kind,3.891100e+00_fp_kind,&
     2.833700e+00_fp_kind,2.063500e+00_fp_kind], 1.4584040913686560e+00_fp_kind,1.098110e+01_fp_kind) /
data ionicXrayFF( 2, 53) / t_ionicFF('I ',  53, -1, &
     [5.400000e+01_fp_kind,5.269170e+01_fp_kind,4.949040e+01_fp_kind,4.573840e+01_fp_kind,4.226830e+01_fp_kind,&
     3.926890e+01_fp_kind,3.664610e+01_fp_kind,3.427660e+01_fp_kind,3.208580e+01_fp_kind,3.004590e+01_fp_kind,&
     2.815600e+01_fp_kind,2.486110e+01_fp_kind,2.223260e+01_fp_kind,2.019540e+01_fp_kind,1.860200e+01_fp_kind,&
     1.729670e+01_fp_kind,1.509520e+01_fp_kind,1.308700e+01_fp_kind,1.121780e+01_fp_kind,9.579400e+00_fp_kind,&
     8.241500e+00_fp_kind,6.144000e+00_fp_kind,5.134000e+00_fp_kind,4.481100e+00_fp_kind,3.890100e+00_fp_kind,&
     2.825700e+00_fp_kind,2.062600e+00_fp_kind], 1.4577422483283371e+00_fp_kind,1.319510e+01_fp_kind) /
data ionicXrayFF( 1, 54) / t_ionicFF('Xe',  54,  0, &
     [5.400000e+01_fp_kind,5.291740e+01_fp_kind,5.012740e+01_fp_kind,4.659160e+01_fp_kind,4.306230e+01_fp_kind,&
     3.997130e+01_fp_kind,3.725590e+01_fp_kind,3.685530e+01_fp_kind,3.267630e+01_fp_kind,3.066210e+01_fp_kind,&
     2.879020e+01_fp_kind,2.547530e+01_fp_kind,2.276240e+01_fp_kind,2.062270e+01_fp_kind,1.894710e+01_fp_kind,&
     1.759580e+01_fp_kind,1.539540e+01_fp_kind,1.343650e+01_fp_kind,1.159540e+01_fp_kind,9.943100e+00_fp_kind,&
     8.558900e+00_fp_kind,6.316900e+00_fp_kind,5.231300e+00_fp_kind,4.568200e+00_fp_kind,3.991600e+00_fp_kind,&
     2.940200e+00_fp_kind,2.147500e+00_fp_kind], 1.4909599344293909e+00_fp_kind,1.076520e+01_fp_kind) /
data ionicXrayFF( 2, 54) / t_ionicFF('Xe',  54,  0, &
     [5.400000e+01_fp_kind,5.291740e+01_fp_kind,5.012740e+01_fp_kind,4.659160e+01_fp_kind,4.306230e+01_fp_kind,&
     3.997130e+01_fp_kind,3.725590e+01_fp_kind,3.685530e+01_fp_kind,3.267630e+01_fp_kind,3.066210e+01_fp_kind,&
     2.879020e+01_fp_kind,2.547530e+01_fp_kind,2.276240e+01_fp_kind,2.062270e+01_fp_kind,1.894710e+01_fp_kind,&
     1.759580e+01_fp_kind,1.539540e+01_fp_kind,1.343650e+01_fp_kind,1.159540e+01_fp_kind,9.943100e+00_fp_kind,&
     8.558900e+00_fp_kind,6.316900e+00_fp_kind,5.231300e+00_fp_kind,4.568200e+00_fp_kind,3.991600e+00_fp_kind,&
     2.940200e+00_fp_kind,2.147500e+00_fp_kind], 1.4909599344293909e+00_fp_kind,1.076520e+01_fp_kind) /
data ionicXrayFF( 1, 55) / t_ionicFF('Cs',  55,  0, &
     [5.500000e+01_fp_kind,5.352750e+01_fp_kind,5.060450e+01_fp_kind,4.729390e+01_fp_kind,4.389250e+01_fp_kind,&
     4.071720e+01_fp_kind,3.790860e+01_fp_kind,3.544500e+01_fp_kind,3.324680e+01_fp_kind,3.124210e+01_fp_kind,&
     2.938790e+01_fp_kind,2.607760e+01_fp_kind,2.330740e+01_fp_kind,2.107640e+01_fp_kind,1.931440e+01_fp_kind,&
     1.790460e+01_fp_kind,1.568110e+01_fp_kind,1.376470e+01_fp_kind,1.195950e+01_fp_kind,1.030530e+01_fp_kind,&
     8.885000e+00_fp_kind,6.505100e+00_fp_kind,5.334300e+00_fp_kind,4.654600e+00_fp_kind,4.086100e+00_fp_kind,&
     3.045100e+00_fp_kind,2.234600e+00_fp_kind], 1.5245899775231491e+00_fp_kind,1.634830e+01_fp_kind) /
data ionicXrayFF( 2, 55) / t_ionicFF('Cs',  55,  1, &
     [5.400000e+01_fp_kind,5.308450e+01_fp_kind,5.063670e+01_fp_kind,4.734850e+01_fp_kind,4.390140e+01_fp_kind,&
     4.070760e+01_fp_kind,3.789820e+01_fp_kind,3.543920e+01_fp_kind,3.324550e+01_fp_kind,3.124360e+01_fp_kind,&
     2.939040e+01_fp_kind,2.607950e+01_fp_kind,2.330770e+01_fp_kind,2.107580e+01_fp_kind,1.931360e+01_fp_kind,&
     1.790420e+01_fp_kind,1.568120e+01_fp_kind,1.376500e+01_fp_kind,1.195970e+01_fp_kind,1.030570e+01_fp_kind,&
     8.885100e+00_fp_kind,6.504200e+00_fp_kind,5.333600e+00_fp_kind,4.653100e+00_fp_kind,4.088000e+00_fp_kind,&
     3.044600e+00_fp_kind,2.235000e+00_fp_kind], 1.5248744432862720e+00_fp_kind,9.015400e+00_fp_kind) /
data ionicXrayFF( 1, 56) / t_ionicFF('Ba',  56,  0, &
     [5.600000e+01_fp_kind,5.434490e+01_fp_kind,5.112400e+01_fp_kind,4.784170e+01_fp_kind,4.459010e+01_fp_kind,&
     4.146040e+01_fp_kind,3.860260e+01_fp_kind,3.606860e+01_fp_kind,3.382340e+01_fp_kind,3.180340e+01_fp_kind,&
     2.995430e+01_fp_kind,2.665790e+01_fp_kind,2.385600e+01_fp_kind,2.155110e+01_fp_kind,1.970510e+01_fp_kind,&
     1.822860e+01_fp_kind,1.595800e+01_fp_kind,1.407300e+01_fp_kind,1.230780e+01_fp_kind,1.066320e+01_fp_kind,&
     9.217200e+00_fp_kind,6.706800e+00_fp_kind,5.443400e+00_fp_kind,4.739800e+00_fp_kind,4.175500e+00_fp_kind,&
     3.151400e+00_fp_kind,2.323200e+00_fp_kind], 1.5581256706808091e+00_fp_kind,1.811230e+01_fp_kind) /
data ionicXrayFF( 2, 56) / t_ionicFF('Ba',  56,  2, &
     [5.400000e+01_fp_kind,5.320340e+01_fp_kind,5.102550e+01_fp_kind,4.798510e+01_fp_kind,4.465440e+01_fp_kind,&
     4.145070e+01_fp_kind,3.857160e+01_fp_kind,3.604470e+01_fp_kind,3.381340e+01_fp_kind,3.180450e+01_fp_kind,&
     2.996110e+01_fp_kind,2.666480e+01_fp_kind,2.385790e+01_fp_kind,2.154940e+01_fp_kind,1.970240e+01_fp_kind,&
     1.822640e+01_fp_kind,1.595790e+01_fp_kind,1.407380e+01_fp_kind,1.230840e+01_fp_kind,1.066350e+01_fp_kind,&
     9.217200e+00_fp_kind,6.705800e+00_fp_kind,5.442400e+00_fp_kind,4.736800e+00_fp_kind,4.179600e+00_fp_kind,&
     3.148600e+00_fp_kind,2.324400e+00_fp_kind], 1.5589653399309920e+00_fp_kind,7.800000e+00_fp_kind) /
data ionicXrayFF( 1, 57) / t_ionicFF('La',  57,  0, &
     [5.700000e+01_fp_kind,5.534910e+01_fp_kind,5.197490e+01_fp_kind,4.851590e+01_fp_kind,4.520970e+01_fp_kind,&
     4.208110e+01_fp_kind,3.921840e+01_fp_kind,3.666660e+01_fp_kind,3.440420e+01_fp_kind,3.237640e+01_fp_kind,&
     3.052840e+01_fp_kind,2.723570e+01_fp_kind,2.440550e+01_fp_kind,2.203550e+01_fp_kind,2.011030e+01_fp_kind,&
     1.856470e+01_fp_kind,1.623090e+01_fp_kind,1.436550e+01_fp_kind,1.264070e+01_fp_kind,1.101410e+01_fp_kind,&
     9.552400e+00_fp_kind,6.922400e+00_fp_kind,5.559800e+00_fp_kind,4.823000e+00_fp_kind,4.262100e+00_fp_kind,&
     3.255800e+00_fp_kind,2.412600e+00_fp_kind], 1.5910924499060140e+00_fp_kind,1.770540e+01_fp_kind) /
data ionicXrayFF( 2, 57) / t_ionicFF('La',  57,  3, &
     [5.400000e+01_fp_kind,5.294300e+00_fp_kind,5.133520e+01_fp_kind,4.852460e+01_fp_kind,4.533970e+01_fp_kind,&
     4.217380e+01_fp_kind,3.925810e+01_fp_kind,3.667040e+01_fp_kind,3.439020e+01_fp_kind,3.235730e+01_fp_kind,&
     3.051120e+01_fp_kind,2.722820e+01_fp_kind,2.440490e+01_fp_kind,2.203730e+01_fp_kind,2.011200e+01_fp_kind,&
     1.856570e+01_fp_kind,1.623070e+01_fp_kind,1.436500e+01_fp_kind,1.264020e+01_fp_kind,1.101370e+01_fp_kind,&
     9.552000e+00_fp_kind,6.921500e+00_fp_kind,5.558800e+00_fp_kind,4.820200e+00_fp_kind,4.267100e+00_fp_kind,&
     3.251500e+00_fp_kind,2.415400e+00_fp_kind], 1.5930207421140781e+00_fp_kind,6.884100e+00_fp_kind) /
data ionicXrayFF( 1, 58) / t_ionicFF('Ce',  58,  0, &
     [5.800000e+01_fp_kind,5.638170e+01_fp_kind,5.303720e+01_fp_kind,4.956820e+01_fp_kind,4.622770e+01_fp_kind,&
     4.304360e+01_fp_kind,4.010910e+01_fp_kind,3.747960e+01_fp_kind,3.514440e+01_fp_kind,3.305500e+01_fp_kind,&
     3.115800e+01_fp_kind,2.779390e+01_fp_kind,2.490440e+01_fp_kind,2.247260e+01_fp_kind,2.048490e+01_fp_kind,&
     1.888540e+01_fp_kind,1.649470e+01_fp_kind,1.463670e+01_fp_kind,1.294380e+01_fp_kind,1.133750e+01_fp_kind,&
     9.869900e+00_fp_kind,7.144100e+00_fp_kind,5.685100e+00_fp_kind,4.908500e+00_fp_kind,4.346600e+00_fp_kind,&
     3.355800e+00_fp_kind,2.503200e+00_fp_kind], 1.6237909212783350e+00_fp_kind,1.729400e+01_fp_kind) /
data ionicXrayFF( 2, 58) / t_ionicFF('Ce',  58,  4, &
     [5.400000e+01_fp_kind,5.336680e+01_fp_kind,5.158940e+01_fp_kind,4.898680e+01_fp_kind,4.595820e+01_fp_kind,&
     4.286230e+01_fp_kind,3.994150e+01_fp_kind,3.730900e+01_fp_kind,3.497880e+01_fp_kind,3.291090e+01_fp_kind,&
     3.104960e+01_fp_kind,2.777050e+01_fp_kind,2.494250e+01_fp_kind,2.253300e+01_fp_kind,2.053960e+01_fp_kind,&
     1.892280e+01_fp_kind,1.650440e+01_fp_kind,1.464170e+01_fp_kind,1.295480e+01_fp_kind,1.135370e+01_fp_kind,&
     9.886500e+00_fp_kind,7.150500e+00_fp_kind,5.684000e+00_fp_kind,4.904200e+00_fp_kind,4.351100e+00_fp_kind,&
     3.352700e+00_fp_kind,2.507600e+00_fp_kind], 1.6267741168159939e+00_fp_kind,6.159800e+00_fp_kind) /
data ionicXrayFF( 1, 59) / t_ionicFF('Pr',  59,  0, &
     [5.900000e+01_fp_kind,5.743820e+01_fp_kind,5.427700e+01_fp_kind,5.094850e+01_fp_kind,4.759830e+01_fp_kind,&
     4.430980e+01_fp_kind,4.124240e+01_fp_kind,3.848330e+01_fp_kind,3.603140e+01_fp_kind,3.384030e+01_fp_kind,&
     3.185620e+01_fp_kind,2.835700e+01_fp_kind,2.537280e+01_fp_kind,2.287120e+01_fp_kind,2.082890e+01_fp_kind,&
     1.918690e+01_fp_kind,1.674990e+01_fp_kind,1.489260e+01_fp_kind,1.322360e+01_fp_kind,1.163700e+01_fp_kind,&
     1.016960e+01_fp_kind,7.368400e+00_fp_kind,5.818300e+00_fp_kind,4.997100e+00_fp_kind,4.429900e+00_fp_kind,&
     3.451900e+00_fp_kind,2.594800e+00_fp_kind], 1.6561026288356400e+00_fp_kind,1.686310e+01_fp_kind) /
data ionicXrayFF( 2, 59) / t_ionicFF('Pr',  59,  3, &
     [5.600000e+01_fp_kind,5.531300e+01_fp_kind,5.339360e+01_fp_kind,5.060580e+01_fp_kind,4.739360e+01_fp_kind,&
     4.414010e+01_fp_kind,4.109030e+01_fp_kind,3.834830e+01_fp_kind,3.591790e+01_fp_kind,3.375410e+01_fp_kind,&
     3.180120e+01_fp_kind,2.836060e+01_fp_kind,2.541130e+01_fp_kind,2.291800e+01_fp_kind,2.086740e+01_fp_kind,&
     1.921160e+01_fp_kind,1.675490e+01_fp_kind,1.489490e+01_fp_kind,1.323080e+01_fp_kind,1.164840e+01_fp_kind,&
     1.018180e+01_fp_kind,7.373500e+00_fp_kind,5.817700e+00_fp_kind,4.993300e+00_fp_kind,4.433600e+00_fp_kind,&
     3.449300e+00_fp_kind,2.597700e+00_fp_kind], 1.6580387679225841e+00_fp_kind,6.690600e+00_fp_kind) /
data ionicXrayFF( 1, 60) / t_ionicFF('Nd',  60,  0, &
     [6.000000e+01_fp_kind,5.846780e+01_fp_kind,5.533610e+01_fp_kind,5.200940e+01_fp_kind,4.864230e+01_fp_kind,&
     4.531500e+01_fp_kind,4.218820e+01_fp_kind,3.935840e+01_fp_kind,3.683540e+01_fp_kind,3.457970e+01_fp_kind,&
     3.254040e+01_fp_kind,2.895470e+01_fp_kind,2.589920e+01_fp_kind,2.332890e+01_fp_kind,2.121960e+01_fp_kind,&
     1.951830e+01_fp_kind,1.700830e+01_fp_kind,1.514220e+01_fp_kind,1.349900e+01_fp_kind,1.193810e+01_fp_kind,&
     1.047690e+01_fp_kind,7.606100e+00_fp_kind,5.960000e+00_fp_kind,5.087400e+00_fp_kind,4.510500e+00_fp_kind,&
     3.547300e+00_fp_kind,2.687000e+00_fp_kind], 1.6877846212901060e+00_fp_kind,1.648970e+01_fp_kind) /
data ionicXrayFF( 2, 60) / t_ionicFF('Nd',  60,  3, &
     [5.700000e+01_fp_kind,5.632290e+01_fp_kind,5.442570e+01_fp_kind,5.165490e+01_fp_kind,4.843810e+01_fp_kind,&
     4.515150e+01_fp_kind,4.204440e+01_fp_kind,3.923120e+01_fp_kind,3.672740e+01_fp_kind,3.449590e+01_fp_kind,&
     3.248460e+01_fp_kind,2.895360e+01_fp_kind,2.593340e+01_fp_kind,2.337370e+01_fp_kind,2.125830e+01_fp_kind,&
     1.954440e+01_fp_kind,1.701420e+01_fp_kind,1.514380e+01_fp_kind,1.350460e+01_fp_kind,1.194810e+01_fp_kind,&
     1.048840e+01_fp_kind,7.611700e+00_fp_kind,5.959900e+00_fp_kind,5.083800e+00_fp_kind,4.514000e+00_fp_kind,&
     3.544600e+00_fp_kind,2.689900e+00_fp_kind], 1.6896916948321490e+00_fp_kind,6.589400e+00_fp_kind) /
data ionicXrayFF( 1, 61) / t_ionicFF('Pm',  61,  0, &
     [6.100000e+01_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind], 0.0000000000000000e+00_fp_kind,0.000000e+00_fp_kind) /
data ionicXrayFF( 2, 61) / t_ionicFF('Pm',  61,  0, &
     [6.100000e+01_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,0.000000e+00_fp_kind,&
     0.000000e+00_fp_kind,0.000000e+00_fp_kind], 0.0000000000000000e+00_fp_kind,0.000000e+00_fp_kind) /
data ionicXrayFF( 1, 62) / t_ionicFF('Sm',  62,  0, &
     [6.200000e+01_fp_kind,6.052380e+01_fp_kind,5.745420e+01_fp_kind,5.414150e+01_fp_kind,5.075540e+01_fp_kind,&
     4.736810e+01_fp_kind,4.413920e+01_fp_kind,4.118010e+01_fp_kind,3.851990e+01_fp_kind,3.613410e+01_fp_kind,&
     3.397880e+01_fp_kind,3.020450e+01_fp_kind,2.699500e+01_fp_kind,2.428120e+01_fp_kind,2.203360e+01_fp_kind,&
     2.020760e+01_fp_kind,1.752850e+01_fp_kind,1.561810e+01_fp_kind,1.401140e+01_fp_kind,1.250380e+01_fp_kind,&
     1.106930e+01_fp_kind,8.099600e+00_fp_kind,6.269200e+00_fp_kind,5.279700e+00_fp_kind,4.669000e+00_fp_kind,&
     3.730300e+00_fp_kind,2.871600e+00_fp_kind], 1.7483578111364440e+00_fp_kind,1.579490e+01_fp_kind) /
data ionicXrayFF( 2, 62) / t_ionicFF('Sm',  62,  3, &
     [5.900000e+01_fp_kind,5.834280e+01_fp_kind,5.649160e+01_fp_kind,5.376120e+01_fp_kind,5.054760e+01_fp_kind,&
     4.721170e+01_fp_kind,4.400700e+01_fp_kind,4.106500e+01_fp_kind,3.842150e+01_fp_kind,3.605530e+01_fp_kind,&
     3.392300e+01_fp_kind,3.019670e+01_fp_kind,2.702200e+01_fp_kind,2.432240e+01_fp_kind,2.207270e+01_fp_kind,&
     2.023640e+01_fp_kind,1.753670e+01_fp_kind,1.561890e+01_fp_kind,1.401450e+01_fp_kind,1.251130e+01_fp_kind,&
     1.107910e+01_fp_kind,8.106000e+00_fp_kind,6.270000e+00_fp_kind,5.276500e+00_fp_kind,4.672300e+00_fp_kind,&
     3.727300e+00_fp_kind,2.874800e+00_fp_kind], 1.7504008443100361e+00_fp_kind,6.388100e+00_fp_kind) /
data ionicXrayFF( 1, 63) / t_ionicFF('Eu',  63,  0, &
     [6.300000e+01_fp_kind,6.155020e+01_fp_kind,5.851600e+01_fp_kind,5.520640e+01_fp_kind,5.181520e+01_fp_kind,&
     4.840510e+01_fp_kind,4.513320e+01_fp_kind,4.211640e+01_fp_kind,3.939220e+01_fp_kind,3.694320e+01_fp_kind,&
     3.472990e+01_fp_kind,3.085790e+01_fp_kind,2.756790e+01_fp_kind,2.478000e+01_fp_kind,2.246120e+01_fp_kind,&
     2.056950e+01_fp_kind,1.779450e+01_fp_kind,1.584830e+01_fp_kind,1.425100e+01_fp_kind,1.276890e+01_fp_kind,&
     1.135270e+01_fp_kind,8.352800e+00_fp_kind,6.436400e+00_fp_kind,5.382900e+00_fp_kind,4.748300e+00_fp_kind,&
     3.817700e+00_fp_kind,2.964000e+00_fp_kind], 1.7773335998401050e+00_fp_kind,1.547000e+01_fp_kind) /
data ionicXrayFF( 2, 63) / t_ionicFF('Eu',  63,  3, &
     [6.000000e+01_fp_kind,4.935250e+01_fp_kind,5.752460e+01_fp_kind,5.481640e+01_fp_kind,5.160920e+01_fp_kind,&
     4.825570e+01_fp_kind,4.500990e+01_fp_kind,4.201020e+01_fp_kind,3.930100e+01_fp_kind,3.686900e+01_fp_kind,&
     3.467570e+01_fp_kind,3.084680e+01_fp_kind,2.758990e+01_fp_kind,2.481700e+01_fp_kind,2.249780e+01_fp_kind,&
     2.059750e+01_fp_kind,1.780310e+01_fp_kind,1.584890e+01_fp_kind,1.425300e+01_fp_kind,1.277500e+01_fp_kind,&
     1.136150e+01_fp_kind,8.359300e+00_fp_kind,6.437600e+00_fp_kind,5.380000e+00_fp_kind,4.751500e+00_fp_kind,&
     3.814600e+00_fp_kind,2.966800e+00_fp_kind], 1.7790955671195330e+00_fp_kind,6.289700e+00_fp_kind) /
data ionicXrayFF( 1, 64) / t_ionicFF('Gd',  64,  0, &
     [6.400000e+01_fp_kind,6.254910e+01_fp_kind,5.940010e+01_fp_kind,5.596690e+01_fp_kind,5.255150e+01_fp_kind,&
     4.917850e+01_fp_kind,4.593880e+01_fp_kind,4.292530e+01_fp_kind,4.018050e+01_fp_kind,3.769980e+01_fp_kind,&
     3.545310e+01_fp_kind,3.152290e+01_fp_kind,2.817960e+01_fp_kind,2.533030e+01_fp_kind,2.293970e+01_fp_kind,&
     2.097330e+01_fp_kind,1.807700e+01_fp_kind,1.607720e+01_fp_kind,1.448310e+01_fp_kind,1.302840e+01_fp_kind,&
     1.163580e+01_fp_kind,8.616600e+00_fp_kind,6.614500e+00_fp_kind,5.491400e+00_fp_kind,4.827400e+00_fp_kind,&
     3.903800e+00_fp_kind,3.055400e+00_fp_kind], 1.8048260223219119e+00_fp_kind,1.525910e+01_fp_kind) /
data ionicXrayFF( 2, 64) / t_ionicFF('Gd',  64,  3, &
     [6.100000e+01_fp_kind,6.036200e+01_fp_kind,5.855680e+01_fp_kind,5.587110e+01_fp_kind,5.267230e+01_fp_kind,&
     4.930490e+01_fp_kind,4.602260e+01_fp_kind,4.296940e+01_fp_kind,4.019810e+01_fp_kind,3.770260e+01_fp_kind,&
     3.544910e+01_fp_kind,3.151710e+01_fp_kind,2.817620e+01_fp_kind,2.532850e+01_fp_kind,2.293840e+01_fp_kind,&
     2.097200e+01_fp_kind,1.807590e+01_fp_kind,1.607630e+01_fp_kind,1.448230e+01_fp_kind,1.302740e+01_fp_kind,&
     1.163460e+01_fp_kind,8.615100e+00_fp_kind,6.613200e+00_fp_kind,5.488900e+00_fp_kind,4.831700e+00_fp_kind,&
     3.899300e+00_fp_kind,3.058300e+00_fp_kind], 1.8066250203062990e+00_fp_kind,6.195400e+00_fp_kind) /
data ionicXrayFF( 1, 65) / t_ionicFF('Tb',  65,  0, &
     [6.500000e+01_fp_kind,6.360140e+01_fp_kind,6.062780e+01_fp_kind,5.735340e+01_fp_kind,5.396530e+01_fp_kind,&
     5.052310e+01_fp_kind,4.717830e+01_fp_kind,4.405630e+01_fp_kind,4.120970e+01_fp_kind,3.863530e+01_fp_kind,&
     3.630260e+01_fp_kind,3.222300e+01_fp_kind,2.875990e+01_fp_kind,2.581640e+01_fp_kind,2.335090e+01_fp_kind,&
     2.132370e+01_fp_kind,1.834140e+01_fp_kind,1.630040e+01_fp_kind,1.470290e+01_fp_kind,1.326480e+01_fp_kind,&
     1.189050e+01_fp_kind,8.864600e+00_fp_kind,6.793900e+00_fp_kind,5.606100e+00_fp_kind,4.910400e+00_fp_kind,&
     3.984800e+00_fp_kind,3.145900e+00_fp_kind], 1.8309602758750001e+00_fp_kind,1.485260e+01_fp_kind) /
data ionicXrayFF( 2, 65) / t_ionicFF('Tb',  65,  3, &
     [6.200000e+01_fp_kind,6.137150e+01_fp_kind,5.958960e+01_fp_kind,5.692810e+01_fp_kind,5.374060e+01_fp_kind,&
     5.036330e+01_fp_kind,4.704870e+01_fp_kind,4.394580e+01_fp_kind,4.111460e+01_fp_kind,3.855610e+01_fp_kind,&
     3.624160e+01_fp_kind,3.220300e+01_fp_kind,2.877410e+01_fp_kind,2.584920e+01_fp_kind,2.338720e+01_fp_kind,&
     2.135420e+01_fp_kind,1.835320e+01_fp_kind,1.630200e+01_fp_kind,1.470370e+01_fp_kind,1.326900e+01_fp_kind,&
     1.189790e+01_fp_kind,8.871700e+00_fp_kind,6.796000e+00_fp_kind,5.603700e+00_fp_kind,4.913500e+00_fp_kind,&
     3.981500e+00_fp_kind,3.148800e+00_fp_kind], 1.8327340455803589e+00_fp_kind,6.099400e+00_fp_kind) /
data ionicXrayFF( 1, 66) / t_ionicFF('Dy',  66,  0, &
     [6.600000e+01_fp_kind,6.462610e+01_fp_kind,6.168540e+01_fp_kind,5.842750e+01_fp_kind,5.504460e+01_fp_kind,&
     5.159220e+01_fp_kind,4.821790e+01_fp_kind,5.404970e+01_fp_kind,4.214690e+01_fp_kind,3.951270e+01_fp_kind,&
     3.712160e+01_fp_kind,3.293700e+01_fp_kind,2.938450e+01_fp_kind,2.636050e+01_fp_kind,2.381930e+01_fp_kind,&
     2.172160e+01_fp_kind,1.862650e+01_fp_kind,1.652610e+01_fp_kind,1.491810e+01_fp_kind,1.349700e+01_fp_kind,&
     1.214470e+01_fp_kind,9.121000e+00_fp_kind,6.983100e+00_fp_kind,5.726600e+00_fp_kind,4.994400e+00_fp_kind,&
     4.064800e+00_fp_kind,3.235300e+00_fp_kind], 1.8556736509534839e+00_fp_kind,1.455770e+01_fp_kind) /
data ionicXrayFF( 2, 66) / t_ionicFF('Dy',  66,  3, &
     [6.300000e+01_fp_kind,6.238100e+01_fp_kind,6.062270e+01_fp_kind,5.798670e+01_fp_kind,5.481350e+01_fp_kind,&
     5.143050e+01_fp_kind,4.808840e+01_fp_kind,4.494040e+01_fp_kind,4.205320e+01_fp_kind,3.943440e+01_fp_kind,&
     3.706050e+01_fp_kind,3.291470e+01_fp_kind,2.939580e+01_fp_kind,2.639130e+01_fp_kind,2.385500e+01_fp_kind,&
     2.175270e+01_fp_kind,1.863980e+01_fp_kind,1.652810e+01_fp_kind,1.491840e+01_fp_kind,1.350030e+01_fp_kind,&
     1.215120e+01_fp_kind,9.128100e+00_fp_kind,6.985700e+00_fp_kind,5.724700e+00_fp_kind,4.997500e+00_fp_kind,&
     4.061300e+00_fp_kind,3.238100e+00_fp_kind], 1.8573625081458689e+00_fp_kind,6.004100e+00_fp_kind) /
data ionicXrayFF( 1, 67) / t_ionicFF('Ho',  67,  0, &
     [6.700000e+01_fp_kind,6.565000e+01_fp_kind,6.274190e+01_fp_kind,5.950150e+01_fp_kind,5.612570e+01_fp_kind,&
     5.266590e+01_fp_kind,4.926580e+01_fp_kind,4.605520e+01_fp_kind,4.309930e+01_fp_kind,4.040760e+01_fp_kind,&
     3.795910e+01_fp_kind,3.366940e+01_fp_kind,3.002620e+01_fp_kind,2.692060e+01_fp_kind,2.430270e+01_fp_kind,&
     2.213330e+01_fp_kind,1.892070e+01_fp_kind,1.675340e+01_fp_kind,1.512790e+01_fp_kind,1.371990e+01_fp_kind,&
     1.238920e+01_fp_kind,9.376200e+00_fp_kind,7.178400e+00_fp_kind,5.853400e+00_fp_kind,5.081100e+00_fp_kind,&
     4.142500e+00_fp_kind,3.323400e+00_fp_kind], 1.8789068511823790e+00_fp_kind,1.427300e+01_fp_kind) /
data ionicXrayFF( 2, 67) / t_ionicFF('Ho',  67,  3, &
     [6.400000e+01_fp_kind,6.339040e+01_fp_kind,6.165540e+01_fp_kind,5.904530e+01_fp_kind,5.588790e+01_fp_kind,&
     5.250190e+01_fp_kind,4.913610e+01_fp_kind,4.594680e+01_fp_kind,4.300670e+01_fp_kind,4.032980e+01_fp_kind,&
     3.789770e+01_fp_kind,3.364480e+01_fp_kind,3.003460e+01_fp_kind,2.694920e+01_fp_kind,2.433760e+01_fp_kind,&
     2.216470e+01_fp_kind,1.893520e+01_fp_kind,1.675610e+01_fp_kind,1.512800e+01_fp_kind,1.372230e+01_fp_kind,&
     1.239480e+01_fp_kind,9.383200e+00_fp_kind,7.181400e+00_fp_kind,5.851700e+00_fp_kind,5.084200e+00_fp_kind,&
     4.139000e+00_fp_kind,3.326200e+00_fp_kind], 1.8805725431810250e+00_fp_kind,5.911100e+00_fp_kind) /
data ionicXrayFF( 1, 68) / t_ionicFF('Er',  68,  0, &
     [6.800000e+01_fp_kind,6.667330e+01_fp_kind,6.379750e+01_fp_kind,6.057520e+01_fp_kind,5.720800e+01_fp_kind,&
     5.374360e+01_fp_kind,5.032110e+01_fp_kind,4.707180e+01_fp_kind,4.406580e+01_fp_kind,4.131860e+01_fp_kind,&
     3.881400e+01_fp_kind,3.441970e+01_fp_kind,3.068490e+01_fp_kind,2.749670e+01_fp_kind,2.480230e+01_fp_kind,&
     2.255890e+01_fp_kind,1.922470e+01_fp_kind,1.698390e+01_fp_kind,1.533390e+01_fp_kind,1.393430e+01_fp_kind,&
     1.262420e+01_fp_kind,9.629100e+00_fp_kind,7.379100e+00_fp_kind,5.986600e+00_fp_kind,5.171100e+00_fp_kind,&
     4.218400e+00_fp_kind,3.410000e+00_fp_kind], 1.9006038086391039e+00_fp_kind,1.399780e+01_fp_kind) /
data ionicXrayFF( 2, 68) / t_ionicFF('Er',  68,  3, &
     [6.500000e+01_fp_kind,6.439940e+01_fp_kind,6.268750e+01_fp_kind,6.010340e+01_fp_kind,5.696310e+01_fp_kind,&
     5.357690e+01_fp_kind,5.019080e+01_fp_kind,4.696370e+01_fp_kind,4.397390e+01_fp_kind,4.124130e+01_fp_kind,&
     3.875220e+01_fp_kind,3.439300e+01_fp_kind,3.069060e+01_fp_kind,2.752310e+01_fp_kind,2.483520e+01_fp_kind,&
     2.259060e+01_fp_kind,1.924050e+01_fp_kind,1.698730e+01_fp_kind,1.533380e+01_fp_kind,1.393610e+01_fp_kind,&
     1.262890e+01_fp_kind,9.636000e+00_fp_kind,7.382400e+00_fp_kind,5.985100e+00_fp_kind,5.174200e+00_fp_kind,&
     4.214800e+00_fp_kind,3.412700e+00_fp_kind], 1.9021882010859701e+00_fp_kind,5.820400e+00_fp_kind) /
data ionicXrayFF( 1, 69) / t_ionicFF('Tm',  69,  0, &
     [6.900000e+01_fp_kind,6.769590e+01_fp_kind,6.485210e+01_fp_kind,6.164850e+01_fp_kind,5.829120e+01_fp_kind,&
     5.482450e+01_fp_kind,5.138290e+01_fp_kind,4.809810e+01_fp_kind,4.504510e+01_fp_kind,4.224480e+01_fp_kind,&
     3.968550e+01_fp_kind,3.518730e+01_fp_kind,3.136040e+01_fp_kind,2.808880e+01_fp_kind,2.531510e+01_fp_kind,&
     2.299900e+01_fp_kind,1.953950e+01_fp_kind,1.721860e+01_fp_kind,1.553740e+01_fp_kind,1.414140e+01_fp_kind,&
     1.285000e+01_fp_kind,9.878900e+00_fp_kind,7.584500e+00_fp_kind,6.126000e+00_fp_kind,5.264700e+00_fp_kind,&
     4.292600e+00_fp_kind,3.494900e+00_fp_kind], 1.9207115171185121e+00_fp_kind,1.373120e+01_fp_kind) /
data ionicXrayFF( 2, 69) / t_ionicFF('Tm',  69,  3, &
     [6.600000e+01_fp_kind,6.769590e+01_fp_kind,6.485210e+01_fp_kind,6.164850e+01_fp_kind,5.829120e+01_fp_kind,&
     5.482450e+01_fp_kind,5.125150e+01_fp_kind,4.799010e+01_fp_kind,4.495360e+01_fp_kind,4.216760e+01_fp_kind,&
     3.962320e+01_fp_kind,3.515860e+01_fp_kind,3.136340e+01_fp_kind,2.811300e+01_fp_kind,2.534790e+01_fp_kind,&
     2.303070e+01_fp_kind,1.955650e+01_fp_kind,1.722290e+01_fp_kind,1.553720e+01_fp_kind,1.414250e+01_fp_kind,&
     1.285390e+01_fp_kind,9.885600e+00_fp_kind,7.588000e+00_fp_kind,6.124800e+00_fp_kind,5.267800e+00_fp_kind,&
     4.288900e+00_fp_kind,3.494900e+00_fp_kind], 1.9207115171185121e+00_fp_kind,5.732000e+00_fp_kind) /
data ionicXrayFF( 1, 70) / t_ionicFF('Yb',  70,  0, &
     [7.000000e+01_fp_kind,6.871800e+01_fp_kind,6.590590e+01_fp_kind,6.272130e+01_fp_kind,5.937500e+01_fp_kind,&
     5.590810e+01_fp_kind,5.245030e+01_fp_kind,4.913350e+01_fp_kind,4.603640e+01_fp_kind,4.318520e+01_fp_kind,&
     4.057270e+01_fp_kind,3.597180e+01_fp_kind,3.205250e+01_fp_kind,2.869690e+01_fp_kind,2.584440e+01_fp_kind,&
     2.345380e+01_fp_kind,1.986560e+01_fp_kind,1.745880e+01_fp_kind,1.573960e+01_fp_kind,1.434210e+01_fp_kind,&
     1.306710e+01_fp_kind,1.012480e+01_fp_kind,7.793600e+00_fp_kind,6.271500e+00_fp_kind,5.362300e+00_fp_kind,&
     4.365500e+00_fp_kind,3.578100e+00_fp_kind], 1.9392941183555481e+00_fp_kind,1.347280e+01_fp_kind) /
data ionicXrayFF( 2, 70) / t_ionicFF('Yb',  70,  3, &
     [6.700000e+01_fp_kind,6.641700e+01_fp_kind,6.474990e+01_fp_kind,6.221830e+01_fp_kind,5.911540e+01_fp_kind,&
     5.573520e+01_fp_kind,5.231770e+01_fp_kind,4.902500e+01_fp_kind,4.594480e+01_fp_kind,4.310790e+01_fp_kind,&
     4.050980e+01_fp_kind,3.594110e+01_fp_kind,3.205270e+01_fp_kind,2.871880e+01_fp_kind,2.587590e+01_fp_kind,&
     2.348530e+01_fp_kind,1.988370e+01_fp_kind,1.746390e+01_fp_kind,1.573960e+01_fp_kind,1.434260e+01_fp_kind,&
     1.307030e+01_fp_kind,1.013130e+01_fp_kind,7.797400e+00_fp_kind,6.270700e+00_fp_kind,5.365500e+00_fp_kind,&
     4.361800e+00_fp_kind,3.580800e+00_fp_kind], 1.9408363846598600e+00_fp_kind,5.645800e+00_fp_kind) /
data ionicXrayFF( 1, 71) / t_ionicFF('Lu',  71,  0, &
     [7.100000e+01_fp_kind,6.970190e+01_fp_kind,6.677480e+01_fp_kind,6.346170e+01_fp_kind,6.009310e+01_fp_kind,&
     5.668090e+01_fp_kind,5.328680e+01_fp_kind,5.000730e+01_fp_kind,4.691650e+01_fp_kind,4.404920e+01_fp_kind,&
     4.140850e+01_fp_kind,3.674200e+01_fp_kind,3.275950e+01_fp_kind,2.934070e+01_fp_kind,2.641980e+01_fp_kind,&
     2.395580e+01_fp_kind,2.022470e+01_fp_kind,1.713100e+01_fp_kind,1.594290e+01_fp_kind,1.453820e+01_fp_kind,&
     1.327940e+01_fp_kind,1.037370e+01_fp_kind,8.011300e+00_fp_kind,6.425000e+00_fp_kind,5.464000e+00_fp_kind,&
     4.438200e+00_fp_kind,3.659200e+00_fp_kind], 1.9561870366850440e+00_fp_kind,1.348340e+01_fp_kind) /
data ionicXrayFF( 2, 71) / t_ionicFF('Lu',  71,  3, &
     [6.800000e+01_fp_kind,6.742550e+01_fp_kind,6.578020e+01_fp_kind,6.327480e+01_fp_kind,6.019190e+01_fp_kind,&
     5.681770e+01_fp_kind,5.338850e+01_fp_kind,5.006760e+01_fp_kind,4.694670e+01_fp_kind,4.406110e+01_fp_kind,&
     4.141110e+01_fp_kind,3.673980e+01_fp_kind,3.275840e+01_fp_kind,2.934030e+01_fp_kind,2.641920e+01_fp_kind,&
     2.395450e+01_fp_kind,2.022270e+01_fp_kind,1.771160e+01_fp_kind,1.594200e+01_fp_kind,1.453760e+01_fp_kind,&
     1.327870e+01_fp_kind,1.037230e+01_fp_kind,8.009700e+00_fp_kind,6.422500e+00_fp_kind,5.467600e+00_fp_kind,&
     4.433600e+00_fp_kind,3.662100e+00_fp_kind], 1.9578216724905320e+00_fp_kind,5.561800e+00_fp_kind) /
data ionicXrayFF( 1, 72) / t_ionicFF('Hf',  72,  0, &
     [7.200000e+01_fp_kind,7.071420e+01_fp_kind,6.772920e+01_fp_kind,6.428970e+01_fp_kind,6.083920e+01_fp_kind,&
     5.741760e+01_fp_kind,5.404170e+01_fp_kind,5.080330e+01_fp_kind,4.772990e+01_fp_kind,4.486400e+01_fp_kind,&
     4.221250e+01_fp_kind,3.750560e+01_fp_kind,3.347330e+01_fp_kind,2.999860e+01_fp_kind,2.701370e+01_fp_kind,&
     2.447870e+01_fp_kind,2.060330e+01_fp_kind,1.797960e+01_fp_kind,1.614940e+01_fp_kind,1.473050e+01_fp_kind,&
     1.348440e+01_fp_kind,1.061910e+01_fp_kind,8.233200e+00_fp_kind,6.585400e+00_fp_kind,5.571200e+00_fp_kind,&
     4.509200e+00_fp_kind,3.738400e+00_fp_kind], 1.9715682023275110e+00_fp_kind,1.322020e+01_fp_kind) /
data ionicXrayFF( 2, 72) / t_ionicFF('Hf',  72,  0, &
     [7.200000e+01_fp_kind,7.071420e+01_fp_kind,6.772920e+01_fp_kind,6.428970e+01_fp_kind,6.083920e+01_fp_kind,&
     5.741760e+01_fp_kind,5.404170e+01_fp_kind,5.080330e+01_fp_kind,4.772990e+01_fp_kind,4.486400e+01_fp_kind,&
     4.221250e+01_fp_kind,3.750560e+01_fp_kind,3.347330e+01_fp_kind,2.999860e+01_fp_kind,2.701370e+01_fp_kind,&
     2.447870e+01_fp_kind,2.060330e+01_fp_kind,1.797960e+01_fp_kind,1.614940e+01_fp_kind,1.473050e+01_fp_kind,&
     1.348440e+01_fp_kind,1.061910e+01_fp_kind,8.233200e+00_fp_kind,6.585400e+00_fp_kind,5.571200e+00_fp_kind,&
     4.509200e+00_fp_kind,3.738400e+00_fp_kind], 1.9715682023275110e+00_fp_kind,1.322020e+01_fp_kind) /
data ionicXrayFF( 1, 73) / t_ionicFF('Ta',  73,  0, &
     [7.300000e+01_fp_kind,7.173250e+01_fp_kind,6.871900e+01_fp_kind,6.517260e+01_fp_kind,6.162160e+01_fp_kind,&
     5.815300e+01_fp_kind,5.478360e+01_fp_kind,5.155000e+01_fp_kind,4.849200e+01_fp_kind,4.563390e+01_fp_kind,&
     4.298160e+01_fp_kind,3.825400e+01_fp_kind,3.418620e+01_fp_kind,3.066500e+01_fp_kind,2.762250e+01_fp_kind,&
     2.502040e+01_fp_kind,2.100190e+01_fp_kind,1.825990e+01_fp_kind,1.636090e+01_fp_kind,1.492020e+01_fp_kind,&
     1.368250e+01_fp_kind,1.085980e+01_fp_kind,8.458500e+00_fp_kind,6.752600e+00_fp_kind,5.684100e+00_fp_kind,&
     4.579500e+00_fp_kind,3.815800e+00_fp_kind], 1.9855516143859400e+00_fp_kind,1.293160e+01_fp_kind) /
data ionicXrayFF( 2, 73) / t_ionicFF('Ta',  73,  0, &
     [7.300000e+01_fp_kind,7.173250e+01_fp_kind,6.871900e+01_fp_kind,6.517260e+01_fp_kind,6.162160e+01_fp_kind,&
     5.815300e+01_fp_kind,5.478360e+01_fp_kind,5.155000e+01_fp_kind,4.849200e+01_fp_kind,4.563390e+01_fp_kind,&
     4.298160e+01_fp_kind,3.825400e+01_fp_kind,3.418620e+01_fp_kind,3.066500e+01_fp_kind,2.762250e+01_fp_kind,&
     2.502040e+01_fp_kind,2.100190e+01_fp_kind,1.825990e+01_fp_kind,1.636090e+01_fp_kind,1.492020e+01_fp_kind,&
     1.368250e+01_fp_kind,1.085980e+01_fp_kind,8.458500e+00_fp_kind,6.752600e+00_fp_kind,5.684100e+00_fp_kind,&
     4.579500e+00_fp_kind,3.815800e+00_fp_kind], 1.9855516143859400e+00_fp_kind,1.293160e+01_fp_kind) /
data ionicXrayFF( 1, 74) / t_ionicFF('W ',  74,  0, &
     [7.400000e+01_fp_kind,7.275310e+01_fp_kind,6.972920e+01_fp_kind,6.609600e+01_fp_kind,6.244090e+01_fp_kind,&
     5.890240e+01_fp_kind,5.550360e+01_fp_kind,5.226640e+01_fp_kind,4.921500e+01_fp_kind,4.636420e+01_fp_kind,&
     4.371570e+01_fp_kind,3.898170e+01_fp_kind,3.489200e+01_fp_kind,3.133480e+01_fp_kind,2.824260e+01_fp_kind,&
     2.557870e+01_fp_kind,2.142030e+01_fp_kind,1.855530e+01_fp_kind,1.657890e+01_fp_kind,1.510860e+01_fp_kind,&
     1.387430e+01_fp_kind,1.109500e+01_fp_kind,8.686100e+00_fp_kind,6.926400e+00_fp_kind,5.802700e+00_fp_kind,&
     4.649500e+00_fp_kind,3.891300e+00_fp_kind], 1.9981371784100901e+00_fp_kind,1.264340e+01_fp_kind) /
data ionicXrayFF( 2, 74) / t_ionicFF('W ',  74,  0, &
     [7.400000e+01_fp_kind,7.275310e+01_fp_kind,6.972920e+01_fp_kind,6.609600e+01_fp_kind,6.244090e+01_fp_kind,&
     5.890240e+01_fp_kind,5.550360e+01_fp_kind,5.226640e+01_fp_kind,4.921500e+01_fp_kind,4.636420e+01_fp_kind,&
     4.371570e+01_fp_kind,3.898170e+01_fp_kind,3.489200e+01_fp_kind,3.133480e+01_fp_kind,2.824260e+01_fp_kind,&
     2.557870e+01_fp_kind,2.142030e+01_fp_kind,1.855530e+01_fp_kind,1.657890e+01_fp_kind,1.510860e+01_fp_kind,&
     1.387430e+01_fp_kind,1.109500e+01_fp_kind,8.686100e+00_fp_kind,6.926400e+00_fp_kind,5.802700e+00_fp_kind,&
     4.649500e+00_fp_kind,3.891300e+00_fp_kind], 1.9981371784100901e+00_fp_kind,1.264340e+01_fp_kind) /
data ionicXrayFF( 1, 75) / t_ionicFF('Re',  75,  0, &
     [7.500000e+01_fp_kind,7.377470e+01_fp_kind,7.075290e+01_fp_kind,6.705040e+01_fp_kind,6.329490e+01_fp_kind,&
     5.967350e+01_fp_kind,5.622530e+01_fp_kind,5.296660e+01_fp_kind,4.991030e+01_fp_kind,4.706190e+01_fp_kind,&
     4.441720e+01_fp_kind,3.968560e+01_fp_kind,3.558560e+01_fp_kind,3.200330e+01_fp_kind,2.887010e+01_fp_kind,&
     2.615090e+01_fp_kind,2.185810e+01_fp_kind,1.886690e+01_fp_kind,1.680540e+01_fp_kind,1.529720e+01_fp_kind,&
     1.406080e+01_fp_kind,1.132410e+01_fp_kind,8.915200e+00_fp_kind,7.106400e+00_fp_kind,5.927300e+00_fp_kind,&
     4.719800e+00_fp_kind,3.964900e+00_fp_kind], 2.0093784621968571e+00_fp_kind,1.236270e+01_fp_kind) /
data ionicXrayFF( 2, 75) / t_ionicFF('Re',  75,  0, &
     [7.500000e+01_fp_kind,7.377470e+01_fp_kind,7.075290e+01_fp_kind,6.705040e+01_fp_kind,6.329490e+01_fp_kind,&
     5.967350e+01_fp_kind,5.622530e+01_fp_kind,5.296660e+01_fp_kind,4.991030e+01_fp_kind,4.706190e+01_fp_kind,&
     4.441720e+01_fp_kind,3.968560e+01_fp_kind,3.558560e+01_fp_kind,3.200330e+01_fp_kind,2.887010e+01_fp_kind,&
     2.615090e+01_fp_kind,2.185810e+01_fp_kind,1.886690e+01_fp_kind,1.680540e+01_fp_kind,1.529720e+01_fp_kind,&
     1.406080e+01_fp_kind,1.132410e+01_fp_kind,8.915200e+00_fp_kind,7.106400e+00_fp_kind,5.927300e+00_fp_kind,&
     4.719800e+00_fp_kind,3.964900e+00_fp_kind], 2.0093784621968571e+00_fp_kind,1.236270e+01_fp_kind) /
data ionicXrayFF( 1, 76) / t_ionicFF('Os',  76,  0, &
     [7.600000e+01_fp_kind,7.479650e+01_fp_kind,7.178540e+01_fp_kind,6.802920e+01_fp_kind,6.418030e+01_fp_kind,&
     6.046950e+01_fp_kind,5.695740e+01_fp_kind,5.366140e+01_fp_kind,5.058770e+01_fp_kind,4.773400e+01_fp_kind,&
     4.509010e+01_fp_kind,4.036430e+01_fp_kind,3.626330e+01_fp_kind,3.266620e+01_fp_kind,2.950130e+01_fp_kind,&
     2.673420e+01_fp_kind,2.231470e+01_fp_kind,1.919570e+01_fp_kind,1.704190e+01_fp_kind,1.548780e+01_fp_kind,&
     1.424300e+01_fp_kind,1.154660e+01_fp_kind,9.144900e+00_fp_kind,7.292200e+00_fp_kind,6.058100e+00_fp_kind,&
     4.790900e+00_fp_kind,4.036800e+00_fp_kind], 2.0194321542121472e+00_fp_kind,1.209170e+01_fp_kind) /
data ionicXrayFF( 2, 76) / t_ionicFF('Os',  76,  0, &
     [7.600000e+01_fp_kind,7.479650e+01_fp_kind,7.178540e+01_fp_kind,6.802920e+01_fp_kind,6.418030e+01_fp_kind,&
     6.046950e+01_fp_kind,5.695740e+01_fp_kind,5.366140e+01_fp_kind,5.058770e+01_fp_kind,4.773400e+01_fp_kind,&
     4.509010e+01_fp_kind,4.036430e+01_fp_kind,3.626330e+01_fp_kind,3.266620e+01_fp_kind,2.950130e+01_fp_kind,&
     2.673420e+01_fp_kind,2.231470e+01_fp_kind,1.919570e+01_fp_kind,1.704190e+01_fp_kind,1.548780e+01_fp_kind,&
     1.424300e+01_fp_kind,1.154660e+01_fp_kind,9.144900e+00_fp_kind,7.292200e+00_fp_kind,6.058100e+00_fp_kind,&
     4.790900e+00_fp_kind,4.036800e+00_fp_kind], 2.0194321542121472e+00_fp_kind,1.209170e+01_fp_kind) /
data ionicXrayFF( 1, 77) / t_ionicFF('Ir',  77,  0, &
     [7.700000e+01_fp_kind,7.582330e+01_fp_kind,7.284040e+01_fp_kind,6.905300e+01_fp_kind,6.512310e+01_fp_kind,&
     6.131890e+01_fp_kind,5.772690e+01_fp_kind,5.437350e+01_fp_kind,5.126360e+01_fp_kind,4.839050e+01_fp_kind,&
     4.573860e+01_fp_kind,4.101430e+01_fp_kind,3.691930e+01_fp_kind,3.331830e+01_fp_kind,3.013280e+01_fp_kind,&
     2.732710e+01_fp_kind,2.279080e+01_fp_kind,1.954340e+01_fp_kind,1.729040e+01_fp_kind,1.568160e+01_fp_kind,&
     1.442180e+01_fp_kind,1.176210e+01_fp_kind,9.374400e+00_fp_kind,7.483700e+00_fp_kind,6.195400e+00_fp_kind,&
     4.863500e+00_fp_kind,4.106900e+00_fp_kind], 2.0282907435683200e+00_fp_kind,1.177930e+01_fp_kind) /
data ionicXrayFF( 2, 77) / t_ionicFF('Ir',  77,  0, &
     [7.700000e+01_fp_kind,7.582330e+01_fp_kind,7.284040e+01_fp_kind,6.905300e+01_fp_kind,6.512310e+01_fp_kind,&
     6.131890e+01_fp_kind,5.772690e+01_fp_kind,5.437350e+01_fp_kind,5.126360e+01_fp_kind,4.839050e+01_fp_kind,&
     4.573860e+01_fp_kind,4.101430e+01_fp_kind,3.691930e+01_fp_kind,3.331830e+01_fp_kind,3.013280e+01_fp_kind,&
     2.732710e+01_fp_kind,2.279080e+01_fp_kind,1.954340e+01_fp_kind,1.729040e+01_fp_kind,1.568160e+01_fp_kind,&
     1.442180e+01_fp_kind,1.176210e+01_fp_kind,9.374400e+00_fp_kind,7.483700e+00_fp_kind,6.195400e+00_fp_kind,&
     4.863500e+00_fp_kind,4.106900e+00_fp_kind], 2.0282907435683200e+00_fp_kind,1.177930e+01_fp_kind) /
data ionicXrayFF( 1, 78) / t_ionicFF('Pt',  78,  0, &
     [7.800000e+01_fp_kind,7.691100e+01_fp_kind,7.407490e+01_fp_kind,7.032480e+01_fp_kind,6.629560e+01_fp_kind,&
     6.232510e+01_fp_kind,5.856860e+01_fp_kind,5.508810e+01_fp_kind,5.189690e+01_fp_kind,4.798180e+01_fp_kind,&
     4.631700e+01_fp_kind,4.161250e+01_fp_kind,3.755100e+01_fp_kind,3.396790e+01_fp_kind,3.077490e+01_fp_kind,&
     2.793690e+01_fp_kind,2.328630e+01_fp_kind,1.990720e+01_fp_kind,1.754920e+01_fp_kind,1.587970e+01_fp_kind,&
     1.459900e+01_fp_kind,1.197140e+01_fp_kind,9.603300e+00_fp_kind,7.680100e+00_fp_kind,6.339700e+00_fp_kind,&
     4.937400e+00_fp_kind,4.175700e+00_fp_kind], 2.0362563546149488e+00_fp_kind,1.082250e+01_fp_kind) /
data ionicXrayFF( 2, 78) / t_ionicFF('Pt',  78,  0, &
     [7.800000e+01_fp_kind,7.691100e+01_fp_kind,7.407490e+01_fp_kind,7.032480e+01_fp_kind,6.629560e+01_fp_kind,&
     6.232510e+01_fp_kind,5.856860e+01_fp_kind,5.508810e+01_fp_kind,5.189690e+01_fp_kind,4.798180e+01_fp_kind,&
     4.631700e+01_fp_kind,4.161250e+01_fp_kind,3.755100e+01_fp_kind,3.396790e+01_fp_kind,3.077490e+01_fp_kind,&
     2.793690e+01_fp_kind,2.328630e+01_fp_kind,1.990720e+01_fp_kind,1.754920e+01_fp_kind,1.587970e+01_fp_kind,&
     1.459900e+01_fp_kind,1.197140e+01_fp_kind,9.603300e+00_fp_kind,7.680100e+00_fp_kind,6.339700e+00_fp_kind,&
     4.937400e+00_fp_kind,4.175700e+00_fp_kind], 2.0362563546149488e+00_fp_kind,1.082250e+01_fp_kind) /
data ionicXrayFF( 1, 79) / t_ionicFF('Au',  79,  0, &
     [7.900000e+01_fp_kind,7.793580e+01_fp_kind,7.513710e+01_fp_kind,7.138370e+01_fp_kind,6.730100e+01_fp_kind,&
     6.324740e+01_fp_kind,5.940180e+01_fp_kind,5.584180e+01_fp_kind,5.258790e+01_fp_kind,4.962860e+01_fp_kind,&
     4.693620e+01_fp_kind,4.221370e+01_fp_kind,3.816040e+01_fp_kind,3.458820e+01_fp_kind,3.139420e+01_fp_kind,&
     2.852650e+01_fp_kind,2.379470e+01_fp_kind,2.029230e+01_fp_kind,1.782520e+01_fp_kind,1.608550e+01_fp_kind,&
     1.477550e+01_fp_kind,1.217220e+01_fp_kind,9.829700e+00_fp_kind,7.880600e+00_fp_kind,6.489700e+00_fp_kind,&
     5.014400e+00_fp_kind,4.243100e+00_fp_kind], 2.0433110522239351e+00_fp_kind,1.054660e+01_fp_kind) /
data ionicXrayFF( 2, 79) / t_ionicFF('Au',  79,  0, &
     [7.900000e+01_fp_kind,7.793580e+01_fp_kind,7.513710e+01_fp_kind,7.138370e+01_fp_kind,6.730100e+01_fp_kind,&
     6.324740e+01_fp_kind,5.940180e+01_fp_kind,5.584180e+01_fp_kind,5.258790e+01_fp_kind,4.962860e+01_fp_kind,&
     4.693620e+01_fp_kind,4.221370e+01_fp_kind,3.816040e+01_fp_kind,3.458820e+01_fp_kind,3.139420e+01_fp_kind,&
     2.852650e+01_fp_kind,2.379470e+01_fp_kind,2.029230e+01_fp_kind,1.782520e+01_fp_kind,1.608550e+01_fp_kind,&
     1.477550e+01_fp_kind,1.217220e+01_fp_kind,9.829700e+00_fp_kind,7.880600e+00_fp_kind,6.489700e+00_fp_kind,&
     5.014400e+00_fp_kind,4.243100e+00_fp_kind], 2.0433110522239351e+00_fp_kind,1.054660e+01_fp_kind) /
data ionicXrayFF( 1, 80) / t_ionicFF('Hg',  80,  0, &
     [8.000000e+01_fp_kind,7.889830e+01_fp_kind,7.601950e+01_fp_kind,7.220200e+01_fp_kind,6.809340e+01_fp_kind,&
     6.403500e+01_fp_kind,6.018360e+01_fp_kind,5.660670e+01_fp_kind,5.332520e+01_fp_kind,5.033320e+01_fp_kind,&
     4.760820e+01_fp_kind,4.283600e+01_fp_kind,3.876000e+01_fp_kind,3.518330e+01_fp_kind,3.198750e+01_fp_kind,&
     2.911870e+01_fp_kind,2.430930e+01_fp_kind,2.069670e+01_fp_kind,1.812000e+01_fp_kind,1.630220e+01_fp_kind,&
     1.495370e+01_fp_kind,1.236450e+01_fp_kind,1.005240e+01_fp_kind,8.084300e+00_fp_kind,6.645100e+00_fp_kind,&
     5.094800e+00_fp_kind,4.309100e+00_fp_kind], 2.0494881154801941e+00_fp_kind,1.093560e+01_fp_kind) /
data ionicXrayFF( 2, 80) / t_ionicFF('Hg',  80,  0, &
     [8.000000e+01_fp_kind,7.889830e+01_fp_kind,7.601950e+01_fp_kind,7.220200e+01_fp_kind,6.809340e+01_fp_kind,&
     6.403500e+01_fp_kind,6.018360e+01_fp_kind,5.660670e+01_fp_kind,5.332520e+01_fp_kind,5.033320e+01_fp_kind,&
     4.760820e+01_fp_kind,4.283600e+01_fp_kind,3.876000e+01_fp_kind,3.518330e+01_fp_kind,3.198750e+01_fp_kind,&
     2.911870e+01_fp_kind,2.430930e+01_fp_kind,2.069670e+01_fp_kind,1.812000e+01_fp_kind,1.630220e+01_fp_kind,&
     1.495370e+01_fp_kind,1.236450e+01_fp_kind,1.005240e+01_fp_kind,8.084300e+00_fp_kind,6.645100e+00_fp_kind,&
     5.094800e+00_fp_kind,4.309100e+00_fp_kind], 2.0494881154801941e+00_fp_kind,1.093560e+01_fp_kind) /
data ionicXrayFF( 1, 81) / t_ionicFF('Tl',  81,  0, &
     [8.100000e+01_fp_kind,7.974250e+01_fp_kind,7.668520e+01_fp_kind,7.286940e+01_fp_kind,6.883480e+01_fp_kind,&
     6.483160e+01_fp_kind,6.099570e+01_fp_kind,5.740160e+01_fp_kind,5.408370e+01_fp_kind,5.104770e+01_fp_kind,&
     4.828010e+01_fp_kind,4.344530e+01_fp_kind,3.934220e+01_fp_kind,3.576240e+01_fp_kind,3.257000e+01_fp_kind,&
     2.969680e+01_fp_kind,2.483250e+01_fp_kind,2.111660e+01_fp_kind,1.842980e+01_fp_kind,1.652860e+01_fp_kind,&
     1.513480e+01_fp_kind,1.255020e+01_fp_kind,1.027140e+01_fp_kind,8.291100e+00_fp_kind,6.805500e+00_fp_kind,&
     5.178700e+00_fp_kind,4.374400e+00_fp_kind], 2.0551669415965388e+00_fp_kind,1.280450e+01_fp_kind) /
data ionicXrayFF( 2, 81) / t_ionicFF('Tl',  81,  0, &
     [8.100000e+01_fp_kind,7.974250e+01_fp_kind,7.668520e+01_fp_kind,7.286940e+01_fp_kind,6.883480e+01_fp_kind,&
     6.483160e+01_fp_kind,6.099570e+01_fp_kind,5.740160e+01_fp_kind,5.408370e+01_fp_kind,5.104770e+01_fp_kind,&
     4.828010e+01_fp_kind,4.344530e+01_fp_kind,3.934220e+01_fp_kind,3.576240e+01_fp_kind,3.257000e+01_fp_kind,&
     2.969680e+01_fp_kind,2.483250e+01_fp_kind,2.111660e+01_fp_kind,1.842980e+01_fp_kind,1.652860e+01_fp_kind,&
     1.513480e+01_fp_kind,1.255020e+01_fp_kind,1.027140e+01_fp_kind,8.291100e+00_fp_kind,6.805500e+00_fp_kind,&
     5.178700e+00_fp_kind,4.374400e+00_fp_kind], 2.0551669415965388e+00_fp_kind,1.280450e+01_fp_kind) /
data ionicXrayFF( 1, 82) / t_ionicFF('Pb',  82,  0, &
     [8.200000e+01_fp_kind,8.066860e+01_fp_kind,7.744710e+01_fp_kind,7.352830e+01_fp_kind,6.949950e+01_fp_kind,&
     6.554940e+01_fp_kind,6.175790e+01_fp_kind,5.817880e+01_fp_kind,5.484730e+01_fp_kind,5.177900e+01_fp_kind,&
     4.897130e+01_fp_kind,4.406330e+01_fp_kind,3.991870e+01_fp_kind,3.632700e+01_fp_kind,3.313630e+01_fp_kind,&
     3.026300e+01_fp_kind,2.535910e+01_fp_kind,2.155100e+01_fp_kind,1.875630e+01_fp_kind,1.676700e+01_fp_kind,&
     1.532080e+01_fp_kind,1.272880e+01_fp_kind,1.048580e+01_fp_kind,8.499200e+00_fp_kind,6.971700e+00_fp_kind,&
     5.266500e+00_fp_kind,4.438800e+00_fp_kind], 2.0602672470255778e+00_fp_kind,1.353440e+01_fp_kind) /
data ionicXrayFF( 2, 82) / t_ionicFF('Pb',  82,  0, &
     [8.200000e+01_fp_kind,8.066860e+01_fp_kind,7.744710e+01_fp_kind,7.352830e+01_fp_kind,6.949950e+01_fp_kind,&
     6.554940e+01_fp_kind,6.175790e+01_fp_kind,5.817880e+01_fp_kind,5.484730e+01_fp_kind,5.177900e+01_fp_kind,&
     4.897130e+01_fp_kind,4.406330e+01_fp_kind,3.991870e+01_fp_kind,3.632700e+01_fp_kind,3.313630e+01_fp_kind,&
     3.026300e+01_fp_kind,2.535910e+01_fp_kind,2.155100e+01_fp_kind,1.875630e+01_fp_kind,1.676700e+01_fp_kind,&
     1.532080e+01_fp_kind,1.272880e+01_fp_kind,1.048580e+01_fp_kind,8.499200e+00_fp_kind,6.971700e+00_fp_kind,&
     5.266500e+00_fp_kind,4.438800e+00_fp_kind], 2.0602672470255778e+00_fp_kind,1.353440e+01_fp_kind) /
data ionicXrayFF( 1, 83) / t_ionicFF('Bi',  83,  0, &
     [8.300000e+01_fp_kind,8.162720e+01_fp_kind,7.827470e+01_fp_kind,7.421740e+01_fp_kind,7.013520e+01_fp_kind,&
     6.620820e+01_fp_kind,6.246490e+01_fp_kind,5.892120e+01_fp_kind,5.559870e+01_fp_kind,5.251500e+01_fp_kind,&
     4.967610e+01_fp_kind,4.469380e+01_fp_kind,4.049550e+01_fp_kind,3.688100e+01_fp_kind,3.368790e+01_fp_kind,&
     3.081620e+01_fp_kind,2.588580e+01_fp_kind,2.199770e+01_fp_kind,1.909920e+01_fp_kind,1.701870e+01_fp_kind,&
     1.551330e+01_fp_kind,1.290120e+01_fp_kind,1.069490e+01_fp_kind,8.708000e+00_fp_kind,7.143400e+00_fp_kind,&
     5.358000e+00_fp_kind,4.503000e+00_fp_kind], 2.0651489865854771e+00_fp_kind,1.389110e+01_fp_kind) /
data ionicXrayFF( 2, 83) / t_ionicFF('Bi',  83,  0, &
     [8.300000e+01_fp_kind,8.162720e+01_fp_kind,7.827470e+01_fp_kind,7.421740e+01_fp_kind,7.013520e+01_fp_kind,&
     6.620820e+01_fp_kind,6.246490e+01_fp_kind,5.892120e+01_fp_kind,5.559870e+01_fp_kind,5.251500e+01_fp_kind,&
     4.967610e+01_fp_kind,4.469380e+01_fp_kind,4.049550e+01_fp_kind,3.688100e+01_fp_kind,3.368790e+01_fp_kind,&
     3.081620e+01_fp_kind,2.588580e+01_fp_kind,2.199770e+01_fp_kind,1.909920e+01_fp_kind,1.701870e+01_fp_kind,&
     1.551330e+01_fp_kind,1.290120e+01_fp_kind,1.069490e+01_fp_kind,8.708000e+00_fp_kind,7.143400e+00_fp_kind,&
     5.358000e+00_fp_kind,4.503000e+00_fp_kind], 2.0651489865854771e+00_fp_kind,1.389110e+01_fp_kind) /
data ionicXrayFF( 1, 84) / t_ionicFF('Po',  84,  0, &
     [8.400000e+01_fp_kind,8.260370e+01_fp_kind,7.914980e+01_fp_kind,7.494670e+01_fp_kind,7.076980e+01_fp_kind,&
     6.682950e+01_fp_kind,6.312180e+01_fp_kind,5.962150e+01_fp_kind,5.632560e+01_fp_kind,5.324420e+01_fp_kind,&
     5.038680e+01_fp_kind,4.533760e+01_fp_kind,4.107730e+01_fp_kind,3.742910e+01_fp_kind,3.422660e+01_fp_kind,&
     3.135600e+01_fp_kind,2.640960e+01_fp_kind,2.245430e+01_fp_kind,1.945780e+01_fp_kind,1.728430e+01_fp_kind,&
     1.571380e+01_fp_kind,1.306840e+01_fp_kind,1.089820e+01_fp_kind,8.916800e+00_fp_kind,7.320000e+00_fp_kind,&
     5.453300e+00_fp_kind,4.567500e+00_fp_kind], 2.0700594844679472e+00_fp_kind,1.406290e+01_fp_kind) /
data ionicXrayFF( 2, 84) / t_ionicFF('Po',  84,  0, &
     [8.400000e+01_fp_kind,8.260370e+01_fp_kind,7.914980e+01_fp_kind,7.494670e+01_fp_kind,7.076980e+01_fp_kind,&
     6.682950e+01_fp_kind,6.312180e+01_fp_kind,5.962150e+01_fp_kind,5.632560e+01_fp_kind,5.324420e+01_fp_kind,&
     5.038680e+01_fp_kind,4.533760e+01_fp_kind,4.107730e+01_fp_kind,3.742910e+01_fp_kind,3.422660e+01_fp_kind,&
     3.135600e+01_fp_kind,2.640960e+01_fp_kind,2.245430e+01_fp_kind,1.945780e+01_fp_kind,1.728430e+01_fp_kind,&
     1.571380e+01_fp_kind,1.306840e+01_fp_kind,1.089820e+01_fp_kind,8.916800e+00_fp_kind,7.320000e+00_fp_kind,&
     5.453300e+00_fp_kind,4.567500e+00_fp_kind], 2.0700594844679472e+00_fp_kind,1.406290e+01_fp_kind) /
data ionicXrayFF( 1, 85) / t_ionicFF('At',  85,  0, &
     [8.500000e+01_fp_kind,8.362540e+01_fp_kind,8.015080e+01_fp_kind,7.582170e+01_fp_kind,7.149370e+01_fp_kind,&
     6.746210e+01_fp_kind,6.373650e+01_fp_kind,6.026200e+01_fp_kind,5.700150e+01_fp_kind,5.394360e+01_fp_kind,&
     5.108990e+01_fp_kind,4.600000e+01_fp_kind,4.167720e+01_fp_kind,3.798200e+01_fp_kind,3.475760e+01_fp_kind,&
     3.188180e+01_fp_kind,2.692580e+01_fp_kind,2.291790e+01_fp_kind,1.983220e+01_fp_kind,1.756600e+01_fp_kind,&
     1.592500e+01_fp_kind,1.323080e+01_fp_kind,1.109520e+01_fp_kind,9.125100e+00_fp_kind,7.501500e+00_fp_kind,&
     5.552500e+00_fp_kind,4.632800e+00_fp_kind], 2.0752346728516140e+00_fp_kind,1.375290e+01_fp_kind) /
data ionicXrayFF( 2, 85) / t_ionicFF('At',  85,  0, &
     [8.500000e+01_fp_kind,8.362540e+01_fp_kind,8.015080e+01_fp_kind,7.582170e+01_fp_kind,7.149370e+01_fp_kind,&
     6.746210e+01_fp_kind,6.373650e+01_fp_kind,6.026200e+01_fp_kind,5.700150e+01_fp_kind,5.394360e+01_fp_kind,&
     5.108990e+01_fp_kind,4.600000e+01_fp_kind,4.167720e+01_fp_kind,3.798200e+01_fp_kind,3.475760e+01_fp_kind,&
     3.188180e+01_fp_kind,2.692580e+01_fp_kind,2.291790e+01_fp_kind,1.983220e+01_fp_kind,1.756600e+01_fp_kind,&
     1.592500e+01_fp_kind,1.323080e+01_fp_kind,1.109520e+01_fp_kind,9.125100e+00_fp_kind,7.501500e+00_fp_kind,&
     5.552500e+00_fp_kind,4.632800e+00_fp_kind], 2.0752346728516140e+00_fp_kind,1.375290e+01_fp_kind) /
data ionicXrayFF( 1, 86) / t_ionicFF('Rn',  86,  0, &
     [8.600000e+01_fp_kind,8.464830e+01_fp_kind,8.117070e+01_fp_kind,7.674100e+01_fp_kind,7.226310e+01_fp_kind,&
     6.811430e+01_fp_kind,6.433870e+01_fp_kind,6.086900e+01_fp_kind,5.763880e+01_fp_kind,5.461150e+01_fp_kind,&
     5.177460e+01_fp_kind,4.666740e+01_fp_kind,4.228910e+01_fp_kind,3.854100e+01_fp_kind,3.528520e+01_fp_kind,&
     3.239780e+01_fp_kind,2.743360e+01_fp_kind,2.338510e+01_fp_kind,2.021980e+01_fp_kind,1.786280e+01_fp_kind,&
     1.614750e+01_fp_kind,1.339010e+01_fp_kind,1.128580e+01_fp_kind,9.332100e+00_fp_kind,7.686700e+00_fp_kind,&
     5.656200e+00_fp_kind,4.699100e+00_fp_kind], 2.0807592535875941e+00_fp_kind,1.345460e+01_fp_kind) /
data ionicXrayFF( 2, 86) / t_ionicFF('Rn',  86,  0, &
     [8.600000e+01_fp_kind,8.464830e+01_fp_kind,8.117070e+01_fp_kind,7.674100e+01_fp_kind,7.226310e+01_fp_kind,&
     6.811430e+01_fp_kind,6.433870e+01_fp_kind,6.086900e+01_fp_kind,5.763880e+01_fp_kind,5.461150e+01_fp_kind,&
     5.177460e+01_fp_kind,4.666740e+01_fp_kind,4.228910e+01_fp_kind,3.854100e+01_fp_kind,3.528520e+01_fp_kind,&
     3.239780e+01_fp_kind,2.743360e+01_fp_kind,2.338510e+01_fp_kind,2.021980e+01_fp_kind,1.786280e+01_fp_kind,&
     1.614750e+01_fp_kind,1.339010e+01_fp_kind,1.128580e+01_fp_kind,9.332100e+00_fp_kind,7.686700e+00_fp_kind,&
     5.656200e+00_fp_kind,4.699100e+00_fp_kind], 2.0807592535875941e+00_fp_kind,1.345460e+01_fp_kind) /
data ionicXrayFF( 1, 87) / t_ionicFF('Fr',  87,  0, &
     [8.700000e+01_fp_kind,8.528730e+01_fp_kind,8.166870e+01_fp_kind,7.743400e+01_fp_kind,7.304040e+01_fp_kind,&
     6.884720e+01_fp_kind,6.500350e+01_fp_kind,6.149670e+01_fp_kind,5.826340e+01_fp_kind,5.525030e+01_fp_kind,&
     5.242800e+01_fp_kind,4.731910e+01_fp_kind,4.290050e+01_fp_kind,3.910380e+01_fp_kind,3.581330e+01_fp_kind,&
     3.290910e+01_fp_kind,2.793470e+01_fp_kind,2.385340e+01_fp_kind,2.061690e+01_fp_kind,1.817230e+01_fp_kind,&
     1.638080e+01_fp_kind,1.354840e+01_fp_kind,1.147080e+01_fp_kind,9.538100e+00_fp_kind,7.873000e+00_fp_kind,&
     5.764100e+00_fp_kind,4.767200e+00_fp_kind], 2.0869920518333340e+00_fp_kind,1.855130e+01_fp_kind) /
data ionicXrayFF( 2, 87) / t_ionicFF('Fr',  87,  0, &
     [8.700000e+01_fp_kind,8.528730e+01_fp_kind,8.166870e+01_fp_kind,7.743400e+01_fp_kind,7.304040e+01_fp_kind,&
     6.884720e+01_fp_kind,6.500350e+01_fp_kind,6.149670e+01_fp_kind,5.826340e+01_fp_kind,5.525030e+01_fp_kind,&
     5.242800e+01_fp_kind,4.731910e+01_fp_kind,4.290050e+01_fp_kind,3.910380e+01_fp_kind,3.581330e+01_fp_kind,&
     3.290910e+01_fp_kind,2.793470e+01_fp_kind,2.385340e+01_fp_kind,2.061690e+01_fp_kind,1.817230e+01_fp_kind,&
     1.638080e+01_fp_kind,1.354840e+01_fp_kind,1.147080e+01_fp_kind,9.538100e+00_fp_kind,7.873000e+00_fp_kind,&
     5.764100e+00_fp_kind,4.767200e+00_fp_kind], 2.0869920518333340e+00_fp_kind,1.855130e+01_fp_kind) /
data ionicXrayFF( 1, 88) / t_ionicFF('Ra',  88,  0, &
     [8.800000e+01_fp_kind,8.610540e+01_fp_kind,8.220570e+01_fp_kind,7.799360e+01_fp_kind,7.373360e+01_fp_kind,&
     6.958280e+01_fp_kind,6.570320e+01_fp_kind,6.214750e+01_fp_kind,4.888690e+01_fp_kind,5.587010e+01_fp_kind,&
     5.305620e+01_fp_kind,4.795660e+01_fp_kind,4.351320e+01_fp_kind,3.967330e+01_fp_kind,3.634450e+01_fp_kind,&
     3.341710e+01_fp_kind,2.842740e+01_fp_kind,2.432060e+01_fp_kind,2.102260e+01_fp_kind,1.849580e+01_fp_kind,&
     1.662670e+01_fp_kind,1.370540e+01_fp_kind,1.164920e+01_fp_kind,9.741500e+00_fp_kind,8.061000e+00_fp_kind,&
     5.879200e+00_fp_kind,4.836300e+00_fp_kind], 2.0935432165716530e+00_fp_kind,2.039260e+01_fp_kind) /
data ionicXrayFF( 2, 88) / t_ionicFF('Ra',  88,  0, &
     [8.800000e+01_fp_kind,8.610540e+01_fp_kind,8.220570e+01_fp_kind,7.799360e+01_fp_kind,7.373360e+01_fp_kind,&
     6.958280e+01_fp_kind,6.570320e+01_fp_kind,6.214750e+01_fp_kind,4.888690e+01_fp_kind,5.587010e+01_fp_kind,&
     5.305620e+01_fp_kind,4.795660e+01_fp_kind,4.351320e+01_fp_kind,3.967330e+01_fp_kind,3.634450e+01_fp_kind,&
     3.341710e+01_fp_kind,2.842740e+01_fp_kind,2.432060e+01_fp_kind,2.102260e+01_fp_kind,1.849580e+01_fp_kind,&
     1.662670e+01_fp_kind,1.370540e+01_fp_kind,1.164920e+01_fp_kind,9.741500e+00_fp_kind,8.061000e+00_fp_kind,&
     5.879200e+00_fp_kind,4.836300e+00_fp_kind], 2.0935432165716530e+00_fp_kind,2.039260e+01_fp_kind) /
data ionicXrayFF( 1, 89) / t_ionicFF('Ac',  89,  0, &
     [8.900000e+01_fp_kind,8.706800e+01_fp_kind,8.296470e+01_fp_kind,7.858970e+01_fp_kind,7.432500e+01_fp_kind,&
     7.021770e+01_fp_kind,6.635830e+01_fp_kind,6.279840e+01_fp_kind,5.952770e+01_fp_kind,5.650540e+01_fp_kind,&
     5.369090e+01_fp_kind,4.858710e+01_fp_kind,4.411830e+01_fp_kind,4.023860e+01_fp_kind,3.687300e+01_fp_kind,&
     3.392210e+01_fp_kind,2.891530e+01_fp_kind,2.478770e+01_fp_kind,2.143550e+01_fp_kind,1.883090e+01_fp_kind,&
     1.688420e+01_fp_kind,1.386270e+01_fp_kind,1.182150e+01_fp_kind,9.940600e+00_fp_kind,8.251600e+00_fp_kind,&
     6.000200e+00_fp_kind,4.906300e+00_fp_kind], 2.1003571016616012e+00_fp_kind,2.046780e+01_fp_kind) /
data ionicXrayFF( 2, 89) / t_ionicFF('Ac',  89,  0, &
     [8.900000e+01_fp_kind,8.706800e+01_fp_kind,8.296470e+01_fp_kind,7.858970e+01_fp_kind,7.432500e+01_fp_kind,&
     7.021770e+01_fp_kind,6.635830e+01_fp_kind,6.279840e+01_fp_kind,5.952770e+01_fp_kind,5.650540e+01_fp_kind,&
     5.369090e+01_fp_kind,4.858710e+01_fp_kind,4.411830e+01_fp_kind,4.023860e+01_fp_kind,3.687300e+01_fp_kind,&
     3.392210e+01_fp_kind,2.891530e+01_fp_kind,2.478770e+01_fp_kind,2.143550e+01_fp_kind,1.883090e+01_fp_kind,&
     1.688420e+01_fp_kind,1.386270e+01_fp_kind,1.182150e+01_fp_kind,9.940600e+00_fp_kind,8.251600e+00_fp_kind,&
     6.000200e+00_fp_kind,4.906300e+00_fp_kind], 2.1003571016616012e+00_fp_kind,2.046780e+01_fp_kind) /
data ionicXrayFF( 1, 90) / t_ionicFF('Th',  90,  0, &
     [9.000000e+01_fp_kind,8.806780e+01_fp_kind,8.382340e+01_fp_kind,7.925200e+01_fp_kind,7.490530e+01_fp_kind,&
     7.080490e+01_fp_kind,6.697080e+01_fp_kind,6.342580e+01_fp_kind,6.015970e+01_fp_kind,5.713820e+01_fp_kind,&
     5.432400e+01_fp_kind,4.921450e+01_fp_kind,4.472180e+01_fp_kind,4.080400e+01_fp_kind,3.740100e+01_fp_kind,&
     3.442360e+01_fp_kind,2.939670e+01_fp_kind,2.525260e+01_fp_kind,2.185420e+01_fp_kind,1.917730e+01_fp_kind,&
     1.715390e+01_fp_kind,1.402180e+01_fp_kind,1.198810e+01_fp_kind,1.013560e+01_fp_kind,8.444100e+00_fp_kind,&
     6.126000e+00_fp_kind,4.978500e+00_fp_kind], 2.1080079744535212e+00_fp_kind,2.020220e+01_fp_kind) /
data ionicXrayFF( 2, 90) / t_ionicFF('Th',  90,  0, &
     [9.000000e+01_fp_kind,8.806780e+01_fp_kind,8.382340e+01_fp_kind,7.925200e+01_fp_kind,7.490530e+01_fp_kind,&
     7.080490e+01_fp_kind,6.697080e+01_fp_kind,6.342580e+01_fp_kind,6.015970e+01_fp_kind,5.713820e+01_fp_kind,&
     5.432400e+01_fp_kind,4.921450e+01_fp_kind,4.472180e+01_fp_kind,4.080400e+01_fp_kind,3.740100e+01_fp_kind,&
     3.442360e+01_fp_kind,2.939670e+01_fp_kind,2.525260e+01_fp_kind,2.185420e+01_fp_kind,1.917730e+01_fp_kind,&
     1.715390e+01_fp_kind,1.402180e+01_fp_kind,1.198810e+01_fp_kind,1.013560e+01_fp_kind,8.444100e+00_fp_kind,&
     6.126000e+00_fp_kind,4.978500e+00_fp_kind], 2.1080079744535212e+00_fp_kind,2.020220e+01_fp_kind) /
data ionicXrayFF( 1, 91) / t_ionicFF('Pa',  91,  0, &
     [9.100000e+01_fp_kind,8.912750e+01_fp_kind,8.502510e+01_fp_kind,8.051980e+01_fp_kind,7.608760e+01_fp_kind,&
     7.181680e+01_fp_kind,6.780710e+01_fp_kind,6.412150e+01_fp_kind,6.075910e+01_fp_kind,5.768090e+01_fp_kind,&
     5.483960e+01_fp_kind,4.972680e+01_fp_kind,4.525210e+01_fp_kind,4.134580e+01_fp_kind,3.794240e+01_fp_kind,&
     3.495680e+01_fp_kind,2.990550e+01_fp_kind,2.572680e+01_fp_kind,2.227390e+01_fp_kind,1.952540e+01_fp_kind,&
     1.742790e+01_fp_kind,1.418090e+01_fp_kind,1.214500e+01_fp_kind,1.031970e+01_fp_kind,8.631800e+00_fp_kind,&
     6.252600e+00_fp_kind,5.054900e+00_fp_kind], 2.1173563123435808e+00_fp_kind,1.964040e+01_fp_kind) /
data ionicXrayFF( 2, 91) / t_ionicFF('Pa',  91,  0, &
     [9.100000e+01_fp_kind,8.912750e+01_fp_kind,8.502510e+01_fp_kind,8.051980e+01_fp_kind,7.608760e+01_fp_kind,&
     7.181680e+01_fp_kind,6.780710e+01_fp_kind,6.412150e+01_fp_kind,6.075910e+01_fp_kind,5.768090e+01_fp_kind,&
     5.483960e+01_fp_kind,4.972680e+01_fp_kind,4.525210e+01_fp_kind,4.134580e+01_fp_kind,3.794240e+01_fp_kind,&
     3.495680e+01_fp_kind,2.990550e+01_fp_kind,2.572680e+01_fp_kind,2.227390e+01_fp_kind,1.952540e+01_fp_kind,&
     1.742790e+01_fp_kind,1.418090e+01_fp_kind,1.214500e+01_fp_kind,1.031970e+01_fp_kind,8.631800e+00_fp_kind,&
     6.252600e+00_fp_kind,5.054900e+00_fp_kind], 2.1173563123435808e+00_fp_kind,1.964040e+01_fp_kind) /
data ionicXrayFF( 1, 92) / t_ionicFF('U ',  92,  0, &
     [9.200000e+01_fp_kind,9.015460e+01_fp_kind,8.606670e+01_fp_kind,8.152180e+01_fp_kind,7.702150e+01_fp_kind,&
     7.267190e+01_fp_kind,6.858170e+01_fp_kind,6.482210e+01_fp_kind,6.139960e+01_fp_kind,5.827840e+01_fp_kind,&
     5.541010e+01_fp_kind,5.027620e+01_fp_kind,4.579640e+01_fp_kind,4.188260e+01_fp_kind,3.846680e+01_fp_kind,&
     3.546820e+01_fp_kind,3.039730e+01_fp_kind,2.619700e+01_fp_kind,2.270110e+01_fp_kind,1.988730e+01_fp_kind,&
     1.771640e+01_fp_kind,1.434490e+01_fp_kind,1.229930e+01_fp_kind,1.050190e+01_fp_kind,8.820800e+00_fp_kind,&
     6.385600e+00_fp_kind,5.133300e+00_fp_kind], 2.1273836809732600e+00_fp_kind,1.929090e+01_fp_kind) /
data ionicXrayFF( 2, 92) / t_ionicFF('U ',  92,  0, &
     [9.200000e+01_fp_kind,9.015460e+01_fp_kind,8.606670e+01_fp_kind,8.152180e+01_fp_kind,7.702150e+01_fp_kind,&
     7.267190e+01_fp_kind,6.858170e+01_fp_kind,6.482210e+01_fp_kind,6.139960e+01_fp_kind,5.827840e+01_fp_kind,&
     5.541010e+01_fp_kind,5.027620e+01_fp_kind,4.579640e+01_fp_kind,4.188260e+01_fp_kind,3.846680e+01_fp_kind,&
     3.546820e+01_fp_kind,3.039730e+01_fp_kind,2.619700e+01_fp_kind,2.270110e+01_fp_kind,1.988730e+01_fp_kind,&
     1.771640e+01_fp_kind,1.434490e+01_fp_kind,1.229930e+01_fp_kind,1.050190e+01_fp_kind,8.820800e+00_fp_kind,&
     6.385600e+00_fp_kind,5.133300e+00_fp_kind], 2.1273836809732600e+00_fp_kind,1.929090e+01_fp_kind) /



type :: t_ionicFF_peng
  character(2) :: element
  integer      :: Z,dZ
  real(fp_kind)    :: a(5),b(5)
end type t_ionicFF_peng

type(t_ionicFF_peng):: ionicFF_Peng(113)

data ionicFF_Peng(1  )/t_ionicFF_peng('H ', 1 ,  -1, &
    [0.140E+0_fp_kind, 0.649E+0_fp_kind, 0.137E+1_fp_kind,  0.337E+0_fp_kind,  0.787E+0_fp_kind]   ,  &
    [0.984E+0_fp_kind,      0.867E+1_fp_kind,   0.389E+2_fp_kind,   0.111E+3_fp_kind,   0.166E+3_fp_kind])/
data ionicFF_Peng(2  )/t_ionicFF_peng('Li', 3 ,  +1, &
    [0.460E-2_fp_kind, 0.165E-1_fp_kind, 0.435E-1_fp_kind,  0.649E-1_fp_kind,  0.270E-1_fp_kind]   ,  &
    [0.358E-1_fp_kind,      0.239E+0_fp_kind,   0.879E+0_fp_kind,   0.264E+1_fp_kind,   0.709E+1_fp_kind])/
data ionicFF_Peng(3  )/t_ionicFF_peng('Be', 4 ,  +2, &
    [0.340E-2_fp_kind, 0.103E-1_fp_kind, 0.233E-1_fp_kind,  0.325E-1_fp_kind,  0.120E-1_fp_kind]   ,  &
    [0.267E-1_fp_kind,      0.162E+0_fp_kind,   0.531E+0_fp_kind,   0.148E+1_fp_kind,   0.388E+1_fp_kind])/
data ionicFF_Peng(4  )/t_ionicFF_peng('O ', 8 ,  -1, &
    [0.205E+0_fp_kind, 0.628E+0_fp_kind, 0.117E+1_fp_kind,  0.103E+1_fp_kind,  0.290E+0_fp_kind]   ,  &
    [0.397E+0_fp_kind,      0.264E+1_fp_kind,   0.880E+1_fp_kind,   0.271E+2_fp_kind,   0.918E+2_fp_kind])/
data ionicFF_Peng(5  )/t_ionicFF_peng('O ', 8 ,  -2, &
    [0.421E-1_fp_kind, 0.210E+0_fp_kind, 0.852E+0_fp_kind,  0.182E+1_fp_kind,  0.117E+1_fp_kind]   ,  &
    [0.609E-1_fp_kind,      0.559E+0_fp_kind,   0.296E+1_fp_kind,   0.115E+2_fp_kind,   0.377E+2_fp_kind])/
data ionicFF_Peng(6  )/t_ionicFF_peng('F ', 9 ,  -1, &
    [0.134E+0_fp_kind, 0.391E+0_fp_kind, 0.814E+0_fp_kind,  0.928E+0_fp_kind,  0.347E+0_fp_kind]   ,  &
    [0.228E+0_fp_kind,      0.147E+1_fp_kind,   0.468E+1_fp_kind,   0.132E+2_fp_kind,   0.360E+2_fp_kind])/
data ionicFF_Peng(7  )/t_ionicFF_peng('Na', 11,  +1, &
    [0.256E-1_fp_kind, 0.919E-1_fp_kind, 0.297E+0_fp_kind,  0.514E+0_fp_kind,  0.199E+0_fp_kind]   ,  &
    [0.397E-1_fp_kind,      0.287E+0_fp_kind,   0.118E+1_fp_kind,   0.375E+1_fp_kind,   0.108E+2_fp_kind])/
data ionicFF_Peng(8  )/t_ionicFF_peng('Mg', 12,  +2, &
    [0.210E-1_fp_kind, 0.672E-1_fp_kind, 0.198E+0_fp_kind,  0.368E+0_fp_kind,  0.174E+0_fp_kind]   ,  &
    [0.331E-1_fp_kind,      0.222E+0_fp_kind,   0.838E+0_fp_kind,   0.248E+1_fp_kind,   0.675E+1_fp_kind])/
data ionicFF_Peng(9  )/t_ionicFF_peng('Al', 13,  +3, &
    [0.192E-1_fp_kind, 0.579E-1_fp_kind, 0.163E+0_fp_kind,  0.284E+0_fp_kind,  0.114E+0_fp_kind]   ,  &
    [0.306E-1_fp_kind,      0.198E+0_fp_kind,   0.713E+0_fp_kind,   0.204E+1_fp_kind,   0.525E+1_fp_kind])/
data ionicFF_Peng(10 )/t_ionicFF_peng('Si', 14,  +4, &
    [0.192E+0_fp_kind, 0.289E+0_fp_kind, 0.100E+0_fp_kind, -0.728E-1_fp_kind,  0.120E-2_fp_kind]   ,  &
    [0.359E+0_fp_kind,      0.196E+1_fp_kind,   0.934E+1_fp_kind,   0.111E+2_fp_kind,   0.134E+2_fp_kind])/
data ionicFF_Peng(11 )/t_ionicFF_peng('Cl', 17,  -1, &
    [0.265E+0_fp_kind, 0.596E+0_fp_kind, 0.160E+1_fp_kind,  0.269E+1_fp_kind,  0.123E+1_fp_kind]   ,  &
    [0.252E+0_fp_kind,      0.156E+1_fp_kind,   0.621E+1_fp_kind,   0.178E+2_fp_kind,   0.478E+2_fp_kind])/
data ionicFF_Peng(12 )/t_ionicFF_peng('K ', 19,  +1, &
    [0.199E+0_fp_kind, 0.396E+0_fp_kind, 0.928E+0_fp_kind,  0.145E+1_fp_kind,  0.450E+0_fp_kind]   ,  &
    [0.192E+0_fp_kind,      0.110E+1_fp_kind,   0.391E+1_fp_kind,   0.975E+1_fp_kind,   0.234E+2_fp_kind])/
data ionicFF_Peng(13 )/t_ionicFF_peng('Ca', 20,  +2, &
    [0.164E+0_fp_kind, 0.327E+0_fp_kind, 0.743E+0_fp_kind,  0.116E+1_fp_kind,  0.307E+0_fp_kind]   ,  &
    [0.157E+0_fp_kind,      0.894E+0_fp_kind,   0.315E+1_fp_kind,   0.767E+1_fp_kind,   0.177E+2_fp_kind])/
data ionicFF_Peng(14 )/t_ionicFF_peng('Sc', 21,  +3, &
    [0.163E+0_fp_kind, 0.307E+0_fp_kind, 0.716E+0_fp_kind,  0.880E+0_fp_kind,  0.139E+0_fp_kind]   ,  &
    [0.157E+0_fp_kind,      0.899E+0_fp_kind,   0.306E+1_fp_kind,   0.705E+1_fp_kind,   0.161E+2_fp_kind])/
data ionicFF_Peng(15 )/t_ionicFF_peng('Ti', 22,  +2, &
    [0.399E+0_fp_kind, 0.104E+1_fp_kind, 0.121E+1_fp_kind, -0.797E-1_fp_kind,  0.352E+0_fp_kind]   ,  &
    [0.376E+0_fp_kind,      0.274E+1_fp_kind,   0.810E+1_fp_kind,   0.142E+2_fp_kind,   0.232E+2_fp_kind])/
data ionicFF_Peng(16 )/t_ionicFF_peng('Ti', 22,  +3, &
    [0.364E+0_fp_kind, 0.919E+0_fp_kind, 0.135E+1_fp_kind, -0.933E+0_fp_kind,  0.589E+0_fp_kind]   ,  &
    [0.364E+0_fp_kind,      0.267E+1_fp_kind,   0.818E+1_fp_kind,   0.118E+2_fp_kind,   0.149E+2_fp_kind])/
data ionicFF_Peng(17 )/t_ionicFF_peng('Ti', 22,  +4, &
    [0.116E+0_fp_kind, 0.256E+0_fp_kind, 0.565E+0_fp_kind,  0.772E+0_fp_kind,  0.132E+0_fp_kind]   ,  &
    [0.108E+0_fp_kind,      0.655E+0_fp_kind,   0.238E+1_fp_kind,   0.551E+1_fp_kind,   0.123E+2_fp_kind])/
data ionicFF_Peng(18 )/t_ionicFF_peng('V ', 23,  +2, &
    [0.317E+0_fp_kind, 0.939E+0_fp_kind, 0.149E+1_fp_kind, -0.131E+1_fp_kind,  0.147E+1_fp_kind]   ,  &
    [0.269E+0_fp_kind,      0.209E+1_fp_kind,   0.722E+1_fp_kind,   0.152E+2_fp_kind,   0.176E+2_fp_kind])/
data ionicFF_Peng(19 )/t_ionicFF_peng('V ', 23,  +3, &
    [0.341E+0_fp_kind, 0.805E+0_fp_kind, 0.942E+0_fp_kind,  0.783E-1_fp_kind,  0.156E+0_fp_kind]   ,  &
    [0.321E+0_fp_kind,      0.223E+1_fp_kind,   0.599E+1_fp_kind,   0.134E+2_fp_kind,   0.169E+2_fp_kind])/
data ionicFF_Peng(20 )/t_ionicFF_peng('V ', 23,  +5, &
    [0.367E-1_fp_kind, 0.124E+0_fp_kind, 0.244E+0_fp_kind,  0.723E+0_fp_kind,  0.435E+0_fp_kind]   ,  &
    [0.330E-1_fp_kind,      0.222E+0_fp_kind,   0.824E+0_fp_kind,   0.280E+1_fp_kind,   0.670E+1_fp_kind])/
data ionicFF_Peng(21 )/t_ionicFF_peng('Cr', 24,  +2, &
    [0.237E+0_fp_kind, 0.634E+0_fp_kind, 0.123E+1_fp_kind,  0.713E+0_fp_kind,  0.859E-1_fp_kind]   ,  &
    [0.177E+0_fp_kind,      0.135E+1_fp_kind,   0.430E+1_fp_kind,   0.122E+2_fp_kind,   0.390E+2_fp_kind])/
data ionicFF_Peng(22 )/t_ionicFF_peng('Cr', 24,  +3, &
    [0.393E+0_fp_kind, 0.105E+1_fp_kind, 0.162E+1_fp_kind, -0.115E+1_fp_kind,  0.407E+0_fp_kind]   ,  &
    [0.359E+0_fp_kind,      0.257E+1_fp_kind,   0.868E+1_fp_kind,   0.110E+2_fp_kind,   0.158E+2_fp_kind])/
data ionicFF_Peng(23 )/t_ionicFF_peng('Cr', 24,  +4, &
    [0.132E+0_fp_kind, 0.292E+0_fp_kind, 0.703E+0_fp_kind,  0.692E+0_fp_kind,  0.959E-1_fp_kind]   ,  &
    [0.109E+0_fp_kind,      0.695E+0_fp_kind,   0.239E+1_fp_kind,   0.565E+1_fp_kind,   0.147E+2_fp_kind])/
data ionicFF_Peng(24 )/t_ionicFF_peng('Mn', 25,  +2, &
    [0.576E-1_fp_kind, 0.210E+0_fp_kind, 0.604E+0_fp_kind,  0.132E+1_fp_kind,  0.659E+0_fp_kind]   ,  &
    [0.398E-1_fp_kind,      0.284E+0_fp_kind,   0.129E+1_fp_kind,   0.423E+1_fp_kind,   0.145E+2_fp_kind])/
data ionicFF_Peng(25 )/t_ionicFF_peng('Mn', 25,  +3, &
    [0.116E+0_fp_kind, 0.523E+0_fp_kind, 0.881E+0_fp_kind,  0.589E+0_fp_kind,  0.214E+0_fp_kind]   ,  &
    [0.117E-1_fp_kind,      0.876E+0_fp_kind,   0.306E+1_fp_kind,   0.644E+1_fp_kind,   0.143E+2_fp_kind])/
data ionicFF_Peng(26 )/t_ionicFF_peng('Mn', 25,  +4, &
    [0.381E+0_fp_kind, 0.183E+1_fp_kind,-0.133E+1_fp_kind,  0.995E+0_fp_kind,  0.618E-1_fp_kind]   ,  &
    [0.354E+0_fp_kind,      0.272E+1_fp_kind,   0.347E+1_fp_kind,   0.547E+1_fp_kind,   0.161E+2_fp_kind])/
data ionicFF_Peng(27 )/t_ionicFF_peng('Fe', 26,  +2, &
    [0.307E+0_fp_kind, 0.838E+0_fp_kind, 0.111E+1_fp_kind,  0.280E+0_fp_kind,  0.277E+0_fp_kind]   ,  &
    [0.230E+0_fp_kind,      0.162E+1_fp_kind,   0.487E+1_fp_kind,   0.107E+2_fp_kind,   0.192E+2_fp_kind])/
data ionicFF_Peng(28 )/t_ionicFF_peng('Fe', 26,  +3, &
    [0.198E+0_fp_kind, 0.387E+0_fp_kind, 0.889E+0_fp_kind,  0.709E+0_fp_kind,  0.117E+0_fp_kind]   ,  &
    [0.154E+0_fp_kind,      0.893E+0_fp_kind,   0.262E+1_fp_kind,   0.665E+1_fp_kind,   0.180E+2_fp_kind])/
data ionicFF_Peng(29 )/t_ionicFF_peng('Co', 27,  +2, &
    [0.213E+0_fp_kind, 0.488E+0_fp_kind, 0.998E+0_fp_kind,  0.828E+0_fp_kind,  0.230E+0_fp_kind]   ,  &
    [0.148E+0_fp_kind,      0.939E+0_fp_kind,   0.278E+1_fp_kind,   0.731E+1_fp_kind,   0.207E+2_fp_kind])/
data ionicFF_Peng(30 )/t_ionicFF_peng('Co', 27,  +3, &
    [0.331E+0_fp_kind, 0.487E+0_fp_kind, 0.729E+0_fp_kind,  0.608E+0_fp_kind,  0.131E+0_fp_kind]   ,  &
    [0.267E+0_fp_kind,      0.141E+1_fp_kind,   0.289E+1_fp_kind,   0.645E+1_fp_kind,   0.158E+2_fp_kind])/
data ionicFF_Peng(31 )/t_ionicFF_peng('Ni', 28,  +2, &
    [0.338E+0_fp_kind, 0.982E+0_fp_kind, 0.132E+1_fp_kind, -0.356E+1_fp_kind,  0.362E+1_fp_kind]   ,  &
    [0.237E+0_fp_kind,      0.167E+1_fp_kind,   0.573E+1_fp_kind,   0.114E+2_fp_kind,   0.121E+2_fp_kind])/
data ionicFF_Peng(32 )/t_ionicFF_peng('Ni', 28,  +3, &
    [0.347E+0_fp_kind, 0.877E+0_fp_kind, 0.790E+0_fp_kind,  0.538E-1_fp_kind,  0.192E+0_fp_kind]   ,  &
    [0.260E+0_fp_kind,      0.171E+1_fp_kind,   0.475E+1_fp_kind,   0.751E+1_fp_kind,   0.130E+2_fp_kind])/
data ionicFF_Peng(33 )/t_ionicFF_peng('Cu', 29,  +1, &
    [0.312E+0_fp_kind, 0.812E+0_fp_kind, 0.111E+1_fp_kind,  0.794E+0_fp_kind,  0.257E+0_fp_kind]   ,  &
    [0.201E+0_fp_kind,      0.131E+1_fp_kind,   0.380E+1_fp_kind,   0.105E+2_fp_kind,   0.282E+2_fp_kind])/
data ionicFF_Peng(34 )/t_ionicFF_peng('Cu', 29,  +2, &
    [0.224E+0_fp_kind, 0.544E+0_fp_kind, 0.970E+0_fp_kind,  0.727E+0_fp_kind,  0.182E+0_fp_kind]   ,  &
    [0.145E+0_fp_kind,      0.933E+0_fp_kind,   0.269E+1_fp_kind,   0.711E+1_fp_kind,   0.194E+2_fp_kind])/
data ionicFF_Peng(35 )/t_ionicFF_peng('Zn', 30,  +2, &
    [0.252E+0_fp_kind, 0.600E+0_fp_kind, 0.917E+0_fp_kind,  0.663E+0_fp_kind,  0.161E+0_fp_kind]   ,  &
    [0.161E+0_fp_kind,      0.101E+1_fp_kind,   0.276E+1_fp_kind,   0.708E+1_fp_kind,   0.190E+2_fp_kind])/
data ionicFF_Peng(36 )/t_ionicFF_peng('Ga', 31,  +3, &
    [0.391E+0_fp_kind, 0.947E+0_fp_kind, 0.690E+0_fp_kind,  0.709E-1_fp_kind,  0.653E-1_fp_kind]   ,  &
    [0.264E+0_fp_kind,      0.165E+1_fp_kind,   0.482E+1_fp_kind,   0.107E+2_fp_kind,   0.152E+2_fp_kind])/
data ionicFF_Peng(37 )/t_ionicFF_peng('Ge', 32,  +4, &
    [0.346E+0_fp_kind, 0.830E+0_fp_kind, 0.599E+0_fp_kind,  0.949E-1_fp_kind, -0.217E-1_fp_kind]   ,  &
    [0.232E+0_fp_kind,      0.145E+1_fp_kind,   0.408E+1_fp_kind,   0.132E+2_fp_kind,   0.295E+2_fp_kind])/
data ionicFF_Peng(38 )/t_ionicFF_peng('Br', 35,  -1, &
    [0.125E+0_fp_kind, 0.563E+0_fp_kind, 0.143E+1_fp_kind,  0.352E+1_fp_kind,  0.322E+1_fp_kind]   ,  &
    [0.530E-1_fp_kind,      0.469E+0_fp_kind,   0.215E+1_fp_kind,   0.111E+2_fp_kind,   0.389E+2_fp_kind])/
data ionicFF_Peng(39 )/t_ionicFF_peng('Br', 37,  +1, &
    [0.368E+0_fp_kind, 0.884E+0_fp_kind, 0.114E+1_fp_kind,  0.226E+1_fp_kind,  0.881E+0_fp_kind]   ,  &
    [0.187E+0_fp_kind,      0.112E+1_fp_kind,   0.398E+1_fp_kind,   0.109E+2_fp_kind,   0.266E+2_fp_kind])/
data ionicFF_Peng(40 )/t_ionicFF_peng('Sr', 38,  +2, &
    [0.346E+0_fp_kind, 0.804E+0_fp_kind, 0.988E+0_fp_kind,  0.189E+1_fp_kind,  0.609E+0_fp_kind]   ,  &
    [0.176E+0_fp_kind,      0.104E+1_fp_kind,   0.359E+1_fp_kind,   0.932E+1_fp_kind,   0.214E+2_fp_kind])/
data ionicFF_Peng(41 )/t_ionicFF_peng('Y ', 39,  +3, &
    [0.465E+0_fp_kind, 0.923E+0_fp_kind, 0.241E+1_fp_kind, -0.231E+1_fp_kind,  0.248E+1_fp_kind]   ,  &
    [0.240E+0_fp_kind,      0.143E+1_fp_kind,   0.645E+1_fp_kind,   0.997E+1_fp_kind,   0.122E+2_fp_kind])/
data ionicFF_Peng(42 )/t_ionicFF_peng('Zr', 40,  +4, &
    [0.234E+0_fp_kind, 0.642E+0_fp_kind, 0.747E+0_fp_kind,  0.147E+1_fp_kind,  0.377E+0_fp_kind]   ,  &
    [0.113E+0_fp_kind,      0.736E+0_fp_kind,   0.254E+1_fp_kind,   0.672E+1_fp_kind,   0.147E+2_fp_kind])/
data ionicFF_Peng(43 )/t_ionicFF_peng('Nb', 41,  +3, &
    [0.377E+0_fp_kind, 0.749E+0_fp_kind, 0.129E+1_fp_kind,  0.161E+1_fp_kind,  0.481E+0_fp_kind]   ,  &
    [0.184E+0_fp_kind,      0.102E+1_fp_kind,   0.380E+1_fp_kind,   0.944E+1_fp_kind,   0.257E+2_fp_kind])/
data ionicFF_Peng(44 )/t_ionicFF_peng('Nb', 41,  +5, &
    [0.828E-1_fp_kind, 0.271E+0_fp_kind, 0.654E+0_fp_kind,  0.124E+1_fp_kind,  0.829E+0_fp_kind]   ,  &
    [0.369E-1_fp_kind,      0.261E+0_fp_kind,   0.957E+0_fp_kind,   0.394E+1_fp_kind,   0.944E+1_fp_kind])/
data ionicFF_Peng(45 )/t_ionicFF_peng('Mo', 42,  +3, &
    [0.401E+0_fp_kind, 0.756E+0_fp_kind, 0.138E+1_fp_kind,  0.158E+1_fp_kind,  0.497E+0_fp_kind]   ,  &
    [0.191E+0_fp_kind,      0.106E+1_fp_kind,   0.384E+1_fp_kind,   0.938E+1_fp_kind,   0.246E+2_fp_kind])/
data ionicFF_Peng(46 )/t_ionicFF_peng('Mo', 42,  +5, &
    [0.479E+0_fp_kind, 0.846E+0_fp_kind, 0.156E+2_fp_kind, -0.152E+2_fp_kind,  0.160E+1_fp_kind]   ,  &
    [0.241E+0_fp_kind,      0.146E+1_fp_kind,   0.679E+1_fp_kind,   0.713E+1_fp_kind,   0.104E+2_fp_kind])/
data ionicFF_Peng(47 )/t_ionicFF_peng('Mo', 42,  +6, &
    [0.203E+0_fp_kind, 0.567E+0_fp_kind, 0.646E+0_fp_kind,  0.116E+1_fp_kind,  0.171E+0_fp_kind]      ,  &
    [0.971E-1_fp_kind,      0.647E+0_fp_kind,   0.228E+1_fp_kind,   0.561E+1_fp_kind,   0.124E+2_fp_kind])/
data ionicFF_Peng(48 )/t_ionicFF_peng('Ru', 44,  +3, &
    [0.428E+0_fp_kind, 0.773E+0_fp_kind, 0.155E+1_fp_kind,  0.146E+1_fp_kind,  0.486E+0_fp_kind]      ,  &
    [0.191E+0_fp_kind,      0.109E+1_fp_kind,   0.382E+1_fp_kind,   0.908E+1_fp_kind,   0.217E+2_fp_kind])/
data ionicFF_Peng(49 )/t_ionicFF_peng('Ru', 44,  +4, &
    [0.282E+0_fp_kind, 0.653E+0_fp_kind, 0.114E+1_fp_kind,  0.153E+1_fp_kind,  0.418E+0_fp_kind]      ,  &
    [0.125E+0_fp_kind,      0.753E+0_fp_kind,   0.285E+1_fp_kind,   0.701E+1_fp_kind,   0.175E+2_fp_kind])/
data ionicFF_Peng(50 )/t_ionicFF_peng('Rh', 45,  +3, &
    [0.352E+0_fp_kind, 0.723E+0_fp_kind, 0.150E+1_fp_kind,  0.163E+1_fp_kind,  0.499E+0_fp_kind]   ,  &
    [0.151E+0_fp_kind,      0.878E+0_fp_kind,   0.328E+1_fp_kind,   0.816E+1_fp_kind,   0.207E+2_fp_kind])/
data ionicFF_Peng(51 )/t_ionicFF_peng('Rh', 45,  +4, &
    [0.397E+0_fp_kind, 0.725E+0_fp_kind, 0.151E+1_fp_kind,  0.119E+1_fp_kind,  0.251E+0_fp_kind]   ,  &
    [0.177E+0_fp_kind,      0.101E+1_fp_kind,   0.362E+1_fp_kind,   0.856E+1_fp_kind,   0.189E+2_fp_kind])/
data ionicFF_Peng(52 )/t_ionicFF_peng('Pd', 46,  +2, &
    [0.935E+0_fp_kind, 0.311E+1_fp_kind, 0.246E+2_fp_kind, -0.436E+2_fp_kind,  0.211E+2_fp_kind]   ,  &
    [0.393E+0_fp_kind,      0.406E+1_fp_kind,   0.431E+2_fp_kind,   0.540E+2_fp_kind,   0.698E+2_fp_kind])/
data ionicFF_Peng(53 )/t_ionicFF_peng('Pd', 46,  +4, &
    [0.348E+0_fp_kind, 0.640E+0_fp_kind, 0.122E+1_fp_kind,  0.145E+1_fp_kind,  0.427E+0_fp_kind]   ,  &
    [0.151E+0_fp_kind,      0.832E+0_fp_kind,   0.285E+1_fp_kind,   0.659E+1_fp_kind,   0.156E+2_fp_kind])/
data ionicFF_Peng(54 )/t_ionicFF_peng('Ag', 47,  +1, &
    [0.503E+0_fp_kind, 0.940E+0_fp_kind, 0.217E+1_fp_kind,  0.199E+1_fp_kind,  0.726E+0_fp_kind]   ,  &
    [0.199E+0_fp_kind,      0.119E+1_fp_kind,   0.405E+1_fp_kind,   0.113E+2_fp_kind,   0.324E+2_fp_kind])/
data ionicFF_Peng(55 )/t_ionicFF_peng('Ag', 47,  +2, &
    [0.431E+0_fp_kind, 0.756E+0_fp_kind, 0.172E+1_fp_kind,  0.178E+1_fp_kind,  0.526E+0_fp_kind]   ,  &
    [0.175E+0_fp_kind,      0.979E+0_fp_kind,   0.330E+1_fp_kind,   0.824E+1_fp_kind,   0.214E+2_fp_kind])/
data ionicFF_Peng(56 )/t_ionicFF_peng('Cd', 48,  +2, &
    [0.425E+0_fp_kind, 0.745E+0_fp_kind, 0.173E+1_fp_kind,  0.174E+1_fp_kind,  0.487E+0_fp_kind]   ,  &
    [0.168E+0_fp_kind,      0.944E+0_fp_kind,   0.314E+1_fp_kind,   0.784E+1_fp_kind,   0.204E+2_fp_kind])/
data ionicFF_Peng(57 )/t_ionicFF_peng('In', 49,  +3, &
    [0.417E+0_fp_kind, 0.755E+0_fp_kind, 0.159E+1_fp_kind,  0.136E+1_fp_kind,  0.451E+0_fp_kind]   ,  &
    [0.164E+0_fp_kind,      0.960E+0_fp_kind,   0.308E+1_fp_kind,   0.703E+1_fp_kind,   0.161E+2_fp_kind])/
data ionicFF_Peng(58 )/t_ionicFF_peng('Sn', 50,  +2, &
    [0.797E+0_fp_kind, 0.213E+1_fp_kind, 0.215E+1_fp_kind, -0.164E+1_fp_kind,  0.272E+1_fp_kind]   ,  &
    [0.317E+0_fp_kind,      0.251E+1_fp_kind,   0.904E+1_fp_kind,   0.242E+2_fp_kind,   0.264E+2_fp_kind])/
data ionicFF_Peng(59 )/t_ionicFF_peng('Sn', 50,  +4, &
    [0.261E+0_fp_kind, 0.642E+0_fp_kind, 0.153E+1_fp_kind,  0.136E+1_fp_kind,  0.177E+0_fp_kind]   ,  &
    [0.957E-1_fp_kind,      0.625E+0_fp_kind,   0.251E+1_fp_kind,   0.631E+1_fp_kind,   0.159E+2_fp_kind])/
data ionicFF_Peng(60 )/t_ionicFF_peng('Sb', 51,  +3, &
    [0.552E+0_fp_kind, 0.114E+1_fp_kind, 0.187E+1_fp_kind,  0.136E+1_fp_kind,  0.414E+0_fp_kind]   ,  &
    [0.212E+0_fp_kind,      0.142E+1_fp_kind,   0.421E+1_fp_kind,   0.125E+2_fp_kind,   0.290E+2_fp_kind])/
data ionicFF_Peng(61 )/t_ionicFF_peng('Sb', 51,  +5, &
    [0.377E+0_fp_kind, 0.588E+0_fp_kind, 0.122E+1_fp_kind,  0.118E+1_fp_kind,  0.244E+0_fp_kind]   ,  &
    [0.151E+0_fp_kind,      0.812E+0_fp_kind,   0.240E+1_fp_kind,   0.527E+1_fp_kind,   0.119E+2_fp_kind])/
data ionicFF_Peng(62 )/t_ionicFF_peng('I ', 53,  -1, &
    [0.901E+0_fp_kind, 0.280E+1_fp_kind, 0.561E+1_fp_kind, -0.869E+1_fp_kind,  0.126E+2_fp_kind]   ,  &
    [0.312E+0_fp_kind,      0.259E+1_fp_kind,   0.141E+2_fp_kind,   0.344E+2_fp_kind,   0.395E+2_fp_kind])/
data ionicFF_Peng(63 )/t_ionicFF_peng('Cs', 55,  +1, &
    [0.587E+0_fp_kind, 0.140E+1_fp_kind, 0.187E+1_fp_kind,  0.348E+1_fp_kind,  0.167E+1_fp_kind]   ,  &
    [0.200E+0_fp_kind,      0.138E+1_fp_kind,   0.412E+1_fp_kind,   0.130E+2_fp_kind,   0.318E+2_fp_kind])/
data ionicFF_Peng(64 )/t_ionicFF_peng('Ba', 56,  +2, &
    [0.733E+0_fp_kind, 0.205E+1_fp_kind, 0.230E+2_fp_kind, -0.152E+3_fp_kind,  0.134E+3_fp_kind]   ,  &
    [0.258E+0_fp_kind,      0.196E+1_fp_kind,   0.118E+2_fp_kind,   0.144E+2_fp_kind,   0.149E+2_fp_kind])/
data ionicFF_Peng(65 )/t_ionicFF_peng('La', 57,  +3, &
    [0.493E+0_fp_kind, 0.110E+1_fp_kind, 0.150E+1_fp_kind,  0.270E+1_fp_kind,  0.108E+1_fp_kind]   ,  &
    [0.167E+0_fp_kind,      0.111E+1_fp_kind,   0.311E+1_fp_kind,   0.961E+1_fp_kind,   0.212E+2_fp_kind])/
data ionicFF_Peng(66 )/t_ionicFF_peng('Ce', 58,  +3, &
    [0.560E+0_fp_kind, 0.135E+1_fp_kind, 0.159E+1_fp_kind,  0.263E+1_fp_kind,  0.706E+0_fp_kind]   ,  &
    [0.190E+0_fp_kind,      0.130E+1_fp_kind,   0.393E+1_fp_kind,   0.107E+2_fp_kind,   0.238E+2_fp_kind])/
data ionicFF_Peng(67 )/t_ionicFF_peng('Ce', 58,  +4, &
    [0.483E+0_fp_kind, 0.109E+1_fp_kind, 0.134E+1_fp_kind,  0.245E+1_fp_kind,  0.797E+0_fp_kind]   ,  &
    [0.165E+0_fp_kind,      0.110E+1_fp_kind,   0.302E+1_fp_kind,   0.885E+1_fp_kind,   0.188E+2_fp_kind])/
data ionicFF_Peng(68 )/t_ionicFF_peng('Pr', 59,  +3, &
    [0.663E+0_fp_kind, 0.173E+1_fp_kind, 0.235E+1_fp_kind,  0.351E+0_fp_kind,  0.159E+1_fp_kind]   ,  &
    [0.226E+0_fp_kind,      0.161E+1_fp_kind,   0.633E+1_fp_kind,   0.110E+2_fp_kind,   0.169E+2_fp_kind])/
data ionicFF_Peng(69 )/t_ionicFF_peng('Pr', 59,  +4, &
    [0.521E+0_fp_kind, 0.119E+1_fp_kind, 0.133E+1_fp_kind,  0.236E+1_fp_kind,  0.690E+0_fp_kind]   ,  &
    [0.177E+0_fp_kind,      0.117E+1_fp_kind,   0.328E+1_fp_kind,   0.894E+1_fp_kind,   0.193E+2_fp_kind])/
data ionicFF_Peng(70 )/t_ionicFF_peng('Nd', 60,  +3, &
    [0.501E+0_fp_kind, 0.118E+1_fp_kind, 0.145E+1_fp_kind,  0.253E+1_fp_kind,  0.920E+0_fp_kind]   ,  &
    [0.162E+0_fp_kind,      0.108E+1_fp_kind,   0.306E+1_fp_kind,   0.880E+1_fp_kind,   0.196E+2_fp_kind])/
data ionicFF_Peng(71 )/t_ionicFF_peng('Pm', 61,  +3, &
    [0.496E+0_fp_kind, 0.120E+1_fp_kind, 0.147E+1_fp_kind,  0.243E+1_fp_kind,  0.943E+0_fp_kind]   ,  &
    [0.156E+0_fp_kind,      0.105E+1_fp_kind,   0.307E+1_fp_kind,   0.856E+1_fp_kind,   0.192E+2_fp_kind])/
data ionicFF_Peng(72 )/t_ionicFF_peng('Sm', 62,  +3, &
    [0.518E+0_fp_kind, 0.124E+1_fp_kind, 0.143E+1_fp_kind,  0.240E+1_fp_kind,  0.781E+0_fp_kind]   ,  &
    [0.163E+0_fp_kind,      0.108E+1_fp_kind,   0.311E+1_fp_kind,   0.852E+1_fp_kind,   0.191E+2_fp_kind])/
data ionicFF_Peng(73 )/t_ionicFF_peng('Eu', 63,  +2, &
    [0.613E+0_fp_kind, 0.153E+1_fp_kind, 0.184E+1_fp_kind,  0.246E+1_fp_kind,  0.714E+0_fp_kind]   ,  &
    [0.190E+0_fp_kind,      0.127E+1_fp_kind,   0.418E+1_fp_kind,   0.107E+2_fp_kind,   0.262E+2_fp_kind])/
data ionicFF_Peng(74 )/t_ionicFF_peng('Eu', 63,  +3, &
    [0.496E+0_fp_kind, 0.121E+1_fp_kind, 0.145E+1_fp_kind,  0.236E+1_fp_kind,  0.774E+0_fp_kind]   ,  &
    [0.152E+0_fp_kind,      0.101E+1_fp_kind,   0.295E+1_fp_kind,   0.818E+1_fp_kind,   0.185E+2_fp_kind])/
data ionicFF_Peng(75 )/t_ionicFF_peng('Gd', 64,  +3, &
    [0.490E+0_fp_kind, 0.119E+1_fp_kind, 0.142E+1_fp_kind,  0.230E+1_fp_kind,  0.795E+0_fp_kind]   ,  &
    [0.148E+0_fp_kind,      0.974E+0_fp_kind,   0.281E+1_fp_kind,   0.778E+1_fp_kind,   0.177E+2_fp_kind])/
data ionicFF_Peng(76 )/t_ionicFF_peng('Tb', 65,  +3, &
    [0.503E+0_fp_kind, 0.122E+1_fp_kind, 0.142E+1_fp_kind,  0.224E+1_fp_kind,  0.710E+0_fp_kind]   ,  &
    [0.150E+0_fp_kind,      0.982E+0_fp_kind,   0.286E+1_fp_kind,   0.777E+1_fp_kind,   0.177E+2_fp_kind])/
data ionicFF_Peng(77 )/t_ionicFF_peng('Dy', 66,  +3, &
    [0.503E+0_fp_kind, 0.124E+1_fp_kind, 0.144E+1_fp_kind,  0.217E+1_fp_kind,  0.643E+0_fp_kind]   ,  &
    [0.148E+0_fp_kind,      0.970E+0_fp_kind,   0.288E+1_fp_kind,   0.773E+1_fp_kind,   0.176E+2_fp_kind])/
data ionicFF_Peng(78 )/t_ionicFF_peng('Ho', 67,  +3, &
    [0.456E+0_fp_kind, 0.117E+1_fp_kind, 0.143E+1_fp_kind,  0.215E+1_fp_kind,  0.692E+0_fp_kind]   ,  &
    [0.129E+0_fp_kind,      0.869E+0_fp_kind,   0.261E+1_fp_kind,   0.724E+1_fp_kind,   0.167E+2_fp_kind])/
data ionicFF_Peng(79 )/t_ionicFF_peng('Er', 68,  +3, &
    [0.522E+0_fp_kind, 0.128E+1_fp_kind, 0.146E+1_fp_kind,  0.205E+1_fp_kind,  0.508E+0_fp_kind]   ,  &
    [0.150E+0_fp_kind,      0.964E+0_fp_kind,   0.293E+1_fp_kind,   0.772E+1_fp_kind,   0.178E+2_fp_kind])/
data ionicFF_Peng(80 )/t_ionicFF_peng('Tm', 69,  +3, &
    [0.475E+0_fp_kind, 0.120E+1_fp_kind, 0.142E+1_fp_kind,  0.205E+1_fp_kind,  0.584E+0_fp_kind]   ,  &
    [0.132E+0_fp_kind,      0.864E+0_fp_kind,   0.260E+1_fp_kind,   0.709E+1_fp_kind,   0.166E+2_fp_kind])/
data ionicFF_Peng(81 )/t_ionicFF_peng('Yb', 70,  +2, &
    [0.508E+0_fp_kind, 0.137E+1_fp_kind, 0.176E+1_fp_kind,  0.223E+1_fp_kind,  0.584E+0_fp_kind]   ,  &
    [0.136E+0_fp_kind,      0.922E+0_fp_kind,   0.312E+1_fp_kind,   0.872E+1_fp_kind,   0.237E+2_fp_kind])/
data ionicFF_Peng(82 )/t_ionicFF_peng('Yb', 70,  +3, &
    [0.498E+0_fp_kind, 0.122E+1_fp_kind, 0.139E+1_fp_kind,  0.197E+1_fp_kind,  0.559E+0_fp_kind]   ,  &
    [0.138E+0_fp_kind,      0.881E+0_fp_kind,   0.263E+1_fp_kind,   0.699E+1_fp_kind,   0.163E+2_fp_kind])/
data ionicFF_Peng(83 )/t_ionicFF_peng('Lu', 71,  +3, &
    [0.483E+0_fp_kind, 0.121E+1_fp_kind, 0.141E+1_fp_kind,  0.194E+1_fp_kind,  0.522E+0_fp_kind]   ,  &
    [0.131E+0_fp_kind,      0.845E+0_fp_kind,   0.257E+1_fp_kind,   0.688E+1_fp_kind,   0.162E+2_fp_kind])/
data ionicFF_Peng(84 )/t_ionicFF_peng('Hf', 72,  +4, &
    [0.522E+0_fp_kind, 0.122E+1_fp_kind, 0.137E+1_fp_kind,  0.168E+1_fp_kind,  0.312E+0_fp_kind]   ,  &
    [0.145E+0_fp_kind,      0.896E+0_fp_kind,   0.274E+1_fp_kind,   0.691E+1_fp_kind,   0.161E+2_fp_kind])/
data ionicFF_Peng(85 )/t_ionicFF_peng('Ta', 73,  +5, &
    [0.569E+0_fp_kind, 0.126E+1_fp_kind, 0.979E+0_fp_kind,  0.129E+1_fp_kind,  0.551E+0_fp_kind]   ,  &
    [0.161E+0_fp_kind,      0.972E+0_fp_kind,   0.276E+1_fp_kind,   0.540E+1_fp_kind,   0.109E+2_fp_kind])/
data ionicFF_Peng(86 )/t_ionicFF_peng('W' , 74,  +6, &
    [0.181E+0_fp_kind, 0.873E+0_fp_kind, 0.118E+1_fp_kind,  0.148E+1_fp_kind,  0.562E+0_fp_kind]   ,  &
    [0.118E-1_fp_kind,      0.442E+0_fp_kind,   0.152E+1_fp_kind,   0.435E+1_fp_kind,   0.942E+1_fp_kind])/
data ionicFF_Peng(87 )/t_ionicFF_peng('Os', 76,  +4, &
    [0.586E+0_fp_kind, 0.131E+1_fp_kind, 0.163E+1_fp_kind,  0.171E+1_fp_kind,  0.540E+0_fp_kind]   ,  &
    [0.155E+0_fp_kind,      0.938E+0_fp_kind,   0.319E+1_fp_kind,   0.784E+1_fp_kind,   0.193E+2_fp_kind])/
data ionicFF_Peng(88 )/t_ionicFF_peng('Ir', 77,  +3, &
    [0.692E+0_fp_kind, 0.137E+1_fp_kind, 0.180E+1_fp_kind,  0.197E+1_fp_kind,  0.804E+0_fp_kind]   ,  &
    [0.182E+0_fp_kind,      0.104E+1_fp_kind,   0.347E+1_fp_kind,   0.851E+1_fp_kind,   0.212E+2_fp_kind])/
data ionicFF_Peng(89 )/t_ionicFF_peng('Ir', 77,  +4, &
    [0.653E+0_fp_kind, 0.129E+1_fp_kind, 0.150E+1_fp_kind,  0.174E+1_fp_kind,  0.683E+0_fp_kind]   ,  &
    [0.174E+0_fp_kind,      0.992E+0_fp_kind,   0.314E+1_fp_kind,   0.722E+1_fp_kind,   0.172E+2_fp_kind])/
data ionicFF_Peng(90 )/t_ionicFF_peng('Pt', 78,  +2, &
    [0.872E+0_fp_kind, 0.168E+1_fp_kind, 0.263E+1_fp_kind,  0.193E+1_fp_kind,  0.475E+0_fp_kind]   ,  &
    [0.223E+0_fp_kind,      0.135E+1_fp_kind,   0.499E+1_fp_kind,   0.136E+2_fp_kind,   0.330E+2_fp_kind])/
data ionicFF_Peng(91 )/t_ionicFF_peng('Pt', 78,  +4, &
    [0.550E+0_fp_kind, 0.121E+1_fp_kind, 0.162E+1_fp_kind,  0.195E+1_fp_kind,  0.610E+0_fp_kind]   ,  &
    [0.142E+0_fp_kind,      0.833E+0_fp_kind,   0.281E+1_fp_kind,   0.721E+1_fp_kind,   0.177E+2_fp_kind])/
data ionicFF_Peng(92 )/t_ionicFF_peng('Au', 79,  +1, &
    [0.811E+0_fp_kind, 0.157E+1_fp_kind, 0.263E+1_fp_kind,  0.268E+1_fp_kind,  0.998E+0_fp_kind]   ,  &
    [0.201E+0_fp_kind,      0.118E+1_fp_kind,   0.425E+1_fp_kind,   0.121E+2_fp_kind,   0.344E+2_fp_kind])/
data ionicFF_Peng(93 )/t_ionicFF_peng('Au', 79,  +3, &
    [0.722E+0_fp_kind, 0.139E+1_fp_kind, 0.194E+1_fp_kind,  0.194E+1_fp_kind,  0.699E+0_fp_kind]   ,  &
    [0.184E+0_fp_kind,      0.106E+1_fp_kind,   0.358E+1_fp_kind,   0.856E+1_fp_kind,   0.204E+2_fp_kind])/
data ionicFF_Peng(94 )/t_ionicFF_peng('Hg', 80,  +1, &
    [0.796E+0_fp_kind, 0.156E+1_fp_kind, 0.272E+1_fp_kind,  0.276E+1_fp_kind,  0.118E+1_fp_kind]   ,  &
    [0.194E+0_fp_kind,      0.114E+1_fp_kind,   0.421E+1_fp_kind,   0.124E+2_fp_kind,   0.362E+2_fp_kind])/
data ionicFF_Peng(95 )/t_ionicFF_peng('Hg', 80,  +2, &
    [0.773E+0_fp_kind, 0.149E+1_fp_kind, 0.245E+1_fp_kind,  0.223E+1_fp_kind,  0.570E+0_fp_kind]   ,  &
    [0.191E+0_fp_kind,      0.112E+1_fp_kind,   0.400E+1_fp_kind,   0.108E+2_fp_kind,   0.276E+2_fp_kind])/
data ionicFF_Peng(96 )/t_ionicFF_peng('Tl', 81,  +1, &
    [0.820E+0_fp_kind, 0.157E+1_fp_kind, 0.278E+1_fp_kind,  0.282E+1_fp_kind,  0.131E+1_fp_kind]   ,  &
    [0.197E+0_fp_kind,      0.116E+1_fp_kind,   0.423E+1_fp_kind,   0.127E+2_fp_kind,   0.357E+2_fp_kind])/
data ionicFF_Peng(97 )/t_ionicFF_peng('Tl', 81,  +3, &
    [0.836E+0_fp_kind, 0.143E+1_fp_kind, 0.394E+0_fp_kind,  0.251E+1_fp_kind,  0.150E+1_fp_kind]   ,  &
    [0.208E+0_fp_kind,      0.120E+1_fp_kind,   0.257E+1_fp_kind,   0.486E+1_fp_kind,   0.135E+2_fp_kind])/
data ionicFF_Peng(98 )/t_ionicFF_peng('Pb', 82,  +2, &
    [0.755E+0_fp_kind, 0.144E+1_fp_kind, 0.248E+1_fp_kind,  0.245E+1_fp_kind,  0.103E+1_fp_kind]   ,  &
    [0.181E+0_fp_kind,      0.105E+1_fp_kind,   0.375E+1_fp_kind,   0.106E+2_fp_kind,   0.279E+2_fp_kind])/
data ionicFF_Peng(99 )/t_ionicFF_peng('Pb', 82,  +4, &
    [0.583E+0_fp_kind, 0.114E+1_fp_kind, 0.160E+1_fp_kind,  0.206E+1_fp_kind,  0.662E+0_fp_kind]   ,  &
    [0.144E+0_fp_kind,      0.796E+0_fp_kind,   0.258E+1_fp_kind,   0.622E+1_fp_kind,   0.148E+2_fp_kind])/
data ionicFF_Peng(100)/t_ionicFF_peng('Bi', 83,  +3, &
    [0.708E+0_fp_kind, 0.135E+1_fp_kind, 0.228E+1_fp_kind,  0.218E+1_fp_kind,  0.797E+0_fp_kind]   ,  &
    [0.170E+0_fp_kind,      0.981E+0_fp_kind,   0.344E+1_fp_kind,   0.941E+1_fp_kind,   0.237E+2_fp_kind])/
data ionicFF_Peng(101)/t_ionicFF_peng('Bi', 83,  +5, &
    [0.654E+0_fp_kind, 0.118E+1_fp_kind, 0.125E+1_fp_kind,  0.166E+1_fp_kind,  0.778E+0_fp_kind]   ,  &
    [0.162E+0_fp_kind,      0.905E+0_fp_kind,   0.268E+1_fp_kind,   0.514E+1_fp_kind,   0.112E+2_fp_kind])/
data ionicFF_Peng(102)/t_ionicFF_peng('Ra', 88,  +2, &
    [0.911E+0_fp_kind, 0.165E+1_fp_kind, 0.253E+1_fp_kind,  0.362E+1_fp_kind,  0.158E+1_fp_kind]   ,  &
    [0.204E+0_fp_kind,      0.126E+1_fp_kind,   0.403E+1_fp_kind,   0.126E+2_fp_kind,   0.300E+2_fp_kind])/
data ionicFF_Peng(103)/t_ionicFF_peng('Ac', 89,  +3, &
    [0.915E+0_fp_kind, 0.164E+1_fp_kind, 0.226E+1_fp_kind,  0.318E+1_fp_kind,  0.125E+1_fp_kind]   ,  &
    [0.205E+0_fp_kind,      0.128E+1_fp_kind,   0.392E+1_fp_kind,   0.113E+2_fp_kind,   0.251E+2_fp_kind])/
data ionicFF_Peng(104)/t_ionicFF_peng('U ', 92,  +3, &
    [0.114E+1_fp_kind, 0.248E+1_fp_kind, 0.361E+1_fp_kind,  0.113E+1_fp_kind,  0.900E+0_fp_kind]   ,  &
    [0.250E+0_fp_kind,      0.184E+1_fp_kind,   0.739E+1_fp_kind,   0.180E+2_fp_kind,   0.227E+2_fp_kind])/
data ionicFF_Peng(105)/t_ionicFF_peng('U ', 92,  +4, &
    [0.109E+1_fp_kind, 0.232E+1_fp_kind, 0.120E+2_fp_kind,  0.911E+1_fp_kind,  0.215E+1_fp_kind]   ,  &
    [0.243E+0_fp_kind,      0.175E+1_fp_kind,   0.779E+1_fp_kind,   0.831E+1_fp_kind,   0.165E+2_fp_kind])/
data ionicFF_Peng(106)/t_ionicFF_peng('U ', 92,  +6, &
    [0.687E+0_fp_kind, 0.114E+1_fp_kind, 0.183E+1_fp_kind,  0.253E+1_fp_kind,  0.957E+0_fp_kind]   ,  &
    [0.154E+0_fp_kind,      0.861E+0_fp_kind,   0.258E+1_fp_kind,   0.770E+1_fp_kind,   0.159E+2_fp_kind])/
    contains




    !--------------------------------------------------------------------------------------
!     This function returns the electron scattering factors in
!     ANGSTROM units. The atom type is given by the index K, and
!     the scattering vector squared (in A ^ (-2)) by S2. This
!     has no relativistic mass correction at this stage. This
!     correction is carried out elsewhere in the program (in CCD
!     or RELM).
!
!     Uses 5 gaussian summation coefficients from Waasmaier and Kirfel
!     modified 4/6/96 for ions
!   elsa_ext() determines the electron scattering factor using the
!   shielded coulomb extrapolation for scattering vectors beyond
!   the limits for which the parameterisation is fitted.
!
!   see Elsa for further descriptions, as they are essentially the
!   same function (AJD)


    function double_elsa_ext(nt,k,atomf,s2)

        implicit none

        integer(4) nt,m,n,k
        real(8) atomf(13,nt),s2,deltak,alpha,total,xbs
        real(8) c1,fxray,s2l
        real(8) double_elsa_ext
        data c1/2.393367d-02/

        ! alpha = inverse yukawa range in A (only used for small s2 in ions)
        alpha=0.02d0

        if(s2.lt.1.0d-03) then
                ! Small scattering vector limit
                total=atomf(11,k)
                do n=1,5
                      total=total+atomf(n,k)
                enddo
                ! if deltak = 0 then we have a neutral atom else deltak = Z - sum a_i - c
                deltak=float(nint(atomf(12,k)-total))


                total = 0.0d0
                do n= 1, 5
                    m = n + 5
                    total = total + atomf(n,k) * atomf(m,k)* (1.0 - atomf(m,k)/2.0d0*s2)
                enddo

                total=total+deltak/(s2+alpha**2d0)
                double_elsa_ext = c1 * total
        else
                if (s2.gt.36.0d0) then
                  ! Large scattering vector limit
                    s2l = s2 + atomf(13,k)
                      double_elsa_ext = c1 * atomf(12,k) / s2l
                  else
                  ! scattering vector in
                  ! parameterisation range
                      fxray = atomf(11,k)
                      do n = 1, 5
                          m = n + 5
                          xbs = - atomf(m,k) * s2
                          xbs = exp(xbs)
                          fxray = fxray + atomf(n,k) * xbs
                      enddo
                      double_elsa_ext = c1 * (atomf(12,k) - fxray) / s2
                endif
        endif

    end function double_elsa_ext



    function single_elsa_ext(nt,k,atomf,s2)

        implicit none

        integer(4) nt,m,n,k
        real(4) atomf(13,nt),s2,deltak,alpha,total,xbs
        real(4) c1,fxray,s2l
        real(4) single_elsa_ext
        data c1/2.393367e-02/

        ! alpha = inverse yukawa range in A (only used for small s2 in ions)
        alpha=0.02

        if(s2.lt.1.0e-03) then
                ! Small scattering vector limit
                total=atomf(11,k)
                do n=1,5
                      total=total+atomf(n,k)
                enddo
                ! if deltak = 0 then we have a neutral atom else deltak = Z - sum a_i - c
                deltak=float(nint(atomf(12,k)-total))


                total = 0.0
                do n= 1, 5
                    m = n + 5
                    total = total + atomf(n,k) * atomf(m,k)* (1.0 - atomf(m,k)/2.0*s2)
                enddo

                total=total+deltak/(s2+alpha**2)
                single_elsa_ext = c1 * total
        else
                if (s2.gt.36.0) then
                  ! Large scattering vector limit
                    s2l = s2 + atomf(13,k)
                      single_elsa_ext = c1 * atomf(12,k) / s2l
                  else
                  ! scattering vector in
                  ! parameterisation range
                      fxray = atomf(11,k)
                      do n = 1, 5
                          m = n + 5
                          xbs = - atomf(m,k) * s2
                          xbs = exp(xbs)
                          fxray = fxray + atomf(n,k) * xbs
                      enddo
                      single_elsa_ext = c1 * (atomf(12,k) - fxray) / s2
                endif
        endif

    end function single_elsa_ext


    function Peng_ionic_FF_integer_ionicity(s2,Z,dZ,cutoff) result(fe)
        use m_precision

        real(fp_kind),intent(in)::s2,cutoff
        integer*4,intent(in):: Z,dZ
        integer*4::Zapprox
        optional:: cutoff

        real(fp_kind)::fe,cutoff_

        integer*4::i

        cutoff_ = 1e-4
        if (present(cutoff)) cutoff_=cutoff

        if (s2<cutoff_) then
            fe = 0
        else
            fe = 0.023934*dZ/s2
        endif

        !First find atom
        do i=1,113
            if (ionicFF_Peng(i)%Z==Z.and.ionicFF_Peng(i)%dZ==dZ) then

            cutoff_ = 1e-4
            if (present(cutoff)) cutoff_=cutoff



            fe = fe+sum(ionicFF_Peng(i)%a(1:5)*exp(-ionicFF_Peng(i)%b(1:5)*s2))
            return
            endif
        enddo
        !write(*,*) 'A parametrization of Z = ',Z,' DeltaZ = ',dZ,'does not seem to be included in Peng (1998).'
        !write(*,*) 'It will be approximated using the ionic scattering factors of Waasmaier and Kirfel (1994).'
        !If no atom, approximate it using the closest Z atom from Waasmeier and Kirfel
        !plus ionic contribtuion

        Zapprox = Z-dZ

        fe = fe + elsa_ext(92,Zapprox,xrayFF(2:14,:),s2)

    end function

    function Peng_ionic_FF_fractional_ionicity(s2,Z,dZ,cutoff) result(fe)
        use m_precision

        real(fp_kind),intent(in)::s2,cutoff,dZ
        integer*4,intent(in):: Z
        integer*4::Zceil,Zfloor
        optional:: cutoff

        real(fp_kind)::fe,frac,cutoff_

        cutoff_ = 1e-4
        if (present(cutoff)) cutoff_=cutoff

        !Find floor and ceiling of fractional ionicity and calculate
        Zceil = ceiling(dZ)
        frac = modulo(dZ,1.0)
        fe = Peng_ionic_FF_integer_ionicity(s2,Z,Zceil,cutoff_)*frac
        Zfloor = floor(dZ)
        fe = fe+Peng_ionic_FF_integer_ionicity(s2,Z,Zfloor,cutoff_)*(1-frac)


    end function

      function wavev(e)
      !  this function returns the wavevector in one over lambda, in a-1,
      !  for an input electron energy e in ev.

      use m_precision

      implicit none

      real(fp_kind) c1,c2,e
      real(fp_kind) wavev
      data c1, c2 / 9.78475598e-07_fp_kind, 12.2642596_fp_kind /

      wavev = sqrt( e + c1 *e ** 2.0_fp_kind ) / c2

      end function



   end module
