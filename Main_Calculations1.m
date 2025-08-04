%Inital function for calculating cation constants from input data

function [T] = Main_Calculations1(T_input)

T = Data_Normalization(T_input);

%Constants From Prohaska, T. et al. (2022) Standard atomic weights of the elements 2021 (IUPAC Technical Report) DOI: 10.1515/pac-2019-0603
atomic_O = 15.999;
atomic_Si = 28.085;     atomic_SiO2 = 60.083;
atomic_Ti = 47.867;     atomic_TiO2 = 79.865;
atomic_Al = 26.982;     atomic_Al2O3 = 101.961;
atomic_Cr = 51.996;     atomic_Cr2O3 = 151.989;
atomic_Fe = 55.845;     atomic_FeO = 71.844;
atomic_Mn = 54.938;     atomic_MnO = 70.937; 
atomic_Mg = 24.305;     atomic_MgO = 40.304;
atomic_Ni = 58.693;     atomic_NiO = 74.692;
atomic_Ca = 40.078;     atomic_CaO = 56.077;
atomic_Na = 22.99;      atomic_Na2O = 61.979;
atomic_K = 39.098;      atomic_K2O = 94.195;
atomic_P = 30.974;      atomic_P2O5 = 141.943;
atomic_S = 32.06;       atomic_SO3 = 80.057;
atomic_Li = 6.94;       atomic_Li2O = 29.879;

atomic_Sr = 87.62;      atomic_SrO = 103.619;
atomic_Y = 88.906;      atomic_Y2O3 = 225.809;
atomic_Zr = 91.224;     atomic_ZrO2 = 123.222;
atomic_Pb = 207.2;      atomic_PbO = 223.199;
atomic_Ag = 107.87;     atomic_Ag2O = 231.739;
atomic_V = 50.942;      atomic_V2O3 = 149.881;
atomic_Ba = 137.33;     atomic_BaO = 153.329;
atomic_Zn = 65.38;      atomic_ZnO = 81.379;

atomic_La = 138.91;     atomic_La2O3 = 325.817;
atomic_Ce = 140.12;     atomic_Ce2O3 = 328.237;
atomic_Pr = 140.91;     atomic_Pr2O3 = 329.817; 
atomic_Nd = 144.24;		atomic_Nd2O3 = 336.477;
atomic_Pm = 145;
atomic_Sm = 150.36;     atomic_Sm2O3 = 348.717;
atomic_Eu = 151.96;     atomic_EuO = 167.959;
atomic_Gd = 157.25;     atomic_Gd2O3 = 362.497;
atomic_Tb = 158.93;     atomic_Tb2O3 = 365.857;
atomic_Dy = 162.5;      atomic_Dy2O3 = 372.997;
atomic_Ho = 164.93;     atomic_Ho2O3 = 377.857;
atomic_Er = 167.26;     atomic_Er2O3 = 382.517;
atomic_Tm = 168.93;     atomic_Tm2O3 = 385.857;
atomic_Yb = 173.05;     atomic_Yb2O3 = 394.097;
atomic_Lu = 174.97;     atomic_Lu2O3 = 397.937;
atomic_B  = 10.81;      atomic_B2O3 = 69.617;
atomic_Sb = 121.76;     atomic_Sb2O3 = 291.517;
atomic_As = 74.922;     atomic_As2O5 = 229.839;
atomic_Cu = 63.546;     atomic_CuO = 79.545;
atomic_Th = 232;        atomic_ThO2 = 263.998;
atomic_Hf = 178.49;     atomic_HfO2 = 210.488;
atomic_Co = 58.933;     atomic_CoO = 74.932;

atomic_F = 18.998;
atomic_Cl = 35.45;

nor_Oxygen = 24; % The normalized Oxygen atomic number is 24.

%Calculating mole #
T.SiO2_Mole = T.SiO2 ./ atomic_SiO2;
T.TiO2_Mole = T.TiO2 ./ atomic_TiO2;
T.Al2O3_Mole = T.Al2O3 ./ atomic_Al2O3;
T.Cr2O3_Mole = T.Cr2O3 ./ atomic_Cr2O3;
T.FeO_Mole = T.FeO ./ atomic_FeO;
T.MnO_Mole = T.MnO ./ atomic_MnO;
T.MgO_Mole = T.MgO ./ atomic_MgO;
T.NiO_Mole = T.NiO ./ atomic_NiO;
T.CaO_Mole = T.CaO ./ atomic_CaO;
T.Na2O_Mole = T.Na2O ./ atomic_Na2O;
T.K2O_Mole = T.K2O ./ atomic_K2O;
T.SO3_Mole = T.SO3 ./ atomic_SO3;
T.P2O5_Mole = T.P2O5 ./ atomic_P2O5;

T.V2O3_Mole = T.V2O3 ./ atomic_V2O3;
T.ZnO_Mole = T.ZnO ./ atomic_ZnO;
T.CoO_Mole = T.CoO ./ atomic_CoO;
T.BaO_Mole = T.BaO ./ atomic_BaO;
T.SrO_Mole = T.SrO ./ atomic_SrO;
T.B2O3_Mole = T.B2O3 ./ atomic_B2O3;
T.PbO_Mole = T.PbO ./ atomic_PbO;
T.CuO_Mole = T.CuO ./ atomic_CuO;
T.Sb2O3_Mole = T.Sb2O3 ./ atomic_Sb2O3;
T.As2O5_Mole = T.As2O5 ./ atomic_As2O5;
T.ThO2_Mole = T.ThO2 ./ atomic_ThO2;
T.ZrO2_Mole = T.ZrO2 ./ atomic_ZrO2;
T.HfO2_Mole = T.HfO2 ./ atomic_HfO2;
T.Ag2O_Mole = T.Ag2O ./ atomic_Ag2O;
T.Y2O3_Mole = T.Y2O3 ./ atomic_Y2O3;
T.La2O3_Mole = T.La2O3 ./ atomic_La2O3;
T.Ce2O3_Mole = T.Ce2O3 ./ atomic_Ce2O3;
T.Nd2O3_Mole = T.Nd2O3 ./ atomic_Nd2O3;
T.Cl_Mole = T.Cl ./ atomic_Cl;
T.F_Mole = T.F ./ atomic_F;

%Calculating oxygens for each oxide
T.SiO2_Oxygen = T.SiO2_Mole * 2;
T.TiO2_Oxygen = T.TiO2_Mole * 2;
T.Al2O3_Oxygen = T.Al2O3_Mole * 3;
T.Cr2O3_Oxygen = T.Cr2O3_Mole * 3;
T.FeO_Oxygen = T.FeO_Mole * 1;
T.MnO_Oxygen = T.MnO_Mole * 1;
T.MgO_Oxygen = T.MgO_Mole * 1;
T.NiO_Oxygen = T.NiO_Mole * 1;
T.CaO_Oxygen = T.CaO_Mole * 1;
T.Na2O_Oxygen = T.Na2O_Mole * 1;
T.K2O_Oxygen = T.K2O_Mole * 1;
T.SO3_Oxygen = T.SO3_Mole * 3;
T.P2O5_Oxygen = T.P2O5_Mole * 5;

T.V2O3_Oxygen = T.V2O3_Mole * 3;
T.ZnO_Oxygen = T.ZnO_Mole * 1;
T.CoO_Oxygen = T.CoO_Mole * 1;
T.BaO_Oxygen = T.BaO_Mole * 1;
T.SrO_Oxygen = T.SrO_Mole * 1;
T.B2O3_Oxygen = T.B2O3_Mole * 3;
T.PbO_Oxygen = T.PbO_Mole * 1;
T.CuO_Oxygen = T.CuO_Mole * 1;
T.Sb2O3_Oxygen = T.Sb2O3_Mole * 3;
T.As2O5_Oxygen = T.As2O5_Mole * 5;
T.ThO2_Oxygen = T.ThO2_Mole * 2;
T.ZrO2_Oxygen = T.ZrO2_Mole * 2;
T.HfO2_Oxygen = T.HfO2_Mole * 2;
T.Ag2O_Oxygen = T.Ag2O_Mole * 2;
T.Y2O3_Oxygen = T.Y2O3_Mole * 3;
T.La2O3_Oxygen = T.La2O3_Mole * 3;
T.Ce2O3_Oxygen = T.Ce2O3_Mole * 3;
T.Nd2O3_Oxygen = T.Nd2O3_Mole * 3;
T.Cl_Oxygen = T.Cl_Mole * 1;
T.F_Oxygen = T.F_Mole * 1;

%Sum of all oxygen
T.Oxygen_sum = T.SiO2_Oxygen...
    + T.TiO2_Oxygen...
    + T.Al2O3_Oxygen...
    + T.Cr2O3_Oxygen...
    + T.FeO_Oxygen...
    + T.MnO_Oxygen...
    + T.MgO_Oxygen...
    + T.NiO_Oxygen...
    + T.CaO_Oxygen...
    + T.Na2O_Oxygen...
    + T.K2O_Oxygen...
    + T.SO3_Oxygen...
    + T.P2O5_Oxygen...
    + T.V2O3_Oxygen...
    + T.ZnO_Oxygen...
    + T.CoO_Oxygen...
    + T.BaO_Oxygen...
    + T.SrO_Oxygen...
    + T.B2O3_Oxygen...
    + T.PbO_Oxygen...
    + T.CuO_Oxygen...
    + T.Sb2O3_Oxygen...
    + T.As2O5_Oxygen...
    + T.ThO2_Oxygen...
    + T.ZrO2_Oxygen...
    + T.HfO2_Oxygen...
    + T.Ag2O_Oxygen...
    + T.Y2O3_Oxygen...
    + T.La2O3_Oxygen...
    + T.Ce2O3_Oxygen...
    + T.Nd2O3_Oxygen;

%Calculating normalized cations
T.Si_Normalized = T.SiO2_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Ti_Normalized = T.TiO2_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Al_Normalized = T.Al2O3_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Cr_Normalized = T.Cr2O3_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Fe_Normalized = T.FeO_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Mn_Normalized = T.MnO_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Mg_Normalized = T.MgO_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Ni_Normalized = T.NiO_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Ca_Normalized = T.CaO_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Na_Normalized = T.Na2O_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.K_Normalized  = T.K2O_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.S_Normalized  = T.SO3_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.P_Normalized  = T.P2O5_Oxygen * nor_Oxygen ./ T.Oxygen_sum;

T.V_Normalized  = T.V2O3_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Zn_Normalized  = T.ZnO_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Co_Normalized  = T.CoO_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Ba_Normalized  = T.BaO_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Sr_Normalized  = T.SrO_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.B_Normalized  = T.B2O3_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Pb_Normalized  = T.PbO_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Cu_Normalized  = T.CuO_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Sb_Normalized  = T.Sb2O3_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.As_Normalized  = T.As2O5_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Th_Normalized  = T.ThO2_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Zr_Normalized  = T.ZrO2_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Hf_Normalized  = T.HfO2_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Ag_Normalized  = T.Ag2O_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Y_Normalized  = T.Y2O3_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.La_Normalized  = T.La2O3_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Ce_Normalized  = T.Ce2O3_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Nd_Normalized  = T.Nd2O3_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.Cl_Normalized  = T.Cl_Oxygen * nor_Oxygen ./ T.Oxygen_sum;
T.F_Normalized  = T.F_Oxygen * nor_Oxygen ./ T.Oxygen_sum;

T.Si_Multi = T.Si_Normalized * 1/2;
T.Ti_Multi = T.Ti_Normalized * 1/2;
T.Al_Multi = T.Al_Normalized * 2/3;
T.Cr_Multi = T.Cr_Normalized * 2/3;
T.Fe_Multi = T.Fe_Normalized * 1/1;
T.Mn_Multi = T.Mn_Normalized * 1/1;
T.Mg_Multi = T.Mg_Normalized * 1/1;
T.Ni_Multi = T.Ni_Normalized * 1/1;
T.Ca_Multi = T.Ca_Normalized * 1/1;
T.Na_Multi = T.Na_Normalized * 2/1;
T.K_Multi = T.K_Normalized * 2/1;
T.S_Multi = T.S_Normalized * 1/3;
T.P_Multi = T.P_Normalized * 2/5;

T.V_Multi = T.V_Normalized * 2/3;
T.Zn_Multi = T.Zn_Normalized * 1;
T.Co_Multi = T.Co_Normalized * 1;
T.Ba_Multi = T.Ba_Normalized * 1;
T.Sr_Multi = T.Sr_Normalized * 1;
T.B_Multi = T.B_Normalized * 2/3;
T.Pb_Multi = T.Pb_Normalized * 1;
T.Cu_Multi = T.Cu_Normalized * 1;
T.Sb_Multi = T.Sb_Normalized * 2/3;
T.As_Multi = T.As_Normalized * 2/5;
T.Th_Multi = T.Th_Normalized * 1/2;
T.Zr_Multi = T.Zr_Normalized * 1/2;
T.Hf_Multi = T.Hf_Normalized * 1/2;
T.Ag_Multi = T.Ag_Normalized * 2/1;
T.Y_Multi = T.Y_Normalized * 2/3;
T.La_Multi = T.La_Normalized * 2/3;
T.Ce_Multi = T.Ce_Normalized * 2/3;
T.Nd_Multi = T.Nd_Normalized * 2/3;

T.Cl_Multi = T.Cl_Normalized * 1 / 1;
T.F_Multi = T.F_Normalized * 1 / 1;

%Calculating A, B, C, D, T, ABCD, and ABCDT
T.A = T.Na_Multi + T.K_Multi;
T.B = T.Ca_Multi + T.Mn_Multi + T.Sr_Multi;
T.C = T.Mg_Multi + T.Fe_Multi +  T.Ni_Multi + T.Zn_Multi + T.Cu_Multi + T.Co_Multi;
T.D = T.Cr_Multi + T.Ti_Multi + T.Al_Multi + T.Zr_Multi + T.V_Multi;
T.T = T.Si_Multi;

T.ABCD = T.A + T.B + T.C + T.D + T.Ce_Multi + T.La_Multi + T.Nd_Multi + T.Y_Multi + T.Hf_Multi + T.Th_Multi;
T.ABCDT = T.A + T.B + T.C + T.D + T.T + T.Ce_Multi + T.La_Multi + T.Nd_Multi + T.Y_Multi + T.Hf_Multi + T.Th_Multi;

T.ABCDS = T.ABCD + T.S_Multi + T.Pb_Multi;
T.ABCDTS = T.ABCDT + T.S_Multi + T.Pb_Multi;

T.ABCDP = T.ABCD + T.P_Multi;
T.ABCDTP = T.ABCDT + T.P_Multi;

end