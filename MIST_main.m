function [] = MIST_main(T_User)

%% "delete" previous export & Read input file

B = table(zeros);

writetable(B, 'table_export.xlsx', 'WriteVariableNames', true, 'Sheet', 'Main Results', 'WriteMode','overwrite');    
writetable(B, 'table_export.xlsx', 'WriteVariableNames', true, 'Sheet','Supplemental Data', 'WriteMode','overwrite');   

T_input = readtable(T_User);

%% Primary calculations

T_input2 = Input_Formatting(T_input);                           % Add missing columns and fill missing data with 0s

T = Main_Calculations1(T_input2); 

%% Set variables

t = 0.1;            %"allowance" value

lowerdelta = 95;    %upper bound for check
upperdelta = 105;   %low

length = size(T,1);

%% 1st Order Flow

for i = 1:length

    % initialize blank columns to avoid errors
    T.group1(i) = {''}; T.group2(i) = {''}; T.group3(i) = {''}; T.group4(i) = {''}; T.species(i) = {''}; T.formula(i) = {''}; T.check(i) = {''};
    T.Fo1(i) = 0; T.Fo(i) = 0; T.Fa(i) = 0; T.An(i) = 0; T.Ab(i) = 0; T.Or(i) = 0; T.Wo(i) = 0; T.En(i) = 0; T.Fs(i) = 0; T.Aeg(i) = 0; T.Jd(i) = 0; T.Di(i) = 0; T.Alm(i) = 0; T.Grs(i) = 0; T.Sps(i) = 0; T.Prp(i) = 0; 
    T.ABCDT_ideal(i) = 0; T.ABCD_ideal(i) = 0; T.T_prime(i) = 0; T.T_ideal(i) = 0;
    T.Delta_ABCDT(i) = 0; T.Delta_ABCD(i) = 0; T.Delta_T(i) = 0;

    if (T.SiO2(i) > 14)
        T.group1(i) = {'Silicate'};

        if(T.SiO2(i) > 85 && T.Na2O(i) < 10 && T.Al2O3(i) < 15)
            T.group1(i) = {'SiO2 Phase'};
        end
        
        if(T.SiO2(i) > 35 && T.SiO2(i) < 38 && T.Al2O3(i) > 60 && T.Al2O3(i) < 65)
            T.group2(i) = {'Al-Silicate'};
        end

    end

    if (T.P2O5(i) > 20)
        T.group1(i) = {'Phosphate'};
        REE_ox = T.La2O3(i) + T.Ce2O3(i) + T.Nd2O3(i);

        if(REE_ox > T.CaO(i))       
            T.group2(i) = {'REE-Phosphate'};
        end

        if(T.Y2O3(i) > T.CaO(i))
            T.group2(i) = {'REE-Phosphate'};
            T.group3(i) = {'Y-Phosphate'};
        end

    end

    if (T.SO3(i) > 12)
        T.group1(i) = {'Sulfate'};
    end

    if(T.Na2O(i) + T.CaO(i) + T.MgO(i) + T.FeO(i) + T.MnO(i) + T.SrO(i) + T.BaO(i) + T.PbO(i) + T.ZnO(i) > 90 && T.SiO2(i) < 1 && T.Al2O3(i) < 0.5 && T.TiO2(i) < 0.5 && T.Cr2O3(i) < 0.5 && T.P2O5(i) < 1 && T.SO3(i) < 4)
        T.group1(i) = {'Carbonate'};

        if(T.FeO(i) > 50 && T.FeO(i) < 99 && T.CaO(i) < 1 && T.CaO(i) > 0.1 && T.Na2O(i) < 1 && T.Al2O3(i) < 0.1 && T.Cr2O3(i) < 0.1)
            T.group2(i) = {'Orthorhombic Carbonate'};
            T.group3(i) = {'Fe-Mg-Sr-Pb-Zn Carbonate'};
            T.group4(i) = {''};
        end

        if(T.CaO(i) + T.MgO(i) + T.MnO(i) + T.SrO(i) + T.BaO(i) > 90 && T.CaO(i) > (T.MgO(i) + T.FeO(i) + T.MnO(i) + T.SrO(i) + T.BaO(i)))
            T.group2(i) = {'Trigonal Carbonate'};
            T.group3(i) = {'Ca-Mg-Fe-Mn Carbonate'};
        end

        if(T.Na2O(i) > 95)
            T.group2(i) = {'Monoclinic Carbonate'};
            T.group3(i) = {'Na-Carbonate'};
        end

    end

    if(T.P2O5(i) < 1 && T.SO3(i) < 1 && T.SiO2(i) < 2 && T.CaO(i) < 0.4 && T.Na2O(i) < 0.1)

        if (strcmp(T.group1(i), 'Carbonate'))
            T.group1(i) = {'Carbonate or Oxide-Hydroxide'};

            if (T.FeO(i) > 95 && T.CaO(i) > 0 && T.CaO(i) < 0.2)
                T.group1(i) = {'Carbonate'};
                T.group2(i) = {'Orthorhombic Carbonate'};
                T.group3(i) = {'Fe-Mg-Sr-Pb-Zn Carbonate'};
                T.group4(i) = {''};
            end
        else
            T.group1(i) = {'Oxide-Hydroxide'};
        end

    end

    if(T.TiO2(i) > 95)
        T.group1(i) = {'Oxide-Hydroxide'};
        T.group2(i) = {'Oxide'};
        T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide'};
        T.group4(i) = {'Ti-Oxide'};
    end

    if(T.Cr2O3(i) > 90)
        T.group1(i) = {'Oxide-Hydroxide'};
        T.group2(i) = {'Oxide'};
        T.group3(i) = {'Other Oxide'};
        T.group4(i) = {'Cr-Oxide'};
    end

    if(T.Al2O3(i) > 96)
        T.group1(i) = {'Oxide-Hydroxide'};
        T.group2(i) = {'Oxide or Hydroxide'};
        T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide'};
        T.group4(i) = {'Al-Oxide or Hydroxide'};
    end

    if(T.MnO(i) > 40)

        if (T.CaO(i) >= 0.4 && T.SiO2(i) < 5)
            T.group1(i) = {'Carbonate or Oxide-Hydroxide'};
            T.group2(i) = {'Orthorhombic Carbonate'};
            T.group3(i) = {'Fe-Mg-Sr-Pb-Zn Carbonate'};
            T.group4(i) = {''};
        end

        if (T.CaO(i) > 0 && T.CaO(i) < 0.4 && T.SiO2(i) < 5)
            T.group1(i) = {'Oxide-Hydroxide'};
            T.group2(i) = {'Oxide or Hydroxide'};
            T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
            T.group4(i) = {'Mn-Oxide or Mn-Hydroxide'};
        end

        if (T.CaO(i) == 0 && T.SiO2(i) < 5)
            T.group1(i) = {'Oxide-Hydroxide'};
            T.group2(i) = {'Oxide or Hydroxide'};
            T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
            T.group4(i) = {'Mn-Oxide or Mn-Hydroxide'};
        end

        if(T.SiO2(i) > 30)
            T.group1(i) = {'Silicate'};
        end

    end

    if (T.Cl(i) > 28)
        T.group1(i) = {'Halide'};
    end

    if (T.F(i) > 18)
        T.group1(i) = {'Halide'};
    end

    if(T.SiO2(i) > 29 && T.SiO2(i) < 33 && T.ZrO2(i) > 40 && T.ZrO2(i) < 67)
        T.group1(i) = {'Silicate'};
        T.group2(i) = {'Nesosilicate'};
        T.group3(i) = {'Zircon'};
    end

    if(T.ZrO2(i) > 90)
        T.group1(i) = {'Oxide-Hydroxide'};
        T.group2(i) = {'Oxide'};
        T.group3(i) = {'Other Oxide'};
        T.group4(i) = {''};
    end

end

%% Start of Suborders

for i = 1:length

    % Set variable names for cations normalized to 24 oxygen formula units (n)

    Ca_n = T.Ca_Multi(i);
    Mg_n = T.Mg_Multi(i);
    Mn_n = T.Mn_Multi(i);
    Fe_n = T.Fe_Multi(i);
    Al_n = T.Al_Multi(i);
    Si_n = T.Si_Multi(i);
    Cr_n = T.Cr_Multi(i);
    Sr_n = T.Sr_Multi(i);
    Ba_n = T.Ba_Multi(i);
    P_n = T.P_Multi(i);
    Na_n = T.Na_Multi(i);
    Ti_n = T.Ti_Multi(i);
    Ag_n = T.Ag_Multi(i);
    S_n =  T.S_Multi(i);
    V_n =  T.V_Multi(i);
    Zn_n = T.Zn_Multi(i);
    Ni_n = T.Ni_Multi(i);
    Cl_n = T.Cl_Multi(i);
    K_n = T.K_Multi(i);
    F_n = T.F_Multi(i);
    Co_n = T.Co_Multi(i);

    Y_n = T.Y_Multi(i);
    La_n = T.La_Multi(i);
    Ce_n = T.Ce_Multi(i);
    Nd_n = T.Nd_Multi(i);
    Pb_n = T.Pb_Multi(i);
    Cu_n = T.Cu_Multi(i);
    Th_n = T.Th_Multi(i);
    Zr_n = T.Zr_Multi(i);
    Hf_n = T.Hf_Multi(i);

    A_n = T.A(i);
    B_n = T.B(i);
    C_n = T.C(i);
    D_n = T.D(i);
    T_n = T.T(i);

    ABCD_n = A_n + B_n + C_n + D_n;
    ABCDT_n = A_n + B_n + C_n + D_n + T_n;

    % Set variable names for cations normalized to 2 oxygen formula units (n2)

    Renorm_factor = 12;

    Si_n2 = T.Si_Multi(i) / Renorm_factor;  
    P_n2 = T.P_Multi(i) / Renorm_factor;
    S_n2 = T.S_Multi(i) / Renorm_factor;
    Ca_n2 = T.Ca_Multi(i) / Renorm_factor;
    Fe_n2 = T.Fe_Multi(i) / Renorm_factor;
    Mg_n2 = T.Mg_Multi(i) / Renorm_factor;
    Mn_n2 = T.Mn_Multi(i) / Renorm_factor;
    Sr_n2 = T.Sr_Multi(i) / Renorm_factor;
    Ba_n2 = T.Ba_Multi(i) / Renorm_factor;
    Pb_n2 = T.Pb_Multi(i) / Renorm_factor;
    Zn_n2 = T.Zn_Multi(i) / Renorm_factor;
    Na_n2 = T.Na_Multi(i) / Renorm_factor;
    Ti_n2 = T.Ti_Multi(i) / Renorm_factor;
    Al_n2 = T.Al_Multi(i) / Renorm_factor;
    K_n2 = T.K_Multi(i) / Renorm_factor;
    Ni_n2 = T.Ni_Multi(i) / Renorm_factor;
    Cr_n2 = T.Cr_Multi(i) / Renorm_factor;
    Cl_n2 = T.Cr_Multi(i) / Renorm_factor;
    F_n2 = T.F_Multi(i) / Renorm_factor;
    Cu_n2 = T.Cu_Multi(i) / Renorm_factor;
    V_n2 = T.V_Multi(i) / Renorm_factor;
    Co_n2 = T.Co_Multi(i) / Renorm_factor;

    T_n2 = Si_n2;

    if(T.MgO(i) >= 95)
        T.group2(i) = {'Carbonate or Oxide-Hydroxide'};
        T.group3(i) = {''};
        T.group4(i) = {''};
        T.species(i) = {'Mg Carbonate or Oxide-Hydroxide'};
        T.formula(i) = {''};
    elseif(T.FeO(i) >= 95)
        T.group2(i) = {'Carbonate or Oxide-Hydroxide'};
        T.group3(i) = {''};
        T.group4(i) = {''};
        T.species(i) = {'Fe Oxide-Hydroxide or Carbonate'};
        T.formula(i) = {''};
    elseif(T.MnO(i) >= 95)
        T.group2(i) = {'Carbonate or Oxide-Hydroxide'};
        T.group3(i) = {''};
        T.group4(i) = {''};
        T.species(i) = {'Mn Oxide-Hydroxide or Carbonate'};
        T.formula(i) = {''};
    elseif(T.Cr2O3(i) >= 94)
        T.group2(i) = {'Carbonate or Oxide-Hydroxide'};
        T.group3(i) = {''};
        T.group4(i) = {'Cr-Oxide'};
        T.species(i) = {'Cr Oxide-Hydroxide'};
        T.formula(i) = {''};
    end

    %% 1. SiO2 Phase
    if(strcmp(T.group1(i), 'SiO2 Phase'))

        if(T.SiO2(i) > 95 && T.Al2O3(i) < 1)
            T.species(i) = {'Quartz (or other SiO2 polymorph)'};
            [sorted_formula, total] = Formula_Output({Si_n2}, {'Si'});
            T.formula(i) = {sorted_formula + 'O2'};
        end

        if(T.SiO2(i) > 95 && T.Al2O3(i) > 1 && T.Al2O3(i) < 5)
            T.species(i) = {'SiO2 glass (or amorphous silica)'};
            [sorted_formula, total] = Formula_Output({Si_n2, Al_n2, Cr_n2, Fe_n2, Mg_n2, Ca_n2, Na_n2, K_n2}, {'Si', 'Al','Cr','Fe','Mg','Ca','Na','K'});
            T.formula(i) = {strcat(sorted_formula,total, 'O2·nH2O')};
        end

        if(T.SiO2(i) > 75 && T.SiO2(i) < 90 && T.Al2O3(i) > 10 && T.Al2O3(i) < 15 && T.Na2O(i) >= 0 && T.Na2O(i) < 10)
            T.species(i) = {'Obsidian glass'};
            [sorted_formula, total] = Formula_Output({Si_n2, Al_n2, Na_n2, K_n2, Ca_n2, Fe_n2, Mg_n2}, {'Si', 'Al','Na','K','Ca','Fe','Mg'});
            T.formula(i) = {strcat(sorted_formula,total, 'O2·nH2O')};
        end

    end %end of SiO2 phase

    %% 2. Carbonates
    if(strcmp(T.group1(i), 'Carbonate') || strcmp(T.group1(i), 'Carbonate or Oxide-Hydroxide'))

        if(Si_n2 > 0.1)
            T.group1(i) = {'Silicate'};
            T.group2(i) = {''};
            T.group3(i) = {''};
            T.group4(i) = {''};
            T.species(i) = {''};
            T.formula(i) = {''};
        end

        if(strcmp(T.group3(i), 'Fe-Mg-Sr-Pb-Zn Carbonate'))
            if(Mn_n2 / (Ca_n2 + Mg_n2 + Fe_n2 + Mn_n2 + Sr_n2 + Ba_n2) > 1.75)
                if(Ca_n2 > 0)
                    T.group2(i) = {'Orthorhombic Carbonate'};
                    T.group3(i) = {'Fe-Mg-Sr-Pb-Zn Carbonate'};
                    T.group4(i) = {''};
                    T.species(i) = {'rhodochrosite'};
                    [sorted_formula, total] = Formula_Output({Mn_n/24, Ca_n/24, Fe_n/24, Sr_n/24, Ba_n/24}, {'Mn','Ca','Fe','Sr','Ba'});
                    T.formula(i) = {strcat(sorted_formula, total, '(CO3)')};
                end
            end

            if(1-(Ca_n2+Mg_n2+K_n2) > 1 - (1)*t && 1-(Ca_n2+Mg_n2+K_n2) < 1 + (1)*t && T.CaO(i) > 1)
                if(Ca_n2/Mn_n2 < 0.5 && T.CaO(i) > 1 && Ca_n/(Na_n+K_n) > 1)
                    T.group1(i) = {'Oxide-Hydroxide'};
                    T.group2(i) = {'Hydroxide'};
                    T.group3(i) = {'Fe-Al-Ti-Mn-Mg Hydroxide'};
                    T.group4(i) = {'Mn-Hydroxide'};
                    T.species(i) = {'takanelite'};
                    [sorted_formula, total] = Formula_Output({Mn_n2*3,Ca_n2*3}, {'Mn2+','Ca'});
                    T.formula(i) = {'(' + sorted_formula + ')_{Σ=' + total + '}(Mn4+)1-xO2·0.7H2O'};
                end

            end

            if(Na_n2 + K_n2 > 0 && Na_n > Ca_n)
                T.group1(i) = {'Oxide-Hydroxide'};
                T.group2(i) = {'Oxide or Hydroxide'};
                T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
                T.group4(i) = {'Mn-Oxide or Mn-Hydroxide'};
                T.species(i) = {'birnessite'};
                Fe3 = Fe_n2*2;
                Mn3 = 2 - (Fe3 + Cr_n2*2 + Ti_n2*2 + Mg_n2*2);
                if(Mn3<0), Mn3=0; end
                Mn4 = Mn_n2*2 - Mn3;
                if(Mn4<0), Mn4=0; end
                [sorted_formula1, total1] = Formula_Output({Na_n2*2, Ca_n2*2, K_n2*2}, {'Na','Ca','K'});
                [sorted_formula2, total2] = Formula_Output({Mn4,Mn3,Cr_n2*2,Ti_n2*2,Mg_n2*2,Fe3}, {'Mn4+','Mn3+','Cr','Ti','Mg','Fe3+'});
                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O4·1.5H2O')};
            end

        end

        %% 2.1 Trigonal Carbonates
        if(Ca_n2 > Fe_n2 && Ca_n2 / (Mg_n2 + Fe_n2 + Mn_n2 + Sr_n2 + Ba_n2) >= 0.95 )   
            T.group2(i) = {'Trigonal Carbonate'};
            T.group3(i) = {'Ca-Mg-Fe-Mn Carbonate'};
            T.group4(i) = {''};
            T.species(i) = {''};
            T.formula(i) = {''};

            if(Mn_n2/(Ca_n2 + Mg_n2 + Fe_n2 + Sr_n2 + Ba_n2) > 0.7 - (0.7)*t && Mn_n2/(Ca_n2 + Mg_n2 + Fe_n2 + Sr_n2 + Ba_n2) < 1 + (1)*t && T.CaO(i) > 1)
                T.species(i) = {'kutnohorite'};
                [sorted_formula1, total1] = Formula_Output({Ca_n2,Na_n2}, {'Ca','Na'});
                [sorted_formula2, total2] = Formula_Output({Mn_n2, Mg_n2, Fe_n2},{'Mn2+','Mg','Fe'});
                T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, '(CO3)2')};
            end

            if(Ca_n2 / (Mg_n2 + Fe_n2 + Mn_n2 + Sr_n2 + Ba_n2) > 1.75)
                T.group4(i) = {'Ca-Carbonate'};
                T.species(i) = {'calcite'};
                [sorted_formula, total] = Formula_Output({Ca_n2/2,Mg_n2/2,Fe_n2/2,Mn_n2/2,Sr_n2/2,Ba_n2/2}, {'Ca', 'Mg', 'Fe', 'Mn','Sr','Ba'});
                T.formula(i) = {strcat(sorted_formula,total,'(CO3)')};

                if(Mg_n2 > 0.1 && Mg_n2 < 0.3)
                    T.species(i) = {'magnesian calcite'};
                    [sorted_formula, total] = Formula_Output({Ca_n2/2,Mg_n2/2,Fe_n2/2,Mn_n2/2}, {'Ca', 'Mg', 'Fe', 'Mn'});
                    T.formula(i) = {strcat(sorted_formula,total,'(CO3)')};
                end

            end

            if(Ca_n2 / (Mg_n2 + Fe_n2 + Mn_n2) > 1 - 1*t && Ca_n2 / (Mg_n2 + Fe_n2 + Mn_n2) < 1 + 1*t && Ca_n2 / Mg_n2 > 1 - 1*t && Ca_n2 / Mg_n2 < 1 + 1*t)
                T.species(i) = {'dolomite'};
                Mgx = 1 - Ca_n2;
                if(Mgx < 0), Mgx = 0; end
                Mgy = Mg_n2 - Mgx;
                [sorted_formula1, total1] = Formula_Output({Ca_n2,Mgx}, {'Ca', 'Mg'});
                [sorted_formula2, total2] = Formula_Output({Mgy, Fe_n2, Mn_n2, Sr_n2, Ba_n2}, {'Mg','Fe','Mn','Sr','Ba'});
                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2, total2,'(CO3)2')};

                if(Mg_n2 > Fe_n2 && Fe_n2 > 0.1)
                    T.species(i) = {'ferroan dolomite'};
                    Cay = Ca_n2 - (Mg_n2 + Fe_n2 + Mn_n2 + Ba_n2 + Sr_n2);
                    Cax = Ca_n2 - Cay;
                    [sorted_formula1, total1] = Formula_Output({Mg_n2 + Fe_n2 + Mn_n2 + Ba_n2 + Sr_n2, Cay},{'Mg','Fe','Mn','Ba','Sr','Ca'});
                    T.formula(i) = {strcat('Ca_', num2str(round(Cax,2)), sorted_formula1, total1, '(CO3)2')};
                end

            end

            if(Fe_n2 > Mg_n2)
                T.species(i) = {'ankerite'};
                Cax = 1;
                Cay = Ca_n2-Cax;
                if(Cay<0), Cay = 0; end
                [sorted_formula, total] = Formula_Output({Cay,Mg_n2,Fe_n2,Mn_n2, Zn_n2}, {'Ca', 'Mg', 'Fe', 'Mn','Zn'});
                T.formula(i) = {strcat('Ca_1',sorted_formula,total,'(CO3)2')};
            end

        end %end of Trigonal Carbonates

        %% 2.2 Orthorhombic Carbonates
      
        if(Mg_n2 > 1.75 && T.MgO(i) < 100)
            T.group2(i) = {'Orthorhombic Carbonate'};
            T.group3(i) = {'Fe-Mg-Sr-Pb-Zn Carbonate'};
            T.group4(i) = {''};
            T.species(i) = {'magnesite'};
            [sorted_formula, total] = Formula_Output({Fe_n2/2,Mn_n2/2,Mg_n2/2,Ca_n2/2,Sr_n2/2,Ba_n2/2}, {'Fe','Mn','Mg','Ca','Sr','Ba'});
            T.formula(i) = {strcat(sorted_formula,total, '(CO3)')};
        end

        if(Sr_n2 > 1.75)
            T.group2(i) = {'Orthorhombic Carbonate'};
            T.group3(i) = {'Fe-Mg-Sr-Pb-Zn Carbonate'};
            T.group4(i) = {''};
            T.species(i) = {'strontianite'};
            [sorted_formula, total] = Formula_Output({Fe_n2/2,Mg_n2/2,Ca_n2/2,Zn_n2/2, Sr_n2/2,Ba_n2/2}, {'Fe','Mg','Ca','Zn','Sr','Ba'});
            T.formula(i) = {strcat(sorted_formula,total, '(CO3)')};
        end

        if(Ba_n2 > 1.75)
            T.group2(i) = {'Orthorhombic Carbonate'};
            T.group3(i) = {'Fe-Mg-Sr-Pb-Zn Carbonate'};
            T.group4(i) = {''};
            T.species(i) = {'witherite'};
            [sorted_formula, total] = Formula_Output({Fe_n2/2,Mg_n2/2,Ca_n2/2,Zn_n2/2, Sr_n2/2,Ba_n2/2}, {'Fe','Mg','Ca','Zn','Sr','Ba'});
            T.formula(i) = {strcat(sorted_formula,total, '(CO3)')};
        end

        if(Pb_n2 > 1.75)
            T.group2(i) = {'Orthorhombic Carbonate'};
            T.group3(i) = {'Fe-Mg-Sr-Pb-Zn Carbonate'};
            T.group4(i) = {''};
            T.species(i) = {'cerussite'};
            [sorted_formula, total] = Formula_Output({Pb_n2*(3/2),Zn_n2*(3/2),Fe_n2*(3/2),Ca_n2*(3/2),Mn_n2*(3/2),Sr_n2*(3/2),Ba_n2*(3/2)}, {'Pb','Zn','Fe','Ca','Mn','Sr','Ba'});
            T.formula(i) = {strcat(sorted_formula, total, '(CO3)')};
        end

        if(Zn_n2 > 1.75)
            T.group2(i) = {'Orthorhombic Carbonate'};
            T.group3(i) = {'Fe-Mg-Sr-Pb-Zn Carbonate'};
            T.group4(i) = {''};
            T.species(i) = {'smithsonite'};
            [sorted_formula, total] = Formula_Output({Pb_n/24,Zn_n/24, Sr_n/24, Ba_n/24, Ca_n/24}, {'Pb','Zn','Sr','Ba','Ca'});
            T.formula(i) = {strcat(sorted_formula, total, '(CO3)')};
        end

        %% 2.3 Monoclinic Carbonates
        if(Na_n/T.ABCDT(i) > 0.95)
            T.group2(i) = {'Monoclinic Carbonate'};
            T.group3(i) = {'Na-Carbonate'};
            T.group4(i) = {''};
            T.species(i) = {'natrite or natron'};
            [sorted_formula1, total1] = Formula_Output({Na_n/24,K_n/24, Ca_n/24, Mg_n/24},{'Na','K','Ca','Mg'});
            T.formula(i) = {strcat(sorted_formula1, total1, '(CO3) or', sorted_formula1, total1,'(CO3)•10H2O')};

            if (T.SO3(i) > 1 && T.SO3(i) < 4)
                T.species(i) = {'natron'};
                T.formula(i) = {strcat(sorted_formula1, total1,'(CO3)•10H2O')};
            end

        end %end of Monoclinic Carbonates

    end

    %% 3. Sulfates
    if(strcmp(T.group1(i), 'Sulfate'))

        if(Fe_n > Mg_n && Fe_n > Al_n && Fe_n > Mn_n && Fe_n >= Cu_n) 
            if((Fe_n + Cu_n) / S_n > 1-1*t && (Fe_n + Cu_n) / S_n < 1.75+1.75*t)  
                if(Cu_n/Fe_n > 0 && Cu_n/Fe_n < 1+1*t)    
                    if((Mg_n + Al_n + Mn_n)/24 > 0 && (Mg_n + Al_n + Mn_n)/24 < 1.53+1.53*t)
                        T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
                        T.group3(i) = {'Fe-Sulfate'};
                        T.group4(i) = {''};
                        T.species(i) = {'valleriite'};
                        [sorted_formula, total] = Formula_Output({Fe_n2/2,Cu_n2/2}, {'Fe','Cu'});
                        T.formula(i) = {'2[(' + sorted_formula + ')S]·1.53[(Mg,Al)(OH)2]'};  
                    end
                end
            end
        end %end of valleriite

        %% 3.1 Ca-Sulfate
        if(Ca_n / (Na_n + Ca_n + Mg_n + Fe_n + Mn_n) > 0.5 && Ca_n > Mn_n && (Fe_n + Al_n)/Ca_n < 0.285 && Na_n/Ca_n < 0.1)
            T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
            T.group3(i) = {'Na-Ca-Sulfate'};
            T.group4(i) = {'Ca-Sulfate'};
            T.species(i) = {'anhydrite/gypsum/bassanite'};
            [sorted_formula1, total1] = Formula_Output({Ca_n2*2,Na_n2*2,K_n2*2,Mg_n2*2,Fe_n2*2,Mn_n2*2},{'Ca','Na','K','Mg','Fe','Mn'});
            [sorted_formula2, total2] = Formula_Output({S_n2*2,P_n2*2},{'S','P'});
            T.formula(i) = {strcat(sorted_formula1,total1,'[',sorted_formula2,total2,'O4] / ·2H2O / ·0.5H2O ')};
        end %end of Ca-Sulfate

        %% 3.2 Na-Sulfate
        if(Na_n / (Na_n + Ca_n + Mg_n + Fe_n + Mn_n) > 0.5 && Na_n > K_n)
            T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
            T.group3(i) = {'Na-Ca-Sulfate'};
            T.group4(i) = {'Na-Sulfate'};

            if(Na_n / S_n > 2-2*t && Na_n / S_n < 2+2*t && Fe_n < 0.5)
                T.species(i) = {'thénardite or mirabilite'};
                [sorted_formula1, total1] = Formula_Output({Na_n/6, Ca_n/6, Mg_n/6, Fe_n/6},{'Na','Ca','Mg', 'Fe'});
                [sorted_formula2, total2] = Formula_Output({S_n/6, P_n/6},{'S','P'});
                T.formula(i) = {strcat(sorted_formula1, total1, '(', sorted_formula2, total2, 'O4)·10H2O')};
            end

            if(Na_n / S_n > (3/2)-(3/2)*t && Na_n / S_n < (3/2)+(3/2)*t && Fe_n < 0.005) 
                T.species(i) = {'ivsite'};
                renorm = 8/24;
                [sorted_formula1, total1] = Formula_Output({Na_n*renorm, K_n*renorm, Ca_n*renorm},{'Na','K','Ca'});
                [sorted_formula2, total2] = Formula_Output({S_n*renorm/2, P_n*renorm/2},{'S','P'});
                T.formula(i) = {strcat(sorted_formula1,total1, 'H(',sorted_formula2,total2,'O4)2')};
            end

        end %end of Na-Sulfate

        %% 3.3 Mg-Sulfate-H2O
        if(Mg_n / (Na_n + Ca_n + Mg_n + Fe_n + Mn_n) > 0.5)
            T.group2(i) = {'Mg-Sulfate-H2O'};

            if((Mg_n+Cu_n+Zn_n) / S_n > 1-1*t && (Mg_n+Cu_n+Zn_n) / S_n < 1+1*t && Cu_n > 0)
                T.species(i) = {'alpersite'};
                [sorted_formula1, total1] = Formula_Output({Mg_n2*2,Fe_n2*2,Mn_n2*2,Ca_n2*2,Na_n2*2},{'Mg','Fe','Mn','Ca','Na'});
                [sorted_formula2, total2] = Formula_Output({S_n2*2,P_n2*2},{'S','P'});
                T.formula(i) = {strcat(sorted_formula1,total1,'[',sorted_formula2,total2,'O4]·7H2O')};
            end

            if(Mg_n / S_n > 1-1*t && Mg_n / S_n < 1+1*t && Cu_n == 0)
                T.species(i) = {'kieserite/epsomite/meridianiite/sanderite'};
                [sorted_formula1, total1] = Formula_Output({Mg_n2*2, Fe_n2*2, Mn_n2*2, Ca_n2*2, Na_n2*2, K_n2*2},{'Mg','Fe','Mn', 'Ca','Na','K'});
                [sorted_formula2, total2] = Formula_Output({S_n2*2, P_n2*2},{'S','P'});
                T.formula(i) = {strcat(sorted_formula1,total1, '(',sorted_formula2,total2,'O4)·nH2O')};
            end

            if(Mg_n / S_n > (7/5)-(7/5)*t && Mg_n / S_n < (7/5)+(7/5)*t)
                T.species(i) = {'caminite'};
                OH = 4 - (Cl_n2*11 + F_n2*11);
                if(OH < 0), OH = 0; end
                [sorted_formula1,total1] = Formula_Output({Mg_n2*11,Ca_n2*11,Mn_n2*11,Fe_n2*11,Na_n2*11},{'Mg','Ca','Mn','Fe','Na'});
                [sorted_formula2,total2] = Formula_Output({(S_n2*11)/5,(P_n2*11)/5},{'S','P'});
                [sorted_formula3,total3] = Formula_Output({OH,Cl_n2*11,F_n2*11},{'OH','Cl','F'});
                T.formula(i) = {strcat(sorted_formula1,total1,'[',sorted_formula2,total2,'O4]5',sorted_formula3,total3,'·H2O')};
            end

            if((Mg_n+Mn_n+Zn_n) / S_n > 7.5-7.5*t && (Mg_n+Mn_n+Zn_n) / S_n < 7.5+7.5*t)
                if(Mg_n/(Mg_n+Mn_n+Zn_n) > 0.5)
                    T.species(i) = {'mooreite'};
                    renorm = 21/24;
                    OH = 26 - (F_n*renorm + Cl_n*renorm);
                    if(OH<=0), OH = 0; end
                    [sorted_formula1, total1] = Formula_Output({Mg_n*renorm, Mn_n*renorm, Zn_n*renorm, Fe_n*renorm},{'Mg','Mn','Zn','Fe'});
                    [sorted_formula2, total2] = Formula_Output({S_n*renorm, P_n*renorm}, {'S','P'});
                    [sorted_formula3, total3] = Formula_Output({OH/2, F_n*renorm/2, Cl_n*renorm/2}, {'OH','F','Cl'}); 
                    T.formula(i) = {strcat(sorted_formula1,total1, '(',sorted_formula2,total2,'O4)',sorted_formula3, total3, '·8H2O')};
                end
            end

        end %end of Mg-Sulfate-H2O

        %% 3.4 Na-Mg-Sulfate
        if(Mg_n / (Ca_n + Mg_n + Fe_n + Mn_n) > 0.5 && Na_n / (Ca_n + Fe_n + Mn_n) > 0.5)
            T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
            T.group3(i) = {'Na-Mg-Sulfate'};
            T.group4(i) = {''};

            if(Na_n / S_n > 2.1-2.1*t && Na_n / S_n < 2.1+2.1*t)
                if((Mg_n + Fe_n)/S_n > 0.1-0.1*t && (Mg_n + Fe_n)/S_n < 0.1+0.1*t)
                    if(Mg_n > Fe_n)
                        if(Cl_n / S_n > (3/10)-(3/10)*0.1 && Cl_n / S_n < (3/10)+(3/10)*0.1)
                            T.species(i) = {'D`ansite'};
                            Renorm_factor = 24/41.5;
                            Na_n41 = T.Na_Multi(i) / Renorm_factor;
                            Ca_n41 = T.Ca_Multi(i) / Renorm_factor;
                            Mg_n41 = T.Mg_Multi(i) / Renorm_factor;
                            Fe_n41 = T.Fe_Multi(i) / Renorm_factor;
                            S_n41 = T.S_Multi(i) / Renorm_factor;
                            P_n41 = T.P_Multi(i) / Renorm_factor;
                            Cl_n41 = T.Cl_Multi(i) / Renorm_factor;
                            F_n41 = T.F_Multi(i) / Renorm_factor;
                            Mn_n41 = T.F_Multi(i) / Renorm_factor;
                            OH = 3 - (Cl_n41+F_n41);
                            [sorted_formula1, total1] = Formula_Output({Na_n41, Ca_n41}, {'Na','Ca'});
                            [sorted_formula2, total2] = Formula_Output({Mg_n41, Fe_n41, Mn_n41}, {'Mg', 'Fe', 'Mn'});
                            [sorted_formula3, total3] = Formula_Output({S_n41/10,P_n41/10}, {'S','P'});
                            [sorted_formula4, total4] = Formula_Output({OH, Cl_n41, F_n41}, {'OH','Cl','F'});
                            T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2,'[',sorted_formula3,total3, 'O4]10',sorted_formula4,total4)};
                        end
                    end
                end
            end

            if((Na_n+Ca_n) / (Mg_n+Fe_n+Mn_n) > 6-6*t && (Na_n+Ca_n) / (Mg_n+Fe_n+Mn_n) < 6+6*t)
                if((Mg_n+Fe_n+Mn_n) / S_n > 0.25-0.25*t && (Mg_n+Fe_n+Mn_n) / S_n < 0.25+0.25*t)
                    T.species(i) = {'vanthoffite'};
                    renorm = 16/24;
                    [sorted_formula1, total1] = Formula_Output({Na_n*renorm, K_n*renorm},{'Na','K'});
                    [sorted_formula2, total2] = Formula_Output({Mg_n*renorm, Ca_n*renorm, Fe_n*renorm, Mn_n*renorm}, {'Mg','Ca','Fe','Mn'});
                    [sorted_formula3, total3] = Formula_Output({S_n*renorm/4, P_n*renorm/4}, {'S','P'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, '(',sorted_formula3, total3, 'O4)4')};
                end
            end

            if((Na_n+Ca_n)/(Mg_n+Fe_n+Mn_n) > 2-2*2*t && (Na_n+Ca_n)/(Mg_n+Fe_n+Mn_n) < 2+2*2*t && Na_n > Ca_n)
                if((Mg_n + Fe_n + Mn_n)/S_n > 1/2-2*(1/2)*t && (Mg_n + Fe_n + Mn_n)/S_n < 1/2 + 2*(1/2)*t && Mg_n/(Mg_n+Fe_n+Mn_n) > 0.5)
                    T.species(i) = {'blödite'};
                    [sorted_formula1,total1] = Formula_Output({Na_n2*4,Ca_n2*4,K_n2*4},{'Na','Ca','K'});
                    [sorted_formula2,total2] = Formula_Output({Mg_n2*4,Fe_n2*4,Mn_n2*4},{'Mg','Fe','Mn'});
                    [sorted_formula3,total3] = Formula_Output({(S_n2*4)/2,(P_n2*4)/2},{'S','P'});
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'[',sorted_formula3,total3,'O4]2·4H2O')};
                end
            end

        end %end of Na-Mg Sulfates

        %% 3.5 Ca-Fe-Sulfate
        if(Fe_n / (Na_n + Mg_n + Fe_n + Mn_n) > 0.5 && Ca_n / (Na_n + Mg_n + Mn_n) > 0.5 && (Fe_n + Al_n)/Ca_n >= 0.3-0.3*t && (Fe_n + Al_n)/Ca_n < 4+4*t)
            T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
            T.group3(i) = {'Ca-Fe-Sulfate'};

            if(Ca_n / (Fe_n + Al_n) > 0.5-0.5*t && Ca_n / (Fe_n + Al_n) < 0.5+0.5*t)
                if(Fe_n / (Ca_n + Al_n) > 0.5-0.5*t && Fe_n / (Ca_n + Al_n) < 0.5+0.5*t)
                    if(Al_n / (Ca_n + Fe_n) > 0.5-0.5*t && Al_n / (Ca_n + Fe_n) < 0.5+0.5*t)
                        if((Ca_n + Fe_n + Al_n) / S_n > 3-3*t && (Ca_n + Fe_n + Al_n) / S_n < 3+3*t)
                            T.species(i) = {'vyalsovite'};
                            renorm = 2.5/24;
                            Fe3 = 1 - (Al_n*renorm + Cr_n*renorm);
                            if(Fe3<=0), Fe3 = 0; end
                            Fe2 = Fe_n*renorm - Fe3;
                            if(Fe2<=0), Fe2=0; end
                            OH = 5 - (F_n*renorm + Cl_n*renorm);
                            if(OH<=0), OH=0; end 
                            [sorted_formula1, total1] = Formula_Output({Ca_n*renorm, Na_n*renorm, K_n*renorm}, {'Ca','Na','K'});
                            [sorted_formula2, total2] = Formula_Output({Fe2, Mg_n*renorm, Mn_n*renorm}, {'Fe2+', 'Mg','Mn'});
                            [sorted_formula3, total3] = Formula_Output({Al_n*renorm, Fe3, Cr_n*renorm}, {'Al','Fe3+','Cr'});
                            [sorted_formula4, total4] = Formula_Output({OH, F_n*renorm, Cl_n*renorm}, {'OH','F','Cl'});
                            T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, 'S', sorted_formula4, total4)};
                        end
                    end
                end
            end

            if(Fe_n / S_n > (2/3)-(2/3)*t && Fe_n / S_n < (2/3)+(2/3)*t)
                if(Ca_n / Fe_n > 0.25-0.25*t && Ca_n / Fe_n < 0.25+0.25*t)
                    T.species(i) = {'calciocopiapite'};
                    Fe3 = 4 - (Al_n*12.5+Cr_n*12.5+V_n*12.5);
                    if(Fe3<0), Fe3 = 0; end
                    OH = 2 - (Cl_n*12.5+F_n*12.5);
                    if(OH<0), OH = 0; end
                    [sorted_formula1,total1] = Formula_Output({Ca_n*12.5,Na_n*12.5,K_n*12.5,Mg_n*12.5},{'Ca','Na','K','Mg'});
                    [sorted_formula2,total2] = Formula_Output({Fe3,Al_n*12.5,Cr_n*12.5,V_n*12.5},{'Fe3+','Al','Cr','V'});
                    [sorted_formula3,total3] = Formula_Output({(S_n*12.5)/6,(P_n*12.5)/6},{'S','P'});
                    [sorted_formula4,total4] = Formula_Output({OH, Cl_n*12.5, F_n*12.5},{'OH','Cl','F'}); 
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'[',sorted_formula3,total3,'O4]6',sorted_formula4,total4,'·20H2O')};
                end
            end

            if(Ca_n / S_n > (6/2.5)-(6/2.5)*t && Ca_n / S_n < (6/2)+(6/2)*t)
                if((Fe_n+Al_n) / Ca_n > (1/3)-2*(1/3)*t && (Fe_n+Al_n) / Ca_n < (1/3)+(1/3)*t)
                    if((Fe_n+Al_n) / S_n > (2/2.5)-(2/2.5)*t && (Fe_n+Al_n) / S_n < 1+1*t)
                        T.species(i) = {'sturmanite'};
                        renorm = 17/24;
                        [sorted_formula1, total1] = Formula_Output({Ca_n*renorm, Na_n*renorm, K_n*renorm},{'Ca','Na','K'});
                        [sorted_formula2, total2] = Formula_Output({Fe_n*renorm, Al_n*renorm, Mg_n*renorm, Mn_n*renorm},{'Fe','Al','Mg','Mn'});
                        [sorted_formula3, total3] = Formula_Output({S_n*renorm/2.5, P_n*renorm/2.5},{'S','P'});
                        T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, '(', sorted_formula3, total3, 'O4)2.5[B(OH)4](OH)12·25H2O')};
                    end
                end
            end

        end %end of Ca-Fe-Sulfates

        %% 3.6 Fe-Sulfate
        if(Fe_n / (Na_n + K_n + Ca_n + Mg_n + Fe_n + Mn_n + Pb_n + Cu_n) > 0.5 && Na_n/Fe_n < 0.1 && Ca_n/Fe_n < 0.05)
            T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
            T.group3(i) = {'Fe-Sulfate'};
            T.group4(i) = {''};

            if(Fe_n / S_n > 0.5-0.5*t && Fe_n / S_n < 0.5+0.5*t)
                T.species(i) = {'rhomboclase'};
                renorm = 8/24;
                [sorted_formula1, total1] = Formula_Output({Fe_n*renorm, Al_n*renorm, Cr_n*renorm, Ti_n*renorm, Mg_n*renorm, Ca_n*renorm},{'Fe','Al','Cr','Ti','Mg','Ca'});
                [sorted_formula2, total2] = Formula_Output({S_n*renorm/2, P_n*renorm/2},{'S','P'});
                T.formula(i) = {strcat('(H5O2)', sorted_formula1, total1, '(', sorted_formula2, total2, 'O4)2·2H2O')};
            end

            if(Al_n/Fe_n > (1/3) - (1/3)*t && Al_n/Fe_n < (1/3) + (1/3)*t)
                Fe3 = 1 - (Al_n+Cr_n+Ti_n);
                if(Fe3<0), Fe3 = 0; end
                Fe2 = Fe_n - Fe3;
                T.species(i) = {'coquimbite'};
                [sorted_formula1,total1] = Formula_Output({Al_n, Fe3, Cr_n, Ti_n},{'Al','Fe3+','Cr','Ti'});
                [sorted_formula2,total2] = Formula_Output({Fe2, Mg_n,Mn_n,Ca_n},{'Fe2+', 'Mg','Mn','Ca'});
                [sorted_formula3,total3] = Formula_Output({S_n/6, P_n/6},{'S','P'});
                [sorted_formula4,total4] = Formula_Output({Fe_n/2, Mn_n/2, Al_n/2}, {'Fe','Mn','Al'});
                [sorted_formula5,total5] = Formula_Output({S_n/2, P_n/2},{'S','P'});   
                T.formula(i) = {strcat(sorted_formula4, total4, '(', sorted_formula5, total5, 'O4)3) or', sorted_formula1,total1,sorted_formula2,total2,'(',sorted_formula3,total3,'O4)6(H2O)12·6H2O')};      
            end
            
            if((Fe_n+Al_n) / S_n > 0.6 && (Fe_n+Al_n) / S_n < 0.73)
                Fe3 = 1 - (Al_n+Cr_n+Ti_n);
                if Fe3<0, Fe3 = 0; end
                Fe2 = Fe_n - Fe3;
                T.species(i) = {'mikasaite/coquimbite'};
                [sorted_formula1,total1] = Formula_Output({Al_n, Fe3, Cr_n, Ti_n},{'Al','Fe3+','Cr','Ti'});
                [sorted_formula2,total2] = Formula_Output({Fe2, Mg_n,Mn_n,Ca_n},{'Fe2+', 'Mg','Mn','Ca'});
                [sorted_formula3,total3] = Formula_Output({S_n/6, P_n/6},{'S','P'});
                [sorted_formula4,total4] = Formula_Output({Fe_n/2, Mn_n/2, Al_n/2}, {'Fe','Mn','Al'});
                [sorted_formula5,total5] = Formula_Output({S_n/2, P_n/2},{'S','P'});
                T.formula(i) = {strcat(sorted_formula4, total4, '(', sorted_formula5, total5, 'O4)3) or', sorted_formula1,total1,sorted_formula2,total2,'(',sorted_formula3,total3,'O4)6(H2O)12·6H2O')};
            end

            if(Fe_n / S_n > 1-1*t && Fe_n / S_n < 1+1*t)
                T.species(i) = {'Fe-Sulfate (Fe/S = 1, many species possible)'};
            end

            if((Fe_n + Mg_n + Na_n) / S_n >= 0.68 && (Fe_n + Mg_n + Na_n) / S_n < 0.8 && Al_n/Fe_n < 1/2)
                T.species(i) = {'römerite/bílinite'};
                Fe2 = 1 -(Mg_n*(16/24) + Mn_n*(16/24));
                if(Fe2<0), Fe2=0; end
                Fe3 = Fe_n*(16/24) - Fe2;
                if(Fe3<0), Fe3=0; end
                [sorted_formula1, total1] = Formula_Output({Fe2,Mg_n*(16/24),Mn_n*(16/24), Ca_n*(16/24)},{'Fe2+','Mg','Mn','Ca'});
                [sorted_formula2, total2] = Formula_Output({Fe3,Cr_n*(16/24),Ti_n*(16/24),V_n*(16/24), Al_n*(16/24)},{'Fe3+','Cr','Ti','V','Al'});
                [sorted_formula3, total3] = Formula_Output({(S_n*(16/24))/4,(P_n*(16/24))/4},{'S','P'});
                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'[',sorted_formula3,total3,'O4]4·14H2O or ·22H2O')};
            end

            if(round(Fe_n / S_n,2) > (5/6)-(5/6)*t && Fe_n / S_n < (5/6)+(5/6)*t)
                Fe2 = 1 -(Mg_n2*12.5 + Mn_n2*12.5 + Ca_n2*12.5 + Na_n2*12.5);
                if(Fe2<0), Fe2=0; end
                Fe3 = Fe_n2*12.5 - Fe2;
                if(Fe3<0), Fe3=0; end
                OH = 2 - (Cl_n2*12.5 + F_n2*12.5);
                if(OH<0), OH=0; end
                T.species(i) = {'copiapite/ferricopiapite'};
                [sorted_formula1,total1] = Formula_Output({Fe2, Mg_n2*12.5, Mn_n2*12.5, Ca_n2*12.5, Na_n2*12.5},{'Fe2+','Mg','Mn','Ca','Na'});
                [sorted_formula2,total2] = Formula_Output({Fe3,Al_n2*12.5,Cr_n2*12.5},{'Fe3+','Al','Cr'});
                [sorted_formula3,total3] = Formula_Output({S_n2*12.5/6, P_n2*12.5/6},{'S','P'});
                [sorted_formula4,total4] = Formula_Output({OH,Cl_n2*12.5,F_n2*12.5},{'OH','Cl','F'});
                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'[',sorted_formula3,total3,'O4]6',sorted_formula4,total4,'·20H2O')};
            end

            if(Fe_n/(Fe_n+Cu_n+Al_n+Mg_n + Mn_n) > 0.5 && Fe_n >= Cu_n)    
                if((Fe_n + Cu_n) / S_n > 1-1*t && (Fe_n + Cu_n) / S_n < 2.5+2.5*t)        
                    if(Cu_n/Fe_n > 0.2-0.2*t && Cu_n/Fe_n < 1+1*t)
                        if(Fe_n > Mg_n && (Mg_n + Al_n + Mn_n) > 0 && (Mg_n + Al_n + Mn_n) < 0.76+0.76*t)
                            T.species(i) = {'ferrovalleriite'};
                            T.formula(i) = {'2[(Fe,Cu)S]·1.53[(Fe,Al,Mg)(OH)2]'}; 
   
                        end
                    end    
                end  
            end

            if(Mg_n/(Cu_n+Zn_n) > 1 && (Na_n + K_n) <= 0.0025)
                T.species(i) = {'melanterite'};
                [sorted_formula1, total1] = Formula_Output({Fe_n2*2,Cu_n2*2,Ni_n2*2, Co_n2*2, Mg_n2*2, Ca_n2*2, Mn_n2*2, Na_n2*2}, {'Fe','Cu','Ni','Co','Mg','Ca','Mn','Na'});
                [sorted_formula2, total2] = Formula_Output({S_n2*2, P_n2*2},{'S','P'});
                T.formula(i) = {strcat(sorted_formula1, total1, '(', sorted_formula2, total2, 'O4)·5H2O')};
            end

            if((Fe_n + Cu_n + Zn_n) / S_n > 1-1*t && (Fe_n + Cu_n + Zn_n) / S_n < 1+1*t)
                if(Fe_n/(Fe_n + Cu_n + Zn_n) > 0.5 && Mg_n/(Cu_n+Zn_n) < 1)
                    T.species(i) = {'siderotil/szomolnokite'};
                    [sorted_formula1, total1] = Formula_Output({Fe_n2*2,Cu_n2*2,Ni_n2*2, Co_n2*2, Mg_n2*2, Ca_n2*2, Mn_n2*2, Na_n2*2}, {'Fe','Cu','Ni','Co','Mg','Ca','Mn','Na'});
                    [sorted_formula2, total2] = Formula_Output({S_n2*2, P_n2*2},{'S','P'});
                    T.formula(i) = {strcat(sorted_formula1, total1, '(', sorted_formula2, total2, 'O4)·5H2O')};
                end
            end

            if((Fe_n + Cu_n) / S_n > (3/2)-(3/2)*t && (Fe_n + Cu_n) / S_n < (3/2)+(3/2)*t)
                T.species(i) = {'hydroniumjarosite'};
                [sorted_formula1, total1] = Formula_Output({Fe_n/2, Al_n/2, Mg_n/2, Ca_n/2, Mn_n/2, Na_n/2},{'Fe3+','Al','Mg','Ca','Mn','Na'});
                [sorted_formula2, total2] = Formula_Output({S_n/2,P_n/2},{'S','P'});
                [sorted_formula3, total3] = Formula_Output({6-(F_n/4 + Cl_n/4), F_n/4, Cl_n/4}, {'OH','Cl','F'});
                T.formula(i) = {strcat('H3O',sorted_formula1,total1,'(',sorted_formula2,total2,'O4)2', sorted_formula3,total3)};
            end

            if((Fe_n + Cu_n) / S_n > (16/3.2)-(16/3.2)*t && (Fe_n + Cu_n) / S_n < (16/3.2)+(16/3.2)*t)
                T.species(i) = {'schwertmannite'};
                renorm = 33.6/24;
                OH = 9.6 - (F_n*renorm + Cl_n*renorm);
                if(OH<=0), OH = 0; end
                [sorted_formula1, total1] = Formula_Output({Fe_n*renorm, Al_n*renorm, Mg_n*renorm, Mn_n*renorm, Ca_n*renorm, Na_n*renorm},{'Fe','Al','Mg','Mn','Ca','Na'});
                [sorted_formula2, total2] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                [sorted_formula3, total3] = Formula_Output({S_n*renorm/3.2, P_n*renorm/3.2},{'S', 'P'});
                T.formula(i) = {strcat(sorted_formula1,total1, 'O16', sorted_formula2, total2, '(',sorted_formula3, total3, 'O4)3.2·10H2O')};
            end

            if((Fe_n + Pb_n) / S_n > 0.5-0.5*t && (Fe_n + Pb_n) / S_n < 0.5+0.5*t && Pb_n > 0)
                T.species(i) = {'viaeneite'};
                [sorted_formula, total] = Formula_Output({Fe_n2/2,Pb_n2/2,Cu_n2/2, Ni_n2/2, Co_n2/2}, {'Fe','Pb','Cu','Ni','Co'});
                T.formula(i) = {strcat(sorted_formula, total, 'S8O')};
            end

        end %end of Fe-Sulfates

        %% 3.7 Na-K-Fe-Sulfate

        if((Fe_n + Al_n) / S_n > (1/8)-(1/8)*t && (Fe_n + Al_n) / S_n < (2/8)+(2/8)*t)
            if((K_n + Na_n + Ca_n + Pb_n) / S_n > (6/8)-(6/8)*t && (K_n + Na_n + Ca_n + Pb_n) / S_n < (6/8)+(6/8)*t)
                if((K_n + Na_n + Pb_n)/(Mg_n + Cu_n) > (4/3)-(4/3)*t && (K_n + Na_n + Pb_n)/(Mg_n + Cu_n) < (4/2)+(4/2)*t)
                    if(Pb_n > 0 && Cu_n > 0)
                        T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
                        T.group3(i) = {'Na-K-Fe-Sulfate'};
                        T.group4(i) = {''};
                        T.species(i) = {'philoxenite'};
                        Nax = 2 - Ca_n2*16;
                        if(Nax<=0), Nax = 0; end
                        Nay = Na_n2*16 - Nax;
                        if(Nay<=0), Nay = 0; end
                        Fe3 = 1 - Al_n2*16;
                        if(Fe3<=0), Fe3 = 0; end
                        Fe2 = Fe_n2*16 - Fe3;
                        if(Fe2<=0), Fe2 = 0; end  
                        [sorted_formula1, total1] = Formula_Output({K_n2*16,Nay,Pb_n2*16}, {'K','Na','Pb'});
                        [sorted_formula2, total2] = Formula_Output({Nax,Ca_n2*16}, {'Na','Ca'});
                        [sorted_formula3, total3] = Formula_Output({Mg_n2*16,Cu_n2*16, Fe2, Zn_n2*16, Mn_n2*16}, {'Mg','Cu','Fe2+','Zn','Mn'});
                        [sorted_formula4, total4] = Formula_Output({Fe3, Al_n2*16}, {'Fe3+','Al'});
                        [sorted_formula5, total5] = Formula_Output({S_n2*16/8, P_n2*16/8}, {'S','P'});
                        T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, sorted_formula3, total3, sorted_formula4,total4,'(',sorted_formula5,total5,'O4)8')};
                    end
                end
            end
        end %end of philoxenite

        if(Fe_n / (Ca_n + Mg_n + Fe_n + Mn_n) > 0.5 && (Na_n + K_n) / (K_n + Na_n + Mg_n + Mn_n) > 0.5 && (Na_n + K_n)/Fe_n > 0.1)
            T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
            T.group3(i) = {'Na-K-Fe-Sulfate'};
            T.group4(i) = {''};

            if(K_n / Fe_n > 0.3-0.3*t && K_n / Fe_n < 0.33+0.33*t)
                if(K_n > Na_n && Fe_n / S_n > 1.5-1.5*t && Fe_n / S_n < 1.5+1.5*t)
                    T.species(i) = {'jarosite'};
                    renorm = 11/24;
                    vac = 1 - (K_n*renorm + Na_n*renorm);
                    if(vac<0), vac =0; end
                    OH = 6 - (F_n*renorm + Cl_n*renorm);
                    if(OH<0), OH =0; end
                    [sorted_formula1, total1] = Formula_Output({K_n*renorm, Na_n*renorm, vac}, {'K','Na','□'});
                    [sorted_formula2, total2] = Formula_Output({Fe_n*renorm, Al_n*renorm, Mg_n*renorm},{'Fe3+','Al','Mg'});
                    [sorted_formula3, total3] = Formula_Output({S_n*renorm, P_n*renorm},{'S','P'});
                    [sorted_formula4, total4] = Formula_Output({OH/2, F_n*renorm/2, Cl_n*renorm/2}, {'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, '(',sorted_formula3,total3,'O4)2',sorted_formula4,total4)};
                end
            end

            if(Na_n / Fe_n > 1-2*1*t && Na_n / Fe_n < 1+1*t)
                if(Fe_n / S_n > 0.5-2*0.5*t && Fe_n / S_n < 0.5+0.5*t)
                    if(Na_n > K_n)
                        T.species(i) = {'eldfellite/amarillite/erdite'};
                        [sorted_formula1, total1] = Formula_Output({Na_n2*4,K_n2*4,Ca_n2*4}, {'Na','K','Ca'});
                        [sorted_formula2, total2] = Formula_Output({Fe_n2*4,Al_n2*4,Mg_n2*4,Mn_n2*4}, {'Fe3+','Al','Mg','Mn'});
                        [sorted_formula3, total3] = Formula_Output({S_n2*2,P_n2*2}, {'S','P'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'[',sorted_formula3,total3,'O4]2·6H2O')};
                        if((Ca_n) > 0 || (Al_n + Ti_n) > 0)
                           T.species(i) = {'eldfellite/amarillite'};
                        end
                    end
                end
            end

            if(Na_n / Fe_n > 3-3*t && Na_n / Fe_n < 3+3*t)
                if(Fe_n / S_n > (1/3)-(1/3)*t && Fe_n / S_n < (1/3)+(1/3)*t)
                    T.species(i) = {'ferrinatrite'};
                    [sorted_formula1, total1] = Formula_Output({Na_n2*6, K_n2*6, Ca_n2*6}, {'Na','K','Ca'});
                    [sorted_formula2, total2] = Formula_Output({Fe_n2*6, Al_n2*6, Cr_n2*6, V_n2*6}, {'Fe3+','Al','Cr','V'});
                    [sorted_formula3, total3] = Formula_Output({S_n2*6/3, P_n2*6/3}, {'S','P'});
                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2,'[',sorted_formula3,total3,'O4]3·3H2O')};
                end
            end

            if(Na_n / (Fe_n+Al_n) > 0.5-0.5*t && Na_n / (Fe_n+Al_n) < 0.5+0.5*t)
                if((Fe_n+Al_n) / S_n > 0.5-0.5*t && (Fe_n+Al_n) / S_n < 0.5+0.5*t)
                    if(Cl_n / S_n > 0.5/4 && Cl_n / S_n < 0.25+0.25*t)
                        Clx = 1;
                        Cly = Clx - Cl_n*17.5/24;
                        if(Cly<0), Cly = 0; end
                        OH = 2 - (Cl_n*17.5/24 + F_n*17.5/24);
                        if(OH<0), OH = 0; end
                        [sorted_formula1, total1] = Formula_Output({Na_n*17.5/24, K_n*17.5/24, Ca_n*17.5/24}, {'Na','K','Ca'});
                        [sorted_formula2, total2] = Formula_Output({Fe_n*17.5/24, Mg_n*17.5/24, Mn_n*17.5/24}, {'Fe','Mg','Mn'});
                        [sorted_formula3, total3] = Formula_Output({(S_n*17.5/24)/4, (P_n*17.5/24)/4}, {'S','P'});
                        [sorted_formula4, total4] = Formula_Output({OH, Cly, F_n*17.5/24}, {'OH','Cl','F'});
                        T.species(i) = {'adranosite-(Fe)'};
                        T.formula(i) = {strcat('(NH4)4', sorted_formula1, total1, sorted_formula2, total2, '[', sorted_formula3,total3,'O4]4Cl_1', sorted_formula4, total4)};
                    end
                end
            end

            if(Na_n / Fe_n > (1/3)-(1/3)*t && Na_n / Fe_n < (1/3)+(1/3)*t)
                if(Fe_n / S_n > (3/5)-(3/5)*t && Fe_n / S_n < (3/5)+(3/5)*t)
                    T.species(i) = {'coyoteite'};
                    T.formula(i) = {'NaFe3S5·2H2O'};
                end
            end

            if(Na_n / Fe_n > 21-21*t && Na_n / Fe_n < 21+21*t)
                if(Fe_n / S_n > 0.1-(2)*0.1*t && Fe_n / S_n < 0.1+0.1*t)
                    if(Cl_n / S_n > (3/10)-(3/10)*t && Cl_n / S_n < (3/10)+(3/10)*t)
                        T.species(i) = {'D`ansite-(Fe)'};
                        Renorm_factor = 24/41.5;
                        Na_n41 = T.Na_Multi(i) / Renorm_factor;
                        Ca_n41 = T.Ca_Multi(i) / Renorm_factor;
                        Mg_n41 = T.Mg_Multi(i) / Renorm_factor;
                        Fe_n41 = T.Fe_Multi(i) / Renorm_factor;
                        S_n41 = T.S_Multi(i) / Renorm_factor;
                        P_n41 = T.P_Multi(i) / Renorm_factor;
                        Cl_n41 = T.Cl_Multi(i) / Renorm_factor;
                        F_n41 = T.F_Multi(i) / Renorm_factor;
                        Mn_n41 = T.F_Multi(i) / Renorm_factor;
                        OH = 3 - (Cl_n41+F_n41);
                        if(OH<0), OH=0; end
                        [sorted_formula1, total1] = Formula_Output({Na_n41, Ca_n41}, {'Na','Ca'});
                        [sorted_formula2, total2] = Formula_Output({Mg_n41, Fe_n41, Mn_n41}, {'Mg', 'Fe', 'Mn'});
                        [sorted_formula3, total3] = Formula_Output({S_n41/10,P_n41/10}, {'S','P'});
                        [sorted_formula4, total4] = Formula_Output({OH, Cl_n41, F_n41}, {'OH','Cl','F'});
                        T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2,'[',sorted_formula3,total3, 'O4]10',sorted_formula4,total4)};
                    end
                end
            end

            if(Na_n / (Fe_n+Mg_n+Mn_n) > 3-3*t && Na_n / (Fe_n+Mg_n+Mn_n) < 3+3*t)
                if((Fe_n+Mg_n+Mn_n) / S_n > 2-2*t && (Fe_n+Mg_n+Mn_n) / S_n < 2+2*t)
                    if((Na_n + Fe_n + Mg_n + Mn_n + S_n)/3 > 9-9*t && (Na_n + Fe_n + Mg_n + Mn_n + S_n)/3 < 9+9*t)
                        T.species(i) = {'ferrotychite'};
                        Renorm = 8/24;
                        [sorted_formula1,total1] = Formula_Output({Na_n*Renorm, K_n*Renorm, Ca_n*Renorm },{'Na','K','Ca'});
                        [sorted_formula2,total2] = Formula_Output({Fe_n*Renorm, Mg_n*Renorm, Mn_n*Renorm, Sr_n*Renorm, Ba_n*Renorm},{'Fe2+','Mg','Mn','Sr','Ba'});
                        [sorted_formula3,total3] = Formula_Output({S_n*Renorm, P_n*Renorm},{'S','P'});  
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'(CO3)4','(',sorted_formula3,total3, 'O4)')};
                    end
                end
            end

            if(Na_n / Fe_n > 2-2*t && Na_n / Fe_n < 2+2*t)
                if(Fe_n / S_n > 0.5-0.5*t && Fe_n / S_n < 0.5+0.5*t)
                    T.species(i) = {'metasideronatrite'};
                    renorm = 8.5/24;
                    OH = 1 - (F_n*renorm + Cl_n*renorm);
                    if(OH<=0), OH = 0; end
                    [sorted_formula1, total1] = Formula_Output({Na_n*renorm, K_n*renorm},{'Na','K'});
                    [sorted_formula2, total2] = Formula_Output({Fe_n*renorm, Al_n*renorm, Cr_n*renorm, Mg_n*renorm},{'Fe','Al','Cr','Mg'});
                    [sorted_formula3, total3] = Formula_Output({S_n*renorm/2, P_n*renorm/2},{'S','P'});
                    [sorted_formula4, total4] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, '(',sorted_formula3, total3,'O4)2',sorted_formula4, total4, '·H2O')};
                end
            end

            if(K_n / Na_n > (1/3)-(1/3)*t && K_n / Na_n < (1/3)+(1/3)*t)
                if(Fe_n / S_n > (7/12)-(7/12)*t && Fe_n / S_n < (7/12)+(7/12)*t)
                    if((K_n + Na_n) / S_n > (2/3)-(2/3)*t && (K_n + Na_n) / S_n < (2/3)+(2/3)*t)
                        T.species(i) = {'metavoltine'};
                        renorm = 50/24;
                        Fe3 = 6 - (Al_n*renorm + Cr_n*renorm);
                        if(Fe3<=0), Fe3=0; end
                        Fe2 = Fe_n*renorm - Fe3;
                        if(Fe2<=0), Fe2=0; end
                        Nay = 6 - Ca_n*renorm;
                        if(Nay<=0), Nay=0; end
                        Nax = Na_n*renorm - Nay;
                        if(Nax<=0), Nax=0; end
                        [sorted_formula1, total1] = Formula_Output({K_n*renorm, Nax},{'K','Na'});
                        [sorted_formula2, total2] = Formula_Output({Nay, Ca_n*renorm},{'Na', 'Ca'});
                        [sorted_formula3, total3] = Formula_Output({Fe2, Mg_n*renorm, Mn_n*renorm},{'Fe2+','Mg','Mn'});
                        [sorted_formula4, total4] = Formula_Output({Fe3, Al_n*renorm, Cr_n*renorm},{'Fe3+','Al','Cr'});
                        [sorted_formula5, total5] = Formula_Output({S_n*renorm/12, P_n*renorm/12},{'S','P'});
                        T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, sorted_formula4, total4, 'O2(', sorted_formula5, total5,'O4)12·18H2O')};
                    end
                end
            end

            if(Na_n / Fe_n > (1/3)-(1/3)*t && Na_n / Fe_n < (1/3)+(1/3)*t)
                if(Fe_n / S_n > (3/2)-(3/2)*t && Fe_n / S_n < (3/2)+(3/2)*t)
                    T.species(i) = {'natrojarosite'};
                    renorm = 11/24;
                    [sorted_formula1, total1] = Formula_Output({Na_n*renorm, K_n*renorm},{'Na','K'});
                    [sorted_formula2, total2] = Formula_Output({Fe_n*renorm, Al_n*renorm},{'Fe','Al'});
                    [sorted_formula3, total3] = Formula_Output({S_n*renorm, P_n*renorm},{'S','P'});
                    [sorted_formula4, total4] = Formula_Output({6 - (F_n*renorm + Cl_n*renorm)/2, F_n*renorm/2, Cl_n*renorm/2},{'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2, total2, '(',sorted_formula3,total3,'O4)2',sorted_formula4,total4)};
                end
            end

            if(Na_n / Fe_n > (1/6)-2*(1/6)*t && Na_n / Fe_n < (1/6)+(1/6)*t)
                if(Al_n / Fe_n > (1/2)-(1/2)*t && Al_n / Fe_n < (1/2)+(1/2)*t)
                    if(Fe_n / S_n > 3-3*t && Fe_n / S_n < 3+3*3*t)
                        T.species(i) = {'nikischerite'};
                        renorm = 17/24;
                        Fe3 = 3 - Al_n*renorm;
                        if(Fe3<=0), Fe3 = 0; end
                        Fe2 = Fe_n*renorm - Fe3;
                        if(Fe2<=0), Fe2 = 0; end
                        OH = 18 - (F_n*renorm + Cl_n*renorm);
                        if(OH<=0), OH = 0; end
                        [sorted_formula1, total1] = Formula_Output({Fe2, Mg_n*renorm, Mn_n*renorm, Ca_n*renorm},{'Fe2+','Mg','Mn','Ca'});
                        [sorted_formula2, total2] = Formula_Output({Al_n*renorm, Fe3},{'Al','Fe3+'});
                        [sorted_formula3, total3] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                        [sorted_formula4, total4] = Formula_Output({Na_n*renorm, K_n*renorm},{'Na','K'});
                        [sorted_formula5, total5] = Formula_Output({S_n*renorm/2, P_n*renorm/2},{'S','P'});
                        T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, sorted_formula3,total3, '[', sorted_formula4,total4, '(H2O)6](',sorted_formula5,total5,'O4)2·6H2O')};
                    end
                end
            end

            if((Fe_n+Al_n) / S_n > (2/3)-(2/3)*t && Fe_n / S_n < (2/3)+(2/3)*t)
                if(Cl_n / S_n > (5/3)-(5/3)*t && Cl_n / S_n < (5/3)+(5/3)*t)
                    if((Na_n+K_n+Mn_n) / S_n > 1-1*t && (Na_n+K_n+Mn_n) / S_n < 1+2*1*t)
                        if(Na_n > K_n)
                            T.species(i) = {'therasiaite'};
                            renorm = 14.5/24;
                            Nax = 1 - K_n*renorm;
                            if(Nax<=0), Nax = 0; end
                            Nay = Na_n*renorm - Nax;
                            if(Nay<=0), Nay = 0; end
                            Fe3 = 3 - (Al_n*renorm+Cr_n*renorm);
                            if(Fe3<=0), Fe3 = 0; end
                            Fe2 = Fe_n*renorm - Fe3;
                            if(Fe2<=0), Fe2 = 0; end
                            OH = 5 - (F_n*renorm + Cl_n*renorm);
                            if(OH<=0), OH = 0; end
                            [sorted_formula1, total1] = Formula_Output({K_n*renorm, Nax},{'K','Na'});
                            [sorted_formula2, total2] = Formula_Output({Nay, Ca_n*renorm},{'Na','Ca'});
                            [sorted_formula3, total3] = Formula_Output({Fe2, Mg_n*renorm,Mn_n*renorm},{'Fe2+','Mg','Mn'});
                            [sorted_formula4, total4] = Formula_Output({Fe3, Al_n*renorm, Cr_n*renorm},{'Fe3+','Al','Cr'});
                            [sorted_formula5, total5] = Formula_Output({S_n*renorm/3, P_n*renorm/3},{'S','P'});
                            [sorted_formula6, total6] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                            T.formula(i) = {strcat('(NH4)3', sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, sorted_formula4, total4, '(', sorted_formula5, total5, 'O4)3', sorted_formula6, total6)};
                        end
                    end
                end
            end

            if(Fe_n / S_n > (1/6)-(1/6)*t && Fe_n / S_n < (1/6)+2*(1/6)*t)
                if((Na_n + K_n) / S_n > 2-2*t && (Na_n + K_n) / S_n < 2+2*t)
                    if(Na_n > K_n)
                        T.species(i) = {'ungemachite or clinoungemachite'};
                        renorm = 25/24;
                        Nax = 1 - K_n*renorm;
                        if(Nax<=0), Nax = 0; end
                        Nay = Na_n*renorm - Nax;
                        if(Nay<=0), Nay = 0; end
                        [sorted_formula1, total1] = Formula_Output({K_n*renorm, Nax},{'K','Na'});
                        [sorted_formula2, total2] = Formula_Output({Nay, Ca_n*renorm, Mn_n*renorm},{'Na','Ca','Mn'});
                        [sorted_formula3, total3] = Formula_Output({Fe_n*renorm, Al_n*renorm,Cr_n*renorm},{'Fe3+','Al','Cr'});
                        [sorted_formula4, total4] = Formula_Output({S_n*renorm/6, P_n*renorm/6},{'S','P'});
                        T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2,total2, sorted_formula3, total3, '(', sorted_formula4, total4, 'O4)6 (NO3)2·6H2O or (OH)2·10H2O')};
                    end
                end
            end

            if((K_n+Na_n) / (Fe_n+Al_n) > 1-1*t && (K_n+Na_n) / (Fe_n+Al_n) < 1+1*t)
                if(K_n > Na_n && Fe_n > Al_n)
                    if(Fe_n / S_n > 0.5-0.5*t && Fe_n / S_n < 0.5+2*(0.5)*t)
                        T.species(i) = {'yavapaiite'};
                        renorm = 8/24;
                        [sorted_formula1, total1] = Formula_Output({K_n*renorm, Na_n*renorm, Ca_n*renorm, Mn_n*renorm},{'K','Na','Ca','Mn'});
                        [sorted_formula2, total2] = Formula_Output({Fe_n*renorm, Mg_n*renorm, Al_n*renorm, Cr_n*renorm}, {'Fe','Mg','Al','Cr'});
                        [sorted_formula3, total3] = Formula_Output({S_n*renorm/2, P_n*renorm/2}, {'S','P'});
                        T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, '(', sorted_formula3, total3, 'O4)2')};
                    end
                end
            end

        end %end of Na-K-Fe-Sulfates

        %% 3.8 K-Mg Sulfate
        if(Cl_n / S_n > 1-1*t && Cl_n / S_n < 1+1*t)
            if((K_n + Na_n) / (Mg_n + Fe_n) > 1-1*t && (K_n + Na_n) / (Mg_n + Fe_n) < 1+1*t)
                if(K_n > Na_n && Mg_n > Fe_n)
                    T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
                    T.group3(i) = {'K-Mg Sulfate'};
                    T.group4(i) = {''};
                    T.species(i) = {'kainite'};
                    renorm = 4.5/24;
                    OH = 1 - (F_n*renorm + Cl_n*renorm);
                    if(OH<=0), OH = 0; end
                    [sorted_formula1, total1] = Formula_Output({K_n*renorm, Na_n*renorm}, {'K','Na'});
                    [sorted_formula2, total2] = Formula_Output({Mg_n*renorm, Ca_n*renorm, Mn_n*renorm, Fe_n*renorm}, {'Mg','Ca','Mn','Fe'});
                    [sorted_formula3, total3] = Formula_Output({S_n*renorm, P_n*renorm},{'S','P'});
                    [sorted_formula4, total4] = Formula_Output({OH, Cl_n*renorm, F_n*renorm},{'OH','Cl','F'});
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'(', sorted_formula3,total3,'O4)',sorted_formula4,total4,'·2.75H2O')};
                end
            end
        end

        %% 3.9 K-Al Sulfate
        if ((K_n + Na_n + Al_n + Fe_n)/S_n > (4/2)-(4/2)*t && (K_n + Na_n + Al_n + Fe_n)/S_n < (4/2)+(4/2)*t)
            if((K_n + Na_n) / (Al_n + Fe_n) > (1/3)-(1/3)*t && (K_n + Na_n) / (Al_n + Fe_n) < (1/3)+(1/3)*2*t)
                if(K_n > Na_n && Al_n > Fe_n)
                    T.group2(i) = {'K-Al Sulfate'};
                    T.species(i) = {'alunite'};
                    OH = 6 - (Cl_n*11/24 + F_n*11/24);
                    [sorted_formula1, total1] = Formula_Output({Na_n*11/24,K_n2*11/24,Ca_n2*11/24}, {'Na','K','Ca'});
                    [sorted_formula2, total2] = Formula_Output({Fe_n*11/24, Al_n*11/24, Ti_n*11/24, Cr_n*11/24, Mg_n*11/24, Mn_n*11/24}, {'Fe3+','Al','Ti','Cr','Mg','Mn'});
                    [sorted_formula3, total3] = Formula_Output({(S_n*11/24)/2,(P_n*11/24)/2}, {'S','P'});
                    [sorted_formula4, total4] = Formula_Output({OH, Cl_n*11/24, F_n*11/24}, {'OH','Cl','F'});
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'[',sorted_formula3,total3,'O4]2', sorted_formula4, total4)};
                end
            end
        end

        %% 3.10 Cu-Pb Sulfate
        if((Cu_n + Pb_n) / S_n > (7/3)-(7/3)*3*t && (Cu_n + Pb_n) / S_n < (7/3)+(7/3)*3*t)
            if(Cu_n / Pb_n > (2/5)-(2/5)*2*t && Cu_n / Pb_n < (2/5)+(2/5)*2*t)
                T.group2(i) = {'Cu-Pb Sulfate'};
                T.species(i) = {'caledonite'};
                Renorm = 18/24;
                OH = 6-(Cl_n*Renorm + F_n*Renorm);
                [sorted_formula1,total1] = Formula_Output({Cu_n*Renorm,Fe_n*Renorm,Zn_n*Renorm},{'Cu','Fe','Zn'});
                [sorted_formula2,total2] = Formula_Output({Pb_n*Renorm},{'Pb'});
                [sorted_formula3,total3] = Formula_Output({(S_n*Renorm)/3,(P_n*Renorm)/3},{'S','P'});
                [sorted_formula4,total4] = Formula_Output({OH,Cl_n*Renorm,F_n*Renorm},{'OH','Cl','F'});
                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'[',sorted_formula3,total3,'O4]3(CO3)',sorted_formula4,total4)};
            end
        end

        %% 3.11 Cu-Sulfate
        if(T.CuO(i) > 40 && T.SO3(i) > 40)
            T.group2(i) = {'Cu-Pb Sulfate'};
            T.group3(i) = {'Cu-Sulfate'};
            T.group4(i) = {''};
            T.species(i) = {'chalcocyanite/dolerophanite or Cu-hydroxyl-Sulfate'};
            [sorted_formula1,total1] = Formula_Output({Cu_n/6,Fe_n/6,Zn_n/6,Mg_n/6,Ca_n/6},{'Cu','Fe','Zn','Mg','Ca'});
            [sorted_formula2,total2] = Formula_Output({S_n/6,P_n/6},{'S','P'});
            T.formula(i) = {strcat(sorted_formula1, total1, '[', sorted_formula2, total2, 'O4]')};
        end

        %% 3.12 Extra Sulfates
        if((K_n+Na_n)/(Ca_n+Mn_n+Fe_n) >= (2/5)-(2/5)*t && (K_n+Na_n)/(Ca_n+Mn_n+Fe_n) < (2/5)+(2/5)*t)
            if((K_n+Na_n+Ca_n+Mn_n+Fe_n)/(S_n+P_n) >= (7/6)-(7/6)*t && (K_n+Na_n+Ca_n+Mn_n+Fe_n)/(S_n+P_n) < (7/6)+(7/6)*t && S_n/(S_n+P_n) > 0.5)
                if(K_n > Na_n && Ca_n > Mn_n && Ca_n > Fe_n)
                    T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
                    T.group3(i) = {'Na-K-Fe-Sulfate'};
                    T.group4(i) = {''};
                    T.species(i) = {'görgeyite'};
                    T.formula(i) = {strcat('K_', num2str(round(K_n2*12,2)), 'Ca_', num2str(round(Ca_n2*12,2)), '(SO4)6·H2O')};
                end
            end
        end
        
        if((K_n+Na_n)/(Ca_n+Mn_n+Fe_n) >= (2/1)-(2/1)*t && (K_n+Na_n)/(Ca_n+Mn_n+Fe_n) < (2/1)+(2/1)*t)
            if((K_n+Na_n+Ca_n+Mn_n+Fe_n)/(S_n+P_n) >= (3/2)-(3/2)*t && (K_n+Na_n+Ca_n+Mn_n+Fe_n)/(S_n+P_n) < (3/2)+(3/2)*t && S_n/(S_n+P_n) > 0.5)
                if(K_n > Na_n && Ca_n > Mn_n && Ca_n > Fe_n)
                    T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
                    T.group3(i) = {'Na-K-Fe-Sulfate'};
                    T.group4(i) = {''};
                    T.species(i) = {'syngenite'};
                    T.formula(i) = {strcat('K_', num2str(round(K_n2*4,2)), 'Ca_', num2str(round(Ca_n2*4,2)), '(SO4)2·H2O')};
                end
            end
        end

        if((Mg_n+Fe_n+Mn_n+Al_n)/(S_n+P_n) >= (5/6)-2*(5/6)*t && (Mg_n+Fe_n+Mn_n+Al_n)/(S_n+P_n) < (5/6)+2*(5/6)*t && S_n/(S_n+P_n) > 0.5)
            if((Mg_n+Mn_n)/(Fe_n+Al_n) >= (1/4)-3*(1/4)*t && (Mg_n+Mn_n)/(Fe_n+Al_n) < (1/4)+3*(1/4)*t)
                if(Fe_n>Al_n && Mg_n>Mn_n)
                    T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};    
                    T.group3(i) = {'Fe-Sulfate'};
                    T.group4(i) = {''};
                    T.species(i) = {'magnesiocopiapite'};
                    T.formula(i) = {strcat('Mg_', num2str(round(Mg_n2*13,2)), 'Fe3+4(SO4)6(OH)2·20H2O')};
                end
            end
        end

        if((Mg_n+Fe_n+Mn_n+Al_n)/(S_n+P_n) >= (1)-2*(1)*t && (Mg_n+Fe_n+Mn_n+Al_n)/(S_n+P_n) < (1)+2*(1)*t && S_n/(S_n+P_n) > 0.5)
            if((Mg_n+Mn_n)/(Fe_n+Al_n) >= (1)-2*(1)*t && (Mg_n+Mn_n)/(Fe_n+Al_n) < (1)+2*(1)*t)
                if(Fe_n>Al_n && Mg_n>Mn_n)
                    T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
                    T.species(i) = {'botryogen'};
                    T.formula(i) = {strcat('Mg_', num2str(round(Mg_n2*4,2)), 'Fe3+(SO4)2(OH)·7H2O')};
                end
            end
        end

        if((Mg_n+Ca_n+Mn_n+Na_n+K_n)/(Fe_n+Al_n) < 0.1 && S_n/(S_n+P_n) > 0.5)
            if((Fe_n+Al_n)/(S_n+P_n) >= 1-1*t && (Fe_n+Al_n)/(S_n+P_n) < 1+1*t && Cl_n/S_n >= 0.5)
                T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
                T.group3(i) = {'Fe-Sulfate'};
                T.group4(i) = {''};
                T.species(i) = {'xitieshanite'};
                T.formula(i) = {strcat('Fe3+_', num2str(round(Fe_n2*2,2)), '(SO4)Cl·6H2O')};
            end
        end

    end %end of Sulfates group


    %% 4. Phosphate Class.
    if(strcmp(T.group1(i), 'Phosphate') && P_n > S_n)

        %% 4.1 Ca-Phosphate
        if(Ca_n / (ABCD_n) > 0.5 && Na_n / Ca_n < 0.1) 
            T.group2(i) = {'Na-Ca-Phosphate'};
            T.group3(i) = {'Ca-Phosphate'};
            T.group4(i) = {''};
            T.species(i) = {''};
            T.formula(i) = {''};
    
            if(Ca_n / (P_n+S_n) > (5/3)-(5/3)*t && Ca_n / (P_n+S_n) < (5/3)+(5/3)*t)
                T.group4(i) = {'Apatite'};
                Renorm_factor = 24 / 12.5;
                F_n1 = T.F_Multi(i) / Renorm_factor;
                P_n1 = T.P_Multi(i) / Renorm_factor;
                Cl_n1 = T.Cl_Multi(i) / Renorm_factor;
                Ca_n1 = T.Ca_Multi(i) / Renorm_factor;
                Na_n1 = T.Na_Multi(i) / Renorm_factor;
                Mg_n1 = T.Mg_Multi(i) / Renorm_factor;
                La_n1 = T.La_Multi(i) / Renorm_factor;
                Ce_n1 = T.Ce_Multi(i) / Renorm_factor;
                Nd_n1 = T.Nd_Multi(i) / Renorm_factor;
                S_n1 = T.S_Multi(i) / Renorm_factor;
    
                if((F_n1 + Cl_n1) < 0.5+0.5*t)
                    T.species(i) = {'hydroxylapatite'};
                    OH = 1 - (F_n1/3 + Cl_n1/3);
                    if(OH<=0), OH=0; end
                    [sorted_formula1, total1] = Formula_Output({Ca_n1, Na_n1, Mg_n1},{'Ca','Na','Mg'});
                    [sorted_formula2, total2] = Formula_Output({P_n1, S_n1},{'P','S'});
                    [sorted_formula3, total3] = Formula_Output({OH, F_n1/3, Cl_n1/3},{'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1,total1, '(',sorted_formula2,total2, 'O4)3',sorted_formula3,total3)};
                end
    
                if(F_n1 > Cl_n1 && (F_n1+Cl_n1)>= 0.5)
                    T.species(i) = {'fluorapatite'};
                    OH = 1 - (F_n1 + Cl_n1);
                    if(OH<=0), OH=0; end
                    [sorted_formula1, total1] = Formula_Output({Ca_n1, Na_n1, Mg_n1},{'Ca','Na','Mg'});
                    [sorted_formula2, total2] = Formula_Output({P_n1/3, S_n1/3},{'P','S'});
                    [sorted_formula3, total3] = Formula_Output({OH, F_n1, Cl_n1},{'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1, total1, '(',sorted_formula2,total2, 'O4)3',sorted_formula3,total3)};
                end
    
                if(Cl_n1 > F_n1 && (F_n1+Cl_n1)>= 0.5)
                    OH = 1-(F_n1+Cl_n1);
                    if(OH<=0), OH=0; end
                    T.species(i) = {'chlorapatite'};
                    [sorted_formula1,total1] = Formula_Output({Ca_n1, Na_n1, Mg_n1, La_n1, Ce_n1, Nd_n1},{'Ca','Na','Mg','La','Ce','Nd'});
                    [sorted_formula2,total2] = Formula_Output({P_n1, S_n1},{'P','S'});
                    [sorted_formula3,total3] = Formula_Output({Cl_n1, F_n1, OH},{'Cl','F','OH'});
                    T.formula(i) = {strcat(sorted_formula1,total1,'[',sorted_formula2,total2,'O4]3',sorted_formula3,total3)};
                end
            end
    
            if((Ca_n+Mg_n+Na_n+K_n+Ni_n+Fe_n) / (P_n+S_n) > (3/2)-(3/2)*t && (Ca_n+Mg_n+Na_n+K_n+Ni_n+Fe_n) / (P_n+S_n) < (3/2)+(3/2)*t && Cl_n < 0.001)
                T.species(i) = {'tuite'};
                renorm = 8/24;
                [sorted_formula1, total1] = Formula_Output({Ca_n*renorm, Mg_n*renorm, Na_n*renorm},{'Ca','Mg','Na'});
                [sorted_formula2, total2] = Formula_Output({P_n*renorm/2, S_n*renorm/2}, {'P','S'});
                T.formula(i) = {strcat(sorted_formula1, total1, '(',sorted_formula2, total2, 'O4)2')};
            end
    
            if((Ca_n+Th_n+Ce_n) / (P_n+S_n) > 1-1*t && (Ca_n+Th_n+Ce_n) / (P_n+S_n) < 1+1*t && (Th_n + Ce_n) > 0.1)
                if(Ca_n > Th_n && Ca_n > Ce_n)
                    T.species(i) = {'brockite'};
                    Renorm_factor = 6;
                    Ca_n4 = T.Ca_Multi(i) / Renorm_factor;
                    P_n4 = T.P_Multi(i) / Renorm_factor;
                    S_n4 = T.S_Multi(i) / Renorm_factor;
                    Th_n4 = T.Th_Multi(i) / Renorm_factor;
                    Ce_n4 = T.Ce_Multi(i) / Renorm_factor;
                    La_n4 = T.La_Multi(i) / Renorm_factor;
                    Nd_n4 = T.Nd_Multi(i) / Renorm_factor;
                    [sorted_formula1,total1] = Formula_Output({Ca_n4,Th_n4,Ce_n4,La_n4,Nd_n4},{'Ca','Th','Ce','La','Nd'});
                    [sorted_formula2,total2] = Formula_Output({P_n4,S_n4},{'P','S'});
                    T.formula(i) = {strcat(sorted_formula1,total1,'[',sorted_formula2,total2,'O4]·H2O')};
                end
            end
    
            if(Ca_n / (P_n+S_n) > 1-1*t && Ca_n / (P_n+S_n) < 1+1*t && (Th_n + Ce_n) == 0)
                Renorm = 3.5/24;
                T.species(i) = {'monetite or brushite'};
                OH = 1 - (F_n*Renorm+Cl_n*Renorm);
                if(OH<=0), OH = 0; end
                [sorted_formula1,total1] = Formula_Output({Ca_n*Renorm,Mn_n*Renorm,Mg_n*Renorm,Fe_n*Renorm,La_n*Renorm,Ce_n*Renorm,Nd_n*Renorm},{'Ca','Mn','Mg','Fe','La','Ce','Nd'});
                [sorted_formula2,total2] = Formula_Output({P_n*Renorm,S_n*Renorm},{'P','S'});
                [sorted_formula3,total3] = Formula_Output({OH, F_n*Renorm, Cl_n*Renorm},{'OH','F','Cl'});
                T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O3', sorted_formula3, total3, 'or', sorted_formula1,total1,'[',sorted_formula2, total2, 'O3OH]·2H2O')};
            end
    
            if(Ca_n / (P_n+S_n) > 2-2*t && Ca_n / (P_n+S_n) < 2+2*t && (Ca_n + P_n + S_n) > 3*(24/4.5)-3*(24/4.5)*t && (Ca_n + P_n + S_n) < 3*(24/4.5)+3*(24/4.5)*t)
                T.species(i) = {'isoclasite'};
                Renorm_factor = 24 / 4.5;
                Ca_n5 = T.Ca_Multi(i) / Renorm_factor;
                P_n5 = T.P_Multi(i) / Renorm_factor;
                Mg_n5 = T.Mg_Multi(i) / Renorm_factor;
                Mn_n5 = T.Mn_Multi(i) / Renorm_factor;
                Na_n5 = T.Na_Multi(i) / Renorm_factor;
                S_n5 = T.S_Multi(i) / Renorm_factor;
                F_n5 = T.F_Multi(i) / Renorm_factor;
                Cl_n5 = T.Cl_Multi(i) / Renorm_factor;
                OH = 1 - (F_n5 + Cl_n5);
                if(OH<=0), OH=0; end
                [sorted_formula1, total1] = Formula_Output({Ca_n5, Mg_n5, Mn_n5, Na_n5},{'Ca','Mg','Mn','Na'});
                [sorted_formula2, total2] = Formula_Output({P_n5, S_n5},{'P','S'});
                [sorted_formula3, total3] = Formula_Output({OH, F_n5, Cl_n5},{'OH','F','Cl'});
                T.formula(i) = {strcat(sorted_formula1,total1, '(',sorted_formula2,total2,'O4)',sorted_formula3,total3,'·2H2O')};
            end
    
        end %end of Ca-Phosphates
    
        %% 4.2 Mg-Phosphate
        if(Mg_n / (Na_n + Ca_n + Mg_n + Fe_n + Mn_n) > 0.5)
            T.group2(i) = {'Mg-Phosphate'};
            T.group3(i) = {''};
            T.group4(i) = {''};
            T.species(i) = {''};
            T.formula(i) = {''};
    
            if((Mg_n + Fe_n + Mn_n) / (P_n+S_n) > (3/2)-2*(3/2)*t && (Mg_n + Fe_n + Mn_n) / (P_n+S_n) < (3/2)+2*(3/2)*t && Mg_n > Fe_n && Mg_n > Mn_n)
                T.species(i) = {'dry: chopinite or ferringtonite, or OH-bearing: barićite, bobierrite, cattiite'};
                [sorted_formula1, total1] = Formula_Output({Mg_n/3,Fe_n/3,Ca_n/3,Mn_n/3,Na_n/3,K_n/3},{'Mg','Fe','Ca','Mn','Na','K'});
                [sorted_formula2, total2]  = Formula_Output({(P_n/3)/2, (S_n/3)/2},{'P','S'});
                T.formula(i) = {strcat(sorted_formula1,total1,'[',sorted_formula2,total2,'O4]2 / ·8H2O / ·22H2O')};
            end
    
            if((Mg_n + Fe_n + Mn_n) / (P_n+S_n) > (12/6)-(12/6)*t && (Mg_n + Fe_n + Mn_n) / (P_n+S_n) < (12/6)+(12/6)*t && Mg_n > Fe_n && Mg_n > Mn_n)
                T.species(i) = {'holtedahlite or hydroxylwagnerite or kovdorskite'};
                renorm = 4.5/24;
                Mg_n4 = T.Mg_Multi(i)*renorm;
                Ca_n4 = T.Ca_Multi(i)*renorm;
                Mn_n4 = T.Mn_Multi(i)*renorm;
                Na_n4 = T.Na_Multi(i)*renorm;
                K_n4 = T.K_Multi(i)*renorm;
                Fe_n4 = T.Fe_Multi(i)*renorm;
                P_n4 = T.P_Multi(i)*renorm;
                S_n4 = T.S_Multi(i)*renorm;
                F_n4 = T.F_Multi(i)*renorm;
                Cl_n4 = T.Cl_Multi(i)*renorm;
                OH = 1 - (F_n4+Cl_n4);
                if(OH<=0), OH=0; end
                [sorted_formula1, total1] = Formula_Output({Mg_n4, Ca_n4, Mn_n4, Na_n4, K_n4, Fe_n4},{'Mg','Ca','Mn','Na','K','Fe'});
                [sorted_formula2, total2] = Formula_Output({P_n4, S_n4}, {'P','S'});
                [sorted_formula3, total3] = Formula_Output({OH, F_n4, Cl_n4},{'OH','F','Cl'});
                T.formula(i) = {strcat(sorted_formula1,total1, '(',sorted_formula2,total2, 'O4)',sorted_formula3,total3, '·nH2O')};
            end
    
            if((Mg_n + Fe_n + Mn_n) / (P_n+S_n) > 1-1*t && (Mg_n + Fe_n + Mn_n) / (P_n+S_n) < 1+1*t && Mg_n > Fe_n && Mg_n > Mn_n)
                T.species(i) = {'newberyite or phosphorrösslerite'};
                renorm = 3.5/24;
                [sorted_formula1,total1] = Formula_Output({Mg_n*renorm, Ca_n*renorm, Fe_n*renorm, Mn_n*renorm},{'Mg','Ca','Fe','Mn'});
                [sorted_formula2,total2] = Formula_Output({S_n*renorm, P_n*renorm}, {'S','P'});
                T.formula(i) = {strcat(sorted_formula1,total1,'(', sorted_formula2,total2, 'O3OH)·3H2O or ·7H2O')};
            end
    
            if((Mg_n + Fe_n + Mn_n) / (P_n+S_n) > (7/2)-(7/2)*t && (Mg_n + Fe_n + Mn_n) / (P_n+S_n) < (7/2)+(7/2)*t && Mg_n > Fe_n && Mg_n > Mn_n)
                T.species(i) = {'raaedite'};
                OH = 8 - (F_n/2 + Cl_n/2);
                if(OH<=0), OH = 0; end
                [sorted_formula1, total1] = Formula_Output({Mg_n/2, Fe_n/2, Mn_n/2, Ca_n/2, Na_n/2},{'Mg','Fe','Mn','Ca','Na'});
                [sorted_formula2, total2] = Formula_Output({P_n/2, S_n/2}, {'P','S'});
                [sorted_formula3, total3] = Formula_Output({OH/2, F_n/2/2, Cl_n/2/2}, {'OH','F','Cl'});
                T.formula(i) = {strcat(sorted_formula1, total1, '(', sorted_formula2, total2, 'O4)2',sorted_formula3, total3)};
            end
    
        end %end of Mg-Phosphates
    
        %% 4.3 Fe-Phosphate
        if(Fe_n / (Na_n + Ca_n + Mg_n + Fe_n + Mn_n) >= 0.49)
            T.group2(i) = {'Fe-Phosphate'};
            T.group3(i) = {''};
            T.group4(i) = {''};
            T.species(i) = {''};
            T.formula(i) = {''};
    
            if((Fe_n + Mn_n + Mg_n+Ca_n+Na_n) / (P_n+S_n) > 1.5-2*1.5*t && (Fe_n + Mn_n + Mg_n+Ca_n+Na_n) / (P_n+S_n) <  1.5+2*1.5*t)
                if(Fe_n > Mn_n && Fe_n > Ca_n && Fe_n > Mg_n)
                    renorm = 8/24;
                    T.species(i) = {'graftonite or sarcopside'};
                    [sorted_formula1, total1] = Formula_Output({Fe_n*renorm, Mg_n*renorm, Mn_n*renorm, Ca_n*renorm, Na_n*renorm},{'Fe','Mg','Mn','Ca','Na'});
                    [sorted_formula2, total2] = Formula_Output({S_n*renorm/2, P_n*renorm/2}, {'S','P'});
                    T.formula(i) = {strcat('Fe2+Fe2+2(PO4)2 or ', sorted_formula1, total1, '(',sorted_formula2, total2, 'O4)2')};
                end
            end
    
            if((Fe_n+Mn_n+Mg_n+Ca_n+Na_n) / (P_n+S_n) > 3-2*3*t && (Fe_n+Mn_n+Mg_n+Ca_n+Na_n) / (P_n+S_n) < 3+2*3*t)
                if(Fe_n/(Fe_n+Mn_n+Mg_n+Ca_n+Na_n) > 0.5)
                    Renorm_factor = 24 / 7;
                    Fe_n7 = T.Fe_Multi(i) / Renorm_factor;
                    P_n7 = T.P_Multi(i) / Renorm_factor;
                    Al_n7 = T.Al_Multi(i) / Renorm_factor; 
                    S_n7 = T.S_Multi(i) / Renorm_factor;
                    T.species(i) = {'grattarolaite'};
                    [sorted_formula1, total1] = Formula_Output({Fe_n7, Al_n7},{'Fe3+', 'Al'});
                    [sorted_formula2, total2] = Formula_Output({P_n7, S_n7},{'P', 'S'});
                    T.formula(i) = {strcat(sorted_formula1,total1, 'O3','(',sorted_formula2,total2,'O4)')};
                end
            end
    
            if((Fe_n+Mn_n+Mg_n+Ca_n+Na_n) / (P_n+S_n) > 1-2*1*t && (Fe_n+Mn_n+Mg_n+Ca_n+Na_n) / (P_n+S_n) < 1+2*1*t)
                if(Fe_n > Mn_n && Fe_n > Ca_n)
                    Renorm_factor = 6;
                    Fe_n4 = T.Fe_Multi(i) / Renorm_factor;
                    P_n4 = T.P_Multi(i) / Renorm_factor;
                    S_n4 = T.S_Multi(i) / Renorm_factor;
                    Mn_n4 = T.Mn_Multi(i) / Renorm_factor;
                    Al_n4 = T.Al_Multi(i) / Renorm_factor;
                    Mg_n4 = T.Mg_Multi(i) / Renorm_factor;
                    Ca_n4 = T.Ca_Multi(i) / Renorm_factor;
                    Na_n4 = T.Na_Multi(i) / Renorm_factor;
                    T.species(i) = {'heterosite or rodolicoite'};
                    [sorted_formula1,total1] = Formula_Output({Fe_n4, Mn_n4, Al_n4},{'Fe3+','Mn3+','Al'});
                    [sorted_formula2,total2] = Formula_Output({P_n4, S_n4},{'P','S'});
                    [sorted_formula3, total3] = Formula_Output({Fe_n4, Mg_n4, Ca_n4, Na_n4}, {'Fe3+','Mg','Ca','Na'});
                    T.formula(i) = {strcat(sorted_formula1, total1, '(',sorted_formula2,total2,'O4) or ', sorted_formula3, total3, '(',sorted_formula2,total2,'O4)')};
                end
            end
    
        end %end of Fe-Phosphates
    
        %% 4.4 Na-Phosphates
        if(Na_n / (K_n+Na_n+Ca_n+Mg_n+Fe_n+Mn_n) > 0.5 && (Na_n+K_n) / (P_n+S_n) > 1.95-2*(1.95)*t && (Na_n+K_n) / (P_n+S_n) < (7/2) + 2*(7/2)*t)
            T.group2(i) = {'Na-Ca-Phosphate'};
            T.group3(i) = {'Na-Phosphate'};
            T.group4(i) = {''};
            T.species(i) = {''};
            T.formula(i) = {''};
    
            if((Na_n+K_n+Ca_n)/(P_n+S_n) > 1.95-1.95*t && (Na_n+K_n+Ca_n)/(P_n+S_n) < 2 + 2*2*t && Na_n/(Na_n+K_n+Ca_n) >= 0.5) 
                T.group4(i) = {'Na-Phosphate-OH'};
                Renorm_factor = 24 / 3.5;
                Na_n4 = T.Na_Multi(i) / Renorm_factor;
                K_n4 = T.K_Multi(i) / Renorm_factor;
                Ca_n4 = T.Ca_Multi(i) / Renorm_factor;
                Mn_n4 = T.Mn_Multi(i) / Renorm_factor;
                Mg_n4 = T.Mg_Multi(i) / Renorm_factor;
                Fe_n4 = T.Fe_Multi(i) / Renorm_factor;
                P_n4 = T.P_Multi(i) / Renorm_factor;
                S_n4 = T.S_Multi(i) / Renorm_factor;
                F_n4 = T.F_Multi(i) / Renorm_factor;
                Cl_n4 = T.Cl_Multi(i) / Renorm_factor;
                T.species(i) = {'catalanoite or dorfmanite or nahpoite'};
                [sorted_formula1,total1] = Formula_Output({Na_n4,K_n4,Ca_n4,Mn_n4,Mg_n4,Fe_n4},{'Na','K','Ca','Mn','Mg','Fe'});
                [sorted_formula2,total2] = Formula_Output({P_n4,S_n4},{'P','S'});
                [sorted_formula3,total3] = Formula_Output({1 - (F_n4+Cl_n1), F_n4, Cl_n4},{'OH','F','Cl'});
                T.formula(i) = {strcat(sorted_formula1,total1,'[H',sorted_formula2,total2,'O4]·8H2O or', sorted_formula1, total1, sorted_formula2,total2, 'O3',sorted_formula3, total3,'·2H2O or',sorted_formula1,total1,sorted_formula2,total2,'O3',sorted_formula4,total4)};
            end
    
            if((Na_n+K_n)/(P_n+S_n) > (7/2)-2*(7/2)*t && (Na_n+K_n)/(P_n+S_n) < (7/2)+2*(7/2)*t && (F_n+Cl_n)/(P_n+S_n) > (0.5/2)-(0.5/2)*t && (F_n+Cl_n)/(P_n+S_n) <= 1 + 2*(1)*t)
                if(F_n > Cl_n)
                    T.group4(i) = {'Na-Phosphate-F'};
                    T.species(i) = {'natrophosphate'};
                    renorm = 8.5/24;
                    [sorted_formula1, total1] = Formula_Output({Na_n*renorm, K_n*renorm, Ca_n*renorm, Mg_n*renorm},{'Na','K', 'Ca','Mg'});
                    [sorted_formula2, total2] = Formula_Output({S_n*renorm, P_n*renorm},{'S','P'});
                    [sorted_formula3, total3] = Formula_Output({1 - (F_n*renorm + Cl_n*renorm)/2, F_n*renorm/2, Cl_n*renorm/2},{'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1,total1, '(',sorted_formula2,total2,'O4)2',sorted_formula3,total3,'·19H2O')};
                end
            end
    
        end %end of Na-Phosphates
    
        %% 4.5 Na-Ca-Phosphate
        if((Na_n+Ca_n) > (Mg_n+Fe_n+Mn_n) && Na_n/Ca_n > 1/9)
    
            if(Na_n / Ca_n > 0.95-0.95*t && Na_n / Ca_n < 1.05+1.05*t)
                if((Na_n + Ca_n) / (P_n+S_n) > 1.95-1.95*t && (Na_n + Ca_n) / (P_n+S_n) < 2.05+2.05*t)
                    T.group2(i) = {'Na-Ca-Phosphate'};
                    T.group3(i) = {''};
                    T.group4(i) = {''};
                    T.species(i) = {''};
                    T.formula(i) = {''};
                    Renorm_factor = 6;
                    Na_n4 = T.Na_Multi(i) / Renorm_factor;
                    K_n4 = T.K_Multi(i) / Renorm_factor;
                    Ca_n4 = T.Ca_Multi(i) / Renorm_factor;
                    Mg_n4 = T.Mg_Multi(i) / Renorm_factor;
                    Fe_n4 = T.Fe_Multi(i) / Renorm_factor;
                    Mn_n4 = T.Mn_Multi(i) / Renorm_factor;
                    P_n4 = T.P_Multi(i) / Renorm_factor;
                    S_n4 = T.S_Multi(i) / Renorm_factor;
                    T.species(i) = {'buchwaldite'};
                    [sorted_formula1, total1] = Formula_Output({Na_n4,K_n4}, {'Na','K'});
                    [sorted_formula2, total2] = Formula_Output({Ca_n4,Mg_n4,Fe_n4,Mn_n4}, {'Ca','Mg','Fe','Mn'});
                    [sorted_formula3, total3] = Formula_Output({P_n4,S_n4}, {'P','S'});  
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'[',sorted_formula3,total3,'O4]')};
                end
            end
    
            if (Na_n/Mg_n > 1-1*t && Na_n/Mg_n < 1+1*t && Na_n/Ca_n > (1/9)-(1/9)*t && Na_n/Ca_n < (1/9)+(1/9)*t)
                T.group2(i) = {'Na-Ca-Phosphate'};
                T.group3(i) = {'Ca-Phosphate'};
                T.group4(i) = {''};
                T.species(i) = {'merrillite'};
                renorm = 28/24;
                Ca_n28 = T.Ca_Multi(i)*renorm;
                Mn_n28 = T.Mn_Multi(i)*renorm;
                Na_n28 = T.Na_Multi(i)*renorm;
                Mg_n28 = T.Mg_Multi(i)*renorm;
                Fe_n28 = T.Fe_Multi(i)*renorm;
                P_n28 = T.P_Multi(i)*renorm;
                S_n28 = T.S_Multi(i)*renorm;
                Cay = 1 - Na_n28;
                if(Cay<=0), Cay = 0; end
                Cax = Ca_n28 - Cay;
                if(Cax<=0), Cax = 0; end
                [sorted_formula1, total1] = Formula_Output({Cax, Mn_n28},{'Ca','Mn'});
                [sorted_formula2, total2] = Formula_Output({Na_n28, Cay},{'Na','Ca'});
                [sorted_formula3, total3] = Formula_Output({Mg_n28, Fe_n28},{'Mg','Fe'});
                [sorted_formula4, total4] = Formula_Output({P_n28/7, S_n28/7},{'P','S'});
                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3, total3, '(',sorted_formula4, total4,'O4)7')};
            end
    
        end
    
        %% 4.6 K-Phosphate
        if(K_n / P_n > 0.2-0.2*t && K_n / P_n < 1+1*t && K_n/(K_n + Na_n + Ca_n + Mg_n + Mn_n + Fe_n) > 0.5-2*0.5*t)
            Renorm_factor = 6;
            K_n4 = T.K_Multi(i) / Renorm_factor;
            Na_n4 = T.Na_Multi(i) / Renorm_factor; 
            Ca_n4 = T.Ca_Multi(i) / Renorm_factor;
            Mn_n4 = T.Mn_Multi(i) / Renorm_factor;
            Mg_n4 = T.Mg_Multi(i) / Renorm_factor;
            P_n4 = T.P_Multi(i) / Renorm_factor;
            S_n4 = T.S_Multi(i) / Renorm_factor;
            T.group2(i) = {'K-Phosphate'};
            T.group3(i) = {''};
            T.group4(i) = {''};
            T.species(i) = {'archerite'};
            [sorted_formula1, total1] = Formula_Output({K_n4,Na_n4,Ca_n4,Mn_n4,Mg_n4},{'K','Na','Ca','Mn','Mg'});
            [sorted_formula2, total2] = Formula_Output({P_n4,S_n4},{'P','S'});
            T.formula(i) = {strcat('H2',sorted_formula1,total1,'[',sorted_formula2,total2,'O4]')};
        end
    
        %% 4.7 REE-Phosphate
        if(strcmp(T.group2(i), 'REE-Phosphate'))
            REE = La_n + Ce_n + Nd_n;
    
            if(Ce_n / REE > 0.5)
                T.species(i) = {'monazite-Ce'};
            end
    
            if(La_n / REE > 0.5)
                T.species(i) = {'monazite-La'};
            end
    
            if(Nd_n / REE > 0.5)
                T.species(i) = {'monazite-Nd'};
            end
    
            Renorm_factor = 6;
            La_n4 = T.La_Multi(i) / Renorm_factor;
            Ce_n4 = T.Ce_Multi(i) / Renorm_factor;
            Nd_n4 = T.Nd_Multi(i) / Renorm_factor;
            P_n4 = T.P_Multi(i) / Renorm_factor;
            S_n4 = T.S_Multi(i) / Renorm_factor;
            [sorted_formula1, total1] = Formula_Output({La_n4, Ce_n4, Nd_n4}, {'La', 'Ce', 'Nd'});
            [sorted_formula2, total2] = Formula_Output({P_n4, S_n4},{'P','S'});
            T.formula(i) = {strcat(sorted_formula1, total1, '(', sorted_formula2, total2, 'O4)')};
        end
    
        %% 4.8 Y-Phosphate

        if(strcmp(T.group3(i), 'Y-Phosphate')) 
            T.group4(i) = {''};
            T.species(i) = {'xenotime-Y'};
            renorm = 4/24;
            [sorted_formula1, total1] = Formula_Output({Y_n*renorm, La_n*renorm, Ce_n*renorm, Nd_n*renorm},{'Y','La','Ce','Nd'});
            [sorted_formula2, total2] = Formula_Output({P_n*renorm, S_n*renorm}, {'P','S'});
            T.formula(i) = {strcat(sorted_formula1, total1, '(',sorted_formula2,total2, 'O4)')};
        end

    end %end of Phosphates group

    %% 5 Oxides and Hydroxides
    if(strcmp(T.group1(i), 'Oxide-Hydroxide') || strcmp(T.group1(i), 'Carbonate or Oxide-Hydroxide'))

        %% 5.1.1 Fe-Al-Ti-Mn-Mg Oxides and Hydroxides

        if(T.CuO(i) > 95)
            T.group2(i) = {'Oxide or Hydroxide'};
            T.group3(i) = {'Other Oxide or Hydroxide'};
            T.group4(i) = {''};
            T.species(i) = {'cuprite/tenorite/spertiniite'};
            [sorted_formula, total] = Formula_Output({Cu_n/24, Fe_n/24, Zn_n/24, Mg_n/24},{'Cu','Fe','Zn', 'Mg'});
            T.formula(i) = {strcat(sorted_formula,total,'O or (OH)2')};
            T.check(i) = {'no check'};
        end

        %% 5.1.1.1 Spinel
        if(Si_n < 0.1 && P_n < 0.1 && S_n < 0.1 && Na_n < 0.1)

            if (T.FeO(i) > 98 && T.FeO(i) ~= 100 && Ca_n == 0)
                T.group2(i) = {'Oxide'};
                T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide'};
                T.group4(i) = {'Fe-Oxide'};
                T.species(i) = {'magnetite/hematite/wüstite/Fe-Hydroxide'};
                T.formula(i) = {''};
                T.check(i) = {'no check'};

                if((T.TiO2(i) > 0 && T.TiO2(i) < 1) || (T.Cr2O3(i) > 0 && T.Cr2O3(i) < 2))
                    T.group4(i) = {'Spinel'};
                    T.species(i) = {'magnetite'};
                    Fe3 = 2 - (Al_n/6 + Cr_n/6 + Ti_n/6 + V_n/6);
                    if(Fe3<=0), Fe3 = 0; end
                    Fe2 = Fe_n/6 - Fe3;
                    if(Fe2<=0), Fe2 = 0; end
                    [sorted_formula1, total1]  = Formula_Output({Fe2, Mg_n/6, Mn_n/6, Ni_n/6}, {'Fe2+', 'Mg', 'Mn', 'Ni'});
                    [sorted_formula2, total2]  = Formula_Output({Fe3, Al_n/6, Cr_n/6, Ti_n/6, V_n/6}, {'Fe3+', 'Al', 'Cr', 'Ti','V'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O4')};
                end

            end

            Renorm_factor = 6;
            Fe_n4 = T.Fe_Multi(i) / Renorm_factor;
            Mg_n4 = T.Mg_Multi(i) / Renorm_factor;
            Mn_n4 = T.Mn_Multi(i) / Renorm_factor;
            Ni_n4 = T.Ni_Multi(i) / Renorm_factor;
            Zn_n4 = T.Zn_Multi(i) / Renorm_factor;
            Al_n4 = T.Al_Multi(i) / Renorm_factor;
            Cr_n4 = T.Cr_Multi(i) / Renorm_factor;
            Ti_n4 = T.Ti_Multi(i) / Renorm_factor;
            Si_n4 = T.Si_Multi(i) / Renorm_factor;
            V_n4 = T.V_Multi(i) / Renorm_factor;
            Co_n4 = T.Co_Multi(i) / Renorm_factor;
            CD_n4 = Fe_n4 + Mg_n4 + Mn_n4 + Ni_n4 + Zn_n4 + Al_n4 + Cr_n4 + Ti_n4 + V_n4 + Si_n4;
            [T, Fe2, Fe3, CD_n4_new] = Fe2O3_Calculation(T, i, CD_n4, Fe_n4, Renorm_factor, 1, 3.000);

            if(CD_n4_new > 2.95-2.95*t && CD_n4_new < 3.5+3.5*t)
                T.group2(i) = {'Oxide'};
                T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide'};
                T.group4(i) = {'Spinel'};

                if(Ti_n4 >= 0.7 && Fe_n4 >= 1.9)
                    T.species(i) = {'ulvöspinel'};
                    Fe3 = 2 - (Al_n/6 + Cr_n/6 + Ti_n/6 + V_n/6);
                    if(Fe3<=0), Fe3 = 0; end
                    Fe2 = Fe_n/6 - Fe3;
                    if(Fe2<=0), Fe2 = 0; end
                    [sorted_formula1, total1]  = Formula_Output({Fe2, Mg_n/6, Mn_n/6, Ni_n/6, Co_n/6}, {'Fe2+', 'Mg', 'Mn', 'Ni', 'Co'});
                    [sorted_formula2, total2]  = Formula_Output({Fe3, Al_n/6, Cr_n/6, Ti_n/6}, {'Fe3+', 'Al', 'Cr', 'Ti'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O4')};
                end

                if (Ti_n4 > 0 && Ti_n4 < 0.2)
                    T.species(i) = {'magnetite'};
                    Fe3 = 2 - (Al_n/6 + Cr_n/6 + Ti_n/6 + V_n/6);
                    if(Fe3<=0), Fe3 = 0; end
                    Fe2 = Fe_n/6 - Fe3;
                    if(Fe2<=0), Fe2 = 0; end
                    [sorted_formula1, total1]  = Formula_Output({Fe2, Mg_n/6, Mn_n/6, Ni_n/6}, {'Fe2+', 'Mg', 'Mn', 'Ni'});
                    [sorted_formula2, total2]  = Formula_Output({Fe3, Al_n/6, Cr_n/6, Ti_n/6, V_n/6}, {'Fe3+', 'Al', 'Cr', 'Ti','V'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O4')};
                end

                if (Ti_n4 >= 0.2 && Ti_n4 < 0.7)
                    T.species(i) = {'Ti-magnetite'};
                    Fe3 = 2 - (Al_n/6 + Cr_n/6 + Ti_n/6 + V_n/6);
                    if(Fe3<=0), Fe3 = 0; end
                    Fe2 = Fe_n/6 - Fe3;
                    if(Fe2<=0), Fe2 = 0; end
                    [sorted_formula1, total1]  = Formula_Output({Fe2, Mg_n/6, Mn_n/6, Ni_n/6}, {'Fe2+', 'Mg', 'Mn', 'Ni'});
                    [sorted_formula2, total2]  = Formula_Output({Fe3, Al_n/6, Cr_n/6, Ti_n/6, V_n/6}, {'Fe3+', 'Al', 'Cr', 'Ti','V'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O4')};
                end

                if((Ti_n4 + Fe3 + Cr_n4 + Al_n4 + V_n4) / (Mg_n4 + Fe2 + Mn_n4 + Ni_n4 + Zn_n4) > 1.95-1.95*t && (Ti_n4 + Fe3 + Cr_n4 + Al_n4 + V_n4) / (Mg_n4 + Fe2 + Mn_n4 + Ni_n4 + Zn_n4) < 2.05+2.05*t)

                    if(Al_n4 / (Al_n4 + Fe3 + Cr_n4 + Ti_n4) > 0.5)
                        if(Mg_n4 / (Mg_n4 + Fe2 + Mn_n4 + Ni_n4 + Zn_n4) > 0.5)
                            T.species(i) = {'spinel'};
                            [sorted_formula1, total1]  = Formula_Output({Fe2, Mg_n4, Mn_n4, Ni_n4, Co_n4}, {'Fe2+', 'Mg', 'Mn', 'Ni', 'Co'});
                            [sorted_formula2, total2]  = Formula_Output({Fe3, Ti_n4, Al_n4, Cr_n4}, {'Fe3+', 'Ti', 'Al', 'Cr'});
                            T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O4')};
                        end
                    end

                    if(Mg_n4 / (Mg_n4 + Fe2) > 0.5)
                        if(Cr_n4 / (Cr_n4 + Fe3 + Al_n4 + Ti_n4) > 0.5)
                            T.species(i) = {'magnesiochromite'};
                            [sorted_formula1, total1]  = Formula_Output({Fe2, Mg_n4, Mn_n4, Ni_n4}, {'Fe2+', 'Mg', 'Mn', 'Ni'});
                            [sorted_formula2, total2]  = Formula_Output({Fe3, Ti_n4, Al_n4, Cr_n4, V_n4}, {'Fe3+', 'Ti', 'Al', 'Cr', 'V'});     
                            T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, 'O4')};
                        end
                    end

                    if(Fe_n4 / (Mg_n4 + Fe2) > 0.5)
                        if(Cr_n4 / (Cr_n4 + Fe3 + Al_n4 + Ti_n4) > 0.5)
                            T.species(i) = {'chromite'};
                            [sorted_formula1, total1]  = Formula_Output({Fe2, Mg_n4, Mn_n4, Ni_n4, Zn_n4}, {'Fe2+', 'Mg', 'Mn', 'Ni', 'Zn'});
                            [sorted_formula2, total2]  = Formula_Output({Fe3, Ti_n4, Al_n4, Cr_n4, V_n4,Si_n4}, {'Fe3+', 'Ti', 'Al', 'Cr', 'V','Si'});
                            T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2,'O4')};
                        end
                    end

                    if(Al_n4 / (Cr_n4 + Fe3 + Al_n4 + Ti_n4 + Zn_n4) > 0.5)
                        if(Zn_n4 / (Mg_n4 + Fe2 + Zn_n4) > 0.5)
                            T.species(i) = {'gahnite'};
                            Fe3 = 2-(Al_n4+Ti_n4+Cr_n4+V_n4);
                            if(Fe3<=0), Fe3 = 0; end
                            Fe2 = Fe_n4 - Fe3;
                            if(Fe2<=0), Fe2 = 0; end
                            [sorted_formula1, total1]  = Formula_Output({Fe2, Mg_n4, Mn_n4, Ni_n4, Zn_n4}, {'Fe2+', 'Mg', 'Mn', 'Ni', 'Zn'});
                            [sorted_formula2, total2]  = Formula_Output({Fe3, Ti_n4, Al_n4, Cr_n4, V_n4}, {'Fe3+', 'Ti', 'Al', 'Cr', 'V'});
                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O4')};
                        end

                        if(Fe_n4 / (Mg_n4 + Fe2 + Zn_n4) > 0.5)
                            T.species(i) = {'hercynite'};
                            Fe3 = 2-(Al_n4+Ti_n4+Cr_n4+V_n4);
                            if(Fe3<=0), Fe3 = 0; end
                            Fe2 = Fe_n4 - Fe3;
                            if(Fe2<=0), Fe2 = 0; end
                            [sorted_formula1, total1]  = Formula_Output({Fe2, Mg_n4, Mn_n4, Ni_n4, Zn_n4}, {'Fe2+', 'Mg', 'Mn', 'Ni', 'Zn'});
                            [sorted_formula2, total2]  = Formula_Output({Fe3, Ti_n4, Al_n4, Cr_n4, V_n4}, {'Fe3+', 'Ti', 'Al', 'Cr', 'V'});
                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O4')};
                        end
                    end

                end

            end

        end %end of Spinels

        %% 5.1.1.2 Ilmenite
        if((T.FeO(i) + T.MgO(i) + T.MnO(i) + T.ZnO(i) + T.NiO(i)) > 30 && T.TiO2(i) > 40)
            if((Mg_n + Fe_n + Mn_n + Zn_n + Ni_n)/Ti_n > 1-2*1*t && (Mg_n + Fe_n + Mn_n + Zn_n + Ni_n)/Ti_n < 1+2*1*t)
                if((Al_n + Cr_n + V_n)/Ti_n < 0.2)
                    T.group2(i) = {'Oxide'};
                    T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide'};
                    T.group4(i) = {'Ilmenite'};
                    Renorm_factor = 8;
                    Fe_n3 = T.Fe_Multi(i) / Renorm_factor;
                    Mg_n3 = T.Mg_Multi(i) / Renorm_factor;
                    Mn_n3 = T.Mn_Multi(i) / Renorm_factor;
                    Ni_n3 = T.Ni_Multi(i) / Renorm_factor;
                    Cu_n3 = T.Cu_Multi(i) / Renorm_factor;
                    Zn_n3 = T.Zn_Multi(i) / Renorm_factor;
                    Al_n3 = T.Al_Multi(i) / Renorm_factor;
                    Cr_n3 = T.Cr_Multi(i) / Renorm_factor;
                    Ti_n3 = T.Ti_Multi(i) / Renorm_factor;
                    Si_n3 = T.Si_Multi(i) / Renorm_factor;
                    Ca_n3 = T.Ca_Multi(i) / Renorm_factor;
                    Fe3 = 1 - (Ti_n3 + Al_n3 + Cr_n3);
                    if(Fe3 < 0), Fe2 = Fe_n3; else, Fe2 = Fe_n3 - Fe3; end

                    if(Fe2 / (Fe2 + Mg_n3 + Mn_n3 + Zn_n3 + Ni_n3) > 0.5)
                        T.species(i) = {'ilmenite'};  
                        [sorted_formula1, total1]  = Formula_Output({Fe_n3, Mg_n3, Mn_n3}, {'Fe', 'Mg', 'Mn'});  
                        [sorted_formula2, total2]  = Formula_Output({Ti_n3, Al_n3, Cr_n3}, {'Ti', 'Al','Cr'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O3')};
                    end

                    if(Zn_n3 / (Zn_n3 + Fe2 + Mg_n3 + Mn_n3) > 0.5)
                        T.species(i) = {'ecandrewsite'};
                        [sorted_formula1, total1]  = Formula_Output({Zn_n3, Fe_n3, Cu_n3}, {'Zn','Fe','Cu'});
                        [sorted_formula2, total2]  = Formula_Output({Ti_n3, Si_n3}, {'Ti','Si'});
                        T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2,'O3')};
                    end

                    if(Mg_n3 / (Mg_n3 + Fe2 + Mn_n3 + Zn_n3 + Ni_n3) > 0.5)
                        T.species(i) = {'geikielite'};
                        [sorted_formula1, total1]  = Formula_Output({Fe2, Mg_n3, Mn_n3}, {'Fe2+', 'Mg', 'Mn'});
                        [sorted_formula2, total2]  = Formula_Output({Fe3, Ti_n3, Al_n3,Cr_n3,Si_n3}, {'Fe3+', 'Ti', 'Al','Cr','Si'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, 'O3')};
                    end

                    if(Mn_n3 / (Mn_n3 + Mg_n3 + Fe2 + Zn_n3 + Ni_n3) > 0.5)
                        T.species(i) = {'pyrophanite'};
                        [sorted_formula1, total1]  = Formula_Output({Fe2, Mg_n3, Mn_n3, Ca_n3}, {'Fe2+', 'Mg', 'Mn', 'Ca'});
                        [sorted_formula2, total2]  = Formula_Output({Fe3, Ti_n3, Al_n3, Cr_n3}, {'Fe3+', 'Ti', 'Al', 'Cr'}); 
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, 'O3')};
                    end

                end
            end
        end %end of Ilmenties

        %% 5.1.1.3 Fe-Oxides
        if(T.FeO(i) > 95 && Fe_n / 24 > 0.95 && Fe_n / 24 < 1.05 && T.CaO(i) < 0.000001)
            if(Mn_n > Mg_n)
                if(round(Al_n + Cr_n + Ti_n + V_n,2) <= 0.01)
                    Renorm_factor = 24;
                    Fe_n1 = T.Fe_Multi(i) / Renorm_factor;
                    Mn_n1 = T.Mn_Multi(i) / Renorm_factor;
                    Mg_n1 = T.Mg_Multi(i) / Renorm_factor;
                    T.group2(i) = {'Oxide or Hydroxide'};
                    T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
                    T.group4(i) = {'Fe-Oxide or Fe-Hydroxide'};
                    T.species(i) = {'wüstite'};
                    [sorted_formula1, total1]  = Formula_Output({Fe_n1, Mn_n1, Mg_n1}, {'Fe', 'Mn', 'Mg'});
                    T.formula(i) = {strcat(sorted_formula1, total1, 'O')};
                end
            end
        end

        %% 5.1.1.4 Mn-Oxides or Hydroxides
        if(T.MnO(i) > 40)

            if(Mn_n / 24 > 3-3*t && Mn_n / 24 < 3+3*t)
                Renorm_factor = 24;
                Cl_n1 = T.Cl_Multi(i) / Renorm_factor;
                Fe_n1 = T.Fe_Multi(i) / Renorm_factor;
                Mn_n1 = T.Mn_Multi(i) / Renorm_factor;
                Ca_n1 = T.Ca_Multi(i) / Renorm_factor;
                Mg_n1 = T.Mg_Multi(i) / Renorm_factor;
                F_n1 =T.F_Multi(i) / Renorm_factor;

                if(Cl_n1 == 0)
                    if(T.CaO(i) > 0 && T.CaO(i) < 0.4)
                        T.group2(i) = {'Oxide or Hydroxide'};
                        T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
                        T.group4(i) = {'Mn-Oxide or Mn-Hydroxide'};
                        T.species(i) = {'manganosite'};
                        [sorted_formula2, total2] = Formula_Output({Mn_n/24, Fe_n/24, Mg_n/24, Zn_n/24, Ca_n/24}, {'Mn','Fe','Mg','Zn','Ca'});
                        T.formula(i) = {strcat(sorted_formula2, total2, 'O')};
                    end
                end

                if(Cl_n1 > 0)
                    T.group2(i) = {'Hydroxide'};
                    T.group3(i) = {'Fe-Al-Ti-Mn-Mg Hydroxide'};
                    T.group4(i) = {'Mn-Hydroxide'};
                    T.species(i) = {'pyrochroite'};
                    [sorted_formula, total]  = Formula_Output({Mn_n1, Fe_n1, Mg_n1, Ca_n1}, {'Mn2+', 'Fe','Mg','Ca'});
                    [sorted_formula2,total2] = Formula_Output({1 - (F_n1+Cl_n1), F_n1, Cl_n1},{'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula,total1,sorted_formula2,total2)};
                end
            end

            if((Mn_n + Fe_n + Al_n)/24 > (2/3)-(2/3)*t && (Mn_n + Fe_n + Al_n)/24 < (2/3)+(2/3)*t && Fe_n > Al_n && Mn_n > Al_n)
                if((Ca_n+Na_n)/(Mn_n+Fe_n+Al_n) < 0.001)
                    Renorm_factor = 8;
                    Mn_n3 = T.Mn_Multi(i) / Renorm_factor;
                    Fe_n3 = T.Fe_Multi(i) / Renorm_factor;
                    Al_n3 = T.Al_Multi(i) / Renorm_factor;
                    Cr_n3 = T.Cr_Multi(i) / Renorm_factor;
                    Ti_n3 = T.Ti_Multi(i) / Renorm_factor;
                    V_n3 = T.V_Multi(i) / Renorm_factor;
                    Mg_n3 = T.Mg_Multi(i) / Renorm_factor;
                    T.group2(i) = {'Oxide or Hydroxide'};
                    T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
                    T.group4(i) = {'Mn-Oxide or Mn-Hydroxide'};
                    T.species(i) = {'pyrolusite/bixbyite'};
                    [sorted_formula, total]  = Formula_Output({Mn_n3, Fe_n3,Al_n3,Cr_n3,Ti_n3,V_n3,Mg_n3}, {'Mn3+', 'Fe3+','Al','Cr','Ti','V','Mg'});
                    T.formula(i) = {strcat(sorted_formula,total,'O3')};
                end    
            end

            if(Mn_n / 24 > 0.5-0.5*t && Mn_n / 24 < 0.5+0.5*t && Al_n == 0 && Mg_n >= 0 && (Ca_n+Na_n)/(Mn_n+Fe_n+Al_n) < 0.001)
                T.group2(i) = {'Oxide or Hydroxide'};
                T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
                T.group4(i) = {'Mn-Oxide or Mn-Hydroxide'};
                T.species(i) = {'pyrolusite/bixbyite'};
                [sorted_formula, total]  = Formula_Output({Mn_n2, Fe_n2,Al_n2,Cr_n2,Ti_n2,V_n2,Mg_n2}, {'Mn', 'Fe','Al','Cr','Ti','V','Mg'});
                T.formula(i) = {strcat(sorted_formula,total,'O2')};
            end

            if((Mn_n + Fe_n)/24 > (18/24)-(18/24)*t && (Mn_n + Fe_n)/24 < (18/24)+(18/24)*t && Mn_n > Fe_n) 
                Renorm_factor = 6;
                Mn_n4 = T.Mn_Multi(i) / Renorm_factor;
                Fe_n4 = T.Fe_Multi(i) / Renorm_factor;
                Al_n4 = T.Al_Multi(i) / Renorm_factor;
                Na_n4 = T.Na_Multi(i) / Renorm_factor;
                K_n4 = T.K_Multi(i) / Renorm_factor;
                Ca_n4 = T.Ca_Multi(i) / Renorm_factor;
                Ba_n4 = T.Ba_Multi(i) / Renorm_factor;

                if(Na_n4 + K_n4 < 0.02 && Na_n+Ca_n < 0.01)
                    T.group2(i) = {'Oxide or Hydroxide'};
                    T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
                    T.group4(i) = {'Mn-Oxide or Mn-Hydroxide'};
                    T.species(i) = {'hausmannite'};
                    Mn3 = 2 - (Fe_n4 + Al_n4);
                    if(Mn3 <=0), Mn3 = 0; end
                    Mn2 = Mn_n4 - Mn3;
                    if(Mn2 <=0), Mn2 = 0; end
                    [sorted_formula1, total1] = Formula_Output({Mn2, Ca_n4, Ba_n4},{'Mn2+','Ca','Ba'});
                    [sorted_formula2, total2] = Formula_Output({Mn3, Fe_n4, Al_n4},{'Mn3+','Fe3+','Al'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O4')};
                end

                if(Na_n4 + K_n4 + Ca_n4 > 0.02 && Na_n4 + K_n4 + Ca_n4 < 1 && Na_n4 > Ca_n4 && Ca_n4 > K_n4)
                    if((Na_n+Ca_n+K_n)/(Mn_n+Fe_n) > 0.6/2-(0.6/2)*t && (Na_n+Ca_n+K_n)/(Mn_n+Fe_n) < (0.6/2)+(0.6/2)*t)
                        if(Mn_n > Fe_n && Mn_n > Mg_n)
                            T.group2(i) = {'Oxide or Hydroxide'};
                            T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
                            T.group4(i) = {'Mn-Oxide or Mn-hydroxide'};
                            T.species(i) = {'birnessite'};
                            Fe3 = Fe_n2*2;
                            Mn3 = 2 - (Fe3 + Cr_n2*2 + Ti_n2*2 + Mg_n2*2);
                            if(Mn3<0), Mn3=0; end
                            Mn4 = Mn_n2*2 - Mn3;
                            if(Mn4<0), Mn4=0; end
                            [sorted_formula1, total1] = Formula_Output({Na_n2*2, Ca_n2*2, K_n2*2}, {'Na','Ca','K'});
                            [sorted_formula2, total2] = Formula_Output({Mn4,Mn3,Cr_n2*2,Ti_n2*2,Mg_n2*2,Fe3}, {'Mn4+','Mn3+','Cr','Ti','Mg','Fe3+'});
                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O4·1.5H2O')};
                        end
                    end
                end

            end
            
        end %end of Mn-Oxides/Hydroxides

        %% 5.1.1.5 Mg-Oxide
        if(T.MgO(i) > 90 && T.MgO(i) < 100)
            if((Fe_n + Mg_n)/24 > 1 - 1*(0.5*t) && (Fe_n + Mg_n)/24 < 1 + 1*(0.5*t))
                if(Mg_n/(Mg_n+Fe_n)>0.75)

                    if(Fe_n > (Ca_n + Mn_n))
                        T.group2(i) = {'Oxide'};
                        T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide'};
                        T.group4(i) = {'Mg-Oxide'};
                        T.species(i) = {'periclase'};
                        [sorted_formula,total] = Formula_Output({Mg_n/24,Fe_n/24,Mn_n/24},{'Mg','Fe','Mn'});
                        T.formula(i) = {strcat(sorted_formula,total, 'O')};
                    else
                        T.group2(i) = {'Orthorhombic Carbonate'};
                        T.group3(i) = {'Fe-Mg-Sr-Pb-Zn Carbonate'};
                        T.group4(i) = {''};
                        T.species(i) = {'magnesite'};
                        [sorted_formula1, total1] = Formula_Output({Fe_n2/2,Mn_n2/2,Mg_n2/2,Ca_n2/2,Sr_n2/2,Ba_n2/2}, {'Fe','Mn','Mg','Ca','Sr','Ba'});
                        T.formula(i) = {strcat(sorted_formula1,total1, '(CO3)')};
                    end

                    Renorm_factor = 24;
                    Fe_n1 = T.Fe_Multi(i) / Renorm_factor;
                    Mg_n1 = T.Mg_Multi(i) / Renorm_factor;
                    Mn_n1 = T.Mn_Multi(i) / Renorm_factor;

                    if(Ca_n == 0)
                        T.group2(i) = {'Oxide or Hydroxide'};
                        T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
                        T.group4(i) = {'Mg-Oxide or Mg-Hydroxide'};
                        T.species(i) = {'periclase or brucite'};
                        [sorted_formula,total] = Formula_Output({Mg_n1,Fe_n1,Mn_n1},{'Mg','Fe','Mn'});
                        T.formula(i) = {strcat(sorted_formula,total, '(OH)2')};
                    end

                    if(Ca_n > 0 || Sr_n > 0 || Ba_n > 0)
                        T.group2(i) = {'Orthorhombic Carbonate'};
                        T.group3(i) = {'Fe-Mg-Sr-Pb-Zn Carbonate'};
                        T.group4(i) = {''};
                        T.species(i) = {'magnesite'};
                        [sorted_formula1, total1] = Formula_Output({Fe_n2/2,Mn_n2/2,Mg_n2/2,Ca_n2/2,Sr_n2/2,Ba_n2/2}, {'Fe','Mn','Mg','Ca','Sr','Ba'}); 
                        T.formula(i) = {strcat(sorted_formula1,total1, '(CO3)')};
                    end
                    
                end
            end
        end %end of Mg-Oxide

        %% 5.1.1.6 Ti-Oxide & Cr-Oxide
        if(strcmp(T.group4(i),'Ti-Oxide'))
            T.species(i) = {'rutile'};
            [sorted_formula, total]  = Formula_Output({Ti_n2, Fe_n2, Al_n2, Cr_n2}, {'Ti','Fe3+','Al','Cr'});
            T.formula(i) = {strcat(sorted_formula, total, 'O2')};
        end

        if(strcmp(T.group4(i), 'Cr-Oxide'))
            T.species(i) = {'eskolaite'};
            [sorted_formula1, total1] = Formula_Output({T.Cr_Multi(i)/8, T.V_Multi(i)/8, T.Fe_Multi(i)/8, T.Mn_Multi(i)/8, T.Mg_Multi(i)/8},{'Cr','V','Fe','Mn','Mg'});
            T.formula(i) = {strcat(sorted_formula1, total1, 'O3')};
        end

        %% 5.1.1.7 Al-Oxides
        if(strcmp(T.group4(i), 'Al-Oxide or Hydroxide'))
            T.species(i) = {'Al Oxide or Hydroxide'};
%                 'Al-Oxide or Al-Hydroxide (corundum, gibbsite, diaspore, or others)']};
            T.check(i) = {'no check'};
            Renorm_factor = 24 / 1.5;
            Fe_n15 = T.Fe_Multi(i) / Renorm_factor;
            Al_n15 = T.Al_Multi(i) / Renorm_factor;
            Mg_n15 = T.Mg_Multi(i) / Renorm_factor;
            Cr_n15 = T.Cr_Multi(i) / Renorm_factor;
            Ti_n15 = T.Ti_Multi(i) / Renorm_factor;
            Cl_n15 = T.Cl_Multi(i) / Renorm_factor;
            F_n15 = T.F_Multi(i) / Renorm_factor;
            OH = 3 - (F_n15 + Cl_n15);
            if(OH<=0), OH=0; end
            [sorted_formula1, total1]  = Formula_Output({Al_n15, Fe_n15, Mg_n15, Cr_n15, Ti_n15}, {'Al', 'Fe3+','Mg','Cr','Ti'});
            [sorted_formula2, total2]  = Formula_Output({OH, F_n15, Cl_n15}, {'OH', 'F','Cl'});
            T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2)};

            if(Cr_n > 0)
                T.group2(i) = {'Oxide'};
                T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide'};
                T.group4(i) = {'Al-Oxide'};
                T.species(i) = {'corundum'};
                Renorm_factor = 8;
                Fe_n3 = T.Fe_Multi(i) / Renorm_factor;
                Al_n3 = T.Al_Multi(i) / Renorm_factor;
                Cr_n3 = T.Cr_Multi(i) / Renorm_factor;
                Ti_n3 = T.Ti_Multi(i) / Renorm_factor;
                Mg_n3 = T.Mg_Multi(i) / Renorm_factor;
                [sorted_formula, total]  = Formula_Output({Al_n3, Fe_n3, Ti_n3, Cr_n3,Mg_n3}, {'Al','Fe3+','Ti','Cr','Mg'});
                T.formula(i) = {strcat(sorted_formula, total, 'O3')};
            end

        end %end of Al-Oxide/Hydroxide

        %% 5.1.2 Other oxides
        if(strcmp(T.group3(i), 'Other Oxide'))
            T.group4(i) = {'Zr Oxide'};
            T.species(i) = {'baddeleyite'};
            Renorm_factor = 12;
            La_n2 = T.La_Multi(i) / Renorm_factor;
            Ce_n2 = T.Ce_Multi(i) / Renorm_factor;
            Nd_n2 = T.Nd_Multi(i) / Renorm_factor;
            Hf_n2 = T.Hf_Multi(i) / Renorm_factor;
            Zr_n2 = T.Zr_Multi(i) / Renorm_factor;
            [sorted_formula, total]  = Formula_Output({Zr_n2,Hf_n2, Nd_n2,La_n2,Ce_n2,}, {'Zr','Hf','Nd','La','Ce'});
            T.formula(i) = {strcat(sorted_formula,total,'O2')};
        end

        %% 5.2 Hydroxides

        %% 5.2.1 Fe-Hydroxides
        if(T.FeO(i) > 80)  
            if(ABCDT_n/24 > (1/1.5) - (1/1.5)*t && ABCDT_n/24 < (1/1.5) + (1/1.5)*t)
                if(Fe_n / ABCDT_n > 1-1*t)
                    T.group2(i) = {'Hydroxide'}; 
                    T.group3(i) = {'Fe-Al-Ti-Mn-Mg Hydroxide'};
                    T.group4(i) = {'Fe-Hydroxide'};
                    T.species(i) = {''};

                    if (S_n > 0.1)
                        T.group2(i) = {'K-Na-Mg-Fe-Ca Sulfate'};
                        T.group3(i) = {'Fe-Sulfate'};
                        T.group4(i) = {''};
                        T.species(i) = {'schwertmannite'};
                        renorm = 33.6/24;
                        OH = 9.6 - (F_n*renorm + Cl_n*renorm);
                        if(OH<=0), OH = 0; end 
                        [sorted_formula1, total1] = Formula_Output({Fe_n*renorm, Al_n*renorm, Mg_n*renorm, Mn_n*renorm, Ca_n*renorm, Na_n*renorm},{'Fe','Al','Mg','Mn','Ca','Na'});
                        [sorted_formula2, total2] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                        [sorted_formula3, total3] = Formula_Output({S_n*renorm/3.2, P_n*renorm/3.2},{'S', 'P'});
                        T.formula(i) = {strcat(sorted_formula1,total1, 'O16', sorted_formula2, total2, '(',sorted_formula3, total3, 'O4)3.2·10H2O')};
                    end

                end
            end
        end %end of Fe-Hydroxides group

        if(T.FeO(i) > 75)
            if(Fe_n/24 > 0.6 && Fe_n/24 < 1+0.5*t && (Fe_n/24 + Mg_n/24 + Mn_n/24) > 1- 1*t && (Fe_n/24 + Mg_n/24 + Mn_n/24) < 1+ 1*t)
                if(Mg_n > Mn_n && (Mn_n+Mg_n)/24 <= 0.45 + 0.45*t && (Mn_n+Mg_n)/24 > 0 && Ca_n/24 < 0.01 && T.TiO2(i) < 1 && (Mg_n+Mn_n)/Fe_n < 0.49)
                    T.group2(i) = {'Oxide or Hydroxide'};
                    T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
                    T.group4(i) = {'Fe-Oxide or Fe-Hydroxide'};
                    T.species(i) = {'amakinite or wüstite'};
                    OH = 2 - (Cl_n/24 + F_n/24);
                    [sorted_formula1, total1]  = Formula_Output({Fe_n/24, Mg_n/24, Mn_n/24, Al_n/24}, {'Fe', 'Mg', 'Mn', 'Al'});
                    [sorted_formula2, total2] = Formula_Output({OH, Cl_n/24, F_n/24}, {'OH','Cl','F'});
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2)};
                end
            end
        end

        %% 5.2.2 Al-Hydroxides
        if(T.Al2O3(i) > 80)
            if(ABCDT_n/24 > (1/1.5) - (1/1.5)*t && ABCDT_n/24 < (1/1.5) + (1/1.5)*t)
                if(Al_n/ABCDT_n > 0.9)
                    T.group2(i) = {'Oxide or Hydroxide'};
                    T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
                    T.group4(i) = {'Al-Oxide or Al-Hydroxide'};
                    Renorm_factor = 24 / 1.5;
                    Fe_n15 = T.Fe_Multi(i) / Renorm_factor;
                    Al_n15 = T.Al_Multi(i) / Renorm_factor;
                    Mg_n15 = T.Mg_Multi(i) / Renorm_factor;
                    Cr_n15 = T.Cr_Multi(i) / Renorm_factor;
                    Ti_n15 = T.Ti_Multi(i) / Renorm_factor;
                    Cl_n15 = T.Cl_Multi(i) / Renorm_factor;
                    F_n15 = T.F_Multi(i) / Renorm_factor;
                    T.species(i) = {'Al Oxide or Hydroxide'};
                    OH = 3 - (F_n15 + Cl_n15);
                    if(OH<=0), OH=0; end
                    [sorted_formula1, total1]  = Formula_Output({Al_n15, Fe_n15, Mg_n15, Cr_n15, Ti_n15}, {'Al', 'Fe3+','Mg','Cr','Ti'});
                    [sorted_formula2, total2]  = Formula_Output({OH, F_n15, Cl_n15}, {'OH', 'F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1, total, sorted_formula2, total2)};
                    T.check(i) = {'no check'};
                end
            end
        end

        %% 5.2.3 Mn-Hydroxide
        if(Ag_n > 0 && Ca_n > 0)
            Renorm_factor = 24 / 7;
            Zn_n7 = T.Zn_Multi(i) / Renorm_factor;
            Ag_n7 = T.Ag_Multi(i) / Renorm_factor;
            Ca_n7 = T.Ca_Multi(i) / Renorm_factor;
            Na_n7 = T.Na_Multi(i) / Renorm_factor;
            Al_n7 = T.Al_Multi(i) / Renorm_factor;
            Cr_n7 = T.Cr_Multi(i) / Renorm_factor;
            Fe_n7 = T.Fe_Multi(i) / Renorm_factor;
            K_n7 =  T.K_Multi(i) / Renorm_factor;
            Mg_n7 = T.Mg_Multi(i) / Renorm_factor;
            Ba_n7 = T.Ba_Multi(i) / Renorm_factor;
            Cu_n7 = T.Cu_Multi(i) / Renorm_factor;
            Pb_n7 = T.Pb_Multi(i) / Renorm_factor;
            Mn_n7 = T.Mn_Multi(i) / Renorm_factor;
            Mn2 = 1 - (Ag_n7 + Ca_n7 + K_n7 + Ba_n7 + Cu_n7 + Pb_n7 + Mg_n7 + Zn_n7);
            if(Mn2 < 0), Mn2 = 0; end
            Mn4 = Mn_n7 - Mn2;
            if(Mn4 < 0), Mn4 = 0; end
            Fe2 = 1 - (Mn2+Mg_n7+Ca_n7+Na_n7);
            if(Fe2<0), Fe2=0; end
            Fe3 = Fe_n7-Fe2;
            if(Fe3<0), Fe3=0; end
            T.group2(i) = {'Hydroxide'};
            T.group3(i) = {'Fe-Al-Ti-Mn-Mg Hydroxide'};
            T.group4(i) = {'Mn-Hydroxide'};
            T.species(i) = {'aurorite'};
            [sorted_formula1, total1] = Formula_Output({Mn2,Fe2,Ag_n7,Ca_n7,Na_n7}, {'Mn2+','Fe2+','Ag','Ca','Na'});
            [sorted_formula2, total2] = Formula_Output({Mn4,Fe3,Al_n7,Cr_n7}, {'Mn4+','Fe3+','Al','Cr'});
            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O7·3H2O')};
        end

        if(Mn_n/24 > (1/1.5)-(1/1.5)*t && Mn_n/24 < (1/1.5)+(1/1.5)*t && Ca_n+Na_n < 0.01)
            T.group2(i) = {'Oxide or Hydroxide'};
            T.group3(i) = {'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'};
            T.group4(i) = {'Mn-Oxide or Mn-Hydroxide'};
            T.species(i) = {'manganosite or manganite'};
            [sorted_formula1, total1] = Formula_Output({Mn_n*(1.5/24), Fe_n*(1.5/24), Al_n*(1.5/24)}, {'Mn3+','Fe3+','Al'});
            [sorted_formula2, total2] = Formula_Output({Mn_n/24, Fe_n/24, Mg_n/24, Zn_n/24, Ca_n/24}, {'Mn','Fe','Mg','Zn','Ca'});
            T.formula(i) = {strcat(sorted_formula2, total2, 'O or', sorted_formula1, total1, 'O(OH)')};
        end

        if(round(Mn_n/24,1) >= (1/1.5)-(1/1.5)*t && Mn_n/24 < (1/1.5)+(1/1.5)*t && Ti_n+Ca_n+Na_n > 0.1 && T.Na2O(i) < 5)
            T.group2(i) = {'Hydroxide'};
            T.group3(i) = {'Fe-Al-Ti-Mn-Mg Hydroxide'};
            T.group4(i) = {'Mn-Hydroxide'};
            T.species(i) = {'vernadite'};
            T.formula(i) = {'(Mn4+,Fe3+,Ca,Na)(O,OH)2·nH2O'};
        end

        if((Mn_n+Fe_n)/24 > (6/11)-(6/11)*t && (Mn_n+Fe_n)/24 < (6/11)+(6/11)*t && Mn_n > Fe_n && Ca_n+Na_n+K_n < 0.1)
            Renorm_factor = 24/11;
            Al_n11 = T.Al_Multi(i) / Renorm_factor;
            Mn_n11 = T.Mn_Multi(i) / Renorm_factor;
            Fe_n11 = T.Fe_Multi(i) / Renorm_factor;
            F_n11 = T.F_Multi(i) / Renorm_factor;
            Cl_n11 = T.Cl_Multi(i) / Renorm_factor;
            T.group2(i) = {'Hydroxide'};
            T.group3(i) = {'Fe-Al-Ti-Mn-Mg Hydroxide'};
            T.group4(i) = {'Mn-Hydroxide'};
            T.species(i) = {'janggunite'};
            OH = 6 - (F_n11 + Cl_n11);
            if(OH<=0), OH=0; end
            [sorted_formula1, total1]  = Formula_Output({Mn_n11*(11/7), Fe_n11, Al_n11}, {'Mn', 'Fe3+', 'Al'});
            [sorted_formula2, total2]  = Formula_Output({OH, F_n11, Cl_n11}, {'OH', 'F', 'Cl'});
            T.formula(i) = {strcat(sorted_formula1,total1, 'O8',sorted_formula2,total2)};
        end

        if((Mn_n+Ca_n)/24 > (5/9)-(5/9)*t && (Mn_n+Ca_n)/24 < (5/9)+(5/9)*t && Ca_n > 0)
            if(Ca_n+Mn_n+Mg_n > Fe_n && Mn_n > Ca_n && Ca_n > Mg_n && Na_n < 0.01 + 3*(0.01)*(t))
                T.group2(i) = {'Hydroxide'};
                T.group3(i) = {'Fe-Al-Ti-Mn-Mg Hydroxide'};
                T.group4(i) = {'Mn-Hydroxide'};
                T.species(i) = {'takanelite'};
                [sorted_formula, total] = Formula_Output({Mn_n2*3,Ca_n2*3}, {'Mn2+','Ca'});
                T.formula(i) = {'(' + sorted_formula + ')_{Σ=' + total + '}(Mn4+)1-xO2·0.7H2O'};
            end
        end

        if(Ba_n/Mn_n > (1/8)-2*(1/8)*t && Ba_n/Mn_n < (2/5)+2*(2/5)*t)  
            T.group2(i) = {'Oxide'};
            T.group3(i) = {'Other Oxide'};
            T.group4(i) = {'Ba-Mn Oxide'};
            T.species(i) = {'hollandite or romanèchite'};
            renorm = 16/24;
            Mn2 = 1 - (Ba_n*renorm + Ca_n*renorm + Sr_n*renorm + Pb_n*renorm + Na_n*renorm + K_n*renorm);
            if(Mn2<=0), Mn2=0; end
            Mn3 = 2 - (Al_n*renorm + Fe_n*renorm);
            if(Mn3<=0), Mn3=0; end
            Mn4 = Mn_n*renorm - (Mn2 + Mn3);
            if(Mn4<=0), Mn4=0; end
            [sorted_formula1, total1] = Formula_Output({Ba_n*renorm, Ca_n*renorm, Sr_n*renorm, Pb_n*renorm, Mn2, Na_n*renorm, K_n*renorm},{'Ba','Ca','Sr','Pb','Mn2+','Na','K'});
            [sorted_formula2, total2] = Formula_Output({Mn3, Al_n*renorm, Fe_n*renorm},{'Mn3+','Al','Fe3+'});
            T.formula(i) = {strcat(sorted_formula1,total1,'(Mn4+_', num2str(round(Mn4,2)), sorted_formula2,total2, ')∑={',num2str(round(Mn4+Mn3+Al_n*renorm+Fe_n*renorm,2)), 'O16')};
        end

        %% 5.2.4 Mg-Hydroxide
        if(T.MgO(i) > 80)
            if(Mg_n/24 > 1-1*t && Mg_n/24 < 1+1*t)
                Renorm_factor = 24;
                Mg_n1 = T.Mg_Multi(i) / Renorm_factor;
                Mn_n1 = T.Mn_Multi(i) / Renorm_factor;
                Fe_n1 = T.Fe_Multi(i) / Renorm_factor;
                Ca_n1 = T.Ca_Multi(i) / Renorm_factor;
                Cl_n1 = T.Cl_Multi(i) / Renorm_factor;

                if(Cl_n1 > 0.0005 && Ca_n1 == 0)
                    T.group2(i) = {'Hydroxide'};
                    T.group3(i) = {'Fe-Al-Ti-Mn-Mg Hydroxide'};
                    T.group4(i) = {'Mg-Hydroxide'};
                    T.species(i) = {'brucite'};
                    [sorted_formula, total]  = Formula_Output({Mg_n1, Fe_n1, Mn_n1}, {'Mg', 'Fe', 'Mn'});
                    T.formula(i) = {strcat('(', sorted_formula, ')Σ=', total, '(OH)2')};
                end

            end
        end

        %% 5.2.5 Other Hydroxides
        if(ABCDT_n/24 > (8/8.75)-(8/8.75)*t && ABCDT_n/24 < (8/8.75)+(8/8.75)*t)
            if((Fe_n + Ni_n + Mg_n)/24 > (8/9)-(8/9)*t && (Fe_n + Ni_n + Mg_n)/24 < (8/9)+2*(8/9)*t)
                if(Cl_n > 0.1 && Fe_n > Ni_n && Ni_n > Mg_n)
                    T.group2(i) = {'Hydroxide'};
                    T.group3(i) = {'Other Hydroxide'};
                    T.group4(i) = {''};
                    T.species(i) = {'akaganeite'};
                    [sorted_formula, total] = Formula_Output({Fe_n*10.5/24, Ni_n*10.5/24, Mg_n*10.5/24, Mn_n*10.5/24, Al_n*10.5/24}, {'Fe3+','Ni2+','Mg','Mn','Al'});
                    T.formula(i) = {strcat(sorted_formula, total, '(OH,O)16Cl_1.25·nH2O')};
                end
            end
        end

        if(strcmp(T.group2(i), 'Oxide or Hydroxide') || strcmp(T.group2(i), 'Oxide') || strcmp(T.group2(i), 'Hydroxide'))
            if(strcmp(T.group3(i), 'Fe-Al-Ti-Mn-Mg Oxide or Hydroxide'))
                if(T.FeO(i) > 98 && T.CaO(i) < 0.01 && T.Cr2O3(i) == 0 && T.TiO2(i) == 0)
                    T.group4(i) = {'Fe-Oxide or Fe-Hydroxide'};
                    T.species(i) = {'magnetite/hematite/wüstite/Fe-Hydroxide'};
                    T.formula(i) = {''};
                    T.check(i) = {'no check'};
                end
            end
        end

    end %end of Oxide-Hydroxide

    %% 6. Halide Class
    if(strcmp(T.group1(i), 'Halide'))

        if(Na_n / (K_n + Ca_n + Mg_n + Fe_n) > 1)
            if(Na_n / Cl_n > 1-1*t && Na_n / Cl_n < 1+1*t)
                if(T.Na2O(i) + T.Cl(i) > 95)
                    Renorm_factor = 48;
                    Na_n05 = T.Na_Multi(i) / Renorm_factor;
                    K_n05 = T.K_Multi(i) / Renorm_factor;
                    Ca_n05 = T.Ca_Multi(i) / Renorm_factor;
                    Mg_n05 = T.Mg_Multi(i) / Renorm_factor;
                    F_n05 = T.F_Multi(i) / Renorm_factor;
                    Cl_n05 = T.Cl_Multi(i) / Renorm_factor;
                    T.species(i) = {'halite'};
                    [sorted_formula1, total1] = Formula_Output({Na_n05, K_n05, Ca_n05, Mg_n05}, {'Na', 'K','Ca', 'Mg'});
                    [sorted_formula2, total2] = Formula_Output({Cl_n05, F_n05},{'Cl','F'});
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2)};
                end
            end
        end

        if(K_n / (Na_n + Ca_n + Mg_n + Fe_n) > 1)
            if(K_n / Cl_n > 1-1*t && K_n / Cl_n < 1+1*t)
                if((T.K2O(i) + T.Cl(i)) > 95)
                    Renorm_factor = 48;
                    Na_n05 = T.Na_Multi(i) / Renorm_factor;
                    K_n05 = T.K_Multi(i) / Renorm_factor;
                    Ca_n05 = T.Ca_Multi(i) / Renorm_factor;
                    F_n05 = T.F_Multi(i) / Renorm_factor;
                    Cl_n05 = T.Cl_Multi(i) / Renorm_factor;
                    OH = 1 - (F_n05 + Cl_n05);
                    if(OH<=0), OH = 0; end
                    T.species(i) = {'sylvite'};
                    [sorted_formula1, total1]  = Formula_Output({Na_n05, K_n05, Ca_n05}, {'Na', 'K', 'Ca'});
                    [sorted_formula2, total2] = Formula_Output({Cl_n05, F_n05, OH}, {'Cl','F','OH'});
                    T.formula(i) = {strcat(sorted_formula1, total, sorted_formula2, total2)};
                end
            end
        end

        if(Mg_n / (Na_n + Ca_n + K_n + Fe_n) > 1)
            if(Mg_n / Cl_n > 0.5-0.5*t && Mg_n / Cl_n < 0.5+0.5*t)
                if((T.MgO(i) + T.Cl(i)) > 95)
                    Renorm_factor = 24;
                    Mg_n1 = T.Mg_Multi(i) / Renorm_factor;
                    Fe_n1 = T.Fe_Multi(i) / Renorm_factor;
                    Ni_n1 = T.Ni_Multi(i) / Renorm_factor;
                    Mn_n1 = T.Mn_Multi(i) / Renorm_factor;
                    Cu_n1 = T.Cu_Multi(i) / Renorm_factor;
                    F_n1 = T.F_Multi(i) / Renorm_factor;
                    Cl_n1 = T.Cl_Multi(i) / Renorm_factor;
                    OH = 2-(F_n1+Cl_n1);
                    T.species(i) = {'chloromagnesite'};
                    [sorted_formula1, total1]  = Formula_Output({Mg_n1, Fe_n1,Ni_n1,Mn_n1,Cu_n1}, {'Mg', 'Fe','Ni','Mn','Cu'});
                    [sorted_formula2, total2]  = Formula_Output({Cl_n1, F_n1, OH}, {'Cl','F','OH'});
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2)};
                end
            end
        end

        if((Fe_n + Ni_n) / (Na_n + Ca_n + K_n + Mg_n) > 1)
            if((Fe_n + Ni_n) / Cl_n > 0.5-0.5*t && (Fe_n + Ni_n) / Cl_n < 0.5+0.5*t)
                if(Fe_n/Ni_n > 1)
                    Renorm_factor = 24;
                    Mg_n1 = T.Mg_Multi(i) / Renorm_factor;
                    Fe_n1 = T.Fe_Multi(i) / Renorm_factor;
                    Mn_n1 = T.Mn_Multi(i) / Renorm_factor;
                    Ca_n1 = T.Ca_Multi(i) / Renorm_factor;
                    Cl_n1 = T.Cl_Multi(i) / Renorm_factor;
                    T.species(i) = {'lawrencite'};
                    [sorted_formula, total]  = Formula_Output({Mg_n1, Fe_n1, Mn_n1, Ca_n1}, {'Mg', 'Fe', 'Mn','Ca'});
                    T.formula(i) = {strcat(sorted_formula, total, 'Cl_',num2str(round(Cl_n1,2)))};
                end
            end
        end

        if(Fe_n / (Na_n + Ca_n + K_n + Mg_n) > 1)
            if(Fe_n / Cl_n > (1/3)-(1/3)*t && Fe_n / Cl_n < (1/3)+(1/3)*t)
                if((T.FeO(i) + T.Cl(i)) > 95)
                    Renorm_factor = 24 / 1.5;
                    OH = 3 - (F_n/Renorm_factor + Cl_n/Renorm_factor);
                    if(OH<=0), OH = 0; end
                    T.species(i) = {'molysite'};
                    [sorted_formula1, total1] = Formula_Output({Fe_n/Renorm_factor, Al_n/Renorm_factor, Cr_n/Renorm_factor},{'Fe','Al','Cr'});
                    [sorted_formula2, total2] = Formula_Output({OH, Cl_n/Renorm_factor, F_n/Renorm_factor},{'OH','Cl','F'});    
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2)};
                end
            end
        end

        if((T.MnO(i) + T.Cl(i)) > 95)
            T.group2(i) = {''};
            T.group3(i) = {''};
            T.group4(i) = {''}; 
            T.species(i) = {'scacchite'};
            [sorted_formula1, total1] = Formula_Output({Mn_n/24, Fe_n/24, Ca_n/24, Na_n/24, K_n/24},{'Mn','Fe','Ca','Na','K'});
            [sorted_formula2, total2] = Formula_Output({2 - (F_n/24 + Cl_n/24), F_n/24, Cl_n/24},{'OH','Cl','F'});
            T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2, total2)};
        end

        if((Ca_n+Na_n)/Cl_n >= (1/2)-2*(1/2)*t && (Ca_n+Na_n)/Cl_n <= (1/2)+2*(1/2)*t)
            if(Ca_n > Na_n && Ca_n/(Ca_n+Na_n+Mg_n+Fe_n+Mn_n) > 0.5)
                T.species(i) = {'antarcticite'};
                OH = 2 - (Cl_n/24 + F_n/24);
                if(OH<=0), OH=0; end
                [sorted_formula1,total1] = Formula_Output({Ca_n/24,Na_n/24,K_n/24,Mg_n/24,Fe_n/24},{'Ca','Na','K','Mg','Fe'});
                [sorted_formula2,total2] = Formula_Output({Cl_n/24, F_n/24, OH},{'Cl','F','OH'});
                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'·6H2O')};
            end
        end

    end %end of Halides

    %% 7. Silicate class
    if(strcmp(T.group1(i), 'Silicate') && ~strcmp(T.group1(i), 'SiO2 Phase'))

        %% 7.1 Staurolite
        if((Fe_n+Mg_n+Mn_n+Zn_n)/Al_n >= 2/9 - (2/9)*3*t && (Fe_n+Mg_n+Mn_n+Zn_n)/Al_n < 2/9 + (2/9)*3*t)
            if(Si_n/Al_n >= 4/9 - 2*(4/9)*t && Si_n/Al_n < 4/9 + 2*(4/9)*t)
                if((Na_n+Ca_n+K_n)/(Fe_n+Mg_n+Mn_n+Zn_n) < 0.15 && Fe_n > Mg_n && S_n <= 0.05 && P_n <= 0.05)
                    T.group2(i) = {'Nesosilicate'};
%                     T.group3(i) = {'Nesosilicate'};
                    T.group3(i) = {'Staurolite'};
                    T.species(i) = {'staurolite'};
                    renorm = 23.5/24;
                    Al_IV = 4 - Si_n*renorm;
                    if(Al_IV < 0), Al_IV = 0; end
                    Al_VI = Al_n*renorm - Al_IV;
                    if(Al_VI < 0), Al_VI = 0; end
                    Fe3 = 9 - (Al_VI + Cr_n*renorm + Ti_n*renorm);
                    if(Fe3 < 0), Fe3 = 0; end
                    Fe2 = Fe_n*renorm - Fe3;
                    if(Fe2 < 0), Fe2 = 0; end
                    [sorted_formula1, total1] = Formula_Output({Fe2, Mg_n*renorm, Mn_n*renorm, Zn_n*renorm}, {'Fe2+','Mg','Mn','Zn'});
                    [sorted_formula2, total2] = Formula_Output({Al_VI, Fe3, Cr_n*renorm, Ti_n*renorm}, {'AlVI','Fe3+','Cr','Ti'});
                    [sorted_formula3, total3] = Formula_Output({Si_n*renorm, Al_IV}, {'Si','AlIV'});
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'[',sorted_formula3,total3, 'O23](OH)')};
                end
            end
        end

        %% 7.2 Cyclosilicate
        if((Fe_n+Mg_n+Mn_n+Zn_n)/(Al_n+Cr_n+Ti_n) >= 2/4 - (2/4)*2*t && (Fe_n+Mg_n+Mn_n+Zn_n)/(Al_n+Cr_n+Ti_n) <= 2/4 + (2/4)*2*t)
            if(Si_n/Al_n >= 5/4 - 2*(5/4)*t && Si_n/Al_n < 5/4 + 2*(5/4)*t)
                if((Na_n+Ca_n+K_n)/(Fe_n+Mg_n+Mn_n+Zn_n) < 0.2 && S_n <= 0.05 && P_n <= 0.05)

                    if(Mg_n > Fe_n)
                        T.group2(i) = {'Cyclosilicate'};
                        T.group3(i) = {'5-Ring Silicate'};
                        T.species(i) = {'cordierite'};                    
                    end

                    if(Fe_n > Mg_n)
                        T.group2(i) = {'Cyclosilicate'};
                        T.group3(i) = {'5-Ring Silicate'};
                        T.species(i) = {'sekaninaite'};  
                    end

                    renorm = 18/24;
                    Al_IV = 5 - Si_n*renorm;
                    if(Al_IV < 0), Al_IV = 0; end
                    Al_VI = Al_n*renorm - Al_IV;
                    if(Al_VI < 0), Al_VI = 0; end
                    Fe3 = 2 - (Al_VI + Cr_n*renorm + Ti_n*renorm);
                    if(Fe3 < 0), Fe3 = 0; end
                    Fe2 = Fe_n*renorm - Fe3;
                    if(Fe2 < 0), Fe2 = 0; end
                    [sorted_formula1, total1] = Formula_Output({Fe2, Mg_n*renorm, Mn_n*renorm, Zn_n*renorm}, {'Fe2+','Mg','Mn','Zn'});
                    [sorted_formula2, total2] = Formula_Output({Al_VI, Fe3, Cr_n*renorm, Ti_n*renorm}, {'AlVI','Fe3+','Cr','Ti'});
                    [sorted_formula3, total3] = Formula_Output({Si_n*renorm, Al_IV}, {'Si','AlIV'});
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'[',sorted_formula3,total3, 'O18]')};
                end
            end
        end
        
        %% 7.3 Framework Silicates
        if((Si_n + Al_n) / 24 > 0.5-0.5*t && (Si_n + Al_n) / 24 < 0.5+0.5*t)
            if((Mg_n + Fe_n + Ti_n + Cr_n) < 0.5 && S_n <= 0.05 && P_n <= 0.05)
                if(ABCDT_n < 15.8 || (ABCDT_n >= 15.8 && K_n > Na_n) || 18-18*t < ABCDT_n)
                    T.group2(i) = {'Tectosilicate'};

                    if( (T.Al2O3(i) > 50 && Ca_n/(K_n+Ca_n+Na_n)>0.5) ...
                            || (T.SiO2(i) < 50 && T.K2O(i) > 9 && T.K2O(i) < 13) ...
                            || (T.SiO2(i) < 50 && T.Na2O(i) > 5 && T.Na2O(i) < 9) )
                        T.group2(i) = {'Phyllosilicate'};
                        T.group3(i) = {'Mica'};
                    end

                    if((Na_n + K_n + Ca_n + Sr_n + Ba_n)/(Si_n + Al_n) < 0.05/6 && Al_n/Si_n >= (1/2)-(1/2)*t && Al_n/Si_n < (1/2)+(1/2)*t)
                        T.group2(i) = {'Phyllosilicate'};
                        T.group3(i) = {'Talc-Pyrophyllite'};    
                        T.group4(i) = {''};
                        T.species(i) = {''};
                        OH = 2 - (F_n11+Cl_n11);
                        if(OH<=0), OH = 0; end
                        [sorted_formula1, total1] = Formula_Output({Fe_n11, Mg_n11, Ni_n11, Mn_n11}, {'Fe', 'Mg', 'Ni','Mn'});
                        [sorted_formula2, total2] = Formula_Output({Si_n11, Al_n11},{'Si','Al'});
                        [sorted_formula3, total3] = Formula_Output({OH, F_n11, Cl_n11},{'OH','F','Cl'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O10',sorted_formula3, total3)};
                    end

                    if( (Na_n + K_n + Ca_n)/(Si_n+Al_n) > (0.2/6)-(0.2/6)*t && (Na_n + K_n + Ca_n)/(Si_n+Al_n) < (0.6/6)+(0.6/6)*t && (Ca_n + Na_n)/(Ca_n+Na_n+K_n) > 0.5 && T.SiO2(i) < 75)
                        T.group2(i) = {'Phyllosilicate'};
                        T.group3(i) = {'Smectite'};
                        T.group4(i) = {'Dioctahedral Smectite'};
                        T.species(i) = {'beidellite'};
                        Fe3 = Fe_n*(11/24);
                        Al_IV = 4 - Si_n*(11/24);
                        if(Al_IV < 0), Al_IV = 0; end
                        Al_VI = Al_n*(11/24) - Al_IV;
                        if(Al_VI < 0), Al_VI = 0; end
                        [sorted_formula1, total1]  = Formula_Output({Na_n*(11/24), Ca_n*(11/24),K_n*(11/24)}, {'Na', 'Ca','K'});
                        [sorted_formula2, total2]  = Formula_Output({Al_VI, Fe3}, {'AlVI', 'Fe3+'});
                        [sorted_formula3, total3]  = Formula_Output({Si_n*(11/24), Al_IV}, {'Si', 'AlIV'});
                        [sorted_formula4, total4]  = Formula_Output({2-(Cl_n*(11/24)+F_n*(11/24)), Cl_n*(11/24),F_n*(11/24)}, {'OH', 'Cl','F'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'[',sorted_formula3,total3,'O10]',sorted_formula4,total4,'·nH2O')};
                    end

                    %% 7.1.1 Feldspar Group
                    if(~strcmp(T.group2(i), 'Phyllosilicate') && (Na_n + K_n + Ca_n + Al_n + Si_n) / 24 > (5/8)-(5/8)*t && (Na_n + K_n + Ca_n + Al_n + Si_n) / 24 < (5/8)+(5/8)*t)
                        if((Na_n + K_n + Ca_n + Al_n + Si_n) > 15-15*t && (Na_n + K_n + Ca_n + Al_n + Si_n) < 15+15*t)
                            if((Si_n+Al_n)/24 > 0.5-0.5*t && (Si_n+Al_n)/24 < 0.5+0.5*t && round(Ca_n/Al_n,2) <= 0.5+0.5*2*t)
                                T.group3(i) = {'Feldspar'};
                                T.group4(i) = {''};
                                T.species(i) = {''};
                                T.formula(i) = {''};
                                Renorm_factor = 3;
                                Na_n8 = T.Na_Multi(i) / Renorm_factor;
                                K_n8 = T.K_Multi(i) / Renorm_factor;
                                Ca_n8 = T.Ca_Multi(i) / Renorm_factor;
                                Al_n8 = T.Al_Multi(i) / Renorm_factor;
                                Fe_n8 = T.Fe_Multi(i) / Renorm_factor;
                                Si_n8 = T.Si_Multi(i) / Renorm_factor;
                                Mg_n8 = T.Mg_Multi(i) / Renorm_factor;
                                Ba_n8 = T.Ba_Multi(i) / Renorm_factor;
                                Sr_n8 = T.Sr_Multi(i) / Renorm_factor;
                                Ab = (Na_n / (Na_n + K_n + Ca_n)) * 100;
                                An = (Ca_n / (Na_n + K_n + Ca_n)) * 100;
                                Or = (K_n / (Na_n + K_n + Ca_n)) * 100;
                                T.An(i) = round(An,2);
                                T.Or(i) = round(Or,2);
                                T.Ab(i) = round(Ab,2);

                                if(An >= 0 && An < 10 && Or >= 40)
                                    T.species(i) = {'K-feldspar'};
                                    [sorted_formula1, total1] = Formula_Output({Na_n8,K_n8,Ca_n8, Mg_n8}, {'Na', 'K', 'Ca', 'Mg'});
                                    [sorted_formula2, total2] = Formula_Output({Al_n8, Fe_n8},{'Al','Fe'});
                                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'Si_',num2str(round(Si_n8,2)),'O8')};
                                    if(An >= 0 && An < 10 && Or >= 80 && Or <= 100)
                                        T.species(i) = {'sanidine/orthoclase/microcline'};
                                    end
                                end

                                if(An < 10 && Or < 10)
                                    T.species(i) = {'albite'};
                                    [sorted_formula1, total1] = Formula_Output({Na_n8,K_n8,Ca_n8}, {'Na', 'K', 'Ca'});
                                    [sorted_formula2, total2] = Formula_Output({Al_n8,Fe_n8}, {'Al', 'Fe3+'});
                                    T.formula(i) = {strcat(sorted_formula1, total1, '[', sorted_formula2,total2, 'Si_', num2str(round(Si_n8,2)), 'O8]')};
                                end

                                if(An > 90 && (Ca_n + Al_n)/Si_n > 3/2 - 2*(3/2)*t && (Ca_n + Al_n)/Si_n < 3/2 + (3/2)*t)
                                    T.species(i) = {'anorthite'};
                                    [sorted_formula, total] = Formula_Output({Na_n8,K_n8,Ca_n8,Ba_n8,Sr_n8}, {'Na', 'K', 'Ca','Ba','Sr'});
                                    T.formula(i) = {strcat(sorted_formula,total,'[Al_', num2str(round(Al_n8)),'Si_',num2str(round(Si_n8)), 'O8]')};
                                end

                                if(An > 10 && An < 90)
                                    T.species(i) = {'plagioclase'};
                                    [sorted_formula1, total1] = Formula_Output({Na_n8,K_n8,Ca_n8}, {'Na', 'K', 'Ca'});
                                    [sorted_formula2, total2] = Formula_Output({Al_n8, Fe_n8, Mg_n8}, {'Al','Fe','Mg'});
                                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'Si_', num2str(round(Si_n8,2)),'O8')};

                                    if(An >= 80 && An < 100 && Or > 0 && Or < 10)
                                        T.species(i) = {'plagioclase (HTLP or High H2O)'};
                                    end

                                    if(An >= 60 && An < 80 && Or > 10 && Or <= 20)
                                        T.species(i) = {'plagioclase (HTMP or High H2O)'};
                                    end

                                    if(An >= 45 && An <= 60 && Or > 20 && Or < 30)
                                        T.species(i) = {'plagioclase (HT & High H2O)'};
                                    end
                                    if(Or > 10)
                                        T.species(i) = {'ternary feldspars'};
                                    end

                                end

                                if((Ab + Or) > An)
                                    if(Ab < 90)
                                        if(Or > An && An <= 10 && Or < 40 && Or > 10)
                                            T.species(i) = {'anorthoclase'};
                                            [sorted_formula1, total1] = Formula_Output({Na_n8,K_n8,Ca_n8,Ba_n8,Sr_n8}, {'Na', 'K', 'Ca','Ba','Sr'});
                                            [sorted_formula2,total2] = Formula_Output({Al_n8, Fe_n8},{'Al','Fe3+'});
                                            T.formula(i) = {strcat(sorted_formula1,total1,'[', sorted_formula2,total2,'Si_',num2str(round(Si_n8)), 'O8]')};
                                        end
                                    end
                                end
                
                                % Feldsapr miscibility gap
                                if(An >= 80 && An < 100 && Or > 10 || ...
                                        An >= 60 && An < 80 && Or > 20 || ...
                                        An > 45 && An < 60 && Or > 30 || ...
                                        Or >= 80 && Or < 100 && An > 10 || ...
                                        Or >= 65 && Or < 80 && An > 20 || ...
                                        Or >= 10 && Or < 65 && An > 30)
                                    T.formula(i) = {''};
                                    T.species(i) = {''};
                                    T.group4(i) = {''};
                                    T.group3(i) = {''};
                                end    

                            end
                        end
                    end %end of Feldspars

                    %% 7.1.2 Feldspathoid Group
                    if(~strcmp(T.group2(i), 'Phyllosilicate') && (Na_n + K_n + Ca_n + Al_n + Si_n) / 24 > 0.66)
                        if((Si_n + Al_n) / 24 > 0.5-0.5*t && (Si_n + Al_n) / 24 < 0.5+0.5*t)

                            if((Na_n+Ca_n+K_n)/(Si_n+Al_n) > (1/2)-2*(1/2)*t && (Na_n+Ca_n+K_n)/(Si_n+Al_n) < (1/2)+2*(1/2)*t)
                                if((Na_n+Ca_n)/K_n > 3*(16/24)-2*(3*(16/24))*t && (Na_n+Ca_n)/K_n < 7*(16/24)+2*(7*(16/24))*t)
                                    if(Na_n > Ca_n && K_n > Ca_n)
                                        T.group3(i) = {'Feldspathoid'};
                                        T.group4(i) = {''};
                                        T.species(i) = {'nepheline'};
                                        Renorm_factor = 24 / 16;
                                        Na_n16 = T.Na_Multi(i) / Renorm_factor;
                                        Ca_n16 = T.Ca_Multi(i) / Renorm_factor;
                                        Al_n16 = T.Al_Multi(i) / Renorm_factor;
                                        Si_n16 = T.Si_Multi(i) / Renorm_factor;

                                        [sorted_formula1,total1] = Formula_Output({Na_n16, Ca_n16},{'Na','Ca'});
                                        [sorted_formula2, total2] = Formula_Output({Al_n16, Si_n16},{'Al','Si'});
                                        T.formula(i) = {strcat(sorted_formula1,total1, 'K_',sorted_formula2,total2, 'O16)')};
                                    end
                                end
                            end

                            if((Na_n+Ca_n+K_n)/(Si_n+Al_n) > (1/2)-2*(1/2)*t && (Na_n+Ca_n+K_n)/(Si_n+Al_n) < (1/2)+2*(1/2)*t)
                                if(K_n/(Ca_n+Na_n) < (0.5/3.5))
                                    if(Na_n > Ca_n && Na_n > K_n)
                                        T.group3(i) = {'Feldspathoid'};
                                        T.group4(i) = {''};
                                        T.species(i) = {'trinepheline'};
                                        Renorm_factor = 24 / 16;
                                        Na_n16 = T.Na_Multi(i) / Renorm_factor;
                                        Ca_n16 = T.Ca_Multi(i) / Renorm_factor;
                                        Al_n16 = T.Al_Multi(i) / Renorm_factor;
                                        Si_n16 = T.Si_Multi(i) / Renorm_factor;

                                        [sorted_formula1,total1] = Formula_Output({Na_n16, Ca_n16},{'Na','Ca'});
                                        [sorted_formula2, total2] = Formula_Output({Al_n16, Si_n16},{'Al','Si'});
                                        T.formula(i) = {strcat(sorted_formula1,total1, 'K_',sorted_formula2,total2, 'O16)')};
                                    end
                                end
                            end

                            if((Na_n + K_n) / (Al_n + Si_n) > (0.666)-(0.666)*t && (Na_n + K_n) / (Al_n + Si_n) < (0.666)+(0.666)*t)
                                if(Na_n > K_n && Cl_n > 1)
                                    T.group3(i) = {'Feldspathoid'};
                                    T.group4(i) = {''};
                                    T.species(i) = {'sodalite'};
                                    Renorm_factor = 24 / 12.5;
                                    Na_n12 = T.Na_Multi(i) / Renorm_factor;
                                    K_n12 = T.K_Multi(i) / Renorm_factor;
                                    Ca_n12 = T.Ca_Multi(i) / Renorm_factor;
                                    Al_n12 = T.Al_Multi(i) / Renorm_factor;
                                    Si_n12 = T.Si_Multi(i) / Renorm_factor;
                                    Cl_n12 = T.Cl_Multi(i) / Renorm_factor;
                                    F_n12 = T.F_Multi(i) / Renorm_factor;

                                    OH = 1 - (Cl_n12+F_n12);
                                    if(OH<=0), OH = 0; end
                                    [sorted_formula, total]  = Formula_Output({Na_n12, Ca_n12, K_n12}, {'Na', 'Ca', 'K'});
                                    [sorted_formula2, total2] = Formula_Output({Si_n12, Al_n12},{'Si','Al'});
                                    [sorted_formula3, total3] = Formula_Output({Cl_n12, OH, F_n12}, {'Cl','OH','F'});
                                    T.formula(i) = {strcat(sorted_formula, total, sorted_formula2, total2, 'O12',sorted_formula3, total3)};
                                end
                            end

                            if((Na_n + K_n) / (Al_n + Si_n) > (0.333)-(0.333)*t && (Na_n + K_n) / (Al_n + Si_n) < (0.333)+(0.333)*t)
                                if(K_n > Na_n)
                                    if(Al_n / Si_n > 0.5-0.5*t && Al_n / Si_n < 0.5+0.5*t)
                                        T.group3(i) = {'Feldspathoid'};
                                        T.group4(i) = {''};
                                        T.species(i) = {'leucite'};
                                        Renorm_factor = 4;
                                        Na_n6 = T.Na_Multi(i) / Renorm_factor;
                                        K_n6 = T.K_Multi(i) / Renorm_factor;
                                        Al_n6 = T.Al_Multi(i) / Renorm_factor;
                                        Si_n6 = T.Si_Multi(i) / Renorm_factor;
                                        [sorted_formula1, total1] = Formula_Output({Na_n6, K_n6}, {'Na', 'K'});
                                        [sorted_formula2, total2] = Formula_Output({Al_n6, Si_n6}, {'Al','Si'});
                                        T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O6')};
                                    end
                                end
                            end

                            if((Na_n + Ca_n) / (Si_n + Al_n) > (0.666)-(0.666)*t && (Na_n + Ca_n) / (Si_n + Al_n) < (0.666)+(0.666)*t)
                                if(Na_n > Ca_n && Ca_n/(Na_n+K_n) > (1/4)-2*(1/4)*t && Ca_n/(Na_n+K_n) < (1/3)+2*(1/3)*t)
                                    if(S_n >= 0)
                                        T.group3(i) = {'Feldspathoid'};
                                        T.group4(i) = {''};
                                        T.species(i) = {'cancrinite'};
                                        Renorm = 26/24;
                                        vac = 8 - (Na_n*Renorm + Ca_n*Renorm + K_n*Renorm);
                                        if(vac < 0), vac = 0; end
                                        [sorted_formula1, total1]  = Formula_Output({Na_n, Ca_n, K_n,vac}, {'Na', 'Ca', 'K','□'});
                                        [sorted_formula2,total2] = Formula_Output({Al_n*Renorm,Fe_n*Renorm},{'Al','Fe3+'});
                                        T.formula(i) = {strcat(sorted_formula1,total1,'[',sorted_formula2,total2,'Si_',num2str(round(Si_n*Renorm)), 'O24](CO3,SO4,PO4)2·2H2O')};
                                    end
                                end
                            end

                        end
                    end %end of Feldspathoids

                    %% 7.1.3 Zeolite group
                    if(strcmp(T.group2(i),'Tectosilicate') && ~strcmp(T.group1(i), 'SiO2 Phase'))
                        if(Al_n/Si_n >= (1/5)-(1/5)*t && Al_n/Si_n <= (1/3)+2*(1/3)*t)
                            if( ((Na_n+K_n+Ca_n)/(Si_n+Al_n) > (1/12)-(1/12)*t && (Na_n+K_n+Ca_n)/(Si_n+Al_n) < (1/4)+(1/4)*t) )

                                Renorm_factor = 24/72;
                                Na_n72 = T.Na_Multi(i) / Renorm_factor;
                                K_n72 = T.K_Multi(i) / Renorm_factor;
                                Ca_n72 = T.Ca_Multi(i) / Renorm_factor;
                                Si_n72 = T.Si_Multi(i) / Renorm_factor;
                                Al_n72 = T.Al_Multi(i) / Renorm_factor;
                                Mg_n72 = T.Mg_Multi(i) / Renorm_factor;
                                Fe_n72 = T.Fe_Multi(i) / Renorm_factor;
                                Mn_n72 = T.Mn_Multi(i) / Renorm_factor;

                                if(Na_n + K_n + Ca_n > 0 && Na_n + K_n + Ca_n < 2.5 + 2.5*t)
                                    T.group3(i) = {'Zeolite'};
                                    T.group4(i) = {''};
                                    T.species(i) = {''};
                                    T.formula(i) = {''};
                                end
                            end
                        end
                    end

                    if(strcmp(T.group3(i), 'Zeolite'))
                        if((Al_n+Fe_n) / Si_n > (6/30)-(6/30)*t && (Al_n+Fe_n) / Si_n < (6.5/30)+(6.5/30)*t)
                            if(Ca_n > Na_n && Ca_n > K_n && Ca_n > Mg_n && (Ca_n+Na_n+K_n)*(72/24) > 2-2*t && Ca_n*(72/24) < 3+3*t)
                                T.species(i) = {'clinoptilolite-Ca'};
                                [sorted_formula1, total1] = Formula_Output({Ca_n72, Na_n72, K_n72},{'Ca','Na','K'});
                                [sorted_formula2, total2] = Formula_Output({Si_n72,Al_n72},{'Si','Al'});
                                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O72·20H2O')};
                            end
                        end

                        if(Al_n / Si_n > (0.2)-(0.2)*t && Al_n / Si_n < (0.22)+(0.22)*t)
                            if(Na_n > K_n && Na_n > Ca_n && (Ca_n+Na_n+K_n)*(72/24) > 2-2*t && Ca_n*(72/24) < 3+3*t)
                                T.species(i) = {'clinoptilolite-Na'};
                                [sorted_formula1, total1] = Formula_Output({Na_n72, Ca_n72, K_n72},{'Na','Ca','K'});
                                [sorted_formula2, total2] = Formula_Output({Si_n72,Al_n72},{'Si','Al'});
                                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O72·20H2O')};
                            end
                        end

                        if(Al_n / Si_n > (0.2)-(0.2)*t && Al_n / Si_n < (0.22)+(0.22)*t)
                            if(K_n > Na_n && K_n > Ca_n && (Ca_n+Na_n+K_n)*(72/24) > 2-2*t && Ca_n*(72/24) < 3+3*t)
                                T.species(i) = {'clinoptilolite-K'};
                                [sorted_formula1, total1] = Formula_Output({K_n72, Na_n72, Ca_n72},{'K','Na','Ca'});
                                [sorted_formula2, total2] = Formula_Output({Si_n72,Al_n72},{'Si','Al'});
                                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O72·20H2O')};
                            end  
                        end

                        if(Al_n / Si_n > (0.333)-2*(0.333)*t && Al_n / Si_n < (0.333)+2*(0.333)*t)
                            if(Ca_n > Na_n && Ca_n > K_n && Ca_n/(Na_n+K_n) < 3)
                                T.species(i) = {'heulandite-Ca'};
                                [sorted_formula1, total1] = Formula_Output({Ca_n72, Na_n72, K_n72, Mg_n72, Mn_n72}, {'Ca','Na','K', 'Mg', 'Mn'});
                                [sorted_formula2, total2] = Formula_Output({Fe_n72, Al_n72},{'Fe3+','Al'});
                                T.formula(i) = {strcat(sorted_formula1,total1, '(Si_', num2str(round(Si_n72,2)), sorted_formula2,total2,')∑={',num2str(round(Si_n72+Fe_n72+Al_n72)),'O72·26H2O')};
                            end
                        end

                        if((Na_n + K_n + Ca_n + Ba_n + Sr_n) / (Si_n + Al_n) > (0.1666)-(0.1666)*t && (Na_n + K_n + Ca_n + Ba_n + Sr_n) / (Si_n + Al_n) < (0.1666)+(0.1666)*t)
                            if(Al_n / Si_n > (1/3)-2*(1/3)*t && Al_n / Si_n < (1/3)+2*(1/3)*t)
                                if(Na_n / (Na_n + K_n + Ca_n) > 0.5-0.5*t && Na_n / (Na_n + K_n + Ca_n) < 0.75+0.75*t && (Na_n*(72/24) + K_n*(72/24) + Ca_n*(72/24)) > 6-6*t && (Na_n*(72/24) + K_n*(72/24) + Ca_n*(72/24)) < 6+6*t )
                                    T.species(i) = {'heulandite-Na'};
                                    [sorted_formula1, total1] = Formula_Output({Ca_n72, Na_n72, K_n72, Mg_n72, Mn_n72}, {'Ca','Na','K', 'Mg', 'Mn'});
                                    [sorted_formula2, total2] = Formula_Output({Si_n72, Fe_n72, Al_n72},{'Si','Fe3+','Al'});
                                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2,'O72·22H2O')};
                                end
                            end
                        end

                        if((Na_n + K_n + Ca_n) / (Si_n + Al_n) > (0.1388)-(0.1388)*t && (Na_n + K_n + Ca_n) / (Si_n + Al_n) < (0.1388)+(0.1388)*t)
                            if(Al_n / Si_n > (0.333)-(0.333)*t && Al_n / Si_n < (0.333)+(0.333)*t)
                                if(K_n / (Na_n + K_n + Ca_n) > 0.5-0.5*t && K_n / (Na_n + K_n + Ca_n) < 0.75+0.75*t)
                                    T.species(i) = {'heulandite-K'};
                                    [sorted_formula1, total1] = Formula_Output({Ca_n72, Na_n72, K_n72, Mg_n72, Mn_n72}, {'Ca','Na','K', 'Mg', 'Mn'});
                                    [sorted_formula2, total2] = Formula_Output({Fe_n72, Al_n72},{'Fe3+','Al'});
                                    T.formula(i) = {strcat(sorted_formula1,total1, '(Si_', num2str(round(Si_n72,2)), sorted_formula2,total2,')∑={',num2str(round(Si_n72+Fe_n72+Al_n72)),'O72·26H2O')};
                                end
                            end
                        end

                        if(Al_n / Si_n > (9/27)-2*(9/27)*t && Al_n / Si_n < (9/27)+ 2*(9/27)*t)
                            if(Ca_n / (Na_n + K_n) > 4 - 2*4*t && Ca_n / (Na_n + K_n) < 4 + 2*4*t)
                                T.species(i) = {'stilbite-Ca'};
                                vac = 1 - (Na_n72 + K_n72);
                                if(vac <= 0), vac = 0; end 
                                [sorted_formula1, total1] = Formula_Output({Na_n72, K_n72, vac}, {'Na','K','□'});
                                [sorted_formula2, total2] = Formula_Output({Ca_n72, Mn_n72, Mg_n72, Fe_n72}, {'Ca','Mn','Mg','Fe'});
                                [sorted_formula3, total3] = Formula_Output({Si_n72,Al_n72},{'Si','Al'});
                                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3, total3, 'O72·28H2O')};
                            end
                        end

                        if((Na_n + K_n + Ca_n) / (Si_n + Al_n) > (9/27)-2*(9/27)*t && (Na_n + K_n + Ca_n) / (Si_n + Al_n) < (9/27)+ 2*(9/27)*t)
                            if(Al_n / Si_n > (9/27)-2*(9/27)*t && Al_n / Si_n < (9/27)+ 2*(9/27)*t)
                                if(Na_n/Ca_n > 1 && Na_n > K_n)
                                    T.species(i) = {'stilbite-Na'};
                                    [sorted_formula1, total1] = Formula_Output({Na_n72, K_n72, Ca_n72}, {'Na','K','Ca'});
                                    [sorted_formula2, total2] = Formula_Output({Si_n72,Al_n72},{'Si','Al'});
                                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O72·20H2O')};
                                end
                            end
                        end
                    end %end of Zeolites

                end
            end
        end %end of Framework Silicates
       
        %% 7.4 Nesosilicate
        %% 7.4.1 Olivine Group
        if(ABCDT_n > 18-18*t && ABCDT_n < 18+18*t && Si_n > 6-6*t && Si_n < 6+6*t)
            if(round((Na_n + K_n) / Si_n,2) <= 0.02)
                if(Al_n / Si_n < 0.02 && Ti_n / Si_n < 0.2)
                    if((Mg_n + Fe_n + Ni_n + Mn_n + Ca_n) / Si_n > 2-2*t && (Mg_n + Fe_n + Ni_n + Mn_n + Ca_n) / Si_n < 2+2*t)
                        if((Mg_n + Fe_n) / (Mg_n + Fe_n + Ni_n + Mn_n + Ca_n) > 0.75)
                            if(Ca_n / Si_n <= 0.1+0.1*t && S_n <= 0.05 && P_n <= 0.05)
                                T.group2(i) = {'Nesosilicate'};
                                T.group3(i) = {'Olivine'};
                                T.group4(i) = {''};
                                T.species(i) = {''};
                                T.formula(i) = {''};
                                Renorm_factor = 6;
                                Mg_n4 = T.Mg_Multi(i) / Renorm_factor;
                                Mn_n4 = T.Mn_Multi(i) / Renorm_factor;
                                Fe_n4 = T.Fe_Multi(i) / Renorm_factor;
                                Ca_n4 = T.Ca_Multi(i) / Renorm_factor;
                                Si_n4 = T.Si_Multi(i) / Renorm_factor;
                                Ni_n4 = T.Ni_Multi(i) / Renorm_factor;
                                
                                Fo1 = 100 * Mg_n / (Mg_n + Fe_n + Mn_n + Ca_n);
                                T.Fo1(i) = round(Fo1,2);
                                Fo = 100 * Mg_n / (Mg_n + Fe_n);
                                T.Fo(i) = round(Fo,2);
                                Fa = 100* Fe_n / (Mg_n + Fe_n);
                                T.Fa(i) = round(Fa,2);

                                if(Mg_n > Fe_n)
                                    T.species(i) = {'forsterite'};
                                    [sorted_formula1, total1]  = Formula_Output({Mg_n4, Fe_n4, Ni_n4, Mn_n4, Ca_n4}, {'Mg','Fe','Ni','Mn','Ca'});
                                    T.formula(i) = {strcat(sorted_formula1,total1,'(Si_',num2str(round(Si_n4,2)), 'O4)')};
                                end

                                if(Mg_n <= Fe_n)
                                    T.species(i) = {'fayalite'};
                                    [sorted_formula1, total1]  = Formula_Output({Fe_n4, Mg_n4, Mn_n4, Ca_n4}, {'Fe','Mg','Mn','Ca'});
                                    T.formula(i) = {strcat(sorted_formula1,total1,'(Si_',num2str(round(Si_n4,2)), 'O4)')};
                                end

                            end
                        end
                    end
                end
            end
        end %end of olivines

        %% 7.4.2 Garnet Group

        Renorm_factor = 2;
        Mg_n12 = T.Mg_Multi(i) / Renorm_factor;
        Mn_n12 = T.Mn_Multi(i) / Renorm_factor;
        Fe_n12 = T.Fe_Multi(i) / Renorm_factor;
        Ca_n12 = T.Ca_Multi(i) / Renorm_factor;
        Si_n12 = T.Si_Multi(i) / Renorm_factor;
        Na_n12 = T.Na_Multi(i) / Renorm_factor;
        K_n12 = T.K_Multi(i) / Renorm_factor;
        Cr_n12 = T.Cr_Multi(i) / Renorm_factor;
        Ti_n12 = T.Ti_Multi(i) / Renorm_factor;
        Al_n12 = T.Al_Multi(i) / Renorm_factor;
        V_n12 = T.V_Multi(i) / Renorm_factor;
        Zr_n12 = T.Zr_Multi(i) / Renorm_factor;
        
        Al_IV = 3 - Si_n12;
        if(Al_IV<0.001), Al_IV = 0; end

        Al_VI = Al_n12 - Al_IV;
        if(Al_VI<0), Al_VI = 0; end

        Fe3 = 2 - (Al_VI + Cr_n12 + Ti_n12 + V_n12);
        if(Fe3<0.001), Fe3 = 0; end

        Fe2 = Fe_n12 - Fe3; 
        if(Fe2<0), Fe2 = 0; end

        % criteria for Zr-rich garnet kimzeyite
        if(~strcmp(T.group2(i),'Tectosilicate') && (Zr_n12 > 1 || Ti_n12 > 1))
            
            Fe3 = Fe_n12; 

            Al_IV = 3 - (Si_n12 + Fe3);
            if(Al_IV<0), Al_IV = 0; end

            Al_VI = Al_n12 - Al_IV; 
            if(Al_VI<0), Al_VI = 0; end

            Mg_b = 2 - (Zr_n12 + Ti_n12 + V_n12 + Al_VI);
            if(Mg_b<0), Mg_b = 0; end

            Mg_a = Mg_n12 - Mg_b;
            if(Mg_a<0), Mg_a = 0; end

            if((Zr_n12+Ti_n12+V_n12+Mg_b+Al_VI)/(Ca_n12 + Mg_a) > (2/3)-(2/3)*2*t && (Zr_n12+Ti_n12+V_n12+Mg_b+Al_VI)/(Ca_n12 + Mg_a) < (2/3)+(2/3)*t) 
                if((Zr_n12+Ti_n12+V_n12+Mg_b)/(Si_n12 + Al_IV + Fe3) > (2/3)-(2/3)*t && (Zr_n12+Ti_n12+V_n12+Mg_b)/(Si_n12 + Al_IV + Fe3) < (2/3)+(2/3)*t && Zr_n12 > Ti_n12)
                    T.group2(i) = {'Nesosilicate'};
                    T.group3(i) = {'Garnet'};
                    T.group4(i) = {'Other Garnet'};
                    T.species(i) = {'kimzeyite'};

                    [sorted_formula1, total1] = Formula_Output({Mg_n12, Ca_n12, Mn_n12}, {'Mg', 'Ca', 'Mn'});
                    [sorted_formula2, total2] = Formula_Output({Al_VI, Ti_n12, Cr_n12, V_n12, Zr_n12}, {'AlVI', 'Ti', 'Cr', 'V', 'Zr'});
                    [sorted_formula3, total3] = Formula_Output({Si_n12, Al_IV, Fe3}, {'Si','AlIV', 'Fe3+'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, 'O12')};

                    if(Si_n12 + Al_IV + Fe3 >= 3)
                        Fe3_b = 3 - (Si_n12 + Al_IV);
                        if(Fe3_b < 0), Fe3_b = 0; end

                        Fe3_a = Fe3 - Fe3_b;
                        if(Fe3_a < 0), Fe3_a = 0; end

                        Mg_b = 2 - (Fe3_a + Zr_n12 + Ti_n12 + V_n12 + Al_VI);
                        if(Mg_b<0), Mg_b = 0; end
        
                        Mg_a = Mg_n12 - Mg_b;
                        if(Mg_a<0), Mg_a = 0; end

                        [sorted_formula1, total1] = Formula_Output({Mg_a, Ca_n12, Mn_n12}, {'Mg', 'Ca', 'Mn'});
                        [sorted_formula2, total2] = Formula_Output({Mg_b, Al_VI, Fe3_a, Ti_n12, Cr_n12, V_n12, Zr_n12}, {'Mg','AlVI', 'Fe3+', 'Ti', 'Cr', 'V', 'Zr'});
                        [sorted_formula3, total3] = Formula_Output({Si_n12, Al_IV, Fe3_b}, {'Si','AlIV', 'Fe3+'});
                        T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, 'O12')};
                    end
                end
                
                if(Ti_n12 > Zr_n12)
                    T.group2(i) = {'Nesosilicate'};
                    T.group3(i) = {'Garnet'};
                    T.group4(i) = {'Other Garnet'};
                    T.species(i) = {'schorlomite'};
                    
                    % Ideal site capacities
                    Z_ideal = 3.00;  % Z site can hold 3 cations (shared by Si, Ti, Fe3, Al)
                    Y_ideal = 2.00;  % Y site can hold 2 cations (Ti, Fe3, Al, Cr, V, Zr)
                    X_ideal = 3.00;  % X site can hold 3 cations (Ca, Fe2+, Mg, etc.)

                    % Initialize variables for Z site (excluding Cr, V, Zr)
                    Si_Z = 0;
                    Fe3_Z = 0;
                    Ti_Z = 0;
                    Al_Z = 0;
                    
                    % Distribute Ti, Fe3, Si, and Al to the Z site first
                    if Si_n12 <= Z_ideal
                        Si_Z = Si_n12;   % All Si goes to Z site
                        remaining_Z = Z_ideal - Si_Z; % Remaining Z capacity after Si
                        Fe3_Z = min(Fe3, remaining_Z); % Fill remaining Z site with Fe3
                        remaining_Z = remaining_Z - Fe3_Z; % Update remaining capacity for Z
                    
                        % Fill remaining Z site with Ti, Al if space allows
                        Ti_Z = min(Ti_n12, remaining_Z);
                        remaining_Z = remaining_Z - Ti_Z;
                        Al_Z = min(Al_n12, remaining_Z);
                    else
                        % If there's not enough room for all Si, fill Z site as much as possible
                        Si_Z = Z_ideal;
                        remaining_Z = Z_ideal - Si_Z;
                        % Distribute Fe3, Ti, Al to Z site
                        Fe3_Z = min(Fe3, remaining_Z);
                        remaining_Z = remaining_Z - Fe3_Z;
                        Ti_Z = min(Ti_n12, remaining_Z);
                        remaining_Z = remaining_Z - Ti_Z;
                        Al_Z = min(Al_n12, remaining_Z);
                    end
                    
                    % Update remaining Ti, Fe3, Al after Z site distribution
                    Ti_remaining = Ti_n12 - Ti_Z;
                    Fe3_remaining = Fe3 - Fe3_Z;
                    Al_remaining = Al_n12 - Al_Z;
                    
                    % Initialize variables for Y site (Cr, V, Zr should be placed here)
                    Ti_Y = 0;
                    Fe3_Y = 0;
                    Al_Y = 0;
                    Cr_Y = 0;
                    V_Y = 0;
                    Zr_Y = 0;
                    
                    % Distribute remaining Ti, Fe3, Al, Cr, V, and Zr to the Y site
                    if Ti_remaining + Fe3_remaining + Al_remaining + Cr_n12 + V_n12 + Zr_n12 <= Y_ideal
                        Ti_Y = Ti_remaining;    % All remaining Ti fits in Y site
                        Fe3_Y = Fe3_remaining;  % All remaining Fe3 fits in Y site
                        Al_Y = Al_remaining;    % All remaining Al fits in Y site
                        Cr_Y = Cr_n12;          % All Cr fits in Y site
                        V_Y = V_n12;            % All V fits in Y site
                        Zr_Y = Zr_n12;          % All Zr fits in Y site
                    else
                        % Fill Y site step by step
                        Ti_Y = min(Ti_remaining, Y_ideal);    % Fill Y site with Ti
                        remaining_Y = Y_ideal - Ti_Y;
                        
                        Fe3_Y = min(Fe3_remaining, remaining_Y);  % Fill Y site with Fe3
                        remaining_Y = remaining_Y - Fe3_Y;
                    
                        Al_Y = min(Al_remaining, remaining_Y);    % Fill Y site with Al
                        remaining_Y = remaining_Y - Al_Y;
                    
                        Cr_Y = min(Cr_n12, remaining_Y);          % Fill Y site with Cr
                        remaining_Y = remaining_Y - Cr_Y;
                    
                        V_Y = min(V_n12, remaining_Y);            % Fill Y site with V
                        remaining_Y = remaining_Y - V_Y;
                    
                        Zr_Y = min(Zr_n12, remaining_Y);          % Fill Y site with Zr
                    end
                    
                    % Update remaining cations after Y site distribution
                    Fe3_excess = Fe3_remaining - Fe3_Y;
                    Al_excess = Al_remaining - Al_Y;
                    Cr_excess = Cr_n12 - Cr_Y;
                    V_excess = V_n12 - V_Y;
                    Zr_excess = Zr_n12 - Zr_Y;
                    
                    % Convert excess Fe3 to Fe2 for the X site
                    if Fe3_excess > 0
                        Fe2 = Fe3_excess;  % Excess Fe3 is reduced to Fe2
                        Fe3_excess = 0;          % Reset excess Fe3 after conversion
                    end

                    [sorted_formula1, total1] = Formula_Output({Fe2, Mg_n12, Ca_n12, Mn_n12, Na_n12, K_n12}, {'Fe2+','Mg', 'Ca', 'Mn', 'Na', 'K'});
                    [sorted_formula2, total2] = Formula_Output({Ti_Y, Fe3_Y, Al_Y, Cr_Y, V_Y, Zr_Y}, {'Ti','Fe3+','AlVI','Cr','V','Zr'});
                    [sorted_formula3, total3] = Formula_Output({Si_Z, Fe3_Z, Ti_Z, Al_Z}, {'Si','Fe3+','Ti','AlIV'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, 'O12')};
                    
                end

            end

        end
           
        if(~strcmp(T.group2(i),'Tectosilicate') && (Si_n12 + Al_IV) > 3-3*t && (Si_n12 + Al_IV) < 3+3*t)

            Prp =  Mg_n12 / (Mg_n12 + Fe2 + Ca_n12 + Mn_n12) * 100;
            Alm = Fe2 / (Mg_n12 + Fe2 + Ca_n12 + Mn_n12) * 100;
            Grs = Ca_n12 / (Mg_n12 + Fe2 + Ca_n12 + Mn_n12) * 100;
            Sps = Mn_n12 / (Mg_n12 + Fe2 + Ca_n12 + Mn_n12) * 100;
            T.Prp(i) = round(Prp,2); T.Alm(i) = round(Alm,2); T.Grs(i) = round(Grs,2); T.Sps(i) = round(Sps,2);

            if(ABCDT_n/24 >= (8/12)-(8/12)*t && ABCDT_n/24 <= (8/12)+(8/12)*1.25*t && (Ca_n12 + Mg_n12 + Mn_n12 + Fe2) / (Fe3 + Al_VI + Cr_n12 + Ti_n12 + Zr_n12 + V_n12) > (3/2)-(3/2)*t && (Ca_n12 + Mg_n12 + Mn_n12 + Fe2) / (Fe3 + Al_VI + Cr_n12 + Ti_n12 + Zr_n12 + V_n12) < (3/2)+(3/2)*2.2*t)
                if(ABCD_n/2 >= 5 - (1/2)*5*t && ABCD_n/2 <= 5+5*1.3*t && Na_n12 + K_n12 >= 0 && Na_n12 + K_n12 < 0.2+0.2*t && S_n <= 0.05 && P_n <= 0.05 && Al_IV + Si_n12 > 3 - 3*t && Al_IV + Si_n12 < 3+3*t && (Na_n12+K_n12) >= 0 && (Na_n12+K_n12) < 0.2+0.2*t)
                   
                    T.group2(i) = {'Nesosilicate'};
                    T.group3(i) = {'Garnet'};
                    T.group4(i) = {''};
                    T.species(i) = {''}; 

                    if((Mg_n + Fe_n + Al_n)/Si_n > 5/3 + 5/3*t && Mg_n > Fe_n ...
                        || (Si_n*(12/24) > 3.3 && Ca_n*(12/24) < 0.3) )
                        T.group2(i) = {'Phyllosilicate'};
                        T.group3(i) = {''}; T.group4(i) = {''}; T.species(i) = {''}; T.formula(i) = {''};
                    end

                    % [sorted_formula1, total1] = Formula_Output({Mg_n12, Fe2, Ca_n12, Mn_n12}, {'Mg', 'Fe2+', 'Ca', 'Mn'});
                    % [sorted_formula2, total2] = Formula_Output({Al_VI, Fe3, Ti_n12, Cr_n12, V_n12, Zr_n12}, {'Al', 'Fe3+', 'Ti', 'Cr', 'V', 'Zr'});
                    % [sorted_formula3, total3] = Formula_Output({Si_n12, Al_IV}, {'Si','AlIV'});
                    % T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, 'O12')};

                    Norm_Factor = (8/ (T.ABCDT(i)/2));

                    Mg_c8 = Mg_n12*Norm_Factor;
                    Mn_c8 = Mn_n12*Norm_Factor;
                    Fe_c8 = Fe_n12*Norm_Factor;
                    Ca_c8 = Ca_n12*Norm_Factor;
                    Si_c8 = Si_n12*Norm_Factor;
                    Na_c8 = Na_n12*Norm_Factor;
                    K_c8 = K_n12*Norm_Factor;
                    Cr_c8 = Cr_n12*Norm_Factor;
                    Ti_c8 = Ti_n12*Norm_Factor;
                    Al_c8 = Al_n12*Norm_Factor;
                    V_c8 = V_n12*Norm_Factor;
                    Zr_c8 = Zr_n12*Norm_Factor;

                    Al_IV = 3 - Si_c8;
                    if(Al_IV<0.001), Al_IV = 0; end
            
                    Al_VI = Al_c8 - Al_IV;
                    if(Al_VI<0), Al_VI = 0; end
            
                    Fe3 = 2 - (Al_VI + Cr_c8 + Ti_c8 + V_c8);
                    if(Fe3<0.001), Fe3 = 0; end
            
                    Fe2 = Fe_c8 - Fe3; 
                    if(Fe2<0), Fe2 = 0; end

                    [sorted_formula1, total1] = Formula_Output({Mg_c8, Fe2, Ca_c8, Mn_c8}, {'Mg', 'Fe2+', 'Ca', 'Mn'});
                    [sorted_formula2, total2] = Formula_Output({Al_VI, Fe3, Ti_c8, Cr_c8, V_c8, Zr_c8}, {'Al', 'Fe3+', 'Ti', 'Cr', 'V', 'Zr'});
                    [sorted_formula3, total3] = Formula_Output({Si_c8, Al_IV}, {'Si','AlIV'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, 'O12')};

                    %% 7.2.2.1. Pyralspite Garnets
                    if(Ca_n12 / (Ca_n12 + Mg_n12 + Fe2 + Mn_n12) < 0.5 && Al_VI/(Al_VI+Ti_n12+Cr_n12+Fe3) > 0.5 && ~strcmp(T.group2(i), 'Phyllosilicate'))
                        T.group4(i) = {'Pyralspite Garnet'}; % updated group
                        [~, maxi] = max([Prp Alm Grs Sps]);
    
                        if(maxi == 1)
                            T.species(i) = {'pyrope'};

                            if(round(Al_n/Si_n,2) < 2/3 && Fe_n > 0.1)
                                T.group2(i) = {'Phyllosilicate'};
                                T.group3(i) = {''}; T.group4(i) = {''}; T.species(i) = {''}; T.formula(i) = {''};
                            end

                        end
    
                        if(maxi == 2)
                            T.species(i) = {'almandine'};
                        end
    
                        if(maxi == 4)
                            T.species(i) = {'spessartine'};
                        end
                    end % end of pyralspite garnets

                    %% 7.2.2.2. Ugrandite garnets
                    if(Ca_n12 / (Ca_n12 + Mg_n12 + Fe2 + Mn_n12) > 0.5 && ~strcmp(T.group2(i), 'Phyllosilicate'))
    
                        T.group4(i) = {'Ugrandite Garnet'};
                        
                        if(Al_VI/(Al_VI+Fe3+Cr_n12+Ti_n12+V_n12+Zr_n12) > 0.5 && Ca_n)
                            T.species(i) = {'grossular'};

                            if(Mn_n > Fe_n)

                                Fe3 = Fe_n12;
                                Mn3 = 2 - (Al_VI + Fe3);
                                if(Al_VI<0), Al_VI = 0; end

                                Mn2_b = 1 - (Mg_n12);
                                if(Mn2_b<0), Mn2_b = 0; end

                                Mn2 = Mn_n12 - Mn3;
                                if(Mn2<0), Mn2 = 0; end

                                Mn2_a = Mn2 - Mn2_b;
                                if(Mn2_a<0), Mn2_a = 0; end

                                if((Mn2_b + Mg_n12)/(Al_VI + Mn3 + Fe3) > 1/2 - (1/2)*t && (Mn2_b + Mg_n12)/(Al_VI + Mn3 + Fe3) < 1/2 + (1/2)*t)
                                    if((Mn2_b+Mg_n12)/(Ca_n12+Mn2_a) > 1/2 - (1/2)*t && (Mn2_b+Mg_n12)/(Ca_n12+Mn2_a) < 1/2 + (1/2)*t)
                                        
                                        Al_IV = 3 - Si_n/2;
                                        if(Al_IV < 0), Al_IV = 0; end
                                        
                                        Al_VI = Al_n/2 - Al_IV;
                                        if(Al_VI<0), Al_VI = 0; end
                        
                                        Si_a = 1;
                                        Si_b = Si_n/2 - Si_a;
                                        if(Si_b < 0), Si_b = 0; end
                        
                                        Al_IV = 2 - Si_b;
                                        if(Al_IV<0), Al_IV = 0; end

                                        T.group2(i) = {'Sorosilicate'};
                                        T.group3(i) = {'Calcic Sorosilicate'};
                                        T.group4(i) = {'Epidote'};
                        
                                        [sorted_formula1, total1] = Formula_Output({Ca_n/2, Na_n/2}, {'Ca','Na'});
                                        [sorted_formula2, total2] = Formula_Output({Al_VI, Fe_n/2, Mg_n/2, Mn_n/2, Cr_n/2}, {'AlVI','Fe3+','Mg','Mn','Cr'});
                                        [sorted_formula3, total3] = Formula_Output({Si_b, Al_IV}, {'Si','AlIV'});
                                        [sorted_formula4, total4] = Formula_Output({Si_a}, {'Si'});
                                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, '[',sorted_formula3,total3,'O7][', sorted_formula4, total4,'O4]O(OH)')};
                                        T.species(i) = {'piemontite or pumpelleyite-Mn'};

                                    end
                                end


                            end

                            if((Fe2+Mg_n12+Mn_n12)/(Al_VI + Fe3 + Ti_n12) > 1/2 - (1/2)*t && (Fe2+Mg_n12+Mn_n12)/(Ca_n12+Na_n12) < 1/2 + (1/2)*1.25*t)

                                Al_IV = 3 - Si_n/2;
                                if(Al_IV < 0), Al_IV = 0; end
                                
                                Al_VI = Al_n/2 - Al_IV;
                                if(Al_VI<0), Al_VI = 0; end
                
                                Si_a = 1;
                                Si_b = Si_n/2 - Si_a;
                                if(Si_b < 0), Si_b = 0; end
                
                                Al_IV = 2 - Si_b;
                                if(Al_IV<0), Al_IV = 0; end

                                T.group2(i) = {'Sorosilicate'};
                                T.group3(i) = {'Calcic Sorosilicate'};
                                T.group4(i) = {'Epidote'};
                
                                [sorted_formula1, total1] = Formula_Output({Ca_n/2, Na_n/2}, {'Ca','Na'});
                                [sorted_formula2, total2] = Formula_Output({Al_VI, Fe_n/2, Mg_n/2, Mn_n/2, Cr_n/2}, {'AlVI','Fe3+','Mg','Mn','Cr'});
                                [sorted_formula3, total3] = Formula_Output({Si_b, Al_IV}, {'Si','AlIV'});
                                [sorted_formula4, total4] = Formula_Output({Si_a}, {'Si'});
                                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, '[',sorted_formula3,total3,'O7][', sorted_formula4, total4,'O4]O(OH)')};

                                if(Fe2 > Mn_n && Fe2 > Mg_n)
                                    T.species(i) = {'epidote or pumpelleyite-Fe'};
                                end

                                if(Mg_n > Mn_n && Mg_n > Fe2)
                                    T.species(i) = {'pumpelleyite-Mg'};
                                end

                            end

                        end
    
                        if(Fe3/(Al_VI+Fe3+Cr_n12+Ti_n12+V_n12+Zr_n12) > 0.5)
                            T.species(i) = {'andradite'};
                        end

                        if(Cr_n12/(Al_VI+Fe3+Cr_n12+Ti_n12+V_n12+Zr_n12) > 0.5)
                            T.species(i) = {'uvarovite'};
                        end

                    end %end of ugrandite garnets

                    if(Mg_n12 / (Ca_n12 + Mg_n12 + Fe2 + Mn_n12) > 0.5 && ~strcmp(T.group2(i), 'Phyllosilicate'))
                        if(Cr_n12/(Al_VI + Fe3 + Ti_n12 + Cr_n12 + V_n12 + Zr_n12) > 0.5)
                            T.group4(i) = {'Other Garnet'};
                            T.species(i) = {'knorringite'};
                        end

                    end

                    if(Ca_n12 / (Ca_n12 + Mg_n12 + Fe2 + Mn_n12) > 0.5 && ~strcmp(T.group2(i), 'Phyllosilicate'))
                        if(V_n12/(Al_VI + Fe3 + Ti_n12 + Cr_n12 + V_n12 + Zr_n12) > 0.5)
                            T.group4(i) = {'Other Garnet'};
                            T.species(i) = {'goldmanite'};
                        end

                        if(round(Ti_n12/(Al_VI + Fe3 + Ti_n12 + Cr_n12 + V_n12 + Zr_n12),2) >= 0.5 && Fe3 > Al_VI && Fe3 > Cr_n12 && Fe3 > V_n12)
                            T.group4(i) = {'Other Garnet'};
                            T.species(i) = {'morimotoite'};
                        end
                    end

                    if(Mn_n12 / (Ca_n12 + Mg_n12 + Fe2 + Mn_n12) > 0.5 && ~strcmp(T.group2(i), 'Phyllosilicate'))
                        if(Fe3/(Al_VI + Fe3 + Ti_n12 + Cr_n12 + V_n12 + Zr_n12) > 0.5)
                            T.group4(i) = {'Other Garnet'};
                            T.species(i) = {'calderite'};
                        end

                        if(V_n12/(Al_VI + Fe3 + Ti_n12 + Cr_n12 + V_n12 + Zr_n12) > 0.5)
                            T.group4(i) = {'Other Garnet'};
                            T.species(i) = {'momoiite'};
                        end
                    end

                    if(Fe2 / (Ca_n12 + Mg_n12 + Fe2 + Mn_n12) > 0.5 && ~strcmp(T.group2(i), 'Phyllosilicate'))
                        if(Fe3/(Al_VI + Fe3 + Ti_n12 + Cr_n12 + V_n12 + Zr_n12) > 0.5)
                            T.group2(i) = {'Nesosilicate'};
                            T.group3(i) = {'Garnet'};
                            T.group4(i) = {'Other Garnet'};
                            T.species(i) = {'skiagite'};
                        end
                    end

                end  % End of garnets
            end
        end 

        %% 7.31 Sorosilicate
        if(strcmp(T.group1(i), 'Silicate'))
                
            Renorm_factor = 24 / 12.5;
            Mg_n12 = T.Mg_Multi(i) / Renorm_factor;
            Mn_n12 = T.Mn_Multi(i) / Renorm_factor;
            Fe_n12 = T.Fe_Multi(i) / Renorm_factor;
            Ca_n12 = T.Ca_Multi(i) / Renorm_factor;
            Na_n12 = T.Na_Multi(i) / Renorm_factor;
            K_n12 = T.K_Multi(i) / Renorm_factor;
            Si_n12 = T.Si_Multi(i) / Renorm_factor;
            Cr_n12 = T.Cr_Multi(i) / Renorm_factor;
            Ti_n12 = T.Ti_Multi(i) / Renorm_factor;
            Al_n12 = T.Al_Multi(i) / Renorm_factor;
            V_n12 = T.V_Multi(i) / Renorm_factor;
            Zr_n12 = T.Zr_Multi(i) / Renorm_factor;
            Sr_n12 = T.Sr_Multi(i) / Renorm_factor;

            Al_IV = 3 - Si_n12;
            if(Al_IV < 0), Al_IV = 0; end
            
            Al_VI = Al_n12 - Al_IV;
            if(Al_VI<0), Al_VI = 0; end
            
            if((Al_VI + Fe_n12 + Mn_n12)/(Si_n12+Al_IV) > 1 - 1*t && (Al_VI + Fe_n12 + Mn_n12)/(Si_n12+Al_IV) < 1 + 1*t)
                if((Ca_n12+Sr_n12)/(Al_VI+Fe_n12+Mn_n12) > 2/3 - (2/3)*1.5*t && (Ca_n12+Sr_n12)/(Al_VI+Fe_n12+Mn_n12) < 2/3 + (2/3)*1.5*t && (Fe_n12+Mg_n12+Mn_n12)/Al_VI >= 0 && (Fe_n12+Mg_n12+Mn_n12)/Al_VI < (1/2) + (1/2)*2*t)
                    if(Al_VI > Fe_n12 && Al_VI > Mg_n12 && Al_VI > Mn_n12 && (Ca_n12+Sr_n12) > Na_n12 && (Ca_n12+Sr_n12) > K_n12)

                        T.group2(i) = {'Sorosilicate'};
                        T.group3(i) = {'Calcic Sorosilicate'};
                        T.group4(i) = {'Epidote'};

                        Si_a = 1;
                        Si_b = Si_n12 - Si_a;
                        if(Si_b < 0), Si_b = 0; end

                        Al_IV = 2 - Si_b;
                        if(Al_IV<0), Al_IV = 0; end

                        [sorted_formula1, total1] = Formula_Output({Ca_n12, Na_n12}, {'Ca','Na'});
                        [sorted_formula2, total2] = Formula_Output({Al_VI, Fe_n12, Mg_n12, Mn_n12, Cr_n12}, {'AlVI','Fe3+','Mg','Mn','Cr'});
                        [sorted_formula3, total3] = Formula_Output({Si_b, Al_IV}, {'Si','AlIV'});
                        [sorted_formula4, total4] = Formula_Output({Si_a}, {'Si'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, '[',sorted_formula3,total3,'O7][', sorted_formula4, total4,'O4]O(OH)')};

                        if((Fe_n12 + Mg_n12 + Mn_n12)/Al_VI >= 0 && (Fe_n12 + Mg_n12 + Mn_n12)/Al_VI < 1/2)
                            if(Na_n/Ca_n < 0.1 )
                                T.species(i) = {'zoisite/clinozoisite or pumpelleyite-Al'};
                            end
    
                                if((Sr_n12+Ca_n12)/(Al_VI+Fe_n12+Ti_n12) >= 2/3 - (2/3)*t && (Sr_n12+Ca_n12)/(Al_VI+Fe_n12+Ti_n12) < 2/3 + (2/3)*t)
                                    if(Sr_n12/Ca_n12 >= 0.5)
                                        T.species(i) = {'niigataite'};
                                    end
                                end
    
                                if((2 - Si_n/4)/ (Si_n/4) >= 1 - 1*t)
                                    T.group1(i) = {'Silicate'};
                                    T.group2(i) = {'Inosilicate'};
                                    T.group3(i) = {'Pyroxene'};
                                    T.group4(i) = {'Clinopyroxne'};
                                    T.species(i) = {'Ca-Tschermaks/kushiroite'};
                                end
                        end

                        if(round((Fe_n12 + Mg_n12 + Mn_n12)/Al_VI,2) >= 0.4 - 0.4*t)
                            if(Fe_n12 > Mg_n12 && Fe_n12 > Mn_n12)
                                T.species(i) = {'epidote or pumpelleyite-Fe'};
                            end
                            if(Mg_n12 > Fe_n12 && Mg_n12 > Mn_n12)
                                T.species(i) = {'pumpelleyite-Mg'};
                            end
                            if(Mn_n12 > Fe_n12 && Mn_n12 > Mg_n12)
                                T.species(i) = {'piemontite or pumpelleyite-Mn'};
                            end
                        end

                    end 
                end
            end 

            if((Al_VI + Fe_n12 + Mn_n12)/(Si_n12+Al_IV) > 1 - 1*t && (Al_VI + Fe_n12 + Mn_n12)/(Si_n12+Al_IV) < 1 + 1*t)
                if((Ca_n12+Sr_n12)/(Al_VI+Fe_n12+Mn_n12) > 2/3 - (2/3)*t && (Ca_n12+Sr_n12)/(Al_VI+Fe_n12+Mn_n12) < 2/3 + (2/3)*t && (Fe_n12+Mg_n12+Mn_n12)/Al_VI > 2 - 2*t && (Fe_n12+Mg_n12+Mn_n12)/Al_VI <  2 + 2*t)
                    if(Al_VI > Fe_n12 && Al_VI > Mg_n12 && Al_VI < Mn_n12 && (Ca_n12+Sr_n12) > Na_n12 && (Ca_n12+Sr_n12) > K_n12)
    
                        T.group2(i) = {'Sorosilicate'};
                        T.group3(i) = {'Calcic Sorosilicate'};
                        T.group4(i) = {'Epidote'};
                       
                        if(Sr_n12/Ca_n12 >= 1 - 2*1*t && Sr_n12/Ca_n12 < 1 + 2*1*t && Mn_n12 > Fe_n12)
                            if(Al_VI/(Mn_n12+Fe_n12) >= 1/2 - (1/2)*2*t && Al_VI/(Mn_n12+Fe_n12) < 1/2 + (1/2)*2*t)
                                T.species(i) = {'tweddillite'};

                                Si_a = 1;
                                Si_b = Si_n12 - Si_a;
                                if(Si_b < 0), Si_b = 0; end
        
                                Al_IV = 2 - Si_b;
                                if(Al_IV<0), Al_IV = 0; end
        
                                [sorted_formula1, total1] = Formula_Output({Ca_n12, Na_n12}, {'Ca','Na'});
                                [sorted_formula2, total2] = Formula_Output({Al_VI, Fe_n12, Mg_n12, Mn_n12, Cr_n12}, {'AlVI','Fe3+','Mg','Mn','Cr'});
                                [sorted_formula3, total3] = Formula_Output({Si_b, Al_IV}, {'Si','AlIV'});
                                [sorted_formula4, total4] = Formula_Output({Si_a}, {'Si'});
                                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, '[',sorted_formula3,total3,'O7][', sorted_formula4, total4,'O4]O(OH)')};
                            end
                        end
    
                    end
                end
             end % end of sorosilicate

        end 

        %% 7.5 Al-Silicate
        if(strcmp(T.group2(i), 'Al-Silicate') && S_n <= 0.05 && P_n <= 0.05) %updated group
            Renorm_factor = 24 / 5;
            Fe_n5 = T.Fe_Multi(i) / Renorm_factor;
            Cr_n5 = T.Cr_Multi(i) / Renorm_factor;
            Si_n5 = T.Si_Multi(i) / Renorm_factor;
            Al_n5 = T.Al_Multi(i) / Renorm_factor;
            Ti_n5 = T.Ti_Multi(i) / Renorm_factor;
            Mg_n5 = T.Mg_Multi(i) / Renorm_factor;
            Fe3 = Fe_n5;
            Al_IV = 1 - Si_n5;
            if(Al_IV < 0), Al_IV = 0; end
            Al_VI = Al_n5 - Al_IV;
            if(Al_VI < 0), Al_VI = 0; end
            T.group2(i) = {'Nesosilicate'};
%             T.group3(i) = {'Nesosilicate'};
            T.group3(i) = {'Al-Silicate'};
            T.species(i) = {'kyanite/sillimanite/andalusite'};
            [sorted_formula1, total1] = Formula_Output({Al_VI, Fe3, Cr_n5, Ti_n5, Mg_n5},{'AlIV','Fe3+','Cr','Ti','Mg'});
            [sorted_formula2, total2] = Formula_Output({Si_n5,Al_IV}, {'Si','AlIV'});
            T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O5')};
        end

        %% 7.6 Zircon
        if(strcmp(T.group3(i), 'Zircon'))
            T.species(i) = {'zircon'};
            renorm = 4/24;
            [sorted_formula1, total1] = Formula_Output({Zr_n*renorm, Hf_n*renorm},{'Zr','Hf'});
            T.formula(i) = {strcat(sorted_formula1, total1, '(Si_', num2str(round(Si_n*renorm, 2)), 'O4)')};
        end

        %% 7.4.1 K-Na Mica Catch
        if(strcmp(T.group1(i), 'Silicate') && ~strcmp(T.group2(i), 'Tectosilicate') && S_n <= 0.05 && P_n <= 0.05)
            if((T.K2O(i) > 8 && T.K2O(i) < 13 && Al_n/Si_n <= 1 && K_n/(K_n+Na_n+Ca_n) > 0.5)...
                    || (T.Na2O(i) > 6*(100/95) && Al_n/Si_n <= 1 && Na_n/(K_n+Na_n+Ca_n) > 0.5)...
                    || (T.Al2O3(i) > 50 && Ca_n/(Si_n+Al_n+Fe_n+Mg_n) <= 1/6 + (1/6)*t && Ca_n/(K_n+Na_n+Ca_n) > 0.5))

                if(Na_n/Si_n <= 1/3)
                    T.group2(i) = {'Phyllosilicate'};
                    T.group3(i) = {'Mica'};
                    T.group4(i) = {'K-Na Mica'};
                    T.species(i) = {''};
                    T.formula(i) = {''};
                end

            end
        end

        % Catch for amphiboles incorrectly in phyllosilicates
        if(strcmp(T.group4(i), 'K-Na Mica') || strcmp(T.group3(i), 'Mica') && S_n <= 0.05 && P_n <= 0.05 )
            Renorm_factor = 24 / 23;
            Mg_n23 = T.Mg_Multi(i) / Renorm_factor;
            Mn_n23 = T.Mn_Multi(i) / Renorm_factor;
            Na_n23 = T.Na_Multi(i) / Renorm_factor;
            Fe_n23 = T.Fe_Multi(i) / Renorm_factor;
            Ca_n23 = T.Ca_Multi(i) / Renorm_factor;
            Si_n23 = T.Si_Multi(i) / Renorm_factor;
            Ti_n23 = T.Ti_Multi(i) / Renorm_factor;
            Cr_n23 = T.Cr_Multi(i) / Renorm_factor;
            Al_n23 = T.Al_Multi(i) / Renorm_factor;
            K_n23 = T.K_Multi(i) / Renorm_factor;
            F_n23 = T.F_Multi(i) / Renorm_factor;
            Cl_n23 = T.Cl_Multi(i) / Renorm_factor;
            Al_IV = 8 - Si_n23;
            if(Al_IV < 0), Al_IV = 0; end
            Al_VI = Al_n23 - Al_IV;
            if(Al_VI < 0), Al_VI = 0; end

            if((Na_n23+Ca_n23+K_n23)/(Fe_n23+Mg_n23+Mn_n23+Al_VI+Cr_n23+Ti_n23) > (3/5) - (3/5)*2*t && (Na_n23+Ca_n23+K_n23)/(Fe_n23+Mg_n23+Mn_n23+Al_VI+Cr_n23+Ti_n23) < 3/5 + (3/5)*t && Na_n23/(Na_n23+Ca_n23+K_n23) > 0.5)
                if(Al_IV/Si_n23 < 1/8)
                    if(Fe_n23/(Fe_n23+Mg_n23+Mn_n23+Al_VI+Cr_n23+Ti_n23) > 0.5)
                        T.group2(i) = {'Inosilicate'};
                        T.group3(i) = {'Amphibole'};
                        T.group4(i) = {'Sodic-Calcic Amphibole'};
                        T.species(i) = {''};
%                         T.species(i) = {'ferro-richterite'};
                        Nax = 1 - K_n23;
                        if(Nax <= 0), Nax = 0; end
                        Nay = Na_n23 - Nax;
                        if(Nay <= 0), Nay = 0; end
                        vac = 2 - (Nay + Ca_n23);
                        if(vac <= 0), vac = 0; end
                        [sorted_formula1, total1] = Formula_Output({Nax, K_n23}, {'Na','K'});
                        [sorted_formula2, total2] = Formula_Output({Nay, Ca_n23, vac}, {'Na', 'Ca','□'});
                        [sorted_formula3, total3] = Formula_Output({Si_n23, Al_IV}, {'Si','AlIV'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, 'Fe2+_', num2str(round(Fe_n23,2)), sorted_formula3,total3,'O22(OH)2')};
                    end

                    if(Mg_n23/(Fe_n23+Mg_n23+Mn_n23+Al_VI+Cr_n23+Ti_n23) > 0.5)
                        T.group2(i) = {'Inosilicate'};
                        T.group3(i) = {'Amphibole'};
                        T.group4(i) = {'Sodic-Calcic Amphibole'};
                        T.species(i) = {''};
%                         T.species(i) = {'richterite'};
                        Nax = 1 - (K_n23+Na_n23);
                        if(Nax <= 0), Nax = 0; end
                        Nay = Na_n23 - Nax;
                        if(Nay <= 0), Nay = 0; end
                        vac = 2 - (Nax + Ca_n23);
                        if(vac <= 0), vac = 0; end
                        OH = 2 - (F_n23+Cl_n23);
                        if(OH <= 0), OH = 0; end
                        [sorted_formula1, total1] = Formula_Output({Nax, K_n23, vac}, {'Na','K','□'});
                        [sorted_formula2, total2] = Formula_Output({Nay, Ca_n23, Mn_n23}, {'Na', 'Ca','Mn'});
                        [sorted_formula3, total3] = Formula_Output({Mg_n23, Fe_n23, Al_VI, Ti_n23, Cr_n23}, {'Mg','Fe','AlVI','Ti','Cr'});
                        [sorted_formula4, total4] = Formula_Output({Si_n23, Al_IV}, {'Si','AlIV'});
                        [sorted_formula5, total5] = Formula_Output({OH, F_n23, Cl_n23}, {'OH','F','Cl'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3, sorted_formula4, total4, 'O22', sorted_formula5, total5)};
                    end

                end   
            end

            if(Al_VI / (Na_n23+Ca_n23+K_n23+Mg_n23+Mn_n23+Fe_n23+Cr_n23+Ti_n23) > 2/5 - 2/5*t && Al_n / (Na_n+Mg_n+Fe_n) < 2/5 + 2/5*t && Mg_n23 >= Fe_n23 && Ca_n23/Na_n23 < 0.3)
                if(Al_IV/Si_n23 < 1/8 && Na_n23/(Na_n23+Ca_n23+K_n23) > 0.5)
                    T.group2(i) = {'Inosilicate'};
                    T.group3(i) = {'Amphibole'};
                    T.group4(i) = {'Non-Calcic Amphibole'};
                    T.species(i) = {''};
%                     T.species(i) = {'glaucophane'};
                    Fe3 = 2 -(Al_VI + Ti_n23 + Cr_n23);
                    if(Fe3 <= 0), Fe3 = 0; end
                    Fe2 = Fe_n23 - Fe3;
                    if(Fe2 <= 0), Fe2 = 0; end
                    [sorted_formula1, total1] = Formula_Output({Na_n23, Ca_n23, K_n23}, {'Na','Ca','K'});
                    [sorted_formula2, total2] = Formula_Output({Mg_n23, Fe2, Mn_n23}, {'Mg', 'Fe2+','Mn'});
                    [sorted_formula3, total3] = Formula_Output({Al_VI, Fe3, Ti_n23, Cr_n23}, {'AlVI','Fe3+','Ti','Cr'});
                    [sorted_formula4, total4] = Formula_Output({Si_n23, Al_IV}, {'Si','AlIV'});
                    [sorted_formula5, total5] = Formula_Output({2 - (F_n23+Cl_n23), F_n23, Cl_n23}, {'OH','F', 'Cl'});
                    T.formula(i) = {strcat('◻', sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3, sorted_formula4, total4, 'O22', sorted_formula5, total5)};
                end
            end

            if((Na_n23+Ca_n23+K_n23) / (Fe_n23+Mg_n23+Mn_n23+Al_VI+Cr_n23+Ti_n23) >= 2/5 - 2/5*t && (Na_n23+Ca_n23+K_n23) / (Fe_n23+Mg_n23+Mn_n23+Al_VI+Cr_n23+Ti_n23) <= 2/5 + 2/5*t)
                if(Fe_n23/(Fe_n23+Mg_n23+Mn_n23+Al_VI+Cr_n23+Ti_n23) > 0.5 && Na_n23/(Na_n23+Ca_n23+K_n23) > 0.5 && Al_IV/Si_n23 < 1/8)
                    T.group2(i) = {'Inosilicate'};
                    T.group3(i) = {'Amphibole'};
                    T.group4(i) = {'Non-Calcic Amphibole'};
                    T.species(i) = {''};
%                     T.species(i) = {'riebeckite'};
                    Nax = 1 - (K_n23+Na_n23);
                    if(Nax <= 0), Nax = 0; end
                    Nay = Na_n23 - Nax;
                    if(Nay <= 0), Nay = 0; end
                    Fe2 = 3 - Mg_n23;
                    if(Fe2 <= 0), Fe2 = 0; end
                    Fe3 = Fe_n23 - Fe2;
                    if(Fe3 <= 0), Fe3 = 0; end
                    vac = 2 - (Nax + Ca_n23);
                    if(vac <= 0), vac = 0; end
                    OH = 2 - (F_n23+Cl_n23);
                    if(OH <= 0), OH = 0; end
                    [sorted_formula1, total1] = Formula_Output({Nax, K_n23, vac}, {'Na','K','□'});
                    [sorted_formula2, total2] = Formula_Output({Nay, Ca_n23, Mn_n23}, {'Na', 'Ca','Mn'});
                    [sorted_formula3, total3] = Formula_Output({Fe2, Mg_n23}, {'Fe2+','Mg'});
                    [sorted_formula4, total4] = Formula_Output({Fe_n23, Al_VI, Ti_n23, Cr_n23}, {'Fe3+','AlVI','Ti','Cr'});
                    [sorted_formula5, total5] = Formula_Output({Si_n23, Al_IV}, {'Si','AlIV'});
                    [sorted_formula6, total6] = Formula_Output({OH, F_n23, Cl_n23}, {'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3, sorted_formula4, total4, sorted_formula5, total5, 'O22', sorted_formula6, total6)};
                end
            end

        end

        %% 7.7 Chain Silicates - Amph-Px Assignment
        if(~strcmp(T.group3(i), 'Garnet') && ~strcmp(T.group4(i), 'K-Na Mica') && ~strcmp(T.group4(i), 'Epidote') && ABCDT_n/24 > 15/23 - 15/23*t && ABCDT_n/24 < 16/23 + 16/23*t && S_n <= 0.05 && P_n <= 0.05)
                
                Al_IV = 8-Si_n;
                if(Al_IV < 0), Al_IV = 0; end
                Al_VI = Al_n - Al_IV;
                if(Al_VI < 0), Al_VI = 0; end

            if(Si_n/24 > 6/23 - 6/23*t && Si_n/24 < 8/23 + 8/23*t && (T.K2O(i) < 6 || T.Na2O(i) < 6))
                if(((Si_n + Al_n)/24 > 8/23 - 8/23*t && (Si_n + Al_n)/24 < 10/23 + 10/23*t)...
                        || ((Na_n + K_n + Ca_n) > 4 - 4*t && (Na_n + K_n + Ca_n) < 4 + 4*t && Na_n/(Na_n + K_n + Ca_n) > 0.5))
                    if( (Al_n/Si_n >= 0 && Al_n/Si_n <= 1/2 + 1/2*t) || (Al_IV/Si_n > 1/7 - (1/7)*t  && Al_IV/Si_n < 1/3 + (1/3)*t) )
                        T.group2(i) = {'Inosilicate'};
                        T.group3(i) = {''}; T.group4(i) = {''};
                        T.species(i) = {''}; T.formula(i) = {''};
                        
                        ABC = Na_n + K_n + Ca_n + Mn_n + Mg_n + Fe_n + Al_VI + Ti_n + Cr_n;

                        if(ABCDT_n / 24 > 16/24 - 16/24*t && ABCDT_n / 24 < 16/24 + 16/24*t)
                            if((Si_n+Al_IV) / 24 > 8/24 - 2*(8/24*t) && (Si_n+Al_IV) / 24 < 8/24 + 8/24*t)
                                if(T.CaO(i) < 7.25 || T.CaO(i) > 15)
                                    if(T.K2O(i) < 0.2)

                                        if(Al_IV < 1.2 && round(ABC,3) >= 7.85)
                                            T.group3(i) = {'Pyroxene'};
                                        end

                                        if(strcmp(T.group3(i), 'Pyroxene'))
                                            if(Al_n < 0.1 && Na_n*(23/24) > 1.5 - 1.5*t && Na_n*(23/24) < 3 + 3*t)
                                                T.group3(i) = {'Amphibole'};
                                                
                                               if(Mg_n/Si_n > 1.5/4 && Mg_n/Si_n < 2.8/4 && Na_n/Si_n > 0.1 && Na_n/Si_n < 0.4)
                                                    if(Mg_n > Fe_n && Na_n/(Na_n+Ca_n+K_n) > 0.5 && Al_n/Si_n < 0.2 && Na_n > Ca_n && Na_n > K_n)
                                                        T.group2(i) = {'Phyllosilicate'};
                                                        T.group3(i) = {'Smectite'};
                                                        T.group4(i) = {'Trioctahedral Smectite'};
                                                        T.species(i) = {''};
                                                    end
                                                end

                                                if(Na_n+Ca_n/Si_n > 0.01 && Na_n+Ca_n/Si_n < 0.8 && Fe_n/(Si_n+Al_n) > 0.5-0.5*t && Fe_n/(Si_n+Al_n) < 0.5+0.5*t)
                                                    if(Al_n/Si_n >= 0.5/3.5 && Al_n/Si_n <= 0.5)
                                                        T.group2(i) = {'Phyllosilicate'};
                                                        T.group3(i) = {'Smectite'};
                                                        T.group4(i) = {'Trioctahedral Smectite'};
                                                        T.species(i) = {''};
                                                    end
                                                end
                                            end
                                        end

                                        if(round(ABC,3) > 7 && round(ABC,3) < 7.85 && (Si_n+Al_IV) / 24 > (8/23)*(24/23) - ((8/23)*(24/23))*t && (Si_n+Al_IV) / 24 < ((8/23)*(24/23)) + (((8/23)*(24/23))*t))
                                            T.group3(i) = {'Amphibole'};

                                            if(Mg_n/Si_n > 1.5/4 && Mg_n/Si_n < 2.8/4 && Na_n/Si_n > 0.1 && Na_n/Si_n < 0.4)
                                                if(Mg_n > Fe_n && Na_n/(Na_n+Ca_n+K_n) > 0.5 && Al_n/Si_n < 0.2 && Na_n>Ca_n && Na_n>K_n)
                                                    T.group2(i) = {'Phyllosilicate'};
                                                    T.group3(i) = {'Smectite'};
                                                    T.group4(i) = {'Trioctahedral Smectite'};
                                                    T.species(i) = {''};
                                                end
                                            end

                                            if(Na_n+Ca_n/Si_n > 0.01 && Na_n+Ca_n/Si_n < 0.8 && Fe_n/(Si_n+Al_n) > 0.5-0.5*0.1 && Fe_n/(Si_n+Al_n) < 0.5+0.5*0.1)
                                                if(Al_n/Si_n >= 0.5/3.5 && Al_n/Si_n <= 0.5)
                                                    T.group2(i) = {'Phyllosilicate'};
                                                    T.group3(i) = {'Smectite'};
                                                    T.group4(i) = {''};
                                                    T.species(i) = {''};
                                                end
                                            end

                                        end

                                    end
                                end
                            end
                        end

                    end
                end
            end
        end %end of Amphibole or Pyroxene assignment

        %% 7.7.1 Amphiboles
        if(strcmp(T.group2(i), 'Inosilicate') && ~strcmp(T.group3(i), 'Pyroxene') && S_n <= 0.05 && P_n <= 0.05)
            Renorm_factor = 24 / 23;
            Mg_n23 = T.Mg_Multi(i) / Renorm_factor;
            Mn_n23 = T.Mn_Multi(i) / Renorm_factor;
            Na_n23 = T.Na_Multi(i) / Renorm_factor;
            Fe_n23 = T.Fe_Multi(i) / Renorm_factor;
            Ca_n23 = T.Ca_Multi(i) / Renorm_factor;
            Si_n23 = T.Si_Multi(i) / Renorm_factor;
            Ni_n23 = T.Ni_Multi(i) / Renorm_factor;
            Ti_n23 = T.Ti_Multi(i) / Renorm_factor;
            Cr_n23 = T.Cr_Multi(i) / Renorm_factor;
            Al_n23 = T.Al_Multi(i) / Renorm_factor;
            K_n23 = T.K_Multi(i) / Renorm_factor;
            Cl_n23 = T.Cl_Multi(i) / Renorm_factor;
            F_n23 = T.F_Multi(i) / Renorm_factor;
            ABCDT_n23 = T.ABCDT(i) / Renorm_factor;
            Al_IV = 8 - Si_n23;
            if(Al_IV < 0), Al_IV = 0; end
            Al_VI = Al_n23 - Al_IV;
            if(Al_VI < 0), Al_VI = 0; end
            A = Na_n23 + K_n23;
            B = Ca_n23 + Mn_n23;

            if(ABCDT_n23 > 15 - 15*t && ABCDT_n23 < 16 + 16*t)
                if(Si_n23 > 6 - 6*t && Si_n23 < 8 + 8*t)
                    if( ((A+B) > 2 - 2*t && (A+B) < 3 + 3*t && T.Na2O(i) + T.K2O(i) < 6)...
                            || ((Mg_n23+Fe_n23) > 5 && (Mg_n23+Fe_n23) < 7 && (A+B) >= 0 && (A+B) < 1.5) )
                        if(Al_n/Si_n >= 0 && Al_n/Si_n <= 2/6 + 2/6*t && S_n <= 0.05 && P_n <= 0.05)
                            T.group3(i) = {'Amphibole'};

                            if(Mg_n/Si_n > 1.5/4 && Mg_n/Si_n < 2.8/4 && Na_n/Si_n > 0.01 && Na_n/Si_n < 0.4)
                                if(Mg_n > Fe_n && Na_n/(Na_n+Ca_n+K_n) > 0.5 && Al_n/Si_n < 0.2 && Na_n>Ca_n && Na_n>K_n)
                                    T.group2(i) = {'Phyllosilicate'};
                                    T.group3(i) = {'Smectite'};
                                    T.group4(i) = {'Trioctahedral Smectite'};
                                    T.species(i) = {''};
                                end
                            end

                            if(Na_n+Ca_n/Si_n > 0.01 && Na_n+Ca_n/Si_n < 0.8 && Fe_n/(Si_n+Al_n) > 0.5-0.5*0.1 && Fe_n/(Si_n+Al_n) < 0.5+0.5*0.1)
                                if(Al_n/Si_n >= 0.5/3.5 && Al_n/Si_n <= 0.5)
                                    T.group2(i) = {'Phyllosilicate'};
                                    T.group3(i) = {'Smectite'};
                                    T.group4(i) = {''};
                                    T.species(i) = {''};                          
                                end
                            end

                        end
                    end
                end
            end

            %% 7.7.1.1. Calcic amphiboles
            if(T.K2O(i) < 5 && T.CaO(i) > 7 && T.CaO(i) < 13+13*t)
                if(Si_n23 + Al_IV > 8-8*t && Si_n23 + Al_IV < 8+8*t)
                    if(round(Mg_n23 + Fe_n23 + Al_VI + Ti_n23 + Cr_n23,1) >= 5-5*t && Mg_n23 + Fe_n23 + Al_VI + Ti_n23 + Cr_n23 < 5+5*t)
                        if(Ca_n23 + Mn_n23 >= 1 && Ca_n23 + Mn_n23 < 2.35)
                            if(Ca_n > Mn_n && S_n <= 0.05 && P_n <= 0.05)
                                T.group2(i) = {'Inosilicate'};
                                T.group3(i) = {'Amphibole'};
                                T.group4(i) = {'Calcic Amphibole'};
                                T.species(i) = {''};
                                T.formula(i) = {''};

                                if(Si_n23 > 7.6 && Mg_n23 > Fe_n23)
                                    if(Mg_n23 + Fe_n23 > 5-5*t && Mg_n23 + Fe_n23 < 5+5*t)
                                        if(Mg_n23 > 4.5-4.5*t && Mg_n23 < 5+5*t)
                                            if(Fe_n23 >= 0 && Fe_n23 < 0.5+0.5*t)
                                                Fe2 = Fe_n23;
                                                T.species(i) = {'tremolite'};
                                                OH = 2 - (Cl_n23+F_n23);
                                                if(OH<0), OH=0; end
                                                vac = 1 - (Na_n23 + K_n23);
                                                if(vac < 0), vac = 0; end
                                                [sorted_formula1, total1] = Formula_Output({Na_n23, K_n23, vac}, {'Na', 'K','□'});
                                                [sorted_formula2, total2] = Formula_Output({Ca_n23, Mn_n23}, {'Ca', 'Mn'});
                                                [sorted_formula3, total3] = Formula_Output({Mg_n23, Fe2, Al_VI}, {'Mg', 'Fe','AlVI'});
                                                [sorted_formula4, total4] = Formula_Output({Si_n23,Al_IV}, {'Si', 'AlIV'});
                                                [sorted_formula5, total5] = Formula_Output({OH, Cl_n23, F_n23}, {'OH','Cl','F'});
                                                T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, sorted_formula4,total4,'O22', sorted_formula5, total5)};
                                            end
                                        end

                                        if(Mg_n23 > 2.5-2.5*t && Mg_n23 < 4.5+4.5*t)
                                            if(Fe_n23 > 0.5-0.5*t && Fe_n23 < 2.5+2.5*t && Al_IV/Si_n23 < 1/8)
                                                if(Ca_n23+Na_n23 > 2-2*t && Ca_n23+Na_n23 < 2+2*(1/2)*t)
                                                    if(Mg_n23 + Fe_n23 + Al_VI + Mn_n23 + Ti_n23 + Cr_n23 > 5 - 5*t && Mg_n23 + Fe_n23 + Al_VI + Mn_n23 + Ti_n23 + Cr_n23 < 5 + 5*t)
                                                        Fe2 = Fe_n23;
                                                        OH = 2 - (Cl_n23+F_n23);
                                                        if(OH<0), OH=0; end
                                                        vac = 1 - (Na_n23 + K_n23);
                                                        if(vac < 0), vac = 0; end
                                                        [sorted_formula1, total1] = Formula_Output({Na_n23, K_n23, vac}, {'Na', 'K','□'});
                                                        [sorted_formula2, total2] = Formula_Output({Ca_n23, Mn_n23}, {'Ca', 'Mn'});
                                                        [sorted_formula3, total3] = Formula_Output({Mg_n23, Fe2, Ni_n23, Al_VI,Ti_n23,Cr_n23}, {'Mg', 'Fe','Ni','AlVI','Ti','Cr'});
                                                        [sorted_formula4, total4] = Formula_Output({Si_n23,Al_IV}, {'Si', 'AlIV'});
                                                        [sorted_formula5, total5] = Formula_Output({OH, Cl_n23, F_n23}, {'OH','Cl','F'});
                                                        T.species(i) = {'actinolite'};
                                                        T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, '[', sorted_formula4,total4,'O22]', sorted_formula5, total5)};
                                                    end
                                                end
                                            end
                                        end

                                    end
                                end

                                 if(Al_IV / Si_n23 >= 0 && Al_IV / Si_n23 < 0.5+0.5*t)
                                    if(Al_IV >= 0 && Al_IV < 0.5)
                                        if((Mg_n23 + Fe_n23+Al_VI) > 5 - 5*t && (Mg_n23 + Fe_n23+Al_VI) < 5 + 5*t)
                                            if(Mg_n > Fe_n)
                                                if(Na_n23 + Ca_n23 +K_n23 > 2 - 2*t && Na_n23 + Ca_n23 +K_n23 < 2 + 2*2*t && Ca_n23 + Na_n23 >= 1.34)
                                                    if(Na_n > K_n && Ca_n > K_n && Ca_n23 < 1.5 && Ca_n23 > 0.5 && Na_n23 > 0.5)
                                                        T.group4(i) = {'Sodic-Calcic Amphibole'};
                                                        T.species(i) = {'winchite'};
                                                        Fe3 = 1 - (Al_VI + Cr_n23 + Ti_n23);
                                                        if(Fe3 <=0), Fe3 = 0; end
                                                        Fe2 = Fe_n23 - Fe3;
                                                        if(Fe2 <=0), Fe2 = 0; end
                                                        Nay = 2 - (Ca_n23+Mn_n23);
                                                        if(Nay <=0), Nay = 0; end
                                                        Nax = Na_n23 - Nay;
                                                        if(Nax <=0), Nax = 0; end 
                                                        vac = 1 - (Na_n23 + K_n23);
                                                        if(vac <= 0), vac = 0; end
                                                        [sorted_formula1, total1] = Formula_Output({vac, Nax,K_n23}, {'□', 'Na','K'});
                                                        [sorted_formula2, total2] = Formula_Output({Ca_n23, Nay, Mn_n23}, {'Ca', 'Na','Mn'});
                                                        [sorted_formula3, total3] = Formula_Output({Fe2, Mg_n23}, {'Fe2+','Mg'});
                                                        [sorted_formula4, total4] = Formula_Output({Fe3, Al_VI,Ti_n23,Cr_n23},{'Fe3+','AlVI','Ti','Cr'});
                                                        [sorted_formula5, total5] = Formula_Output({Si_n23,Al_IV},{'Si','AlIV'});
                                                        [sorted_formula6, total6] = Formula_Output({2-(Cl_n23+F_n23), Cl_n23, F_n23},{'OH','Cl','F'});
                                                        T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, sorted_formula4, total4, sorted_formula5, total5, 'O22', sorted_formula6, total6)};
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end

                                if(Al_IV / Si_n23 > 0.2-0.2*t && Al_IV / Si_n23 < 0.333+0.333*t)
                                    if(Al_VI >= 0 && Al_VI < 1 && (Mg_n23 + Fe_n23 + Ti_n23 + Al_VI) > 5 - 2*5*(0.1) && (Mg_n23 + Fe_n23 + Ti_n23 + Al_VI) < 5 + 2*5*(0.1))
                                        if(Na_n23 + K_n23 > 1-2*1*t && Na_n23 + K_n23 < 1.5+(1/2)*1.5*t && Na_n > K_n)
                                            if((Na_n + K_n) / (Ca_n + Mn_n) >= 0.5 - 3*0.5*t && (Na_n + K_n) / (Ca_n + Mn_n) < 1 + 1*t)
                                                
                                                Fe3 = 1-(Al_VI+Cr_n23);
                                                if(Fe3<=0), Fe3 = 0; end
                                                Fe2 = Fe_n23 - Fe3;
                                                if(Fe2<=0), Fe2 = 0; end
                                                OH = 2 - (Cl_n23+F_n23);
                                                if(OH<0), OH=0; end
                                                vac = 1 - (Na_n23 + K_n23);
                                                if(vac <= 0), vac = 0; end
                                                [sorted_formula1, total1] = Formula_Output({Na_n23, K_n23, vac}, {'Na', 'K','□'});
                                                [sorted_formula2, total2] = Formula_Output({Ca_n23, Mn_n23}, {'Ca', 'Mn'});
                                                [sorted_formula3, total3] = Formula_Output({Mg_n23, Fe2}, {'Mg', 'Fe2+'});
                                                [sorted_formula4, total4] = Formula_Output({Fe3, Al_VI, Cr_n23}, {'Fe3+', 'Al','Cr'});
                                                [sorted_formula5, total5] = Formula_Output({Si_n23,Al_IV}, {'Si', 'AlIV'});
                                                [sorted_formula6, total6] = Formula_Output({OH, Cl_n23, F_n23}, {'OH','Cl','F'});
                                                T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3,sorted_formula4, total4, sorted_formula5, total5,sorted_formula6, total6)};

                                                if(Fe_n / (Ti_n + Al_VI + Mg_n) >= 1 && (Fe2+Mg_n) > Fe3 && (Fe2+Mg_n)/Fe3 >= 4 - 4*t && Fe3 > 0.5 && Fe3 > Al_VI)
                                                    T.species(i) = {'hastingsite'}; 
                                                end
                                                
                                                if(Mg_n / (Ti_n + Al_VI + Fe_n) >= 1 && Fe3 > Fe2)
                                                    T.species(i) = {'magnesio-hastingsite'}; 
                                                end
                                            end
                                        end
                                    end
                                end

                                if(Al_IV / Si_n23 > 0.6/7.4 - 2*(0.6/7.4)*t && Al_IV / Si_n23 < 1/7 + 4*(1/7)*t)
                                    if(Al_IV > 0.55-0.55*t && Al_IV < 1.5+1.5*t)
                                        if(Al_VI > 0.5-0.5*t && Al_VI < 1+1*t && (Fe_n23 + Al_VI + Mg_n23 + Ti_n23 + Cr_n23 + Mn_n23)/(Si_n23+Al_IV) > 5/8-(5/8)*t && (Fe_n23 + Al_VI + Mg_n23 + Ti_n23 + Cr_n23 + Mn_n23)/(Si_n23+Al_IV) > 5/8-(5/8)*t)
                                            if(Na_n23 + K_n23 < 0.65 + 0.65*t && Ca_n23 > 1.5 - 1.5*t && Ca_n23+K_n23+Na_n23 > 2-2*t && Ca_n23+K_n23+Na_n23 < 2+2*2*t)
                                                if(Ca_n23 > (Na_n23+K_n23) && Al_IV + Si_n23 > 8 - 8*t && Al_IV + Si_n23 < 8 + 8*t)
                                                    if(Mg_n23 + Fe_n23 > 4-4*t && round(Mg_n23 + Fe_n23,1) <= 4+4*t)

                                                        if(Mg_n > Fe_n)
                                                            Fe2 = Fe_n23;
                                                            T.species(i) = {'magnesio-hornblende'};
                                                            OH = 2 - (Cl_n23+F_n23);
                                                            if(OH<0), OH=0; end
                                                            vac = 1 - (Na_n23 + K_n23);
                                                            if(vac <= 0), vac = 0; end
                                                            [sorted_formula1, total1] = Formula_Output({Na_n23, K_n23, vac}, {'Na', 'K','□'});
                                                            [sorted_formula2, total2] = Formula_Output({Ca_n23, Mn_n23}, {'Ca', 'Mn'});
                                                            [sorted_formula3, total3] = Formula_Output({Mg_n23, Fe2}, {'Mg', 'Fe2+'});
                                                            [sorted_formula4, total4] = Formula_Output({Fe2, Al_VI, Cr_n23}, {'Fe3+', 'Al','Cr'});
                                                            [sorted_formula5, total5] = Formula_Output({Si_n23,Al_IV}, {'Si', 'AlIV'});
                                                            T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3,sorted_formula4, total4, sorted_formula5, total5, 'O22','OH_', num2str(round(OH,2)) )};
                                                        end

                                                        if(Fe_n > Mg_n)
                                                            Fe3 = 1 - (Al_VI + Cr_n23 + Ti_n23);
                                                            if(Fe3 <=0), Fe3 = 0; end
                                                            Fe2 = Fe_n23 - Fe3;
                                                            if(Fe2 <=0), Fe2 = 0; end
                                                            Nay = 2 - Ca_n23;
                                                            if(Nay <=0), Nay = 0; end
                                                            Nax = Na_n23 - Nay;
                                                            if(Nax <=0), Nax = 0; end
                                                            vac = 1 - (Nax + K_n23);
                                                            if(vac <= 0), vac = 0; end 
                                                            [sorted_formula1, total1] = Formula_Output({vac, Nax}, {'□', 'Na'});
                                                            [sorted_formula2, total2] = Formula_Output({Ca_n23, Nay}, {'Ca', 'Na'});
                                                            [sorted_formula3, total3] = Formula_Output({Fe2, Mg_n23, Mn_n23, Ni_n23}, {'Fe2+','Mg','Mn','Ni'});
                                                            [sorted_formula4, total4] = Formula_Output({Fe3, Al_VI,Ti_n23,Cr_n23},{'Fe3+','AlVI','Ti','Cr'});
                                                            [sorted_formula5, total5] = Formula_Output({Si_n23,Al_IV},{'Si','AlIV'});
                                                            [sorted_formula6, total6] = Formula_Output({2-(Cl_n23+F_n23), Cl_n23, F_n23},{'OH','Cl','F'});
                                                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3,total3,sorted_formula4,total4,sorted_formula5,total5,'O22',sorted_formula6,total6)};
                                                            T.species(i) = {'ferro-hornblende'};

                                                            if(Fe3 > Al_VI)
                                                                T.species(i) = {'ferro-ferri-hornblende'};
                                                            end

                                                        end

                                                    end
                                                end
                                            end
                                        end
                                    end
                                end

                                if(Al_IV / Si_n23 > (1.3/6.7)-(1.3/6.7)*t && Al_IV / Si_n23 < 0.333+0.333*t)
                                    if(Al_IV > 1.3 - 1.3*t && Al_IV < 2+2*t)
                                        if(Al_VI > 0.5 - 0.5*t && Al_VI < 1+1*t)
                                            if(Na_n23 + K_n23 > 0.5 && Na_n23 + K_n23 < 1+1*t && Na_n > K_n)
                                                if( (Na_n+K_n+Ca_n)/(Mg_n+Fe_n+Ti_n+Al_VI) >= (2.5/5) - (2.5/5)*t && (Na_n+K_n+Ca_n)/(Mg_n+Fe_n+Ti_n+Al_VI) <= (3/5) + (3/5)*t )
                                                    if(Mg_n23 + Fe_n23 > 4-4*t && Mg_n23 + Fe_n23 < 4.5 + 4.5*t)
                                                        if(Mg_n > Fe_n)
                                                            T.species(i) = {'pargasite'};
                                                            OH = 2 - (Cl_n23+F_n23);
                                                            if(OH <=0), OH = 0; end
                                                            Fe3 = 1 - (Al_VI + Cr_n23 + Ti_n23);
                                                            if(Fe3 <= 0), Fe3 = 0; end
                                                            Fe2 = Fe_n23 - Fe3;
                                                            if(Fe2 <=0), Fe2 = 0; end
                                                            vac = 1 - (Na_n23 + K_n23);
                                                            if(vac <= 0), vac = 0; end
                                                            [sorted_formula1, total1] = Formula_Output({Na_n23, K_n23, vac}, {'Na', 'K','□'});
                                                            [sorted_formula2, total2] = Formula_Output({Ca_n23, Mn_n23}, {'Ca', 'Mn'});
                                                            [sorted_formula3, total3] = Formula_Output({Mg_n23, Fe2}, {'Mg', 'Fe2+'});
                                                            [sorted_formula4, total4] = Formula_Output({Fe3, Al_VI, Cr_n23, Ti_n23}, {'Fe3+', 'Al_VI','Cr','Ti'});
                                                            [sorted_formula5, total5] = Formula_Output({Si_n23,Al_IV}, {'Si', 'AlIV'});
                                                            [sorted_formula6, total6] = Formula_Output({OH,F_n23, Cl_n23}, {'OH', 'F','Cl'});
                                                            T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3,sorted_formula4, total4, sorted_formula5, total5, 'O22',sorted_formula6, total6)};
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end

                                if(Al_IV > 1-1*t && Al_IV < 1+1*t)
                                     if(Na_n23 + Ca_n23 > 3-3*t && Na_n23 + Ca_n23 < 3+3*t && Na_n23 > 1 - 1*t && Na_n23 < 1 + 1*t && Ca_n23 > 2-2*t && Ca_n23 < 2+2*t)
                                        if(Ca_n > Na_n && Si_n23 + Al_IV > 8 - 8*(1/2)*t && Si_n23 + Al_IV < 8 + 8*(1/2)*t)
                                            if((Mg_n23 + Fe_n23 + Al_VI + Cr_n23 + Ti_n23 + Mn_n23) > 5-5*t && (Mg_n23 + Fe_n23 + Al_VI + Cr_n23 + Ti_n23 + Mn_n23) < 5+5*t)

                                                if(Mg_n23/(Mg_n23 + Fe_n23 + Al_VI + Cr_n23 + Ti_n23 + Mn_n23) > 0.5)
                                                    T.species(i) = {'edenite'};
                                                    Nax = 1 - K_n23;
                                                    if(Nax < 0), Nax = 0; end
                                                    Nay = Na_n23 - Nax;
                                                    OH = 2 - (Cl_n23+F_n23);
                                                    if(OH<0), OH=0; end
                                                    [sorted_formula1, total1] = Formula_Output({Nax, K_n23}, {'Nax','K'});
                                                    [sorted_formula2, total2] = Formula_Output({Ca_n23, Nay}, {'Ca','Nay'});
                                                    [sorted_formula3, total3] = Formula_Output({Mg_n23, Fe_n23, Mn_n23, Cr_n23, Ti_n23,Al_VI}, {'Mg','Fe','Mn','Cr','Ti','ALVI'});
                                                    [sorted_formula4, total4] = Formula_Output({Si_n23, Al_IV}, {'Si','AlIV'});
                                                    [sorted_formula5, total5] = Formula_Output({OH, Cl_n23, F_n23}, {'OH','Cl','F'});
                                                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3,total3,sorted_formula4,total4,'O22',sorted_formula5,total5)};
                                                end

                                                if(Fe_n23/(Mg_n23 + Fe_n23 + Al_VI + Cr_n23 + Ti_n23 + Mn_n23) > 0.5)
                                                    T.species(i) = {'ferro-edenite'};
                                                    Cay = 2 - Mn_n23;
                                                    if(Cay<0), Cay=0; end
                                                    Cax = Ca_n23 - Cay;
                                                    if(Cax<0), Cax=0; end
                                                    [sorted_formula1, total1] = Formula_Output({Na_n23, K_n23, Cax},{'Na','K','Ca'});
                                                    [sorted_formula2, total2] = Formula_Output({Cay, Mn_n23},{'Ca','Mn'});
                                                    [sorted_formula3, total3] = Formula_Output({Si_n23, Al_n23},{'Si','Al'});
                                                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, 'Fe_', num2str(round(Fe_n23, 2)), sorted_formula3, total3, 'O22(OH)2')};
                                                end

                                            end
                                        end
                                    end
                                end

                                if(Ti_n23 > 0.5-0.5*t && Ti_n23 < 1+1*t)
                                    if(Al_IV * 24 / 23 > 2 - 2*2*t && Al_IV * 24 / 23 < 2 + 2*2*t)
                                        if(Na_n + Ca_n > 3-3*t && Na_n + Ca_n < 3+3*t)
                                            if(Ca_n > Na_n)
                                                if(Mg_n + Fe_n + Al_VI*24/23 + Ti_n > 5-5*t && Mg_n + Fe_n + Al_VI*24/23 + Ti_n < 5+5*t)
                                                    if(Mg_n > Fe_n)
                                                        T.species(i) = {'kaersutite'};
                                                        OH = 2 - (F_n+Cl_n);
                                                        if(OH<=0), OH = 0; end
                                                        vac = 1 - (Na_n+K_n);
                                                        if(vac<=0), vac = 0; end
                                                        Fe3 = 2 - (Al_VI+Cr_n+Ti_n);
                                                        if(Fe3<=0), Fe3 = 0; end
                                                        Fe2 = Fe_n - Fe3;
                                                        if(Fe2<=0), Fe2 = 0; end
                                                        [sorted_formula1, total1] = Formula_Output({K_n, Na_n, vac}, {'K','Na','□'});
                                                        [sorted_formula2, total2] = Formula_Output({Ca_n, Mn_n},{'Ca', 'Mn'});
                                                        [sorted_formula3, total3] = Formula_Output({Mg_n, Fe2},{'Mg','Fe2+'});
                                                        [sorted_formula4, total4] = Formula_Output({Al_VI, Fe3, Ti_n, Cr_n},{'AlVI','Fe3+','Ti','Cr'});
                                                        [sorted_formula5, total5] = Formula_Output({Si_n, Al_IV},{'Si','AlIV'});
                                                        [sorted_formula6, total6] = Formula_Output({OH, F_n, Cl_n},{'OH','F','Cl'});
                                                        T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, sorted_formula3,total3, sorted_formula4,total4,sorted_formula5,total5,'O22',sorted_formula6,total6)};
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end

                            end
                        end
                    end
                end
            end

            %% 7.7.1.2. Non-Calcic amphiboles
            if(~strcmp(T.group4(i), 'Calcic Amphibole') && S_n <= 0.05 && P_n <= 0.05)

                Mgy = 5 - (Al_VI + Cr_n23 + Ti_n23 + Fe_n23);
                if(Mgy <= 0), Mgy = 0; end

                Fex = 2 - (Ca_n23+Na_n23+K_n23);
                if(Fex <= 0), Fex = 0; end
                Fey = Fe_n23 - Fex;
                if(Fey <= 0), Fey = 0; end

                if((Mg_n23 + Fe_n23 + Mn_n23 + Cr_n23 + Al_VI + Ti_n23)/(Si_n23 + Al_IV) > (5/8)-2*(5/8)*t && (Mg_n23 + Fe_n23 + Mn_n23 + Cr_n23 + Al_VI + Ti_n23)/(Si_n23 + Al_IV) < (5/8)+2*(5/8)*t)
                    if((K_n23+Na_n23+Ca_n23)/(Si_n23+Al_IV) > (3/8)-(3/8)*2*t && (K_n23+Na_n23+Ca_n23)/(Si_n23+Al_IV) < (3/8)+(3/8)*2*t)
                        if((K_n23+Ca_n23+Na_n23) > 3 - 3*2*t && (K_n23+Ca_n23+Na_n23) < 3 + 3*2*t  && (Al_IV/Si_n23) < (1/8) + (1/8)*2*t)
                            if(T.SiO2(i) >= 55)
                                T.group3(i) = {'Amphibole'};
                                T.group4(i) = {'Sodic-Calcic Amphibole'};
                                T.species(i) = {'potassic-richterite'};
                                Nax = 1 - (K_n23+Na_n23);
                                if(Nax <= 0), Nax = 0; end
                                Nay = Na_n23 - Nax;
                                if(Nay <= 0), Nay = 0; end
                                vac = 2 - (Nax + Ca_n23);
                                if(vac <= 0), vac = 0; end
                                OH = 2 - (F_n23+Cl_n23);
                                if(OH <= 0), OH = 0; end
                                [sorted_formula1, total1] = Formula_Output({Nax, K_n23, vac}, {'Na','K','□'});
                                [sorted_formula2, total2] = Formula_Output({Nay, Ca_n23, Mn_n23}, {'Na', 'Ca','Mn'});
                                [sorted_formula3, total3] = Formula_Output({Mg_n23, Fe_n23, Al_VI, Ti_n23, Cr_n23}, {'Mg','Fe','AlVI','Ti','Cr'});
                                [sorted_formula4, total4] = Formula_Output({Si_n23, Al_IV}, {'Si','AlIV'});
                                [sorted_formula5, total5] = Formula_Output({OH, F_n23, Cl_n23}, {'OH','F','Cl'});
                                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3, sorted_formula4, total4, 'O22', sorted_formula5, total5)};
                            end
                        end
                    end
                end

                if(T.CaO(i) < 2 && S_n <= 0.05 && P_n <= 0.05)
                    if(Al_n / Si_n >= 0 && Al_n / Si_n < (2/3)+(2/3)*t)
                        if(Fe_n23 < 5)
                            if((Mgy + Fe_n23 + Al_VI + Cr_n23 + Ti_n23 + Mn_n23) > 5-5*t && (Mgy + Fe_n23 + Al_VI + Cr_n23 + Ti_n23 + Mn_n23) < 5+5*t)
                                T.group2(i) = {'Inosilicate'};
                                T.group3(i) = {'Amphibole'};
                                T.group4(i) = {''};
                                T.species(i) = {''};
                            end
                        elseif(Fe_n23 >= 5)
                            if((Mgy + Fey + Al_VI + Cr_n23 + Ti_n23 + Mn_n23) > 5-5*t && (Mgy + Fey + Al_VI + Cr_n23 + Ti_n23 + Mn_n23) < 5+5*t)
                                T.group2(i) = {'Inosilicate'};
                                T.group3(i) = {'Amphibole'};
                                T.group4(i) = {''};
                                T.species(i) = {''};
                            end
                        end

                        if(Ca_n/(Na_n+K_n) <= 0.25 + 0.25*t && (K_n/Na_n)/(Mg_n+Fe_n+Al_VI+Ti_n) < 1/5 - 2*(1/5)*t && (K_n/Na_n)/(Mg_n+Fe_n+Al_VI+Ti_n) <= 2/5 - 2*(2/5)*t)
                            T.group4(i) = {'Non-Calcic Alk Amphibole'};
                        else
                            T.group4(i) = {'Non-Calcic Non-Alk Amphibole'};
                        end

                            if(Al_n23 + Si_n23 +Al_IV > 10-10*t && Al_n23 + Si_n23+Al_IV < 10+10*t)
                                if(Mg_n23 + Fe_n23+Al_VI > 5-5*t && Mg_n23 + Fe_n23+Al_VI < 5+5*t)
                                    if(Mg_n > Fe_n && Mg_n23 + Fe_n23 > 3-3*t && Mg_n23 + Fe_n23 < 3+3*t)
                                        if(Na_n23 + K_n23 +Ca_n23 > 2-2*t && Na_n23 + K_n23+Ca_n23 < 2+2*t)
                                            if(Na_n > K_n && Na_n > Ca_n && Al_IV/Si_n23 < 1/8)
                                                T.species(i) = {'glaucophane'};
                                                Fe3 = 2 -(Al_VI + Ti_n23 + Cr_n23);
                                                if(Fe3 <= 0), Fe3 = 0; end   
                                                Fe2 = Fe_n23 - Fe3;
                                                if(Fe2 <= 0), Fe2 = 0; end
                                                [sorted_formula1, total1] = Formula_Output({Na_n23, Ca_n23, K_n23}, {'Na','Ca','K'});
                                                [sorted_formula2, total2] = Formula_Output({Mg_n23, Fe2, Mn_n23}, {'Mg', 'Fe2+','Mn'});
                                                [sorted_formula3, total3] = Formula_Output({Al_VI, Fe3, Ti_n23, Cr_n23}, {'AlVI','Fe3+','Ti','Cr'});
                                                [sorted_formula4, total4] = Formula_Output({Si_n23, Al_IV}, {'Si','AlIV'});
                                                [sorted_formula5, total5] = Formula_Output({2 - (F_n23+Cl_n23), F_n23, Cl_n23}, {'OH','F', 'Cl'}); 
                                                T.formula(i) = {strcat('◻', sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3, sorted_formula4, total4, 'O22', sorted_formula5, total5)};
                                            end
                                        end
                                    end
                                end
                            end

                            if(Si_n23 + Al_n23 > 8-8*t && Si_n23 + Al_n23 < 8+8*t)
                                if(Al_n23 < 0.05 && Ti_n23 < 0.05 && Cr_n23 < 0.05 && Si_n23 >= 7 && Mg_n23 >= 6 && Na_n23 + Ca_n23 + K_n23 <= 0.15)
                                    if(Mg_n23 + Fe_n23 + Ca_n23 + Mn_n23 + Na_n23 > 7-7*t && Mg_n23 + Fe_n23 + Ca_n23 + Mn_n23 + Na_n23 < 7+7*t)
                                        if(Mg_n > Fe_n && Ca_n23 < 0.02)
                                            T.species(i) = {'anthophyllite/cummingtonite'};     
                                            Mgy = 5 - (Fe_n23+Ni_n23+Al_VI+Ti_n23+Cr_n23);
                                            Mgx = Mg_n23-Mgy;
                                            if(Mgx < 0), Mgx = 0; end
                                            vac = 8 - (Na_n23 + K_n23);
                                            if(vac < 0), vac = 0; end
                                            [sorted_formula1, total1] = Formula_Output({Na_n23, K_n23, vac},{'Na','K','□'});
                                            [sorted_formula2, total2] = Formula_Output({Mgx,Ca_n23,Mn_n23},{'Mg','Ca','Mn'});
                                            [sorted_formula3, total3] = Formula_Output({Mgy,Fe_n23,Ni_n23,Al_VI,Ti_n23,Cr_n23},{'Mg','Fe','Ni','AlVI','Ti','Cr'});
                                            [sorted_formula4, total4] = Formula_Output({Si_n23,Al_IV},{'Si','AlIV'});
                                            [sorted_formula5, total5] = Formula_Output({2-(Cl_n23+F_n23), Cl_n23, F_n23},{'OH','Cl','F'});
                                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3,total3,sorted_formula4,total4,'O22',sorted_formula5,total5)};
                                        end
                                    end
                                end


                                if(Al_n23 < 0.05 && Si_n23 >= 7 && Fe_n23 >= 5 && Na_n23 + Ca_n23 + K_n23 <= 0.15)
                                    if(Mg_n23 + Fe_n23 + Ca_n23 + Mn_n23 + Na_n23 > 7-7*t && Mg_n23 + Fe_n23 + Ca_n23 + Mn_n23 + Na_n23 < 7+7*t)
                                        if(Mg_n < Fe_n)
                                            T.species(i) = {'ferro-anthophyllite/grunerite'};
                                            vac = 1      - (Na_n23 + K_n23);
                                            if(vac < 0), vac = 0; end
                                            [sorted_formula1, total1] = Formula_Output({Na_n23, K_n23, vac},{'Na','K','□'});
                                            [sorted_formula2, total2] = Formula_Output({Mg_n23,Fe_n23,Mn_n23},{'Mg','Fe','Mn'});
                                            [sorted_formula3, total3] = Formula_Output({Fe_n23, Ni_n23,Al_VI,Ti_n23,Cr_n23},{'Fe', 'Ni','AlVI','Ti','Cr'});
                                            [sorted_formula4, total4] = Formula_Output({Si_n23,Al_IV},{'Si','AlIV'});
                                            [sorted_formula5, total5] = Formula_Output({2-(Cl_n23+F_n23), Cl_n23, F_n23},{'OH','Cl','F'});
                                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3,total3,sorted_formula4,total4,'O22',sorted_formula5,total5)};
                                        end

                                    end
                                end

                            end

                            if(Al_IV / Si_n > (1.5/6)-(1.5/6)*t && Al_IV / Si_n < (2/6)+(2/6)*t)
                                if(Mg_n23 + Fe_n23 + Mn_n23 + Ti_n23 + Cr_n23 + Al_VI > 7-7*t && Mg_n23 + Fe_n23 + Mn_n23 + Ti_n23 + Cr_n23 + Al_VI  < 7+7*t)
                                    if(Fe_n23 < 4 && Si_n23 <= 7 && Al_n23 > 2 - 2*t && Al_n23 < 4 + 4*t)
                                        T.species(i) = {'gedrite'};
                                        Fe3 = 2 - (Al_VI + Cr_n23 + Ti_n23);
                                        if(Fe3 <= 0), Fe3 = 0; end
                                        Fe2 = Fe_n23 - Fe3;
                                        if(Fe2 <= 0), Fe2 = 0; end
                                        Mga = 2 - (Fe_n23 + Mn_n23);
                                        if(Mga <= 0), Mga = 0; end
                                        Mgb = Mg_n23 - Mga;
                                        if(Mgb <= 0), Mgb = 0; end
                                        OH = 4 - (F_n23 + Cl_n23);
                                        if(OH <= 0), OH = 0; end
                                        vac = 1 - (Na_n23 + K_n23 + Ca_n23);
                                        if(vac < 0), vac = 0; end
                                        [sorted_formula1, total1] = Formula_Output({vac, Ca_n23, Na_n23, K_n23}, {'□', 'Ca', 'Na','K'});
                                        [sorted_formula2, total2] = Formula_Output({Fe2, Mga, Mn_n23}, {'Fe2+','Mg','Mn'});
                                        [sorted_formula3, total3] = Formula_Output({Fe3, Al_VI,Ti_n23,Cr_n23},{'Fe2+','AlVI','Ti','Cr'});
                                        [sorted_formula4, total4] = Formula_Output({Si_n23,Al_IV},{'Si','AlIV'});
                                        [sorted_formula5, total5] = Formula_Output({OH, Cl_n23, F_n23},{'OH','Cl','F'});
                                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, '(Mg_',num2str(round(Mgb,2)),sorted_formula3,total3, sorted_formula4,total4, 'O22',sorted_formula5,total5)};
                                    end

                                    if(Fe_n23 >= 4)
                                        T.species(i) = {'ferro-gedrite'};
                                        Fex = 5 - (Al_VI + Cr_n23 + Ti_n23);
                                        if(Fex <= 0), Fex = 0; end
                                        Fey = Fe_n23 - Fex;
                                        if(Fey <=0), Fey = 0; end
                                        vac = 8 - (Na_n23 + Ca_n23);
                                        if(vac < 0), vac = 0; end
                                        [sorted_formula1, total1] = Formula_Output({vac, Ca_n23, Na_n23}, {'□', 'Ca', 'Na'});
                                        [sorted_formula2, total2] = Formula_Output({Fey, Mg_n23, Mn_n23, Ni_n23}, {'Fe2+','Mg','Mn','Ni'});
                                        [sorted_formula3, total3] = Formula_Output({Fex, Al_VI,Ti_n23,Cr_n23},{'Fe2+','AlVI','Ti','Cr'});
                                        [sorted_formula4, total4] = Formula_Output({Si_n23,Al_IV},{'Si','AlIV'});
                                        [sorted_formula5, total5] = Formula_Output({2-(Cl_n23+F_n23), Cl_n23, F_n23},{'OH','Cl','F'});
                                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3,total3,sorted_formula4,total4,'O22',sorted_formula5,total5)};
                                    end

                                end
                            end

                            if(Si_n23 + Al_IV > 8-8*t && Si_n23 + Al_IV < 8+8*t)
                                if(Al_n23 < 0.5)
                                    if(Mg_n23 + Fe_n23 > 5-5*t && Mg_n23 + Fe_n23 < 5+5*t && Na_n23 > K_n23)

                                        if(Na_n23 + K_n23 > 2-2*t && Na_n23 + K_n23 < 2+2*t)
                                            T.species(i) = {'riebeckite'};
                                            Nax = 1 - (K_n23+Na_n23);
                                            if(Nax <= 0), Nax = 0; end
                                            Nay = Na_n23 - Nax;
                                            if(Nay <= 0), Nay = 0; end
                                            Fe2 = 3 - Mg_n23;
                                            if(Fe2 <= 0), Fe2 = 0; end
                                            Fe3 = Fe_n23 - Fe2;
                                            if(Fe3 <= 0), Fe3 = 0; end
                                            vac = 2 - (Nax + Ca_n23);
                                            if(vac <= 0), vac = 0; end
                                            OH = 2 - (F_n23+Cl_n23);
                                            if(OH <= 0), OH = 0; end
                                            [sorted_formula1, total1] = Formula_Output({Nax, K_n23, vac}, {'Na','K','□'});
                                            [sorted_formula2, total2] = Formula_Output({Nay, Ca_n23, Mn_n23}, {'Na', 'Ca','Mn'});
                                            [sorted_formula3, total3] = Formula_Output({Fe2, Mg_n23}, {'Fe2+','Mg'});
                                            [sorted_formula4, total4] = Formula_Output({Fe_n23, Al_VI, Ti_n23, Cr_n23}, {'Fe3+','AlVI','Ti','Cr'});
                                            [sorted_formula5, total5] = Formula_Output({Si_n23, Al_IV}, {'Si','AlIV'});
                                            [sorted_formula6, total6] = Formula_Output({OH, F_n23, Cl_n23}, {'OH','F','Cl'});
                                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3, sorted_formula4, total4, sorted_formula5, total5, 'O22', sorted_formula6, total6)};
                                        end

                                        if(Na_n23 + K_n23 > 3-3*t && Na_n23 + K_n23 < 3+3*t)
                                            T.species(i) = {'arfvedsonite'};
                                            Nax = 1 - K_n23;
                                            if(Nax<0), Nax = 0; end
                                            Nay = Na_n23 - Nax;
                                            if(Nay<0), Nay = 0; end
                                            Fe2 = 4 - Mg_n23;
                                            if(Fe2<0), Fe2=0; end
                                            Fe3 = Fe_n23-Fe2;
                                            if(Fe3<0), Fe3=0; end
                                            [sorted_formula1, total1] = Formula_Output({Nax,K_n23},{'Na','K'});
                                            [sorted_formula2, total2] = Formula_Output({Ca_n23,Mn_n23,Nay},{'Ca','Mn','Na'});
                                            [sorted_formula3, total3] = Formula_Output({Fe2,Mg_n23},{'Fe2+','Mg'});
                                            [sorted_formula4, total4] = Formula_Output({Fe3,Al_VI,Ti_n23,Cr_n23},{'Fe3+','AlVI','Ti','Cr'});
                                            [sorted_formula5, total5] = Formula_Output({Si_n23,Al_IV},{'Si','AlIV'});
                                            [sorted_formula6, total6] = Formula_Output({2-(Cl_n23+F_n23), Cl_n23,F_n23},{'OH','Cl','F'});
                                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3,total3,sorted_formula4,total4,'[',sorted_formula5,total5,'O22]',sorted_formula6,total6)};
                                        end

                                    end                                     
                                end
                            end
                    end
                end

                if(Na_n > Ca_n && Na_n > K_n && Ca_n > K_n && Ca_n/Na_n > (1/2) - (1/2)*2*t && Ca_n/Na_n <= (1/2) + (1/2)*2*t)
                    if(Na_n23 + Ca_n23 > 3-3*t && Na_n23 + Ca_n23 < 3+3*t)
                        if(Mg_n23 + Fe_n23 > 5-5*t && Mg_n23 + Fe_n23 < 5+5*t)
                            if(Na_n > Ca_n && Al_n23 < 0.5)
                                if(Mg_n > Fe_n)
                                    T.group4(i) = {'Sodic-Calcic Amphibole'};
                                    T.species(i) = {'richterite'};
                                    Al_IV = 8 - Si_n23;
                                    if(Al_IV < 0), Al_IV = 0; end
                                    Al_VI = Al_n23 - Al_IV;
                                    if(Al_VI < 0), Al_VI = 0; end
                                    Nax = 1 - (K_n23+Na_n23);
                                    if(Nax <= 0), Nax = 0; end
                                    Nay = Na_n23 - Nax;
                                    if(Nay <= 0), Nay = 0; end
                                    vac = 2 - (Nax + Ca_n23);
                                    if(vac <= 0), vac = 0; end
                                    OH = 2 - (F_n23+Cl_n23);
                                    if(OH <= 0), OH = 0; end
                                    [sorted_formula1, total1] = Formula_Output({Nax, K_n23, vac}, {'Na','K','□'});
                                    [sorted_formula2, total2] = Formula_Output({Nay, Ca_n23, Mn_n23}, {'Na', 'Ca','Mn'});
                                    [sorted_formula3, total3] = Formula_Output({Mg_n23, Fe_n23, Al_VI, Ti_n23, Cr_n23}, {'Mg','Fe','AlVI','Ti','Cr'});
                                    [sorted_formula4, total4] = Formula_Output({Si_n23, Al_IV}, {'Si','AlIV'});
                                    [sorted_formula5, total5] = Formula_Output({OH, F_n23, Cl_n23}, {'OH','F','Cl'});
                                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3, sorted_formula4, total4, 'O22', sorted_formula5, total5)};
                                end

                                if(Mg_n <= Fe_n)
                                    T.group4(i) = {'Sodic-Calcic Amphibole'};
                                    T.species(i) = {'ferro-richterite'};
                                    Nax = 1 - K_n23;
                                    if(Nax <= 0), Nax = 0; end
                                    Nay = Na_n23 - Nax;
                                    if(Nay <= 0), Nay = 0; end
                                    [sorted_formula1, total1] = Formula_Output({Nax, K_n23}, {'Na','K'});
                                    [sorted_formula2, total2] = Formula_Output({Nay, Ca_n23, vac}, {'Na', 'Ca','□'});
                                    [sorted_formula3, total3] = Formula_Output({Si_n23, Al_IV}, {'Si','AlIV'});
                                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, 'Fe2+_', num2str(round(Fe_n23,2)), sorted_formula3,total3,'O22(OH)2')};
                                end

                            end
                        end
                    end
                end
            end
        end %end of Amphiboles

        % Catch
        if(strcmp(T.group3(i), 'Amphibole'))
            Renorm_factor = 24 / 11;
            Mg_n11 = T.Mg_Multi(i) / Renorm_factor;
            Mn_n11 = T.Mn_Multi(i) / Renorm_factor;
            Fe_n11 = T.Fe_Multi(i) / Renorm_factor;
            Ca_n11 = T.Ca_Multi(i) / Renorm_factor;
            Si_n11 = T.Si_Multi(i) / Renorm_factor;
            Ni_n11 = T.Ni_Multi(i) / Renorm_factor;
            Ti_n11 = T.Ti_Multi(i) / Renorm_factor;
            Cr_n11 = T.Cr_Multi(i) / Renorm_factor;
            Na_n11 = T.Na_Multi(i) / Renorm_factor;
            Al_n11 = T.Al_Multi(i) / Renorm_factor;
            Cl_n11 = T.Cl_Multi(i) / Renorm_factor;
            F_n11 = T.F_Multi(i) / Renorm_factor;
            Zn_n11 = T.Zn_Multi(i) / Renorm_factor;
            V_n11 = T.V_Multi(i) / Renorm_factor;
            K_n11 = T.K_Multi(i) / Renorm_factor;

            Al_IV = 4 - Si_n11;
            if(Al_IV < 0), Al_IV = 0; end
            Al_VI = Al_n11 - Al_IV;
            if(Al_VI < 0), Al_VI = 0; end

            if((Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11)/(Si_n11+Al_IV) >= (3/4)-(3/4)*1.5*t && (Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11)/(Si_n11+Al_IV) <= (3/4)+(3/4)*t)
                if(Fe_n11/(Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11) > 0.5 && Al_IV/Si_n11 < 0.1 && (Na_n11+K_n11+Ca_n11)/(Si_n11+Al_IV) < 0.1)
                    if(ABCDT_n/24 > 7/11 - 2*(7/11)*t && ABCDT_n/24 < 7/11 + (7/11)*t)
                        T.group2(i) = {'Phyllosilicate'};
                        T.group3(i) = {'Talc-Pyrophyllite'};
                        T.group4(i) = {''};
%                         T.species(i) = {'minnesotaite'};
                        OH = 2 - (F_n11+Cl_n11);
                        if(OH<=0), OH = 0; end
                        [sorted_formula1, total1] = Formula_Output({Fe_n11, Mg_n11, Ni_n11, Mn_n11}, {'Fe', 'Mg', 'Ni','Mn'});
                        [sorted_formula2, total2] = Formula_Output({Si_n11, Al_n11},{'Si','Al'});
                        [sorted_formula3, total3] = Formula_Output({OH, F_n11, Cl_n11},{'OH','F','Cl'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O10',sorted_formula3, total3)};

                    end
                end
            end

            if((Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11+Ti_n11)/(Si_n11+Al_IV) >= 1.5/4 - (1.5/4)*t && (Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11+Ti_n11)/(Si_n11+Al_IV) <= 3/4 + 3/4*t)
                if((Na_n11+Ca_n11+K_n11)/(Si_n11+Al_IV) > 0.1/4 && (Na_n11+Ca_n11+K_n11)/(Si_n11+Al_IV) < 0.4/4 && Na_n/(Na_n+K_n+Ca_n) > 0.5 && Mg_n11 > Fe_n11 && Al_IV/Si_n11 < 0.1 && Na_n11 > Ca_n11 && Na_n11 > K_n11)
                    T.group2(i) = {'Phyllosilicate'};
                    T.group3(i) = {'Smectite'};
                    T.group4(i) = {''};
%                         T.species(i) = {'hectorite'};
                    renorm = 11/24;
                    Al_IV = 4 - Si_n*renorm;
                    if(Al_IV <=0), Al_IV = 0; end
                    Al_VI = Al_n*renorm - Al_IV;
                    if(Al_VI <=0), Al_VI = 0; end
                    Li = 3 - (Mg_n*renorm + Fe_n*renorm + Al_VI);
                    if(Li <=0), Li = 0; end
                    OH = 2 - (F_n*renorm + Cl_n*renorm);
                    if(OH <=0), OH = 0; end
                    [sorted_formula1, total1] = Formula_Output({Na_n*renorm, Ca_n*renorm, K_n*renorm},{'Na','Ca','K'});
                    [sorted_formula2, total2] = Formula_Output({Mg_n*renorm, Li, Fe_n*renorm, Al_VI, Mn_n*renorm},{'Mg','Li','Fe','AlVI','Mn'});
                    [sorted_formula3, total3] = Formula_Output({Si_n*renorm, Al_IV},{'Si','AlIV'});
                    [sorted_formula4, total4] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2,total2,sorted_formula3,total3,'O10',sorted_formula4,total4,'·nH2O')};
                end
            end

            if((Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11+Ti_n11)/(Si_n11+Al_IV) > 2/4 - 2/4*t && (Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11+Ti_n11)/(Si_n11+Al_IV) < 2/4 + 2/4*t)
                if(Si_n11 > Al_IV && Mg_n11 < Fe_n11 && Ca_n11/(Si_n11+Al_IV) < 0.05 && Na_n11 > Ca_n11 && Na_n11 > K_n11 && (Na_n11+Ca_n11+K_n11)/(Si_n11+Al_IV) > (0.3/4)-(0.3/4)*t && (Na_n11+Ca_n11+K_n11)/(Si_n11+Al_IV) < (0.3/4)+(0.3/4)*t)
                    T.group2(i) = {'Phyllosilicate'};
                    T.group3(i) = {'Smectite'};
                    T.group4(i) = {''};
%                     T.species(i) = {'nontronite'};
                    renorm = 11/24;
                    Al_IV = 4 - Si_n*renorm;
                    if(Al_IV<=0), Al_IV = 0; end
                    Al_VI = Al_n*renorm - Al_IV;
                    if(Al_VI<=0), Al_VI = 0; end
                    OH = 2 - (F_n*renorm + Cl_n*renorm);
                    if(OH<=0), OH = 0; end      
                    [sorted_formula1, total1] = Formula_Output({Na_n*renorm, K_n*renorm, Ca_n*renorm},{'Na','K','Ca'});
                    [sorted_formula2, total2] = Formula_Output({Fe_n*renorm, Al_VI, Mg_n*renorm, Mn_n*renorm},{'Fe','AlVI','Mg','Mn'});
                    [sorted_formula3, total3] = Formula_Output({Si_n*renorm, Al_IV},{'Si','AlIV'});
                    [sorted_formula4, total4] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, sorted_formula3, total3, 'O10', sorted_formula4,total4,'·nH2O')};
                end
            end

            if((Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11)/(Si_n11+Al_IV) > 3/4 - (3/4)*t && (Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11)/(Si_n11+Al_IV) < 3/4 + (3/4)*t)    
                if(Si_n11 > Al_IV && Mg_n11/(Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11) >= 0.5 && (Na_n11+Ca_n11+K_n11)/(Si_n11+Al_IV) < 0.5 + 0.5*t && (Ca_n11+Na_n11) > K_n11 && Ca_n > Na_n)
                    if(Mg_n11 > Fe_n11 && Fe_n11 > Al_VI)
                        T.group2(i) = {'Phyllosilicate'};
                        T.group3(i) = {'Smectite'};
                        T.group4(i) = {''};
%                         T.species(i) = {'saponite'};
                        renorm = 11/24;
                        Al_IV = 4 - Si_n*renorm;
                        if(Al_IV<=0), Al_IV = 0; end 
                        Al_VI = Al_n*renorm - Al_IV;
                        if(Al_VI<=0), Al_VI = 0; end            
                        OH = 2 - (F_n*renorm + Cl_n*renorm);
                        if(OH<=0), OH = 0; end
                        [sorted_formula1, total1] = Formula_Output({Na_n*renorm, K_n*renorm, Ca_n*renorm},{'Na','K','Ca'});
                        [sorted_formula2, total2] = Formula_Output({Fe_n*renorm, Al_VI, Mg_n*renorm, Mn_n*renorm, Ti_n*renorm, Cr_n*renorm},{'Fe','AlVI','Mg','Mn', 'Ti','Cr'});
                        [sorted_formula3, total3] = Formula_Output({Si_n*renorm, Al_IV},{'Si','AlIV'});
                        [sorted_formula4, total4] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                        T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, sorted_formula3, total3, 'O10', sorted_formula4,total4,'·4H2O')};
                    end
                end
            end

            if( (Zn_n*(11/24))/(Si_n11+Al_IV) > 3/4 - (3/4)*2*t && (Zn_n*(11/24))/(Si_n11+Al_IV) < 3/4 + (3/4)*2*t)
                if((Na_n11+Ca_n11+K_n11)/(Si_n11+Al_IV) > (0.3/4)-(0.3/4)*0.1 && (Na_n11+Ca_n11+K_n11)/(Si_n11+Al_IV) < (0.3/4)+(0.3/4)*0.1 && Na_n11 > Ca_n11 && Na_n11 > K_n11)
                    T.group2(i) = {'Phyllosilicate'};
                    T.group3(i) = {'Smectite'};
                    T.group4(i) = {''};
%                     T.species(i) = {'sauconite'};
                    renorm = 11/24;
                    Al_IV = 4 - Si_n*renorm;
                    if(Al_IV<=0), Al_IV = 0; end  
                    Al_VI = Al_n*renorm - Al_IV;
                    if(Al_VI<=0), Al_VI = 0; end
                    OH = 2 - (F_n*renorm + Cl_n*renorm);
                    if(OH<=0), OH = 0; end
                    [sorted_formula1, total1] = Formula_Output({Na_n*renorm, K_n*renorm, Ca_n*renorm},{'Na','K','Ca'});
                    [sorted_formula2, total2] = Formula_Output({Fe_n*renorm, Al_VI, Mg_n*renorm, Mn_n*renorm, Cu_n*renorm, Zn_n*renorm},{'Fe','AlVI','Mg','Mn', 'Cu','Zn'});
                    [sorted_formula3, total3] = Formula_Output({Si_n*renorm, Al_IV},{'Si','AlIV'});
                    [sorted_formula4, total4] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'}); 
                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, sorted_formula3, total3, 'O10', sorted_formula4,total4,'·4H2O')};
                end
            end
            
            Renorm_factor = 24 / 16;
            Mg_n16 = T.Mg_Multi(i) / Renorm_factor;
            Mn_n16 = T.Mn_Multi(i) / Renorm_factor;
            Fe_n16 = T.Fe_Multi(i) / Renorm_factor;
            Ca_n16 = T.Ca_Multi(i) / Renorm_factor;
            Si_n16 = T.Si_Multi(i) / Renorm_factor;
            Ti_n16 = T.Ti_Multi(i) / Renorm_factor;
            Cr_n16 = T.Cr_Multi(i) / Renorm_factor;
            Na_n16 = T.Na_Multi(i) / Renorm_factor;
            Al_n16 = T.Al_Multi(i) / Renorm_factor;
            F_n16 = T.F_Multi(i) / Renorm_factor;
            Cl_n16 = T.Cl_Multi(i) / Renorm_factor;
            Zn_n16 = T.Zn_Multi(i) / Renorm_factor;
            V_n16 = T.V_Multi(i) / Renorm_factor;
            K_n16 = T.K_Multi(i) / Renorm_factor;

            Al_IV = 6 - Si_n16;
            if(Al_IV < 0), Al_IV = 0; end
            Al_VI = Al_n16 - Al_IV;
            if(Al_VI < 0), Al_VI = 0; end

            if((Mg_n16+Mn_n16+Fe_n16+Al_VI+Zn_n16+Cr_n16+V_n16+Ti_n16)/(Si_n16+Al_IV) > 4/6 - 4/6*t && (Mg_n16+Mn_n16+Fe_n16+Al_VI+Zn_n16+Cr_n16+V_n16+Ti_n16)/(Si_n16+Al_IV) < 4/6 + 4/6*t)
                if(Mg_n16 >= Fe_n16 && Al_IV/Si_n16 < 0.1 && (Na_n16+K_n16+Ca_n16)/(Si_n16+Al_IV) < 0.025 && Al_VI < Mg_n16 && Al_VI <= Fe_n16)
                    T.group2(i) = {'Phyllosilicate'};
                    T.group3(i) = {'Palygorskite-Sepiolite'};
                    T.group4(i) = {''};
%                     T.species(i) = {'sepiolite'};
                    OH = 2 - (F_n16+Cl_n16);
                    if(OH<=0), OH=0; end
                    [sorted_formula1, total1] = Formula_Output({Mg_n16, Fe_n16, Mn_n16, Ca_n16, Al_VI, Cr_n16, Ti_n16},{'Mg','Fe','Mn','Ca','AlVI','Cr','Ti'});
                    [sorted_formula2, total2] = Formula_Output({Si_n16, Al_IV},{'Si','AlIV'});
                    [sorted_formula3, total3] = Formula_Output({OH, F_n16, Cl_n16},{'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O15', sorted_formula3, total3, '·6H2O')};
                end
            end

        end %end of Catch

        %% 7.7.2 Pyroxenes
        if(strcmp(T.group3(i), 'Pyroxene') && S_n <= 0.05 && P_n <= 0.05)
            Renorm_factor = 4;
            Mg_n6 = T.Mg_Multi(i) / Renorm_factor;
            Mn_n6 = T.Mn_Multi(i) / Renorm_factor;
            Fe_n6 = T.Fe_Multi(i) / Renorm_factor;
            Ca_n6 = T.Ca_Multi(i) / Renorm_factor;
            Si_n6 = T.Si_Multi(i) / Renorm_factor;
            Ti_n6 = T.Ti_Multi(i) / Renorm_factor;
            Cr_n6 = T.Cr_Multi(i) / Renorm_factor;
            Na_n6 = T.Na_Multi(i) / Renorm_factor;
            Al_n6 = T.Al_Multi(i) / Renorm_factor;
            ABCD_n6 = T.ABCD(i) / Renorm_factor;
            Al_IV = 2 - Si_n6;
            if(Al_IV < 0), Al_IV = 0; end
            Al_VI = Al_n6 - Al_IV;
            if(Al_VI < 0), Al_VI = 0; end
            delta = (Na_n6 + Ca_n6 + Mg_n6 + Fe_n6 + Cr_n6 + Ti_n6 + Mn_n6 + Al_VI) - 2;
            if(delta > 0), Fe3 = Fe_n6 - delta; end
            if(delta <= 0), Fe3 = 0; end
            Fe2 = Fe_n6 - Fe3;
            if(Fe2 < 0), Fe2 = 0; end

            %% 7.7.2.1. Clinopyroxenes & Na-rich Pyroxenes
            if((Ca_n / Si_n >= 0.04 && Ca_n / Si_n < 1/2 + 1/2*t) || (Na_n/Si_n > 1/2 - 2*(1/2)*t && Na_n/Si_n < 1/2 + (1/2)*t))                   
                T.group4(i) = {'Clinopyroxene'};
                Wo = Ca_n6 / (Ca_n6 + Mg_n6 + Fe_n6) * 100;
                En = Mg_n6 / (Ca_n6 + Mg_n6 + Fe_n6) * 100;
                Fs = Fe_n6 / (Ca_n6 + Mg_n6 + Fe_n6) * 100;

                Aeg = ((Fe3+Na_n6) / ((Fe3+Na_n6) + (Al_VI+Na_n6) + (Ca_n6+Fe2+Mg_n6))) * 100;
                Jd = ((Al_VI+Na_n6) / ((Fe3+Na_n6) + (Al_VI+Na_n6) + (Ca_n6+Fe2+Mg_n6))) * 100;
                Di = ((Ca_n6+Fe2+Mg_n6) / ((Fe3+Na_n6) + (Al_VI+Na_n6) + (Ca_n6+Fe2+Mg_n6))) * 100;

                T.Wo(i) = round(Wo,2);
                T.En(i) = round(En,2);
                T.Fs(i) = round(Fs,2);
                T.Aeg(i) = round(Aeg,2);
                T.Jd(i) = round(Jd,2);
                T.Di(i) = round(Di,2);

                if(Na_n > Ca_n)
                    T.group4(i) = {'Na-rich Pyroxene'};

                    if(Fe3 > Al_VI && round(Aeg/Jd) >= 1 && round(Di) <= 20 && Mg_n < 1 && Ca_n < 1 && Al_n < 1)
                        T.species(i) = {'aegirine'};
                        [sorted_formula1, total1] = Formula_Output({Na_n6, Ca_n6}, {'Na', 'Ca'});
                        [sorted_formula2, total2] = Formula_Output({Fe3,Fe2,Mg_n6,Mn_n6,Al_VI,Cr_n6,Ti_n6},{'Fe3+','Fe2+', 'Mg','Mn','AlVI','Cr','Ti'});
                        [sorted_formula3, total3] = Formula_Output({Si_n6,Al_IV}, {'Si','AlIV'});
                        T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, '[', sorted_formula3, total3, 'O6]')};
                    end

                    if(Al_VI > Fe3 && T.SiO2(i) > 58-1*t && T.SiO2(i) <= 61.5+0.5*t && round(Aeg/Jd) <= 1 && round(Di) <= 20)
                        T.species(i) = {'jadeite'};
                        [sorted_formula1, total1] = Formula_Output({Na_n6, Ca_n6, Mn_n6}, {'Na', 'Ca', 'Mn'});
                        [sorted_formula2, total2] = Formula_Output({Fe3,Al_VI,Cr_n6,Ti_n6},{'Fe3+', 'AlVI','Cr','Ti'});
                        [sorted_formula3, total3] = Formula_Output({Si_n6,Al_IV}, {'Si','AlIV'});
                        T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, '[', sorted_formula3, total3, 'O6]')};
                    end
                end

                if(Na_n >= 0 && Na_n < 0.4)
                    if(Mg_n > Fe_n && Ca_n > Fe_n && Al_IV/Si_n < 1/2 + (1/2)*t && round(Wo) <= 45 && round(Wo) >= 25 && round(En/Fs) >= 1)
                        T.species(i) = {'augite'};
                        [sorted_formula1, total1] = Formula_Output({Na_n6,Ca_n6,Mg_n6,Fe2,Mn_n6,Al_VI,Fe3,Cr_n6,Ti_n6},{'Na','Ca','Mg','Fe2+','Mn','AlVI','Fe3+','Cr','Ti'});
                        [sorted_formula2, total2] = Formula_Output({Si_n6,Al_IV},{'Si','AlIV'});
                        T.formula(i) = {strcat(sorted_formula1,total1,'[',sorted_formula2,total2,'O6]')};
                    end
                end

                if(Na_n >= 0 && Na_n < 0.4)
                    if(Mg_n < Fe_n && Ca_n > Mg_n && Al_IV/Si_n < 1/2 + (1/2)*t && round(Wo) <= 45 && round(Wo) >= 25 && round(En/Fs) <= 1)
                        T.species(i) = {'augite'}; % was ferroaugite
                        [sorted_formula1, total1] = Formula_Output({Na_n6,Ca_n6,Mg_n6,Fe2,Mn_n6,Al_VI,Fe3,Cr_n6,Ti_n6},{'Na','Ca','Mg','Fe2+','Mn','AlVI','Fe3+','Cr','Ti'});
                        [sorted_formula2, total2] = Formula_Output({Si_n6,Al_IV},{'Si','AlIV'});
                        T.formula(i) = {strcat(sorted_formula1,total1,'[',sorted_formula2,total2,'O6]')};
                    end
                end

                if(Ca_n > Na_n)

                    if(round((Ca_n6+Na_n6)/(Mg_n6+Fe_n6+Al_VI+Ti_n6+Cr_n6),2) > 1 - 2*1*t &&  round((Ca_n6+Na_n6)/(Mg_n6+Fe_n6+Al_VI+Ti_n6+Cr_n6),2) < 1 + 2*1*t)

                        if(Mg_n > Fe_n && Al_IV/Si_n < 1/2 + (1/2)*t && round(Wo) >= 45 && round(Wo) <= 50 && round(En/Fs) >= 1)
                            T.species(i) = {'diopside'};
                            [sorted_formula1, total1] = Formula_Output({Ca_n6,Na_n6,Mn_n6}, {'Ca','Na','Mn'});
                            [sorted_formula2, total2] = Formula_Output({Mg_n6,Fe_n6,Al_VI,Cr_n6,Ti_n6},{'Mg','Fe','AlVI','Cr','Ti'});
                            [sorted_formula3, total3] = Formula_Output({Si_n6, Al_IV}, {'Si','AlIV'});

                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3,total3,'O6')};
                        end

                        if(Fe_n > Mg_n && Al_IV/Si_n < 1/2 + (1/2)*t && round(Wo) >= 45 && round(Wo) <= 50 && round(En/Fs) <= 1)
                            T.species(i) = {'hedenbergite'};

                            [sorted_formula1, total1] = Formula_Output({Ca_n6, Na_n6, Mn_n6},{'Ca','Na','Mn'});
                            [sorted_formula2, total2] = Formula_Output({Fe_n6, Mg_n6, Al_VI, Cr_n6, Ti_n6},{'Fe','Mg','AlVI','Cr','Ti'});
                            [sorted_formula3, total3] = Formula_Output({Si_n6, Al_IV},{'Si','AlIV'});

                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3,'O6')};
                        end

                    end
                    
                    if(Ca_n + Na_n > 4-4*t && Ca_n + Na_n < 4+4*t && Fe3 > 0 && Fe_n > Al_n && Fe3 > Al_VI && T.Na2O(i) > 2 && Fe3 > Fe2 && round(Aeg/Jd) >= 1 && round(Di,2) <= 80 && round(Di,2) > 20)
                        T.group4(i) = {'Na-rich Pyroxene'};
                        T.species(i) = {'aegirine-augite'};
                        [sorted_formula1, total1] = Formula_Output({Na_n6, Ca_n6}, {'Na', 'Ca'});
                        [sorted_formula2, total2] = Formula_Output({Fe3,Fe2,Mg_n6,Mn_n6,Al_VI,Cr_n6,Ti_n6},{'Fe3+','Fe2+', 'Mg','Mn','AlVI','Cr','Ti'});
                        [sorted_formula3, total3] = Formula_Output({Si_n6,Al_IV}, {'Si','AlIV'});
                        T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, '[', sorted_formula3, total3, 'O6]')};
                    end

                end

                if(T.Na2O(i) >= 4 && T.Na2O(i) < 8 && Ca_n > Na_n && Al_VI > Fe3 && round(Aeg/Jd) <= 1 && round(Di,2) <= 80 && round(Di,2) > 20)
                    if(Na_n6+Ca_n6 > 1-1*t && Na_n6+Ca_n6 <= 1+1*t && Mg_n6+Fe_n6+Al_VI+Ti_n6+Cr_n6 > 1-1*t && Mg_n6+Fe_n6+Al_VI+Ti_n6+Cr_n6 <= 1+1*t && Si_n6+Al_IV > 2-1*t && Si_n6+Al_IV < 2+1*t)
                        T.group4(i) = {'Na-rich Pyroxene'};
                        T.species(i) = {'omphacite'};
                        Al_IV = 2 - Si_n/4;
                        if(Al_IV<=0), Al_IV = 0; end
                        Al_VI = Al_n/4 - Al_IV;
                        if(Al_VI<=0), Al_VI = 0; end
                        [sorted_formula1, total1] = Formula_Output({Ca_n/4, Na_n/4},{'Ca','Na'});
                        [sorted_formula2, total2] = Formula_Output({Mg_n/4, Fe_n/4, Mn_n/4, Al_VI},{'Mg','Fe','Mn','AlVI'});
                        [sorted_formula3, total3] = Formula_Output({Si_n/4, Al_IV},{'Si','AlIV'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3,'O6')};
                    end
                end

                if(T.CaO(i) > 2 - 2*t && T.CaO(i) < 8 + 8*t && Ca_n > Na_n && Mg_n > Ca_n && Fe_n > Ca_n && T.Al2O3(i) < 8)
                    if(Mg_n + Fe_n + Ca_n > 8-8*t && Mg_n + Fe_n + Ca_n < 8+8*t)
                        T.group4(i) = {'Clinopyroxene'};
                        T.species(i) = {'pigeonite'};
                        Al_IV = 2 - Si_n/4;
                        if(Al_IV<=0), Al_IV = 0; end
                        Al_VI = Al_n/4 - Al_IV;
                        if(Al_VI<=0), Al_VI = 0; end
                        [sorted_formula1, total1] = Formula_Output({Mg_n/4, Fe_n/4, Mn_n/4, Al_VI, Ca_n/4, Na_n/4},{'Mg','Fe','Mn','AlVI','Ca','Na'});
                        [sorted_formula2, total2] = Formula_Output({Si_n/4, Al_IV},{'Si','AlIV'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O6')};
                    end
                end
            
            end %end of clinopyroxene

            %% 7.7.2.2. Orthopyroxenes
            if(T.CaO(i) < 2.2 && Na_n6 < 0.1 && (Mg_n6+Fe_n6+Mn_n6+Ca_n6+Al_n6)/Si_n6 > 1-1*t && round( ((Mg_n6+Fe_n6+Mn_n6+Ca_n6+Al_n6)/Si_n6),1) <= 1+1*t && (Mg_n6+Fe_n6+Mn_n6)>(Ca_n6+Al_n6))
                T.group4(i) = {'Orthopyroxene'};
                T.species(i) = {''};
                T.formula(i) = {''};
                Wo = Ca_n6 / (Ca_n6 + Mg_n6 + Fe_n6) * 100;
                En = Mg_n6 / (Ca_n6 + Mg_n6 + Fe_n6) * 100;
                Fs = Fe_n6 / (Ca_n6 + Mg_n6 + Fe_n6) * 100;
                T.Wo(i) = round(Wo,2);
                T.En(i) = round(En,2);
                T.Fs(i) = round(Fs,2);

                if(round(En,2) > 10 && round(En,2) < 90)
                    T.species(i) = {'orthopyroxene'};
                    [sorted_formula1, total1] = Formula_Output({Mg_n6, Fe_n6, Mn_n6, Ca_n6, Al_VI, Cr_n6, Ti_n6}, {'Mg','Fe','Mn','Ca','AlVI','Cr','Ti'});
                    [sorted_formula2, total2] = Formula_Output({Si_n6, Al_IV},{'Si','AlIV'});
                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, 'O6')};
                end

                if(round(En,2) >= 90)
                    T.species(i) = {'enstatite'};
                    [sorted_formula1, total1] = Formula_Output({Mg_n6, Fe_n6, Mn_n6, Ca_n6, Al_VI, Cr_n6, Ti_n6}, {'Mg','Fe','Mn','Ca','AlVI','Cr','Ti'});
                    [sorted_formula2, total2] = Formula_Output({Si_n6, Al_IV},{'Si','AlIV'});
                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, 'O6')};
                end

                if(round(Fs,2) > 50)
                    T.species(i) = {'ferrosilite'};
                    [sorted_formula1, total1] = Formula_Output({Fe_n6, Mg_n6, Mn_n6, Ca_n6, Na_n6, Cr_n6, Al_VI, Ti_n6},{'Fe2+','Mg','Mn','Ca','Na','Cr','AlVI','Ti'});
                    [sorted_formula2, total2] = Formula_Output({Si_n6, Al_IV},{'Si','AlIV'});
                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2,'O6')};
                end

                if(round(En,2) > 50 && round(En,2) < 90)
                    T.species(i) = {'hypersthene'};
                    [sorted_formula1, total1] = Formula_Output({Mg_n6, Fe_n6, Ca_n6, Mn_n6, Ti_n6 Al_VI, Cr_n6}, {'Mg','Fe','Ca','Mn','Ti','AlVI','Cr'});
                    [sorted_formula2, total2] = Formula_Output({Si_n6, Al_IV},{'Si','AlIV'});
                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, 'O6')};
                end

                if(round(Wo,2) > 90)
                    T.group4(i) = {'Pyroxenoid'};
                    T.species(i) = {'wollastonite'};
                    renorm = 3/24;
                    [sorted_formula1, total1] = Formula_Output({Ca_n*renorm, Mg_n*renrom, Mn_n*renorm, Na_n*renorm});
                    [sorted_formula2, total2] = Formula_Output({Si_n*renorm, Al_n*renrom});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O3')};

                    if(Al_n >= Si_n)

                        T.group1(i) = {'Silicate'};
                        T.group2(i) = {'Chain Silicate'};
                        T.group3(i) = {'Pyroxene'};
                        T.group4(i) = {'Clinopyroxne'};
                        T.species(i) = {'Ca-Tschermaks/kushiroite'};

                        Al_IV = 3 - Si_n/2;
                        if(Al_IV < 0), Al_IV = 0; end
                        
                        Al_VI = Al_n/2 - Al_IV;
                        if(Al_VI<0), Al_VI = 0; end

                        Si_a = 1;
                        Si_b = Si_n/2 - Si_a;
                        if(Si_b < 0), Si_b = 0; end

                        Al_IV = 2 - Si_b;
                        if(Al_IV<0), Al_IV = 0; end

                        [sorted_formula1, total1] = Formula_Output({Ca_n/2, Na_n/2}, {'Ca','Na'});
                        [sorted_formula2, total2] = Formula_Output({Al_VI, Fe_n/2, Mg_n/2, Mn_n/2, Cr_n/2}, {'AlVI','Fe3+','Mg','Mn','Cr'});
                        [sorted_formula3, total3] = Formula_Output({Si_b, Al_IV}, {'Si','AlIV'});
                        [sorted_formula4, total4] = Formula_Output({Si_a}, {'Si'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, '[',sorted_formula3,total3,'O7][', sorted_formula4, total4,'O4]O(OH)')};

                    end
                end

                if((Mg_n+Fe_n)/(Si_n+Al_n) <= 7/8 && (Mg_n+Fe_n)/(Si_n+Al_n) > 7/8 - (7/8)*t && Al_n/Si_n < 0.15 && Mg_n/(Mg_n+Fe_n) >= 0.5)
                    T.group3(i) = {''};
                    T.group4(i) = {''};
                    T.species(i) = {''}; 
                end

            end %end of orthopyroxene

            %% 7.7.2.3. Pyroxenoids
            Wo = Ca_n6 / (Ca_n6 + Mg_n6 + Fe_n6) * 100;
            En = Mg_n6 / (Ca_n6 + Mg_n6 + Fe_n6) * 100;
            Fs = Fe_n6 / (Ca_n6 + Mg_n6 + Fe_n6) * 100;
            T.Wo(i) = round(Wo,2);
            T.En(i) = round(En,2);
            T.Fs(i) = round(Fs,2);

            if(Wo >= 50 && Ca_n/Si_n > 0.7 && Ca_n/Si_n < 1.1)
                T.group4(i) = {'Pyroxenoid'};
                T.species(i) = {'wollastonite'};
                renorm = 3/24;
                [sorted_formula1, total1] = Formula_Output({Ca_n*renorm, Mg_n*renorm, Mn_n*renorm, Na_n*renorm}, {'Ca','Mg','Mn','Na'});
                [sorted_formula2, total2] = Formula_Output({Si_n*renorm, Al_n*renorm}, {'Si','Al'});
                T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O3')};

                if(Al_n >= Si_n)

                    T.group1(i) = {'Silicate'};
                    T.group2(i) = {'Chain Silicate'};
                    T.group3(i) = {'Pyroxene'};
                    T.group4(i) = {'Clinopyroxne'};
                    T.species(i) = {'Ca-Tschermaks/kushiroite'};
    
                    Al_IV = 3 - Si_n/2;
                    if(Al_IV < 0), Al_IV = 0; end
                    
                    Al_VI = Al_n/2 - Al_IV;
                    if(Al_VI<0), Al_VI = 0; end
    
                    Si_a = 1;
                    Si_b = Si_n/2 - Si_a;
                    if(Si_b < 0), Si_b = 0; end
    
                    Al_IV = 2 - Si_b;
                    if(Al_IV<0), Al_IV = 0; end
    
                    [sorted_formula1, total1] = Formula_Output({Ca_n/2, Na_n/2}, {'Ca','Na'});
                    [sorted_formula2, total2] = Formula_Output({Al_VI, Fe_n/2, Mg_n/2, Mn_n/2, Cr_n/2}, {'AlVI','Fe3+','Mg','Mn','Cr'});
                    [sorted_formula3, total3] = Formula_Output({Si_b, Al_IV}, {'Si','AlIV'});
                    [sorted_formula4, total4] = Formula_Output({Si_a}, {'Si'});
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, '[',sorted_formula3,total3,'O7][', sorted_formula4, total4,'O4]O(OH)')};
                end
            end

            if(ABCD_n6 / 3 > (2/3)-(2/3)*t && ABCD_n6 / 3 < (2/3)+(2/3)*t)
                if(Mn_n / ABCD_n > 0.9)
                    if(Mn_n / Si_n >= 1-1*2*t && Mn_n / Si_n < 1+1*t)
                        T.group4(i) = {'Pyroxenoid'};
                        T.species(i) = {'rhodonite'};
                        renorm = 15/24;
                        Al_IV = 5 - Si_n*renorm;
                        if(Al_IV < 0), Al_IV = 0; end
                        Al_VI = Al_n*renorm - Al_IV;
                        if(Al_VI < 0), Al_VI = 0; end
                        Mny = 1 - (Al_VI+Ti_n*renorm+Cr_n*renorm);
                        if(Mny < 0), Mny = 0; end
                        Mnx = Mn_n*renorm - Mny;
                        if(Mnx < 0), Mnx = 0; end
                        [sorted_formula1, total1] = Formula_Output({Na_n*renorm, Ca_n*renorm},{'Na','Ca'});
                        [sorted_formula2, total2] = Formula_Output({Mnx,Fe_n*renorm, Mg_n*renorm},{'Mn','Fe','Mg'});
                        [sorted_formula3, total3] = Formula_Output({Mny,Al_VI,Ti_n*renorm,Cr_n*renorm},{'Mn','AlVI','Ti','Cr'});
                        [sorted_formula4, total4] = Formula_Output({Si_n*renorm,Al_IV},{'Si','AlIV'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3,total3,sorted_formula4,total4,'O15')};
                    end
                end
            end

        end %end of pyroxenes

        % Catch
        if(strcmp(T.group2(i), 'Inosilicate') && ~strcmp(T.group3(i), 'Amphibole') && ~strcmp(T.group3(i), 'Pyroxene'))
            if(Na_n >= Ca_n && Na_n >= K_n && S_n <= 0.05 && P_n <= 0.05)
                Renorm_factor = 24 / 23;
                Mg_n23 = T.Mg_Multi(i) / Renorm_factor;
                Mn_n23 = T.Mn_Multi(i) / Renorm_factor;
                Fe_n23 = T.Fe_Multi(i) / Renorm_factor;
                Si_n23 = T.Si_Multi(i) / Renorm_factor;
                Ti_n23 = T.Ti_Multi(i) / Renorm_factor;
                Cr_n23 = T.Cr_Multi(i) / Renorm_factor;
                Al_n23 = T.Al_Multi(i) / Renorm_factor;
                Al_IV = 8 - Si_n23;
                if(Al_IV < 0), Al_IV = 0; end
                Al_VI = Al_n23 - Al_IV;
                if(Al_VI < 0), Al_VI = 0; end

                if((Mg_n23 + Fe_n23 + Mn_n23 + Cr_n23 + Al_VI + Ti_n23)/(Si_n23 + Al_IV) > (5/8)-2*(5/8)*t && (Mg_n23 + Fe_n23 + Mn_n23 + Cr_n23 + Al_VI + Ti_n23)/(Si_n23 + Al_IV) < (5/8)+2*(5/8)*t)
                    if((K_n23+Na_n23+Ca_n23)/(Si_n23+Al_IV) > (3/8)-(3/8)*2*t && (K_n23+Na_n23+Ca_n23)/(Si_n23+Al_IV) < (3/8)+(3/8)*2*t)
                        if((K_n23+Ca_n23+Na_n23) > 3 - 3*2*t && (K_n23+Ca_n23+Na_n23) < 3 + 3*2*t  && (Al_IV/Si_n23) < (1/8) + (1/8)*2*t)
                            if(T.SiO2(i) >= 55)
                                T.group3(i) = {'Amphibole'};
%                                 T.group4(i) = {'Sodic-Calcic Amphibole'};
%                                 T.species(i) = {'potassic-richterite'};
                                Nax = 1 - (K_n23+Na_n23);
                                if(Nax <= 0), Nax = 0; end
                                Nay = Na_n23 - Nax;
                                if(Nay <= 0), Nay = 0; end
                                vac = 2 - (Nax + Ca_n23);
                                if(vac <= 0), vac = 0; end
                                OH = 2 - (F_n23+Cl_n23);
                                if(OH <= 0), OH = 0; end
                                [sorted_formula1, total1] = Formula_Output({Nax, K_n23, vac}, {'Na','K','□'});
                                [sorted_formula2, total2] = Formula_Output({Nay, Ca_n23, Mn_n23}, {'Na', 'Ca','Mn'});
                                [sorted_formula3, total3] = Formula_Output({Mg_n23, Fe_n23, Al_VI, Ti_n23, Cr_n23}, {'Mg','Fe','AlVI','Ti','Cr'});
                                [sorted_formula4, total4] = Formula_Output({Si_n23, Al_IV}, {'Si','AlIV'});
                                [sorted_formula5, total5] = Formula_Output({OH, F_n23, Cl_n23}, {'OH','F','Cl'});
                                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3, sorted_formula4, total4, 'O22', sorted_formula5, total5)};
                            end
                        end
                    end
                end

            end
        end

        %% 7.8 Phyllosilicates

        %% 7.81 Ca Mica & Prehnite
        if(strcmp(T.group1(i), 'Silicate') && Ca_n > Na_n && Ca_n > K_n)

            Renorm_factor = 24 / 11;
            Mg_n11 = T.Mg_Multi(i) / Renorm_factor;
            Mn_n11 = T.Mn_Multi(i) / Renorm_factor;
            Fe_n11 = T.Fe_Multi(i) / Renorm_factor;
            Ca_n11 = T.Ca_Multi(i) / Renorm_factor;
            Si_n11 = T.Si_Multi(i) / Renorm_factor;
            Ni_n11 = T.Ni_Multi(i) / Renorm_factor;
            Ti_n11 = T.Ti_Multi(i) / Renorm_factor;
            Cr_n11 = T.Cr_Multi(i) / Renorm_factor;
            Na_n11 = T.Na_Multi(i) / Renorm_factor;
            Al_n11 = T.Al_Multi(i) / Renorm_factor;
            K_n11 = T.K_Multi(i) / Renorm_factor;
            F_n11 = T.F_Multi(i) / Renorm_factor;
            Cl_n11 = T.Cl_Multi(i) / Renorm_factor;

            Al_IV = 4- Si_n11;
            if(Al_IV<0), Al_IV=0; end

            Al_VI = Al_n11 - Al_IV;
            if(Al_VI<0), Al_VI=0; end

            OH = 2 - (F_n11 + Cl_n11);
            if(OH<0), OH = 0; end

            if((Ca_n11+Na_n11)/(Al_n11+Fe_n11+Mg_n11+Mn_n11) >= 1/4 - (1/4)*t && (Ca_n11+Na_n11)/(Al_n11+Fe_n11+Mg_n11+Mn_n11) < 1/4 + (1/4)*t && (Al_VI+Fe_n11+Mg_n11+Mn_n11)/(Si_n11+Al_IV) >= 2/4 - (2/4)*t && (Ca_n11+Na_n11)/(Al_VI+Fe_n11+Mg_n11+Mn_n11) < 2/4 + (2/4)*t)
                if(Al_VI > Fe_n11 && Al_VI > Mg_n11 && Al_VI > Mg_n11)
                    T.group2(i) = {'Phyllosilicate'};
                    T.group3(i) = {'Mica'};
                    T.group4(i) = {'Ca-Mica'};
                    T.species(i) = {'margarite'};

                    [sorted_formula1, total1] = Formula_Output({Ca_n11, Na_n11, K_n11}, {'Ca','Na','K'});
                    [sorted_formula2, total2] = Formula_Output({Al_VI, Fe_n11, Cr_n11, Ti_n11, Mg_n11, Mn_n11}, {'AlVI','Fe','Cr','Ti','Mg','Mn'});
                    [sorted_formula3, total3] = Formula_Output({Si_n11, Al_IV},{'Si','AlIV'});
                    [sorted_formula4, total4] = Formula_Output({OH, F_n11, Cl_n11}, {'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, 'O10', sorted_formula4, total4)};
                end
            end

            if((Ca_n11+Na_n11)/(Al_n11+Fe_n11+Mg_n11+Mn_n11) >= 1 - 1*t && (Ca_n11+Na_n11)/(Al_n11+Fe_n11+Mg_n11+Mn_n11) < 1 + 1*t && (Al_VI+Fe_n11+Mg_n11+Mn_n11)/(Si_n11+Al_IV) >= 1/4 - (1/4)*t && (Al_VI+Fe_n11+Mg_n11+Mn_n11)/(Si_n11+Al_IV) < 1.5/4 + (1.5/4)*t)
                if(Al_VI > Fe_n11 && Al_VI > Mg_n11 && Al_VI > Mg_n11)
                    T.group2(i) = {'Inosilicate'};
                    T.group3(i) = {'Prehnite'};
                    T.group4(i) = {''};
                    T.species(i) = {'prehnite'};
                end
                if(Fe_n11 > Mg_n11 && Fe_n11 > Mn_n11 && Fe_n11 > Al_VI)
                    if((Fe_n11+Al_VI)/(Ca_n11+Na_n11) >= 1/2 - (1/2)*t && (Fe_n11+Al_VI)/(Ca_n11+Na_n11) < 1/2 + (1/2)*2*t)
                        T.group2(i) = {'Inosilicate'};
                        T.group3(i) = {'Prehnite'};
                        T.group4(i) = {''};
                        T.species(i) = {'ferriprehnite'};
                    end
                end

                [sorted_formula1, total1] = Formula_Output({Ca_n11, Na_n11, K_n11}, {'Ca','Na','K'});
                [sorted_formula2, total2] = Formula_Output({Al_VI, Fe_n11, Cr_n11, Ti_n11, Mg_n11, Mn_n11}, {'AlVI','Fe','Cr','Ti','Mg','Mn'});
                [sorted_formula3, total3] = Formula_Output({Si_n11, Al_IV},{'Si','AlIV'});
                [sorted_formula4, total4] = Formula_Output({OH, F_n11, Cl_n11}, {'OH','F','Cl'});
                T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, 'O10', sorted_formula4, total4)};

            end
        end

        %% 7.8.1 Di & Tri Micas Assignment
        if(strcmp(T.group1(i), 'Silicate') && ~strcmp(T.group3(i), 'Amphibole') && ~strcmp(T.group3(i), 'Pyroxene') && ~strcmp(T.group3(i), 'Prehnite') && ~strcmp(T.group2(i), 'Tectosilicate') && ~strcmp(T.group4(i), 'Ca-Mica') && S_n <= 0.05 && P_n <= 0.05)
            if((T.K2O(i) > 6 && T.K2O(i) < 13 && round(Al_n/Si_n,2) <= 1.51 && K_n/(K_n+Na_n+Ca_n) > 0.5)...
                    || (T.Na2O(i) > 6*(100/95) && round(Al_n/Si_n,1) <= 1 && Na_n/(K_n+Na_n+Ca_n) > 0.5)...
                    || (T.Al2O3(i) > 50 && Ca_n/(Si_n+Al_n+Fe_n+Mg_n) <= (1/6)+(1/6)*t && Ca_n/(K_n+Na_n+Ca_n) > 0.5))

                T.group2(i) = {'Phyllosilicate'};
                T.group3(i) = {'Mica'};
                T.group4(i) = {''};
                T.species(i) = {''};
                T.formula(i) = {''};

                if(Ca_n/Si_n > (1/8)-2*(1/8)*t && Ca_n/Si_n < (1/8)+2*(1/8)*t && Al_n23 < 0.5)
                    if((Mg_n+Fe_n)/Si_n > (5/8)-(5/8)*t && (Mg_n+Fe_n)/Si_n < (5/8)+(5/8)*t && Mg_n > Fe_n && Na_n > Ca_n)
                        T.group2(i) = {'Inosilicate'};
                        T.group3(i) = {'Amphibole'};
                        T.group4(i) = {'Sodic-Calcic Amphibole'};
                        T.species(i) = {'richterite'};
                        Al_IV = 8 - Si_n23;
                        if(Al_IV < 0), Al_IV = 0; end
                        Al_VI = Al_n23 - Al_IV;
                        if(Al_VI < 0), Al_VI = 0; end
                        Nax = 1 - (K_n23+Na_n23);
                        if(Nax <= 0), Nax = 0; end
                        Nay = Na_n23 - Nax;
                        if(Nay <= 0), Nay = 0; end
                        vac = 2 - (Nax + Ca_n23);
                        if(vac <= 0), vac = 0; end
                        OH = 2 - (F_n23+Cl_n23);
                        if(OH <= 0), OH = 0; end
                        [sorted_formula1, total1] = Formula_Output({Nax, K_n23, vac}, {'Na','K','□'});
                        [sorted_formula2, total2] = Formula_Output({Nay, Ca_n23, Mn_n23}, {'Na', 'Ca','Mn'});
                        [sorted_formula3, total3] = Formula_Output({Mg_n23, Fe_n23, Al_VI, Ti_n23, Cr_n23}, {'Mg','Fe','AlVI','Ti','Cr'});
                        [sorted_formula4, total4] = Formula_Output({Si_n23, Al_IV}, {'Si','AlIV'});
                        [sorted_formula5, total5] = Formula_Output({OH, F_n23, Cl_n23}, {'OH','F','Cl'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3, sorted_formula4, total4, 'O22', sorted_formula5, total5)};
                    end
                end

                Renorm_factor = 24 / 11;
                Mg_n11 = T.Mg_Multi(i) / Renorm_factor;
                Mn_n11 = T.Mn_Multi(i) / Renorm_factor;
                Fe_n11 = T.Fe_Multi(i) / Renorm_factor;
                Ca_n11 = T.Ca_Multi(i) / Renorm_factor;
                Si_n11 = T.Si_Multi(i) / Renorm_factor;
                Ni_n11 = T.Ni_Multi(i) / Renorm_factor;
                Ti_n11 = T.Ti_Multi(i) / Renorm_factor;
                Cr_n11 = T.Cr_Multi(i) / Renorm_factor;
                Na_n11 = T.Na_Multi(i) / Renorm_factor;
                Al_n11 = T.Al_Multi(i) / Renorm_factor;
                K_n11 = T.K_Multi(i) / Renorm_factor;
                F_n11 = T.F_Multi(i) / Renorm_factor;
                Cl_n11 = T.Cl_Multi(i) / Renorm_factor;
                V_n11 = T.V_Multi(i) / Renorm_factor;
                Ba_n11 = T.Ba_Multi(i) / Renorm_factor;

                if(K_n > (Na_n + Ca_n) && (Mg_n + Fe_n) > Al_n)
                    T.group4(i) = {'Trioctahedral Mica'};
                end

                if(K_n > (Na_n + Ca_n) || Na_n > (K_n + Ca_n))
                    if(round((Mg_n + Fe_n),2) <= round(Al_n,2) || Al_n/Si_n < 0.2)
                        T.group4(i) = {'Dioctahedral Mica'};
                    end
                end

            end
        end

        %% 7.8.1.1 Trioctahedral Micas
        if(strcmp(T.group4(i), 'Trioctahedral Mica') && S_n <= 0.05 && P_n <= 0.05)
            Al_IV = 4 - Si_n11;
            if(Al_IV < 0), Al_IV = 0; end
            Al_VI = Al_n11 - Al_IV;
            if(Al_VI < 0), Al_VI = 0; end

            if((K_n + Na_n) / (Mg_n + Fe_n) > 0.25-0.25*t && (K_n + Na_n) / (Mg_n + Fe_n) <= 1+1*t)
                if( ((K_n + Na_n) / (Si_n+Al_IV) > 0.25-0.25*t && (K_n + Na_n) / (Si_n+Al_IV) < 0.5+0.5*t) || (Al_n / Si_n < 0.5+0.5*t) )
                    T.group4(i) = {'Trioctahedral Mica (biotite)'};

                    if(K_n11 > Na_n11)

                        if((Na_n11+K_n11) > 0.9-0.9*t && (Na_n11+K_n11) < 1+1*t && Si_n11/Al_IV > 2.5-2.5*t && Si_n11/Al_IV < 3+3*t)
                            if(Fe_n11 + Mg_n11 + Mn_n11 + Al_VI + Cr_n11 + Ti_n11 > 3-3*t && Fe_n11 + Mg_n11 + Mn_n11 + Al_VI + Cr_n11 + Ti_n11 < 3+3*t)

                                if(Mg_n11 >= Fe_n11)
                                    Fe2 = Fe_n11;
                                    T.species(i) = {'biotite (phlogopite)'};
                                    OH = 2 - (Cl_n11 + F_n11);
                                    if(OH<=0), OH=0; end
                                    vac = 1 - (Na_n11 + K_n11);
                                    if(vac <= 0), vac = 0; end
                                    [sorted_formula1, total1] = Formula_Output({Na_n11, K_n11,vac},{'Na','K','□'});
                                    [sorted_formula2, total2] = Formula_Output({Fe2,Mg_n11,Mn_n11,Al_VI,Cr_n11,Ti_n11},{'Fe2+','Mg','Mn','AlVI','Cr','Ti'});
                                    [sorted_formula3, total3] = Formula_Output({Al_IV,Si_n11},{'AlIV','Si'});
                                    [sorted_formula4, total4] = Formula_Output({OH, Cl_n11, F_n11}, {'OH','Cl','F'});
                                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3,'O10',sorted_formula4, total4)};
                                end
                                if(Fe_n11 > Mg_n11)
                                    Fe2 = Fe_n11;
                                    T.species(i) = {'biotite (annite)'};
                                    OH = 2 - (Cl_n11 + F_n11);
                                    if(OH<=0), OH=0; end
                                    [sorted_formula1, total1] = Formula_Output({K_n11,Na_n11},{'K','Na'});
                                    [sorted_formula2, total2] = Formula_Output({Fe2,Mg_n11,Mn_n11,Al_VI,Cr_n11,Ti_n11},{'Fe2+','Mg','Mn','AlVI','Cr','Ti'});
                                    [sorted_formula3, total3] = Formula_Output({Al_IV,Si_n11},{'AlIV','Si'});
                                    [sorted_formula4, total4] = Formula_Output({OH, Cl_n11, F_n11}, {'OH','Cl','F'});
                                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'(',sorted_formula3,total3,'O10)',sorted_formula4, total4)};
                                end

                            end
                        end
                        
                        Fe3 = 1 - (Al_VI + Ti_n11 + Cr_n11);
                        if(Fe3 < 0), Fe3 = 0; end
                        Fe2 = Fe_n11 - Fe3;
                        if(Fe2 < 0), Fe2 = 0; end

                        if(Fe_n11 > Mg_n11)
                            if((Fe2 + Mg_n11 + Mn_n11)/(Al_n11 + Fe3 + Cr_n11 + Ti_n11) >= (1.5/3)-(1.5/3)*2*t && (Fe2 + Mg_n11 + Mn_n11)/(Al_n11 + Fe3 + Cr_n11 + Ti_n11) <= (2/3)+(2/3)*3*t)
                                if(Si_n11 / (Al_n11 + Fe3 + Cr_n11 + Ti_n11 + Mg_n11 + Mn_n11) > (2/3)-(2/3)*2*t && Si_n11 / (Al_n11 + Fe3 + Cr_n11 + Ti_n11 + Mg_n11 + Mn_n11) < 2+2*t)
                                    T.species(i) = {'siderophyllite'};
                                    OH = 2 - (Cl_n11 + F_n11);
                                    if(OH<=0), OH=0; end
                                    vac = 1 - (Na_n11 + K_n11 + Ca_n11);
                                    if(vac <= 0), vac = 0; end
                                    [sorted_formula1, total1] = Formula_Output({Na_n11, K_n11, Ca_n11, vac},{'Na','K','Ca', '□'});
                                    [sorted_formula2, total2] = Formula_Output({Fe2,Mg_n11,Mn_n11},{'Fe2+','Mg','Mn'});
                                    [sorted_formula3, total3] = Formula_Output({Al_VI,Cr_n11,Ti_n11},{'AlVI','Cr','Ti'});
                                    [sorted_formula4, total4] = Formula_Output({Al_IV,Si_n11},{'AlIV','Si'});
                                    [sorted_formula5, total5] = Formula_Output({OH, Cl_n11, F_n11}, {'OH','Cl','F'});
                                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3, total3,'(',sorted_formula4,total4,'O10)',sorted_formula5, total5)};
                                end
                            end
                        end

                    end  

                end
            end
        end %end of Trioctahderal Micas

        %% 7.8.1.2 Dioctahedral Micas
        if(strcmp(T.group4(i), 'Dioctahedral Mica') && S_n <= 0.05 && P_n <= 0.05)
            Al_IV = 4 - Si_n11;
            if(Al_IV < 0), Al_IV = 0;end
            Al_VI = Al_n11 - Al_IV;
            if(Al_VI < 0), Al_VI = 0; end

            if(K_n11/(K_n11 + Na_n11 + Ca_n11) > 0.5 && Al_VI/(Fe_n + Cr_n + Mg_n + Mn_n) > 0.5)
                T.species(i) = {'muscovite'};
                OH = 2 - (Cl_n11 + F_n11);
                if(OH<=0), OH = 0; end
                vac = 1 - (Na_n11 + K_n11);
                if(vac <= 0), vac = 0; end
                [sorted_formula1, total1] = Formula_Output({Na_n11, K_n11,vac},{'Na','K','□'});
                [sorted_formula2,total2] = Formula_Output({Mg_n11, Fe_n11, Al_VI, Ti_n11, Cr_n11},{'Mg','Fe','AlVI','Ti','Cr'});
                [sorted_formula3,total3] = Formula_Output({Si_n11, Al_IV},{'Si','AlIV'});
                [sorted_formula4,total4] = Formula_Output({OH,Cl_n11,F_n11},{'OH','Cl','F'});
                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3,total3,'O10',sorted_formula4,total4)};

                if(Si_n11 / Al_IV > 3 + 3*t && Al_n/Si_n > 0.3)
                    if(Mg_n11 + Fe_n11 > 0)
                        T.species(i) = {'illite'};
                        T.formula(i) = {'KAl2(Si3Al)O10(OH)2'};
                    end
                end

            end

            if((Na_n11+K_n11+Ca_n11)/Al_VI > (1/2)-2*(1/2)*t && (Na_n11+K_n11+Ca_n11)/Al_VI < (1/2)+(1/2)*t)
                if(Na_n11 > K_n11)
                    if((Al_IV+Al_VI)/Si_n11 > 1-2*1*t && (Al_IV+Al_VI)/Si_n11 < 1+2*1*t)
                        if((Mg_n11+Fe_n11) < Al_VI)
                            T.species(i) = {'paragonite'};
                            renorm = 11/24;
                            Al_IV = 4 - Si_n*renorm;
                            if(Al_IV < 0), Al_IV = 0; end
                            Al_VI = Al_n*renorm - Al_IV;
                            if(Al_VI < 0), Al_VI = 0; end
                            OH = 1 - (F_n*renorm + Cl_n*renorm);
                            if(OH<=0), OH = 0; end
                            [sorted_formula1, total1] = Formula_Output({Ca_n*renorm, Na_n*renorm, K_n*renorm},{'Ca','Na','K'});
                            [sorted_formula2, total2] = Formula_Output({Mg_n*renorm, Fe_n*renorm, Mn_n*renorm, Al_VI, Ti_n*renorm, Cr_n*renorm},{'Mg','Fe','Mn','AlVI','Ti','Cr'});
                            [sorted_formula3, total3] = Formula_Output({Si_n*renorm, Al_IV},{'Si','AlIV'});
                            [sorted_formula4, total4] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                            T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2, total2, sorted_formula3, total3, 'O10',sorted_formula4, total4)};
                        end
                    end
                end
            end

            if(Ca_n11 > Na_n11 + K_n11)
                if(Si_n11 / Al_n11 > (1/2)-(1/2)*t && Si_n11 / Al_n11 < (1/2)+(1/2)*t)
                    if(Al_IV / Al_VI > 1-1*t && Al_IV / Al_VI < 1+1*t)
                        T.species(i) = {'margarite'};
                        Al_VI = 2 - (Fe_n11+Cr_n11+Ti_n11+Mg_n11);
                        if(Al_VI<=0), Al_VI = 0; end
                        Al_IV = Al_n11 - Al_VI;
                        if(Al_IV<=0), Al_IV = 0; end
                        OH = 2 - (F_n11 + Cl_n11);
                        if(OH<=0), OH = 0; end
                        [sorted_formula1, total1] = Formula_Output({Ca_n11, Na_n11, K_n11}, {'Ca','Na','K'});
                        [sorted_formula2, total2] = Formula_Output({Al_VI, Fe_n11, Cr_n11, Ti_n11, Mg_n11}, {'AlVI','Fe','Cr','Ti','Mg'});
                        [sorted_formula3, total3] = Formula_Output({Si_n11, Al_IV},{'Si','AlIV'});
                        [sorted_formula4, total4] = Formula_Output({OH, F_n11, Cl_n11}, {'OH','F','Cl'});
                        T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, 'O10', sorted_formula4, total4)};

                    end
                end
            end

            if(Al_n11/Si_n11 <= 0.2)
                if(Mg_n11/Fe_n11 > 0.9-3*0.9*t && Mg_n11/Fe_n11 < 1+2*1*t) 
                    T.species(i) = {'celadonite'};
                    OH = 2 - (Cl_n11 + F_n11);
                    Fe2 = 1 - (Mg_n11+Mn_n11);
                    if(Fe2<0), Fe2 = 0; end
                    Fe3 = Fe_n11 - Fe2;
                    if(Fe3<0), Fe3 = 0; end
                    [sorted_formula1,total1] = Formula_Output({K_n11,Na_n11,Ba_n11},{'K','Na','Ba'});
                    [sorted_formula2,total2] = Formula_Output({Mg_n11, Fe2, Mn_n11},{'Mg','Fe2+','Mn'});
                    [sorted_formula3,total3] = Formula_Output({Fe3,Al_VI,Ti_n11,V_n11,Cr_n11},{'Fe3+','AlVI','Ti','V','Cr'});
                    [sorted_formula4,total4] = Formula_Output({Si_n11, Al_IV},{'Si','AlIV'});
                    [sorted_formula5,total5] = Formula_Output({OH,Cl_n11,F_n11},{'OH','Cl','F'});
                    T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3,total3,'[',sorted_formula4,total4,'O10]',sorted_formula5,total5)};
                end
            end

            if(Al_IV < 0.5+0.5*t && Al_VI / Si_n11 > (1/4)-(1/4)*t && Al_VI / Si_n11 < (1/4)+(1/4)*t)

                if(Fe_n11 < Mg_n11)
                    T.species(i) = {'aluminoceladonite'};
                        Al_VIx = 1 - (Fe_n11+Mg_n11+Mn_n11+Ni_n11);
                        Al_VIy = Al_VI - Al_VIx;
                        [sorted_formula1, total1] = Formula_Output({K_n11, Na_n11, Ca_n11},{'K','Na','Ca'});
                        [sorted_formula2, total2] = Formula_Output({Fe_n11, Mg_n11,Mn_n11,Ni_n11, Al_VIx},{'Fe','Mg','Mn','Ni','AlVI'});
                        [sorted_formula3, total3] = Formula_Output({Al_VIy,Ti_n11, Cr_n11, V_n11},{'AlVI','Ti','Cr','V'});
                        [sorted_formula4, total4] = Formula_Output({Si_n11, Al_IV},{'Si','AlIV'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3, sorted_formula4,total4, 'O10(OH)2')};
                end

                if(Fe_n11 > Mg_n11)
                    Fe3 = 1 - (Al_VI + Ti_n11 + Cr_n11);
                    if(Fe3 < 0), Fe3 = 0; end

                    if(Al_VI > (Fe3 + Cr_n11 + Ti_n11))
                        T.species(i) = {'ferroaluminoceladonite'};
                        Al_VIx = 1 - (Fe_n11+Mg_n11+Mn_n11+Ni_n11);
                        Al_VIy = Al_VI - Al_VIx;
                        [sorted_formula1, total1] = Formula_Output({K_n11, Na_n11, Ca_n11},{'K','Na','Ca'});
                        [sorted_formula2, total2] = Formula_Output({Fe_n11, Mg_n11,Mn_n11,Ni_n11, Al_VIx},{'Fe','Mg','Mn','Ni','AlVI'});
                        [sorted_formula3, total3] = Formula_Output({Al_VIy,Ti_n11, Cr_n11, V_n11},{'AlVI','Ti','Cr','V'});
                        [sorted_formula4, total4] = Formula_Output({Si_n11, Al_IV},{'Si','AlIV'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, sorted_formula3,total3, sorted_formula4,total4, 'O10(OH)2')};
                    end

                end
            end

        end %end of Dioctahdedral Micas
        
        % Catch
        Al_IV = 4*(24/11) - Si_n;
        if(Al_IV < 0), Al_IV = 0; end
        Al_VI = Al_n - Al_IV;
        if(Al_VI < 0), Al_VI = 0; end
        Fe3 = 1 - (Al_VI + Ti_n + Cr_n);
        if(Fe3 < 0), Fe3 = 0; end
        Fe2 = Fe_n - Fe3;
        if(Fe2 < 0), Fe2 = 0; end

        if(strcmp(T.group3(i), 'Mica'))
            if(strcmp(T.group4(i), 'Dioctahedral Mica'))
                if(Al_n > Si_n && K_n/(Ca_n+Na_n) > 1 && Fe_n > Mg_n)
                    if((Mg_n + Fe2 + Mn_n)/(Al_n + Fe3 + Cr_n + Ti_n) > (1/3)-2*(1/3)*t && (Mg_n + Fe2 + Mn_n)/(Al_n + Fe3 + Cr_n + Ti_n) < (2/3)+3*(2/3)*t)
                        T.group4(i) = {'Trioctahedral Mica (biotite)'};
                        T.species(i) = {'siderophyllite'};
                        OH = 2 - (Cl_n11 + F_n11);
                        if(OH<=0), OH=0; end
                        vac = 1 - (Na_n11 + K_n11 + Ca_n11);
                        if(vac <= 0), vac = 0; end
                        [sorted_formula1, total1] = Formula_Output({Na_n11, K_n11, Ca_n11, vac},{'Na','K','Ca', '□'});
                        [sorted_formula2, total2] = Formula_Output({Fe2,Mg_n11,Mn_n11},{'Fe2+','Mg','Mn'});
                        [sorted_formula3, total3] = Formula_Output({Al_VI,Cr_n11,Ti_n11},{'AlVI','Cr','Ti'});
                        [sorted_formula4, total4] = Formula_Output({Al_IV,Si_n11},{'AlIV','Si'});
                        [sorted_formula5, total5] = Formula_Output({OH, Cl_n11, F_n11}, {'OH','Cl','F'});
                        T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,sorted_formula3, total3,'(',sorted_formula4,total4,'O10)',sorted_formula5, total5)};
                    end
                end
            end
        end

        %% 7.8.2 Talc-Pyrophyllite group
        %% 7.8.2.1 Talc 
        if(strcmp(T.group1(i), 'Silicate') && S_n <= 0.05 && P_n <= 0.05)

            if((Si_n + Al_n) / 24 > (4/11)-(4/11)*t && (Si_n + Al_n) / 24 < (4/11)+(4/11)*t)
                if(Al_n / Si_n < 0.1 && Fe_n / Mg_n < 0.3 && Ca_n / Mg_n < 0.1)
                    if((Mg_n + Fe_n + Ca_n) / (Si_n + Al_n) > 0.75-0.75*t && (Mg_n + Fe_n + Ca_n) / (Si_n + Al_n) < 0.75+0.75*t)
                        if(ABCDT_n > 15.27-15.27*t && ABCDT_n < 15.27+15.27*t)
                            Renorm_factor = 24 / 11;
                            Mg_n11 = T.Mg_Multi(i) / Renorm_factor;
                            Mn_n11 = T.Mn_Multi(i) / Renorm_factor;
                            Fe_n11 = T.Fe_Multi(i) / Renorm_factor;
                            Ca_n11 = T.Ca_Multi(i) / Renorm_factor;
                            Si_n11 = T.Si_Multi(i) / Renorm_factor;
                            Ni_n11 = T.Ni_Multi(i) / Renorm_factor;
                            Ti_n11 = T.Ti_Multi(i) / Renorm_factor;
                            Cr_n11 = T.Cr_Multi(i) / Renorm_factor;
                            Na_n11 = T.Na_Multi(i) / Renorm_factor;
                            Al_n11 = T.Al_Multi(i) / Renorm_factor;
                            T.group2(i) = {'Phyllosilicate'};
                            T.group3(i) = {'Talc-Pyrophyllite'};
                            T.group4(i) = {''};
                            T.species(i) = {'talc'};
                            OH = 2 - (F_n11+Cl_n11);
                            if(OH<=0), OH = 0; end
                            [sorted_formula1, total1] = Formula_Output({Fe_n11, Mg_n11, Ni_n11, Mn_n11}, {'Fe', 'Mg', 'Ni','Mn'});
                            [sorted_formula2, total2] = Formula_Output({Si_n11, Al_n11},{'Si','Al'});
                            [sorted_formula3, total3] = Formula_Output({OH, F_n11, Cl_n11},{'OH','F','Cl'});
                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O10',sorted_formula3, total3)};

                            if(Na_n/Si_n > 0.15/4 && Mg_n/Si_n < 3/4 + (3/4)*t && Mg_n/Si_n > 3/4 - (3/4)*t && Na_n > Ca_n && Na_n > K_n)
                                T.group3(i) = {'Smectite'};
                                T.group4(i) = {'Trioctahedral Smectite'};
                                T.species(i) = {''};
                            end

                        end
                    end
                end
            end

            if((Si_n + Al_n) / 24 > (4/11)-(4/11)*t && (Si_n + Al_n) / 24 < (4/11)+(4/11)*t)
                if(Al_n / Si_n < 0.1 && Fe_n > Mg_n && Fe_n > Mn_n && (Na_n + K_n + Ca_n) / (Mg_n + Fe_n + Mn_n) < 0.1)
                    if((Mg_n + Fe_n + Mn_n + Ca_n + Na_n + K_n + Cr_n + Ti_n) / (Si_n + Al_n) > 0.75-0.75*t && (Mg_n + Fe_n + Mn_n + Ca_n + Na_n + K_n + Cr_n + Ti_n) / (Si_n + Al_n) < 0.75+0.75*t)
                        if(ABCDT_n > 15.27-15.27*t && ABCDT_n < 15.27+15.27*t)
                                Renorm_factor = 24 / 11;
                                Mg_n11 = T.Mg_Multi(i) / Renorm_factor;
                                Mn_n11 = T.Mn_Multi(i) / Renorm_factor;
                                Fe_n11 = T.Fe_Multi(i) / Renorm_factor;
                                Ca_n11 = T.Ca_Multi(i) / Renorm_factor;
                                Si_n11 = T.Si_Multi(i) / Renorm_factor;
                                Ni_n11 = T.Ni_Multi(i) / Renorm_factor;
                                Ti_n11 = T.Ti_Multi(i) / Renorm_factor;
                                Cr_n11 = T.Cr_Multi(i) / Renorm_factor;
                                Na_n11 = T.Na_Multi(i) / Renorm_factor;
                                Al_n11 = T.Al_Multi(i) / Renorm_factor;
                                T.group2(i) = {'Phyllosilicate'};
                                T.group3(i) = {'Talc-Pyrophyllite'};
                                T.group4(i) = {''};
                                T.species(i) = {'minnesotaite'};
                                OH = 2 - (F_n11+Cl_n11);
                                if(OH<=0), OH = 0; end
                                [sorted_formula1, total1] = Formula_Output({Fe_n11, Mg_n11, Ni_n11, Mn_n11}, {'Fe', 'Mg', 'Ni','Mn'});
                                [sorted_formula2, total2] = Formula_Output({Si_n11, Al_n11},{'Si','Al'});
                                [sorted_formula3, total3] = Formula_Output({OH, F_n11, Cl_n11},{'OH','F','Cl'});
                                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O10',sorted_formula3, total3)};

                        end
                    end
                end
            end

        end

        %% 7.8.2.2 Pyrophyllite
        if(strcmp(T.group1(i), 'Silicate'))

            if((Mg_n + Fe_n + Ca_n) / (Si_n + Al_n) >= 0 && (Mg_n + Fe_n + Ca_n) / (Si_n + Al_n) < 0.1+0.1*t && S_n <= 0.05 && P_n <= 0.05)
                if(Al_n / Si_n > 0.5-0.5*t && Al_n / Si_n < 0.5+0.5*t && (Na_n+Ca_n+K_n)/Al_n < 0.25/2)
                    Renorm_factor = 24 / 11;
                    Mg_n11 = T.Mg_Multi(i) / Renorm_factor;
                    Mn_n11 = T.Mn_Multi(i) / Renorm_factor;
                    Fe_n11 = T.Fe_Multi(i) / Renorm_factor;
                    Ca_n11 = T.Ca_Multi(i) / Renorm_factor;
                    K_n11 = T.K_Multi(i) / Renorm_factor;
                    Si_n11 = T.Si_Multi(i) / Renorm_factor;
                    Ni_n11 = T.Ni_Multi(i) / Renorm_factor;
                    Ti_n11 = T.Ti_Multi(i) / Renorm_factor;
                    Cr_n11 = T.Cr_Multi(i) / Renorm_factor;
                    Na_n11 = T.Na_Multi(i) / Renorm_factor;
                    Al_n11 = T.Al_Multi(i) / Renorm_factor;
                    F_n11 = T.F_Multi(i) / Renorm_factor;
                    Cl_n11 = T.Cl_Multi(i) / Renorm_factor;
                    Al_IV = 4 - Si_n11;
                    if(Al_IV < 0), Al_IV = 0; end
                    Al_VI = Al_n11 - Al_IV;
                    if(Al_VI < 0), Al_VI = 0; end
                    Fe3 = Fe_n11;
                    OH = 2 - (F_n11+Cl_n11);
                    if(OH<=0), OH=0; end
                    T.group2(i) = {'Phyllosilicate'};
                    T.group3(i) = {'Talc-Pyrophyllite'};
                    T.group4(i) = {''};
                    T.species(i) = {'pyrophyllite'};
                    [sorted_formula1, total1] = Formula_Output({Al_VI, Fe3, Mg_n11, Cr_n11, Ti_n11},{'AlVI','Fe3+','Mg','Cr','Ti'});
                    [sorted_formula2, total2] = Formula_Output({Si_n11, Al_IV},{'Si','AlIV'});
                    [sorted_formula3, total3] = Formula_Output({OH, F_n11, Cl_n11},{'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O10',sorted_formula3,total3)};

                    if((Na_n+Ca_n)/Si_n > 0.15/4 && Al_n/Mg_n > 1)
                        T.group3(i) = {'Smectite'};
                    end

                end
            end

        end

        %% 7.8.3 Chlorite group
        if(strcmp(T.group1(i), 'Silicate') && ~strcmp(T.group2(i), 'Nesosilicate') && S_n <= 0.05 && P_n <= 0.05)
            if(ABCDT_n*(14/24) >= 10-2*10*t && ABCDT_n*(14/24) <= 10+2*10*t)
                if((Si_n+Al_n)/(C_n + D_n) >= (4/6)-(4/6)*t && (Si_n+Al_n)/(T.C(i) + T.D(i)) <= (4/6)+5*(4/6)*t)
                    if( (Na_n + K_n + Ca_n) < 0.3 && Al_n/Si_n >= (2/3)-5*(2/3)*t && Al_n/Si_n <= (2/3)+5*(2/3)*t)
                        T.group2(i) = {'Phyllosilicate'};
                        T.group3(i) = {'Chlorite'};
                        T.group4(i) = {''};
                        Renorm_factor = 24 / 14;
                        Mg_n14 = T.Mg_Multi(i) / Renorm_factor;
                        Mn_n14 = T.Mn_Multi(i) / Renorm_factor;
                        Fe_n14 = T.Fe_Multi(i) / Renorm_factor;
                        Si_n14 = T.Si_Multi(i) / Renorm_factor;
                        Ni_n14 = T.Ni_Multi(i) / Renorm_factor;
                        Ti_n14 = T.Ti_Multi(i) / Renorm_factor;
                        Cr_n14 = T.Cr_Multi(i) / Renorm_factor;
                        Al_n14 = T.Al_Multi(i) / Renorm_factor;
                        V_n14 = T.V_Multi(i) / Renorm_factor;
                        Cl_n14 = T.Cl_Multi(i) / Renorm_factor;
                        F_n14 = T.F_Multi(i) / Renorm_factor;
                        Al_IV = 4 - Si_n14;
                        if(Al_IV < 0), Al_IV = 0; end
                        Al_VI = Al_n14 - Al_IV;
                        if(Al_VI < 0), Al_VI = 0; end
                        Fe3 = 1 - Al_VI;
                        if(Fe3 < 0), Fe3 = 0; end
                        Fe2 = Fe_n14 - Fe3;
                        if(Fe2 < 0), Fe2 = 0; end

                        if(Mg_n > Fe_n)
                            T.species(i) = {'clinochlore'};
                            OH = 8 - (F_n14 + Cl_n14);
                            [sorted_formula1, total1] = Formula_Output({Mg_n14, Fe2, Ni_n14, Mn_n14}, {'Mg','Fe2+','Ni','Mn'});
                            [sorted_formula2, total2] = Formula_Output({Al_VI, Ti_n14,Cr_n14,V_n14}, {'AlVI','Ti','Cr','V'});
                            [sorted_formula3, total3] = Formula_Output({Al_IV, Si_n14}, {'AlIV','Si'});
                            [sorted_formula4, total4] = Formula_Output({OH, F_n14, Cl_n14} , {'OH','F','Cl'});
                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, '[',sorted_formula3,total3,'O10]', sorted_formula4,total4)};
                        end

                        if(Fe_n > Mg_n)
                            T.species(i) = {'chamosite'};
                            [sorted_formula1,total1] = Formula_Output({Fe2, Mg_n14, Mn_n14, Cr_n14, Ti_n14, Al_VI, Fe3},{'Fe2+','Mg','Mn','Cr','Ti','AlVI','Fe3+'});
                            [sorted_formula2,total2] = Formula_Output({Si_n14, Al_IV},{'Si','AlIV'});
                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O10(OH,O)8')};
                        end

                    end
                end
            end 
        end

        %% 7.8.4 Kaolinite-Serpentine group
        %% 7.8.4.1 Serpentine
        if(strcmp(T.group1(i), 'Silicate') && S_n <= 0.05 && P_n <= 0.05)
            if(ABCDT_n / 24 > (5/7)-(5/7)*t && ABCDT_n / 24 < (5/7)+(5/7)*t)
                if(Mg_n / Si_n > (3/2)-(3/2)*t && Mg_n / Si_n < (3/2)+(3/2)*t)
                    if(Si_n / 24 > (2/7)-(2/7)*t && Si_n / 24 < (2/7)+(2/7)*t)
                        if(Al_n / Si_n < 0.2 && (Na_n + K_n + Ca_n) / Si_n < 0.05)
                            if(Mg_n/(Mg_n+Fe_n) > 0.8)
                                T.group2(i) = {'Phyllosilicate'};
                                T.group3(i) = {'Kaolinite-Serpentine'};
                                T.group4(i) = {'Serpentine'};
                                Renorm_factor = 24 / 7;
                                Mg_n7 = T.Mg_Multi(i) / Renorm_factor;
                                Mn_n7 = T.Mn_Multi(i) / Renorm_factor;
                                Fe_n7 = T.Fe_Multi(i) / Renorm_factor;
                                Al_n7 = T.Al_Multi(i) / Renorm_factor; 
                                Ni_n7 = T.Ni_Multi(i) / Renorm_factor;
                                Si_n7 = T.Si_Multi(i) / Renorm_factor;
                                Cl_n7 = T.Cl_Multi(i) / Renorm_factor;
                                F_n7 = T.F_Multi(i) / Renorm_factor;
                                Al_IV = 2 - Si_n7;
                                if(Al_IV < 0), Al_IV = 0; end
                                Al_VI = Al_n7 - Al_IV;
                                if(Al_VI < 0), Al_VI = 0; end
                                T.species(i) = {'antigorite/lizardite/chrysotile'};
                                [sorted_formula1, total1] = Formula_Output({Mg_n7,Fe_n7,Mn_n7,Ni_n7,Al_VI},{'Mg','Fe','Mn','Ni','AlVI'});
                                [sorted_formula2, total2] = Formula_Output({Si_n7,Al_IV},{'Si','AlIV'});
                                [sorted_formula3, total3] = Formula_Output({4-(Cl_n7+F_n7),Cl_n7,F_n7},{'OH','Cl','F'});
                                T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'O5',sorted_formula3,total3)};
                            end
                        end
                    end
                end
            end
        end

        %% 7.8.4.2 Kaolinite
        if(strcmp(T.group1(i), 'Silicate') && S_n <= 0.05 && P_n <= 0.05)
            Renorm_factor = 24 / 7;
            Mg_n7 = T.Mg_Multi(i) / Renorm_factor;
            Si_n7 = T.Si_Multi(i) / Renorm_factor;
            Fe_n7 = T.Fe_Multi(i) / Renorm_factor;
            Ca_n7 = T.Ca_Multi(i) / Renorm_factor;
            Al_n7 = T.Al_Multi(i) / Renorm_factor;
            F_n7 = T.F_Multi(i) / Renorm_factor;
            Cl_n7 = T.Cl_Multi(i) / Renorm_factor;
            Cr_n7 = T.Cr_Multi(i) / Renorm_factor;
            Ti_n7 = T.Ti_Multi(i) / Renorm_factor;
            Na_n7 = T.Na_Multi(i)/Renorm_factor;
            K_n7 = T.K_Multi(i)/Renorm_factor;
            
            Al_IV = 2 - Si_n7;
            if(Al_IV<=0), Al_IV = 0; end
            Al_VI = Al_n7 - Al_IV;
            if(Al_VI<=0), Al_VI = 0; end

            if(ABCDT_n / 24 > (4/7)-(4/7)*t && ABCDT_n / 24 < (4/7)+(4/7)*t)
                if((Al_VI+Fe_n7+Cr_n7+Ti_n7+Mg_n7+Ca_n7) / Si_n7 > 1-1*t && (Al_VI+Fe_n7+Cr_n7+Ti_n7+Mg_n7+Ca_n7) / Si_n7 < 1+1*t)
                    if((Si_n7+Al_IV) > (2)-(2)*t && (Si_n7+Al_IV) < (2)+(2)*t)
                        if(Mg_n7 / (Si_n7+Al_IV) < 0.15 && (Na_n7 + K_n7 + Ca_n7) / (Si_n7+Al_IV) < 0.15 && Al_VI > Fe_n7 && Al_VI > Cr_n7 && Al_VI > Ti_n7)
                            T.group2(i) = {'Phyllosilicate'};
                            T.group3(i) = {'Kaolinite-Serpentine'};
                            T.group4(i) = {'Kaolinite'};
                            T.species(i) = {'kaolinite'};
                            
                            OH = 4 - (F_n7+Cl_n7);
                            if(OH<=0), OH = 0; end
                            Fe3 = Fe_n7;
                            
                            [sorted_formula1, total1] = Formula_Output({Al_VI, Fe3, Mg_n7, Ca_n7}, {'AlVI','Fe3+','Mg','Ca'});
                            [sorted_formula2, total2] = Formula_Output({Si_n7, Al_IV},{'Si','AlIV'});
                            [sorted_formula3, total3] = Formula_Output({OH, F_n7, Cl_n7}, {'OH','F','Cl'});
                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2, 'O5',sorted_formula3, total3)};
                        end
                    end
                end
            end
        end

        %% 7.8.4.3 Greenalite
        if(strcmp(T.group1(i), 'Silicate') && S_n <= 0.05 && P_n <= 0.05)
            Renorm_factor = 24 / 7;
            Mg_n7 = T.Mg_Multi(i) / Renorm_factor;
            Fe_n7 = T.Fe_Multi(i) / Renorm_factor;
            Ca_n7 = T.Ca_Multi(i) / Renorm_factor;
            Al_n7 = T.Al_Multi(i) / Renorm_factor;
            Si_n7 = T.Si_Multi(i) / Renorm_factor;
            Na_n7 = T.Na_Multi(i) / Renorm_factor;
            Mn_n7 = T.Mn_Multi(i) / Renorm_factor;
            K_n7 = T.K_Multi(i) / Renorm_factor;
            F_n7 = T.F_Multi(i) / Renorm_factor;
            Cl_n7 = T.Cl_Multi(i) / Renorm_factor;
            V_n7 = T.V_Multi(i) / Renorm_factor;
            Cr_n7 = T.Cr_Multi(i) / Renorm_factor;
            ABCDT_n7 = T.ABCDT(i) / Renorm_factor;
            ABCD_n7 = T.ABCD(i) / Renorm_factor;
            Al_IV = 2 - Si_n7;
            if(Al_IV < 0), Al_IV = 0; end
            Al_VI = Al_n7 - Al_IV;
            if(Al_VI < 0), Al_VI = 0; end

            if((Fe_n7+Mg_n7+Mn_n7+Al_VI)/(Si_n7+Al_IV) > 2/2 - (2/2)*t && (Fe_n7+Mg_n7+Mn_n7+Al_VI)/(Si_n7+Al_IV) < 3/2 + (3/2)*t)
                if(Fe_n7/(Fe_n7+Mg_n7+Mn_n7+Al_VI) > 0.75-0.75*t && round(Fe_n7,2) <= 3.00)
                    if(Al_IV/Si_n7 >= 0 && Al_IV/Si_n7 <= 0.5+0.5*t && Al_VI/Si_n7 >= 0 && Al_VI/Si_n7 <= 1/3)
                        if(Na_n7+Ca_n7+K_n7 < 0.1 && Fe_n/Al_n > 3 && Fe_n/Mg_n > 1 && Fe_n > Mn_n)
                            if(ABCDT_n7 > 4 - 4*t && ABCDT_n7 < 5 + 5*t && round(ABCD_n7,2) <= 3.00)
                                T.group2(i) = {'Phyllosilicate'};
                                T.group3(i) = {'Kaolinite-Serpentine'};
                                T.group4(i) = {'Serpentine'};
                                T.species(i) = {'greenalite'};
                                [sorted_formula1, total1] = Formula_Output({Fe_n7, Mg_n7, Mn_n7, Al_VI},{'Fe','Mg','Mn','AlVI'});
                                [sorted_formula2, total2] = Formula_Output({Si_n7, Al_IV},{'Si','AlIV'});
                                [sorted_formula3, total3] = Formula_Output({4 - (F_n7+Cl_n7), F_n7, Cl_n7},{'OH','F','Cl'});
                                T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, 'O5', sorted_formula3, total3)};

                                if((Fe_n7+Mg_n7+Mn_n7+Al_VI)/(Si_n7+Al_IV) > 1 - 1*(0.5)*t && (Fe_n7+Mg_n7+Mn_n7+Al_VI)/(Si_n7+Al_IV) < 1 + 1*(0.5)*t)
                                    T.species(i) = {'hisingerite'};
                                    [sorted_formula1, total1] = Formula_Output({Fe_n7, Mg_n7, Mn_n7, Al_VI, Cr_n7, V_n7},{'Fe3+','Mg','Mn','AlVI','Cr','V'});
                                    [sorted_formula2, total2] = Formula_Output({Si_n7, Al_IV},{'Si','AlIV'});
                                    [sorted_formula3, total3] = Formula_Output({4 - (F_n7+Cl_n7), F_n7, Cl_n7},{'OH','F','Cl'});
                                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, 'O5', sorted_formula3, total3, '·2H2O')};
                                end

                            end
                        end
                    end
                end
            end

        end

        %% 7.8.5 Smectites
        if(strcmp(T.group1(i), 'Silicate') && ~strcmp(T.group3(i), 'Amphibole') && ~strcmp(T.group3(i), 'Pyroxene') && ~strcmp(T.group2(i), 'Tectosilicate') && ~strcmp(T.group2(i), 'Nesosilicate') ...
                && ~strcmp(T.group3(i), 'Mica') && ~strcmp(T.group3(i), 'Talc-Pyrophyllite') && ~strcmp(T.group3(i), 'Kaolinite-Serpentine') && ~strcmp(T.group3(i), 'Palygorskite-Sepiolite') && ~strcmp(T.group3(i), 'Chlorite') && ~strcmp(T.group2(i), 'Cyclosilicate') && ~strcmp(T.group2(i), 'Nesoosilicate')...
                && S_n <= 0.05 && P_n <= 0.05)

            Renorm_factor = 24 / 11;
            Mg_n11 = T.Mg_Multi(i) / Renorm_factor;
            Mn_n11 = T.Mn_Multi(i) / Renorm_factor;
            Fe_n11 = T.Fe_Multi(i) / Renorm_factor;
            Ca_n11 = T.Ca_Multi(i) / Renorm_factor;
            Si_n11 = T.Si_Multi(i) / Renorm_factor;
            Ni_n11 = T.Ni_Multi(i) / Renorm_factor;
            Ti_n11 = T.Ti_Multi(i) / Renorm_factor;
            Cr_n11 = T.Cr_Multi(i) / Renorm_factor;
            Na_n11 = T.Na_Multi(i) / Renorm_factor;
            Al_n11 = T.Al_Multi(i) / Renorm_factor;
            Cl_n11 = T.Cl_Multi(i) / Renorm_factor;
            F_n11 = T.F_Multi(i) / Renorm_factor;
            Zn_n11 = T.Zn_Multi(i) / Renorm_factor;
            V_n11 = T.V_Multi(i) / Renorm_factor;
            K_n11 = T.K_Multi(i) / Renorm_factor;
            ABCDT_n11 = T.ABCDT(i) / Renorm_factor;
            Al_IV = 4 - Si_n11;
            if(Al_IV < 0),Al_IV = 0; end
            Al_VI = Al_n11 - Al_IV;
            if(Al_VI < 0), Al_VI = 0; end

            if( (ABCDT_n11 > 6.3 - 6.3*t && ABCDT_n11 < 7.4 + 7.4*t && ...
                    (Na_n11 + Ca_n11+K_n11) / (Si_n11 + Al_IV) <= 0.3/4 + 4*(0.3/4)*t && K_n11/(Na_n11+Ca_n11) < 0.2 && ...
                    Al_IV/Si_n11 >= 0 && Al_IV/Si_n11 < 0.5 + 0.5*t && (Si_n11+Al_IV) >= 4 - 4*t && (Si_n11+Al_IV) <= 4 + 4*t) || ...
                    ...
                    (Zn_n11/(Si_n11+Al_IV) >= 0.3 - 0.3*t && Zn_n11/(Si_n11+Al_IV) <= (3/4)+(3/4)*t && (Na_n11+K_n11) > 0 && (Na_n11+K_n11) < 0.5) || ...
                    ...
                    ((Mg_n11+Al_n11+Fe_n11+Al_VI)/(Si_n11+Al_IV) > 1/2 - 2*(1/2)*t && (Mg_n11+Al_n11+Fe_n11+Al_VI)/(Si_n11+Al_IV) < 1/2 + 3*(1/2)*t && (Na_n11+Ca_n11+K_n11)/Si_n11 < 0.005) )

                T.group2(i) = {'Phyllosilicate'};
                T.group3(i) = {'Smectite'};

                if(Al_VI > Mg_n11 + Fe_n11 + Zn_n11 + Cr_n11 + V_n11)
                    T.group4(i) = {'Dioctahedral Smectite'};
                    T.species(i) = {''};
                    T.formula(i) = {''};
                end

                if(Al_VI <= Mg_n11 + Fe_n11 + Zn_n11 + Cr_n11 + V_n11)
                    T.group4(i) = {'Trioctahedral Smectite'};
                    T.species(i) = {''};
                    T.formula(i) = {''};
                end

                %% 7.8.5.1. Dioctahedral smectites
                if(strcmp(T.group4(i), 'Dioctahedral Smectite'))

                    Renorm_factor = 24 / 11;
                    Mg_n11 = T.Mg_Multi(i) / Renorm_factor;
                    Mn_n11 = T.Mn_Multi(i) / Renorm_factor;
                    Fe_n11 = T.Fe_Multi(i) / Renorm_factor;
                    Ca_n11 = T.Ca_Multi(i) / Renorm_factor;
                    Si_n11 = T.Si_Multi(i) / Renorm_factor;
                    Ni_n11 = T.Ni_Multi(i) / Renorm_factor;
                    Ti_n11 = T.Ti_Multi(i) / Renorm_factor;
                    Cr_n11 = T.Cr_Multi(i) / Renorm_factor;
                    Na_n11 = T.Na_Multi(i) / Renorm_factor;
                    Al_n11 = T.Al_Multi(i) / Renorm_factor;
                    Cl_n11 = T.Cl_Multi(i) / Renorm_factor;
                    F_n11 = T.F_Multi(i) / Renorm_factor;
                    Zn_n11 = T.Zn_Multi(i) / Renorm_factor;
                    V_n11 = T.V_Multi(i) / Renorm_factor;
                    K_n11 = T.K_Multi(i) / Renorm_factor;
                    Al_IV = 4 - Si_n11;
                    if(Al_IV < 0),Al_IV = 0; end
                    Al_VI = Al_n11 - Al_IV;
                    if(Al_VI < 0), Al_VI = 0; end

                    if(Al_VI / (Mg_n11 + Al_VI + Fe_n11 + Zn_n11 + Cr_n11 + V_n11) > 0.5)
                        if( (Al_VI + Mg_n11 + Fe_n11 + Zn_n11 + Cr_n11 + V_n11 + Ti_n11) > 2-2*t && (Al_VI + Mg_n11 + Fe_n11 + Zn_n11 + Cr_n11 + V_n11 + Ti_n11) < 2+2*t&& Na_n > Ca_n && Al_IV < 0.5)
                            if(Na_n11 + Ca_n11 + K_n11 > 0.3-0.3*t && Na_n11 + Ca_n11 + K_n11 < 0.3+0.3*t)
                                T.species(i) = {'montmorillonite'};
                                OH = 2 - (F_n11+Cl_n11);
                                if(OH <=0), OH = 0; end
                                [sorted_formula1, total1] = Formula_Output({Na_n11, Ca_n11, K_n11},{'Na','Ca','K'});
                                [sorted_formula2, total2] = Formula_Output({Al_VI, Mg_n11, Fe_n11, Mn_n11},{'AlVI','Mg','Fe','Mn'});
                                [sorted_formula3, total3] = Formula_Output({Si_n11, Al_IV},{'Si','AlIV'});
                                [sorted_formula4, total4] = Formula_Output({OH, F_n11, Cl_n11},{'OH','F','Cl'});
                                T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, sorted_formula3, total3, 'O10',sorted_formula4,total4,'·nH2O')};
                            end
                        end
                    end

                    if(Al_VI / (Mg_n11 + Fe_n11 + Zn_n11 + Cr_n11 + V_n11) > 0.75)
                        if((Ca_n + Na_n)/(Ca_n+Na_n+K_n) > 0.5 && Al_IV > 0.5 && Na_n11 + Ca_n11 + K_n11 > 0.3-2*0.3*t && Na_n11 + Ca_n11 + K_n11 < 0.3+0.3*t)
                            OH = 2 - (Cl_n11+F_n11);
                            if(OH < 0), OH = 0; end
                            T.species(i) = {'beidellite'};
                            Fe3 = Fe_n*(11/24);
                            [sorted_formula1, total1]  = Formula_Output({Na_n11, Ca_n11,K_n11}, {'Na', 'Ca','K'});
                            [sorted_formula2, total2]  = Formula_Output({Al_VI, Fe3}, {'AlVi', 'Fe3+'});
                            [sorted_formula3, total3]  = Formula_Output({Si_n11, Al_IV}, {'Si', 'AlIV'});
                            [sorted_formula4, total4]  = Formula_Output({OH, Cl_n11,F_n11}, {'OH', 'Cl','F'});
                            T.formula(i) = {strcat(sorted_formula1,total1,sorted_formula2,total2,'[',sorted_formula3,total3,'O10]',sorted_formula4,total4,'·nH2O')};
                        end
                    end

                end

                %% 7.8.5.2. Trioctahedral smectites
                if(strcmp(T.group4(i), 'Trioctahedral Smectite'))

                    Renorm_factor = 24 / 11;
                    Mg_n11 = T.Mg_Multi(i) / Renorm_factor;
                    Mn_n11 = T.Mn_Multi(i) / Renorm_factor;
                    Fe_n11 = T.Fe_Multi(i) / Renorm_factor;
                    Ca_n11 = T.Ca_Multi(i) / Renorm_factor;
                    Si_n11 = T.Si_Multi(i) / Renorm_factor;
                    Ni_n11 = T.Ni_Multi(i) / Renorm_factor;
                    Ti_n11 = T.Ti_Multi(i) / Renorm_factor;
                    Cr_n11 = T.Cr_Multi(i) / Renorm_factor;
                    Na_n11 = T.Na_Multi(i) / Renorm_factor;
                    Al_n11 = T.Al_Multi(i) / Renorm_factor;
                    Cl_n11 = T.Cl_Multi(i) / Renorm_factor;
                    F_n11 = T.F_Multi(i) / Renorm_factor;
                    Zn_n11 = T.Zn_Multi(i) / Renorm_factor;
                    V_n11 = T.V_Multi(i) / Renorm_factor;
                    K_n11 = T.K_Multi(i) / Renorm_factor;
                    Al_IV = 4 - Si_n11;
                    if(Al_IV < 0),Al_IV = 0; end
                    Al_VI = Al_n11 - Al_IV;
                    if(Al_VI < 0), Al_VI = 0; end

                    if(Fe_n11 / (Al_VI + Mg_n11 + Fe_n11 + Zn_n11 + Cr_n11 + V_n11) > 0.5 && Fe_n11 > Al_VI)
                        if((Ca_n11+Na_n11) > K_n11 && (Na_n + Ca_n + K_n)/(Si_n + Al_IV) >= 0.3/4 - (0.3/4)*t && (Na_n + Ca_n + K_n)/(Si_n + Al_IV) < 0.3/4 + 3*(0.3/4)*t)
                            if((Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11)/(Si_n11+Al_IV) > 1.5/4 - (1.5/4)*t && (Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11)/(Si_n11+Al_IV) < 2.5/4 + 3*(2.5/4)*t)
                                Fe3 = Fe_n11;
                                T.species(i) = {'nontronite'};
                                renorm = 11/24;
    
                                Al_IV = 4 - Si_n*renorm;
                                if(Al_IV<=0), Al_IV = 0; end
                                Al_VI = Al_n*renorm - Al_IV;
                                if(Al_VI<=0), Al_VI = 0; end
                                OH = 2 - (F_n*renorm + Cl_n*renorm);
                                if(OH<=0), OH = 0; end
                                [sorted_formula1, total1] = Formula_Output({Na_n*renorm, K_n*renorm, Ca_n*renorm},{'Na','K','Ca'});
                                [sorted_formula2, total2] = Formula_Output({Fe_n*renorm, Al_VI, Mg_n*renorm, Mn_n*renorm},{'Fe','AlVI','Mg','Mn'});
                                [sorted_formula3, total3] = Formula_Output({Si_n*renorm, Al_IV},{'Si','AlIV'});
                                [sorted_formula4, total4] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                                T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, sorted_formula3, total3, 'O10', sorted_formula4,total4,'·nH2O')};
                            end
                        end
                    end

                    if((Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11+Ti_n11)/(Si_n11+Al_IV) >= (3/4)-2*(3/4)*t && (Mg_n11+Mn_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11+Ti_n11)/(Si_n11+Al_IV) <= (3/4)+(3/4)*t)
                        if((Na_n11+Ca_n11+K_n11)/(Si_n11+Al_IV) > 0.3/4 - (0.3/4)*t && (Na_n11+Ca_n11+K_n11)/(Si_n11+Al_IV) < 0.4/4 + (0.4/4)*t && Na_n/(Na_n+K_n+Ca_n) > 0.5 && Mg_n11 > Fe_n11 && Al_IV/Si_n11 < 0.1 && Na_n11 > Ca_n11 && Na_n11 > K_n11)
                            T.species(i) = {'hectorite'};
                            renorm = 11/24;
                            Al_IV = 4 - Si_n*renorm;
                            if(Al_IV <=0), Al_IV = 0; end
                            Al_VI = Al_n*renorm - Al_IV;
                            if(Al_VI <=0), Al_VI = 0; end
                            Li = 3 - (Mg_n*renorm + Fe_n*renorm + Al_VI);
                            if(Li < 0), Li = 0; end
                            OH = 2 - (F_n*renorm + Cl_n*renorm);
                            if(OH < 0), Li = 0; end
                            [sorted_formula1, total1] = Formula_Output({Na_n*renorm, Ca_n*renorm, K_n*renorm},{'Na','Ca','K'});
                            [sorted_formula2, total2] = Formula_Output({Mg_n*renorm, Li, Fe_n*renorm, Al_VI, Mn_n*renorm},{'Mg','Li','Fe','AlVI','Mn'});
                            [sorted_formula3, total3] = Formula_Output({Si_n*renorm, Al_IV},{'Si','AlIV'});
                            [sorted_formula4, total4] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                            T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2,total2,sorted_formula3,total3,'O10',sorted_formula4,total4,'·nH2O')};
                        end
                    end

                    if((Ca_n11+Na_n11+K_n11)/(Mg_n11+Fe_n11+Al_n11+Zn_n11+Cr_n11+V_n11+Si_n11) > 0.275/7 - 2*(0.275/7)*t && (Ca_n11+Na_n11+K_n11)/(Mg_n11+Fe_n11+Al_n11+Zn_n11+Cr_n11+V_n11+Si_n11) < 0.275/7 + 2*(0.275/7)*t)
                        if(Ca_n11 > Na_n11 && Na_n11 > K_n11)
                            if((Mg_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11)/(Si_n11+Al_IV) > 3/4 - (3/4)*t && (Mg_n11+Fe_n11+Al_VI+Zn_n11+Cr_n11+V_n11)/(Si_n11+Al_IV) < 3/4 + (3/4)*t)
                                if(Mg_n11 > Fe_n11 && Fe_n11 > Al_VI && Si_n11 > Al_IV)
                                    T.species(i) = {'saponite'};
                                    renorm = 11/24;
                                    Al_IV = 4 - Si_n*renorm;
                                    if(Al_IV<=0), Al_IV = 0; end
                                    Al_VI = Al_n*renorm - Al_IV;
                                    if(Al_VI<=0), Al_VI = 0; end
                                    OH = 2 - (F_n*renorm + Cl_n*renorm);
                                    if(OH<=0), OH = 0; end
                                    [sorted_formula1, total1] = Formula_Output({Na_n*renorm, K_n*renorm, Ca_n*renorm},{'Na','K','Ca'});
                                    [sorted_formula2, total2] = Formula_Output({Fe_n*renorm, Al_VI, Mg_n*renorm, Mn_n*renorm, Ti_n*renorm, Cr_n*renorm},{'Fe','AlVI','Mg','Mn', 'Ti','Cr'});
                                    [sorted_formula3, total3] = Formula_Output({Si_n*renorm, Al_IV},{'Si','AlIV'});
                                    [sorted_formula4, total4] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                                    T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, sorted_formula3, total3, 'O10', sorted_formula4,total4,'·4H2O')};
                                end
                            end
                        end
                    end

                    if(Zn_n11 / (Mg_n11 + Fe_n11 + Al_VI + Zn_n11 + Cr_n11 + V_n11) > 0.5)
                        if(Na_n11 >=0 && Na_n11 < 0.5)
                            T.species(i) = {'sauconite'};
                            renorm = 11/24;
                            Al_IV = 4 - Si_n*renorm;
                            if(Al_IV<=0), Al_IV = 0; end
                            Al_VI = Al_n*renorm - Al_IV;
                            if(Al_VI<=0), Al_VI = 0; end
                            OH = 2 - (F_n*renorm + Cl_n*renorm);
                            [sorted_formula1, total1] = Formula_Output({Na_n*renorm, K_n*renorm, Ca_n*renorm},{'Na','K','Ca'});
                            [sorted_formula2, total2] = Formula_Output({Fe_n*renorm, Al_VI, Mg_n*renorm, Mn_n*renorm, Cu_n*renorm, Zn_n*renorm},{'Fe','AlVI','Mg','Mn', 'Cu','Zn'});
                            [sorted_formula3, total3] = Formula_Output({Si_n*renorm, Al_IV},{'Si','AlIV'});
                            [sorted_formula4, total4] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                            T.formula(i) = {strcat(sorted_formula1,total1, sorted_formula2,total2, sorted_formula3, total3, 'O10', sorted_formula4,total4,'·4H2O')};
                        end
                    end

                    if(Mg_n11 /(Mg_n11 + Fe_n11 + Al_VI + Zn_n11 + Cr_n11 + V_n11) > 0.5)
                        if((Mg_n11 + Fe_n11 + Al_VI) > 3)
                            Renorm_factor = 24 / 22;
                            Mg_n22 = T.Mg_Multi(i) / Renorm_factor;
                            Fe_n22 = T.Fe_Multi(i) / Renorm_factor;
                            Ca_n22 = T.Ca_Multi(i) / Renorm_factor;
                            Si_n22 = T.Si_Multi(i) / Renorm_factor;
                            Ti_n22 = T.Ti_Multi(i) / Renorm_factor;
                            Cr_n22 = T.Cr_Multi(i) / Renorm_factor;
                            Na_n22 = T.Na_Multi(i) / Renorm_factor;
                            Al_n22 = T.Al_Multi(i) / Renorm_factor;
                            Cl_n22 = T.Cl_Multi(i) / Renorm_factor;
                            Zn_n22 = T.Zn_Multi(i) / Renorm_factor;
                            V_n22 = T.V_Multi(i) / Renorm_factor;
                            F_n22 = T.F_Multi(i) / Renorm_factor;

                            Al_IV = 8 - Si_n22;
                            if(Al_IV < 0), Al_IV = 0; end
                            Al_VI = Al_n22 - Al_IV;
                            if(Al_VI < 0), Al_VI = 0; end
                            Mgy = 6 - (Al_VI + Fe_n22 + Zn_n22 + Cr_n22 + V_n22);
                            if(Mgy < 0), Mgy = 0; end
                            Mgx = Mg_n22 - Mgy;
                            if(Mgx < 0), Mgx = 0; end
                            OH = 4 - (F_n22+Cl_n22);
                            if(OH < 0), OH = 0; end

                            if((Mgx + Ca_n22 + Na_n22)/(Mgy + Fe_n22 + Al_VI +Zn_n22 + Cr_n22 + V_n22) > 0.7/6 - (0.7/6)*t && ...
                                    (Mgx + Ca_n22 + Na_n22)/(Mgy + Fe_n22 + Al_VI +Zn_n22 + Cr_n22 + V_n22) < 0.85/6 + (0.85/6)*t)
                                if((Mgy + Fe_n22 + Al_VI +Zn_n22 + Cr_n22 + V_n22)/(Si_n22+Al_IV) > 6/8 - (6/8)*t && (Mgy + Fe_n22 + Al_VI +Zn_n22 + Cr_n22 + V_n22)/(Si_n22+Al_IV) < 6/8 + (6/8)*t)
                                    if(Mgx > 0.4-0.4*t && Mgx <= 0.7+0.7*t)
                                        T.species(i) = {'vermiculite'};
                                        [sorted_formula1, total1] = Formula_Output({Mgy, Fe_n22, Al_VI, Cr_n22, Ti_n22}, {'Mg','Fe','AlVI','Cr','Ti'});
                                        [sorted_formula2, total2] = Formula_Output({Si_n22, Al_IV},{'Si','AlIV'});
                                        [sorted_formula3, total3] = Formula_Output({OH, F_n22, Cl_n22},{'OH','F','Cl'});
                                        T.formula(i) = {strcat('Mg_', num2str(round(Mgx)), sorted_formula1, total1, sorted_formula2, total2, 'O22', sorted_formula3, total3,'·8H2O')};
                                    end
                                end
                            end

                        end
                    end

                end %end of Trioctahderal smectites

            end 
        end   %end of smectites

        %% 7.8.6 Palygorskite-Sepiolite group
        if(ABCDT_n/24 > 6/10.5 - (6/10.5)*t && ABCDT_n/24 < 6/10.5 + (6/10.5)*t && S_n <= 0.05 && P_n <= 0.05)
            if((Mg_n + Fe_n + Al_n + Mn_n + Cr_n + Ti_n + Ca_n + Na_n + K_n) / Si_n > 0.5-0.5*t && (Mg_n + Fe_n + Al_n + Mn_n + Cr_n + Ti_n + Ca_n + Na_n + K_n) / Si_n < 0.5+0.5*t)
                if(Mg_n > Fe_n && Mg_n > Al_n && Al_n > Fe_n && S_n+P_n+Ca_n < 0.1)
                    renorm = 10.5 / 24;
                    T.group2(i) = {'Phyllosilicate'};
                    T.group3(i) = {'Palygorskite-Sepiolite'};
                    T.group4(i) = {''};
                    T.species(i) = {'palygorskite'};
                    Al_IV = 4 - Si_n*renorm;
                    if(Al_IV < 0), Al_IV = 0; end
                    Al_VI = Al_n*renorm - Al_IV;
                    if(Al_VI < 0), Al_VI = 0; end
                    OH = 1 - (F_n*renorm + Cl_n*renorm);
                    if(OH<=0), OH = 0; end
                    [sorted_formula1, total1] = Formula_Output({Mg_n*renorm, Fe_n*renorm, Mn_n*renorm, Al_VI, Ca_n*renorm, Na_n*renorm, K_n*renorm},{'Mg','Fe','Mn','AlVI','Ca','Na','K'});
                    [sorted_formula2, total2] = Formula_Output({Si_n*renorm, Al_IV},{'Si','AlIV'});
                    [sorted_formula3, total3] = Formula_Output({OH, F_n*renorm, Cl_n*renorm},{'OH','F','Cl'});
                    T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O10', sorted_formula3, total3, '·4H2O')};
                end
            end
        end

        if(ABCDT_n > 10*(24/16) - 10*(24/16)*t && ABCDT_n < 10*(24/16) + 10*(24/16)*t && S_n <= 0.05 && P_n <= 0.05)
            if(Si_n / 24 > 6/16 - (6/16)*t && Si_n / 24 < 6/16 + (6/16)*t)
                if((Na_n + Ca_n + K_n) / 24 < 0.02)
                    Renorm_factor = 16/24;
                    Al_IV = 6 - Si_n*Renorm_factor;
                    if(Al_IV < 0), Al_IV = 0; end
                    Al_VI = Al_n*Renorm_factor - Al_IV;
                    if(Al_VI < 0), Al_VI = 0; end

                    if((Mg_n*Renorm_factor + Fe_n*Renorm_factor + Mn_n*Renorm_factor + Cr_n*Renorm_factor + Ti_n*Renorm_factor + Al_VI) / (Si_n*Renorm_factor + Al_IV) > (2/3)-(2/3)*t && (Mg_n*Renorm_factor + Fe_n*Renorm_factor + Mn_n*Renorm_factor + Cr_n*Renorm_factor + Ti_n*Renorm_factor + Al_VI) / (Si_n*Renorm_factor + Al_IV) < (2/3)+(2/3)*t)
                        if(Al_IV / Si_n < 0.1 && (Fe_n+Mg_n)/Si_n > 4/6 - (4/6)*t && (Fe_n+Mg_n)/Si_n < 4/6 + (4/6)*t)
                            if(Mg_n > Fe_n)
                                T.group2(i) = {'Phyllosilicate'};
                                T.group3(i) = {'Palygorskite-Sepiolite'};
                                T.group4(i) = {''};
                                T.species(i) = {'sepiolite'};

                                OH = 2 - (F_n*Renorm_factor+Cl_n*Renorm_factor);
                                if(OH<=0), OH=0; end
                                [sorted_formula1, total1] = Formula_Output({Mg_n*Renorm_factor, Fe_n*Renorm_factor, Mn_n*Renorm_factor, Al_VI, Cr_n*Renorm_factor, Ti_n*Renorm_factor},{'Mg','Fe','Mn','AlVI','Cr','Ti'});
                                [sorted_formula2, total2] = Formula_Output({Si_n*Renorm_factor, Al_IV},{'Si','AlIV'});
                                [sorted_formula3, total3] = Formula_Output({OH, F_n*Renorm_factor, Cl_n*Renorm_factor},{'OH','F','Cl'});
                                T.formula(i) = {strcat(sorted_formula1, total1, sorted_formula2, total2, 'O15', sorted_formula3, total3, '·6H2O')};
                            end
                            if(Fe_n > Mg_n && Mn_n < 0.05)
                                T.group2(i) = {'Phyllosilicate'};
                                T.group3(i) = {'Palygorskite-Sepiolite'};
                                T.group4(i) = {''};
                                T.species(i) = {'ferrisepiolite'};
                                T.formula(i) = {'(Fe3+,Fe2+,Mg)4[(Si,Fe3+)6O15](O,OH)2·6H2O'};
                            end
                           
                        end
                    end
                end
            end
        end

    end %end of Silicates

end %end of Group and Species assignments

%% Check Results
for i = 1:size(T,1)
    %SiO2 Phase check
    if(strcmp(T.group1(i), 'SiO2 Phase'))
        T.ABCDT_ideal(i) = 1;
        T.ABCD_ideal(i) = 0;
        T.T_ideal(i) = 1;
        T.T_prime(i) = (T.Si_Multi(i)/12 + T.Al_Multi(i)/12 + T.Ti_Multi(i)/12 + T.Cr_Multi(i)/12 + T.Fe_Multi(i)/12 + T.Mg_Multi(i)/12);
        if(strcmp(T.species(i), 'Obsidian glass'))
            T.ABCDT_ideal(i) = 1.05;
            T.T_ideal(i) = 1.05;
        end
        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)/12))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)/12))*100;
        T.Delta_T(i) = (T.T_ideal(i)/T.T_prime(i))*100;
    end

    %Carbonates check
    if(strcmp(T.group2(i), 'Monoclinic Carbonate')||strcmp(T.group2(i), 'Trigonal Carbonate')|| strcmp(T.group2(i), 'Orthorhombic Carbonate'))
        if(strcmp(T.group2(i), 'Trigonal Carbonate') || strcmp(T.group2(i), 'Orthorhombic Carbonate'))
            T.ABCDT_ideal(i) = 2;   T.ABCD_ideal(i) = 2;
            if(strcmp(T.species(i), 'witherite'))
                T.ABCDTS(i) = T.ABCDTS(i) + T.Ba_Multi(i); T.ABCDS(i) = T.ABCDS(i) + T.Ba_Multi(i);
                T.ABCDT_ideal(i) = 1*(2/1);   T.ABCD_ideal(i) = 1*(2/1);
            end
        end
        if(strcmp(T.group2(i), 'Monoclinic Carbonate'))
            T.ABCDT_ideal(i) = 2*2;   T.ABCD_ideal(i) = 2*2;
            if(strcmp(T.species(i), 'natron'))
                T.ABCDT_ideal(i) = 1.8*2;   T.ABCD_ideal(i) = 1.8*2;
            end
        end
        T.T_ideal(i) = 0;
        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/((T.ABCDTS(i)/12)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/((T.ABCDS(i)/12)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/T_n2)*100;
    end

    %Sulfates check
    if(strcmp(T.group1(i), 'Sulfate'))
        if(strcmp(T.group4(i), 'Ca-Sulfate'))
            T.ABCDT_ideal(i) = 2;   T.ABCD_ideal(i) = 2;
        elseif(strcmp(T.species(i), 'ivsite'))
            T.ABCDT_ideal(i) = 5*(4/8);   T.ABCD_ideal(i) =5*(4/8);
        elseif(strcmp(T.group4(i), 'Na-Sulfate'))
            T.ABCDT_ideal(i) = 3;   T.ABCD_ideal(i) = 3;
        elseif(strcmp(T.group2(i), 'Mg-Sulfate-H2O'))
            T.ABCDT_ideal(i) = 2;   T.ABCD_ideal(i) = 2;
            if(strcmp(T.species(i), 'caminite'))
                T.ABCDT_ideal(i) = 12*(4/22);   T.ABCD_ideal(i) = 12*(4/22);
            elseif(strcmp(T.species(i), 'mooreite'))
                T.ABCDT_ideal(i) = 17*(4/21);   T.ABCD_ideal(i) = 17*(4/21);
            end
        elseif(strcmp(T.species(i), 'D`ansite'))
            T.ABCDT_ideal(i) = 32*(4/41.5);   T.ABCD_ideal(i) = 32*(4/41.5);
        elseif(strcmp(T.species(i), 'vanthoffite'))
            T.ABCDT_ideal(i) = 11*(4/16);   T.ABCD_ideal(i) = 11*(4/16);
        elseif(strcmp(T.species(i), 'blödite'))
            T.ABCDT_ideal(i) = 5*(4/8);   T.ABCD_ideal(i) = 5*(4/8);
        elseif(strcmp(T.group3(i), 'Fe-Sulfate'))
            if(strcmp(T.species(i), 'rhomboclase'))
                T.ABCDT_ideal(i) = 3.5*(4/8);   T.ABCD_ideal(i) = 3.5*(4/8);
            elseif(strcmp(T.species(i), 'mikasaite/coquimbite'))
                T.ABCDT_ideal(i) = 5.5*(4/12);   T.ABCD_ideal(i) = 5.5*(4/12);
            elseif(strcmp(T.species(i), 'Fe-Sulfate (Fe/S = 1, many species possible)') || strcmp(T.species(i), 'siderotil/szomolnokite') || strcmp(T.species(i), 'melanterite'))
                T.ABCDT_ideal(i) = 2;   T.ABCD_ideal(i) = 2;
            elseif(strcmp(T.species(i), 'römerite/bílinite'))
                T.ABCDT_ideal(i) = 7.5*(4/16);   T.ABCD_ideal(i) = 7.5*(4/16);
            elseif(strcmp(T.species(i), 'copiapite/ferricopiapite') || strcmp(T.species(i), 'magnesiocopiapite'))
                T.ABCDT_ideal(i) = 11.67*(4/25);   T.ABCD_ideal(i) = 11.67*(4/25);
            elseif(strcmp(T.species(i), 'hydroniumjarosite'))
                T.ABCDT_ideal(i) = 6*(4/11);   T.ABCD_ideal(i) = 6*(4/11);
            elseif(strcmp(T.species(i), 'schwertmannite'))
                T.ABCDT_ideal(i) = 24.5*(4/33.6);   T.ABCD_ideal(i) = 24.5*(4/33.6);
            elseif(strcmp(T.species(i), 'viaeneite'))
                T.ABCDT_ideal(i) = 12*(4/27);   T.ABCD_ideal(i) = 12*(4/27);
            elseif(strcmp(T.species(i), 'valleriite') || strcmp(T.species(i), 'ferrovalleriite'))
                T.ABCDT_ideal(i) = 3.53*(4/6.225);   T.ABCD_ideal(i) = 3.53*(4/6.225);
            elseif(strcmp(T.species(i), 'botryogen'))
                T.ABCDT_ideal(i) = 4.25*(4/8.5);   T.ABCD_ideal(i) = 4.25*(4/8.5);
            end
%         elseif(strcmp(T.group4(i), 'Fe-Sulfide'))
%             if(strcmp(T.species(i), 'pyrrhotite'))
%                 T.ABCDT_ideal(i) = 2*(4/4);   T.ABCD_ideal(i) = 2*(4/4);
%             elseif(strcmp(T.species(i), 'pyrite'))
%                 T.ABCDT_ideal(i) = 3.5*(4/8);   T.ABCD_ideal(i) = 3.5*(4/8);
%             elseif(strcmp(T.species(i), 'chalcopyrite'))
%                 T.ABCDT_ideal(i) = 4*(4/8);   T.ABCD_ideal(i) = 4*(4/8);
%             end
        elseif(strcmp(T.species(i), 'philoxenite'))
            T.ABCDT_ideal(i) = 18.25*(4/32);   T.ABCD_ideal(i) = 18.25*(4/32);
        elseif(strcmp(T.species(i), 'jarosite'))
            T.ABCDT_ideal(i) = 7*(4/11);   T.ABCD_ideal(i) = 7*(4/11);
        elseif(strcmp(T.group3(i), 'Ca-Fe-Sulfate'))
            if(strcmp(T.species(i), 'vyalsovite'))
                T.ABCDT_ideal(i) = 4*(4/6.51);   T.ABCD_ideal(i) = 4*(4/6.51);
            elseif(strcmp(T.species(i), 'sturmanite'))
                T.ABCDT_ideal(i) = 10.50*(4/16.5);   T.ABCD_ideal(i) = 10.50*(4/16.5);
            elseif(strcmp(T.species(i), 'calciocopiapite'))
                T.ABCDT_ideal(i) = 11.75*(4/25);   T.ABCD_ideal(i) = 11.75*(4/25);
            end
        elseif(strcmp(T.group3(i), 'Na-K-Fe-Sulfate'))
            if(strcmp(T.species(i), 'eldfellite/amarillite') || strcmp(T.species(i), 'erdite') || strcmp(T.species(i), 'eldfellite/amarillite/erdite'))
                T.ABCDT_ideal(i) = 4.25*(4/8);   T.ABCD_ideal(i) = 4.25*(4/8);
            elseif(strcmp(T.species(i), 'ferrinatrite'))
                T.ABCDT_ideal(i) = 7.25*(4/12);   T.ABCD_ideal(i) = 7.25*(4/12);
            elseif(strcmp(T.species(i), 'adranosite-(Fe)'))
                T.ABCDT_ideal(i) = 7*(4/15);   T.ABCD_ideal(i) = 7*(4/15);
            elseif(strcmp(T.species(i), 'coyoteite'))
                T.ABCDT_ideal(i) = 9*(4/10)*(10/18);   T.ABCD_ideal(i) = 9*(4/10)*(10/18);
            elseif(strcmp(T.species(i), 'D`ansite-(Fe)'))
                T.ABCDT_ideal(i) = 32*(4/41.5);   T.ABCD_ideal(i) = 32*(4/41.5);
            elseif(strcmp(T.species(i), 'ferrotychite'))
                T.ABCDT_ideal(i) = 9*(4/8);   T.ABCD_ideal(i) = 9*(4/8);
            elseif(strcmp(T.species(i), 'metasideronatrite'))
                T.ABCDT_ideal(i) = 5.25*(4/8.5);   T.ABCD_ideal(i) = 5.25*(4/8.5);
            elseif(strcmp(T.species(i), 'metavoltine'))
                T.ABCDT_ideal(i) = 30*(4/50);   T.ABCD_ideal(i) = 30*(4/50);
            elseif(strcmp(T.species(i), 'natrojarosite'))
                T.ABCDT_ideal(i) = 6.75*(4/11);   T.ABCD_ideal(i) = 6.75*(4/11);
            elseif(strcmp(T.species(i), 'nikischerite'))
                T.ABCDT_ideal(i) = 12*(4/17);   T.ABCD_ideal(i) = 12*(4/17);
            elseif(strcmp(T.species(i), 'therasiaite'))
                T.ABCDT_ideal(i) = 8.25*(4/13);   T.ABCD_ideal(i) = 8.25*(4/13);
            elseif(strcmp(T.species(i), 'ungemachite or clinoungemachite'))
                T.ABCDT_ideal(i) = 18.25*(4/25);   T.ABCD_ideal(i) = 18.25*(4/25);
            elseif(strcmp(T.species(i), 'yavapaiite'))
                T.ABCDT_ideal(i) = 4.25*(4/8);   T.ABCD_ideal(i) = 4.25*(4/8);
            elseif(strcmp(T.species(i), 'syngenite'))
                T.ABCDT_ideal(i) = 5*(4/8);   T.ABCD_ideal(i) = 5*(4/8);
            elseif(strcmp(T.species(i), 'görgeyite'))
                T.ABCDT_ideal(i) = 13*(4/24);   T.ABCD_ideal(i) = 13*(4/24);
            end
        elseif(strcmp(T.species(i), 'kainite'))
            T.ABCDT_ideal(i) = 3*(4/4.5);   T.ABCD_ideal(i) = 3*(4/4.5);
        elseif(strcmp(T.species(i), 'alunite'))
            T.ABCDT_ideal(i) = 6*(4/11);   T.ABCD_ideal(i) = 6*(4/11);
        elseif(strcmp(T.species(i), 'caledonite'))
            T.ABCDT_ideal(i) = 10*(4/16);   T.ABCD_ideal(i) = 10*(4/16);
        elseif(strcmp(T.group3(i), 'Cu-Sulfate'))
            T.ABCDT_ideal(i) = 2*(4/4);   T.ABCD_ideal(i) =  2*(4/4);
        end

        T.T_ideal(i) = 0;
        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDTS(i)/6))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCDS(i)/6))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T(i)/6))*100;
    end

    %Phosphate check
    if(strcmp(T.group1(i), 'Phosphate'))
        if(strcmp(T.group4(i), 'Apatite'))
            T.ABCDT_ideal(i) = 8*(4/12.5);   T.ABCD_ideal(i) = 8*(4/12.5);
        elseif(strcmp(T.species(i), 'tuite'))
            T.ABCDT_ideal(i) = 5*(4/8);   T.ABCD_ideal(i) = 5*(4/8);
        elseif(strcmp(T.species(i), 'brockite'))
            T.ABCDT_ideal(i) = 2*(4/4);   T.ABCD_ideal(i) = 2*(4/4);
        elseif(strcmp(T.species(i), 'monetite or brushite'))
            T.ABCDT_ideal(i) = 2*(4/3.5);   T.ABCD_ideal(i) = 2*(4/3.5);
        elseif(strcmp(T.species(i), 'isoclasite'))
            T.ABCDT_ideal(i) = 3*(4/4.5);   T.ABCD_ideal(i) = 3*(4/4.5);
        elseif(strcmp(T.group4(i), 'Na-Phosphate-OH'))
            T.ABCDT_ideal(i) = 3*(4/3.5);   T.ABCD_ideal(i) = 3*(4/3.5);
        elseif(strcmp(T.group4(i), 'Na-Phosphate-F'))
            T.ABCDT_ideal(i) = 9*(4/8.5);   T.ABCD_ideal(i) = 9*(4/8.5);
        elseif(strcmp(T.species(i), 'buchwaldite'))
            T.ABCDT_ideal(i) = 3*(4/4);   T.ABCD_ideal(i) = 3*(4/4);
        elseif(strcmp(T.species(i), 'merrillite'))
            T.ABCDT_ideal(i) = 18*(4/28);   T.ABCD_ideal(i) = 18*(4/28);
        elseif(strcmp(T.group2(i), 'Mg-Phosphate'))
            if(strcmp(T.species(i), 'dry: chopinite or ferringtonite, or OH-bearing: barićite, bobierrite, cattiite'))
                T.ABCDT_ideal(i) = 5*(4/8);   T.ABCD_ideal(i) = 5*(4/8);
            elseif(strcmp(T.species(i), 'holtedahlite or hydroxylwagnerite or kovdorskite'))
                T.ABCDT_ideal(i) = 3*(4/4.5);   T.ABCD_ideal(i) = 3*(4/4.5);
            elseif(strcmp(T.species(i), 'newberyite or phosphorrösslerite'))
                T.ABCDT_ideal(i) = 2*(4/3.5);   T.ABCD_ideal(i) = 2*(4/3.5);
            elseif(strcmp(T.species(i), 'raaedite'))
                T.ABCDT_ideal(i) = 9*(4/12);   T.ABCD_ideal(i) = 9*(4/12);
            end
        elseif(strcmp(T.group2(i), 'Fe-Phosphate'))
            if(strcmp(T.species(i), 'graftonite or sarcopside'))
                T.ABCDT_ideal(i) = 5*(4/8);   T.ABCD_ideal(i) = 5*(4/8);
            elseif(strcmp(T.species(i), 'grattarolaite'))
                T.ABCDT_ideal(i) = 5*(4/7);   T.ABCD_ideal(i) = 5*(4/7);
            elseif(strcmp(T.species(i), 'heterosite or rodolicoite'))
                T.ABCDT_ideal(i) = 2.25*(4/4);   T.ABCD_ideal(i) = 2.25*(4/4);
            end
        elseif(strcmp(T.species(i), 'archerite'))
            T.ABCDT_ideal(i) = 1.85*(4/3);   T.ABCD_ideal(i) = 1.85*(4/3);
        elseif( (strcmp(T.group2(i), 'REE-Phosphate') || strcmp(T.species(i), 'xenotime-Y')) && ~strcmp(T.species(i), 'monazite-Nd'))
            T.ABCDT_ideal(i) = 2*(4/4);   T.ABCD_ideal(i) = 2*(4/4);
        elseif(strcmp(T.species(i), 'monazite-Nd'))
            T.ABCDT_ideal(i) = 2*(4/4);   T.ABCD_ideal(i) = 1.8*(4/4);
        end

        T.T_ideal(i) = 0;
        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDTP(i)/6))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCDP(i)/6))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T(i)/6))*100;
    end


    %%Spinels Check
    if(strcmp(T.group4(i), 'Spinel'))
        T.ABCDT_ideal(i) = 3.02*(4/4);   T.ABCD_ideal(i) = 3.02*(4/4);
        if(strcmp(T.species(i), 'Ti-magnetite'))
            T.ABCDT_ideal(i) = 3.5*(4/4);   T.ABCD_ideal(i) = 3.5*(4/4);
        end
        if(strcmp(T.species(i), 'magnetite'))
            T.ABCDT_ideal(i) = 3.9*(4/4);   T.ABCD_ideal(i) = 3.9*(4/4);
        end
        if(strcmp(T.species(i), 'eskolaite'))
            T.ABCDT_ideal(i) = 1.95*(4/3);   T.ABCD_ideal(i) = 1.95*(4/3);
        end
        if(strcmp(T.species(i), 'ulvöspinel'))
            T.ABCDT_ideal(i) = 3.1*(4/4);   T.ABCD_ideal(i) = 3.1*(4/4);
        end

        T.T_ideal(i) = 0;
        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)/6))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)/6))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T(i)/6))*100;
    end

    %%Ilmentie check
    if(strcmp(T.group4(i), 'Ilmenite'))
        T.ABCDT_ideal(i) = 2*(3/3);   T.ABCD_ideal(i) = 2*(3/3);
        T.T_ideal(i) = 0;
        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)/8))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)/8))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T(i)/6))*100;
    end

    %Oxide-Hydroxide check (minus spinels and ilmenties)
    if( (strcmp(T.group2(i), 'Oxide') || strcmp(T.group2(i), 'Hydroxide') || strcmp(T.group2(i), 'Oxide or Hydroxide')) && ~strcmp(T.group4(i), 'Spinel') && ~strcmp(T.group4(i), 'Ilmenite'))
        if(strcmp(T.species(i), 'manganosite or manganite') || strcmp(T.species(i), 'pyrochroite') || strcmp(T.group4(i), 'Mg-Oxide or Mg-Hydroxide') || strcmp(T.species(i), 'brucite') || strcmp(T.species(i), 'periclase')|| strcmp(T.species(i), 'periclase or brucite')|| strcmp(T.species(i), 'amakinite or wüstite') || strcmp(T.species(i), 'wüstite' ))
            T.ABCDT_ideal(i) = 1*(24/1);   T.ABCD_ideal(i) = 1*(24/1);
        elseif(strcmp(T.species(i), 'bixbyite'))
            T.ABCDT_ideal(i) = 2.5*(24/3);   T.ABCD_ideal(i) = 2.5*(24/3);
        elseif(strcmp(T.species(i), 'pyrolusite') || strcmp(T.species(i), 'pyrolusite/bixbyite') )
            T.ABCDT_ideal(i) = 2*(24/2);   T.ABCD_ideal(i) = 2*(24/2);
        elseif(strcmp(T.species(i), 'hausmannite'))
            T.ABCDT_ideal(i) = 2.5*(24/4);   T.ABCD_ideal(i) = 2.5*(24/4);
        elseif(strcmp(T.species(i), 'birnessite'))
            T.ABCDT_ideal(i) = 3.95*(24/4);   T.ABCD_ideal(i) = 3.95*(24/4);
        elseif(strcmp(T.group4(i), 'Ti-Oxide'))
            T.ABCDT_ideal(i) = 1*(24/2);   T.ABCD_ideal(i) = 1*(24/2);
        elseif(strcmp(T.species(i), 'corundum'))
            T.ABCDT_ideal(i) = 2*(24/3);   T.ABCD_ideal(i) = 2*(24/3);
        elseif(strcmp(T.species(i), 'aurorite'))
            T.ABCDT_ideal(i) = 6.5*(24/7);   T.ABCD_ideal(i) = 6.5*(24/7);
        elseif(strcmp(T.species(i), 'vernadite'))
            T.ABCDT_ideal(i) = 1.03*(24/1);   T.ABCD_ideal(i) = 1.02*(24/1);
        elseif(strcmp(T.species(i), 'janggunite'))
            T.ABCDT_ideal(i) = 8*(24/11);   T.ABCD_ideal(i) = 8*(24/11);
        elseif(strcmp(T.species(i), 'takanelite'))
            T.ABCDT_ideal(i) = 8.5*(24/9);   T.ABCD_ideal(i) = 8.25*(24/9);
        elseif(strcmp(T.species(i), 'hollandite or romanèchite'))
            T.ABCDT_ideal(i) = 13.75*(24/16);   T.ABCD_ideal(i) = 13.75*(24/16);
        elseif(strcmp(T.species(i), 'akaganeite'))
            T.ABCDT_ideal(i) = 8.75*(24/8.75);   T.ABCD_ideal(i) = 8.75*(24/8.75);
        end

        T.T_ideal(i) = 0;
        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T(i)))*100;
    end

    %%Halides check
    if(strcmp(T.group1(i), 'Halide'))
        if(strcmp(T.species(i), 'halite') || strcmp(T.species(i), 'sylvite'))
            T.ABCDT_ideal(i) = 1*(24/0.5);   T.ABCD_ideal(i) = 1*(24/0.5);
        elseif(strcmp(T.species(i), 'chloromagnesite') || strcmp(T.species(i), 'lawrencite') || strcmp(T.species(i), 'antarcticite') || strcmp(T.species(i), 'scacchite'))
            T.ABCDT_ideal(i) = 1*(24/1);   T.ABCD_ideal(i) = 1*(24/1);
        elseif(strcmp(T.species(i), 'molysite'))
            T.ABCDT_ideal(i) = 1.5*(24/1.5);   T.ABCD_ideal(i) = 1.5*(24/1.5);
        end

        T.T_ideal(i) = 0;
        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T(i)))*100;
    end

    %%Feldspars check
    if(strcmp(T.group3(i), 'Feldspar'))
        T.T_prime(i) = T.Si_Multi(i)/3 + T.Al_Multi(i)/3 + T.Fe_Multi(i)/3  + T.Mg_Multi(i)/3 ;
        T.ABCDT_ideal(i) = 5*(24/8);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 1;

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Feldspathoids check
    if(strcmp(T.group3(i), 'Feldspathoid'))
        if(strcmp(T.species(i), 'sodalite'))
            T.T_prime(i) = T.Si_Multi(i)*(12.5/24) + T.Al_Multi(i)*(12.5/24) + T.Fe_Multi(i)*(12.5/24) + T.Mg_Multi(i)*(12.5/24);
            T.ABCDT_ideal(i) = 10*(24/12.5);   T.T_ideal(i) = 6;  T.ABCD_ideal(i) = 4;
        end

        if(strcmp(T.species(i), 'leucite'))
            T.T_prime(i) = T.Si_Multi(i)/4 + T.Al_Multi(i)/4 + T.Fe_Multi(i)/4 + T.Mg_Multi(i)/4;
            T.ABCDT_ideal(i) = 4*(24/6);   T.T_ideal(i) = 3;  T.ABCD_ideal(i) = 1;
        end
        if(strcmp(T.species(i), 'cancrinite'))
            T.T_prime(i) = T.Si_Multi(i) + T.Al_Multi(i) + T.Fe_Multi(i) + T.Mg_Multi(i);
            T.ABCDT_ideal(i) = 18.5*(24/24);   T.T_ideal(i) = 11.5;  T.ABCD_ideal(i) = 8;
        end
        if(strcmp(T.species(i), 'nepheline'))
            T.T_prime(i) = T.Si_Multi(i)*(16/24) + T.Al_Multi(i)*(16/24) + T.Fe_Multi(i)*(16/24) + T.Mg_Multi(i)*(16/24);
            T.ABCDT_ideal(i) = 11.8*(24/16);   T.T_ideal(i) = 8;  T.ABCD_ideal(i) = 4;
        end
        if(strcmp(T.species(i), 'trinepheline'))
            T.T_prime(i) = T.Si_Multi(i)*(16/24) + T.Al_Multi(i)*(16/24) + T.Fe_Multi(i)*(16/24) + T.Mg_Multi(i)*(16/24);
            T.ABCDT_ideal(i) = 3*(24/4);   T.T_ideal(i) = 8;  T.ABCD_ideal(i) = 3;
        end

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Zeolites check
    if(strcmp(T.group3(i), 'Zeolite'))
        T.T_prime(i) = T.Si_Multi(i)*(72/24) + T.Al_Multi(i)*(72/24) + T.Fe_Multi(i)*(72/24) + T.Mg_Multi(i)*(72/24);

        if(strcmp(T.species(i), 'clinoptilolite-Ca'))
            T.ABCDT_ideal(i) = 39*(24/72);   T.T_ideal(i) = 36;  T.ABCD_ideal(i) = 3;
        end

        if(strcmp(T.species(i), 'clinoptilolite-Na') || strcmp(T.species(i), 'clinoptilolite-K') || strcmp(T.species(i), 'heulandite-Na'))
            T.ABCDT_ideal(i) = 42*(24/72);   T.T_ideal(i) = 36;  T.ABCD_ideal(i) = 6;
        end

        if(strcmp(T.species(i), 'heulandite-Ca') || strcmp(T.species(i), 'heulandite-K') || strcmp(T.species(i), 'stilbite-Ca') )
            T.ABCDT_ideal(i) = 41*(24/72);   T.T_ideal(i) = 36;  T.ABCD_ideal(i) = 5;
        end

        if(strcmp(T.species(i), 'stilbite-Na'))
            T.ABCDT_ideal(i) = 45*(24/72);   T.T_ideal(i) = 36;  T.ABCD_ideal(i) = 9;
        end

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Olivines check
    if(strcmp(T.group3(i), 'Olivine'))
        T.T_prime(i) = T.Si_Multi(i)/6 + T.Al_Multi(i)/6 + T.Ti_Multi(i)/6 + T.Cr_Multi(i)/6;
        
        T.ABCDT_ideal(i) = 3*(24/4);   T.T_ideal(i) = 1;  T.ABCD_ideal(i) = 2;

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Garnets check
    if(strcmp(T.group3(i), 'Garnet'))
        Al_IV = 3 - T.Si_Multi(i)/2; 
        if(Al_IV < 0), Al_IV = 0; end       

        T.T_prime(i) = T.Si_Multi(i)/2 + Al_IV;
        T.ABCDT_ideal(i) = 8*(24/12);   T.T_ideal(i) = 3;  T.ABCD_ideal(i) = 5;
        if(strcmp(T.group4(i), 'Pyralspite Garnet') || strcmp(T.species(i), 'grossular') || strcmp(T.species(i), 'uvarovite'))
            T.T_prime(i) = T.Si_Multi(i)/2 + Al_IV;
            T.ABCDT_ideal(i) = 8*(24/12);   T.T_ideal(i) = 3;  T.ABCD_ideal(i) = 5;
        end
        if(strcmp(T.species(i), 'andradite') || strcmp(T.species(i), 'skiagite') || strcmp(T.species(i), 'calderite'))
            T.T_prime(i) = T.Si_Multi(i)/2 + Al_IV;
            T.ABCDT_ideal(i) = 8.5*(24/12);   T.T_ideal(i) = 3.2;  T.ABCD_ideal(i) = 5;
        end

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Epidotes check
    if(strcmp(T.group4(i), 'Epidote'))
        Al_IV = 3 - T.Si_Multi(i)*(12.5/24); 
        if(Al_IV < 0), Al_IV = 0; end 

        T.T_prime(i) = T.Si_Multi(i)*(12.5/24) + Al_IV;

        T.ABCDT_ideal(i) = 8*(24/12.5); T.T_ideal(i) = 3;

        if(strcmp(T.species(i), 'epidote or pumpelleyite-Fe') || strcmp(T.species(i), 'tweddillite') || strcmp(T.species(i), 'piemontite or pumpelleyite-Mn'))
            T.ABCDT_ideal(i) = 8.4*(24/12.5); T.T_ideal(i) = 3.2;
        end

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Al-silicates check
    if(strcmp(T.group3(i), 'Al-Silicate'))
        T.T_prime(i) = T.Si_Multi(i)*(5/24) + Al_IV;
        T.ABCDT_ideal(i) = 3*(24/5);   T.T_ideal(i) = 1;  T.ABCD_ideal(i) = 2;

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Zircon check
    if(strcmp(T.species(i), 'zircon'))
        Al_IV = 1 - T.Si_Multi(i)/6; 
        if(Al_IV < 0), Al_IV = 0; end
     
        T.T_prime(i) = T.Si_Multi(i)/6 + Al_IV;
        T.ABCDT_ideal(i) = 2*(24/4);   T.T_ideal(i) = 1;  T.ABCD_ideal(i) = 2;

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    if(strcmp(T.species(i), 'baddeleyite'))
        T.ABCDT_ideal(i) = 1*(24/2);   T.T_ideal(i) = 0;  T.ABCD_ideal(i) = 1*(24/2);
        
        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Pyroxenes check
    if(strcmp(T.group3(i), 'Pyroxene'))
        Al_IV = 2 - T.Si_Multi(i)/4;
        if(Al_IV < 0), Al_IV = 0; end
 
        T.T_prime(i) = T.Si_Multi(i)/4 + Al_IV;
        T.ABCDT_ideal(i) = 4*(24/6);   T.T_ideal(i) = 2;  T.ABCD_ideal(i) = 2;

        if(strcmp(T.species(i), 'aegirine'))
            T.ABCDT_ideal(i) = 4.3*(24/6);   T.T_ideal(i) = 2.15;  T.ABCD_ideal(i) = 4.25-2;
        end

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Amphiboles check
    if(strcmp(T.group3(i), 'Amphibole'))

        Renorm_factor = 24 / 23;

        Si_n23 = T.Si_Multi(i) / Renorm_factor;

        Al_IV = 8 - Si_n23; 
        if(Al_IV < 0), Al_IV = 0; end        

        T.T_prime(i) = Si_n23 + Al_IV;
        T.ABCDT_ideal(i) = 15.5*(24/23);   T.T_ideal(i) = 8; 

        if(strcmp(T.species(i), 'potassic-richterite') || strcmp(T.species(i), 'richterite') || strcmp(T.species(i), 'ferro-richterite'))
            T.ABCDT_ideal(i) = 16*(24/23);   T.T_ideal(i) = 8;
        end
        if(strcmp(T.species(i), 'magnesio-hastingsite') || strcmp(T.species(i), 'arfvedsonite'))
            T.ABCDT_ideal(i) = 15.75*(24/23);   T.T_ideal(i) = 8;  
        end
        if(strcmp(T.species(i), 'edenite') || strcmp(T.species(i), 'ferro-edenite') || strcmp(T.species(i), 'hastingsite'))
            T.ABCDT_ideal(i) = 16*(24/23);   T.T_ideal(i) = 8;
        end
        if(strcmp(T.species(i), 'glaucophane'))
            T.ABCDT_ideal(i) = 15.25*(24/23);   T.T_ideal(i) = 8;
        end
        if(strcmp(T.species(i), 'riebeckite'))
            T.ABCDT_ideal(i) = 15.75*(24/23);   T.T_ideal(i) = 8.3;
        end
        if(strcmp(T.species(i), 'pargasite'))
            T.ABCDT_ideal(i) = 15.8*(24/23);   T.T_ideal(i) = 8;
        end

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Micas check
    if(strcmp(T.group3(i), 'Mica'))
        Al_IV = 4 - T.Si_Multi(i)*(11/24); 
        if(Al_IV < 0), Al_IV = 0; end        

        T.T_prime(i) = T.Si_Multi(i)*(11/24) + Al_IV;
        T.T_ideal(i) = 4;

        if(strcmp(T.species(i), 'biotite (annite)') || strcmp(T.species(i), 'biotite (phlogopite)') || strcmp(T.species(i), 'siderophyllite'))
            T.ABCDT_ideal(i) = 8*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 4;
        end
        if(strcmp(T.species(i), 'ferroaluminoceladonite') || strcmp(T.species(i), 'aluminoceladonite')|| strcmp(T.species(i), 'muscovite') || strcmp(T.species(i), 'paragonite') || strcmp(T.species(i), 'margarite'))
            T.ABCDT_ideal(i) = 7*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 4;
        end
        if(strcmp(T.species(i), 'celadonite'))
            T.ABCDT_ideal(i) = 7.25*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 4;
        end
        if(strcmp(T.species(i), 'illite'))
            T.ABCDT_ideal(i) = 6.75*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 4;
        end

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Prehnite check
    if(strcmp(T.group3(i), 'Prehnite'))
        Al_IV = 4 - T.Si_Multi(i)*(11/24); 
        if(Al_IV < 0), Al_IV = 0; end

        T.T_prime(i) = T.Si_Multi(i)*(11/24) + Al_IV;

        if(strcmp(T.species(i), 'prehnite'))
            T.ABCDT_ideal(i) = 7*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 4;
        end
        if(strcmp(T.species(i), 'ferriprehnite'))
            T.ABCDT_ideal(i) = 7.25*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 4;
        end

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;

    end

    %%Talc-Pyrophyllite check
    if(strcmp(T.group3(i), 'Talc-Pyrophyllite'))
        Al_IV = 4 - T.Si_Multi(i)*(11/24); 
        if(Al_IV < 0), Al_IV = 0; end        
        
        T.T_prime(i) = T.Si_Multi(i)*(11/24) + Al_IV;

        if(strcmp(T.species(i), 'talc'))
            T.ABCDT_ideal(i) = 7*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 4;
        end
        if(strcmp(T.species(i), 'minnesotaite'))
            T.ABCDT_ideal(i) = 7*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 4;
        end
        if(strcmp(T.species(i), 'pyrophyllite'))
            T.ABCDT_ideal(i) = 6*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 4;
        end

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Chlorites check
    if(strcmp(T.group3(i), 'Chlorite'))
        Al_IV = 4 - T.Si_Multi(i)*(14/24); 
        if(Al_IV < 0), Al_IV = 0; end

        T.T_prime(i) = T.Si_Multi(i)*(14/24) + Al_IV;
        T.ABCDT_ideal(i) = 10*(24/14);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 6;

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Kaolinite-Serpentine check
    if(strcmp(T.group3(i), 'Kaolinite-Serpentine'))
        Si_n7 = T.Si_Multi(i)*(7/24);

        Al_IV = 2 - Si_n7; 
        if(Al_IV < 0), Al_IV = 0; end

        T.T_prime(i) = Si_n7 + Al_IV;

        if(strcmp(T.group4(i), 'Serpentine'))
            T.ABCDT_ideal(i) = 5*(24/7);   T.T_ideal(i) = 2;
        end
        if(strcmp(T.group4(i), 'Kaolinite'))
            T.ABCDT_ideal(i) = 4*(24/7);   T.T_ideal(i) = 2;
        end
        if(strcmp(T.species(i), 'greenalite'))
            T.ABCDT_ideal(i) = 4.85*(24/7);   T.T_ideal(i) = 2.07;
        end
        if(strcmp(T.species(i), 'hisingerite'))
            T.ABCDT_ideal(i) = 4.7*(24/7);   T.T_ideal(i) = 2.3;
        end

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Smectites check
    if(strcmp(T.group3(i), 'Smectite'))
        Al_IV = 4 - T.Si_Multi(i)*(11/24); 
        if(Al_IV < 0), Al_IV = 0; end

        T.T_prime(i) = T.Si_Multi(i)*(11/24) + Al_IV;

        if(strcmp(T.species(i), 'montmorillonite') || strcmp(T.species(i), 'beidellite'))
            T.ABCDT_ideal(i) = 6.3*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 3;
        end
        if(strcmp(T.species(i), 'hectorite')|| strcmp(T.species(i), 'saponite') || strcmp(T.species(i), 'sauconite') || strcmp(T.species(i), 'nontronite'))
            T.ABCDT_ideal(i) = 7.3*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 2;
        end
        if(strcmp(T.species(i), 'vermiculite'))
            Al_IV = 8 - T.Si_Multi(i)*(22/24); 
            if(Al_IV < 0), Al_IV = 0; end

            T.T_prime(i) = T.Si_Multi(i)*(22/24) + Al_IV;
            T.ABCDT_ideal(i) = 14.7*(24/22);   T.T_ideal(i) = 8;  T.ABCD_ideal(i) = 2;
        end

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    %%Palygorskite-Sepiolite check
    if(strcmp(T.group3(i), 'Palygorskite-Sepiolite'))
        if(strcmp(T.species(i), 'palygorskite'))
            Al_IV = 4 - T.Si_Multi(i)*(24/10.5); 
            if(Al_IV < 0), Al_IV = 0; end

            T.T_prime(i) = T.Si_Multi(i) + Al_IV;
            T.ABCDT_ideal(i) = 6*(24/10.5);   T.T_ideal(i) = 4*(24/10.5);  T.ABCD_ideal(i) = 2;
        end
        if(strcmp(T.species(i), 'sepiolite') || strcmp(T.species(i), 'ferrisepiolite'))
            Si_n16 = T.Si_Multi(i)*(16/24);

            Al_IV = 6 - Si_n16;
            if(Al_IV < 0), Al_IV = 0; end

            T.T_prime(i) = Si_n16 + Al_IV;
            T.ABCDT_ideal(i) = 10*(24/16);   T.T_ideal(i) = 6;  T.ABCD_ideal(i) = 6;
        end

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    if(strcmp(T.species(i), 'staurolite'))
        T.ABCDT_ideal(i) = 15*(24/23.5);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 11*(24/23.5);

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    if(strcmp(T.species(i), 'cordierite') || strcmp(T.species(i), 'sekaninaite'))
        T.ABCDT_ideal(i) = 11*(24/18);   T.T_ideal(i) = 5;  T.ABCD_ideal(i) = 6*(24/18);

        T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
        T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
        T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;
    end

    ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
    ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
    T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;
  
    if(strcmp(T.species(i), 'Mg Carbonate or Oxide-Hydroxide') || strcmp(T.species(i), 'Fe Oxide-Hydroxide or Carbonate') || strcmp(T.species(i), 'Mn Oxide-Hydroxide or Carbonate') || strcmp(T.species(i), 'Cr Oxide-Hydroxide'))
        T.check(i) = {'no check'};
    elseif( strcmp(T.group2(i), 'Monoclinic Carbonate') || strcmp(T.group2(i), 'Trigonal Carbonate') || strcmp(T.group2(i), 'Orthorhombic Carbonate') ||...
        strcmp(T.group1(i), 'Sulfate')|| strcmp(T.group1(i), 'Phosphate') || strcmp(T.group4(i), 'Spinel') || strcmp(T.group4(i), 'Ilmenite') || ...
        strcmp(T.group2(i), 'Oxide') || strcmp(T.group2(i), 'Hydroxide') || strcmp(T.group2(i), 'Oxide or Hydroxide') ||...
        strcmp(T.group1(i), 'Halide') || strcmp(T.species(i), 'staurolite') || strcmp(T.species(i), 'cordierite') || strcmp(T.species(i), 'sekaninaite'))
        
        if(ABCDT_range == 0 && ABCD_range == 0)
            T.check(i) = {'TRUE'};
        elseif(ABCDT_range == 1 && ABCD_range == 1)
            T.check(i) = {'FALSE'};
        elseif(ABCDT_range == 0 && ABCD_range == 1)
            T.check(i) = {'MAYBE SIMILAR TO'};
        elseif(ABCDT_range == 1 && ABCD_range == 0)
            T.check(i) = {'MAYBE RELATED'};
        end
    
    elseif( strcmp(T.group1(i), 'SiO2 Phase') || strcmp(T.group3(i), 'Feldspar') || strcmp(T.group3(i   ), 'Feldspathoid') || strcmp(T.group3(i), 'Zeolite') ||...
            strcmp(T.group3(i), 'Olivine') || strcmp(T.group3(i), 'Garnet') || strcmp(T.group3(i), 'Al-Silicate') || strcmp(T.group4(i), 'Epidote') || strcmp(T.species(i), 'zircon') ||...
            strcmp(T.group3(i), 'Pyroxene') || strcmp(T.group3(i), 'Amphibole') || strcmp(T.group3(i), 'Mica') || strcmp(T.group3(i), 'Prehnite') || strcmp(T.group3(i), 'Talc-Pyrophyllite') ||...
            strcmp(T.group3(i), 'Chlorite') || strcmp(T.group3(i), 'Kaolinite-Serpentine') || strcmp(T.group3(i), 'Smectite') || strcmp(T.group3(i), 'Palygorskite-Sepiolite') )
    
            if(ABCDT_range == 0 && T_range == 0)
                T.check(i) = {'TRUE'};
            elseif(ABCDT_range == 1 && T_range == 1)
                T.check(i) = {'FALSE'};
            elseif(ABCDT_range == 0 && T_range == 1)
                T.check(i) = {'MAYBE SIMILAR TO'};
            elseif(ABCDT_range == 1 && T_range == 0)
                T.check(i) = {'MAYBE RELATED'};
            end
    
    end

    if(strcmp(T.species(i), ''))
        T.check(i) = {'no species to check'};
    end

    if(strcmp(T.species(i), 'Al Oxide or Hydroxide') || strcmp(T.species(i), 'magnetite/hematite/wüstite/Fe-Hydroxide') || strcmp(T.species(i), 'hisingerite or greenalite') || strcmp(T.species(i), 'cuprite/tenorite/spertiniite'))

        T.check(i) = {'no check'};  
    end 

end

%% Organize Output

T_export = table;

    variablesToExport = {'Point','Input_SiO2','Input_TiO2','Input_Al2O3','Input_Cr2O3','Input_FeO','Input_NiO','Input_MnO','Input_MgO','Input_CaO','Input_Na2O','Input_K2O','Input_P2O5','Input_SO3','Input_F','Input_Cl','Input_V2O3','Input_ZnO','Input_CoO','Input_BaO','Input_SrO','Input_B2O3','Input_PbO','Input_CuO','Input_Sb2O3','Input_As2O5','Input_ThO2','Input_ZrO2','Input_HfO2','Input_Ag2O','Input_Y2O3','Input_La2O3','Input_Ce2O3','Input_Nd2O3','Original_Total',...
        'SiO2','TiO2','Al2O3','Cr2O3','FeO','NiO','MnO','MgO','CaO','Na2O','K2O','P2O5','SO3','F','Cl','V2O3','ZnO','CoO','BaO','SrO','B2O3','PbO','CuO','Sb2O3','As2O5','ThO2','ZrO2','HfO2','Ag2O','Y2O3','La2O3','Ce2O3','Nd2O3','New_Total',...
        'Si_Multi','Ti_Multi','Al_Multi','Cr_Multi','Fe_Multi','Mn_Multi','Mg_Multi','Ni_Multi','Ca_Multi','Na_Multi','K_Multi','S_Multi','P_Multi','V_Multi','Zn_Multi','Co_Multi','Ba_Multi','Sr_Multi','B_Multi','Pb_Multi','Cu_Multi','Sb_Multi','As_Multi','Th_Multi','Zr_Multi','Hf_Multi','Ag_Multi','Y_Multi','La_Multi','Ce_Multi','Nd_Multi',...
        'group1','group2','group3','group4','species','formula','check','Fo1','Fo2','An','Ab','Or','Wo','En','Fs','Aeg','Jd','Di','Alm','Grs','Sps','Prp'};

    columnNames = T.Properties.VariableNames;

    for i = 1:numel(variablesToExport)
        columnName = variablesToExport{i};
        
        % Check if the specified column exists in the table
        if(strcmp(columnName, 'formula'))
            T_export.(columnName) = T.formula;
        elseif(strcmp(columnName, 'check'))
            T_export.(columnName) = T.check;
        elseif ismember(columnName, columnNames)
            columnData = T.(columnName); % Get the column data
            
            if isnumeric(columnData)
                % Check if there is at least one non-zero value in the numeric column
                if any(columnData ~= 0)
                    % Add the column to the export table
                    T_export.(columnName) = columnData;
                end
            elseif iscell(columnData) && all(cellfun(@ischar, columnData))
                % Check if there is at least one non-empty value in the text column
                if any(~cellfun(@isempty, columnData))
                    % Add the column to the export table
                    T_export.(columnName) = columnData;
                end
            end
        else
            
        end
        
    end

    % Remove answers with false and maybe check results
    for i = 1:size(T_export,1)
        if(strcmp(T_export.check(i), 'FALSE') || strcmp(T_export.check(i), 'MAYBE RELATED') || strcmp(T_export.check(i), 'MAYBE SIMILAR TO'))
            T_export.group2{i} = '';
            T_export.group3{i} = '';
            T_export.group4{i} = '';
            T_export.species{i} = '';
            T_export.formula{i} = '';
        end
    end

%% Export final table
    writetable(T_export, 'table_export.xlsx', 'WriteVariableNames', true, 'Sheet', 'Main Results', 'WriteMode','overwrite');    
    writetable(T, 'table_export.xlsx', 'WriteVariableNames', true, 'Sheet','Supplemental Data', 'WriteMode','overwrite');   

end