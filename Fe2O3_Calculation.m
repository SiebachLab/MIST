% Function to incorporate Fe2O3 and recalculate Fe2+ and Fe3+ cations
% Inputs are:
    %1: the original table of data (T)
    %2: the current sample number 'i' (s)
    %3: the current CD (sum of cations) value (CD_initial)
    %4: the current Fe cations (Fe_initial)
    %5: Renormalization factor to be able to normalize to however many
    %oxygens desiered (24/x = y where Renorm_factor = y)
    %6: number to decide which cations to use in recalculation of CD (see
    %bottom of code). 1 = spinels, etc. (CD_cations)
    %7: Ideal_cations (what we want ABCetc. to add up to)

function [T, Fe2, Fe3, CD_renorm_new] = Fe2O3_Calculation(T, s, CD_initial, Fe_initial, Renorm_factor, CD_cations, Ideal_cations)

%If CD < 3, Fe2 will return as previous Fe_n4 value, Fe3 will be 0, and CD
%will be the same as previous
Fe2 = Fe_initial;
Fe3 = 0;
CD_renorm_new = CD_initial;

% Set constants for FeO and Fe2O3
atomic_FeO = 71.846;
atomic_Fe2O3 = 159.6882;

nor_Oxygen = 24;

%If CD >= 3, indication of Fe2 AND Fe3 present and this loop will calculate
%the cations Fe2 and Fe3 based on adding Fe2O3 to our oxide list
if(round(CD_initial,3) > Ideal_cations && T.FeO(s) ~= 0)
    
    FeO = T.FeO(s); %FeO starts as a constant
    w = 0:0.05:FeO; %w ( = Fe2O3) will iteratively increase 
    
    for i = 1:length(w)
        
        T.FeO_Calc(s) = FeO - w(i);
        T.Fe2O3_Calc(s) = (FeO - T.FeO_Calc(s))*1.1; 
        
        FeO_Mole = T.FeO_Calc(s) / atomic_FeO;      %Recalculate FeO moles based on new FeO%
        FeO_Oxygen = FeO_Mole * 1;                  %Realculate oxygen contributed by FeO

        Fe2O3_Mole = T.Fe2O3_Calc(s) / atomic_Fe2O3;    %Calculate Fe2O3 moles based on Fe2O3%
        Fe2O3_Oxygen = Fe2O3_Mole * 3;                  %Calculate oxygen contributed by Fe2O3
        
        %Inital sum of oxygens (FeO and and Fe2O3 will be added in loop)
        O2_sum = T.SiO2_Oxygen(s)...
                    + T.TiO2_Oxygen(s)...
                    + T.Al2O3_Oxygen(s)...
                    + T.Cr2O3_Oxygen(s)...
                    + T.MnO_Oxygen(s)...
                    + T.MgO_Oxygen(s)...
                    + T.NiO_Oxygen(s)...
                    + T.CaO_Oxygen(s)...
                    + T.Na2O_Oxygen(s)...
                    + T.K2O_Oxygen(s)...
                    + T.SO3_Oxygen(s)...
                    + T.P2O5_Oxygen(s)...
                    + T.V2O3_Oxygen(s)...
                    + T.ZnO_Oxygen(s)...
                    + T.CoO_Oxygen(s)...
                    + T.BaO_Oxygen(s)...
                    + T.SrO_Oxygen(s)...
                    + T.B2O3_Oxygen(s)...
                    + T.PbO_Oxygen(s)...
                    + T.CuO_Oxygen(s)...
                    + T.Sb2O3_Oxygen(s)...
                    + T.As2O5_Oxygen(s)...
                    + T.ThO2_Oxygen(s)...
                    + T.ZrO2_Oxygen(s)...
                    + T.HfO2_Oxygen(s)...
                    + T.Ag2O_Oxygen(s)...
                    + T.Y2O3_Oxygen(s)...
                    + T.La2O3_Oxygen(s)...
                    + T.Ce2O3_Oxygen(s)...
                    + T.Nd2O3_Oxygen(s)...
                    + FeO_Oxygen...             %changed to refer to the 'new' FeO (above)
                    + Fe2O3_Oxygen;             %added
      
                %Recalculate X_Normalized and X_Multi for ALL oxides (new
                %O2_sum)
                Si_Norm_new = T.SiO2_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Ti_Norm_new = T.TiO2_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Al_Norm_new = T.Al2O3_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Cr_Norm_new = T.Cr2O3_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Mn_Norm_new = T.MnO_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Mg_Norm_new = T.MgO_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Ni_Norm_new = T.NiO_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Ca_Norm_new = T.CaO_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Na_Norm_new = T.Na2O_Oxygen(s) * nor_Oxygen ./ O2_sum;
                K_Norm_new  = T.K2O_Oxygen(s) * nor_Oxygen ./ O2_sum;
                S_Norm_new  = T.SO3_Oxygen(s) * nor_Oxygen ./ O2_sum;
                P_Norm_new  = T.P2O5_Oxygen(s) * nor_Oxygen ./ O2_sum;
                V_Norm_new  = T.V2O3_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Zn_Norm_new  = T.ZnO_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Co_Norm_new  = T.CoO_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Ba_Norm_new  = T.BaO_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Sr_Norm_new  = T.SrO_Oxygen(s) * nor_Oxygen ./ O2_sum;
                B_Norm_new  = T.B2O3_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Pb_Norm_new  = T.PbO_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Cu_Norm_new  = T.CuO_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Sb_Norm_new  = T.Sb2O3_Oxygen(s) * nor_Oxygen ./ O2_sum;
                As_Norm_new  = T.As2O5_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Th_Norm_new  = T.ThO2_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Zr_Norm_new  = T.ZrO2_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Hf_Norm_new  = T.HfO2_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Ag_Norm_new  = T.Ag2O_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Y_Norm_new  = T.Y2O3_Oxygen(s) * nor_Oxygen ./ O2_sum;
                La_Norm_new  = T.La2O3_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Ce_Norm_new  = T.Ce2O3_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Nd_Norm_new  = T.Nd2O3_Oxygen(s) * nor_Oxygen ./ O2_sum;
                Cl_Norm_new  = T.Cl_Oxygen(s) * nor_Oxygen ./ O2_sum;
                F_Norm_new  = T.F_Oxygen(s) * nor_Oxygen ./ O2_sum;
                
                Si_Multi_new = Si_Norm_new * 1/2;
                Ti_Multi_new = Ti_Norm_new * 1/2;
                Al_Multi_new = Al_Norm_new * 2/3;
                Cr_Multi_new = Cr_Norm_new * 2/3;
                Mn_Multi_new = Mn_Norm_new * 1/1;
                Mg_Multi_new = Mg_Norm_new * 1/1;
                Ni_Multi_new = Ni_Norm_new * 1/1;
                Ca_Multi_new = Ca_Norm_new * 1/1;
                Na_Multi_new = Na_Norm_new * 2/1;
                K_Multi_new = K_Norm_new * 2/1;
                S_Multi_new = S_Norm_new * 1/3;
                P_Multi_new = P_Norm_new * 2/5;
                V_Multi_new = V_Norm_new * 2/3;
                Zn_Multi_new = Zn_Norm_new * 1;
                Co_Multi_new = Co_Norm_new * 1;
                Ba_Multi_new = Ba_Norm_new * 1;
                Sr_Multi_new = Sr_Norm_new * 1;
                B_Multi_new = B_Norm_new * 2/3;
                Pb_Multi_new = Pb_Norm_new * 1;
                Cu_Multi_new = Cu_Norm_new * 1;
                Sb_Multi_new = Sb_Norm_new * 2/3;
                As_Multi_new = As_Norm_new * 2/5;
                Th_Multi_new = Th_Norm_new * 1/2;
                Zr_Multi_new = Zr_Norm_new * 1/2;
                Hf_Multi_new = Hf_Norm_new * 1/2;
                Ag_Multi_new = Ag_Norm_new * 2/1;
                Y_Multi_new = Y_Norm_new * 2/3;
                La_Multi_new = La_Norm_new * 2/3;
                Ce_Multi_new = Ce_Norm_new * 2/3;
                Nd_Multi_new = Nd_Norm_new * 2/3;
                Cl_Multi_new = Cl_Norm_new * 1 / 1;
                F_Multi_new = F_Norm_new * 1 / 1;         
                
                %Calculate renormalized cations from all oxides
                Si_renorm = Si_Multi_new / Renorm_factor;
                Ti_renorm = Ti_Multi_new / Renorm_factor;
                Al_renorm = Al_Multi_new / Renorm_factor;
                Cr_renorm = Cr_Multi_new / Renorm_factor;
                Mn_renorm = Mn_Multi_new / Renorm_factor;
                Mg_renorm = Mg_Multi_new / Renorm_factor;
                Ni_renorm = Ni_Multi_new / Renorm_factor;
                Ca_renorm = Ca_Multi_new / Renorm_factor;
                Na_renorm = Na_Multi_new / Renorm_factor;
                K_renorm = K_Multi_new / Renorm_factor;
                S_renorm = S_Multi_new / Renorm_factor;
                P_renorm = P_Multi_new / Renorm_factor;
                V_renorm = V_Multi_new / Renorm_factor;
                Zn_renorm = Zn_Multi_new / Renorm_factor;
                Co_renorm = Co_Multi_new / Renorm_factor;
                Ba_renorm = Ba_Multi_new / Renorm_factor;
                Sr_renorm = Sr_Multi_new / Renorm_factor;
                B_renorm = B_Multi_new / Renorm_factor;
                Pb_renorm = Pb_Multi_new / Renorm_factor;
                Cu_renorm = Cu_Multi_new / Renorm_factor;
                Sb_renorm = Sb_Multi_new / Renorm_factor;
                As_renorm = As_Multi_new / Renorm_factor;
                Th_renorm = Th_Multi_new / Renorm_factor;
                Zr_renorm = Zr_Multi_new / Renorm_factor;
                Hf_renorm = Hf_Multi_new / Renorm_factor;
                Ag_renorm = Ag_Multi_new / Renorm_factor;
                Y_renorm = Y_Multi_new / Renorm_factor;
                La_renorm = La_Multi_new / Renorm_factor;
                Ce_renorm = Ce_Multi_new / Renorm_factor;
                Nd_renorm = Nd_Multi_new / Renorm_factor;
                Cl_renorm = Cl_Multi_new / Renorm_factor;
                F_renorm = F_Multi_new / Renorm_factor;
                
                
        %Recalculate Normalize, Multi, and cations for FeO and Fe2O3
        Fe_Norm_new = FeO_Oxygen * nor_Oxygen / O2_sum;     %Calculate normalized Fe from FeO
        Fe_Multi_new = Fe_Norm_new * 1/1;                   %Calculate Fe2 cations from FeO
        Fe2_renorm = Fe_Multi_new / Renorm_factor;          %Renormalize
        
        Fe2O3_Normalized = Fe2O3_Oxygen * nor_Oxygen / O2_sum;      %Calculate normalized Fe from Fe2O3
        Fe2O3_Multi = Fe2O3_Normalized * 2/3;                       %Calculate Fe3 cations from Fe2O3
        Fe3_renorm = Fe2O3_Multi / Renorm_factor;                   %Renormalize
        
        if(CD_cations == 1) %spinels (ABCD)
            CD_renorm_new = Fe2_renorm + Fe3_renorm + Mg_renorm + Mn_renorm + Ni_renorm + Zn_renorm + Al_renorm + Cr_renorm + Ti_renorm + V_renorm + Si_renorm; %Calcualte new sum of cations
        end
        
        if(CD_cations == 2) %pyroxenes (ABDCT)
            CD_renorm_new = Si_renorm + Ti_renorm + Al_renorm + Cr_renorm + Mn_renorm + Mg_renorm + Ni_renorm + Ca_renorm + Na_renorm + K_renorm + S_renorm + P_renorm + V_renorm + Zn_renorm + Co_renorm + Ba_renorm + Sr_renorm + B_renorm + Pb_renorm + Cu_renorm + Sb_renorm + As_renorm + Th_renorm + Zr_renorm + Hf_renorm + Ag_renorm + Y_renorm + La_renorm + Ce_renorm + Nd_renorm + Cl_renorm + F_renorm; %Calcualte new sum of cations
        end
        
        Fe2 = Fe2_renorm;
        T.Fe2_Calc(s) = Fe2;
        Fe3 = Fe3_renorm;
        T.Fe3_Calc(s) = Fe3;
        
        if(round(CD_renorm_new,3) <= Ideal_cations)
            %Once the Fe2 and Fe3 have 'equilibrated' to get CD <= 3 (to 3 decimals), end
            %the loop.
            break
        end  
    end   
end