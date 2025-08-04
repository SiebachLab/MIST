%Normalization of input data to 100%

function [T] = Data_Normalization(T_input)

T = T_input;



    for i = 1:height(T_input)
        
        T.Original_Total(i) = T_input.Input_SiO2(i)+T_input.Input_TiO2(i)+T_input.Input_Al2O3(i)+T_input.Input_Cr2O3(i)+...
            T_input.Input_FeO(i)+T_input.Input_NiO(i)+T_input.Input_MnO(i)+T_input.Input_MgO(i)+T_input.Input_CaO(i)+...
            T_input.Input_Na2O(i)+T_input.Input_K2O(i)+T_input.Input_P2O5(i)+T_input.Input_SO3(i)+T_input.Input_F(i)+T_input.Input_Cl(i)+...
            T_input.Input_V2O3(i)+T_input.Input_ZnO(i)+T_input.Input_CoO(i)+T_input.Input_BaO(i)+T_input.Input_SrO(i)+T_input.Input_B2O3(i)+T_input.Input_PbO(i)+...
            T_input.Input_CuO(i)+T_input.Input_Sb2O3(i)+T_input.Input_As2O5(i)+T_input.Input_ThO2(i)+T_input.Input_ZrO2(i)+T_input.Input_HfO2(i)+...
            T_input.Input_Ag2O(i)+T_input.Input_Y2O3(i)+T_input.Input_La2O3(i)+T_input.Input_Ce2O3(i)+T_input.Input_Nd2O3(i);

        T.SiO2(i) = T_input.Input_SiO2(i) * 100/T.Original_Total(i);	
        T.TiO2(i) = T_input.Input_TiO2(i) * 100/T.Original_Total(i);
        T.Al2O3(i) = T_input.Input_Al2O3(i) * 100/T.Original_Total(i);
        T.Cr2O3(i) = T_input.Input_Cr2O3(i) * 100/T.Original_Total(i);
        T.FeO(i) = T_input.Input_FeO(i) * 100/T.Original_Total(i);
        T.NiO(i) = T_input.Input_NiO(i) * 100/T.Original_Total(i);
        T.MnO(i) = T_input.Input_MnO(i) * 100/T.Original_Total(i);
        T.MgO(i) = T_input.Input_MgO(i) * 100/T.Original_Total(i);
        T.CaO(i) = T_input.Input_CaO(i) * 100/T.Original_Total(i);
        T.Na2O(i) = T_input.Input_Na2O(i) * 100/T.Original_Total(i);	
        T.K2O(i) = T_input.Input_K2O(i) * 100/T.Original_Total(i);
        T.P2O5(i) = T_input.Input_P2O5(i) * 100/T.Original_Total(i);
        T.SO3(i) = T_input.Input_SO3(i) * 100/T.Original_Total(i);
        T.F(i) = T_input.Input_F(i) * 100/T.Original_Total(i);										
        T.Cl(i) = T_input.Input_Cl(i) * 100/T.Original_Total(i);	
        T.V2O3(i) = T_input.Input_V2O3(i) * 100/T.Original_Total(i);	
        T.ZnO(i) = T_input.Input_ZnO(i) * 100/T.Original_Total(i);
        T.CoO(i) = T_input.Input_CoO(i) * 100/T.Original_Total(i);
        T.BaO(i) = T_input.Input_BaO(i) * 100/T.Original_Total(i);
        T.SrO(i) = T_input.Input_SrO(i) * 100/T.Original_Total(i);
        T.B2O3(i) = T_input.Input_B2O3(i) * 100/T.Original_Total(i);
        T.PbO(i) = T_input.Input_PbO(i) * 100/T.Original_Total(i);	
        T.CuO(i) = T_input.Input_CuO(i) * 100/T.Original_Total(i);	
        T.Sb2O3(i) = T_input.Input_Sb2O3(i) * 100/T.Original_Total(i);
        T.As2O5(i) = T_input.Input_As2O5(i) * 100/T.Original_Total(i);
        T.ThO2(i) = T_input.Input_ThO2(i) * 100/T.Original_Total(i);	
        T.ZrO2(i) = T_input.Input_ZrO2(i) * 100/T.Original_Total(i);	
        T.HfO2(i) = T_input.Input_HfO2(i) * 100/T.Original_Total(i);	
        T.Ag2O(i) = T_input.Input_Ag2O(i) * 100/T.Original_Total(i);	
        T.Y2O3(i) = T_input.Input_Y2O3(i) * 100/T.Original_Total(i);	
        T.La2O3(i) = T_input.Input_La2O3(i) * 100/T.Original_Total(i);	
        T.Ce2O3(i) = T_input.Input_Ce2O3(i) * 100/T.Original_Total(i);
        T.Nd2O3(i) = T_input.Input_Nd2O3(i) * 100/T.Original_Total(i);

        T.New_Total(i) = round(T.SiO2(i)+T.TiO2(i)+T.Al2O3(i)+T.Cr2O3(i)+T.FeO(i)+...
            T.NiO(i)+T.MnO(i)+T.MgO(i)+T.CaO(i)+T.Na2O(i)+T.K2O(i)+T.P2O5(i)+T.SO3(i)+T.F(i)+...
            T.Cl(i)+T.V2O3(i)+T.ZnO(i)+T.CoO(i)+T.BaO(i)+T.SrO(i)+T.B2O3(i)+T.PbO(i)+T.CuO(i)+...
            T.Sb2O3(i)+T.As2O5(i)+T.ThO2(i)+T.ZrO2(i)+T.HfO2(i)+T.Ag2O(i)+T.Y2O3(i)+T.La2O3(i)+T.Ce2O3(i)+T.Nd2O3(i));
        
    end 
end
     