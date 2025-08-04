function [T_output] = Check_Function(T)
    for i = 1:height(T)
        
        %SiO2 Phase check
        if(strcmp(T.group2(i), 'SiO2 Phase'))
            T.ABCDT_ideal(i) = 1;
            T.ABCD_ideal(i) = 0;
            T.T_ideal(i) = 1;
            T.T_prime(i) = (Si_n2 + Al_n2 + Ti_n2 + Cr_n2 + Fe_n2 + Mg_n2);
            if(strcmp(T.species(i), 'Obsidian glass'))
                T.ABCDT_ideal(i) = 1.03;
                T.T_ideal(i) = 1.03;
            end    
            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/ABCDT_n2)*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/ABCD_n2)*100;
            T.Delta_T(i) = (T.T_ideal(i)/T.T_prime(i))*100;

            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

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
        
        %Carbonates check
        if(strcmp(T.group2(i), 'Monoclinic Carbonate')||strcmp(T.group2(i), 'Trigonal Carbonate')|| strcmp(T.group2(i), 'Orthorhombic Carbonate'))
            if(strcmp(T.group2(i), 'Trigonal Carbonate') || strcmp(T.group2(i), 'Orthorhombic Carbonate'))
                T.ABCDT_ideal(i) = 2;   T.ABCD_ideal(i) = 2;
            elseif(strcmp(T.group2(i), 'Monoclinic Carbonate')) 
                T.ABCDT_ideal(i) = 4;   T.ABCD_ideal(i) = 4;
            end     
        
            T.T_ideal(i) = 0;
            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(ABCDT_n2+S_n2))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(ABCD_n2+S_n2))*100;
            T.Delta_T(i) = (T.T_ideal(i)/T_n2)*100;

            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

            if(strcmp(T.species(i), 'Pure magnetite/hematite/wustite/Fe-Hydroxide'))
                    T.check(i) = {'no check'};
            elseif(strcmp(T.group2(i), 'Pure Carbonate or Oxide-Hydroxide') || strcmp(T.group3(i), 'Pure Carbonate or Oxide-Hydroxide') || strcmp(T.group4(i), 'Pure Carbonate or Oxide-Hydroxide') || strcmp(T.species(i), 'Pure Carbonate or Oxide-Hydroxide'))
                    T.check(i) = {'no check'};
            elseif(strcmp(T.species(i), 'Manganosite/Pyrolusite/Hausmannite/Birnessite/Bixbyite/etc.'))
                T.check(i) = {'no check'};
            elseif(ABCDT_range == 0 && ABCD_range == 0)
                T.check(i) = {'TRUE'};
            elseif(ABCDT_range == 1 && ABCD_range == 1)
                T.check(i) = {'FALSE'};
            elseif(ABCDT_range == 0 && ABCD_range == 1)
                T.check(i) = {'MAYBE SIMILAR TO'};
            elseif(ABCDT_range == 1 && ABCD_range == 0)
                T.check(i) = {'MAYBE RELATED'};     
            end
        else
            T.check(i) = {'no check'};    
        end

        %Sulfates check
        if(strcmp(T.group1(i), 'Sulfate'))    
            if(strcmp(T.group2(i), ''))
                T.ABCDT_ideal(i) = nan;   T.ABCD_ideal(i) = nan;
            elseif(strcmp(T.group4(i), 'Ca-Sulfate'))
                T.ABCDT_ideal(i) = 2;   T.ABCD_ideal(i) = 2;
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
            elseif(strcmp(T.species(i), 'blodite'))
                T.ABCDT_ideal(i) = 5*(4/8);   T.ABCD_ideal(i) = 5*(4/8);
            elseif(strcmp(T.group4(i), 'Fe-Sulfate'))
                if(strcmp(T.species(i), 'rhomboclase'))
                    T.ABCDT_ideal(i) = 3.5*(4/8);   T.ABCD_ideal(i) = 3.5*(4/8);
                elseif(strcmp(T.species(i), 'mikasaite/coquimbite'))
                    T.ABCDT_ideal(i) = 5.5*(4/12);   T.ABCD_ideal(i) = 5.5*(4/12);
                elseif(strcmp(T.species(i), 'Fe-Sulfate (Fe/S = 1)') || strcmp(T.species(i), 'siderotile'))
                    T.ABCDT_ideal(i) = 2;   T.ABCD_ideal(i) = 2;
                elseif(strcmp(T.species(i), 'romerite/bilinite'))
                    T.ABCDT_ideal(i) = 7.5*(4/16);   T.ABCD_ideal(i) = 7.5*(4/16);
                elseif(strcmp(T.species(i), 'copiapite/ferricopiapite'))
                    T.ABCDT_ideal(i) = 11.67*(4/25);   T.ABCD_ideal(i) = 11.67*(4/25);
                elseif(strcmp(T.species(i), 'hydroniumjarosite'))
                    T.ABCDT_ideal(i) = 6*(4/11);   T.ABCD_ideal(i) = 6*(4/11);
                elseif(strcmp(T.species(i), 'schwertmannite'))
                    T.ABCDT_ideal(i) = 24.5*(4/33.6);   T.ABCD_ideal(i) = 24.5*(4/33.6);   
                elseif(strcmp(T.species(i), 'viaeneite'))
                    T.ABCDT_ideal(i) = 12*(4/27);   T.ABCD_ideal(i) = 12*(4/27);  
                end    
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
            elseif(strcmp(T.group4(i), 'Na-K-Fe-Sulfate'))
                if(strcmp(T.species(i), 'eldfellite/amarillite/erdite'))
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
                end  
            elseif(strcmp(T.species(i), 'kainite'))
                T.ABCDT_ideal(i) = 3*(4/4.5);   T.ABCD_ideal(i) = 3*(4/4.5);
            elseif(strcmp(T.species(i), 'Alunite'))
                T.ABCDT_ideal(i) = 6*(4/11);   T.ABCD_ideal(i) = 6*(4/11); 
            elseif(strcmp(T.species(i), 'caledonite'))
                T.ABCDT_ideal(i) = 10*(4/16);   T.ABCD_ideal(i) = 10*(4/16);    
            end

            T.T_ideal(i) = 0;
            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDTS(i)/6))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCDS(i)/6))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T(i)/6))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

            if(ABCDT_range == 0 && ABCD_range == 0)
                T.check(i) = {'TRUE'};
            elseif(ABCDT_range == 1 && ABCD_range == 1)
                T.check(i) = {'FALSE'};
            elseif(ABCDT_range == 0 && ABCD_range == 1)
                T.check(i) = {'MAYBE SIMILAR TO'};
            elseif(ABCDT_range == 1 && ABCD_range == 0)
                T.check(i) = {'MAYBE RELATED'};
            end
        else
            T.check(i) = {'no check'};    
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
               if(strcmp(T.species(i), 'dry: chopinite or ferringtonite, or OH-bearing: baricite, bobierrite, cattiite'))
                   T.ABCDT_ideal(i) = 5*(4/8);   T.ABCD_ideal(i) = 5*(4/8);
               elseif(strcmp(T.species(i), 'holtedahlite or hydroxylwagnerite or kovdorskite'))
                   T.ABCDT_ideal(i) = 3*(4/4.5);   T.ABCD_ideal(i) = 3*(4/4.5);
               elseif(strcmp(T.species(i), 'newberyite or phospharosslerite'))
                   T.ABCDT_ideal(i) = 2*(4/3.5);   T.ABCD_ideal(i) = 2*(4/3.5);   
               elseif(strcmp(T.species(i), 'raaedite'))
                   T.ABCDT_ideal(i) = 9*(4/12);   T.ABCD_ideal(i) = 9*(4/12);    
               end 
            elseif(strcmp(T.group2(i), 'Fe-Phosphate'))
                if(strcmp(T.species(i), 'graftonite or sarcopside'))
                   T.ABCDT_ideal(i) = 5*(4/8);   T.ABCD_ideal(i) = 5*(4/8);
                elseif(strcmp(T.species(i), 'grattarolaite'))
                   T.ABCDT_ideal(i) = 5*(4/7);   T.ABCD_ideal(i) = 5*(4/7);
                elseif(strcmp(T.species(i), 'heteroisite or rodolicoite'))
                   T.ABCDT_ideal(i) = 2.25*(4/4);   T.ABCD_ideal(i) = 2.25*(4/4);  
                end 
            elseif(strcmp(T.species(i), 'archerite'))
               T.ABCDT_ideal(i) = 1.85*(4/3);   T.ABCD_ideal(i) = 1.85*(4/3);
            elseif(strcmp(T.group2(i), 'REE-Phosphate') || strcmp(T.species(i), 'Xenotime-Y'))  
               T.ABCDT_ideal(i) = 2*(4/4);   T.ABCD_ideal(i) = 2*(4/4); 
            end

            T.T_ideal(i) = 0;
            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDTP(i)/6))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCDP(i)/6))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T(i)/6))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

            if(ABCDT_range == 0 && ABCD_range == 0)
                T.check(i) = {'TRUE'};
            elseif(ABCDT_range == 1 && ABCD_range == 1)
                T.check(i) = {'FALSE'};
            elseif(ABCDT_range == 0 && ABCD_range == 1)
                T.check(i) = {'MAYBE SIMILAR TO'};
            elseif(ABCDT_range == 1 && ABCD_range == 0)
                T.check(i) = {'MAYBE RELATED'};
            end
        else
            T.check(i) = {'no check'};    
        end
        
        
        %%Spinels Check
        if(strcmp(T.group4(i), 'Spinel'))  
           T.ABCDT_ideal(i) = 3.02*(4/4);   T.ABCD_ideal(i) = 3.02*(4/4); 
           if(strcmp(T.species(i), 'Ti-Magnetite'))
               T.ABCDT_ideal(i) = 3.5*(4/4);   T.ABCD_ideal(i) = 3.5*(4/4);
           end  
           if(strcmp(T.species(i), 'magnetite'))
               T.ABCDT_ideal(i) = 4*(4/4);   T.ABCD_ideal(i) = 4*(4/4);
           end  
           if(strcmp(T.species(i), 'Eskolaite'))
               T.ABCDT_ideal(i) = 2*(4/3);   T.ABCD_ideal(i) = 2*(4/3);
           end 

            T.T_ideal(i) = 0;
            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)/6))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)/6))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T(i)/6))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

            if(strcmp(T.species(i), 'Pure magnetite/hematite/wustite/Fe-Hydroxide'))
                T.check(i) = {'no check'};
            elseif(strcmp(T.group2(i), 'Pure Carbonate or Oxide-Hydroxide') || strcmp(T.group3(i), 'Pure Carbonate or Oxide-Hydroxide') || strcmp(T.group4(i), 'Pure Carbonate or Oxide-Hydroxide') || strcmp(T.species(i), 'Pure Carbonate or Oxide-Hydroxide'))
                T.check(i) = {'no check'};
            elseif(strcmp(T.species(i), 'Manganosite/Pyrolusite/Hausmannite/Birnessite/Bixbyite/etc.'))
                T.check(i) = {'no check'};
            elseif(ABCDT_range == 0 && ABCD_range == 0)
                T.check(i) = {'TRUE'};
            elseif(ABCDT_range == 1 && ABCD_range == 1)
                T.check(i) = {'FALSE'};
            elseif(ABCDT_range == 0 && ABCD_range == 1)
                T.check(i) = {'MAYBE SIMILAR TO'};
            elseif(ABCDT_range == 1 && ABCD_range == 0)
                T.check(i) = {'MAYBE RELATED'};
            end
        end
        
        %%Illmenties check
        if(strcmp(T.group4(i), 'Ilmenite'))  
           T.ABCDT_ideal(i) = 2*(3/3);   T.ABCD_ideal(i) = 2*(3/3);    

            T.T_ideal(i) = 0;
            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)/8))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)/8))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T(i)/6))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

            if(strcmp(T.species(i), 'Pure magnetite/hematite/wustite/Fe-Hydroxide'))
                T.check(i) = {'no check'};
            elseif(strcmp(T.group2(i), 'Pure Carbonate or Oxide-Hydroxide') || strcmp(T.group3(i), 'Pure Carbonate or Oxide-Hydroxide') || strcmp(T.group4(i), 'Pure Carbonate or Oxide-Hydroxide') || strcmp(T.species(i), 'Pure Carbonate or Oxide-Hydroxide'))
                T.check(i) = {'no check'};
            elseif(strcmp(T.species(i), 'Manganosite/Pyrolusite/Hausmannite/Birnessite/Bixbyite/etc.'))
                T.check(i) = {'no check'};                            
            elseif(ABCDT_range == 0 && ABCD_range == 0)
                T.check(i) = {'TRUE'};
            elseif(ABCDT_range == 1 && ABCD_range == 1)
                T.check(i) = {'FALSE'};
            elseif(ABCDT_range == 0 && ABCD_range == 1)
                T.check(i) = {'MAYBE SIMILAR TO'};
            elseif(ABCDT_range == 1 && ABCD_range == 0)
                T.check(i) = {'MAYBE RELATED'};
            end
        end

        %Oxide-Hydroxide check (minus spinels and ilmenties)
        if( (strcmp(T.group2(i), 'Oxide') || strcmp(T.group2(i), 'Hydroxide') || strcmp(T.group2(i), 'Oxide or Hydroxide')) && ~strcmp(T.group4(i), 'Spinel') && ~strcmp(T.group4(i), 'Ilmenite'))
            if(strcmp(T.species(i), 'Manganosite or manganite') || strcmp(T.species(i), 'Pyrochroite') || strcmp(T.group4(i), 'Mg-Oxide or Mg-Hydroxide') || strcmp(T.species(i), 'brucite') || strcmp(T.species(i), 'amakinite or wustite'))  
                T.ABCDT_ideal(i) = 1*(24/1);   T.ABCD_ideal(i) = 1*(24/1);     
            elseif(strcmp(T.species(i), 'Bixbyite'))
                T.ABCDT_ideal(i) = 2.5*(24/3);   T.ABCD_ideal(i) = 2.5*(24/3);
            elseif(strcmp(T.species(i), 'Pyrolusite'))
                T.ABCDT_ideal(i) = 2*(24/2);   T.ABCD_ideal(i) = 2*(24/2);
            elseif(strcmp(T.species(i), 'Hausmannite'))
                T.ABCDT_ideal(i) = 2.5*(24/4);   T.ABCD_ideal(i) = 2.5*(24/4);   
            elseif(strcmp(T.species(i), 'Birnessite'))
                T.ABCDT_ideal(i) = 3.95*(24/4);   T.ABCD_ideal(i) = 3.95*(24/4);
            elseif(strcmp(T.group4(i), 'Ti-Oxide'))
                T.ABCDT_ideal(i) = 1*(24/2);   T.ABCD_ideal(i) = 1*(24/2);
            elseif(strcmp(T.species(i), 'corundum'))
                T.ABCDT_ideal(i) = 2*(24/3);   T.ABCD_ideal(i) = 2*(24/3);        
            elseif(strcmp(T.species(i), 'Aurorite'))
                T.ABCDT_ideal(i) = 6.5*(24/7);   T.ABCD_ideal(i) = 6.5*(24/7);  
            elseif(strcmp(T.species(i), 'Vernadite'))
                T.ABCDT_ideal(i) = 2.1*(24/2);   T.ABCD_ideal(i) = 2.1*(24/2); 
            elseif(strcmp(T.species(i), 'Janggunite'))
                T.ABCDT_ideal(i) = 8*(24/11);   T.ABCD_ideal(i) = 8*(24/11);  
            elseif(strcmp(T.species(i), 'Takanelite'))
                T.ABCDT_ideal(i) = 8.75*(24/9);   T.ABCD_ideal(i) = 8.75*(24/9);   
            elseif(strcmp(T.species(i), 'Hollandite or Romanechite'))
                T.ABCDT_ideal(i) = 13.75*(24/16);   T.ABCD_ideal(i) = 13.75*(24/16);   
            elseif(strcmp(T.species(i), 'akaganeite'))
                T.ABCDT_ideal(i) = 8.75*(24/8.75);   T.ABCD_ideal(i) = 8.75*(24/8.75);
            end 

            T.T_ideal(i) = 0;
            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

            if(strcmp(T.species(i), 'Pure magnetite/hematite/wustite/Fe-Hydroxide'))
                T.check(i) = {'no check'};
            elseif(strcmp(T.group2(i), 'Pure Carbonate or Oxide-Hydroxide') || strcmp(T.group3(i), 'Pure Carbonate or Oxide-Hydroxide') || strcmp(T.group4(i), 'Pure Carbonate or Oxide-Hydroxide') || strcmp(T.species(i), 'Pure Carbonate or Oxide-Hydroxide'))
                T.check(i) = {'no check'};
            elseif(strcmp(T.species(i), 'Manganosite/Pyrolusite/Hausmannite/Birnessite/Bixbyite/etc.'))
                T.check(i) = {'no check'};
            elseif(strcmp(T.species(i), 'Al-Oxide or Al-Hydroxide (corundum, gibbsite, diaspore, or others)'))
                T.check(i) = {'no check'};       
            elseif(ABCDT_range == 0 && ABCD_range == 0)
                T.check(i) = {'TRUE'};
            elseif(ABCDT_range == 1 && ABCD_range == 1)
                T.check(i) = {'FALSE'};
            elseif(ABCDT_range == 0 && ABCD_range == 1)
                T.check(i) = {'MAYBE SIMILAR TO'};
            elseif(ABCDT_range == 1 && ABCD_range == 0)
                T.check(i) = {'MAYBE RELATED'}; 
            end      
        else
            T.check(i) = {'no check'};
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
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

            if(ABCDT_range == 0 && ABCD_range == 0)
                T.check(i) = {'TRUE'};
            elseif(ABCDT_range == 1 && ABCD_range == 1)
                T.check(i) = {'FALSE'};
            elseif(ABCDT_range == 0 && ABCD_range == 1)
                T.check(i) = {'MAYBE SIMILAR TO'};
            elseif(ABCDT_range == 1 && ABCD_range == 0)
                T.check(i) = {'MAYBE RELATED'};
            end

            if(strcmp(T.species(i), ''))
                T.check(i) = {'no check'};
            end
        else
            T.check(i) = {'no check'};
        end
        
        %%Feldspars check
        if(strcmp(T.group3(i), 'Feldspar'))
            T.T_prime(i) = Si_n8 + Al_n8 + Fe_n8 + Mg_n8;
            T.ABCDT_ideal(i) = 5*(24/8);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 1;

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

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
        
        %%Feldspathoids check
        if(strcmp(T.group3(i), 'Feldspathoid'))
            if(strcmp(T.species(i), 'sodalite'))
                T.T_prime(i) = Si_n12 + Al_n12 + Fe_n12 + Mg_n12;
                T.ABCDT_ideal(i) = 10*(24/12.5);   T.T_ideal(i) = 6;  T.ABCD_ideal(i) = 4;
            end  

            if(strcmp(T.species(i), 'leucite'))
                T.T_prime(i) = Si_n6 + Al_n6 + Fe_n6 + Mg_n6;
                T.ABCDT_ideal(i) = 4*(24/6);   T.T_ideal(i) = 3;  T.ABCD_ideal(i) = 1;
            end

            if(strcmp(T.species(i), 'cancrinite'))
                T.T_prime(i) = Si_n + Al_n + Fe_n + Mg_n;
                T.ABCDT_ideal(i) = 18.5*(24/24);   T.T_ideal(i) = 11.5;  T.ABCD_ideal(i) = 8;
            end

            if(strcmp(T.species(i), 'nepheline'))
                T.T_prime(i) = Si_n16 + Al_n16 + Fe_n16 + Mg_n16;
                T.ABCDT_ideal(i) = 12*(24/16);   T.T_ideal(i) = 8;  T.ABCD_ideal(i) = 4;
            end

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

            if(ABCDT_range == 0 && T_range == 0)
                T.check(i) = {'TRUE'};
            elseif(ABCDT_range == 1 && T_range == 1)
                T.check(i) = {'FALSE'};
            elseif(ABCDT_range == 0 && T_range == 1)
                T.check(i) = {'MAYBE SIMILAR TO'};
            elseif(ABCDT_range == 1 && T_range == 0)
                T.check(i) = {'MAYBE RELATED'};
            end
        else
            T.check(i) = {'no check'};         
        end
        
        %%Zeolites check
        if(strcmp(T.group3(i), 'Zeolite'))
            if(strcmp(T.species(i), 'Clinoptilolite-Ca'))
                T.T_prime(i) = Si_n72 + Al_n72 + Fe_n72 + Mg_n72;
                T.ABCDT_ideal(i) = 39*(24/72);   T.T_ideal(i) = 36;  T.ABCD_ideal(i) = 3;
            end  

            if(strcmp(T.species(i), 'Clinoptilolite-Na') || strcmp(T.species(i), 'Clinoptilolite-K') || strcmp(T.species(i), 'Heulandite-Na'))
                T.T_prime(i) = Si_n72 + Al_n72 + Fe_n72 + Mg_n72;
                T.ABCDT_ideal(i) = 42*(24/72);   T.T_ideal(i) = 36;  T.ABCD_ideal(i) = 6;
            end

            if(strcmp(T.species(i), 'Heulandite-Ca') || strcmp(T.species(i), 'Heulandite-K') || strcmp(T.species(i), 'Stilbite-Ca') )
                T.T_prime(i) = Si_n72 + Al_n72 + Fe_n72 + Mg_n72;
                T.ABCDT_ideal(i) = 41*(24/72);   T.T_ideal(i) = 36;  T.ABCD_ideal(i) = 5;
            end

            if(strcmp(T.species(i), 'Stilbite-Na'))
                T.T_prime(i) = Si_n72 + Al_n72 + Fe_n72 + Mg_n72;
                T.ABCDT_ideal(i) = 45*(24/72);   T.T_ideal(i) = 36;  T.ABCD_ideal(i) = 9;
            end

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

            if(ABCDT_range == 0 && T_range == 0)
                T.check(i) = {'TRUE'};
            elseif(ABCDT_range == 1 && T_range == 1)
                T.check(i) = {'FALSE'};
            elseif(ABCDT_range == 0 && T_range == 1)
                T.check(i) = {'MAYBE SIMILAR TO'};
            elseif(ABCDT_range == 1 && T_range == 0)
                T.check(i) = {'MAYBE RELATED'};
            end
        else
            T.check(i) = {'no check'};             
        end
        
        %%Olivines check
        if(strcmp(T.group3(i), 'Olivine'))

            T.T_prime(i) = Si_n4 + Al_n4 + Ti_n4 + Cr_n4;
            T.ABCDT_ideal(i) = 3*(24/4);   T.T_ideal(i) = 1;  T.ABCD_ideal(i) = 2;                        

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

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
        
        %%Garnets check
        if(strcmp(T.group3(i), 'Garnet')) 
            if(strcmp(T.group4(i), 'Pyralspite Garnet') || strcmp(T.species(i), 'grossular') || strcmp(T.species(i), 'uvarovite'))
                Al_IV = 3 - Si_n12; Al_VI = Al_n12 - Al_IV;
                T.T_prime(i) = Si_n12 + Al_IV;
                T.ABCDT_ideal(i) = 8*(24/12);   T.T_ideal(i) = 3;  T.ABCD_ideal(i) = 5;                        
            end
            if(strcmp(T.species(i), 'andradite'))
                Al_IV = 3 - Si_n12; Al_VI = Al_n12 - Al_IV;
                T.T_prime(i) = Si_n12 + Al_IV;
                T.ABCDT_ideal(i) = 8.5*(24/12);   T.T_ideal(i) = 3;  T.ABCD_ideal(i) = 5;                        
            end

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

            if(ABCDT_range == 0 && T_range == 0)
                T.check(i) = {'TRUE'};
            elseif(ABCDT_range == 1 && T_range == 1)
                T.check(i) = {'FALSE'};
            elseif(ABCDT_range == 0 && T_range == 1)
                T.check(i) = {'MAYBE SIMILAR TO'};
            elseif(ABCDT_range == 1 && T_range == 0)
                T.check(i) = {'MAYBE RELATED'};
            end
        else
            T.check(i) = {'no check'};        
        end
        
        %%Al-silicates check
        if(strcmp(T.group3(i), 'Al-Silicate')) 
            T.T_prime(i) = Si_n5 + Al_IV;
            T.ABCDT_ideal(i) = 3*(24/5);   T.T_ideal(i) = 1;  T.ABCD_ideal(i) = 2;                        

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

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
        
        %%Zircon check
        if(strcmp(T.species(i), 'zircon')) 
            Al_IV = 1 - Si_n4; Al_VI = Al_n4 - Al_IV;
            T.T_prime(i) = Si_n4 + Al_IV;
            T.ABCDT_ideal(i) = 2*(24/4);   T.T_ideal(i) = 1;  T.ABCD_ideal(i) = 2;                        

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

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
        
        %%Pyroxenes check
        if(strcmp(T.group3(i), 'Pyroxene')) 
            Al_IV = 2 - Si_n6; Al_VI = Al_n6 - Al_IV;
            T.T_prime(i) = Si_n6 + Al_IV;
            T.ABCDT_ideal(i) = 4*(24/6);   T.T_ideal(i) = 2;  T.ABCD_ideal(i) = 2;

            if(strcmp(T.species(i), 'Aegirine')) 
                Al_IV = 2 - Si_n6; Al_VI = Al_n6 - Al_IV;
                T.T_prime(i) = Si_n6 + Al_IV;
                T.ABCDT_ideal(i) = 4.25*(24/6);   T.T_ideal(i) = 2;  T.ABCD_ideal(i) = 4.25-2;
            end

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

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
        
        %%Amphiboles check
        if(strcmp(T.group3(i), 'Amphibole')) 
            Al_IV = 8 - Si_n23; Al_VI = Al_n23 - Al_IV;
            T.T_prime(i) = Si_n23 + Al_IV;
            T.ABCDT_ideal(i) = 15.5*(24/23);   T.T_ideal(i) = 8;  T.ABCD_ideal(i) = 15.5-8;

            if(strcmp(T.species(i), 'magnesio-hastingsite') || strcmp(T.species(i), 'arfvedsonite'))
                T.ABCDT_ideal(i) = 15.75*(24/23);   T.T_ideal(i) = 8;  T.ABCD_ideal(i) = 15.5-8;
            end

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

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
               
        %%Micas check
        if(strcmp(T.group3(i), 'Mica')) 
            Al_IV = 4 - Si_n11; Al_VI = Al_n11 - Al_IV;
            T.T_prime(i) = Si_n11 + Al_IV;

            if(strcmp(T.species(i), 'biotite (annite)') || strcmp(T.species(i), 'biotite (phlogopite)') || strcmp(T.species(i), 'siderophyllite'))
                T.ABCDT_ideal(i) = 8*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 4;
            end    
            if(strcmp(T.species(i), 'ferroaluminoceladonite') || strcmp(T.species(i), 'muscovite') || strcmp(T.species(i), 'paragonite') || strcmp(T.species(i), 'margarite'))
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
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

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
        
        %%Talc-Pyrophyllite check
        if(strcmp(T.group3(i), 'Talc-Pyrophyllite')) 
            Al_IV = 4 - Si_n11; Al_VI = Al_n11 - Al_IV;
            T.T_prime(i) = Si_n11 + Al_IV;

            if(strcmp(T.species(i), 'talc') || strcmp(T.species(i), 'minnesotaite'))
                T.ABCDT_ideal(i) = 7*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 4;   
            end
            if(strcmp(T.species(i), 'Pyrophyllite'))
                T.ABCDT_ideal(i) = 6*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 4;   
            end

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

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
        
        %%Chlorites check
        if(strcmp(T.group3(i), 'Chlorite')) 
            Al_IV = 4 - Si_n14; Al_VI = Al_n14 - Al_IV;
            T.T_prime(i) = Si_n14 + Al_IV;
            T.ABCDT_ideal(i) = 10*(24/14);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 6;

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

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
        
        %%Kaolinite-Serpentine check
        if(strcmp(T.group3(i), 'Kaolinite-Serpentine'))
            Al_IV = 2 - Si_n7; Al_VI = Al_n7 - Al_IV;
            T.T_prime(i) = Si_n7 + Al_IV;

            if(strcmp(T.group4(i), 'Serpentine'))
                T.ABCDT_ideal(i) = 5*(24/7);   T.T_ideal(i) = 2;  T.ABCD_ideal(i) = 3;
            end
            if(strcmp(T.group4(i), 'Kaolinite'))
                T.ABCDT_ideal(i) = 4*(24/7);   T.T_ideal(i) = 2;  T.ABCD_ideal(i) = 2;
            end

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

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
        
        %%Smectites check
        if(strcmp(T.group3(i), 'Smectite'))
            Al_IV = 4 - Si_n11; Al_VI = Al_n11 - Al_IV;
            T.T_prime(i) = Si_n11 + Al_IV;

            if(strcmp(T.species(i), 'montmorillonite') || strcmp(T.species(i), 'Beidellite'))
                T.ABCDT_ideal(i) = 6.3*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 3;
            end
            if(strcmp(T.species(i), 'hectorite') || strcmp(T.species(i), 'Hectorite') || strcmp(T.species(i), 'saponite') || strcmp(T.species(i), 'sauconite') || strcmp(T.species(i), 'Nontronite'))
                T.ABCDT_ideal(i) = 7.3*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 2;
            end
            if(strcmp(T.species(i), 'vermiculite'))
                T.ABCDT_ideal(i) = 7.35*(24/11);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 2;
            end

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

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
        
        %%Palygorskite-Sepiolite check
        if(strcmp(T.group3(i), 'Palygorskite-Sepiolite'))

            if(strcmp(T.species(i), 'palygorskite'))
                Al_IV = 4 - Si_n10; Al_VI = Al_n10 - Al_IV;
                T.T_prime(i) = Si_n10 + Al_IV;
                T.ABCDT_ideal(i) = 6*(24/10.5);   T.T_ideal(i) = 4;  T.ABCD_ideal(i) = 2;
            end
            if(strcmp(T.species(i), 'sepiolite'))
                Al_IV = 6 - Si_n16; Al_VI = Al_n16 - Al_IV;
                T.T_prime(i) = Si_n16 + Al_IV;
                T.ABCDT_ideal(i) = 10*(24/16);   T.T_ideal(i) = 6;  T.ABCD_ideal(i) = 6;
            end

            T.Delta_ABCDT(i) = (T.ABCDT_ideal(i)/(T.ABCDT(i)))*100;
            T.Delta_ABCD(i) = (T.ABCD_ideal(i)/(T.ABCD(i)))*100;
            T.Delta_T(i) = (T.T_ideal(i)/(T.T_prime(i)))*100;    
            ABCDT_range = T.Delta_ABCDT(i) <= lowerdelta || T.Delta_ABCDT(i) >= upperdelta;
            ABCD_range = T.Delta_ABCD(i) <= lowerdelta || T.Delta_ABCD(i) >= upperdelta;
            T_range = T.Delta_T(i) <= lowerdelta || T.Delta_T(i) >= upperdelta;

            if(ABCDT_range == 0 && T_range == 0)
                T.check(i) = {'TRUE'};
            elseif(ABCDT_range == 1 && T_range == 1)
                T.check(i) = {'FALSE'};
            elseif(ABCDT_range == 0 && T_range == 1)
                T.check(i) = {'MAYBE SIMILAR TO'};
            elseif(ABCDT_range == 1 && T_range == 0)
                T.check(i) = {'MAYBE RELATED'};
            end
        else    
            T.check(i) = {'no check'};
        end
        
        
        
    end
    
    T_output = T;
    
end