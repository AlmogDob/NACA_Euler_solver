clc; clear;

db_file = 'NACA.db';
db_table = sqlread(sqlite(db_file),'NACA_data');

%%
clc;
NACA_list_fat = fetch(sqlite(db_file),sprintf('SELECT NACA, Mach_inf, angle_of_attack_deg FROM NACA_data;'));
NACA_list = [];

NACA_list_fat = sortrows(NACA_list_fat);

for i = 1:length(NACA_list_fat.NACA)
    current_NACA = NACA_list_fat.NACA(i);
    current_Mach_inf = NACA_list_fat.Mach_inf(i);
    current_alpha_deg = NACA_list_fat.angle_of_attack_deg(i);

    NACA_index = get_index_of_NACA_from_NACA_list(NACA_list, current_NACA);
    if -1 == NACA_index
        NACA_list(end+1).NACA = NACA_list_fat.NACA(i);
        NACA_list(end).max_Mach_inf = current_Mach_inf;
        NACA_list(end).min_Mach_inf = current_Mach_inf;
        NACA_list(end).max_alpha_deg = current_alpha_deg;
        NACA_list(end).min_alpha_deg = current_alpha_deg;
    else
        if NACA_list(NACA_index).max_Mach_inf < current_Mach_inf
            NACA_list(NACA_index).max_Mach_inf = current_Mach_inf;
        end
        if NACA_list(NACA_index).min_Mach_inf > current_Mach_inf
            NACA_list(NACA_index).min_Mach_inf = current_Mach_inf;
        end
        if NACA_list(NACA_index).max_alpha_deg < current_alpha_deg
            NACA_list(NACA_index).max_alpha_deg = current_alpha_deg;
        end
        if NACA_list(NACA_index).min_alpha_deg > current_alpha_deg
            NACA_list(NACA_index).min_alpha_deg = current_alpha_deg;
        end
    end
    
end

for i = 1:length(NACA_list)
    fprintf('%6s: Mach ∈ [%3g,%3g] | alpha ∈ [%3g,%3g]\n', NACA_list(i).NACA, NACA_list(i).min_Mach_inf, NACA_list(i).max_Mach_inf, NACA_list(i).min_alpha_deg, NACA_list(i).max_alpha_deg)
end









function index = get_index_of_NACA_from_NACA_list(NACA_list, wanted_NACA)
    
    for i = 1:length(NACA_list)
        if NACA_list(i).NACA == wanted_NACA
            index = i;
            return
        end
    end
    index = -1;
end
