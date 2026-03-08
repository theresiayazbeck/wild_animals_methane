%Get mammal biomass estimation from
%https://gitlab.com/milo-lab-public/mammal_biomass/ and https://gitlab.com/milo-lab-public/mammal_biomass_since_1850
wild_animals_data = readtable('\mammal_biomass-master\results\wild_land_mammal_biomass_inluding_populations.csv');
SpeciesSummary = readtable('\mammal_biomass-master\data\all_species_list.csv');
wild_terrestrial_biomass = readtable('\mammal_biomass_since_1850-main\results\wild_land\wild_terrestrial_biomass.csv');

% Number of Monte Carlo simulations
nIter = 10000;

%Get ruminants
wild_animals_data_formethane = innerjoin(wild_animals_data,SpeciesSummary);
cetartiodactyl = strcmp(wild_animals_data_formethane.Order,'CETARTIODACTYLA');

swine_genera = {'Babyrousa';'Catagonus';'Choeropsis';'Hippopotamus';'Hylochoerus';...
    'Pecari';'Phacochoerus';'Porcula';'Potamochoerus';'Sus';'Tayassu'};
genus = extractBefore(wild_animals_data_formethane.binomial, " ");
swine = ismember(genus,swine_genera);

%Specific species to exclude for error generated later in R code
ruminants = (cetartiodactyl & ~swine);
ruminants_data = wild_animals_data_formethane(ruminants,:);

final_binomials = readtable('final_binomials.csv');
final_binomials.Properties.VariableNames(1) = "binomial";
ruminants_data = innerjoin(ruminants_data,final_binomials);

nRuminants = height(ruminants_data);

%Get relative biomass compared to reference biomass for yearly estimates
wild_terrestrial_biomass.relative_biomass = wild_terrestrial_biomass.biomassBaseline_Mt_/(sum(wild_animals_data_formethane.biomass_g)*10^(-12));

%Initialize methane emissions arrays
total_emissions_ruminants = zeros(nRuminants,25);
CI_low = zeros(nRuminants,25);
CI_high = zeros(nRuminants,25);

%Define regression parameters for Pérez-Barbería (2017)
%Intercept
a_estimate = -0.573;
a_se = 0.036;

%logBM
b_estimate = 1.075;
b_se = 0.017;

% Assuming log-normal biomass distribution using 2.5% and 97.5% CI
for Y = 2000:2024

    scaling_factor = wild_terrestrial_biomass.relative_biomass(wild_terrestrial_biomass.Year == Y);
    pop_estimate = ruminants_data.estimated_population*scaling_factor;
    pop_ub = ruminants_data.estimated_population_97pt5*scaling_factor;
    pop_lb = ruminants_data.estimated_population_2pt5*scaling_factor;
    pop_se = (pop_ub-pop_lb)/(2*1.96);

    pop_mu = (log(pop_lb) + log(pop_ub)) / 2;
    pop_sigma = (log(pop_ub) - log(pop_lb)) / (2*1.96);

    for i = 1:nRuminants
                       
        %Get parameter samples
        a = normrnd(a_estimate,a_se,nIter,1);
        b = normrnd(b_estimate,b_se,nIter,1);
        pop = lognrnd(pop_mu(i),pop_sigma(i),nIter,1);

        %Pérez-Barbería (2017)
        emission_sample = pop.*(exp(a + b .* log(ruminants_data.AdultBodyMass_g(i)*0.001)) * 1e-12 * 365);%in Tg

        % Total emissions for this ruminant
        total_emissions_ruminants(i,Y-1999) = mean(emission_sample);
        CI_low(i,Y-1999) = prctile(emission_sample, 2.5);
        CI_high(i,Y-1999) = prctile(emission_sample, 97.5);

    end

end

methane_per_species_2024 = total_emissions_ruminants(:,end);
methane_per_species_2024 = table(ruminants_data.binomial,methane_per_species_2024*10^12,'VariableNames',...
    ["binomial","methane_g_yr"]);

disp(['2024 Total Methane = ' num2str(sum(methane_per_species_2024.methane_g_yr)*10^(-12)) 'Tg']);

%Write all years into one table
years = 2000:2024;
yearVarNames = strcat("y_", string(years));

methane_per_species = array2table(total_emissions_ruminants*10^(12),'VariableNames', yearVarNames);

methane_per_species = addvars(methane_per_species,ruminants_data.binomial,'Before', 1, ...
    'NewVariableNames', 'binomial');

disp(['Total Methane = ' num2str(sum(total_emissions_ruminants)) 'Tg']);

writetable(methane_per_species,'methane_per_species_ruminants_perez2017_10000Iter.csv')

%Write total methane emissions with uncertainty in one file
methane_total = table((2000:1:2024)',sum(total_emissions_ruminants)',sum(CI_low)',sum(CI_high)',...
    'VariableNames',{'year','methane_Tg','CI_lb','CI_ub'});

writetable(methane_total,'total_methane_Tg.csv')



