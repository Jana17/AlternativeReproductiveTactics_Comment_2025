
#ifndef ART_sim_h
#define ART_sim_h

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <cfenv>
#include <climits>
#include <fstream>
#include <vector>
#include <string>
#include <string>
#include <sstream>
#include <algorithm>
#include "rndutils.hpp"



namespace ART_sim_code {

    using reng_type = rndutils::default_engine;
    enum sex { male, female };
    enum maturity { juvenile, mature };
    enum class morph_type { minor, major, female };
    enum vitality { alive, dead };
    enum class gene_type { env, threshold };
    enum inheritance_mode { blending, haploid, diploid };


    struct Param {

        inheritance_mode InheritanceMode = blending;
        bool evolving_threshold = true; //can be set to true or false
        const int max_time = 5001;//maximum runtime
        const int save_interval = 1;//which timesteps to output (here: every timestep)
        size_t used_seed;
        const size_t start_seed = 0;
        const size_t end_seed = 30;//number of replicates

        const double initial_env = 1.0;//the initial value of the environment
        const double initial_env_trait = 1.0;//the value of the "environmental genotype" in the initial population (the noise is added later)
        double initi_threshold = 0.25; //the value of the threshold in the initial population (noise is added later), or the global threshold if the threshold does not evolve. Note that if init_threshold=0 and evolving_threshold=false, there are no ARTs (all males are major males).

        double mut_rate = 1.0;//the mutation rate in the Mendelian inheritance scenario(s)
        double mut_stdev = 0.05;//the standard deviation of the distribution of mutational step sizes

        const double pop_size_K = 1000;//the carrying capacity
        const int initial_nr_m = 100; //initial number of males
        const int initial_nr_f = 100; //initial number of females
        const double cost_display = 2; //cost of ornamentation
        const int max_offspring = 6; //constant for calculating the max nr of offspring
        const double env_effect = 1; //constant to weigh how much condition is affected by the environment
        int group_size = 6; //size of the group of competing males
        double beta = 2; //strength of the female preference
        const double ART_success = 0.5; //probability that sneaking will succeed
        const double env_change_mean = 0.005;//mean environmental change
        const double env_change_sd = 0.005;//standard deviation of the distribution of environmental change step sizes
        const double resource_min = 0.0;//lower bound of the resource distribution
        const double resource_max = 1.0;//upper bound of the resource distribution
        const int AgeOfMaturity = 2; //individuals above this age (i.e. from age 3 onwards) are mature and can reproduce
        const double probBirth_Constant = 1; //constant for calculating the impact of density dependence on reproduction; default value 1
    };


    class rnd_j {
    public:
        reng_type rndgen;

        rnd_j() {
            std::random_device rd;
            reng_type rndgen_t(rd());
            rndgen = rndgen_t;
        }

        rnd_j(const Param& P) {
            rndgen = reng_type(P.used_seed);
            set_resource_dist(P.resource_min, P.resource_max);
        }

        // true or false with 50/50 probability:
        bool flip_coin() {
            return coin_flip(rndgen);
        }

        bool bernouilli(double p) {
            return std::bernoulli_distribution(p)(rndgen);
        }

        int poisson(double o) {
            return std::poisson_distribution<int>(o)(rndgen);
        }

        double uniform_real(double min, double max) {
            return std::uniform_real_distribution<double>(min, max)(rndgen);
        }

        double normalDist(double m, double stdev) {
            return std::normal_distribution<double>(m, stdev) (rndgen);
        }

        void set_resource_dist(double min, double max) {
            resource_dist = std::uniform_real_distribution<double>(min, max);
        }

        double generate_resources() {
            return resource_dist(rndgen);
        }

        sex get_random_sex() {
            if (flip_coin()) {
                return female;
            }
            return male;
        }

        double mutate_bounded_trait(const double& old_trait_value,
            const double& MRate,
            const double& stdev_mut,
            const double& lowbound,
            const double& highbound) {
            if (bernouilli(MRate)) {
                auto new_trait_value = normalDist(old_trait_value, stdev_mut);
                while (new_trait_value < lowbound || new_trait_value > highbound) {
                    new_trait_value = normalDist(old_trait_value, stdev_mut);
                }
                return new_trait_value;
            }
            return old_trait_value;
        }

        double mutate_unbounded_trait(const double& old_trait_value,
            const double& MRate,
            const double& stdev_mut) {
            if (bernouilli(MRate)) {
                auto new_trait_value = normalDist(old_trait_value, stdev_mut);
                return new_trait_value;
            }
            return old_trait_value;
        }

        // picks a random number in [0, n-1], useful when picking randomly from a vector.
        size_t random_number(size_t n) {
            if (n == 1) return 0;
            return std::uniform_int_distribution<>(0, static_cast<int>(n) - 1)(rndgen);
        }

    private:
        std::bernoulli_distribution coin_flip = std::bernoulli_distribution(0.5); // we can keep this fixed.
        std::uniform_real_distribution<double> unif = std::uniform_real_distribution<double>(0, 1.0);
        std::uniform_real_distribution<double> resource_dist;
    };


    struct gene {
        gene() {
            chromo = { 0, 0 };
            gene_type_ = gene_type::env;
            evolving_threshold = false;
        }

        gene(double init_value, gene_type gt,
            bool et, inheritance_mode ih, rnd_j& rnd,
            const Param& P) {
            chromo = { init_value, init_value };
            gene_type_ = gt;
            evolving_threshold = et;
            inheritance_mode_ = static_cast<inheritance_mode>(ih);

            if (gene_type_ == gene_type::env) {        // mutation for an environmental gene
                chromo[0] = rnd.mutate_unbounded_trait(chromo[0], 1.0, 0.25);
                chromo[1] = rnd.mutate_unbounded_trait(chromo[1], 1.0, 0.25);
            }
            else {
                if (evolving_threshold) {
                    chromo[0] = rnd.mutate_bounded_trait(chromo[0], 1.0, 0.05, 0.0, 1.0);
                    chromo[1] = rnd.mutate_bounded_trait(chromo[1], 1.0, 0.05, 0.0, 1.0);
                }
                check_bounds();
            }

            sd_mutate = P.mut_stdev;
            p_mutate = P.mut_rate;
        }

        gene(const gene& parent1,
            const gene& parent2,
            rnd_j& rnd) {

            evolving_threshold = parent1.evolving_threshold;
            gene_type_ = parent1.gene_type_;
            inheritance_mode_ = parent1.inheritance_mode_;
            p_mutate = parent1.p_mutate;
            sd_mutate = parent1.sd_mutate;

            if (inheritance_mode_ == blending) { // blending
                chromo[0] = 0.5 * (parent1.chromo[0] + parent2.chromo[0]);
            }
            else if (inheritance_mode_ == haploid) { // haploid
                chromo[0] = rnd.flip_coin() ? parent1.chromo[0] : parent2.chromo[0];
            }
            else if (inheritance_mode_ == diploid) { // diploid
                chromo[0] = rnd.flip_coin() ? parent1.chromo[0] : parent1.chromo[1];
                chromo[1] = rnd.flip_coin() ? parent2.chromo[0] : parent2.chromo[1];
            }
            else {
                throw "not matching inheritance mode";
            }

            if (gene_type_ == gene_type::threshold && evolving_threshold == true) {
                chromo[0] = rnd.mutate_bounded_trait(chromo[0], p_mutate, sd_mutate, 0.0, 1.0);
                if (inheritance_mode_ == diploid) {
                    chromo[1] = rnd.mutate_bounded_trait(chromo[1], p_mutate, sd_mutate, 0.0, 1.0);
                }
            }
            if (gene_type_ == gene_type::env) {
                chromo[0] = rnd.mutate_unbounded_trait(chromo[0], p_mutate, sd_mutate);
                if (inheritance_mode_ == diploid) {
                    chromo[1] = rnd.mutate_unbounded_trait(chromo[1], p_mutate, sd_mutate);
                }
            }
            if (gene_type_ == gene_type::threshold) check_bounds();
        }

        double get_trait_value() {
            if (inheritance_mode_ == blending) {
                return chromo[0];
            }
            if (inheritance_mode_ == haploid) {
                return chromo[0];
            }
            if (inheritance_mode_ == diploid) {
                return 0.5 * (chromo[0] + chromo[1]);
            }
            return 0;
        }

        void check_bounds() {
            if (gene_type_ == gene_type::threshold) {
                chromo[0] = chromo[0] > 1.0 ? 1.0 : chromo[0];
                chromo[1] = chromo[1] > 1.0 ? 1.0 : chromo[1];
                chromo[0] = chromo[0] < 0.0 ? 0.0 : chromo[0];
                chromo[1] = chromo[1] < 0.0 ? 0.0 : chromo[1];
            }
        }

        std::array<double, 2> chromo;
        gene_type gene_type_;
        bool evolving_threshold;
        inheritance_mode inheritance_mode_;

        double p_mutate;
        double sd_mutate;
    };

    struct Individual {

        sex S;
        maturity IsMature;

        //evolving gene loci:
        gene env_gene;
        gene threshold_gene;

        //traits
        double env_trait;
        double threshold;
        morph_type morph;
        size_t age;
        double res;
        double condition;
        double display_trait;
        vitality status;
        size_t rank;

        void set_other_properties(rnd_j& rnd,
            const Param& P,
            const double current_env) {
            if (age > P.AgeOfMaturity) { IsMature = mature; }
            else { IsMature = juvenile; }

            res = rnd.generate_resources();
            condition = res - P.env_effect * std::abs(env_trait - current_env);
            display_trait = 0;

            if (S == male) {
                if (condition < threshold) { morph = morph_type::minor; }
                else { morph = morph_type::major; }
            }
            else { morph = morph_type::female; }
            if (S == male && morph == morph_type::major && IsMature == mature) {
                display_trait = condition;
            }

            status = alive;
        }

        // default constructor with default values (e.g. individuals at initialisation):
        Individual(const Param& parameters,
            sex initial_sex, rnd_j& rnd) : S(initial_sex) {

            env_gene = gene(parameters.initial_env_trait,
                gene_type::env,
                parameters.evolving_threshold,
                parameters.InheritanceMode,
                rnd,
                parameters);

            threshold_gene = gene(parameters.initi_threshold,
                gene_type::threshold,
                parameters.evolving_threshold,
                parameters.InheritanceMode,
                rnd,
                parameters);

            env_trait = env_gene.get_trait_value();
            threshold = threshold_gene.get_trait_value();

            age = rnd.random_number(10) + 1; //random_nr(10) picks from [0,9], +1 is to avoid age=0

            set_other_properties(rnd, parameters, parameters.initial_env);
        }

        // constructor, creating Individual based on a mating event
        Individual(const Individual& parent1, 
            const Individual& parent2, 
            const Param& P,
            sex random_sex,
            double current_env,
            rnd_j& rnd) : S(random_sex) {


            env_gene = gene(parent1.env_gene, parent2.env_gene, rnd);
            env_trait = env_gene.get_trait_value();

            threshold_gene = gene(parent1.threshold_gene, parent2.threshold_gene, rnd);
            threshold = threshold_gene.get_trait_value();


            age = 0;
            set_other_properties(rnd, P, current_env);
        }


        void determine_morph() {
            if (morph != morph_type::female) {
                if (condition < threshold) {
                    morph = morph_type::minor;
                }
                else {
                    morph = morph_type::major;
                }
            }
        }

        void update(const double& AgeOfMaturity,
            const double& env_effect,
            const double& current_env) {
            age++;
            if (age > AgeOfMaturity) { IsMature = mature; }

            condition = res - env_effect * std::abs(env_trait - current_env);

            if (age == 1) {
                determine_morph();
            }

            if (morph == morph_type::major && IsMature == mature) {
                display_trait = condition;
            }
        }
    };


    void initialise(const Param& parameters,
        std::vector<Individual>& major_males,
        std::vector<Individual>& minor_males,
        std::vector<Individual>& females,
        std::vector<Individual>& juvenile_males,
        rnd_j& rng)
    {
        rng = rnd_j(parameters);
        major_males.clear();
        minor_males.clear();
        females.clear();
        juvenile_males.clear();
        sex focal_sexF = female;
        sex focal_sexM = male;

        //add inidividuals
        //remove individuals with a condition < 0
        for (size_t i = 0; i < parameters.initial_nr_f; i++) {
            Individual NewIndiv = Individual(parameters, focal_sexF, rng);
            if (NewIndiv.condition >= 0) { females.push_back(NewIndiv); }
        }
        for (size_t i = 0; i < parameters.initial_nr_m; i++) {
            Individual NewIndiv = Individual(parameters, focal_sexM, rng);
            if (NewIndiv.condition >= 0) {
                if (NewIndiv.IsMature == mature) {
                    if (NewIndiv.morph == morph_type::minor) {
                        minor_males.push_back(NewIndiv);
                    }
                    else {
                        major_males.push_back(NewIndiv);
                    }
                }
                else {
                    juvenile_males.push_back(NewIndiv);
                }
            }
        }
    }

    void update_pop(const Param& parameters,
        rnd_j& rng,
        const double& current_env,
        std::vector<Individual>& population) {

        ///Update the individual age, condition & display trait

        for (int i = 0; i < population.size(); i++) {
            population[i].update(parameters.AgeOfMaturity,
                parameters.env_effect,
                current_env);
        }

        std::vector<Individual> new_pop;
        new_pop.reserve(population.size());
        for (const auto& i : population) {
            if (i.condition >= 0.0) {//remove individuals whose condition has dropped under zero
                new_pop.push_back(i);
            }
        }
        population = new_pop;
    }

    void death_pop(const Param& parameters,
        rnd_j& rng,
        std::vector<Individual>& population,
        double current_density) {

        for (int i = 0; i < population.size(); i++) {
            double prob_death = 0.08 * (current_density + population[i].display_trait * parameters.cost_display +
                0.0154 * std::pow(population[i].age, 2) - 0.169 * population[i].age + 0.46);
            //NOTE: this is based on the R code - the equation in the methods section is different
            if (prob_death < 0) { prob_death = 0; }
            if (prob_death > 1) { prob_death = 1; } 
            if (rng.bernouilli(prob_death)) { population[i].status = dead; }
        }

        std::vector<Individual> new_pop;
        new_pop.reserve(population.size());
        for (const auto& i : population) {
            if (i.status == alive) {//remove individuals whose condition has dropped under zero
                new_pop.push_back(i);
            }
        }
        population = new_pop;
    }

    std::vector<size_t> get_mature_ids(const std::vector<Individual>& v) {
        std::vector<size_t> ids;
        ids.reserve(v.size()); // we reserve the maximum required size, but don't yet fill it
        for (size_t i = 0; i < v.size(); ++i) {
            if (v[i].IsMature == mature) ids.push_back(i);
        }
        return ids;
    }

    void reproduce(std::vector<Individual>& pop_females,
        std::vector<Individual>& pop_major_males,
        std::vector<Individual>& pop_minor_males,
        std::vector<Individual>& pop_juvenile_males,
        double current_env,
        const Param& parameters,
        rnd_j& rng) {

        //Step 1: preparing some variables outside of the big for-loop
        auto mature_major_male_IDs = get_mature_ids(pop_major_males);
        auto mature_minor_male_IDs = get_mature_ids(pop_minor_males);

        size_t NrMatureMales = mature_minor_male_IDs.size() + mature_major_male_IDs.size();
        if (NrMatureMales == 0) {//if there are no mature males, mating cannot take place, and we exit the function
            return;
        }

        std::vector<Individual> MatedMales;
        size_t current_group_size = parameters.group_size;
        if (current_group_size > mature_major_male_IDs.size()) {
            current_group_size = mature_major_male_IDs.size();
        } //can't have more males competing than the total number of mature major males

        //Alternative version of the mating system
        double ratio_minor_males = 1;
        if (mature_major_male_IDs.size() != 0) {
            ratio_minor_males= static_cast<double>(mature_minor_male_IDs.size()) / static_cast<double>(mature_major_male_IDs.size());
        }
        double expected_minor_matings = ratio_minor_males * parameters.ART_success;
        if (mature_major_male_IDs.size() == 0) { expected_minor_matings = 1; }//if no major males around, sneaking always successful

        //calculate density-dependent probability of reproduction
        double current_pop_size = static_cast<double>(pop_females.size()) + static_cast<double>(pop_major_males.size()) + static_cast<double>(pop_minor_males.size()) + static_cast<double>(pop_juvenile_males.size());
        double current_pop_density = current_pop_size / parameters.pop_size_K;
        double prob_birth = parameters.probBirth_Constant * (1 - current_pop_density);
        if (prob_birth > 1) { prob_birth = 1; }
        if (prob_birth < 0) { prob_birth = 0; }

        //Preparing the mating groups (which stay constant for the entire mating season
        size_t groupNr; //how many groups do we have
        if (current_group_size == 0) { groupNr = 0; }
        else { groupNr = mature_major_male_IDs.size() / current_group_size; }

        std::shuffle(mature_major_male_IDs.begin(), mature_major_male_IDs.end(), rng.rndgen);

        std::vector<std::vector<size_t>> MatingGroups; //contain positions of males in the vector of mature individuals
        for (int j = 0; j < groupNr; j++) {
            std::vector<size_t>NewGroup;
            for (size_t MaleNr = 0; MaleNr < current_group_size; MaleNr++) {
                NewGroup.push_back(mature_major_male_IDs[j * current_group_size + MaleNr]);
            }
            MatingGroups.push_back(NewGroup);
        }

        //let's check if there are males left over
        if (current_group_size != 0) {
            size_t leftOverMales = mature_major_male_IDs.size() % current_group_size;
            std::vector<size_t>NewGroup;
            for (size_t MaleNr = 0; MaleNr < leftOverMales; MaleNr++) {
                NewGroup.push_back(mature_major_male_IDs[groupNr * current_group_size + MaleNr]);
            }
            auto additionalMalesNr = current_group_size - leftOverMales; //fill up the last group with random major males, so that all groups are the same size
            std::vector<size_t>additionalMales;
            std::sample(mature_major_male_IDs.begin(), mature_major_male_IDs.end(), std::back_inserter(additionalMales),
                additionalMalesNr, rng.rndgen);
            NewGroup.insert(std::end(NewGroup), std::begin(additionalMales), std::end(additionalMales));
            if (NewGroup.size() != current_group_size) { std::cout << "ERROR - last group wrong size" << std::endl; }
            MatingGroups.push_back(NewGroup);
        }

        //Done with preparing all the variables and vectors!
        //Now we are ready to loop through all mature females and let them mate

        for (size_t i = 0; i < pop_females.size(); i++) {
            if (pop_females[i].IsMature) {//only mature females

                MatedMales.clear();

                //Step 2: mating with major males
                if (!MatingGroups.empty()) {//only if there are major males

                    //First, assign a mating group to this female at random
                    size_t GroupID = rng.random_number(MatingGroups.size());
                    auto ThisMatingGroup = MatingGroups[GroupID];

                    //Then, sort & rank males
                    std::vector<std::pair<double, int>>DisplayTraits;
                    DisplayTraits.reserve(ThisMatingGroup.size());
                    for (size_t j = 0; j < ThisMatingGroup.size(); j++) {
                        std::pair<double, size_t> NextMale;
                        NextMale.first = pop_major_males[ThisMatingGroup[j]].display_trait;
                        NextMale.second = j;//position in the current mating group vector
                        DisplayTraits.push_back(NextMale);
                    }
                    std::sort(DisplayTraits.begin(), DisplayTraits.end());
                    for (size_t j = 0; j < DisplayTraits.size(); j++) {
                        //male with the largest ornament needs to have rank=1,
                        //male with the smallest orament needs to have rank=group size
                        pop_major_males[ThisMatingGroup[DisplayTraits[j].second]].rank = current_group_size - j;
                    }

                    std::vector<double>attractiveness; //a measure of likelihood of being chosen
                    attractiveness.reserve(ThisMatingGroup.size());
                    int ConstSum = 0;
                    for (int j = 0; j < ThisMatingGroup.size(); j++) {
                        ConstSum = ConstSum + (1 / pow((j + 1), parameters.beta));
                    }
                    for (int j = 0; j < ThisMatingGroup.size(); j++) {
                        double new_attractiveness = (1 / pow(pop_major_males[ThisMatingGroup[j]].rank, parameters.beta)) / ConstSum;
                        attractiveness.push_back(new_attractiveness);
                    }
                    std::discrete_distribution<int> weightedLottery(attractiveness.begin(), attractiveness.end());
                    int ChosenMalePosition = weightedLottery(rng.rndgen);//choose a male from the lottery (weighted by attractiveness)
                    MatedMales.push_back(pop_major_males[ThisMatingGroup[ChosenMalePosition]]);//add this chosen major male to MatedMales
                    if (pop_major_males[ThisMatingGroup[ChosenMalePosition]].IsMature == juvenile) { std::cout << "ERROR! Mated male is not mature. " << std::endl; }
                }

                //Step 3: mating with minor males
                //Alternative versions of the mating system
                int NrMinorMates = rng.poisson(expected_minor_matings);
                for (int SneakerMaleCount = 0; SneakerMaleCount < NrMinorMates; SneakerMaleCount++) {
                    int MinorMale_Pos = rng.random_number(mature_minor_male_IDs.size());
                    int MinorMale_ID = mature_minor_male_IDs[MinorMale_Pos];
                    MatedMales.push_back(pop_minor_males[MinorMale_ID]);//add the chosen minor male to MatedMales
                    if (pop_minor_males[MinorMale_ID].IsMature == juvenile) { std::cout << "ERROR! Mated male is not mature. " << std::endl; }
                }

                //the female has now chosen her mate(s)!


                //Step 4: produce offspring
                if (MatedMales.size() != 0) {//only produce offspring if mating has taken place

                    int MaxIndividualOffsprNr = std::round(parameters.max_offspring * pop_females[i].condition);//fecundity depends on condition
                    int OffspCounter = 0;
                    for (int Offsp_i = 0; Offsp_i < MaxIndividualOffsprNr; Offsp_i++) {
                        if (rng.bernouilli(prob_birth)) { OffspCounter++; }//fecundity depends on population density
                    }
                    int ActualOffsprNr = OffspCounter;

                    for (int j = 0; j < ActualOffsprNr; j++) {//now, loop through each of the offspring that this female is producing to add it to the population
                        size_t WhichFather = rng.random_number(MatedMales.size());//chose one of the possible fathers at random (if the female has mated more than once)
                        sex random_sex = rng.get_random_sex();
                        Individual NewOffspring(pop_females[i], MatedMales[WhichFather], parameters, random_sex, current_env, rng);//produce the offspring individual
                        if (NewOffspring.condition >= 0) {//if the new offspring has a condition that is not below zero, add it to the population

                            if (NewOffspring.morph == morph_type::female) {
                                pop_females.push_back(NewOffspring);
                            }
                            else {
                                pop_juvenile_males.push_back(NewOffspring);
                            }
                        }
                    }
                }
            }
        }
    }

    void check_major_males(const std::vector<Individual>& v) {
        for (const auto& i : v) {
            if (i.morph == morph_type::minor) {
                std::cout << "error! found minor male in major vector";
            }
        }
    }

    void update_juveniles(std::vector<Individual>& pop_juvenile_males,
        std::vector<Individual>& pop_major_males,
        std::vector<Individual>& pop_minor_males) {
        std::vector<Individual> still_juvenile;
        for (auto& i : pop_juvenile_males) {
            if (i.IsMature == mature) {
                if (i.morph == morph_type::major) {
                    pop_major_males.push_back(i);
                }
                else if (i.morph == morph_type::minor) {
                    pop_minor_males.push_back(i);
                }
                else {
                    std::cout << "this can not happen!\n";
                }

            }
            else {
                still_juvenile.push_back(i);
            }
        }

        pop_juvenile_males = still_juvenile; // these are still juvenile
    }

    void run() {
        Param parameters;

        //prepare the output file
        std::string TEvoString = "NoTEvo";
        if (parameters.evolving_threshold == true) { TEvoString = "TEvo"; }
        std::ofstream ofs1("ART_Fig1_" + TEvoString + "_initT" + std::to_string(parameters.initi_threshold) + "_AlternativeMatingSystem.csv");
        ofs1 << "Time" << "," << "Rep" << "," << "PopSizeF" << "," << "PopSizeM" << ","
            << "PopSizeM_Maj" << "," << "PopSizeM_Min" << "," << "PopSizeM_Juv" << "," << "PopSize" << ","
            << "CurrentEnv" << "," << "avg_env" << "," << "avg_thr" << "," << "avg_res" << ","
            << "avg_condition" << "," << "avg_majM_condition" << "," << "avg_minM_condition" << ","
            << "avg_juvM_condition" << "," << "avg_F_condition" << "," << "avg_majM_display" << ","
            << "iniT" << "," << "MRate" << "," << "MStdDev" << "," << "Extinction";
        ofs1 << '\n';



        std::vector<Individual> pop_females;
        std::vector<Individual> pop_major_males;
        std::vector<Individual> pop_minor_males;
        std::vector<Individual> pop_juvenile_males;
        rnd_j rng;


        for (auto current_seed = parameters.start_seed; current_seed < parameters.end_seed; ++current_seed) {//looping through the replicates
            parameters.used_seed = current_seed;
            std::cout << "Seed " << parameters.used_seed << std::endl;
            initialise(parameters, pop_major_males, pop_minor_males, pop_females, pop_juvenile_males, rng);//initialise the population
            double current_env = parameters.initial_env;//initialise the environment

            for (int time = 0; time < parameters.max_time; time++) {//loop through all timesteps
                ///Update the environment & the population
                if (time < 125) { current_env = parameters.initial_env + rng.uniform_real(-0.02, 0.02); }
                else { current_env = current_env + rng.normalDist(parameters.env_change_mean, parameters.env_change_sd); }

                update_pop(parameters, rng, current_env, pop_females);
                update_pop(parameters, rng, current_env, pop_major_males);
                update_pop(parameters, rng, current_env, pop_minor_males);
                update_pop(parameters, rng, current_env, pop_juvenile_males);
                update_juveniles(pop_juvenile_males, pop_major_males, pop_minor_males);

                ///Death - depending on age & display costs
                double current_pop_size = static_cast<double>(pop_females.size()) + static_cast<double>(pop_major_males.size()) + static_cast<double>(pop_minor_males.size()) + static_cast<double>(pop_juvenile_males.size());
                double current_density = current_pop_size / parameters.pop_size_K;
                death_pop(parameters, rng, pop_females, current_density);
                death_pop(parameters, rng, pop_major_males, current_density);
                death_pop(parameters, rng, pop_minor_males, current_density);
                death_pop(parameters, rng, pop_juvenile_males, current_density);

                ///Finding a mate & reproducing
                reproduce(pop_females, pop_major_males, pop_minor_males, pop_juvenile_males, current_env, parameters, rng);
                check_major_males(pop_major_males);

                ////Save the data
                bool Extinction = false;
                if (pop_females.size() < 1) { Extinction = true; }
                int NrMalesAlive = pop_major_males.size() + pop_minor_males.size() + pop_juvenile_males.size();
                if (NrMalesAlive < 1) { Extinction = true; }
                if (time % parameters.save_interval == 0 || Extinction) {
                    //if (time != 0) {
                        double avg_env = 0;
                        double avg_thr = 0;
                        double avg_res = 0;
                        double avg_condition = 0;
                        double avg_majM_condition = 0;
                        double avg_minM_condition = 0;
                        double avg_juvM_condition = 0;
                        double avg_F_condition = 0;
                        double avg_majM_display = 0;

                        for (const auto& i : pop_females) {
                            avg_env += i.env_trait;
                            avg_thr += i.threshold;
                            avg_res += i.res;
                            avg_condition += i.condition;
                            avg_F_condition += i.condition;
                        }

                        for (const auto& i : pop_major_males) {
                            avg_env += i.env_trait;
                            avg_thr += i.threshold;
                            avg_res += i.res;
                            avg_condition += i.condition;
                            avg_majM_condition += i.condition;
                            avg_majM_display += i.display_trait;
                        }

                        for (const auto& i : pop_minor_males) {
                            avg_env += i.env_trait;
                            avg_thr += i.threshold;
                            avg_res += i.res;
                            avg_condition += i.condition;
                            avg_minM_condition += i.condition;
                        }

                        for (const auto& i : pop_juvenile_males) {
                            avg_env += i.env_trait;
                            avg_thr += i.threshold;
                            avg_res += i.res;
                            avg_condition += i.condition;
                            avg_juvM_condition += i.condition;
                        }


                        double popsize = static_cast<double>(pop_females.size()) + static_cast<double>(pop_major_males.size()) + static_cast<double>(pop_minor_males.size()) + static_cast<double>(pop_juvenile_males.size());
                        double popsize_f = static_cast<double>(pop_females.size());
                        double popsize_m = static_cast<double>(pop_major_males.size()) + static_cast<double>(pop_minor_males.size());
                        double popsize_majM = static_cast<double>(pop_major_males.size());
                        double popsize_minM = static_cast<double>(pop_minor_males.size());
                        double popsize_juvM = static_cast<double>(pop_juvenile_males.size());

                        double popsize_mult = 1.0 / (static_cast<double>(pop_females.size()) + static_cast<double>(pop_major_males.size()) + static_cast<double>(pop_minor_males.size()) + static_cast<double>(pop_juvenile_males.size()));
                        double popsize_f_mult = 1.0 / (static_cast<double>(pop_females.size()));
                        double popsize_majM_mult = 1.0 / (static_cast<double>(pop_major_males.size()));
                        double popsize_minM_mult = 1.0 / (static_cast<double>(pop_minor_males.size()));
                        double popsize_juvM_mult = 1.0 / (static_cast<double>(pop_juvenile_males.size()));

                        avg_env *= popsize_mult;
                        avg_thr *= popsize_mult;
                        avg_res *= popsize_mult;
                        avg_condition *= popsize_mult;
                        avg_F_condition *= popsize_f_mult;
                        avg_majM_condition *= popsize_majM_mult;
                        avg_minM_condition *= popsize_minM_mult;
                        avg_juvM_condition *= popsize_juvM_mult;
                        avg_majM_display *= popsize_majM_mult;

                        ofs1 << time << "," << current_seed << "," << popsize_f << "," << popsize_m << ","
                            << popsize_majM << "," << popsize_minM << "," << popsize_juvM << "," << popsize << ","
                            << current_env << "," << avg_env << "," << avg_thr << "," << avg_res << ","
                            << avg_condition << "," << avg_majM_condition << "," << avg_minM_condition << ","
                            << avg_juvM_condition << "," << avg_F_condition << "," << avg_majM_display << ","
                            << parameters.initi_threshold << "," << parameters.mut_rate << "," << parameters.mut_stdev << "," << Extinction
                            << '\n';
                    //}
                }

                ///account for extinctions
                if (Extinction) { time = parameters.max_time + 1; }
            }
        }

        ofs1.close();

    }


}   // end namespace new_code



#endif /* ART_sim_h */
