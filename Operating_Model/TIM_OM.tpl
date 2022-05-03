
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Spatial Operating Model built by Daniel Goethel (daniel.goethel@noaa.gov)
// Contributed to by Katelyn Bosley, Jon Deroba, Dana Hanselman, Amy Schueller, Brian Langseth, Aaron Berger
// Part of the Spatial Processes and Stock Assessment Methods (SPASAM) Working Group
// Multistage Stock-Recruit Relationships and IBM Functionaliyt Added as Part of the
// Spatially Integrated Life Cycle (SILC) Working Group
// Code and Supporting Materials Can Be Downloaded from https://github.com/dgoethel/Spatially-Integrated-Life-Cycle-SILC-Model
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// This model can be used by itself to simulate various population structures/dynamics and/or search for F_MSY
// OR the resulting .rep file can be read directly into the associated TIM_EM.tpl estimation model as a .dat file through an R wrapper
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// SEE README IN DAT FILE FOR LIMITATIONS ON ESTIMATION MODEL IN REGARDS TO AXIS OF ESTIMATION/VARIATION IN SPATIAL PARAMETERS
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

GLOBALS_SECTION
  #include "admodel.h"
  #include "qfclib.h"
  #define EOUT(var) cout << #var <<" "<<var<<endl;

TOP_OF_MAIN_SECTION
  arrmblsize=500000000;
  gradient_structure::set_MAX_NVAR_OFFSET(5000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000000);

DATA_SECTION

  init_number phase_F
   // phase for estimating F
   // must be turned on (==1) if F_type==3 AND MSY_model_type_switch==3 (performing F_MSY search)
    
  init_number F_guess
   // starting value for F search if F_type==3
    
  init_number phase_dummy
   // phase for dummy parameter in the OM, used when OM does no estimation (i.e., all phases<0) because ADMB must have at least 1 active parameter to run
   // DEFAULT: set ==1 because not actually estimating any parameters in OM
   // must be turned on (==1) if F_type!=3
   // must be turned off (==(-1)) if F_type=3

  init_int ph_dummy_EM
   // phase for dummy parameter for the EM, used when EM does no actual estimation (i.e., all phases<0) because ADMB must have at least 1 active parameter to run
   // DEFAULT: set ==(-1) because EM will be estimating parameters so dummy variable unnecessary
   
////////////////////////////////////////////////////////////////////////////////////
///// OM MODEL STRUCTURE INPUTS /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

  init_number OM_structure
   // Determines the overarching population structure assumed in the operating model (OM)
   //==0 the OM is panmictic
   //==1 the OM is metamictic/spatial heterogeneity (ie 1 population with multiple areas/subpopulations)
   //==2 the OM is metapop (multiple subpopulations each with their own vital rates, but reproductive mixing, such that fish immediately adopt new vital rates when join new population)
   //==3 the OM is natal homing/overlap/home fidelity (vital rates are based on natal population; all spawners can return to spawn or only fraction...based on input overlap parameters)

  init_number natal_homing_switch
   // determines how SSB is tallied (i.e., follow natal homing assumption that fish must return to natal population to add to SSB), mainly check on OM_structure setting
   //==0 no natal homing, use if OM_structure!=3 (SSB is sum of SSB in population regardless of natal origin; weight/mat/fecund/ are based on current population not natal population)
   //==1 do natal homing, use if OM_structure==3 (a fish only adds to SSB if it is in its natal population at spawning time; weight/mat/fecund/ are based on natal population)
      // natal homing assumes genetic based life history and contribution to SSB (i.e., natal homing and no demographic mixing), natal_homing_switch==0 assumes demographic mixing (e.g. metapopulations where life history is more location based)

  init_number spawn_return_switch
   // determines whether an instantaneous return spawning migration to natal population occurs IF natal_homing_switch==1 (i.e., when OM_structure==3 and natal homing assumed)
   //==0 if natal_homing_switch==1, then ONLY fish that are currently in natal population at time of spawning add to SSB
   //==1 if natal_homing_switch==1, then a fraction of fish return to natal population to spawn (instantaneous migration to natal population and back at time of spawning) based on spawn_return_prob; weight/mat/fecund/ are based on natal population

  init_int nages
   // number of ages (last age is a plus group)
   
  init_int nyrs
   // number of years in the model
   // NOTE: for reference point sims, MSY_model_type_switch>0,  (e.g., F_MSY search), should set this sufficiently high (>100) to approximate equilibrium
 
  init_int npops
   // number of populations
   // NOTE: if OM_structure==0 (panmictic), npops should ==1
   // NOTE: if OM_structure==1 (spatial heterogeneity), npops typically ==1 (but could have multiple populations with no movement among pops; typically would just use metapop structure though)
   // NOTE: if OM_structure>1 (natal homing or metapopulation), npops should >1, otherwise have either panmictic or spatial heterogeneity

  !! int np=npops;
    // turn npops into an integer that can be used as an index for following arrays
  
  init_ivector nregions(1,np)
   // number of regions per population
   // SEE README IN DAT FOR POTENTIAL MISSPECIFICATION OF EM IF DIFFERING NUMBER OF REGIONS IN EACH POPULATION
   // NOTE: if OM_structure==0 (panmictic), nregs should ==1, otherwise have spatial heterogeneity
   // NOTE: if OM_structure==1 (spatial heterogeneity), nregs should be >1, otherwise have panmictic
   // NOTE: if OM_structure>1 (natal homing or metapopulation), nregs can be =>1
 
  init_ivector nfleets(1,np)
   // number of fishing fleets in each region for each population
   // SEE README IN DAT FOR POTENTIAL MISSPECIFICATION OF EM IF DIFFERING NUMBER OF FLEETS IN EACH POPULATION
   
  init_ivector nfleets_survey(1,np)
   // number of survey fleets in each region for each population
   // SEE README IN DAT FOR POTENTIAL MISSPECIFICATION OF EM IF DIFFERING NUMBER OF FLEETS IN EACH POPULATION

////////////////////////////////////////////////////////////////////////////////////
///// INDICES FOR Random Number Generators /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

   // use these to input the MAX values across ALL runs/scenarios for a simulation experiment so that always have consistent RNG values
   // ie these may be greater than the actual dimensions above for any given run, if the sim experiment is exploring variation in model/population structure
   // by using MAX value across all sim scenarios, ensures that associated error/stochasticity in parameter values or observed data is the same for that pop/reg/fleet/age/data source across all runs/scenarios
   // ie, data/parameters in pop 2, reg 2 will always have same set of RNGs, but reg 3 will be different from reg 2 (even though there may not be a reg 2 in this sim scenario)
   // just a way to ensure consistentcy across sim scenarios/runs
  init_int max_pops
  init_int max_regs
  init_int max_ages
  init_int max_yrs
  init_int max_flts
  init_int max_surv_flts
  init_int max_tag_yrs

////////////////////////////////////////////////////////////////////////////////////
///// INDICES for ragged arrays /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

   //NOTE: issue with using type 'init' for indices of ragged arrays in ADMB so need to assign each index value as an integer
  int ny
  !! int ny=nyrs;
  !! int na=nages;
  !! ivector nreg=nregions;
  !! ivector nf=nfleets;
  !! ivector nfs=nfleets_survey;

////////////////////////////////////////////////////////////////////////////////////
////////////// OM Model SWITCHES ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

  init_number larval_move_switch
   // Changes the type of larval movement pattern (sets age class 0 movements)
   //==(-1) larval (pre-age-1) movement is determined in same way as other ages (e.g., based on move_switch) 
   //==0 no movement
   //==1 input movement based on adult movement matrix
   //==2 movement within population only based on population/region specific input residency (symmetric)
   //==3 symmetric movement but only allow movement within a population (ie regions within a population) not across populations
   //==4 symmetric movement across all populations and regions
   //==5 allow movement across all regions and populations, based on population/region specific input residency (symmetric off-diag)
   //==6 movement based on input IBM connectivity matrices (NOT IMPLEMENTED YET)

  init_number first_post_settle_move_switch
   // allows different movement patterns for newly settled fish up until age_first_move
   // assume no movement until age_first_move
   // NOTE: NOT IMPLEMENTED YET, in development as part of multistage SRR for SILC project
   //==(-1) first movement is determined in same way as other ages (e.g., based on move_switch) 
   //==0 no movement
   //==1 input first movement rates using first_move_input
   //==2 return to natal population based on first_return parameter (proportion of fish that return during first post-settlement movement)

  init_number move_switch
   // Sets the type of adult movement pattern (sets age class>0 movements)
   //==0 no movement
   //==1 input age-based movement
   //==2 movement within population only (e.g., among regions within population), based on input pop,reg,age specific residency (symmetric off diagonal)
   //==3 symmetric movement but only allow movement within a population (ie regions within a population) not across populations
   //==4 symmetric movement across all populations and regions
   //==5 allow movement across all regions and populations, based on input pop,reg,age specific residency term (symmetric off diagonal)
   //==6 natal return (for natal_homing_switch==1 ONLY), no movement of fish until return_age when a certain fraction, return_probability, of fish make return migration to natal population (eg, an ontogenetic migration); all fish remain in given population for remainder of lifepsan after return_age
   //==7 larvae stay in the population that they settle in (i.e., for natal homing/overlap, do not return to natal population)
      // if adult movement==0 for natal homing would return to natal population because natal residency is 100% and use natal movement rates (not current population movement rates like with metapopulation/random movement)
   //==8 density dependent movement based on relative biomass among potential destination population/regions, partitions (1-input_residency) based on a logistic function of biomass in current population/region and 'suitability' of destination population/regions
      // uses use_input_Bstar switch
      // DD MOVEMENT CAN BE AGE BASED OR CONSTANT ACROSS AGES...FOR AGE BASED MAKE SURE DD_move_age_switch==1, FOR AGE-INVARIANT DD_move_age_switch==0
   //==9 use input T_year to allow T to vary by year
   //==10 use T_FULL_Input to input time and age varying T

  init_number rand_move_larval
   // adjust  movement to include random variation based on input sigma and lognormal random variable
   // IF larval_move_switch==0, then NO RANDOMNESS IS EMPLOYED (I.E., RESIDENCY STAYS AT 100%)
   //==0 no randomness, just use movement from T calcs
   //==1 add randomness to T (bounded so movement proportion cannot exceed 1 or go below 0)

  init_number rand_move
   // adjust  movement to include random variation based on input sigma and lognormal random variable
   // IF move_switch==0, then NO RANDOMNESS IS EMPLOYED (I.E., RESIDENCY STAYS AT 100%)
   //==0 no randomness, just use movement from T calcs
   //==1 add randomness to T (bounded so movement proportion cannot exceed 1 or go below 0)
  
  init_number DD_move_age_switch
   // Allow age-based movement when using DD movement (Y/N)
   //==0 no age-based DD movement (assumes movement based on total biomass relative to Bstar) 
   //==1 DD movement is age-based (assumes movement based on age-specfic biomass relative to age-specific Bstar)

  init_number use_input_Bstar
   //==0 set Bstar for DD movement equal to SSB_zero*SSB_zero_appor  #######NOT FULLY IMPLEMENTED YET##########
      // (if nreg==1, Bstar=SSB0), **NOTE** for Rec_type!=2 (not BH), defaults to using input_Bstar since no SSB0 calcs are done 
   //==1 use input_Bstar
      // This will be typical approach for implementing DD movement since SSB0 calcs are not well defined for spatial models

  init_number select_switch
   // determine how fishery selectivity is simulated
   //==0 input selectivity
   //==1 logistic selectivity based on input sel_beta1 and sel_beta2
   //==2 double logistic selectivity based on input sel_beta1, sel_beta2, sel_beta3 and sel_beta4  // NOT TEsted YET

  init_number select_switch_survey
   // determine how survey selectivity is simulated
   //==0 input selectivity
   //==1 logistic selectivity based on input sel_beta1 and sel_beta2
   //==2 double logistic selectivity based on input sel_beta1, sel_beta2, sel_beta3 and sel_beta4  // NOT Tested  YET

  init_number F_switch
   // determine how F is simulated
   //==1 input_F by pop, reg, fleets
   //==2 input_F_MSY
   //==3 Estimate F_MSY
   //==4 input_F_const is split evenly among populations (each fleet uses population F)
   //==5 input_F_const is is split evenly among all regions (each fleet uses region F)
   //==6 input_F_const is split evenly among fleets
   //==7 F devs about input_F based on sigma_F
   //==8 random walk in F  //NOT Tested YET
   //==9 dunce cap F with random devs; F increases to a peak half way through time series then decrease back to min over remainder of time series
  
  init_number init_abund_switch
   // determines how the abundance is calculated in the first year of the simulation
   //==0 input init abundance
   //==1 decay from R_ave

  init_number maturity_switch_equil
   // MORE WORK IS NEEDED TO REFINE THESE CALCULATIONS AND DEAL WITH SPATIAL REFERENCE POITNS!
   // SSB0 must be calculated to determine stock-recruit function (if only know steepness and R0 for the population)
   // Use equilibrium SPR calcs to get SSB0, but to do so requires vital rates (maturity, weight), which are typically constant across a population
   // With multiple regions within a pop each with different vitals, must make assumption regarding the proportional contribution of each region's demograhics to equil SSB
   //==0  equal by region, assume equal (average) contributions to SSB0 by each region
   //==1 weighted average, use input equil_ssb_apportion to determine proportional contribution to equil vital rates by region

  init_number SSB_type
   // units of spawning stock biomass
   //==1 fecundity based SSB
   //==2 weight based SSB

  init_number Rec_type
   // form of the stock-recruit relationship
   //==1 stock-recruit relationship assumes an average value based on R_ave
   //==2 Beverton-Holt population-recruit functions based on population-specific input steepness, R0 (R_ave), M, and weight
         // NOTE: SRR DOES NOT TAKE INTO ACCOUNT SPATIAL DYNAMICS (IE, MOVEMENT AMONG POPULATIONS OR REGIONS); SEE README TOPIC AT TOP OF DAT FILE
   //==3 environmental recruitment - sine fucntion based on amplitude and frequency

  init_number recruit_devs_switch
   // determine whether recruitment deviations are incorporated
   //==0 use stock-recruit relationphip directly
   //==1 allow lognormal error around SR curve (i.e., include randomness based on input sigma_recruit)

  init_number recruit_randwalk_switch
   //==0 no random walk recruitment deviations (DEFAULT)
   //==1 have random walk lognormal recruitment deviations (requires recruit_devs_switch==1)....HAS NOT BEEN FULLY IMPLEMENTED OR TESTED

  init_number apportionment_type
   // determines how recruits are apportioned to regions within a population
   // because stock-recruit relationships are assumed only at the population level, if have more than 1 region in a population, must make assumption about how to apportion/assign recruits to each region
   // typically used with OM_structure==1, but also applies if nreg>1 in any pop for OM_structure==2 AND 3
   //==-1 no recruitment apportionment to regions within a population (each region within a population gets full amount of recruits from SR curve); WOULD NOT SUGGEST USING Because leads to sum(recruits)>SRR
   //==0 apportionment to each region is based on relative SSB in region compared to population SSB
   //==1 input apportionment (DEFAULT for OM_structur==2 AND 3)
   //==2 recruits are apportioned equally to each region within a population
   //==3 recruits are apportioned in a completely random manner with uniform equilibrium distribution
   //==4 recruits are apportioned stochastically with normal error surrounding the input proportions...uses the multivariate logistic random variables (Cox and Krunland 2008, FIsheries Research); SUGGESTED METHOD FOR RANDOM APPORTIONMENT
   //==5 recruits are approtioned based on theoretical enviormental phase shift....HAS NOT BEEN FULLY IMPLEMENTED OR TESTED
  
  init_number use_stock_comp_info_survey
   // Determines whether it is assumed that info (stock composition data) is available to determine natal origin for survey age composition data
   //==0 calc OBS survey age comps by pop/area (summed across natal population); assumes no stock compisition data is available; typical for most situations when no genetic/otolith analysis is available
   //==1 calc OBS survey age comps by natal population within each area; assumes stock composition data is available
      // currently assumes that observation error is just that associated with typical observation error based on assumed distribution/CV of given data source
      // in future should incorporate observation error that mimics error more typically associated with genetic/otolith analysis

  init_number use_stock_comp_info_catch
   // Determines whether it is assumed that info (stock composition data) is available to determine natal origin for catch age composition data
   //==0 calc OBS survey age comps by pop/area (summed across natal population); assumes no stock compisition data is available; typical for most situations when no genetic/otolith analysis is available
   //==1 calc OBS survey age comps by natal population within each area; assumes stock composition data is available
      // currently assumes that observation error is just that associated with typical observation error based on assumed distribution/CV of given data source
      // in future should incorporate observation error that mimics error more typically associated with genetic/otolith analysis

////////////////////////////////////////////////////////////////////////////////////
////////////// Tagging Data SWITCHES ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

  init_number do_tag
   // Turn off tagging in OM to save computation time if not going to use tagging data in sim/est models
   //==0 do not calculate tagging data
   //==1 calculate tagging data

  init_number number_tags_switch
   // determines how total tags are determined and how they are distributed across regions
   //==(-2) fixed tags distributed evenly,  use input_total_tags and distribute evenly across regions
   //==(-1) fixed tagging protocol distributed according to survey abundnance, use input_total tags and distribute across regions according to survey abundance
   //==0 fixed tagging protocol, use input tags by year, pop, reg
   //==1 tagging based on total abundance, total tags are based on the input fract of total abundance tagged, frac_abund_tagged, and is distributed evenly across regions
   //==2 tagging based on regional abundnace, regional tags are based on the input fract of regional abund tagged, frac_abund_tagged
   //==3 tagging based on total abundance and distributed according to regional survey biomass, tags are based on the input of total abund tagged, frac_abund_tagged, and is distributed based on survey prop of abundance in each region
   //==4 opportunitstic tagging based on target min/max tags, regional tags are randomly (uniformly) distributed based on input max and min number of tags per region
   //==5 opportunistic tagging based target min/max but with no tagging in some years AND areas, regional tags are randomly (uniformly) distributed based on input max and min number of tags per region AND prob that tagging occurs in a given year/region is based on opport_tag_prob and a uniform distributed random number
   //==6 fixed opportunistic tagging distributed by regional survey biomass but with no tagging in some years, tags are input and distributed according to regional survey abundance AND
      // prob that tagging occurs in a given year is based on opport_tag_year_prob and a uniform distributed random number (random tagging by year)
   //==7 opportunistic tagging based on min/max tags and no tagging in some years and regions, yearly tags are randomly (uniformly) distributed based on max and min number of tags per region AND
      // prob that tagging occurs in a given year is based on opport_tag_year_prob and a uniform distributed random number AND
      // prob that tagging occurs in a given region is based on opport_tag_prob and a uniform distributed random number (random tagging by area and year)
      
  init_number age_dist_tags
   // determines how tags are distributed across ages
   //==0 regional tags distributed evely across all ages
   //==1 regional tags distributed according to survey proportions at age
   //==2 regional tags distributed according to regional catch proportions at age
   
  init_number tag_fit_ages_switch_OM
   // #DOES NOT WORK WITH TAG MIXING!!!   
   // determines whether OM should use same dynamics as EM when fitting tags by cohort, ie match/mismatch with tag_fit_ages_switch
   //==0, OM maintains age-based dynamics; matches EM if tag_fit_ages_switch==0
   //==1, uses only fully selected F to calculate tag recaps (use tag_age_sel to define age of full selection); matches EM if tag_fit_ages_switch==1
      //OM tag dynamics use only fully selected F (instead of maintaining age-based dynamics then summing across ages)

  init_number sim_tag_mixing_switch
   // determines whether tags have a different F or T in first year of release compared to rest of population (i.e., if there is incomplete mixing)
   //==0 F and T same as rest of pop (complete mixing)
   //==1 F and/or T have are different from rest of population (incomplete mixing);
      // used in combo with sim_tag_mixing_T_switch to determine whether applies to F+T, just F, or just T
      // if sim_tag_mixing_T_switch==0, tag mixing only applies to only to F
      // if sim_tag_mixing_T_switch==1, tag mixing applies to both F and T
   
  init_number sim_tag_mixing_T_switch
   // determines whether tags have a different T in first year of release compared to rest of population (i.e., if there is incomplete mixing)
   //==0 T same as rest of pop (complete mixing)
   //==1 T have are different from rest of population (incomplete mixing)
      // used in combo with sim_tag_mixing_switch
      // if sim_tag_mixing_switch==0, then no incomplete mixing for F or T
      // if sim_tag_mixing_switch==1, then incomplete mixing for both F AND T
      

 //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 //##############################################################################
 //##############################################################################
 //### Inputs for MSY/reference point simulations, set MSY_model_type_switch==0 to do normal F based sims (e.g., dunce cap F) and ignore other inputs for this section ###################################################
 //###########################################################################
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   //############################################################################################################################################################################# 
   //#*** Can use following inputs to either calculate the MSY from the system (by pop/reg) or to explore impact of implementing harvest that does not equal TRUE MSY/uMSY
   //#*** To determine MSY use MSY_model_type_switch==3 AND F_switch==3 AND make sure phase_F>0, this does search algorithm to determine F that maximizes equilibrium system-wide catch
   //#** Options 1 and 2 let you input a TAC or harvest rate, solve for an F that achieves it, and see resulting impact on yield, Bio, etc...
   //#** Can also use to test different approaches for spatial catch allocation (e.g., based on survey biomass, equal) 
   //#** Mainly useful for inputting MSY values from incorrect assumptions regarding pop structure or movement 
   //#** NOTE: should use an appropriately long time frame for sims (e.g., >100 years) to approximate equilibrium 
   //#** NOTE: these methods have not been tested with time-varying movement so equilibrium may not be feasible
   //################################################################################################################################################
   //#** Examples of how the reference point models and mismatches in population structure can be used can be found in Goethel and Berger (2017); https://cdnsciencepub.com/doi/10.1139/cjfas-2016-0290
   //#** Examples of how the reference point models can be used to allocate panmictic TAC to areas can be found in Bosley et al. (2019); https://www.sciencedirect.com/science/article/abs/pii/S0165783619301997?via%3Dihub
   //###################################################################################################################################

  init_number MSY_model_type_switch
   // Changes the type of harvest model (i.e., how F is treated) for simulations
   // DEFAULT SHOULD BE ==0 UNLESS ARE EXPLORING MSY OR ASSOCIATED MISMATCHES IN REFERENCE POINTS
   // Only set >0 if want to fix F at a specific value based on TAC or harvest rate or find F_MSY
   // NOTE: IF >0, Suggest making # years >100 to simulate equilibrium
   //==0 use to simulate catch based on input F (DEFAULT)
   //==1 use input TAC to set F OR calculate TAC from input harvest rate; USE THIS IF calc_TAC_from_uMSY==1; (use if want to input a desired TAC level and use Newton-Raphson search to determine associated F)
   //==2 use input harvest rate (uMSY) to set F; (use if want to input a desired harvest rate level and use Newton-Raphson search to determine associated F)
   //==3 Search for F_MSY (maximize equilibrium total system catch), ENSURE F_Switch==3 AND phase_F>0

  init_number nyrs_quasi_equil
   //number of years over which to average yield to determine 'equilibrium' yield when estimating FMSY (used when MSY_model_type_switch==3)
   
  init_number F_max
   //max F allowed in FMSY search
    
 //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 //##############################################################################
 //##############################################################################
 //### Following used when MSY_model_type_switch>0 and <3 (i.e., for reference point simulations ###################################################
 //###########################################################################
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  init_number parse_TAC
   // method for allocating the TAC to population/region
   // NOTE: IF PARSE input_u directly, when MSY_model_type_switch==2, (instead of basing TAC on u*bio and then allocating TAC; MSY_model_type_switch==1), then sum(u) unlikely to equal input_u because of differences in population sizes (ie, applying less than the full uMSY to each area/region)
   //==0 do not alter the input TAC or harvest rate
   //==1 use observed data source to parse TAC or harvest rate (used when allocating harvest but pop structure unknown; eg, Bosley et al. (2017))

  init_number calc_TAC_from_uMSY
   // determines how TAC is determined when MSY_model_type_switch==1
   //==0 just use input TAC, do not calculate TAC from input harvest rate
   //==1 calculate TAC based on input_u and population biomass; uMSY(input)*biomass_population(j,y) to get a yearly TAC that obtains MSY in equil without crashing the stock
   //==2 calculate TAC based on input_u and regional biomass;  uMSY(input)*biomass_region(j,r,y) to get a yearly TAC that obtains MSY in equil without crashing the stock
      //if there are regions within populations and want to match true MSY, then need to set ==2 otherwise will be based on the total pop bio and TAC will be >MSY
     
  init_number parse_TAC_source
   // if allocating the TAC or harvest rate (parse_TAC==1), then this determines data source to use for parsing to population/region
   //==0 use recruitment index_BM, assume index occurs at tspawn so always have 1 year timelag in using rec index to parse TAC, use equal allocation in first year (no data)
   //==1 use recruitment index_AM, assume index occurs at tspawn so always have 1 year timelag in using rec index to parse TAC, use equal allocation in first year (no data)
   //==2 use survey biomass, if tsurvey==0 can use survey in current year or assume timelag; use equal allocation if year<timelag (no data)
   //==3 use equal apportionment across all fleets in a given region
  
  init_number TAC_survey_parse_timelag_switch
   // if parse_TAC_source==2, determine whether or not to implement a timelag in availability of survey data
   //==0 no timelag, use survey apportionment in current year (if tsurvey==0) or previous year (if tsurvey>0)
   //==1 use timelag, use survey apportionment from y-TAC_survey_parse_timelag, assume equal apportionment of TAC among fleets in first years<timelag
  
  init_number TAC_survey_parse_timelag
   // whole number value timelag in years for use with survey apportionment
   // NOTE: only used if parse_TAC_source==2 AND TAC_survey_parse_timelag_switch==1
  
 //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  //##############################################################################
  //##### End reference point inputs
  //##############################################################################
 //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

///////////////////////////////////////////////////////////////////////////////
//////// ADDITIONAL OM PARAMETERS FROM DAT FILE //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//#############################################################################################################################################
//#/////tagging data parameters
//#############################################################################################################################################

  init_int nyrs_release                                                     // number of years with tag release events
  !! int ny_rel=nyrs_release;
  init_vector yrs_releases(1,ny_rel)                                        // model years with releases
  init_3darray input_ntags(1,np,1,nreg,1,ny_rel)                            // number tag releases per year
  init_vector input_total_tags(1,ny_rel)                                    // total tags releases
  init_vector frac_abund_tagged(1,ny_rel)                                   // proportion of abundance that is tagged in each release year
  init_matrix max_tags(1,np,1,nreg)                                         // maximum tags for uniform distribution if tags_switch==4 or 5
  init_matrix min_tags(1,np,1,nreg)                                         // minimum tags for uniform distribution if tags_switch==4 or 5
  init_number opport_tag_prob                                               // cutoff value defining whether uniform RNG represents no tagging(<tag_prob) or tagging (>tag_prob) in a given AREA
  init_number opport_tag_prob_year                                          // cutoff value defining whether uniform RNG represents no tagging(<tag_prob) or tagging (>tag_prob) in a given YEAR
  init_int max_life_tags                                                   // number of years that tag recaptures will be tallied for after release (assume proportional to longevity of the species)...use this to avoid calculating tag recaptures for all remaining model years after release since # recaptures are often extremely limited after a few years after release
  init_int tag_age_sel                                                      // age used to compute tag dynamics (F and T) when fitting by cohort and OM tag dynamics match EM
  
  init_3darray report_rate_TRUE(1,np,1,ny_rel,1,nreg)                       // tag reporting rate (assume constant for all recaptures within a given release cohort, but can be variable across populations or regions)...could switch to allow variation across fleets instead
                                                                               //reporting rate is assumed to be function of release event and recap location (not a function of recap year...could expand to this, but not priority at momement)

  init_vector F_tag_scalar(1,ny_rel)                                       // scalar applied to F when incomplete tag mixing assumed
  init_vector T_tag_res(1,ny_rel)                                           // movement residence term when incomplete tag mixing assumed

//#############################################################################################################################################
//######################### Movement Parameter Inputs #########################################################################################
//#############################################################################################################################################

  init_vector sigma_T(1,np)                                                 // variance term for random variation in movement
  
// First Post-settlement movement parameters
  init_vector age_first_move(1,np)                                          // age at which first movement occurs
  
// Density-dependent movement parameters
  init_3darray input_Bstar(1,np,1,nreg,1,na)                                // carrying capacity term for density dependent movement, used with move_switch==8
  init_matrix SSB_zero_appor(1,np,1,nreg)                                   // SSB apportionment for density depenent movement, used with move_switch==8 (NOT CURRENTLY IMPLEMENTED)
  init_3darray A(1,np,1,nreg,1,na)                                          // DD movement logistic parameter
  init_3darray DD_residency(1,np,1,nreg,1,na)                               // DD movement residency term

// Natal homing/overlap movement parameters
  init_number return_age                                                    // age of return to natal population, used if move_swith ==6
  init_vector return_probability(1,np)                                      // probability of return to natal population, used if move_swith==6
  init_vector spawn_return_prob(1,np)                                       // probability of performing yearly spawning migration to natal population, used if spawn_return_switch==1

//input residency and full movement arrays
  init_matrix input_residency_larval(1,np,1,nreg)                           // larval residency probability for larval_move_switch==2 or 5
  init_3darray input_residency(1,np,1,nreg,1,na)                            // adult residency probability  for move_switch=2 or 5
  init_5darray input_T(1,np,1,nreg,1,na,1,np,1,nreg)                        // age-based movement array for move_switch==1
  init_5darray input_T_year(1,np,1,nreg,1,ny,1,np,1,nreg)                   // yearly movement array for move_switch==9
  !! int xn=na*ny;
  init_5darray T_Full_Input(1,np,1,nreg,1,xn,1,np,1,nreg)                   // full year+age movement array for move_switch==10; can't input 6darray, so need to compress age*year into single index
  
//#############################################################################################################################################
//######################### Recruitment Parameter Inputs ######################################################################################
//#############################################################################################################################################

  init_vector tspawn(1,np)                                                  // time of spawning in proportion of year (0-1)
  init_vector steep(1,np)                                                   // Beverton-Holt steepness parameter
  init_vector ln_R_ave(1,np)                                                // Average Recruitment or virgin recruitment (R0) for B-H S-R curve
  init_matrix input_Rec_prop(1,np,1,nreg)                                   // recruit apportionment to each region within a population
  init_vector amplitude(1,np)                                               // amplitude of periodic recruitment in % of R_ave, for rec_type==3 
  init_vector freq(1,np)                                                    // frequency of periodic recruitment in years (ie, 10 for peak every 10 years), for rec_type==3 
  init_vector sigma_recruit(1,np)                                           // recruitment deviations variance term, used when recruit_devs_switch==1
  init_vector sigma_rec_prop(1,np)                                          // recruitment apportionment variance term, used when apportionment_type==4

//#############################################################################################################################################
//######################### Fishing Fleet Parameter Inputs ####################################################################################
//#############################################################################################################################################

  init_3darray input_TAC(1,np,1,nreg,1,nf)                                  // TAC to find associated F for, if MSY_model_type_switch==1 and calc_TAC_from_uMSY==0
  init_3darray input_u(1,np,1,nreg,1,nf)                                    // harvest rate to find associated F for, if MSY_model_type_switch==2, OR if MSY_model_type_switch==1 and calc_TAC_from_uMSY>0

//Newton-Raphson parameters
  init_number max_Fnew                                                      // max F for Newton-Raphson search, used when MSY_model_type_switch>0 && MSY_model_type_switch<3
  init_number Fnew_start                                                    // starting F for Newton-Raphson search, used when MSY_model_type_switch>0 && MSY_model_type_switch<3
  init_number NR_iterations                                                 // number of iterations to undergo for Newton-Raphson search, used when MSY_model_type_switch>0 && MSY_model_type_switch<3
  init_number NR_dev                                                        // step size for Newton-Raphson search, used when MSY_model_type_switch>0 && MSY_model_type_switch<3
  
  init_3darray input_F_MSY(1,np,1,nreg,1,nf)                                // specify F_MSY directly to simulate MSY dynamics, used with F_switch==2
  init_3darray input_F(1,np,1,nreg,1,nf)                                    // specify F for simulations, used with F_switch==1,7,8
  init_number input_F_const                                                 // specify F for simulations, used with F_switch==4,5,6
  init_3darray dunce_F(1,np,1,nreg,1,4)                                     // dunce cap F parameters (fishery start year, min F start, max F, min F end), used with F_switch==9
  init_3darray F_rho(1,np,1,nreg,1,nf)                                      // degree of autocorrelation (0-1) for  random walk in F, used with F_switch==8; NOT YET TESTED
  init_3darray sigma_F(1,np,1,nreg,1,nf)                                    // F deviations variance term, used when F_switch==7,8,9

  init_3darray sel_beta1(1,np,1,nreg,1,nf)                                  // fishery selectivity slope parameter for logistic selectivity/double logistic, used when select_switch==1,2
  init_3darray sel_beta2(1,np,1,nreg,1,nf)                                  // fishery selectivity inflection parameter for logistic selectivity/double logistic, used when select_switch==1,2
  init_3darray sel_beta3(1,np,1,nreg,1,nf)                                  // fishery selectivity descending slope parameter for double selectivity, used when select_switch==2
  init_3darray sel_beta4(1,np,1,nreg,1,nf)                                  // fishery selectivity descending inflection parameter for double logistic selectivity, used when select_switch==2
  init_4darray input_selectivity(1,np,1,nreg,1,na,1,nf)                     // fishery selectivity by pop/region/age/fleet, used when select_switch==0

//#############################################################################################################################################
//######################### Survey Fleet Parameter Inputs ######################################################################################
//#############################################################################################################################################

//#******************************************************************************************************************************************************************
//# NOTE: EM assumes that survey parameters (q+selectivity) are constant across regions, therefore to avoid bias make selectivity parameters constant across REGIONS
//#******************************************************************************************************************************************************************

  init_matrix tsurvey(1,np,1,nreg)                                          // time of survey in proportion of year (0,1)
  init_3darray q_survey(1,np,1,nreg,1,nfs)                                  // survey catchability scalar
  init_3darray sel_beta1_survey(1,np,1,nreg,1,nfs)                          // survey selectivity slope parameter for logistic selectivity/double logistic, used when select_switch_survey==1,2
  init_3darray sel_beta2_survey(1,np,1,nreg,1,nfs)                          // survey selectivity inflection parameter for logistic selectivity/double logistic, used when select_switch_survey==1,2
  init_3darray sel_beta3_survey(1,np,1,nreg,1,nfs)                          // survey selectivity descending slope parameter for double selectivity, used when select_switch_survey==2
  init_3darray sel_beta4_survey(1,np,1,nreg,1,nfs)                          // survey selectivity descending inflection parameter for double logistic selectivity, used when select_switch_survey==2
  init_4darray input_survey_selectivity(1,np,1,nreg,1,na,1,nfs)             // survey selectivity by pop/region/age/fleet, used when select_switch_survey==0

//#############################################################################################################################################
//######################### Biological Parameter Inputs ########################################################################################
//#############################################################################################################################################

//#**************************************************************************************************************************************************************
//# NOTE: FOR NATAL HOMING THERE IS NO ACCOUNTING OF REGIONAL DIFFERENCES IN VITAL RATES ACROSS REGIONS WITHIN A POPULATION
//#       BECAUSE IT IS ASSUMED THAT GENETICS DEFINE VITAL RATES, WEIGHT, MATURITY, FECUNDITY ARE ALL TREATED AS CONSTANT ACROSS REGIONS WITHIN A POPULATION
//#       FOR OM_structure==3, ENSURE THAT VITAL RATES ARE INPUT AS CONSTANT ACROSS REGIONS, BECAUSE NATAL REGION WILL BE IGNORED (e.g, IN SSB CALCS) 
//#***************************************************************************************************************************************************************

  init_3darray input_weight(1,np,1,nreg,1,na)                               // weight-at-age
  init_3darray input_catch_weight(1,np,1,nreg,1,na)                         // catch weight-at-age
  init_3darray fecundity(1,np,1,nreg,1,na)                                  // fecundity-at-age, used for SSB_type==1 (also uses maturity)
  init_3darray maturity(1,np,1,nreg,1,na)                                   // maturity-at-age, used for SSB_type==1 AND 2 
  init_matrix prop_fem(1,np,1,nreg)                                         // proportion of population assumed to be female for SSB calcs (typically use 0.5)

//#############################################################################################################################################
//######################### Other Demographic Parameter Inputs ################################################################################
//#############################################################################################################################################

  init_matrix input_M_TRUE(1,np,1,na)                                       // natural mortality-at-age
  init_matrix equil_ssb_apport(1,np,1,nreg)                                 // weights given to each population for the weighted average of maturity- and weight-at-age for spawner-per-recruit (SPR) calcs when maturity_switch_equil==1
  init_3darray frac_natal(1,np,1,np,1,nreg)                                 // distribution of all natal populations across populations and regions in first year, used with init_abund_switch==1
  init_4darray input_init_abund(1,np,1,np,1,nreg,1,na)                      // distribution of all natal populations across populations and regions in first year, used with init_abund_switch==0

//#############################################################################################################################################
//######################### Observation Error Parameter Inputs ###############################################################################
//#############################################################################################################################################

  init_matrix rec_index_sigma(1,np,1,nreg)                                  // variance for recruitment index lognormal observation error
  init_3darray sigma_survey(1,np,1,nreg,1,nfs)                              // variance for survey biomass index lognormal observation error
  init_4darray sigma_survey_overlap(1,np,1,np,1,nreg,1,nfs)                 // variance for survey biomass index when natal_homing_switch==1 (natal homing/overlap) lognormal observation error
  init_3darray sigma_catch(1,np,1,nreg,1,nf)                                // variance for catch data lognormal observation error
  init_4darray sigma_catch_overlap(1,np,1,np,1,nreg,1,nf)                   // variance for catch data when natal_homing_switch==1 (natal homing/overlap) lognormal observation error
  init_3darray SIM_ncatch(1,np,1,nreg,1,nf)                                 // multinomial effective sample size for catch age composition data
  init_4darray SIM_ncatch_overlap(1,np,1,np,1,nreg,1,nf)                    // multinomial effective sample size for catch age composition data when use_stock_comp_info_survey==1 (natal homing/overlap)
  init_3darray SIM_nsurvey(1,np,1,nreg,1,nfs)                               // multinomial effective sample size for survey age composition data
  init_4darray SIM_nsurvey_overlap(1,np,1,np,1,nreg,1,nfs)                  // multinomial effective sample size for survey age composition data when use_stock_comp_info_survey==1 (natal homing/overlap)
  init_number SIM_ntag                                                      // multinomial effective sample size for tagging data

//****************************************************************************************************************************************************************************************************************
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//################################################################################################################################################################################################################
//###############################################################################################################################################################################################################
//###############################################################################################################################################################################################################
//# Estimation Model Inputs
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//# These inputs define the structure of the EM, which can either match or mismatch the OM leading to misspecification along various axes (e.g., spatial structure or otherwise)
//# Those values not defined for the EM explicitly here likely match the OM assumptions
//# The EM has some options/functionality not explicitly incorporated in the OM (e.g., areas-as-fleets structure)
//# Conversely, the OM has a lot of structure that cannot be accounted for in the EM (as noted in the README at the top of this file)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//####################################################################################################################################################################################################################
//###############################################################################################################################################################################################################
//#################################################################################################################################################################################################################
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//******************************************************************************************************************************************************************************************************************



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// ESTIMATION MODEL (EM) STRUCTURE+SWITCHES INPUTS ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_number EM_structure //
   // define population structure in EM to determine level of data aggregation that needs to occur
   // mainly used when OM_structure>0 (not panmictic), but EM is panmictic (i.e., aggregating data to a single area)
   //==(-1) the EM is areas-as-fleets (AAF) (aka, fleets-as-areas, FAA), models each region in OM as a different fleet
   //==0 the EM is panmictic (1 area well mixed)
   //==1 the EM is metamictic (spatial heterogeneity, one population with multiple regions)
   //==2 the EM is metapop (multiple subpopulations with unique demographics/vital rates, but reproductive mixing such that fish immediately mix when change location and take on characteristics of new population)
   //==3 the EM is natal homing/overlap/home fidelity (vital rates are based on natal population; all spawners can return to spawn or only fraction...based on input overlap parameters)

  init_number natal_homing_switch_EM
   // determines how SSB is tallied (i.e., follow natal homing assumption that fish must return to natal population to add to SSB), mainly check on OM_structure setting
   //==0 no natal homing, use if EM_structure!=3 (SSB is sum of SSB in population regardless of natal origin; weight/mat/fecund/ are based on current population not natal population)
   //==1 do natal homing, use if EM_structure==3 (a fish only adds to SSB if it is in its natal population at spawning time; weight/mat/fecund/ are based on natal population)
      // natal homing assumes genetic based life history and contribution to SSB (i.e., natal homing and no demographic mixing), natal_homing_switch_EM==0 assumes demographic mixing (e.g. metapopulations where life history is more location based)

  init_number spawn_return_switch_EM
   // determines whether an instantaneous return spawning migration to natal population occurs IF natal_homing_switch_EM==1
   //==0 if natal_homing_switch_EM==1, then ONLY fish that are currently in natal population at time of spawning add to SSB
   //==1 if natal_homing_switch_EM==1, then a fraction of fish return to natal population to spawn (instantaneous migration to natal population and back at time of spawning) based on spawn_return_prob; weight/mat/fecund/ are based on natal population

  init_number npops_EM
   // number of populations in the EM
   // NOTE: if EM_structure==(-1) OR 0, npops should be ==1
   // NOTE: if EM_structure==1 (spatial heterogeneity), npops typically ==1 (but could have multiple populations with no movement among pops; typically would just use metapop structure though)
   // NOTE: if EM_structure>1 (natal homing or metapopulation), npops should >1, otherwise have either panmictic or spatial heterogeneity

  !! int np_em=npops_EM;
   //can't be type init so need to establish as integer
   
  init_ivector nregions_EM(1,np_em)
   // number of regions per population in EM
   // SEE README IN DAT FOR POTENTIAL MISSPECIFICATION OF EM IF DIFFERING NUMBER OF REGIONS IN EACH POPULATION
   // NOTE: if EM_structure==(-1) (AAF) OR 0 (panmictic), nregs should ==1, otherwise have spatial heterogeneity
   // NOTE: if EM_structure==1 (spatial heterogeneity), nregs should be >1, otherwise have panmictic
   // NOTE: if EM_structure>1 (natal homing or metapopulation), nregs can be =>1

  init_ivector nfleets_EM(1,np_em)
   // number of fishing fleets in each region for each population in EM
   // SEE README IN DAT FOR POTENTIAL MISSPECIFICATION OF EM IF DIFFERING NUMBER OF FLEETS IN EACH POPULATION

  init_ivector nfleets_survey_EM(1,np_em)
   // number of survey fleets in each region for each population in EM
   // SEE README IN DAT FOR POTENTIAL MISSPECIFICATION OF EM IF DIFFERING NUMBER OF FLEETS IN EACH POPULATION

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// INDICES for ragged arrays /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //can't be type init so need to establish as integer
  !! ivector nreg_em=nregions_EM;
  !! ivector nf_em=nfleets_EM;
  !! ivector nfs_em=nfleets_survey_EM;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////// EM Model SWITCHES /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_number diagnostics_switch
   // determines whether want to do a 'diagnostic' run without observation error (i.e., fit 'true'/expected value of each data source)
   // essentially just reports the expected value of each data source as the 'observed' value when set to TRUE(==1)
   //==0 fit data with observation error (DEFAULT)
   //==1 fit expected value (true) data (no observation error)

  init_number move_switch_EM
   // Sets the EM parametrization of adult movement pattern (sets age class>0 movements)
   // NOTE: MUST CHANGE THE PHASES OF THE MOVEMENT PARAMETERS (T) TO DETERMINE TYPE OF MOVEMENT (constant, Time varying, age-varying, both time and age varying); SEE BELOW FOR DESCRIPTION
   // NOTE: IF PANMICTIC, NO MOVE, OR INPUT T, THEN  MAKE SURE ALL T PAR PHASES ARE SET TO NEGATIVE NUMBER
   //==(-2) no movement because panmictic, SET T PAR PHASES==(-1); residency set equal to 1.0 in EM (to ensure mathematical consistency in calcs)
   //==(-1) input age based movement using input_T_EM, SET T PAR PHASES==(-1); used when want fixed movement that differs from true movement
   //==0    no movement, SET T PAR PHASES==(-1); residency set equal to 1.0 in EM (to ensure mathematical consistency in calcs)
   //==1    use true movement from OM, SET T PAR PHASES==(-1)
   //==2    estimate movement, MUST SET DESIRED T PARAMETER PHASE>0; see DAT file for full list of movement estimation options by phase input
   //==6    natal return (for natal_homing_switch==1 ONLY), no movement of fish until return_age when a certain fraction, return_probability, of fish make return migration to natal population (eg, an ontogenetic migration); all fish remain in given population for remainder of lifepsan after return_age
   //==7    larvae stay in the population that they settle in (i.e., for natal homing/overlap, do not return to natal population)
            // if adult movement==0 for natal homing would return to natal population because natal residency is 100% and use natal movement rates (not current population movement rates like with metapopulation/random movement)
   //################################################################################################################################################################
   //#** Examples of how different movement parametrizations can be used can be found in Goethel et al. (2020); https://onlinelibrary.wiley.com/doi/full/10.1111/faf.12510
   //###############################################################################################################################################################
   
 //################################################################################################################################################################################################################################################################################################################################################################################################
 //# NOTE: following phase parameters for movement estimation used when move_switch_EM==2
 //###################################################################################################################################################################################################################################################################################################################################################################################
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_int phase_T_CNST                                                     // estimate a time+age-invariant movement rate
  init_int phase_T_CNST_AGE                                                 // estimate an age-varying and time-invariant movement rate
  init_int phase_T_CNST_AGE_no_AG1                                          // estimate an age-varying and time-invariant movement rate, but fix age-1 movement at 100% residency; may be useful when little data to estimate age-1 movement
  init_int phase_T_YR                                                       // estimate a yearly and age-invariant movement rate 
  init_int phase_T_YR_ALT_FREQ                                              // estimate a yearly and age-invariant movement rate, but with time blocks (interval width) based on T_est_freq (number of year parameters is nyrs/T_est_freq, rounded up)
  init_int phase_T_YR_AGE                                                   // estimate a yearly and age-varying movement rate (i.e., a movement parameter for every year AND every age); this is likely to cause model instability given high number of freely estimated parameters
  init_int phase_T_YR_AGE_no_AG1                                            // estimate a yearly and age-varying movement rate, but fix age-1 movement at 100% residency; may be useful when little data to estimate age-1 movement
  init_int phase_T_YR_AGE_ALT_FREQ                                          // estimate an age+time-varying movement rate, but with age+time blocks (interval width) based on T_est_age_freq and T_est_freq
  init_int phase_T_YR_AGE_ALT_FREQ_no_AG1                                   // estimate an age+time-varying movement rate, but with age+time blocks (interval width) based on T_est_age_freq and T_est_freq, AND fix age-1 movement at 100% residency
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_number select_switch_EM
   // determine how fishery selectivity is estimated
   //#==(-1) fix selectivity at TRUE value from OM
   //==0    input time-invariant selectivity based on input_selectivity_EM
   //==1    estimated logistic selectivity based on estimated sel_beta1 and sel_beta2 parameters; ph_sel_log must be >0 else fixed at starting values for each parameter
   //==2    estimated double logistic selectivity based on estimated sel_beta1, sel_beta2, sel_beta3, sel_beta4 parameters; ph_sel_log must be >0 AND ph_sel_dubl must be >0 else fixed at starting values for each parameter
            // NOTE: NOT TESTED YET

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_int ph_sel_log                                                       // phase for fishery logistic OR double logistic selectivity estimation, used with select_switch_EM==1 OR 2
  init_int ph_sel_dubl                                                      // phase for fishery double logistic selectivity estimation, used with select_switch_EM==2
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_number select_switch_survey_EM
   // determine how survey selectivity is estimated
   //==(-1) fix selectivity at TRUE value from OM
   //==0    input time-invariant selectivity based on input_survey_selectivity_EM
   //==1    estimated logistic selectivity based on estimated sel_beta1 and sel_beta2 parameters; ph_sel_log_surv must be >0 else fixed at starting values for each parameter
   //==2    estimated double logistic selectivity based on estimated sel_beta1, sel_beta2, sel_beta3, sel_beta4 parameters; ph_sel_log_surv must be >0 AND ph_sel_dubl_surv must be >0 else fixed at starting values for each parameter
            // NOTE: NOT TESTED YET

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_int ph_q                                                             // phase for survey catchability scalar estimation, used with select_switch_survey_EM==1 OR 2
  init_int ph_sel_log_surv                                                  // phase for survey logistic OR double logistic selectivity estimation, used with select_switch_survey_EM==1 OR 2
  init_int ph_sel_dubl_surv                                                 // phase for survey double logistic selectivity estimation, used with select_switch_survey_EM==2
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_number F_switch_EM
   // determine how fishing mortality is estimated
   //==0 fix F at TRUE value from OM
   //==1 estimate yearly F by population, region, fleet; ph_F must be >0 else F is fixed at input value in EM
   //==2 estimate F as a random walk; NOTE: NOT IMPLEMENTED DO NOT USE

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_int ph_F                                                             // phase for fishing mortality estimationn, used with F_switch_EM==1
  init_int ph_F_rho                                                         // (NOT IMPLEMENTED) phase for fishing mortality random walk estimationn, used with F_switch_EM==2
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_number M_switch_EM
   // determine how natural mortality is estimated
   // NOTE: all M approaches assume M is spatially-invariant within a population (i.e., constant across regions within a population)
   //==(-1) fix M at input_M_EM (can differ from TRUE M)
   //==0    fix M at TRUE value from OM
   //==1    estimate constant M (const across pop and age)
   //==2    estimate population-based M (const across ages)
   //==3    estimate age-based M (const across pop)
   //==4    estimate age- and population-varying M

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_int ph_M_CNST                                                        // phase for estimating age+time+space-invariant constant M (const across pop, age, time), used when M_switch_EM==1
  init_int ph_M_pop_CNST                                                    // phase for estimating population-specific, age+time-invariant M (const across ages, time), used when M_switch_EM==2
  init_int ph_M_age_CNST                                                    // phase for estimating age-varying, pop+time-invariant M (const across pop+time), used when M_switch_EM==3
  init_int ph_M_pop_age                                                     // phase for estimating age+population-varying, time-invariant M (const across time), used when M_switch_EM==4
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_number init_abund_switch_EM
   // determine whether initial abundance at age in first year is estimated, fixed, or based on exponential decay from Rave
   //==(-2) fix initial abundance at init_abund_EM by pop, reg, age (can differ from TRUE distribution in OM)
   //==(-1) fix initial abundance by pop, reg, age at TRUE values from OM
   //==0    estimate initial abundance at age by population, where R_ave is the age-1 abundance in year 1 (to avoid overparametrization of recruitment); ph_init_abund_no_ag1 MUST BE>0 else fixed at starting values for each parameter
   //==1    estimate initial abundance at age by population where all ages estimated; ph_init_abund MUST BE>0 else fixed at starting values for each parameter
   //==2    initial abundance at age is based on exponential decay with age starting from R_ave (R0) at age-1 and assuming no fishing mortality (i.e., decay based on M only)
      // NOTE: all estimation or derivation (all options except -2 or -1) of init_abund requires distributing init_abund to regions using the est_dist_init_abund_EM switch

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_int ph_init_abund                                                    // phase for initial abundance estimation of all ages, used when init_abund_switch_EM==1
  init_int ph_init_abund_no_ag1                                             // phase for initial abundance estimation of all ages except age-1, used when init_abund_switch_EM==0
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_number est_dist_init_abund_EM
   // determine how initial abundance at age from a natal population is distributed across all populations and regions (natal_homing_switch_EM==1) or across regions within a population (natal_homing_switch_EM==0)
   //==(-2) fix the initial spatial distribution at input_dist_init_abund_EM (can differ from TRUE distribution in OM) (CONSTANT ACROSS AGES)
   //==(-1) equally distribute across regions in a population (no fish start outside natal population) (CONSTANT ACROSS AGES)
   //==0    fix the initial spatial distribution at TRUE values from OM averaged across ages (CONSTANT ACROSS AGES)
   //==1    fix the initial spatial distribution at TRUE values from OM by age (VARIES ACROSS AGES)
   //==2    estimate the initial spatial distribution of init abund (constant across ages), ph_non_natal_init OR ph_reg_init MUST BE >0 (CONSTANT ACROSS AGES)
   //==3    estimate the spatial distribution of init abundance by age, need to make ph_non_natal_age_init OR ph_reg_age_init non-negative depending on spatial structure used (VARIES BY AGE)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_int ph_reg_init                                                      // phase for estimation of initial distribution of abundance of each population across all regions within that population averaged across ages (non-natal homing models), used when est_dist_init_abund_EM==2 AND natal_homing_switch_EM==0 (CONSTANT BY AGE)
  init_int ph_reg_age_init                                                  // phase for estimation of initial distribution of abundance of each population across all regions within that population by age (non-natal homing models), used when est_dist_init_abund_EM==3 AND natal_homing_switch_EM==0 (VARIES BY AGE)
  init_int ph_non_natal_init                                                // phase for estimation of initial distribution of abundance of each natal population across all other populations and regions averaged across ages (including regions within natal population; for natal homing models), used when est_dist_init_abund_EM==2 AND natal_homing_switch_EM==1  (CONSTANT BY AGE)
  init_int ph_non_natal_age_init                                            // phase for estimation of initial distribution of abundance of each natal population across all other populations and regions by age (including regions within natal population; for natal homing models), used when est_dist_init_abund_EM==3 AND natal_homing_switch_EM==1  (VARIES BY AGE)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_number maturity_switch_equil_EM
   // MORE WORK IS NEEDED TO REFINE THESE CALCULATIONS AND DEAL WITH SPATIAL REFERENCE POITNS!
   // SSB0 must be calculated to determine stock-recruit function
   // Use equilibrium SPR calcs to get SSB0, but to do so requires vital rates (maturity, weight), which are typically constant across a population
   // With multiple regions within a pop each with different vitals, must make assumption regarding the proportional contribution of each region's demograhics to equil SSB
   //==0  equal by region, assume equal (average) contributions to SSB0 by each region
   //==1 weighted average, use input equil_ssb_apportion to determine proportional contribution to equil vital rates by region

  init_number SSB_type_EM
   // units of spawning stock biomass
   //==1 fecundity based SSB
   //==2 weight based SSB
  
  init_number Rec_type_EM
   // form of the stock-recruit relationship
   //==1 stock-recruit relationship assumes an average value based on R_ave
   //==2 Beverton-Holt population-recruit functions based on population-specific estimated steepness, R0 (R_ave)
         // NOTE: SRR DOES NOT TAKE INTO ACCOUNT SPATIAL DYNAMICS (IE, MOVEMENT AMONG POPULATIONS OR REGIONS); SEE README TOPIC AT TOP OF DAT FILE
   //==3 environmental recruitment - sine fucntion based on amplitude and frequency
   //==4 (NOT IMPLEMENTED...USE OPTION 2 for BH) Beverton-Holt population-recruit functions based on population-specific fixed (at true value) steepness, estimated R0 (R_ave)
   //==5 (NOT IMPLEMENTED...USE OPTION 2 for BH) Beverton-Holt stock-recruit functions based on stock-specific fixed (at true value) steepness and R0 (R_ave) (i.e., both SRR parameters fixed at true values)
   //==6 (NOT IMPLEMENTED...USE OPTION 1 for AVE R) stock-recruit relationship assumes an average value based on fixed (at true value) R_ave
   
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_int ph_lmr                                                           // phase for stock-recruit relationship log_mean_recruitment (LMR; or R0, virgin recruitment) estimation, used when Rec_type_EM==1 OR 2 OR 4
  init_int ph_steep                                                         // phase for stock-recruit relationship steepness estimation, used when Rec_type_EM==2
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_number recruit_devs_switch_EM
   // determine whether recruitment deviations are estimated
   //==(-1) fix recruit deviations at TRUE value from OM; make sure to set ph_rec==(-1)
   //==0    use stock-recruit relationship directly (no deviations); make sure to set ph_rec==(-1)
   //==1    estimate recruit deviations assuming lognormal variance around SR curve based on input sigma_recruit_EM; make sure to set ph_rec >0

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_int ph_rec                                                           // phase for stock-recruit relationship deviations estimation, used when recruit_devs_switch_EM==1
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_number recruit_randwalk_switch_EM
   //==0 no random walk recruitment deviations (DEFAULT)
   //==1 have random walk lognormal recruitment deviations (requires recruit_devs_switch_EM==1)....HAS NOT BEEN FULLY IMPLEMENTED OR TESTED

  init_number apportionment_type_EM
   // determines how recruits are apportioned to regions within a population
   // because stock-recruit relationships are assumed only at the population level, if have more than 1 region in a population, must make assumption about how to apportion/assign recruits to each region
   // typically used with EM_structure==1, but also applies if nreg>1 in any pop for EM_structure==2 AND 3
   //==(-2) fix at TRUE values from OM
   //==(-1) fix with no recruitment apportionment to regions within a stock (each region within a stock gets full amount of recruits from SR curve); WOULD NOT SUGGEST USING Because leads to sum(recruits)>SRR
   //==0    fix with apportionment to each region based on relative SSB in region compared to pop SSB
   //==1    fix with yearly apportionment based on input_rec_prop_EM
   //==2    fix with equal apportionment to each region within a population
   //==3    estimate time-invariant apportionment; ph_rec_app_CNST MUST BE >0
   //==4    estimate time-varying apportionment; ph_rec_app_YR MUST BE >0
            // NOTE: THIS HAS NOT BEEN THOROUGHLY TESTED, MAY NOT BE IMPLEMENTED CORRECTLY

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_int ph_rec_app_CNST                                                  // phase for time-invariant recruitment apportionment (to region within a population) estimation, used when apportionment_type_EM==3
  init_int ph_rec_app_YR                                                    // phase for time-varying recruitment apportionment (to region within a population) estimation, used when apportionment_type_EM==4
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_number use_stock_comp_info_survey_EM 
   // Determines whether to fit survey age composition by natal area for maximum likelihood calcs
   //==0 fit OBS survey age comps by pop/area (summed across natal population); assumes no stock compisition data is available; typical for most situations when no genetic/otolith analysis is available
   //==1 fit OBS survey age comps by natal population within each area; assumes stock composition data is available
         // NOTE: NOT YET IMPLEMENTED IN MLE CALCS
         
  init_number use_stock_comp_info_catch_EM
   // Determines whether to fit catch age composition by natal area for maximum likelihood calcs
   //==0 fit OBS catch age comps by pop/area (summed across natal population); assumes no stock compisition data is available; typical for most situations when no genetic/otolith analysis is available
   //==1 fit OBS catch age comps by natal population within each area; assumes stock composition data is available
         // NOTE: NOT YET IMPLEMENTED IN MLE CALCS

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////// Tagging Data SWITCHES //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_int do_tag_EM
   // Turn off tagging calcs in EM to save computation time if not going to use tagging data
   //==0 do not calculate tagging data
   //==1 calculate tagging data

  init_number tag_fit_ages_switch
   // determines whether tags are fit by age cohorts or by region-only cohorts in the ESTIMATION MODEL
   // the latter is for situations where age is assumed to be unknown for tags and so assume no tag age dynamics in the EM and that all fish fully selected
   // the tag_fit_ages_switch_OM switch determines whether tags are simulated with matching or mismatching age dynamics as EM
   //==0, EM fits tags assuming age-based cohorts
   //==1, EM fits tags assuming region-based cohorts
      // if tag_fit_ages_switch_OM==0, OM assumes normal tagging age-based dynamics, but EM ignores age structure in tags and assumes all tagged fish are fully selected
      // therefore creates inherent process error in tag dynamics because EM uses T, F, M dynamics of fully selected age (error magnified if OM or EM has age-based M or T)

  init_number est_tag_mixing_switch
   // determines whether EM actually estimates different F or T for tagged fish compared to rest of population in first year of release (i.e., estimate F and/or T to account for  incomplete mixing)
   // leads to match/mismatch between EM and OM depending on sim_tag_mixing_switch AND sim_tag_mixing_T_switch
   // if sim_tag_mixing_T_switch==1 then at least F is simulated assuming incomplete mixing of tagged fish
   // NOTE: EM only setup to est an F scalar between 0 and 1 (i.e., F_tag can only be less than F for incomplete mixing)
   // NOTE: EM not setup to do tag mixing and fit by region based cohorts (i.e., can't use with tag_fit_ages==1)
   //==0 F and T same as rest of pop (assume complete mixing)
      // if sim_tag_mixing_T_switch==0, then this assumes match between OM and EM tag mixing assumptions (no incomplete mixing in either)
      // if sim_tag_mixing_T_switch==1, then this assumes mismatch between OM ane EM in tag mixing assumptions (ie, OM assumes incomplete mixing in at least F, EM assumes complete mixing)
   //==1 F is estimated different from rest of population (incomplete mixing F only) in first year of release
   //==2 T is estimated different from rest of population (incomplete mixing T only) in first year of release
   //==3 F AND T are estimated different from rest of population (incomplete mixing F AND T) in first year of release
   
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_number ph_T_tag                                                      // phase for movement tag mixing estimation, used with est_tag_mixing_switch==2 OR 3
  init_number ph_F_tag                                                      // phase for fishing mortality tag mixing estimation, used with est_tag_mixing_switch==1 OR 3
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  init_int do_tag_mult
   // determines assumed MLE distribution for fitting tagging data
   //==0 assume negative binomial distribution for tagging data
     // NOTE: NOT YET IMPLEMENTED/TESTED
   //==1 assume multinomial distribution (same as OM) (DEFAULT)

  init_number report_rate_switch_EM
   // determines how reporting rate is estimated for tagging data
   // NOTE: if tagging data is not fit in the model make sure all reporting rate phases are negative (not estimated)
   //==(-1) fix based on input_report_rate_EM (reporting rate may differ from TRUE report rate)
   //==0    fix at TRUE values from OM
   //==1    estimate population+region-varying, time-invariant reporting rate; ph_rep_rate_CNST MUST BE >0
   //==2    estimate time-, pop-, reg-varying reporting rate; ph_rep_rate_YR MUST BE >0
   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_int phase_rep_rate_CNST                                              // phase for reporting rate estimation by population+region, used with report_rate_switch_EM==1
  init_int phase_rep_rate_YR                                                // phase for reporting rate estimation by year, used with report_rate_switch_EM==2
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// ADDITIONAL EM PARAMETERS FROM DAT FILE ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#############################################################################################################################################
//#  tagging data parameters
//#############################################################################################################################################

  init_number lb_B                                                          // lower bound on reporting rate estimate in logit transform space, initial estimate is B_start
  init_number ub_B                                                          // upper bound on reporting rate estimate in logit transform space, initial estimate is B_star
  init_number B_start                                                       // initial estimate for reporting rate in logit transform space
  init_number lb_scalar_T                                                   // lower bound on tag mixing movement rate scalar estimate in logit transform space, initial estimate is T_scalar_start
  init_number ub_scalar_T                                                   // upper bound on tag mixing movement rate scalar estimate in logit transform space, initial estimate is T_scalar_start
  init_number scalar_T_start                                                // initial estimate for tag mixing movement rate scalar in logit transform space
  init_number lb_scalar_F                                                   // lower bound on tag mixing fishing mortality scalar estimate in logit transform space, initial estimate is F_scalar_start
  init_number ub_scalar_F                                                   // upper bound on tag mixing fishing mortality scalar estimate in logit transform space, initial estimate is F_scalar_start
  init_number scalar_F_start                                                // initial estimate for tag mixing fishing mortality rate scalar in logit transform space
  init_matrix input_report_rate_EM(1,np_em,1,nreg_em)                       // input tag reporting rate for EM which can differ from true reporting rate, not time-varying; used with report_rate_switch_EM==(-1)

//#############################################################################################################################################
//######################### Movement Parameter Inputs #########################################################################################
//#############################################################################################################################################

  init_vector return_probability_EM(1,np_em)                                // probability of return to natal population, used if move_swith_EM==6
  init_vector spawn_return_prob_EM(1,np_em)                                 // probability of return to natal population for spawning, used if spawn_return_switch_EM==1
  init_int T_est_freq                                                       // defines interval width (ie, yearly time blocks) for movement estimation (number of year movement parameters is nyrs/T_est_freq, rounded up)
  init_int T_est_age_freq                                                   // defines interval width (ie, age time blocks) for movement estimation (number of age parameters is nages/T_age_freq, rounded up)
  init_int juv_age                                                          // age at/below which all Ts are set to the estimated movement rate for age-1 in that timeblock
  init_number lb_T                                                          // lower bound on movement rate estimate in logit transform space, initial estimate is T_start
  init_number ub_T                                                          // upper bound on movement rate estimate in logit transform space, initial estimate is T_start
  init_number T_start                                                       // initial estimate for movement parameters in logit transform space
  init_5darray input_T_EM(1,np_em,1,nreg_em,1,na,1,np_em,1,nreg_em)         // age based movement for EM, used when move_switch_EM==(-1); for fixed movement that differs from true movement

//#############################################################################################################################################
//######################### Recruitment Parameter Inputs ######################################################################################
//#############################################################################################################################################

  init_vector tspawn_EM(1,np_em)                                            // time of spawning for EM in proportion of year (0-1); 0 is Jan 1st  
  init_number lb_steep                                                      // lower bound on steepness parameter estimate, initial estimate is steep_start 
  init_number ub_steep                                                      // upper bound on steepness parameter estimate, initial estimate is steep_start 
  init_number steep_start                                                   // initial estimate for steepness parameter
  init_number lb_R_ave                                                      // lower bound on LMR (R0) parameter estimate in log-space, initial estimate is Rave_start
  init_number ub_R_ave                                                      // upper bound on LMR (R0) parameter estimate in log-space, initial estimate is Rave_start
  init_number Rave_start                                                    // initial estimate for LMR (R0) parameter in log-space 
  init_number lb_rec_devs                                                   // lower bound on recruitment deviations estimates in log-space, initial estimate is Rdevs_start
  init_number ub_rec_devs                                                   // upper bound on recruitment deviations estimates in log-space, initial estimate is Rdevs_start 
  init_number Rdevs_start                                                   // initial estimate for recruitment deviations parameters in log-space
  init_number lb_rec_app                                                    // lower bound on recruitment apportionment estimates in logit transform space, initial estimate is Rapp_start 
  init_number ub_rec_app                                                    // upper bound on recruitment apportionment estimates in logit transform space, initial estimate is Rapp_start
  init_number Rapp_start                                                    // initial estimate for recruitment apportionment parameters in logit transform space 
  init_3darray input_rec_prop_EM(1,np_em,1,nreg_em,1,ny)                    // recruit apportionment for EM to each region within a population by year, used with apportionment_type_EM==1
  init_vector sigma_recruit_EM(1,np_em)                                     // recruitment variance term defining strength of recruitment deviations when recruit_devs_switch_EM==1

//#############################################################################################################################################
//######################### Fishing Fleet Parameter Inputs ####################################################################################
//#############################################################################################################################################

  init_number lb_F                                                          // lower bound on fishing mortality parameter estimate, initial estimate is F_start
  init_number ub_F                                                          // upper bound on fishing mortality parameter estimate, initial estimate is F_start
  init_number F_start                                                       // initial estimate for fishing mortality parameter in log-space
  init_number lb_F_rho                                                      // lower bound on fishing mortality random walk parameter estimate, initial estimate is Frho_start
  init_number ub_F_rho                                                      // upper bound on fishing mortality random walk parameter estimate, initial estimate is Frho_start
  init_number Frho_start                                                    // initial estimate for fishing mortality random walk parameter in log-space
  init_number lb_sel_beta1                                                  // lower bound on fishery selectivity beta1 parameter estimate, initial estimate is sel_beta1_start
  init_number ub_sel_beta1                                                  // upper bound on fishery selectivity beta1 parameter estimate, initial estimate is sel_beta1_start
  init_number sel_beta1_start                                               // initial estimate for fishery selectivity beta1 parameter in log-space
  init_number lb_sel_beta2                                                  // lower bound on fishery selectivity beta2 parameter estimate, initial estimate is sel_beta2_start
  init_number ub_sel_beta2                                                  // upper bound on fishery selectivity beta2 parameter estimate, initial estimate is sel_beta2_start
  init_number sel_beta2_start                                               // initial estimate for fishery selectivity beta2 parameter in log-space
  init_number lb_sel_beta3                                                  // lower bound on fishery selectivity beta3 parameter estimate, initial estimate is sel_beta3_start
  init_number ub_sel_beta3                                                  // upper bound on fishery selectivity beta3 parameter estimate, initial estimate is sel_beta3_start
  init_number sel_beta3_start                                               // initial estimate for fishery selectivity beta3 parameter in log-space
  init_number lb_sel_beta4                                                  // lower bound on fishery selectivity beta4 parameter estimate, initial estimate is sel_beta4_start
  init_number ub_sel_beta4                                                  // upper bound on fishery selectivity beta4 parameter estimate, initial estimate is sel_beta4_start
  init_number sel_beta4_start                                               // initial estimate for fishery selectivity beta4 parameter in log-space
  init_4darray input_selectivity_EM(1,np,1,nreg,1,na,1,nf)                  // fishery selectivity for EM for each population, region, age, and fleet; used with select_switch_EM==0

//#############################################################################################################################################
//######################### Survey Fleet Parameter Inputs ######################################################################################
//#############################################################################################################################################

  init_matrix tsurvey_EM(1,np_em,1,nreg_em)                                 // time of survey for EM in proportion of year (0-1); 0 is Jan 1st  
  init_number lb_q                                                          // lower bound on survey catchability parameter estimate, initial estimate is q_start
  init_number ub_q                                                          // upper bound on survey catchability parameter estimate, initial estimate is q_start
  init_number q_start                                                       // initial estimate for survey catchability parameter in log-space
  init_number lb_sel_beta1_surv                                             // lower bound on survey selectivity beta1 parameter estimate, initial estimate is sel_beta1_surv _start
  init_number ub_sel_beta1_surv                                             // upper bound on survey selectivity beta1 parameter estimate, initial estimate is sel_beta1_surv _start
  init_number sel_beta1_surv_start                                          // initial estimate for survey selectivity beta1 parameter in log-space
  init_number lb_sel_beta2_surv                                             // lower bound on survey selectivity beta2 parameter estimate, initial estimate is sel_beta2_surv _start
  init_number ub_sel_beta2_surv                                             // upper bound on survey selectivity beta2 parameter estimate, initial estimate is sel_beta2_surv _start
  init_number sel_beta2_surv_start                                          // initial estimate for survey selectivity beta2 parameter in log-space
  init_number lb_sel_beta3_surv                                             // lower bound on survey selectivity beta3 parameter estimate, initial estimate is sel_beta3_surv _start
  init_number ub_sel_beta3_surv                                             // upper bound on survey selectivity beta3 parameter estimate, initial estimate is sel_beta3_surv _start
  init_number sel_beta3_surv_start                                          // initial estimate for survey selectivity beta3 parameter in log-space
  init_number lb_sel_beta4_surv                                             // lower bound on survey selectivity beta4 parameter estimate, initial estimate is sel_beta4_surv _start
  init_number ub_sel_beta4_surv                                             // upper bound on survey selectivity beta4 parameter estimate, initial estimate is sel_beta4_surv _start
  init_number sel_beta4_surv_start                                          // initial estimate for survey selectivity beta4 parameter in log-space
  init_4darray input_survey_selectivity_EM(1,np,1,nreg,1,na,1,nfs)          // survey selectivity for EM for each population, region, age, and survey fleet; used with select_switch_survey_EM==0

//#############################################################################################################################################
//######################### Natural Mortality Parameter Inputs ########################################################################################
//#############################################################################################################################################

  init_number lb_M                                                          // lower bound on natural mortality parameter estimate, initial estimate is M_start
  init_number ub_M                                                          // upper bound on natural mortality parameter estimate, initial estimate is M_start
  init_number M_start                                                       // initial estimate for natural mortality parameter in log-space
  init_matrix input_M_EM(1,np_em,1,na)                                      // natural mortality for EM for each population and age; used with M_switch_EM==(-1)

//#############################################################################################################################################
//######################### Other Demographic Parameter Inputs ################################################################################
//#############################################################################################################################################

  init_number lb_init_abund                                                 // lower bound on initial abundance parameter estimates, initial estimate is N_start
  init_number ub_init_abund                                                 // upper bound on initial abundance parameter estimates, initial estimate is N_start
  init_int N_start                                                          // initial estimate for initial abundance parameters in log-space
  init_number lb_init_dist                                                  // lower bound on initial abundance distribution parameter estimates, initial estimate is init_dist_start
  init_number ub_init_dist                                                  // upper bound on initial abundance distribution parameter estimates, initial estimate is init_dist_start
  init_int init_dist_start                                                  // initial estimate for initial abundance distribution parameters in logit transform space
  init_4darray init_abund_EM(1,np_em,1,np_em,1,nreg_em,1,na)                // initial abundance for EM for each natal population across all population and regions by age; used with init_abund_switch_EM==(-2)
  init_3darray input_dist_init_abund_EM(1,np_em,1,np_em,1,nreg_em)          // initial abundance distribution for EM for each  natal population across all population and regions by age; used when est_dist_init_abund_EM==(-2)

//#############################################################################################################################################
//######################### Likelihood Function Weights and Penalties Inputs ###############################################################################
//#############################################################################################################################################

 // Likelihood weights /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_number wt_surv                                                       // maximum likelihood added weight for survey biomass data
  init_number wt_catch                                                      // maximum likelihood added weight for fishery yield data
  init_number wt_fish_age                                                   // maximum likelihood added weight for fishery age composition data
  init_number wt_srv_age                                                    // maximum likelihood added weight for survey age composition data
  init_number wt_tag                                                        // maximum likelihood added weight for tag recapture data

 // Penalty Function Inputs /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  init_number wt_rec                                                        // maximum likelihood penalty function weight for recruitment deviations from the stock-recruit relationship
  init_number Rave_pen_switch                                               // determines whether to implement a penalty for R_ave (R0) estimation
  init_number wt_Rave_pen                                                   // maximum likelihood penalty function weight for recruitment R_ave estimation; used with Rave_pen_switch==1
  init_number Rave_mean                                                     // value of R_ave (R0) against which to penalize deviations from in log-space; used with Rave_pen_switch==1 AND wt_Rave_pen>0
  init_number R_app_pen_switch                                               // determines whether to implement a penalty for recruit apportionment estimation
  init_number wt_R_app_pen                                                   // maximum likelihood penalty function weight for recruit apportionment estimation; used with R_app_pen_switch==1
  init_number R_app_pen                                                      // value of recruit apportionment against which to penalize deviations from in logit transform-space; used with R_app_pen_switch==1 AND wt_R_app_pen>0  
  init_number abund_pen_switch                                              // determines whether to implement a penalty for initial abundance estimation
  init_number wt_abund_pen                                                  // maximum likelihood penalty function weight for initial abdundance estimation; used with abund_pen_switch==1
  init_number Mean_N                                                        // value of initial abundance against which to penalize deviations from in log-space; used with abund_pen_switch==1 AND wt_abund_pen>0
  init_number abund_dist_pen_switch                                         // determines whether to implement a penalty for initial abundance distribution estimation
  init_number wt_abund_dist_pen                                             // maximum likelihood penalty function weight for initial abdundance distribution estimation; used with abund_dist_pen_switch==1
  init_number Mean_N_dist                                                   // value of initial abundance distribution against which to penalize deviations from in log-space; used with abund_dist_pen_switch==1 AND wt_abund_dist_pen>0
  init_number move_pen_switch                                               // determines whether to implement a penalty for movement estimation
  init_number wt_move_pen                                                   // maximum likelihood penalty function weight for movement estimation; used with move_pen_switch==1 OR 2
  init_number Tpen                                                          // value of movement against which to penalize deviations from in logit transform-space; used with move_pen_switch==1 OR 2 AND wt_move_pen>0
  init_number sigma_Tpen_EM                                                 // movement variance term that controls deviations from Tpen; used when move_pen_switch==2 AND wt_move_pen>0
  init_number wt_F_pen                                                      // maximum likelihood penalty function weight for yearly fishing mortality estimates
  init_number wt_M_pen                                                      // maximum likelihood penalty function weight for natural mortality estimates
  init_number wt_B_pen                                                      // maximum likelihood penalty function weight for tag reporting rate estimates
  init_number report_rate_ave                                               // value of tag reporting rate against which to penalize deviations from; used with wt_B_pen >0
  init_number report_rate_sigma                                             // tag reporting rate variance term that controls deviations from report_rate_ave; used when wt_B_pen >0

//#############################################################################################################################################
//######################### Observation Error Parameter Inputs ###############################################################################
//#############################################################################################################################################

  init_4darray OBS_survey_fleet_bio_se_EM(1,np_em,1,nreg_em,1,ny,1,nfs_em)  // assumed lognormal variance for the survey biomass index data
  init_4darray OBS_yield_fleet_se_EM(1,np_em,1,nreg_em,1,ny,1,nf_em)        // assumed lognormal variance for the fishery yield data
  init_4darray OBS_survey_prop_N_EM(1,np_em,1,nreg_em,1,ny,1,nfs_em)        // assumed multinomial effective sample size for the survey age composition data
  init_4darray OBS_catch_prop_N_EM(1,np_em,1,nreg_em,1,ny,1,nf_em)          // assumed multinomial effective sample size for the fishery age composition data
  init_4darray tag_N_EM(1,np_em,1,nreg_em,1,ny_rel,1,na)                    // assumed multinomial effective sample size for the tagging data (effective N applies to a given cohort)

//#############################################################################################################################################
//######################### Random Number Generator Seed Inputs ###############################################################################
//#############################################################################################################################################

  init_number myseed_yield                                                  // RNG seed for fishery yield observation error
  init_number myseed_survey                                                 // RNG seed for survey biomass observation error
  init_number myseed_F                                                      // RNG seed for fishing mortality yearly deviations (stochasticity)
  init_number myseed_rec_devs                                               // RNG seed for recruitment yearly deviations (stochasticity)
  init_number myseed_rec_apport                                             // RNG seed for recruitment apportionment yearly deviations (stochasticity)
  init_number myseed_rec_index                                              // RNG seed for recruitment index observation error
  init_number myseed_survey_age                                             // RNG seed for survey age compositions observation error
  init_number myseed_catch_age                                              // RNG seed for fishery age compositions observation error
  init_number myseed_tag                                                    // RNG seed for tag recaptures observation error
  init_number myseed_ntags                                                  // RNG seed for uniform distribution of number of tags; used with number_tags_switch==4:7
  init_number myseed_prob_tag                                               // RNG seed for uniform distribution of probability of tagging by area; used with number_tags_switch==5 OR 7
  init_number myseed_prob_tag_year                                          // RNG seed for uniform distribution of prob of tagging by year; used with number_tags_switch==5:7
  init_number myseed_T                                                      // RNG seed for movement yearly deviations (stochasticity)

//#############################################################################################################################################
//######################### Debug Code to Check that Everything Read In Correctly #############################################################
//#############################################################################################################################################

  init_number debug                                                         // debug number used below to double check that parameters have been read in correctly

  !! cout << "debug = " << debug << endl;
  !! cout << "If debug != 1541 then .dat file not setup correctly" << endl;
  !! cout << "input read" << endl;
  
//#############################################################################################################################################
//######################### Counters and Indices ##############################################################################################
//#############################################################################################################################################

  vector years(1,ny)
  !!years.fill_seqadd(double(1),1.0);
  
    init_imatrix nregions_temp(1,np,1,np)                                   //used to fill tag_recap matrices

  !! for(int j=1;j<=npops;j++)                                              //recap stock
  !! {
  !!  for (int r=1;r<=npops;r++)                                            //recap region
  !!  {
  !!    if(j<=r)
  !!     {
  !!     nregions_temp(j,r)=0;
  !!     }
  !!    if(j>r)
  !!     {
  !!     nregions_temp(j,r)=nreg(r);                                        //create temp matrix that holds the number of regions that exist in all previous populations (so can sum for use in calcs below)
  !!     }
  !!   }
  !!  }
  ivector nreg_temp(1,np)

  int a
  int y
  int z
  int k
  int j
  int i
  int s
  int r
  int n
  int w
  int p
  int v
  int x
  int u
  int d
  int xx
  int region_counter                                                        // counter to enumerate the regions for the report out section for 6d arrays

INITIALIZATION_SECTION  //set initial values
     
  log_F_est F_guess;
  
PARAMETER_SECTION

  !! ivector nr=nregions;
  !! int nps=npops;
  !! int nyr=nyrs;
  !! int nag=nages;
  !! ivector nfl=nfleets;
  !! ivector nfls=nfleets_survey;
  
 init_number dummy(phase_dummy)
 init_3darray log_F_est(1,nps,1,nr,1,nfl,phase_F)
 3darray F_est(1,nps,1,nr,1,nfl)

 //For dunce cap F
 matrix Fstartyr(1,nps,1,nr)
 matrix minF_start(1,nps,1,nr)
 matrix minF_end(1,nps,1,nr)
 matrix maxF(1,nps,1,nr)
 matrix stepF_up(1,nps,1,nr)
 matrix stepF_down(1,nps,1,nr)
 vector R_ave(1,nps)
 
 // vitals
 6darray T(1,nps,1,nr,1,nyr,1,nag,1,nps,1,nr)
 5darray T_year(1,nps,1,nr,1,nyr,1,nps,1,nr)
 6darray rel_bio(1,nps,1,nr,1,nyr,1,nag,1,nps,1,nr)
 3darray Bstar(1,nps,1,nr,1,nag)
 3darray c(1,nps,1,nr,1,nag)
 4darray Fract_Move_DD(1,nps,1,nr,1,nyr,1,nag)
 5darray selectivity(1,nps,1,nr,1,nyr,1,nag,1,nfl)
 4darray selectivity_age(1,nps,1,nr,1,nag,1,nfl)
 4darray F_year(1,nps,1,nr,1,nyr,1,nfl)
 5darray F_fleet(1,nps,1,nr,1,nyr,1,nag,1,nfl)
 4darray F(1,nps,1,nr,1,nyr,1,nag)
 4darray M(1,nps,1,nr,1,nyr,1,nag)
 matrix rec_devs(1,nps,1,nyr)
 matrix rec_devs_randwalk(1,nps,1,nyr)
 4darray init_abund(1,nps,1,nps,1,nr,1,nag)
 4darray weight_population(1,nps,1,nr,1,nyr,1,nag)
 4darray weight_catch(1,nps,1,nr,1,nyr,1,nag)
 3darray wt_mat_mult(1,nps,1,nyr,1,nag)
 4darray wt_mat_mult_reg(1,nps,1,nr,1,nyr,1,nag)
 3darray ave_mat_temp(1,nps,1,nag,1,nr) //to calc average maturity
 matrix ave_mat(1,nps,1,nag) //to calc average maturity
 matrix SPR_N(1,nps,1,nag)
 matrix SPR_SSB(1,nps,1,nag)
 vector SPR(1,nps)
 vector SSB_zero(1,nps)
 vector alpha(1,nps)
 vector beta(1,nps)

//recruitment 
 3darray recruits_BM(1,nps,1,nr,1,nyr)
 3darray recruits_AM(1,nps,1,nr,1,nyr)
 3darray Rec_Prop(1,nps,1,nr,1,nyr-1)
 3darray Rec_prop_temp1(1,nps,1,nyr-1,1,nr)
 3darray Rec_prop_temp2(1,nps,1,nyr-1,1,nr)

 3darray rec_index_BM(1,nps,1,nr,1,nyr)
 3darray rec_index_prop_BM(1,nps,1,nr,1,nyr)
 3darray rec_index_BM_temp(1,nps,1,nyr,1,nr)
 3darray rec_index_AM(1,nps,1,nr,1,nyr)
 3darray rec_index_prop_AM(1,nps,1,nr,1,nyr)
 3darray rec_index_AM_temp(1,nps,1,nyr,1,nr)

 vector env_rec(1,nyr)

//abundance 
 4darray abundance_at_age_BM(1,nps,1,nr,1,nyr,1,nag)
 4darray abundance_at_age_AM(1,nps,1,nr,1,nyr,1,nag)
 4darray abundance_in(1,nps,1,nr,1,nyr,1,nag)
 4darray abundance_res(1,nps,1,nr,1,nyr,1,nag)
 4darray abundance_leave(1,nps,1,nr,1,nyr,1,nag)
 4darray abundance_spawn(1,nps,1,nr,1,nyr,1,nag)

//biomass
 4darray biomass_BM_age(1,nps,1,nr,1,nyr,1,nag)
 4darray biomass_AM_age(1,nps,1,nr,1,nyr,1,nag)
 3darray biomass_BM(1,nps,1,nr,1,nyr)
 3darray biomass_AM(1,nps,1,nr,1,nyr)
 4darray bio_in(1,nps,1,nr,1,nyr,1,nag)
 4darray bio_res(1,nps,1,nr,1,nyr,1,nag)
 4darray bio_leave(1,nps,1,nr,1,nyr,1,nag)

 //tagging data
  !! int nyr_rel=nyrs_release;
  !! ivector xy(1,nyr_rel);
  !! ivector nt(1,nyr_rel);
  !! ivector nt2(1,nyr_rel);
  !! int nt3=max_life_tags*sum(nregions)+1;
  !! ivector tag_age(1,nyr_rel);

  !!  for(int x=1; x<=nyrs_release; x++)
  !!   {
  !!    xx=yrs_releases(x);
  !!    xy(x)=min(max_life_tags,nyr-xx+1);
  !!    nt(x)=xy(x)*sum(nregions)+1;
  !!    nt2(x)=nt(x)-1;
  !!    tag_age(x)=xy(x);
  !!   }

 7darray tags_avail(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr)
 7darray recaps(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr)
 7darray tag_prop(1,nps,1,nr,1,nyr_rel,1,nag,1,tag_age,1,nps,1,nr)
 5darray tag_prop_final(1,nps,1,nr,1,nyr_rel,1,nag,1,nt)
 5darray SIM_tag_prop(1,nps,1,nr,1,nyr_rel,1,nag,1,nt)
 5darray OBS_tag_prop_final(1,nps,1,nr,1,nyr_rel,1,nag,1,nt)
 3darray total_recap_temp(1,nps,1,nr,1,tag_age)
 vector rand_tag_prop_temp2(1,nt3)  
 5darray tag_prop_temp2(1,nps,1,nr,1,nyr_rel,1,nag,1,nt2)
 4darray tag_prop_temp2_no_age(1,nps,1,nr,1,nyr_rel,1,nt2)
 5darray rand_tag_prop_temp(1,nps,1,nr,1,nyr_rel,1,nag,1,2000) //should make function of max(ncatch) but had issues making an index
 4darray rand_tag_prop_temp_no_age(1,nps,1,nr,1,nyr_rel,1,2000) //should make function of max(ncatch) but had issues making an index

 matrix tags_avail_temp(1,nps,1,nr)
 3darray tag_prop_temp(1,nps,1,max_life_tags,1,nr)

 vector ntags_total(1,nyr_rel)
 matrix ntags_total_temp(1,nps,1,nr)
 4darray ntags(1,nps,1,nr,1,nyr_rel,1,nag)
 3darray ntags_region(1,nps,1,nr,1,nyr_rel)
 3darray prob_tag(1,nps,1,nr,1,nyr_rel)
 vector prob_tag_year(1,nyr_rel)

 3darray prob_tag_RN(1,nps,1,nr,1,nyr_rel)
 3darray ntags_RN(1,nps,1,nr,1,nyr_rel)
 vector prob_tag_year_RN(1,nyr_rel)
 4darray total_rec(1,nps,1,nr,1,nyr_rel,1,nag)
 4darray not_rec(1,nps,1,nr,1,nyr_rel,1,nag)
 4darray tag_prop_not_rec(1,nps,1,nr,1,nyr_rel,1,nag)

  3darray total_rec_no_age(1,nps,1,nr,1,nyr_rel)
  3darray not_rec_no_age(1,nps,1,nr,1,nyr_rel)
  3darray ntags_no_age(1,nps,1,nr,1,nyr_rel)
  4darray tag_prop_final_no_age(1,nps,1,nr,1,nyr_rel,1,nt)
  6darray tag_prop_no_age(1,nps,1,nr,1,nyr_rel,1,tag_age,1,nps,1,nr)
  3darray tag_prop_not_rec_no_age(1,nps,1,nr,1,nyr_rel)
  7darray tag_recap_no_age_temp(1,nps,1,nr,1,nyr_rel,1,tag_age,1,nps,1,nr,1,nag) //recaps
  4darray SIM_tag_prop_no_age(1,nps,1,nr,1,nyr_rel,1,nt)
  4darray OBS_tag_prop_final_no_age(1,nps,1,nr,1,nyr_rel,1,nt)
  3darray age_full_selection(1,nps,1,nr,1,nyr)
  vector age_full_selection_temp(1,nag)
  4darray F_tag(1,nps,1,nr,1,nyr_rel,1,nag)
  6darray T_tag(1,nps,1,nr,1,nyr_rel,1,nag,1,nps,1,nr)
  
 //survey index
 5darray survey_selectivity(1,nps,1,nr,1,nyr,1,nag,1,nfls)
 4darray survey_selectivity_age(1,nps,1,nr,1,nag,1,nfls)
 5darray survey_selectivity_temp(1,nps,1,nr,1,nyr,1,nfls,1,nag)
 6darray true_survey_fleet_overlap_age(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,nag)
 6darray survey_at_age_region_fleet_overlap_prop(1,nps,1,nps,1,nr,1,nfls,1,nyr,1,nag)
 6darray SIM_survey_prop_overlap(1,nps,1,nps,1,nr,1,nfls,1,nyr,1,nag)
 6darray OBS_survey_prop_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,nag) 

 6darray true_survey_fleet_overlap_age_bio(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,nag)
 5darray true_survey_fleet_bio_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfls)
 4darray true_survey_region_bio_overlap(1,nps,1,nps,1,nyr,1,nr)
 3darray true_survey_population_bio_overlap(1,nps,1,nyr,1,nps)
 matrix true_survey_natal_bio_overlap(1,nyr,1,nps)
 vector true_survey_total_bio_overlap(1,nyr)
 5darray true_survey_fleet_age(1,nps,1,nr,1,nyr,1,nfls,1,nag)
 
 // true survey abundances for tag releases
 5darray true_survey_fleet_age_temp(1,nps,1,nr,1,nyr,1,nag,1,nfls)
 4darray true_survey_region_abundance(1,nps,1,nyr,1,nr,1,nag)
 4darray true_survey_region_abundance_temp(1,nps,1,nyr,1,nag,1,nr) 
 3darray true_survey_population_abundance_temp(1,nyr,1,nag,1,nps)
 3darray true_survey_population_abundance(1,nyr,1,nps,1,nag)
 3darray true_survey_total_abundance_temp(1,nyr,1,nag,1,nps)
 matrix true_survey_total_abundance(1,nyr,1,nag)

 5darray survey_at_age_fleet_prop(1,nps,1,nr,1,nyr,1,nfls,1,nag)
 5darray SIM_survey_prop(1,nps,1,nr,1,nfls,1,nyr,1,nag)
 5darray OBS_survey_prop(1,nps,1,nr,1,nyr,1,nfls,1,nag)

 5darray true_survey_fleet_age_bio(1,nps,1,nr,1,nyr,1,nfls,1,nag)
 4darray true_survey_fleet_bio(1,nps,1,nr,1,nyr,1,nfls)
 3darray true_survey_region_bio(1,nps,1,nyr,1,nr)
 matrix true_survey_population_bio(1,nyr,1,nps)
 vector true_survey_total_bio(1,nyr)
 5darray OBS_survey_fleet_bio_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfls)
 4darray OBS_survey_region_bio_overlap(1,nps,1,nps,1,nyr,1,nr)
 3darray OBS_survey_population_bio_overlap(1,nps,1,nyr,1,nps)
 matrix OBS_survey_natal_bio_overlap(1,nyr,1,nps)
 vector OBS_survey_total_bio_overlap(1,nyr)
 4darray OBS_survey_fleet_bio(1,nps,1,nr,1,nyr,1,nfls)
 3darray OBS_survey_region_bio(1,nps,1,nyr,1,nr)
 matrix OBS_survey_population_bio(1,nyr,1,nps)
 vector OBS_survey_total_bio(1,nyr)
 3darray apport_region_survey_biomass(1,nps,1,nr,1,nyr)

 //catch, abundance, biomass
 5darray catch_at_age_fleet(1,nps,1,nr,1,nyr,1,nag,1,nfl)
 5darray catch_at_age_fleet_prop(1,nps,1,nr,1,nyr,1,nfl,1,nag)
 5darray SIM_catch_prop(1,nps,1,nr,1,nfl,1,nyr,1,nag)
 5darray OBS_catch_prop(1,nps,1,nr,1,nyr,1,nfl,1,nag)

 4darray yield_fleet(1,nps,1,nr,1,nyr,1,nfl)
 4darray catch_at_age_region(1,nps,1,nr,1,nyr,1,nag)
 4darray catch_at_age_region_prop(1,nps,1,nr,1,nyr,1,nag)
 3darray yield_region(1,nps,1,nr,1,nyr)
 3darray catch_at_age_population(1,nps,1,nyr,1,nag)
 3darray catch_at_age_population_prop(1,nps,1,nyr,1,nag)
 matrix yield_population(1,nps,1,nyr)
 3darray SSB_region(1,nps,1,nr,1,nyr)
 matrix SSB_population(1,nps,1,nyr)
 vector SSB_total(1,nyr)
 3darray abundance_population(1,nps,1,nyr,1,nag)
 matrix abundance_total(1,nyr,1,nag)
 matrix biomass_population(1,nps,1,nyr)
 vector biomass_total(1,nyr)
 matrix catch_at_age_total(1,nyr,1,nag)
 matrix catch_at_age_total_prop(1,nyr,1,nag)
 vector yield_total(1,nyr)
 4darray harvest_rate_region_num(1,nps,1,nr,1,nyr,1,nag)
 3darray harvest_rate_population_num(1,nps,1,nyr,1,nag)
 matrix harvest_rate_total_num(1,nyr,1,nag)
 3darray harvest_rate_region_bio(1,nps,1,nr,1,nyr)
 matrix harvest_rate_population_bio(1,nps,1,nyr)
 vector harvest_rate_total_bio(1,nyr)
 3darray depletion_region(1,nps,1,nr,1,nyr)
 matrix depletion_population(1,nps,1,nyr)
 vector depletion_total(1,nyr)
 5darray abundance_at_age_BM_overlap_region(1,nps,1,nps,1,nyr,1,nag,1,nr)
 4darray abundance_at_age_BM_overlap_population(1,nps,1,nps,1,nyr,1,nag)
 5darray abundance_at_age_AM_overlap_region(1,nps,1,nps,1,nyr,1,nag,1,nr)
 4darray abundance_at_age_AM_overlap_population(1,nps,1,nps,1,nyr,1,nag)
 4darray abundance_AM_overlap_region_all_natal(1,nps,1,nr,1,nyr,1,nag)
 5darray abundance_spawn_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 4darray SSB_region_overlap(1,nps,1,nps,1,nr,1,nyr)
 3darray SSB_population_overlap(1,nps,1,nps,1,nyr)
 matrix SSB_natal_overlap(1,nps,1,nyr)

 6darray catch_at_age_region_fleet_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)
 6darray catch_at_age_region_fleet_overlap_prop(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)
 6darray SIM_catch_prop_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)
 6darray OBS_catch_prop_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)
 5darray catch_at_age_region_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 5darray catch_at_age_region_overlap_prop(1,nps,1,nps,1,nr,1,nyr,1,nag)
 5darray yield_region_fleet_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr)
 4darray yield_region_overlap(1,nps,1,nps,1,nr,1,nyr)
 4darray catch_at_age_population_overlap(1,nps,1,nps,1,nyr,1,nag)
 4darray catch_at_age_population_overlap_prop(1,nps,1,nps,1,nyr,1,nag)
 3darray yield_population_overlap(1,nps,1,nps,1,nyr)
 3darray abundance_natal_overlap(1,nps,1,nyr,1,nag)
 5darray biomass_BM_age_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 5darray biomass_AM_age_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 4darray biomass_BM_overlap_region(1,nps,1,nps,1,nr,1,nyr)
 4darray biomass_AM_overlap_region(1,nps,1,nps,1,nr,1,nyr)
 5darray biomass_AM_overlap_region_all_natal_temp(1,nps,1,nr,1,nyr,1,nag,1,nps)
 4darray biomass_AM_overlap_age_region_all_natal(1,nps,1,nr,1,nyr,1,nag)
 3darray biomass_AM_overlap_region_all_natal(1,nps,1,nr,1,nyr)
 3darray biomass_population_overlap(1,nps,1,nps,1,nyr)
 matrix biomass_natal_overlap(1,nps,1,nyr)
 3darray catch_at_age_natal_overlap(1,nps,1,nyr,1,nag)
 3darray catch_at_age_natal_overlap_prop(1,nps,1,nyr,1,nag)
 matrix yield_natal_overlap(1,nps,1,nyr)
 5darray harvest_rate_region_fleet_bio_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr)
 4darray harvest_rate_region_bio_overlap(1,nps,1,nps,1,nr,1,nyr)
 3darray harvest_rate_population_bio_overlap(1,nps,1,nps,1,nyr)
 matrix harvest_rate_natal_bio_overlap(1,nps,1,nyr)
 4darray depletion_region_overlap(1,nps,1,nps,1,nr,1,nyr)
 3darray depletion_population_overlap(1,nps,1,nps,1,nyr)
 matrix depletion_natal_overlap(1,nps,1,nyr)
 3darray Bratio_population_overlap(1,nps,1,nps,1,nyr)
 matrix Bratio_natal_overlap(1,nps,1,nyr)
 matrix Bratio_population(1,nps,1,nyr)
 vector Bratio_total(1,nyr)

 //Observed Yield
 5darray OBS_yield_region_fleet_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfl)
 4darray OBS_yield_region_overlap(1,nps,1,nps,1,nyr,1,nr)
 3darray OBS_yield_population_overlap(1,nps,1,nyr,1,nps)
 matrix OBS_yield_natal_overlap(1,nyr,1,nps)
 vector OBS_yield_total_overlap(1,nyr)
 4darray OBS_yield_fleet(1,nps,1,nr,1,nyr,1,nfl)
 3darray OBS_yield_region(1,nps,1,nyr,1,nr)
 matrix OBS_yield_population(1,nyr,1,nps)
 vector OBS_yield_total(1,nyr)
 3darray apport_yield_region(1,nps,1,nr,1,nyr)

 matrix biomass_BM_temp(1,nps,1,nr)
 number biomass_BM_temp2
 5darray biomass_BM_overlap_temp(1,nps,1,nr,1,nyr,1,nag,1,nps)
 4darray init_abund_temp(1,nps,1,nr,1,nag,1,nps)
 5darray rand_SIM_survey_prop_temp(1,nps,1,nr,1,nyr,1,nfls,1,2000) //should make function of max(ncatch) but had issues making an index, used 2000 as placeholder since nsurvey unlikely to exceed 2000
 6darray rand_SIM_survey_prop_temp_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfls,1,2000) //should make function of max(ncatch) but had issues making an index, used 2000 as placeholder since nsurvey unlikely to exceed 2000
 5darray rand_SIM_catch_prop_temp(1,nps,1,nr,1,nyr,1,nfl,1,2000) //should make function of max(ncatch) but had issues making an index
 6darray rand_SIM_catch_prop_temp_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfl,1,2000) //should make function of max(ncatch) but had issues making an index
 vector rand_SIM_survey_prop_temp2(1,nages)
 vector rand_SIM_catch_prop_temp2(1,nages)
 5darray OBS_survey_fleet_bio_temp(1,nps,1,nr,1,nyr,1,nfls,1,nps)
 5darray true_survey_fleet_bio_overlap_temp(1,nps,1,nr,1,nyr,1,nfls,1,nps)
 5darray catch_at_age_fleet_prop_temp(1,nps,1,nr,1,nyr,1,nfl,1,nag)
 matrix abundance_move_temp(1,nps,1,nr)
 matrix bio_move_temp(1,nps,1,nr)
 5darray yield_fleet_temp(1,nps,1,nr,1,nyr,1,nfl,1,nag)
 4darray yield_region_temp(1,nps,1,nr,1,nyr,1,nag)
 3darray yield_population_temp(1,nps,1,nyr,1,nag)
 4darray SSB_region_temp(1,nps,1,nr,1,nyr,1,nag)
 matrix SSB_total_temp(1,nyr,1,nps)
 4darray abundance_population_temp(1,nps,1,nyr,1,nag,1,nr)
 3darray abundance_total_temp(1,nyr,1,nag,1,nps)
 3darray biomass_population_temp(1,nps,1,nyr,1,nr)
 matrix biomass_total_temp(1,nyr,1,nps)
 3darray catch_at_age_total_temp(1,nyr,1,nag,1,nps)
 matrix yield_total_temp(1,nyr,1,nps)
 4darray catch_at_age_population_temp(1,nps,1,nyr,1,nag,1,nr)
 5darray SSB_region_temp_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 matrix abundance_move_overlap_temp(1,nps,1,nr)
 5darray OBS_yield_fleet_temp(1,nps,1,nr,1,nyr,1,nfl,1,nps)
 6darray yield_region_fleet_temp_overlap(1,nps,1,nps,1,nr,1,nfl,1,nyr,1,nag)
 5darray yield_region_temp_overlap(1,nps,1,nps,1,nr,1,nyr,1,nag)
 4darray yield_population_temp_overlap(1,nps,1,nps,1,nyr,1,nag)
 4darray abundance_natal_temp_overlap(1,nps,1,nyr,1,nag,1,nps)
 4darray biomass_population_temp_overlap(1,nps,1,nps,1,nyr,1,nr)
 3darray biomass_natal_temp_overlap(1,nps,1,nyr,1,nps)
 4darray catch_at_age_natal_temp_overlap(1,nps,1,nyr,1,nag,1,nps)
 3darray yield_natal_temp_overlap(1,nps,1,nyr,1,nps)
 5darray catch_at_age_population_temp_overlap(1,nps,1,nps,1,nyr,1,nag,1,nr)
 3darray SSB_natal_overlap_temp(1,nps,1,nyr,1,nps)
 matrix SSB_overlap_natal(1,nps,1,nr)
 5darray abundance_AM_overlap_region_all_natal_temp(1,nps,1,nr,1,nyr,1,nag,1,nps)
 3darray SSB_population_temp(1,nps,1,nyr,1,nr)
 4darray SSB_population_temp_overlap(1,nps,1,nps,1,nyr,1,nr)

 //random numbers and NR vectors
 4darray res_TAC(1,nps,1,nr,1,nfl,1,nyr)
 3darray res_u(1,nps,1,nr,1,nyr)
 number Fnew
 number delt
 number fofF
 number fprimeF
 vector fofFvect(1,nag)
 vector fprimeFhigh(1,nag)
 vector fprimeFlow(1,nag)
 4darray TAC(1,nps,1,nr,1,nfl,1,nyr)
 3darray u(1,nps,1,nr,1,nfl)
 4darray yield_RN(1,nps,1,nr,1,nyr,1,nfl)
 5darray yield_RN_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfl)
 4darray survey_RN(1,nps,1,nr,1,nyr,1,nfls)
 5darray survey_RN_overlap(1,nps,1,nps,1,nr,1,nyr,1,nfls)
 4darray F_RN(1,nps,1,nr,1,nyr,1,nfl)
 4darray T_RN(1,nps,1,nr,1,nyr,1,nag)
 matrix rec_devs_RN(1,nps,1,nyr)
 3darray Rec_apport_RN(1,nps,1,nyr-1,1,nr)
 3darray rec_index_RN(1,nps,1,nr,1,nyr)

 4darray yield_RN_temp(1,max_pops,1,max_regs,1,max_yrs,1,max_flts)
 5darray yield_RN_temp_overlap(1,max_pops,1,max_pops,1,max_regs,1,max_yrs,1,max_flts)
 4darray survey_RN_temp(1,max_pops,1,max_regs,1,max_yrs,1,max_surv_flts)
 5darray survey_RN_temp_overlap(1,max_pops,1,max_pops,1,max_regs,1,max_yrs,1,max_surv_flts)
 4darray F_RN_temp(1,max_pops,1,max_regs,1,max_yrs,1,max_flts)
 4darray T_RN_temp(1,max_pops,1,max_regs,1,max_yrs,1,max_ages)
 matrix rec_devs_RN_temp(1,max_pops,1,max_yrs)
 3darray Rec_apport_RN_temp(1,max_pops,1,max_yrs,1,max_regs)
 3darray rec_index_RN_temp(1,max_pops,1,max_regs,1,max_yrs)
 3darray prob_tag_RN_temp(1,max_pops,1,max_regs,1,max_tag_yrs)
 3darray ntags_RN_temp(1,max_pops,1,max_regs,1,max_tag_yrs)
 vector prob_tag_year_RN_temp(1,max_tag_yrs)
 
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///Setting up parameters for export to mismatch EM model (i.e., pop structure mismatches OM structure) and for running panmictic or FAA EM
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////
 4darray abund_frac_age_region(1,nps,1,nr,1,nyr,1,nag)
 3darray abund_frac_region_year(1,nps,1,nr,1,nyr)
 matrix abund_frac_region(1,nps,1,nr)

 3darray input_weight_region_temp(1,nps,1,nag,1,nr)
 matrix input_weight_region(1,nps,1,nag)
 matrix input_weight_population_temp(1,nag,1,nps)
 vector input_weight_population(1,nag)
 3darray input_catch_weight_region_temp(1,nps,1,nag,1,nr)
 matrix input_catch_weight_region(1,nps,1,nag)
 matrix input_catch_weight_population_temp(1,nag,1,nps)
 vector input_catch_weight_population(1,nag)
 3darray fecundity_region_temp(1,nps,1,nag,1,nr)
 matrix fecundity_region(1,nps,1,nag)
 matrix fecundity_population_temp(1,nag,1,nps)
 vector fecundity_population(1,nag)
 3darray maturity_region_temp(1,nps,1,nag,1,nr)
 matrix maturity_region(1,nps,1,nag)
 matrix maturity_population_temp(1,nag,1,nps)
 vector maturity_population(1,nag)
 
 number prop_fem_pan
 matrix prop_fem_temp(1,nps,1,nr)

 matrix rec_index_BM_population(1,nps,1,nyr)
 vector rec_index_pan(1,nyr)
 3darray rec_index_temp(1,nps,1,nyr,1,nr)
 matrix rec_index_temp2(1,nyr,1,nps)

 4darray selectivity_age_temp(1,nps,1,nag,1,nfl,1,nr)
 3darray selectivity_age_pop(1,nps,1,nag,1,nfl)
 4darray survey_selectivity_age_temp(1,nps,1,nag,1,nfl,1,nr)
 3darray survey_selectivity_age_pop(1,nps,1,nag,1,nfl)

//weighted average of age comps
 5darray OBS_survey_prop_temp(1,nps,1,nr,1,nyr,1,nag,1,nfls)
 4darray OBS_survey_prop_temp2(1,nps,1,nr,1,nyr,1,nag)
 4darray OBS_survey_prop_temp3(1,nps,1,nr,1,nyr,1,nag)
 4darray OBS_survey_prop_temp4(1,nps,1,nyr,1,nag,1,nr)
 3darray OBS_survey_prop_population(1,nps,1,nyr,1,nag)
 3darray OBS_survey_prop_pan_temp(1,nyr,1,nag,1,nps)
 matrix  OBS_survey_prop_pan(1,nyr,1,nag)

 5darray OBS_catch_prop_temp(1,nps,1,nr,1,nyr,1,nag,1,nfls)
 4darray OBS_catch_prop_temp2(1,nps,1,nr,1,nyr,1,nag)
 4darray OBS_catch_prop_temp3(1,nps,1,nr,1,nyr,1,nag)
 4darray OBS_catch_prop_temp4(1,nps,1,nyr,1,nag,1,nr)
 3darray OBS_catch_prop_population(1,nps,1,nyr,1,nag)
 3darray OBS_catch_prop_pan_temp(1,nyr,1,nag,1,nps)
 matrix  OBS_catch_prop_pan(1,nyr,1,nag)


//summing ntags
 4darray ntags_temp1(1,nps,1,nyr_rel,1,nag,1,nr)
 3darray ntags_population(1,nps,1,nyr_rel,1,nag)
 3darray ntags_pan_temp(1,nyr_rel,1,nag,1,nps)
 matrix  ntags_pan(1,nyr_rel,1,nag)


//summing OBS tag prop 
 5darray OBS_tag_prop_population_temp(1,nps,1,nyr_rel,1,nag,1,nt,1,nr)
 4darray OBS_tag_prop_population_temp2(1,nps,1,nyr_rel,1,nag,1,nt)
 4darray OBS_tag_prop_population_final(1,nps,1,nyr_rel,1,nag,1,nt)
 4darray OBS_tag_prop_pan_temp(1,nyr_rel,1,nag,1,nt,1,nps)
 3darray OBS_tag_prop_pan_temp2(1,nyr_rel,1,nag,1,nt)
 3darray OBS_tag_prop_pan_final_temp(1,nyr_rel,1,nag,1,max_life_tags+1)//set up a new array for summed proportion
 3darray OBS_tag_prop_pan_final(1,nyr_rel,1,nag,1,max_life_tags+1)//set up a new array for summed proportion

 //calculating true fraction natal
 3darray init_abund_reg_temp(1,nps,1,nps,1,nr)
 4darray init_abund_age_temp(1,nps,1,nag,1,nps,1,nr)
 3darray init_abund_reg_age_temp(1,nps,1,nag,1,nps)
 matrix init_abund_pop_temp(1,nps,1,nps)
 3darray frac_natal_true(1,nps,1,nps,1,nr)
 4darray frac_natal_true_age(1,nps,1,nps,1,nr,1,nag)
 matrix T_temp(1,nps,1,nr)
 number T_temp_sum
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 //temp arrays for F_MSY search
 vector yield_MSY_temp(1,nyrs_quasi_equil)


  objective_function_value f

 !! cout << "parameters set...PARAMETER SECTION FINISHED" << endl;
 
  
PROCEDURE_SECTION

  get_random_numbers();

  get_movement();

  get_selectivity();

  get_F_age();

  get_vitals();

  get_SPR();

  get_env_Rec();

  get_DD_move_parameters();

  get_abundance();

  get_rand_survey_CAA_prop();

  get_rand_CAA_prop();

  get_tag_recaptures();

  get_observed_tag_recaptures();

  evaluate_the_objective_function();

FUNCTION get_random_numbers

 // NOTE THAT IT APPEARS THE RNG JUST CONTINUES A SINGLE SET OF RANDOM NUMBERS WITHIN A DATA OBJECT
 // THUS FOR EACH DATA SET OR SIMULATED PARAMETER THE VALUES WILL NOT REPEAT SO DON'T NEED MULTIPLE RANDOM
 // NUMBERS FOR EACH DATA SET, ETC...

   random_number_generator myrand_yield(myseed_yield);
   random_number_generator myrand_survey(myseed_survey);
   random_number_generator myrand_F(myseed_F);
   random_number_generator myrand_T(myseed_T);
   random_number_generator myrand_rec_devs(myseed_rec_devs);
   random_number_generator myrand_rec_apport(myseed_rec_apport);
   random_number_generator myrand_rec_index(myseed_rec_index);
   random_number_generator myrand_ntags(myseed_ntags);
   random_number_generator myrand_prob_tag(myseed_prob_tag);
   random_number_generator myrand_prob_tag_year(myseed_prob_tag_year);

 //to avoid issues of varying RNGs across scenarios/runs (due to different input index lengths) we generate the RNG
 //and put them into temp arrays that using fixed index legnths meant to be greater than the maximum possible
 //values that would be used for a given index (note that there is no error protection built in so that if someone
 //uses a higher index than input as the max value, then RNGs will repeat)
 //the values from the temp arrays are then read into the actual RN arrays
 //this way if one run has a different index it will have the same RNGs up until the new dimension (e.g., if one run
 //has 3 fleets and another has 4, the one with 4 fleets with have the same RNGs for the first 3 fleets, then will
 //keep reading the RNG string for the 4th fleet)
  for (int p=1;p<=max_pops;p++)
   {
    for (int j=1;j<=max_pops;j++)
     {  
      for (int r=1;r<=max_regs;r++)   
       {       
        for (int y=1;y<=max_yrs;y++)
         {
            for (int x=1;x<=max_surv_flts;x++)
             {
               survey_RN_temp_overlap(p,j,r,y,x)=randn(myrand_survey);
             }
           }
          }
         }
        }

  for (int p=1;p<=max_pops;p++)
   {
    for (int j=1;j<=max_pops;j++)
     {  
      for (int r=1;r<=max_regs;r++)   
       {       
        for (int y=1;y<=max_yrs;y++)
         {
          for (int z=1;z<=max_flts;z++)
           {
               yield_RN_temp_overlap(p,j,r,y,z)=randn(myrand_yield);
             }
            }
           }
          }
         }

    for (int j=1;j<=max_pops;j++)
     {  
      for (int r=1;r<=max_regs;r++)   
       {       
        for (int y=1;y<=max_yrs;y++)
         {
            for (int x=1;x<=max_surv_flts;x++)
             {
                 survey_RN_temp(j,r,y,x)=randn(myrand_survey);
             }
           }
          }
         }

    for (int j=1;j<=max_pops;j++)
     {  
      for (int r=1;r<=max_regs;r++)   
       {       
        for (int y=1;y<=max_yrs;y++)
         {
          for (int z=1;z<=max_flts;z++)
           {
                 yield_RN_temp(j,r,y,z)=randn(myrand_yield);
                 F_RN_temp(j,r,y,z)=randn(myrand_F);
           }
          }
         }
        }

    for (int j=1;j<=max_pops;j++)
     {  
      for (int r=1;r<=max_regs;r++)   
       {  
        for(int s=1; s<=max_tag_yrs; s++)
          {
             ntags_RN_temp(j,r,s)=randu(myrand_ntags);
             prob_tag_RN_temp(j,r,s)=randu(myrand_prob_tag);
          }
         }
        }

        for(int s=1; s<=max_tag_yrs; s++)
          {
             prob_tag_year_RN_temp(s)=randu(myrand_prob_tag_year);
          }


    for (int j=1;j<=max_pops;j++)
     {  
      for (int r=1;r<=max_regs;r++)   
       {       
        for (int y=1;y<=max_yrs;y++)
         {
          for (int a=1;a<=max_ages;a++)
           {
             T_RN_temp(j,r,y,a)=randn(myrand_T);
           }
          }
         }
        }

    for (int j=1;j<=max_pops;j++)
     {  
      for (int r=1;r<=max_regs;r++)   
       {       
        for (int y=1;y<=max_yrs;y++)
         {
                 rec_index_RN_temp(j,r,y)=randn(myrand_rec_index);

               if(apportionment_type==3)//completely random apportionment
                {
                 Rec_apport_RN_temp(j,y,r)=randu(myrand_rec_apport);//generate a positive random number bw 0-1
                }
               if(apportionment_type==4)//completely random apportionment
                {
                 Rec_apport_RN_temp(j,y,r)=randn(myrand_rec_apport);//generate a positive random number bw 0-1
               }
       }
      }
     }
    for (int j=1;j<=max_pops;j++)
     {
       for (int y=1;y<=max_yrs;y++)
         {
             rec_devs_RN_temp(j,y)=randn(myrand_rec_devs);
         }
      }


 //once RNGs are created they are now input into the actual RN arrays used in the calculations based on the true
 //index lengths (
 //*******NOTE*********
 //if index length > input associated max value used in previous loop then will have issues with repeating RNGs
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {  
      for (int r=1;r<=nregions(j);r++)   
       {       
        for (int y=1;y<=nyrs;y++)
         {
          for (int z=1;z<=nfleets(j);z++)
           {
                 yield_RN_overlap(p,j,r,y,z)=yield_RN_temp_overlap(p,j,r,y,z);
            }
           }
          }
         }
        }


  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {  
      for (int r=1;r<=nregions(j);r++)   
       {       
        for (int y=1;y<=nyrs;y++)
         {
            for (int x=1;x<=nfleets_survey(j);x++)
             {
                 survey_RN_overlap(p,j,r,y,x)=survey_RN_temp_overlap(p,j,r,y,x);
             }
           }
          }
         }
        }


    for (int j=1;j<=npops;j++)
     {  
      for (int r=1;r<=nregions(j);r++)   
       {       
        for (int y=1;y<=nyrs;y++)
         {
          for (int x=1;x<=nfleets_survey(j);x++)
             {
                 survey_RN(j,r,y,x)=survey_RN_temp(j,r,y,x);
             }
           }
          }
         }

    for (int j=1;j<=npops;j++)
     {  
      for (int r=1;r<=nregions(j);r++)   
       {       
        for (int y=1;y<=nyrs;y++)
         {
          for (int z=1;z<=nfleets(j);z++)
           {
            yield_RN(j,r,y,z)=yield_RN_temp(j,r,y,z);
            F_RN(j,r,y,z)=F_RN_temp(j,r,y,z);
           }
          }
         }
        }

    for (int j=1;j<=npops;j++)
     {  
      for (int r=1;r<=nregions(j);r++)   
       {       
        for (int y=1;y<=nyrs;y++)
         {
            rec_index_RN(j,r,y)=rec_index_RN_temp(j,r,y);
             if(y>1) //rec apport is of length nyrs-1, so need to make sure y>1
              {
               if(apportionment_type==3)//completely random apportionment
                {
                 Rec_apport_RN(j,y-1,r)=Rec_apport_RN_temp(j,y-1,r);//generate a positive random number bw 0-1
                }
               if(apportionment_type==4)//completely random apportionment
                {
                 Rec_apport_RN(j,y-1,r)=Rec_apport_RN_temp(j,y-1,r);//generate a positive random number bw 0-1
               }
              }
            }
           }
          }

    for (int j=1;j<=npops;j++)
     {
       for (int y=1;y<=nyrs;y++)
         {
           rec_devs_RN(j,y)=rec_devs_RN_temp(j,y);
         }
       }

    for (int j=1;j<=npops;j++)
     {  
      for (int r=1;r<=nregions(j);r++)   
       {       
         for(int s=1; s<=nyrs_release; s++)
               {
                 ntags_RN(j,r,s)=ntags_RN_temp(j,r,s);
                 prob_tag_RN(j,r,s)=prob_tag_RN_temp(j,r,s);
               }
             }
            }

              for(int s=1; s<=nyrs_release; s++)
               {
                 prob_tag_year_RN(s)=prob_tag_year_RN_temp(s);
               }


    for (int j=1;j<=npops;j++)
     {  
      for (int r=1;r<=nregions(j);r++)   
       {       
        for (int y=1;y<=nyrs;y++)
         {
          for (int a=1;a<=nages;a++)
             {
              T_RN(j,r,y,a)=T_RN_temp(j,r,y,a);
             }
            }
           }
          }

///////BUILD MOVEMENT MATRIX////////
FUNCTION get_movement

//POSSIBLE Future ADDITIONS:
  //new functional forms
  //random walk
  //preference functions! (yea right)
  
  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
                 if(move_switch==0)  //fix at no movement
                  {
                   if(j==k && r==n)
                    {
                     T(j,r,y,a,k,n)=1;
                    }
                   else
                    {
                     T(j,r,y,a,k,n)=0;
                    }
                  }
                 if(move_switch==1) // use input movement
                  {
                   T(j,r,y,a,k,n)=input_T(j,r,a,k,n);
                  }
                 if(move_switch==2) // only allow movement within a population (ie regions within a population) not across populations based on input residency term
                  {
                   if(j==k && r==n)
                   {
                    T(j,r,y,a,k,n)=input_residency(j,r,a);
                   }
                   if(j==k && r!=n)
                   {
                    T(j,r,y,a,k,n)=(1-input_residency(j,r,a))/(nregions(j)-1);
                   }
                   if(j!=k)
                   {
                    T(j,r,y,a,k,n)=0;
                   }
                  }
                 if(move_switch==3) //symmetric movement but only allow movement within a population (ie regionp within a population) not across populations based on input residency term
                  {
                   if(j==k)
                   {
                    T(j,r,y,a,k,n)=1/(nregions(j));
                   }
                   if(j!=k)
                   {
                    T(j,r,y,a,k,n)=0;
                   }
                  }
                 if(move_switch==4) //symmetric movement across all populations and regions
                  {
                   T(j,r,y,a,k,n)=1/sum(nregions);
                  }
                 if(move_switch==5) // allow movement across all regions and populations, based on population/region specific residency
                  {
                   if(j==k && r==n)
                   {
                    T(j,r,y,a,k,n)=input_residency(j,r,a);
                   }
                   else
                   {
                    T(j,r,y,a,k,n)=(1-input_residency(j,r,a))/(sum(nregions)-1);
                   }
                  }
                 if(move_switch==9) // use input yearly  movement
                  {
                   T(j,r,y,a,k,n)=input_T_year(j,r,y,k,n);
                  }
                 if(move_switch==10) // use input year and age  movement
                  {
                   T(j,r,y,a,k,n)=T_Full_Input(j,r,a+(y-1)*nages,k,n);  //ADMB does not read in 6d arrays,  need to read in carefully
                  }


           /// Movement for ages <= age at first post-settlement movement
             if(a<=age_first_move(j) && first_post_settle_move_switch>(-1)) // allow different movement from adults
              {
               if(first_post_settle_move_switch==0) // No movement
                {
                   if(j==k && r==n)
                    {
                     T(j,r,y,a,k,n)=1;
                    }
                   else
                    {
                     T(j,r,y,a,k,n)=0;
                    }
                }
   /*           if(first_post_settle_move_switch==1) // use input movement
                {
                 if(a<age_first_move(j)) // No movement if age is less than the age of first movement
                  {
                   if(j==k && r==n)
                    {
                     T(j,r,y,a,k,n)=1;
                    }
                   else
                    {
                     T(j,r,y,a,k,n)=0;
                    }
                  }
                 if(a==age_first_move(j))
                  {
                   T(j,r,y,a,k,n)=first_move_input(j,r,y,k,n);
                  }
                }
     */
     /*          if(first_post_settle_move_switch==2) // if natal homing pop stucture, assume a given fraction of fish move back to their natal pop
                {
                 if(a<age_first_move(j)) // No movement if age is less than the age of first movement
                  {
                   if(j==k && r==n)
                    {
                     T(j,r,y,a,k,n)=1;
                    }
                   else
                    {
                     T(j,r,y,a,k,n)=0;
                    }
                  }
                 if(a==age_first_move(j))
                  {
                   if(j==k && r==n)
                    {
                     T(j,r,y,a,k,n)=1-first_return(j,r,y);
                    }
                   if()
                    {
                     T(j,r,y,a,k,n)=;                    
                    }
                  }
      */          }
  


            if(larval_move_switch>(-1)) //fill in age-1 movement if age-1 movement is different than other ages
              {
               if(a==1 && larval_move_switch==0) // allow different movement from adults
                {
                   if(j==k && r==n)
                    {
                     T(j,r,y,a,k,n)=1;
                    }
                   else
                    {
                     T(j,r,y,a,k,n)=0;
                    }
                  }
                 if(a==1 && larval_move_switch==1) // use input movement
                  {
                   T(j,r,y,a,k,n)=input_T(j,r,a,k,n);
                  }
                 if(a==1 && larval_move_switch==2) // only allow movement within a population (ie regionp within a population) not across populations based on input residency term
                  {
                   if(j==k && r==n)
                   {
                    T(j,r,y,a,k,n)=input_residency_larval(j,r);
                   }
                   if(j==k && r!=n)
                   {
                    T(j,r,y,a,k,n)=(1-input_residency_larval(j,r))/(nregions(j)-1);
                   }
                   if(j!=k)
                   {
                    T(j,r,y,a,k,n)=0;
                   }
                  }
                 if(a==1 && larval_move_switch==3) //symmetric movement but only allow movement within a population (ie regionp within a population) not across populations
                  {
                   if(j==k)
                   {
                    T(j,r,y,a,k,n)=1/(nregions(j));
                   }
                   if(j!=k)
                   {
                    T(j,r,y,a,k,n)=0;
                   }
                  }
                 if(a==1 && larval_move_switch==4) //symmetric movement across all populations and regionp
                  {
                   T(j,r,y,a,k,n)=1/sum(nregions);
                  }
                 if(larval_move_switch==5) // allow movement across all regions and populations, based on population/region specific residency
                  {
                   if(j==k && r==n)
                   {
                    T(j,r,y,a,k,n)=input_residency_larval(j,r);
                   }
                   else
                   {
                    T(j,r,y,a,k,n)=(1-input_residency_larval(j,r))/(sum(nregions)-1);
                   }
                  }
                 }
                }
               }
              }
             }
            }
           }
         
  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
               if(T(j,r,y,a,k,n)<0)
                {
                   T(j,r,y,a,k,n)=-T(j,r,y,a,k,n);
                }
               T_temp(k,n)=T(j,r,y,a,k,n);
             }
           }          
            T_temp_sum=sum(T_temp);
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              if(T_temp_sum!=1) //ensure T doesn't exceed bounds
                {
                 T(j,r,y,a,k,n)=T(j,r,y,a,k,n)/T_temp_sum;
                }
              if(T_temp_sum==0) //ensure T doesn't exceed bounds
                {
                 if(j==k && r==n)
                  {
                   T(j,r,y,a,k,n)=1.0;
                  }
                 if(j!=k || r!=n)
                  {
                   T(j,r,y,a,k,n)=0.0;
                  }
                }
             }
            }
       }
      } 
     }
    }

 if(rand_move==1)
  {
  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              if((a>1 && move_switch!=0) || (a==1 && larval_move_switch>0) || (a==1 && larval_move_switch==(-1) &&  move_switch!=0)) 
              {
               if(j==k && r==n) //add random variation onto residency term then rescale movements below so still some to 1
                {
                 T(j,r,y,a,k,n)*=mfexp(T_RN(j,r,y,a)*sigma_T(j)-0.5*square(sigma_T(j)));
                }
              }
               if(move_switch==0)
                {
                 T(j,r,y,a,k,n)=T(j,r,y,a,k,n);
                }
               if(a==1 && larval_move_switch==0)
                {
                 T(j,r,y,a,k,n)=T(j,r,y,a,k,n);
                }
               if(T(j,r,y,a,k,n)<0)
                {
                   T(j,r,y,a,k,n)=-T(j,r,y,a,k,n);
                }
               T_temp(k,n)=T(j,r,y,a,k,n);
             }
           }         
            T_temp_sum=sum(T_temp);
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              if(T_temp_sum!=1) //ensure T doesn't exceed bounds
                {
                 T(j,r,y,a,k,n)=T(j,r,y,a,k,n)/T_temp_sum;
                }
              if(T_temp_sum==0) //ensure T doesn't exceed bounds
                {
                 if(j==k && r==n)
                  {
                   T(j,r,y,a,k,n)=1.0;
                  }
                 if(j!=k || r!=n)
                  {
                   T(j,r,y,a,k,n)=0.0;
                  }
                }
             }
            }
              }
             }
            }
           }
          }

///////SELECTIVITY CALCULATIONS///////
FUNCTION get_selectivity
//POSSIBLE ADDITIONS:
  //yearly selectivity
  
//fishery selectivity
  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
        {
        age_full_selection_temp=0;
         for (int a=1;a<=nages;a++)
           {
            for (int z=1;z<=nfleets(j);z++)
              {
               if(select_switch==2) //4 parameter double logistic selectivity
                {
                 selectivity(j,r,y,a,z)=1/((1+mfexp(-sel_beta1(j,r,z)*(a-sel_beta2(j,r,z))))*(1+mfexp(-sel_beta3(j,r,z)*(a-sel_beta4(j,r,z)))));
                }
                if(select_switch==1) //two parameter logistic selectivity
                {
                 selectivity(j,r,y,a,z)=1/(1+mfexp(-sel_beta1(j,r,z)*(a-sel_beta2(j,r,z))));
                //selectivity(j,r,y,a,z)=1/(1+mfexp(-log(19)*(a-(sel_beta1(j,r,z)))/(sel_beta2(j,r,z)))); 
                }
                if(select_switch==0) //input selectivity at age constant by year
                {
                 selectivity(j,r,y,a,z)=input_selectivity(j,r,a,z);
                }
               // age_full_selection_temp(a)=selectivity(j,r,y,a,z);
               }
              }
   //    if(tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
     //    {
       //  for (int a=1;a<=nages;a++)
         //  {
           // for (int z=1;z<=nfleets(j);z++)
             // {
               //if(selectivity(j,r,y,a,z)==max(age_full_selection_temp)); //really only getting max age of last fleet but good enough for now
                //{
                 age_full_selection(j,r,y)=tag_age_sel;
                }
               }
              }


//survey selectivity 
 for (int j=1;j<=npops;j++)
    {
     for (int r=1;r<=nregions(j);r++)
      {
       for (int y=1;y<=nyrs;y++)
         {
          for (int a=1;a<=nages;a++)
            {
             for (int z=1;z<=nfleets_survey(j);z++)
               {

               if(select_switch_survey==2) //4 parameter double logistic selectivity
                {
                 survey_selectivity(j,r,y,a,z)=1/((1+mfexp(-sel_beta1_survey(j,r,z)*(a-sel_beta2_survey(j,r,z))))*(1+mfexp(-sel_beta3_survey(j,r,z)*(a-sel_beta4_survey(j,r,z)))));
                }
                if(select_switch_survey==1) //two parameter logistic selectivity
                {
                 survey_selectivity(j,r,y,a,z)=1/(1+mfexp(-sel_beta1_survey(j,r,z)*(a-sel_beta2_survey(j,r,z))));
                //selectivity(j,r,y,a,z)=1/(1+mfexp(-log(19)*(a-(sel_beta1(j,r,z)))/(sel_beta2(j,r,z)))); 
                }
                if(select_switch_survey==0) //input selectivity
                {
                 survey_selectivity(j,r,y,a,z)=input_survey_selectivity(j,r,a,z);
                }
              }
             }
            }
          }
        }

  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
          for (int a=1;a<=nages;a++)
            {
            for (int z=1;z<=nfleets(j);z++)
              {
                selectivity_age(j,r,a,z)=selectivity(j,r,nyrs,a,z);
              }
             for (int z=1;z<=nfleets_survey(j);z++)
               {
                survey_selectivity_age(j,r,a,z)=survey_selectivity(j,r,nyrs,a,z);
               }
              }
             }
            }


 
///////FISHING MORTALITY CALCULATIONS///////
FUNCTION get_F_age

  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for (int y=1;y<=nyrs;y++)
       {
        for (int a=1;a<=nages;a++)
         {
          for (int z=1;z<=nfleets(j);z++)
           { 
             if(F_switch==1) //input F directly
              {
               F_year(j,r,y,z)=input_F(j,r,z);
              }
             if(F_switch==2) //input F_MSY
              {
               F_year(j,r,y,z)=input_F_MSY(j,r,z);
              }
             if(F_switch==3) //estimate F that achieves desired Total F (i.e., F_MSY)
              {
               F_est(j,r,z)=F_max*mfexp(log_F_est(j,r,z))/(1+mfexp(log_F_est(j,r,z)));  //multinomial logit transform to bound Fs between 0 and 5 (i.e., avoid neg Fs and large Fs)
               F_year(j,r,y,z)=F_est(j,r,z);
              }
             if(F_switch==4) //split constant input_F_const by npops
              {
               F_year(j,r,y,z)=input_F_const/npops;
              }
             if(F_switch==5) //split constant input_F_const by nregions
              {
               F_year(j,r,y,z)=input_F_const/sum(nregions);
              }
             if(F_switch==6) //split constant input_F_const by nfleets
              {
               F_year(j,r,y,z)=input_F_const/sum(nfleets);
              }
             if(F_switch==7) //F devs about input F based on sigma_F
              {
               F_year(j,r,y,z)=input_F(j,r,z)*mfexp(F_RN(j,r,y,z)*sigma_F(j,r,z)-0.5*square(sigma_F(j,r,z)));
              }
             if(F_switch==8) //random walk in F
              {
               if(y==1)
                {
                 F_year(j,r,y,z)=input_F(j,r,z);
                }
               if(y>1)
                {
                 F_year(j,r,y,z)=F_rho(j,r,z)*F_year(j,r,y-1,z)*mfexp(F_RN(j,r,y,z)*sigma_F(j,r,z)-0.5*square(sigma_F(j,r,z)));
                }
              }
             if(F_switch==9)  //Dunce cap F
              {
               Fstartyr(j,r)=dunce_F(j,r,1);
               minF_start(j,r)=dunce_F(j,r,2);
               maxF(j,r)=dunce_F(j,r,3);
               minF_end(j,r)=dunce_F(j,r,4);
               stepF_up(j,r)=(maxF(j,r)-minF_start(j,r))/((nyrs-Fstartyr(j,r))/2);
               stepF_down(j,r)=(maxF(j,r)-minF_end(j,r))/((nyrs-Fstartyr(j,r))/2);

              if(y<Fstartyr(j,r))
               {
                F_year(j,r,y,z)=0;
               }
               if(y>=Fstartyr(j,r))
               {
               if(y<((nyrs-Fstartyr(j,r))/2+Fstartyr(j,r)))
                {
                 F_year(j,r,y,z)=(minF_start(j,r)+(y-Fstartyr(j,r))*stepF_up(j,r))*mfexp(F_RN(j,r,y,z)*sigma_F(j,r,z)-0.5*square(sigma_F(j,r,z)));
                }
               }
               if(y>=((nyrs-Fstartyr(j,r))/2+Fstartyr(j,r)))
                {
                 F_year(j,r,y,z)=(maxF(j,r)-((y-Fstartyr(j,r))-((nyrs-Fstartyr(j,r))/2))*stepF_down(j,r))*mfexp(F_RN(j,r,y,z)*sigma_F(j,r,z)-0.5*square(sigma_F(j,r,z)));
                  if(F_year(j,r,y,z)<minF_end(j,r)) //needed because the stepF decrease can be randomly be greater than preceding F, and so F goes negative
                  {
                  F_year(j,r,y,z)=minF_end(j,r)*mfexp(F_RN(j,r,y,z)*sigma_F(j,r,z)-0.5*square(sigma_F(j,r,z)));
                  }
                }
              }
             F_fleet(j,r,y,a,z)=F_year(j,r,y,z)*selectivity(j,r,y,a,z);
             F(j,r,y,a)=sum(F_fleet(j,r,y,a)); 
             M(j,r,y,a)=input_M_TRUE(j,a);
           }
          }
         }
        }
       }


FUNCTION get_vitals
//POSSIBLE ADDITIONS:
  //random walk in apportionment or random to give time-varying
  //switch for input recruitment devs by year to recreate a given population trajectory

//need to add all_natal calcs for all values that might want to compare across an area, because regular non-natal
//homing calcs do not account for different weights across natal populations so if if any part of calc involves weight
//it needs to use the biomass_all_natal value not just biomass_AM


 R_ave=ln_R_ave; ///this is annoying...//a quick fix
 
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {  
      for (int r=1;r<=nregions(j);r++)   
       {       
        for (int y=1;y<=nyrs;y++)
         {
          for (int a=1;a<=nages;a++)
           {
            for (int z=1;z<=nfleets(j);z++)
             {
              weight_population(j,r,y,a)=input_weight(j,r,a);
              weight_catch(j,r,y,a)=input_catch_weight(j,r,a);

              if(maturity_switch_equil==0) // for SPR calculations when maturity across areas is equal or if want a straight average of maturity across areas
               {
                if(SSB_type==1) //fecundity based SSB
                 {
                  ave_mat_temp(j,a,r)=prop_fem(j,r)*fecundity(j,r,a)*maturity(j,r,a);//rearranging for summing
                  ave_mat(j,a) = sum(ave_mat_temp(j,a))/nregions(j); //average maturity across regions
                  wt_mat_mult(j,y,a)=ave_mat(j,a);//for SPR calcs
                 }
               if(SSB_type==2) //weight based SSB
                {
                  ave_mat_temp(j,a,r)=prop_fem(j,r)*weight_population(j,r,y,a)*maturity(j,r,a);//rearranging for summing
                  ave_mat(j,a) = sum(ave_mat_temp(j,a))/nregions(j); //average maturity across regions
                  wt_mat_mult(j,y,a)=ave_mat(j,a);//for SPR calcs
                }
               }
              if(maturity_switch_equil==1)
               {// calculates the weighted average matruity based on equilibrium apportionment of SSB - allows for unequal influence of maturity/weight
                if(SSB_type==1) //fecundity based SSB
                 {
                  ave_mat_temp(j,a,r)=prop_fem(j,r)*fecundity(j,r,a)*maturity(j,r,a)*equil_ssb_apport(j,r);//rearranging for summing
                  ave_mat(j,a) = sum(ave_mat_temp(j,a))/nregions(j); //average maturity across regions
                  wt_mat_mult(j,y,a)=ave_mat(j,a);//for SPR calcs
                 }
               if(SSB_type==2) //weight based SSB
                {
                  ave_mat_temp(j,a,r)=prop_fem(j,r)*weight_population(j,r,y,a)*maturity(j,r,a);//rearranging for summing
                  ave_mat(j,a) = sum(ave_mat_temp(j,a))/nregions(j); //average maturity across regions
                  wt_mat_mult(j,y,a)=ave_mat(j,a);//for SPR calcs
                }
               }                        
               if(recruit_devs_switch==0)  //use population recruit relationship directly
                {
                 rec_devs(j,y)=1;
                }
               if(recruit_devs_switch==1)  // allow lognormal error around SR curve
                {
                 rec_devs(j,y)=mfexp(rec_devs_RN(j,y)*sigma_recruit(j)-.5*square(sigma_recruit(j)));

                 if(recruit_randwalk_switch==1)
                 {
                  rec_devs_randwalk(j,y)=rec_devs(j,y);
                  if(y!=1)
                   {
                    rec_devs(j,y)=rec_devs(j,y-1)*rec_devs_randwalk(j,y); //is this correct?
                   }
                 }
                }
               if(SSB_type==1) //fecundity based SSB
                {
                 wt_mat_mult_reg(j,r,y,a)=prop_fem(j,r)*fecundity(j,r,a)*maturity(j,r,a);// for yearly SSB calcs
                }
               if(SSB_type==2) //weight based SSB
                {
                 wt_mat_mult_reg(j,r,y,a)=prop_fem(j,r)*weight_population(j,r,y,a)*maturity(j,r,a);
                }
             }
            }
           }

        for (int y=1;y<=nyrs-1;y++)
         {
               if(apportionment_type==1) //input recruitment apportionment directly by population and region
                {
                 Rec_Prop(j,r,y)=input_Rec_prop(j,r);
                }
               if(apportionment_type==2) //equal apportionment by nregions
                {
               Rec_Prop(j,r,y)=1.0/nregions(j);
                }
               if(apportionment_type==(-1)) // no apportionment
                {
                 Rec_Prop(j,r,y)=1;
                }
               if(apportionment_type==3)//completely random apportionment
                {
                Rec_prop_temp1(j,y,r)=Rec_apport_RN(j,y,r); //generate a positive random number bw 0-1
                }                 
               if(apportionment_type==4) //add input observation error to input recruit proportions following Schnute and Richards 1995 (as implemented in Cox and Kronlund 2008)
                {
                 Rec_prop_temp1(j,y,r)=log(input_Rec_prop(j,r))+sigma_rec_prop(j)*Rec_apport_RN(j,y,r);//applying the additive error in log space; this equals "log(u) + error" in Cox and Kronlund Table 1
                }
               }   
             }
        
  for (int r=1;r<=nregions(j);r++)
       {
        for (int y=1;y<=nyrs-1;y++)
         {
          if(apportionment_type==3)
           { //need to standardize year by region matrix to ensure new proportions sum to one
            Rec_Prop(j,r,y)=Rec_prop_temp1(j,y,r)/sum(Rec_prop_temp1(j,y));
           }
          if(apportionment_type==4)
           { //need to run through region by year matrix to calculate second half of Schnute and Richards 1995 random mult equations and to standardize randomized apportioments to total one
            Rec_prop_temp2(j,y,r)=Rec_prop_temp1(j,y,r)-(sum(Rec_prop_temp1(j,y))/nregions(j)); //finish equation T1.21 in Cox and Kronlund Table 1 (note that here nregions is the same as A (number of ages) in Table 1 paper
            Rec_Prop(j,r,y)=mfexp(Rec_prop_temp2(j,y,r))/sum(mfexp(Rec_prop_temp2(j,y))); // back-transform and standardize
           }
         }
        }
       }
     }
//SPR calcs are done with either  average maturity/weight across all the regions within a population or assuming an input population fraction at equilibrium
// while the full SSB calcs use the region specific maturity/weight
FUNCTION get_SPR
     for (int k=1;k<=npops;k++)
     {
      for (int n=1;n<=nages;n++)
       {
        if(n==1)
         {
          SPR_N(k,n)=1000;
          SPR_SSB(k,n)=wt_mat_mult(k,1,n)*SPR_N(k,n);
         }
        if(n>1 && n<nages)
         {
          SPR_N(k,n)=SPR_N(k,n-1)*mfexp(-M(k,1,1,n-1));
          SPR_SSB(k,n)=wt_mat_mult(k,1,n)*SPR_N(k,n);
         }
        if(n==nages)
         {
          SPR_N(k,n)=SPR_N(k,n-1)*mfexp(-M(k,1,1,n))*(1/(1-mfexp(-M(k,1,1,n))));
          SPR_SSB(k,n)=wt_mat_mult(k,1,n)*SPR_N(k,n);
         }
       }
     SPR(k)=sum(SPR_SSB(k))/1000;
     SSB_zero(k)=SPR(k)*R_ave(k);
      if(Rec_type==2) //BH recruitment
      {
      alpha(k)=(SSB_zero(k)/R_ave(k))*((1-steep(k))/(4*steep(k)));//alternate parameterization
      beta(k)=(5*steep(k)-1)/(4*steep(k)*R_ave(k));
      }
    }


FUNCTION get_env_Rec // calculate autocorrelated recruitment - input period and amplitude
       for (int p=1;p<=npops;p++)
        {
        for (int j=1;j<=npops;j++)
         {
         for (int r=1;r<=nregions(j);r++)
          {
          for (int y=1;y<=nyrs;y++)
           {
           if(y==1)
             {
              env_rec(y)=input_init_abund(p,j,r,1);
             }           
           if(y>1)
            {
             env_rec(y) = R_ave(j)+(R_ave(j)*amplitude)*sin(((2*M_PI)/freq)*years(y)+freq);
            }
           }
         }
        }
       }
FUNCTION get_DD_move_parameters
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////DD Movement///////////////////////////////////////////////////////////////////////////////////////////////////
           /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
            /// because want to weight T and make more fish move to area with lower bio compared to Bstar
            /// idea is that if B/Bstar<1 then not densely populated and more fish would move there           
      if(move_switch==8)
       {
        for (int s=1;s<=npops;s++)
         {
         for (int a=1;a<=nages;a++)
          {
          for (int n=1;n<=nregions(s);n++)
           {       
            if(use_input_Bstar==0 && Rec_type==2) //set Bstar==SSB0 # THIS IS NOT TESTED, SSB_ZERO CALCS ARE NOT WELL DEFINED IN SPATIAL CONTEXT
             {
              if(nregions(s)==1)
               {
                Bstar(s,n,a)=SSB_zero(s);
               }
              if(nregions(s)>1)
               {
                Bstar(s,n,a)=SSB_zero(s)*SSB_zero_appor(s,n);
               }
             }
            if(use_input_Bstar==1 || Rec_type!=2) // THIS IS TYPICALLY USED FOR DD MOVEMENT
             {
              Bstar(s,n,a)=input_Bstar(s,n,a);
             }
            if(DD_move_age_switch==0) //DD movement is not age based (based on total bio not bio at age); MAKE SURE A is constant across ages
             {
              c(s,n,a)=-log(A(s,n,a))/sum(Bstar(s,n));  /// Make sure input A is constant across ages, so c is constant; Bstar is the point of inflection in T
             }
            if(DD_move_age_switch==1) //DD movement is  age based (based on bio at age)
             {
              c(s,n,a)=-log(A(s,n,a))/Bstar(s,n,a);  ///Bstar is the point of inflection in T
             }
          }
         }
        }
       }
FUNCTION get_abundance
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // Be Aware of Difference Between Natal and Non-natal (Metapop) Calcs
   // regular non-natal homing calcs do not account for different weights/maturity across natal populationps (only account for biology of current pop)
   // so if if any part of calc involves weight/maturity it needs to use the biomass_all_natal value not just biomass_AM 
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   for (int y=1;y<=nyrs;y++)
    {
     for (int p=1;p<=npops;p++)
      {
       for (int j=1;j<=npops;j++)
        {
         for (int r=1;r<=nregions(j);r++)
          {
           for (int a=1;a<=nages;a++)
            {
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////// Fill in initial abundance in first year /////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              if(y==1)
               {
                if(init_abund_switch==1)  //initial abundance is exponential decay from Rave or R0
                 {
                  init_abund(p,j,r,a)=R_ave(p)*pow(mfexp(-(M(p,r,y,a))),a-1)*frac_natal(p,j,r);  //distributes natal abundance across all pop and regions
                 }                
                if(init_abund_switch==0)  //input initial abundance
                 {
                  init_abund(p,j,r,a)=input_init_abund(p,j,r,a);
                 }
               //Natal homing calculations
                   abundance_at_age_BM_overlap_region(p,j,y,a,r)=init_abund(p,j,r,a);

               //NON-NATAL Homing (metapopulation) calcs /////////////////////////////////////////////////////////////// //////////////////////////////////////////////////////////////////////////////    
                   init_abund_temp(j,r,a,p)=init_abund(p,j,r,a);
                   abundance_at_age_BM(j,r,y,a)=sum(init_abund_temp(j,r,a));  //sum across natal populations in a given population unit j

////////////////////// ADD A natal recruit matrix so can track recruit by natal source and a self-recruitment index
                   recruits_BM(j,r,y)=abundance_at_age_BM(j,r,y,1);
////////////////////////////////////////////////////////////////////////////////
                 } // end Y==1 loop

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////// Fill in  abundance in subsequent years /////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               if(y>1)
                {
                 if(a==1) //because year loop doesn't close until after SSB calcs, SSB matrix will be filled in for previous year and can be called in recruit eqns below
                  {
////////////////////// ADD A natal recruit matrix so can track recruit by natal source and a self-recruitment index
                 
                   if(natal_homing_switch>0) //natal homing occurs, use natal pop demographics and only recruit to natal pop
                    {
                     if(p==j) //recruits only go to natal populations (i.e., no recruits assigned to non-natal pops)
                      {
                       if(Rec_type==1) //average recruitment
                        {
                         if(apportionment_type!=0) //use Rec_Prop calculated in Get_Vitals section to apportion recruitment among regions within a population
                          {
                           recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*Rec_Prop(j,r,y-1);  //(no recruit apportionment in first year so need y-1)
                          }
                         if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                          {
                           recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1))); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
                          }
                        }

                       if(Rec_type==2) //BH recruitment  
                        {
                         if(apportionment_type!=0)  //use prespecified Rec_Prop to apportion recruitment among regions within a population
                          {
                           recruits_BM(j,r,y)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)))*rec_devs(j,y)*Rec_Prop(j,r,y-1);
                          }
                         if(apportionment_type==0) //use relative SSB to apportion recruitment among regions within a population
                          {
                           recruits_BM(j,r,y)=((SSB_population_overlap(p,j,y-1))/(alpha(j)+beta(j)*SSB_population_overlap(p,j,y-1)))*rec_devs(j,y)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1)));
                          }
                        }
                  
                       if(Rec_type==3) //environmental recruitment //HAS NOT BEEN TESTED
                        {
                         if(apportionment_type!=0) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                          {
                           recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y)*Rec_Prop(j,r,y-1);
                          }
                         if(apportionment_type==0)  //use relative SSB to apportion recruitment among regions within a population
                          {
                           recruits_BM(j,r,y)=env_rec(y)*rec_devs(j,y)*(SSB_region_overlap(p,j,r,y-1)/sum(SSB_region_overlap(p,j,y-1))); //assume with natal homing that fish aren't from a particular region so when apportion them use relative SSB in each region (but don't account for SSB that moves back to the region to spawn ie fish that move back add to population SSB but not region SSB)
                          }
                        }
                      }
                    }

                   if(natal_homing_switch==0) //for non-natal homing population structures (i.e., metapopulation), use demographics of current pop not natal pop
                    {
                     if(Rec_type==1) //average recruitment
                      {
                       if(apportionment_type!=0) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                        {
                         recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*Rec_Prop(j,r,y-1);
                        }
                       if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                        {
                         recruits_BM(j,r,y)=R_ave(j)*rec_devs(j,y)*(SSB_region(j,r,y-1)/SSB_population(j,y-1));
                        }
                      }
                      
                     if(Rec_type==2) //BH recruitment
                      {
                       if(apportionment_type!=0)  //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                        {
                         recruits_BM(j,r,y)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)))*rec_devs(j,y)*Rec_Prop(j,r,y-1);
                        }
                       if(apportionment_type==0) //use relative SSB to apportion recruitment among regionp within a population
                        {
                         recruits_BM(j,r,y)=((SSB_population(j,y-1))/(alpha(j)+beta(j)*SSB_population(j,y-1)))*rec_devs(j,y)*(SSB_region(j,r,y-1)/SSB_population(j,y-1));
                        }
                      }

                     if(Rec_type==3) //average recruitment
                      {
                       if(apportionment_type!=0) //use prespecified Rec_Prop to apportion recruitment among regionp within a population
                        {
                         recruits_BM(j,r,y)=env_rec(y)* rec_devs(j,y)*Rec_Prop(j,r,y-1);
                        }
                       if(apportionment_type==0)  //use relative SSB to apportion recruitment among regionp within a population
                        {
                         recruits_BM(j,r,y)=env_rec(y)* R_ave(j)*rec_devs(j,y)*(SSB_region(j,r,y-1)/SSB_population(j,y-1));
                        }
                      }
                    }

                  /////a==1 overlap calcs/////////////////////////////////////////////////////////////////////////////////
                   if(p==j)
                    {
                     abundance_at_age_BM_overlap_region(p,j,y,a,r)=recruits_BM(j,r,y); //only recruits in  natal pop count towards recruitment
                    }
                   if(p!=j)
                    {
                     abundance_at_age_BM_overlap_region(p,j,y,a,r)=0;
                    }

                   //// metapop calcs /////////////////////////////////////////////////////////////////////////
                    abundance_at_age_BM(j,r,y,a)=recruits_BM(j,r,y);
                   /////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  }  //end age 1 loop
       
                 if(a>1 && a<nages)
                  {
                  ///// overlap calcs/////////////////////////////////////////////////////////////////////////////////
                   abundance_at_age_BM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)));
                  
                  /////////metapop/////////////////////////////////////////////////////////
                   abundance_at_age_BM(j,r,y,a)=abundance_at_age_AM(j,r,y-1,a-1)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)));
                  /////////////////////////////////////////////////////////////////////////////////////////
                  }
        
                 if(a==nages) //account for fish already in plus group
                  {
                  ///// overlap calcs/////////////////////////////////////////////////////////////////////////////////
                   abundance_at_age_BM_overlap_region(p,j,y,a,r)=abundance_at_age_AM_overlap_region(p,j,y-1,a-1,r)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+abundance_at_age_AM_overlap_region(p,j,y-1,a,r)*mfexp(-(M(j,r,y-1,a)+F(j,r,y-1,a)));
              
                  ///////metapop/////////////////////////////////////////////////////////
                   abundance_at_age_BM(j,r,y,a)=abundance_at_age_AM(j,r,y-1,a-1)*mfexp(-(M(j,r,y-1,a-1)+F(j,r,y-1,a-1)))+abundance_at_age_AM(j,r,y-1,a)*mfexp(-(M(j,r,y-1,a)+F(j,r,y-1,a)));
                  /////////////////////////////////////////////////////////////////////////////////////////
                  }
                 } // end if y=1 loop, rest of calcs for all years
                 
                  ///// overlap calcs/////////////////////////////////////////////////////////////////////////////////
                  abundance_at_age_BM_overlap_population(p,j,y,a)=sum(abundance_at_age_BM_overlap_region(p,j,y,a));
                  biomass_BM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_BM_overlap_region(p,j,y,a,r); //maintain demographics of natal population p

                 if(natal_homing_switch==0) // no natal homing, biomass summed across all natal populations
                  {
                   biomass_BM_age(j,r,y,a)=weight_population(j,r,y,a)*abundance_at_age_BM(j,r,y,a);  //demographics based on current pop unit j
                  }
                 if(natal_homing_switch>0) // natal homing,  biomass summed across natal population (based on demographics of natal population) into biomass_BM matrix
                  {
                   biomass_BM_overlap_temp(j,r,y,a,p)=biomass_BM_age_overlap(p,j,r,y,a);
                   biomass_BM_age(j,r,y,a)=sum(biomass_BM_overlap_temp(j,r,y,a));
                  }

 ///get age one recruitment index
//////////////////////////////////////////////////////////////////
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/////////////////Make age-1 settlement index by natal pop p, too so can identify self-recruitment, etc.////////////////////////////////////////////////////////////////

               rec_index_BM(j,r,y) = recruits_BM(j,r,y)*mfexp(rec_index_RN(j,r,y)*rec_index_sigma(j,r)-0.5*square(rec_index_sigma(j,r)));
               rec_index_BM_temp(j,y,r)=rec_index_BM(j,r,y);
               rec_index_prop_BM(j,r,y)=rec_index_BM(j,r,y)/sum(rec_index_BM_temp(j,y));


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
  
            } //end age loop

             biomass_BM_overlap_region(p,j,r,y)=sum(biomass_BM_age_overlap(p,j,r,y)); //all non-age based calcs moved out of age loop so don't repeat unnecessarily
             biomass_BM(j,r,y)=sum(biomass_BM_age(j,r,y));
          }
        }
      } //close loops so have full biomass vectors filled in at start of DD movement calcs
          
     for (int a=1;a<=nages;a++)
      {
       for (int p=1;p<=npops;p++) //natal population
        {
         for (int j=1;j<=npops;j++) //current population
          {
           for (int r=1;r<=nregions(j);r++) // current region within pop j
            {
             ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
             ///////Density-Dependent Movement Calcs///////////////////////////////////////////////////////////////////////////////////////////////////
             /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /// rel_bio gives ratio of biomass in each area weighted by ratio to Bstar, use 1-ratio of bio
                /// because want to weight T and make more fish move to area with lower bio compared to Bstar
                /// idea is that if B/Bstar<1 then not densely populated and more fish would move there
                /// all calcs are independent of natal populations because assuming everything is inherent physical properties of the geographic area within j or r, not the biological demographics
                /// in other words carrying capacity simply refers to fish that area within j can support regardless of pop source or biology of natal population p found within j
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                
             if(move_switch==8 && DD_move_age_switch==1) //DD movement is age based (based on age based biomass not total bio)
              {
               Fract_Move_DD(j,r,y,a)=(1-DD_residency(j,r,a))/(1+A(j,r,a)*mfexp(c(j,r,a)*biomass_BM_age(j,r,y,a)));  // fraction of population (i.e., T) that moves based on exponential function of biomass; c is calculated previosuly based on Bstar and A
               
                if(rand_move==1) //incorporate stochasticity
                 {
                  Fract_Move_DD(j,r,y,a)*=mfexp(T_RN(j,r,y,a)*sigma_T(j)-0.5*square(sigma_T(j))); //add random variation to proportion moving
                 }
                 
                if(Fract_Move_DD(j,r,y,a)>1) //ensure that movement is not >1
                 {
                  Fract_Move_DD(j,r,y,a)=1;
                 }
                if(Fract_Move_DD(j,r,y,a)<0)  //or movement is not<0
                 {
                  Fract_Move_DD(j,r,y,a)=0;
                 }
                 
                  biomass_BM_temp=0;
                  biomass_BM_temp2=0;
                for (int s=1;s<=npops;s++)
                 {
                  for (int n=1;n<=nregions(s);n++)
                   {       
                    biomass_BM_temp(s,n)=biomass_BM_age(s,n,y,a)/Bstar(s,n,a);  //tally the bio in each pop, reg compared to associated carrying capacity by age
                   }
                 }
                 
                   biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM_age(j,r,y,a)/Bstar(j,r,a)); //do not include population/region moving from in summation since this only applies to the prop of pop that has already been tallied as moving out of current pop
                   
                for (int s=1;s<=npops;s++)
                 {
                  for (int n=1;n<=nregions(s);n++)
                   {
                    if(j==s && r==n) //do not include current pop/region in calcs (residency is already defined)
                     {
                      rel_bio(j,r,y,a,s,n)=0;  //rel bio defined for movement out of pop j region r as suitability of destination population s, region n
                     }
                    if(j!=s || r!=n) //if movement out of current population, define relative suitability by 1-B/carrying capacity (i.e., suitability is higher for lower relative densities)
                     {
                      rel_bio(j,r,y,a,s,n)=1-((biomass_BM_temp(s,n))/(biomass_BM_temp2));
                     }
                     
                    if(sum(nregions)==2)  //if there are only two regions (i.e., 2 populations, 1 region each or 1 pop with 2 regions), then set rel_bio to one so all fish that move go to opposite region (don't stay resident)
                     {
                      rel_bio(j,r,y,a,s,n)=1;
                     }
                   }
                 }
                 
                for (int s=1;s<=npops;s++)
                 {
                  for (int n=1;n<=nregions(s);n++)
                   {
                    if(a>1 || (a==1 && first_post_settle_move_switch!=0)) //replace movement from earlier movement section with DD movement here
                     {
                      if(j==s && r==n) //if no movement
                       {
                        T(j,r,y,a,s,n)=1-Fract_Move_DD(j,r,y,a); // proportion that remain resident
                       }
                      if(j!=s || r!=n) //if movement occurs
                       {
                        T(j,r,y,a,s,n)=rel_bio(j,r,y,a,s,n)*Fract_Move_DD(j,r,y,a);  //movement is the fraction of fish moving mult by the suitability of the pop, reg
                       }
                     }
                     
                    if(a==1 && first_post_settle_move_switch==0) //maintain input age-1 movement in movement section earlier (i.e., likely no movement)
                     {
                      T(j,r,y,a,s,n)=T(j,r,y,a,s,n);
                     }
                   }
                 }
              }

             if(move_switch==8 && DD_move_age_switch==0) //DD movement is not age based (based on total bio not bio at age) //same calcs as above, but now movement is constant across ages
              {
               Fract_Move_DD(j,r,y,a)=(1-DD_residency(j,r,a))/(1+A(j,r,a)*mfexp(c(j,r,a)*biomass_BM(j,r,y)));  //MAKE SURE A and DD_residency are constant across ages otherwise T will vary with age; c is calculated previosuly based on Bstar and A
               
                if(rand_move==1)
                 {
                  Fract_Move_DD(j,r,y,a)*=mfexp(T_RN(j,r,y,1)*sigma_T(j)-0.5*square(sigma_T(j))); //only use RN from one age to ensure constant values across ages
                 }
                 
                if(Fract_Move_DD(j,r,y,a)>1)
                 {
                  Fract_Move_DD(j,r,y,a)=1;
                 }
                if(Fract_Move_DD(j,r,y,a)<0)
                 {
                  Fract_Move_DD(j,r,y,a)=0;
                 }
                 
                  biomass_BM_temp=0;
                  biomass_BM_temp2=0;
                for (int s=1;s<=npops;s++)
                 {
                  for (int n=1;n<=nregions(s);n++)
                   {       
                    biomass_BM_temp(s,n)=biomass_BM(s,n,y)/sum(Bstar(s,n)); //sum Bstar so only one value (no variation by age)
                   }
                 }
                 
                    biomass_BM_temp2=sum(biomass_BM_temp)-(biomass_BM(j,r,y)/sum(Bstar(j,r))); //do not include population/region moving from in summation
                    
                for (int s=1;s<=npops;s++)
                 {
                  for (int n=1;n<=nregions(s);n++)
                   {
                    if(j==s && r==n) //do not include current pop/region in calcs (residency is already defined)
                     {
                      rel_bio(j,r,y,a,s,n)=0;  //rel bio defined for movement out of pop j region r as suitability of destination population s, region n
                     }
                    if(j!=s || r!=n) 
                     {
                      rel_bio(j,r,y,a,s,n)=1-((biomass_BM_temp(s,n))/(biomass_BM_temp2));
                     }
                     
                    if(sum(nregions)==2)
                     {
                      rel_bio(j,r,y,a,s,n)=1;
                     }
                   }
                 }
                 
                for (int s=1;s<=npops;s++)
                 {
                  for (int n=1;n<=nregions(s);n++)
                   {
                    if(a>1 || (a==1 && first_post_settle_move_switch!=0))
                     {
                      if(j==s && r==n) 
                       {
                        T(j,r,y,a,s,n)=1-Fract_Move_DD(j,r,y,a);
                       }
                      if(j!=s || r!=n) 
                       {
                        T(j,r,y,a,s,n)=rel_bio(j,r,y,a,s,n)*Fract_Move_DD(j,r,y,a);
                       }
                     }
                     
                    if(a==1 && first_post_settle_move_switch==0)
                     {
                      T(j,r,y,a,s,n)=T(j,r,y,a,s,n);
                     }
                   }
                 }
              }
            }
          }
        }
      } //close DD calcs loops so T matrix is fully filled for movement calcs
       /////////// END DD CALCS ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
    for (int p=1;p<=npops;p++) //natal population
     {
      for (int j=1;j<=npops;j++) //current population
       {
        for (int r=1;r<=nregions(j);r++) // current region within pop j
         {
          for (int a=1;a<=nages;a++)
           {
          ///// overlap movement calcs/////////////////////////////////////////////////////////////////////////////////           
              abundance_move_overlap_temp=0;
             for (int k=1;k<=npops;k++)
              {
               for (int n=1;n<=nregions(k);n++)  //calc all fish that are moving out of all pops k, regs n into pop j, reg r
                {              
                 if(move_switch!=6 || move_switch!=7)  //all movement scenarios not involving natal return or no movement
                  {
                   if(natal_homing_switch>0) //for natal homing pop sturctures assume movement is based on natal pop
                    {
                     abundance_move_overlap_temp(k,n)=abundance_at_age_BM_overlap_region(p,k,y,a,n)*T(p,n,y,a,j,r); 
                    }
                   if(natal_homing_switch==0) //for non-natal homing pop structures assume movement is based on current population //do this calc to ensure even when no natal homing that still filling overlap matrices with appropriate natal specific values
                    {
                     abundance_move_overlap_temp(k,n)=abundance_at_age_BM_overlap_region(p,k,y,a,n)*T(k,n,y,a,j,r);
                    }
                  }
                 if(move_switch==6) //assumes natal return occurs at return_age so that prop of pop returns at that age, otherwise no movement for other ages
                  {
                   if(a==return_age && p==j && p==k && j==k)
                    {
                     abundance_move_overlap_temp(k,n)=0; // fish already in natal pop are already accounted for
                    }
                   if(a==return_age && p==j && j!=k)
                    {
                     abundance_move_overlap_temp(k,n)=abundance_at_age_BM_overlap_region(p,k,y,a,n)*return_probability(p); // if not in natal pop then certain proprotion return
                    }
                  }
                }
              }

          ///// Update abundance matrices with movement from above calcs /////////////////////////////////////////////////////////////////////////////////           
             if(move_switch!=6 || move_switch!=7)  //regular movement, just add all fish moving to the population,reg they are moving to 
              {
               abundance_at_age_AM_overlap_region(p,j,y,a,r)=sum(abundance_move_overlap_temp);
              }
             if(move_switch==7)  //all fish stay where they were (i.e., no return migration)
              {
               abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_BM_overlap_region(p,j,y,a,r);
              }
             if(move_switch==6) // perform one time return migration to natal population
              {
               if(a<return_age || a>return_age) //if before return age or after return age then no movement
                {
                 abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_BM_overlap_region(p,j,y,a,r);                    
                }
               if(a==return_age && p==j) //add returning fish to natal pop
                {
                 abundance_at_age_AM_overlap_region(p,j,y,a,r)=abundance_at_age_BM_overlap_region(p,j,y,a,r)+(sum(abundance_move_overlap_temp)/nregions(p)); //fish that remain in natal pop plus those that return to natal pop; returning fish split evenly among regions within natal pop
                }
               if(a==return_age && p!=j) //if not in natal pop then fish staying put is 1-return_prob
                {
                 abundance_at_age_AM_overlap_region(p,j,y,a,r)=(1-return_probability(p))*abundance_at_age_BM_overlap_region(p,j,y,a,r);
                }
              }

          ///// overlap abundance/biomass calcs/////////////////////////////////////////////////////////////////////////////////           
                abundance_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r);
                
                 if(natal_homing_switch>0) //for natal homing pop sturctures assume movement is based on natal pop
                  {
                   biomass_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*weight_population(p,r,y,a);  //weight is based on natal pop
                  }
                 if(natal_homing_switch==0) //for non-natal homing pop structures assume demographics are based on current population //do this calc to ensure even when no natal homing that still filling overlap matrices with appropriate natal specific values
                  {
                   biomass_AM_overlap_region_all_natal_temp(j,r,y,a,p)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*weight_population(j,r,y,a);  //weight is based on current pop
                  }

                abundance_AM_overlap_region_all_natal(j,r,y,a)=sum(abundance_AM_overlap_region_all_natal_temp(j,r,y,a)); //abundance across all natal pop in current pop unit
                biomass_AM_overlap_age_region_all_natal(j,r,y,a)=sum(biomass_AM_overlap_region_all_natal_temp(j,r,y,a)); //biomass across all natal pop in current pop unit
                abundance_at_age_AM_overlap_population(p,j,y,a)=sum(abundance_at_age_AM_overlap_region(p,j,y,a));        //abundance from a given natal population in a given pop unit summed across all regions in j pop unit
                biomass_AM_age_overlap(p,j,r,y,a)=weight_population(p,r,y,a)*abundance_at_age_AM_overlap_region(p,j,y,a,r); //biomass from a given natal population in a given pop unit summed across all regions in j pop uni

        ///////////MIGHT WANT TO DO FOLLOWING CALCS FOR NATAL HOMING AS WELL (DO METAPOP TYPE CALCS BELOW)
              //  abundance_in abundance_res abundance_leave  bio_in  bio_res  bio_leave                
        //////////////////////////////////////////////////////////////////////////////

         //////////////////////////////////////////////////////////////////////////////
         ////NON-NATAL Homing movement calcs and putting natal homing abundance into area abundance/////////////////////////////////////////////////////////////// //////////////////////////////////////////////////////////////////////////////
         //////////////////////////////////////////////////////////////////////////////
    
                abundance_move_temp=0;
                bio_move_temp=0;

             for (int k=1;k<=npops;k++)
              {
               for (int n=1;n<=nregions(k);n++) //calc all fish that are moving out of all pops k, regs n into pop j, reg r
                {
                 abundance_move_temp(k,n)=abundance_at_age_BM(k,n,y,a)*T(k,n,y,a,j,r);
                 bio_move_temp(k,n)=abundance_move_temp(k,n)*weight_population(k,n,y,a);                   
                }
              }
              
             if(natal_homing_switch>0) //if natal homing put abundance summed across natal population by region into abundance at age AM
              {
               abundance_at_age_AM(j,r,y,a)=abundance_AM_overlap_region_all_natal(j,r,y,a);
               biomass_AM(j,r,y)= biomass_AM_overlap_region_all_natal(j,r,y);
              }
             if(natal_homing_switch==0) //no natal homing
              {
               abundance_at_age_AM(j,r,y,a)=sum(abundance_move_temp);
               biomass_AM_age(j,r,y,a)=weight_population(j,r,y,a)*abundance_at_age_AM(j,r,y,a);
               biomass_AM(j,r,y)=sum(biomass_AM_age(j,r,y));
              }
                   
                abundance_population_temp(j,y,a,r)=abundance_at_age_AM(j,r,y,a);
                abundance_population(j,y,a)=sum(abundance_population_temp(j,y,a));
                abundance_total_temp(y,a,j)=abundance_population(j,y,a);
                abundance_total(y,a)=sum(abundance_total_temp(y,a));

                abundance_in(j,r,y,a)=sum(abundance_move_temp)-abundance_move_temp(j,r);
                abundance_res(j,r,y,a)=abundance_move_temp(j,r);
                abundance_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)-abundance_res(j,r,y,a);
                bio_in(j,r,y,a)=sum(bio_move_temp)-bio_move_temp(j,r);
                bio_res(j,r,y,a)=bio_move_temp(j,r);
                bio_leave(j,r,y,a)=abundance_at_age_BM(j,r,y,a)*weight_population(j,r,y,a)-bio_res(j,r,y,a);

               //calculate age-1 recruit index with observation error
                recruits_AM(j,r,y)=abundance_at_age_AM(j,r,y,1);                
                rec_index_AM(j,r,y)=recruits_AM(j,r,y)*mfexp(rec_index_RN(j,r,y)*rec_index_sigma(j,r)-0.5*square(rec_index_sigma(j,r)));
                rec_index_AM_temp(j,y,r)=rec_index_AM(j,r,y);
                rec_index_prop_AM(j,r,y)=rec_index_AM(j,r,y)/sum(rec_index_AM_temp(j,y));


             for (int z=1;z<=nfleets_survey(j);z++)    /// survey index //THESE SURVEY CALCS ARE NEEDED HERE FOR THE MSY CALCS THAT FOLLOW
              {
               if(tsurvey(j,r)==0) //if survey at beggining of year, do calcs without temporal adjustment for mortality
                {
                ///// overlap survey calcs/////////////////////////////////////////////////////////////////////////////////           
                 true_survey_fleet_overlap_age(p,j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM_overlap_region(p,j,y,a,r)*q_survey(j,r,z);
                 true_survey_fleet_overlap_age_bio(p,j,r,y,z,a)=true_survey_fleet_overlap_age(p,j,r,y,z,a)*weight_population(p,r,y,a);
                 true_survey_fleet_bio_overlap(p,j,r,y,z)=sum(true_survey_fleet_overlap_age_bio(p,j,r,y,z));  
                 true_survey_fleet_bio_overlap_temp(j,r,y,z,p)=true_survey_fleet_bio_overlap(p,j,r,y,z);

                //add observation error for observed survey index //USE THIS DATA IF ASSUME HAVE GENETIC DATA TO DIFFERENTIATE SURVEY DATA BY NATAL POPULATION
                 OBS_survey_fleet_bio_overlap(p,j,r,y,z)=true_survey_fleet_bio_overlap(p,j,r,y,z)*mfexp(survey_RN_overlap(p,j,r,y,z)*sigma_survey_overlap(p,j,r,z)-.5*square(sigma_survey_overlap(p,j,r,z)));
                 OBS_survey_fleet_bio_temp(j,r,y,z,p)=OBS_survey_fleet_bio_overlap(p,j,r,y,z);

               if(natal_homing_switch==0) //for non-natal homing pop structures
                {
                 true_survey_fleet_age(j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM(j,r,y,a)*q_survey(j,r,z);
                 true_survey_fleet_age_bio(j,r,y,z,a)=true_survey_fleet_age(j,r,y,z,a)*weight_population(j,r,y,a);                  
                 true_survey_fleet_bio(j,r,y,z)=sum(true_survey_fleet_age_bio(j,r,y,z));
                 OBS_survey_fleet_bio(j,r,y,z)=true_survey_fleet_bio(j,r,y,z)*mfexp(survey_RN(j,r,y,z)*sigma_survey(j,r,z)-.5*square(sigma_survey(j,r,z)));
                }
               if(natal_homing_switch==1) // for natal homing structures OBS survey is summed over natal populations in a given pop unit j //WOULD FIT THIS DATA IF ASSUME NO GENETIC DATA TO DIFFERENTIATE THEN FIT BY NATAL POP
                {
                 true_survey_fleet_bio(j,r,y,z)=sum(true_survey_fleet_bio_overlap_temp(j,r,y,z));
                 OBS_survey_fleet_bio(j,r,y,z)=sum(OBS_survey_fleet_bio_temp(j,r,y,z));  
                }
                         
                 // calculate true survey abundance for tagging simulation **dh**
                 true_survey_fleet_age_temp(j,r,y,a,z)=true_survey_fleet_age(j,r,y,z,a);
                 true_survey_region_abundance(j,y,r,a)=sum(true_survey_fleet_age_temp(j,r,y,a));
                 true_survey_region_abundance_temp(j,y,a,r)=true_survey_region_abundance(j,y,r,a);
                 true_survey_population_abundance(y,j,a)=sum(true_survey_region_abundance_temp(j,y,a));
                 true_survey_population_abundance_temp(y,a,j)=true_survey_population_abundance(y,j,a);
                 true_survey_total_abundance(y,a)=sum(true_survey_population_abundance_temp(y,a));

                }  //tsurvey==0
              } //end survey_fleets
           } //end age loop
           
               ///// overlap bio calcs/////////////////////////////////////////////////////////////////////////////////           
                biomass_AM_overlap_region_all_natal(j,r,y)=sum(biomass_AM_overlap_age_region_all_natal(j,r,y));
                biomass_AM_overlap_region(p,j,r,y)=sum(biomass_AM_age_overlap(p,j,r,y));
                biomass_population_temp_overlap(p,j,y,r)=biomass_AM_overlap_region(p,j,r,y);
                biomass_population_overlap(p,j,y)=sum(biomass_population_temp_overlap(p,j,y));
                biomass_natal_temp_overlap(p,y,j)=biomass_population_overlap(p,j,y);
                biomass_natal_overlap(p,y)=sum(biomass_natal_temp_overlap(p,y));

               ///// metapop biomass calcs/////////////////////////////////////////////////////////////////////////////////           
                biomass_population_temp(j,y,r)=biomass_AM(j,r,y);
                biomass_population(j,y)=sum(biomass_population_temp(j,y));
                biomass_total_temp(y,j)=biomass_population(j,y);
                biomass_total(y)=sum(biomass_total_temp(y));

               if(tsurvey(j,r)==0) //if survey at beggining of year, do calcs without temporal adjustment for mortality
                {
               ///// overlap survey calcs/////////////////////////////////////////////////////////////////////////////////           
                 true_survey_region_bio_overlap(p,j,y,r)=sum(true_survey_fleet_bio_overlap(p,j,r,y));               
                 true_survey_population_bio_overlap(p,y,j)=sum(true_survey_region_bio_overlap(p,j,y));               
                 true_survey_natal_bio_overlap(y,p)=sum(true_survey_population_bio_overlap(p,y));               
                 true_survey_total_bio_overlap(y)=sum(true_survey_natal_bio_overlap(y));
                 OBS_survey_region_bio_overlap(p,j,y,r)=sum(OBS_survey_fleet_bio_overlap(p,j,r,y));
                 OBS_survey_population_bio_overlap(p,y,j)=sum(OBS_survey_region_bio_overlap(p,j,y));
                 OBS_survey_natal_bio_overlap(y,p)=sum(OBS_survey_population_bio_overlap(p,y));
                 OBS_survey_total_bio_overlap(y)=sum(OBS_survey_natal_bio_overlap(y));

               ///// metapop survey calcs/////////////////////////////////////////////////////////////////////////////////           
                 true_survey_region_bio(j,y,r)=sum(true_survey_fleet_bio(j,r,y));
                 true_survey_population_bio(y,j)=sum(true_survey_region_bio(j,y));
                 true_survey_total_bio(y)=sum(true_survey_population_bio(y));
                 OBS_survey_region_bio(j,y,r)=sum(OBS_survey_fleet_bio(j,r,y));
                 OBS_survey_population_bio(y,j)=sum(OBS_survey_region_bio(j,y));
                 OBS_survey_total_bio(y)=sum(OBS_survey_population_bio(y));

                 apport_region_survey_biomass(j,r,y)=OBS_survey_region_bio(j,y,r)/OBS_survey_population_bio(y,j); //use this if doing MSY search so need these survey values before MSY
                } 
         }
       }
     } //end pop loop

 ///////////////////////////////////////////////////////////////////////////////////
 ////////////////NEWTON-RAPHSON Calcs for Determining F corresponding to desired TAC Values ////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////////////////////////////
    ///***** These Calcs are used to determine the F that achieves a desired TAC or harvest rate (u) //////////////////////////////////
    // ** Only used if the MSY search switches are turned on /////////////////////////////////////
    // These calcs are needed in the middle of abundance calcs, because they rely on current year biomass or survey //////
    // they then do search for the F rate that achieves the input TAC and use this in remaining abundance calcs, etc... //////////
    // So new F needs to be determined mid cycle and then inserted into later calculations //////////////////////////////////
    // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // These sections were used in the Goethel and Berger 2017; Bosley et al., 2020 papers to compare reference points across spatial structure and movement scenarios ////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int p=1;p<=npops;p++)
     {
      for (int j=1;j<=npops;j++)
       {
        for (int r=1;r<=nregions(j);r++)
         {
          if(MSY_model_type_switch>0 && MSY_model_type_switch<3)  //Find the fully selected F that would result in the TAC
           {
            for (int a=1;a<=nages;a++)
             {
              for(int x=1;x<=nfleets(j);x++)
               {
                if(MSY_model_type_switch==1) // use TAC to set F
                 {                   
                    if(parse_TAC==0) //use input TAC and DO NOT parse it based on data
                     {
                      if(calc_TAC_from_uMSY==0) // just use input TAC
                       {
                        TAC(j,r,x,y)=input_TAC(j,r,x);
                       }
                      if(calc_TAC_from_uMSY==1) //calc the TAC from input harvest rate and pop biomass
                       {
                        TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y);
                       }
                      if(calc_TAC_from_uMSY==2) //calc TAC from input harvest rate and regional biomass
                       {
                        TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y);
                       }
                     }
                     
                    if(parse_TAC==1) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                     {
                      if(parse_TAC_source==0) //parse using recruit index before movement, but assume occurs during spawning so inherent 1 year lag
                       {
                        if(calc_TAC_from_uMSY==0) //parse INPUT TAC
                         {
                          if(y==1) //no data in first year for partition so split evenly among regions
                           {
                            TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                           }
                          if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                           {
                            TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_TAC(j,r,x)/nfleets(j); //partition regional quota evenly by fleets in region
                           }
                         }
                        if(calc_TAC_from_uMSY==1) //calc TAC from input harvest and pop biomass
                         {
                          if(y==1)
                           {
                            TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                           }
                          if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                           {
                            TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j); //partition regional quota evenly by fleets in region
                           }
                         }
                        if(calc_TAC_from_uMSY==2) // calc TAC from input harvest rate and reg biomass
                         {
                          if(y==1)
                           {
                            TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                           }
                          if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                           {
                            TAC(j,r,x,y)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j); //partition regional quota evenly by fleets in region
                           }
                         }
                       }
                       
                      if(parse_TAC_source==1) //parse TAC using recruitment index after movement
                       {
                        if(calc_TAC_from_uMSY==0) //parse INPUT TAC
                         {
                          if(y==1)
                           {
                            TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                           }
                          if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                           {
                            TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_TAC(j,r,x)/nfleets(j); //partition regional quota evenly by fleets in region
                           }
                         }
                        if(calc_TAC_from_uMSY==1) // calc TAC from input u and pop biomass
                         {
                          if(y==1)
                           {
                            TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                           }
                          if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                           {
                            TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j); //partition regional quota evenly by fleets in region
                           }
                         }
                        if(calc_TAC_from_uMSY==2) //calc TAC from input u and reg biomass
                         {
                          if(y==1)
                           {
                            TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j)); //partition evenly by total fleets in area
                           }
                          if(y>1) //assume rec index occurs at time of spawning so never have available in current year
                           {
                            TAC(j,r,x,y)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j); //partition regional quota evenly by fleets in region
                           }
                         }
                       }
                       
                      if(parse_TAC_source==2) //parse TAC using survey biomass
                       {
                        if(calc_TAC_from_uMSY==0) //use input TAC
                         {
                          if(TAC_survey_parse_timelag_switch==1) //assume timelag between survey year and when available to parse TAC
                           {
                            if(y<=TAC_survey_parse_timelag) //apportion TAC equally among fleets if y<timelag
                             {                                                     
                              TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                             }
                            if(y>TAC_survey_parse_timelag) // use survey data from timelag year (usually assume 1 year lag)
                             {                                                     
                              TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_TAC(j,r,x)/nfleets(j);
                             }
                           }
                          if(TAC_survey_parse_timelag_switch==0) //no timelag
                           {
                            if(tsurvey(j,r)==0) //if survey at beginning of year then use current year relative survey biomass to parse
                             {
                              TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_TAC(j,r,x)/nfleets(j);
                             }
                            if(tsurvey(j,r)>0) //if survey after start of year then in first year assume equal parsing among regions
                             {
                              if(y==1) //first year apportion TAC equally among fleets 
                               {                                                     
                                TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                               }
                              if(y>1) //use relative survey biomass from previous year
                               {                                                     
                                TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_TAC(j,r,x)/nfleets(j);
                               }
                             }
                           }
                         }
                         
                        if(calc_TAC_from_uMSY==1) //parse TAC using input u and population biomass
                         {
                          if(TAC_survey_parse_timelag_switch==1) //assume timelag when survey data available
                           {
                            if(y<=TAC_survey_parse_timelag) //for years less than timelag assume equal parsing among regions, fleets
                             {                                                     
                              TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                             }
                            if(y>TAC_survey_parse_timelag) //parse by survey in timelag year
                             {                                                     
                              TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                             }
                           }
                          if(TAC_survey_parse_timelag_switch==0) //no timelag
                           {
                            if(tsurvey(j,r)==0) //if survey at beginning of year use current year relative biomass 
                             {
                              TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                             }
                            if(tsurvey(j,r)>0)  //if survey not at start of year use equal parsing for first year, then prev year rel survey biomass to parse
                             {
                              if(y==1) //first year apportion TAC equally among fleets
                               {                                                     
                                TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                               }
                              if(y>1) //use prev year rel survey bio
                               {                                                     
                                TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)*biomass_population(j,y)/nfleets(j);
                               }
                             }
                           }
                         }
                         
                        if(calc_TAC_from_uMSY==2) //use input u and regional biomass to parse TAC
                         {
                          if(TAC_survey_parse_timelag_switch==1) //assume timelag in survey data availability
                           {
                            if(y<=TAC_survey_parse_timelag) //in years before survey data available assume equal parsing
                             {                                                     
                              TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                             }
                            if(y>TAC_survey_parse_timelag) //in years after timelag use survey data from timelag year
                             {                                                     
                              TAC(j,r,x,y)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                             }
                           }
                          if(TAC_survey_parse_timelag_switch==0) //no timelag
                           {
                            if(tsurvey(j,r)==0) //if survey at start of year use current year rel survey bio to parse
                             {
                              TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                             }
                            if(tsurvey(j,r)>0) //if survey in mid year, then assume equal parsing for first year and prev year rel survey bio for subsequent years
                             {
                              if(y==1) //first year apportion TAC equally among fleets
                               {                                                     
                                TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                               }
                              if(y>1) //use prev year rel survey bio
                               {                                                     
                                TAC(j,r,x,y)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)*biomass_AM(j,r,y)/nfleets(j);
                               }
                             }
                           }
                         }
                       }
                       
                      if(parse_TAC_source==3) //equal parsing across all pop, reg, fleets
                       {
                        if(calc_TAC_from_uMSY==0) //use input TAC
                         {
                          TAC(j,r,x,y)=input_TAC(j,r,x)/(nregions(j)*nfleets(j));
                         }
                        if(calc_TAC_from_uMSY==1) //use input u and pop bio
                         {
                          TAC(j,r,x,y)=input_u(j,r,x)*biomass_population(j,y)/(nregions(j)*nfleets(j));
                         }
                        if(calc_TAC_from_uMSY==2) //use input u and reg bio
                         {
                          TAC(j,r,x,y)=input_u(j,r,x)*biomass_AM(j,r,y)/(nregions(j)*nfleets(j));
                         }
                       }
                     }
                     
                  //Use Newton-Raphson search to find new F that achieves the TAC     
                   if(TAC(j,r,x,y)==0) //iteration may have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                    {
                     Fnew=0;
                    }
                   if(TAC(j,r,x,y)>0) 
                    {
                     Fnew=Fnew_start; //F searching for is set to initial guess

                     for(int i=1;i<=NR_iterations;i++)  //newton-raphson iteration
                      {
                       delt=Fnew*NR_dev;  //step size,  NR_dev~0.001
                       
                        for(int s=1;s<=nages;s++) //calc catch (using baranov's catch eqn) for the new F, and a step above and below new F
                         {
                          if(natal_homing_switch>0) //for natal homing pop structures use biomass summed over all natal regions in a given pop unit j
                           {
                            fofFvect(s)=((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                            fprimeFhigh(s)=(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                            fprimeFlow(s)=(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                           }
                          if(natal_homing_switch==0) //no natal homing, use regular abundance
                           {
                            fofFvect(s)=weight_catch(j,r,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                            fprimeFhigh(s)=weight_catch(j,r,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                            fprimeFlow(s)=weight_catch(j,r,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s))));
                           }
                         }
                      }
                         
                        fofF=sum(fofFvect)-TAC(j,r,x,y); //diff between catch at new F and desired TAC
                        fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt); //derivative at new F
                        Fnew=Fnew-(fofF/fprimeF); //set new F to old F minus diff/deriv
                        
                        if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF)); //if step size brings neg reduce step
                        if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                      }
                    } 
                                                                 
                if(MSY_model_type_switch==2) //do same calcs, but based on input harvest rate (u) instead of TAC
                 {
                  if(parse_TAC==0) //do not parse the u, just use input values (use this if are inputting uMSY directly)
                   {
                    u(j,r,x)=input_u(j,r,x);
                   }
                  if(parse_TAC==1) //use observed data to parse input_TAC by region // if multiple fleets then evenly distribute across fleets
                   {
                    if(parse_TAC_source==0) //parse using recruit index before movement
                     {
                      if(y==1)
                       {
                        u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                       }
                      if(y>1)
                       {
                        u(j,r,x)=rec_index_prop_BM(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                       }
                     }
                    if(parse_TAC_source==1) //parse using recruit index after movement
                     {
                      if(y==1)
                       {
                        u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                       }
                      if(y>1)
                       {
                        u(j,r,x)=rec_index_prop_AM(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                       }
                     }
                    if(parse_TAC_source==2) //parse using survey biomass
                     {
                      if(TAC_survey_parse_timelag_switch==1) //use timelag
                       {
                        if(y<=TAC_survey_parse_timelag) //first year apportion TAC equally among fleets
                         {                                                     
                          u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                         }
                        if(y>TAC_survey_parse_timelag)
                         {                                                     
                          u(j,r,x)=apport_region_survey_biomass(j,r,(y-TAC_survey_parse_timelag))*input_u(j,r,x)/nfleets(j);
                         }
                       }
                      if(TAC_survey_parse_timelag_switch==0) //no timelag
                       {
                        if(tsurvey(j,r)==0)
                         {
                          u(j,r,x)=apport_region_survey_biomass(j,r,y)*input_u(j,r,x)/nfleets(j);
                         }
                        if(tsurvey(j,r)>0)
                         {
                          if(y==1) //first year apportion TAC equally among fleets
                           {                                                     
                            u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                           }
                          if(y>1)
                           {                                                     
                            u(j,r,x)=apport_region_survey_biomass(j,r,y-1)*input_u(j,r,x)/nfleets(j);
                           }
                         }
                       }
                     }
                     
                    if(parse_TAC_source==3) //parse equally across regions and fleets
                     {
                      u(j,r,x)=input_u(j,r,x)/(nregions(j)*nfleets(j));
                     }
                   }
                   
                  if(u(j,r,x)==0) //iterations have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
                   {
                    Fnew=0;
                   }
                  if(u(j,r,x)>0) 
                   {
                    Fnew=Fnew_start;
                    
                    for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterations
                     {
                      delt=Fnew*NR_dev;  // NR_dev~0.001
                      for(int s=1;s<=nages;s++)
                       {
                        if(natal_homing_switch>0) //for natal homing pop structures use abundnace summed across all natal pops in a given pop unit j
                         {
                          fofFvect(s)=(((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                          fprimeFhigh(s)=((((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                          fprimeFlow(s)=((((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*biomass_AM_overlap_age_region_all_natal(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM_overlap_region_all_natal(j,r,y);
                         }
                        if(natal_homing_switch==0) //for non-natal homing pop structures just use total abundance in a given pop unit j
                         {
                          fofFvect(s)=(weight_catch(j,r,y,s)*((Fnew*selectivity(j,r,y,s,x))/(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*(Fnew*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                          fprimeFhigh(s)=(weight_catch(j,r,y,s)*(((Fnew+delt)*selectivity(j,r,y,s,x))/((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew+delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                          fprimeFlow(s)=(weight_catch(j,r,y,s)*(((Fnew-delt)*selectivity(j,r,y,s,x))/((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))*abundance_at_age_AM(j,r,y,s)*(1-mfexp(-1*((Fnew-delt)*selectivity(j,r,y,s,x)+M(j,r,y,s)))))/biomass_AM(j,r,y);
                         }
                       } 
                        fofF=sum(fofFvect)-u(j,r,x);
                        fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                        Fnew=Fnew-(fofF/fprimeF);
                        if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                        if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                     }
                   } 
                 }
                 F_fleet(j,r,y,a,x)=Fnew*selectivity(j,r,y,a,x); //replace fleet F with new F that achieves TAC or u                                        
               }
               F(j,r,y,a)=sum(F_fleet(j,r,y,a)); //replace total F with new Fleet Fs that achieved TAC or u summed across fleets
             }
           }
         }
       }                     
     } //finish TAC calcs

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int p=1;p<=npops;p++)
     {
      for (int j=1;j<=npops;j++)
       {
        for (int r=1;r<=nregions(j);r++)
         {
          for (int a=1;a<=nages;a++)
           {
             ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
             ///////////////SSB calcs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
             ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
             //for natal homing pop structure we don't adjust SSB for fish not in natal population area because SR calculations do this automatically by only using /////////////////
             //SSB that is in natal population area (i.e., p==j) ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
             //for spawning migration scenarios, the SSB is adjusted outside the loop to remove SSB of fish that actually returned (to ensure SSB is correctly reported at end) //////////////////////
             // to natal population (i.e., remove this SSB from the non-natal areas...doesn't impact SR calcs so can do this outside loops without conpequence to model /////////////////////////////////////
             ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

             ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
             ////////////////////// NATAL HOMING CALCS ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
             ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

               abundance_spawn_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(p));
               SSB_region_temp_overlap(p,j,r,y,a)=abundance_spawn_overlap(p,j,r,y,a)*wt_mat_mult_reg(p,r,y,a); //use  mat of natal pop

               abundance_natal_temp_overlap(p,y,a,j)=abundance_at_age_AM_overlap_population(p,j,y,a);
               abundance_natal_overlap(p,y,a)=sum(abundance_natal_temp_overlap(p,y,a));

               catch_at_age_region_overlap(p,j,r,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F(j,r,y,a)+M(j,r,y,a))))*(F(j,r,y,a)/(F(j,r,y,a)+M(j,r,y,a))); //
               catch_at_age_population_temp_overlap(p,j,y,a,r)=catch_at_age_region_overlap(p,j,r,y,a);
               catch_at_age_population_overlap(p,j,y,a)=sum(catch_at_age_population_temp_overlap(p,j,y,a));
               catch_at_age_natal_temp_overlap(p,y,a,j)=catch_at_age_population_overlap(p,j,y,a);
               catch_at_age_natal_overlap(p,y,a)=sum(catch_at_age_natal_temp_overlap(p,y,a));
               
               yield_population_temp_overlap(p,j,y,a)=weight_catch(p,r,y,a)*catch_at_age_population_overlap(p,j,y,a); //use natal pop weight
               yield_region_temp_overlap(p,j,r,y,a)=weight_catch(p,r,y,a)*catch_at_age_region_overlap(p,j,r,y,a);     //use natal pop weight
               
              for (int z=1;z<=nfleets(j);z++)
               {
                catch_at_age_region_fleet_overlap(p,j,r,z,y,a)=abundance_at_age_AM_overlap_region(p,j,y,a,r)*(1.0-mfexp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z)/(F_fleet(j,r,y,a,z)+M(j,r,y,a)));              

                yield_region_fleet_temp_overlap(p,j,r,z,y,a)=weight_catch(p,r,y,a)*catch_at_age_region_fleet_overlap(p,j,r,z,y,a);  //use natal pop weight
               }

             ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
             ////////////////////// NON-NATAL HOMING CALCS ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
             ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                abundance_spawn(j,r,y,a)=abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tspawn(j));

               if(natal_homing_switch==0)  //no natal homing use current pop demographics
                {
                 SSB_region_temp(j,r,y,a)=abundance_spawn(j,r,y,a)*wt_mat_mult_reg(j,r,y,a); // use mat of current pop
                }

               for (int z=1;z<=nfleets(j);z++)
                {
                 catch_at_age_fleet(j,r,y,a,z)=abundance_at_age_AM(j,r,y,a)*(1.0-exp(-(F_fleet(j,r,y,a,z)+M(j,r,y,a))))*(F_fleet(j,r,y,a,z))/(F(j,r,y,a)+M(j,r,y,a));

                 yield_fleet_temp(j,r,y,z,a)=weight_catch(j,r,y,a)*catch_at_age_fleet(j,r,y,a,z); // use mat of current pop
                }

                 catch_at_age_region(j,r,y,a)=sum(catch_at_age_fleet(j,r,y,a));
                 catch_at_age_population_temp(j,y,a,r)=catch_at_age_region(j,r,y,a);
                 catch_at_age_population(j,y,a)=sum(catch_at_age_population_temp(j,y,a));
                 catch_at_age_total_temp(y,a,j)=catch_at_age_population(j,y,a);
                 catch_at_age_total(y,a)=sum(catch_at_age_total_temp(y,a));
                 
                 yield_region_temp(j,r,y,a)=weight_catch(j,r,y,a)*catch_at_age_region(j,r,y,a); // use mat of current pop
                 yield_population_temp(j,y,a)=weight_catch(j,r,y,a)*catch_at_age_population(j,y,a); // use mat of current pop

                 harvest_rate_region_num(j,r,y,a)=catch_at_age_region(j,r,y,a)/abundance_at_age_AM(j,r,y,a);
                 harvest_rate_population_num(j,y,a)=catch_at_age_population(j,y,a)/abundance_population(j,y,a);
                 harvest_rate_total_num(y,a)=catch_at_age_total(y,a)/abundance_total(y,a);
                 
               if(tsurvey(j,r)>0) //if survey not at beginning of year do time adjustments for mortality at time of survey
                {               
                 for (int z=1;z<=nfleets_survey(j);z++)    /// survey index 
                  {
                  ///// overlap survey calcs/////////////////////////////////////////////////////////////////////////////////           
                   true_survey_fleet_overlap_age(p,j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM_overlap_region(p,j,y,a,r)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tsurvey(j,r))*q_survey(j,r,z);
                   true_survey_fleet_overlap_age_bio(p,j,r,y,z,a)=true_survey_fleet_overlap_age(p,j,r,y,z,a)*weight_population(p,r,y,a);

                   if(natal_homing_switch==0) //for non-natal homing pop structures
                    {
                     true_survey_fleet_age(j,r,y,z,a)=survey_selectivity(j,r,y,a,z)*abundance_at_age_AM(j,r,y,a)*mfexp(-(M(j,r,y,a)+F(j,r,y,a))*tsurvey(j,r))*q_survey(j,r,z);
                     true_survey_fleet_age_bio(j,r,y,z,a)=true_survey_fleet_age(j,r,y,z,a)*weight_population(j,r,y,a);                  
                    }
                
                   // calculate true survey abundance for tagging simulation **dh**
                   true_survey_fleet_age_temp(j,r,y,a,z)=true_survey_fleet_age(j,r,y,z,a);
                  }
                  
                 true_survey_region_abundance(j,y,r,a)=sum(true_survey_fleet_age_temp(j,r,y,a));
                 true_survey_region_abundance_temp(j,y,a,r)=true_survey_region_abundance(j,y,r,a);
                 true_survey_population_abundance(y,j,a)=sum(true_survey_region_abundance_temp(j,y,a));
                 true_survey_population_abundance_temp(y,a,j)=true_survey_population_abundance(y,j,a);
                 true_survey_total_abundance(y,a)=sum(true_survey_population_abundance_temp(y,a));
                }
                
           } //end age loop

             ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
             ////////////////////// NATAL HOMING CALCS ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
             ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

               SSB_region_overlap(p,j,r,y)=sum(SSB_region_temp_overlap(p,j,r,y));

                SSB_overlap_natal=0;
               if(natal_homing_switch==1 && spawn_return_switch==1) //if natal homing and return spawning migration occurs each year, then add SSB of fish not in natal pop, but that do instantaneous spawn migration
                {
                 for(int k=1;k<=npops;k++)
                  {
                   for (int n=1;n<=nregions(k);n++)
                    {
                     if(p==k && j==k)
                      {
                       SSB_overlap_natal(k,n)=0;  // already account for all fish already in natal population in calc below
                      }   
                     if(p==j && j!=k) //if j is natal pop, but not currently in j then prop of SSB does instantaneous spawn migration
                      {
                       SSB_overlap_natal(k,n)=spawn_return_prob(p)*SSB_region_overlap(p,k,n,y);
                      } 
                    }
                  } 
                 if(p==j)  //here just add SSB that 'returns' to natal pop SSB, SRR uses all SSB, but then overlap recruit calcs only keep recruits that were from SSB in natal pop
                  {
                   SSB_region_overlap(p,j,r,y)=SSB_region_overlap(p,j,r,y)+(sum(SSB_overlap_natal)/nregions(p));  //reutrning SSB is split evenly across natal regions
                  }
                }
                   
                SSB_population_temp_overlap(p,j,y,r)=SSB_region_overlap(p,j,r,y); 
                SSB_population_overlap(p,j,y)=sum(SSB_population_temp_overlap(p,j,y));
                SSB_natal_overlap_temp(p,y,j)=SSB_population_overlap(p,j,y);
                SSB_natal_overlap(p,y)=sum(SSB_natal_overlap_temp(p,y));

              if(natal_homing_switch>0) //natal homing occurs
               {
                if(p==j)
                 {
                   SSB_region(j,r,y)=SSB_region_overlap(p,j,r,y);  //if natal homing only account for SSB that is in its natal populations area, don't sum across natal populations
                  }
               }
            
                Bratio_population_overlap(p,j,y)=SSB_population_overlap(p,j,y)/SSB_zero(p);
                Bratio_natal_overlap(p,y)=SSB_natal_overlap(p,y)/SSB_zero(p);
                
                yield_region_overlap(p,j,r,y)=sum(yield_region_temp_overlap(p,j,r,y));
                yield_population_overlap(p,j,y)=sum(yield_population_temp_overlap(p,j,y));
                yield_natal_temp_overlap(p,y,j)=yield_population_overlap(p,j,y);
                yield_natal_overlap(p,y)=sum(yield_natal_temp_overlap(p,y));
                
                harvest_rate_region_bio_overlap(p,j,r,y)=yield_region_overlap(p,j,r,y)/biomass_AM_overlap_region(p,j,r,y);
                harvest_rate_population_bio_overlap(p,j,y)=yield_population_overlap(p,j,y)/biomass_population_overlap(p,j,y);
                harvest_rate_natal_bio_overlap(p,y)=yield_natal_overlap(p,y)/biomass_natal_overlap(p,y);
               
                depletion_region_overlap(p,j,r,y)=biomass_AM_overlap_region(p,j,r,y)/biomass_AM_overlap_region(p,j,r,1); //relative to year 1 biomass
                depletion_population_overlap(p,j,y)=biomass_population_overlap(p,j,y)/biomass_population_overlap(p,j,1);
                depletion_natal_overlap(p,y)=biomass_natal_overlap(p,y)/biomass_natal_overlap(p,1);
                
               for (int z=1;z<=nfleets(j);z++)
                {
                 harvest_rate_region_fleet_bio_overlap(p,j,r,z,y)=yield_region_fleet_overlap(p,j,r,z,y)/biomass_AM_overlap_region(p,j,r,y);

                 yield_region_fleet_overlap(p,j,r,z,y)=sum(yield_region_fleet_temp_overlap(p,j,r,z,y));
                 yield_fleet(j,r,y,z)=sum(yield_fleet_temp(j,r,y,z));

                //YIELD observation error
                 OBS_yield_region_fleet_overlap(p,j,r,y,z)=yield_region_fleet_overlap(p,j,r,z,y)*mfexp(yield_RN_overlap(p,j,r,y,z)*sigma_catch_overlap(p,j,r,z)-.5*square(sigma_catch_overlap(p,j,r,z)));
                 OBS_yield_fleet_temp(j,r,y,z,p)=OBS_yield_region_fleet_overlap(p,j,r,y,z);
                 
                 if(natal_homing_switch==0) //no natal homing
                  {
                   OBS_yield_fleet(j,r,y,z)=yield_fleet(j,r,y,z)*mfexp(yield_RN(j,r,y,z)*sigma_catch(j,r,z)-.5*square(sigma_catch(j,r,z)));
                  }
                 if(natal_homing_switch==1) //natal homing
                  {
                   OBS_yield_fleet(j,r,y,z)=sum(OBS_yield_fleet_temp(j,r,y,z));  
                  }
                } //end fleets loop

                OBS_yield_region_overlap(p,j,y,r)=sum(OBS_yield_region_fleet_overlap(p,j,r,y));
                OBS_yield_population_overlap(p,y,j)=sum(OBS_yield_region_overlap(p,j,y));
                OBS_yield_natal_overlap(y,p)=sum(OBS_yield_population_overlap(p,y));
                OBS_yield_total_overlap(y)=sum(OBS_yield_natal_overlap(y));

           ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
           ////////////////////// NON-NATAL HOMING CALCS ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
           ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                SSB_region(j,r,y)=sum(SSB_region_temp(j,r,y));
                SSB_population_temp(j,y,r)=SSB_region(j,r,y); 
                SSB_population(j,y)=sum(SSB_population_temp(j,y)); 
                SSB_total_temp(y,j)=SSB_population(j,y);
                SSB_total(y)=sum(SSB_total_temp(y));

                Bratio_population(j,y)=SSB_population(j,y)/SSB_zero(j);
                Bratio_total(y)=SSB_total(y)/sum(SSB_zero);
                

                yield_region(j,r,y)=sum(yield_region_temp(j,r,y));
                yield_population(j,y)=sum(yield_population_temp(j,y));
                yield_total_temp(y,j)=yield_population(j,y);
                yield_total(y)=sum(yield_total_temp(y));
  
                harvest_rate_region_bio(j,r,y)=yield_region(j,r,y)/biomass_AM(j,r,y);
                harvest_rate_population_bio(j,y)=yield_population(j,y)/biomass_population(j,y);
                harvest_rate_total_bio(y)=yield_total(y)/biomass_total(y);
                
                depletion_region(j,r,y)=biomass_AM(j,r,y)/biomass_AM(j,r,1); //relative to year 1 biomass
                depletion_population(j,y)=biomass_population(j,y)/biomass_population(j,1);
                depletion_total(y)=biomass_total(y)/biomass_total(1);

                OBS_yield_region(j,y,r)=sum(OBS_yield_fleet(j,r,y));
                OBS_yield_population(y,j)=sum(OBS_yield_region(j,y));
                OBS_yield_total(y)=sum(OBS_yield_population(y));
                
                //apportion variables
                apport_yield_region(j,r,y)=OBS_yield_region(j,y,r)/OBS_yield_population(y,j);

               if(tsurvey(j,r)>0) //if survey not at beginning of year do time adjustments for mortality at time of survey
                {               
                 for (int z=1;z<=nfleets_survey(j);z++)    /// survey index 
                  {
                  ///// overlap survey calcs/////////////////////////////////////////////////////////////////////////////////
                   true_survey_fleet_bio_overlap(p,j,r,y,z)=sum(true_survey_fleet_overlap_age_bio(p,j,r,y,z));  
                   true_survey_fleet_bio_overlap_temp(j,r,y,z,p)=true_survey_fleet_bio_overlap(p,j,r,y,z);

                  //add observation error for observed survey index //USE THIS DATA IF ASSUME HAVE GENETIC DATA TO DIFFERENTIATE SURVEY DATA BY NATAL POPULATION
                   OBS_survey_fleet_bio_overlap(p,j,r,y,z)=true_survey_fleet_bio_overlap(p,j,r,y,z)*mfexp(survey_RN_overlap(p,j,r,y,z)*sigma_survey_overlap(p,j,r,z)-.5*square(sigma_survey_overlap(p,j,r,z)));
                   OBS_survey_fleet_bio_temp(j,r,y,z,p)=OBS_survey_fleet_bio_overlap(p,j,r,y,z);
                   
                   if(natal_homing_switch==0) //for non-natal homing pop structures
                    {
                     true_survey_fleet_bio(j,r,y,z)=sum(true_survey_fleet_age_bio(j,r,y,z));
                     OBS_survey_fleet_bio(j,r,y,z)=true_survey_fleet_bio(j,r,y,z)*mfexp(survey_RN(j,r,y,z)*sigma_survey(j,r,z)-.5*square(sigma_survey(j,r,z)));
                    }
                   if(natal_homing_switch==1) // for natal homing structures OBS survey is summed over natal populations in a given pop unit j //WOULD FIT THIS DATA IF ASSUME NO GENETIC DATA TO DIFFERENTIATE THEN FIT BY NATAL POP
                    {
                     true_survey_fleet_bio(j,r,y,z)=sum(true_survey_fleet_bio_overlap_temp(j,r,y,z));
                     OBS_survey_fleet_bio(j,r,y,z)=sum(OBS_survey_fleet_bio_temp(j,r,y,z));  
                    }
                  } //end fleets loop
                    
                 true_survey_region_bio_overlap(p,j,y,r)=sum(true_survey_fleet_bio_overlap(p,j,r,y));               
                 true_survey_population_bio_overlap(p,y,j)=sum(true_survey_region_bio_overlap(p,j,y));               
                 true_survey_natal_bio_overlap(y,p)=sum(true_survey_population_bio_overlap(p,y));               
                 true_survey_total_bio_overlap(y)=sum(true_survey_natal_bio_overlap(y));
                 OBS_survey_region_bio_overlap(p,j,y,r)=sum(OBS_survey_fleet_bio_overlap(p,j,r,y));
                 OBS_survey_population_bio_overlap(p,y,j)=sum(OBS_survey_region_bio_overlap(p,j,y));
                 OBS_survey_natal_bio_overlap(y,p)=sum(OBS_survey_population_bio_overlap(p,y));
                 OBS_survey_total_bio_overlap(y)=sum(OBS_survey_natal_bio_overlap(y));
                  
                 true_survey_region_bio(j,y,r)=sum(true_survey_fleet_bio(j,r,y));
                 true_survey_population_bio(y,j)=sum(true_survey_region_bio(j,y));
                 true_survey_total_bio(y)=sum(true_survey_population_bio(y));
                 OBS_survey_region_bio(j,y,r)=sum(OBS_survey_fleet_bio(j,r,y));
                 OBS_survey_population_bio(y,j)=sum(OBS_survey_region_bio(j,y));
                 OBS_survey_total_bio(y)=sum(OBS_survey_population_bio(y));
                
                 apport_region_survey_biomass(j,r,y)=OBS_survey_region_bio(j,y,r)/OBS_survey_population_bio(y,j);
         
                } //tsurvey>0

         }
       }
     }
   } //end year loop

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////now adjusting natal homing ssb for spawning migration scenarios where fish left non-natal area
       for (int p=1;p<=npops;p++)
        {
         for (int j=1;j<=npops;j++)
          {
           for (int r=1;r<=nregions(p);r++)
            {
             for (int t=1;t<=nregions(j);t++)
              {
               for (int y=1;y<=nyrs;y++)
                {
                   SSB_overlap_natal=0;
                  if(natal_homing_switch==1 && spawn_return_switch==1)
                   {
                    if(p!=j)  //update SSB that doesn't spawn (ie doesn't return to natal population)  //not setting SSB to 0, but just updating so SSB in non-natal areas is discounted for fish that actually migrated back to natal area to spawn
                     {
                      SSB_region_overlap(p,j,r,y)=(1-spawn_return_prob(p))*SSB_region_overlap(p,j,r,y);
                     }
                      SSB_population_temp_overlap(p,j,y,r)=SSB_region_overlap(p,j,r,y); 
                      SSB_population_overlap(p,j,y)=sum(SSB_population_temp_overlap(p,j,y));
                      SSB_natal_overlap_temp(p,y,j)=SSB_population_overlap(p,j,y);
                      SSB_natal_overlap(p,y)=sum(SSB_natal_overlap_temp(p,y));
                      Bratio_population_overlap(p,j,y)=SSB_population_overlap(p,j,y)/SSB_zero(p);
                      Bratio_natal_overlap(p,y)=SSB_natal_overlap(p,y)/SSB_zero(p);
                      Bratio_population(j,y)=SSB_population(j,y)/SSB_zero(j);
                      Bratio_total(y)=SSB_total(y)/sum(SSB_zero);
                    }
                   }
                  }
                 }
                }
               }

  for (int j=1;j<=npops;j++)
   {
    for (int r=1;r<=nregions(j);r++)
     {
      for(int y=1;y<=nyrs;y++)
       {
          for (int k=1;k<=npops;k++)
           {
            for (int n=1;n<=nregions(k);n++)
             {
              T_year(j,r,y,k,n)=T(j,r,y,4,k,n);     //because can't report out 6d arrays properly report out a yearly T for certain graphics       
             } 
            }
           }
          }
         }

 ///manually calculating the true distribution of initial abundance for use in EM (so don't have to estimate)
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
    {
     for (int r=1;r<=nregions(j);r++)
      {
        init_abund_reg_temp(p,j,r)=sum(init_abund(p,j,r));
         for (int a=1;a<=nages;a++)
          {
           init_abund_age_temp(p,a,j,r)=init_abund(p,j,r,a);
          }
      }
      init_abund_pop_temp(p,j)=sum(init_abund_reg_temp(p,j));
       for (int a=1;a<=nages;a++)
        {
         init_abund_reg_age_temp(p,a,j)=sum(init_abund_age_temp(p,a,j));
        }
     }
    }
    
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
    {
     for (int r=1;r<=nregions(j);r++)
      {
       frac_natal_true(p,j,r)=sum(init_abund(p,j,r))/sum(init_abund_pop_temp(p));
        for (int a=1;a<=nages;a++)
         {
          frac_natal_true_age(p,j,r,a)=init_abund(p,j,r,a)/sum(init_abund_reg_age_temp(p,a));
         }
      }
     }
    }

FUNCTION get_rand_survey_CAA_prop
 random_number_generator myrand_survey_age(myseed_survey_age);
 
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int z=1;z<=nfleets_survey(j);z++)
         {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
            for (int a=1;a<=nages;a++)
             {
              survey_at_age_region_fleet_overlap_prop(p,j,r,z,y,a)=true_survey_fleet_overlap_age(p,j,r,y,z,a)/sum(true_survey_fleet_overlap_age(p,j,r,y,z));              
              survey_at_age_fleet_prop(j,r,y,z,a)=true_survey_fleet_age(j,r,y,z,a)/sum(true_survey_fleet_age(j,r,y,z));
             }
            }
           }
          }
         }
        }

  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int z=1;z<=nfleets_survey(j);z++)
         {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
            if(use_stock_comp_info_survey==0)
             {
              rand_SIM_survey_prop_temp(j,r,y,z).fill_multinomial(myrand_survey_age,value(survey_at_age_fleet_prop(j,r,y,z)));
             }
            if(use_stock_comp_info_survey==1)
             {
              rand_SIM_survey_prop_temp_overlap(p,j,r,y,z).fill_multinomial(myrand_survey_age,value(survey_at_age_region_fleet_overlap_prop(p,j,r,z,y)));
             }
           }
         }
       }
     }
   }

  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int z=1;z<=nfleets_survey(j);z++)
         {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
             rand_SIM_survey_prop_temp2=0;
            if(use_stock_comp_info_survey==0)
             {
               for(int n=1;n<=SIM_nsurvey(j,r,z);n++) 
                {
                 rand_SIM_survey_prop_temp2(value(rand_SIM_survey_prop_temp(j,r,y,z,n)))+= 1.0;
                }
               SIM_survey_prop(j,r,z,y)=rand_SIM_survey_prop_temp2;
             }
            if(use_stock_comp_info_survey==1)
             {              
               for(int n=1;n<=SIM_nsurvey_overlap(p,j,r,z);n++) /// look into changing this so can have ncatch change by year (ie different sample sizes for beginning and end of timeseries)
                {
                 rand_SIM_survey_prop_temp2(value(rand_SIM_survey_prop_temp_overlap(p,j,r,y,z,n)))+= 1.0;
                }
               SIM_survey_prop_overlap(p,j,r,z,y)=rand_SIM_survey_prop_temp2;
             }
             
        for(int a=1;a<=nages;a++)
         {
          if(use_stock_comp_info_survey==0)
           {
            OBS_survey_prop(j,r,y,z,a)=SIM_survey_prop(j,r,z,y,a)/SIM_nsurvey(j,r,z);
           }
          if(use_stock_comp_info_survey==1)
           {
            OBS_survey_prop_overlap(p,j,r,y,z,a)=SIM_survey_prop_overlap(p,j,r,z,y,a)/SIM_nsurvey_overlap(p,j,r,z);
           }
         }
       }
      }
     }
    }
   }

FUNCTION get_rand_CAA_prop
 random_number_generator myrand_catch_age(myseed_catch_age);
 
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int z=1;z<=nfleets(j);z++)
         {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
            for (int a=1;a<=nages;a++)
             {
                 catch_at_age_region_fleet_overlap_prop(p,j,r,z,y,a)=catch_at_age_region_fleet_overlap(p,j,r,z,y,a)/sum(catch_at_age_region_fleet_overlap(p,j,r,z,y));              
                 catch_at_age_region_overlap_prop(p,j,r,y,a)=catch_at_age_region_overlap(p,j,r,y,a)/sum(catch_at_age_region_overlap(p,j,r,y));
                 catch_at_age_population_overlap_prop(p,j,y,a)=catch_at_age_population_overlap(p,j,y,a)/sum(catch_at_age_population_overlap(p,j,y));
                 catch_at_age_natal_overlap_prop(p,y,a)=catch_at_age_natal_overlap(p,y,a)/sum(catch_at_age_natal_overlap(p,y));

                 catch_at_age_fleet_prop_temp(j,r,y,z,a)=catch_at_age_fleet(j,r,y,a,z);
                 catch_at_age_fleet_prop(j,r,y,z,a)=catch_at_age_fleet_prop_temp(j,r,y,z,a)/sum(catch_at_age_fleet_prop_temp(j,r,y,z));
                 catch_at_age_region_prop(j,r,y,a)=catch_at_age_region(j,r,y,a)/sum(catch_at_age_region(j,r,y));
                 catch_at_age_population_prop(j,y,a)=catch_at_age_population(j,y,a)/sum(catch_at_age_population(j,y));
                 catch_at_age_total_prop(y,a)=catch_at_age_total(y,a)/sum(catch_at_age_total(y));
              }
            }
           }
          }
         }
        }
        
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int z=1;z<=nfleets(j);z++)
         {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
            if(use_stock_comp_info_catch==0)
             {
              rand_SIM_catch_prop_temp(j,r,y,z).fill_multinomial(myrand_catch_age,value(catch_at_age_fleet_prop(j,r,y,z)));
             }
            if(use_stock_comp_info_catch==1)
             {
              rand_SIM_catch_prop_temp_overlap(p,j,r,y,z).fill_multinomial(myrand_catch_age,value(catch_at_age_region_fleet_overlap_prop(p,j,r,z,y)));
             }
           }
         }
       }
     }
   }
            
  for (int p=1;p<=npops;p++)
   {
    for (int j=1;j<=npops;j++)
     {
      for (int r=1;r<=nregions(j);r++)
       {
        for (int z=1;z<=nfleets(j);z++)
         {
          for (int y=1;y<=nyrs;y++) //need to alter to fit number of years of catch data
           {
             rand_SIM_catch_prop_temp2=0;

            if(use_stock_comp_info_catch==0)
             {
               for(int n=1;n<=SIM_ncatch(j,r,z);n++) /// look into changing this so can have ncatch change by year (ie different sample sizes for beginning and end of timeseries)
                {
                 rand_SIM_catch_prop_temp2(value(rand_SIM_catch_prop_temp(j,r,y,z,n)))+= 1.0;
                }
               SIM_catch_prop(j,r,z,y)=rand_SIM_catch_prop_temp2;
             }
            if(use_stock_comp_info_catch==1)
             {              
               for(int n=1;n<=SIM_ncatch_overlap(p,j,r,z);n++) /// look into changing this so can have ncatch change by year (ie different sample sizes for beginning and end of timeseries)
                {
                 rand_SIM_catch_prop_temp2(value(rand_SIM_catch_prop_temp_overlap(p,j,r,y,z,n)))+= 1.0;
                }
               SIM_catch_prop_overlap(p,j,r,z,y)=rand_SIM_catch_prop_temp2;
             }
             
        for(int a=1;a<=nages;a++)
         {
          if(use_stock_comp_info_catch==0)
           {
            OBS_catch_prop(j,r,y,z,a)=SIM_catch_prop(j,r,z,y,a)/SIM_ncatch(j,r,z);
           }
          if(use_stock_comp_info_catch==1)
           {
            OBS_catch_prop_overlap(p,j,r,y,z,a)=SIM_catch_prop_overlap(p,j,r,z,y,a)/SIM_ncatch_overlap(p,j,r,z);
           }
         }
       }
      }
     }
    }
   }

FUNCTION get_tag_recaptures
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  // In some calcs (mortality and movement) need to account for true age (age+time at large),
  // because multiple cohorts this means need to factor in release year
  // so true age becomes age+year-release year
  // using subscript notation this equals= (a+y-xx)
  // similarly because account for release age, do not need to worry about plus group calcs as carry recaptures out to max_life_tags and never assume plus group (just use plus group mortality and movement values in calcs where a>=max_age)
  /////////////////////////////////////////////////////////////////////////////
  //reporting rate is assumed to be function of release event and recap location (not a function of recap year...could expand to this, but not priority at momement)
  ///////////////////////////////////////////////////////////////////////////////////
 if(do_tag==1)
  {
   for(int x=1; x<=nyrs_release; x++)
    {
     xx=yrs_releases(x);
     ntags_total_temp=0;
      for (int i=1;i<=npops;i++)
       {
         for (int n=1;n<=nregions(i);n++)
          {
           if(prob_tag_RN(i,n,x)>=opport_tag_prob) //set prob tag (for opportunistic tagging) to 1 or 0 (occurs/doesn't occur) based on uniform random number, if opport_Tag_prop==0.5 then 50/50 chance of tagging
           {
            prob_tag(i,n,x)=1.0;
           }
            if(prob_tag_RN(i,n,x)<opport_tag_prob)
           {
            prob_tag(i,n,x)=0.0;
           }
           if(prob_tag_year_RN(x)>=opport_tag_prob_year) //set prob tag (for opportunistic tagging) to 1 or 0 (occurs/doesn't occur) based on uniform random number, if opport_Tag_prop==0.5 then 50/50 chance of tagging
           {
            prob_tag_year(x)=1.0;
           }
            if(prob_tag_year_RN(x)<opport_tag_prob_year)
           {
            prob_tag_year(x)=0.0;
           }
         if(number_tags_switch==(-2)) //input total tags is distributed evenly across regions
          {
            ntags_region(i,n,x)=input_total_tags(x)/(sum(nregions)); //sum abundance across ages and all areas/pops
          }
         if(number_tags_switch==(-1)) //input total tags is distributed across regions according to survey abund
          {
            ntags_region(i,n,x)=input_total_tags(x)*(sum(true_survey_population_abundance(xx,i))/sum(true_survey_total_abundance(xx)))*(sum(true_survey_region_abundance(i,xx,n))/sum(true_survey_population_abundance(xx,i))); //sum abundance across ages and all areas/pops; //sum abundance across ages by region
          }
         if(number_tags_switch==0) // tags by area is a fixed number based on input # tags
          {
            ntags_region(i,n,x)=input_ntags(i,n,x); 
          }
         if(number_tags_switch==1) //total tags is a fraction of total abundance and tags assigned evely to each region
          {
            ntags_region(i,n,x)=frac_abund_tagged(x)*sum(abundance_total(xx))/(sum(nregions)); //sum abundance across ages and all areas/pops
          }
         if(number_tags_switch==2) //total tags is a fraction of pop/reg abundance
          {
            ntags_region(i,n,x)=frac_abund_tagged(x)*sum(abundance_at_age_AM(i,n,xx)); //sum abundance across ages by region
          }
         if(number_tags_switch==3) //total tags is a fraction of total abundance then assigned to region based on regional survey abundance
          {
            ntags_region(i,n,x)=frac_abund_tagged(x)*sum(abundance_total(xx))*(sum(true_survey_population_abundance(xx,i))/sum(true_survey_total_abundance(xx)))*(sum(true_survey_region_abundance(i,xx,n))/sum(true_survey_population_abundance(xx,i))); //sum abundance across ages and all areas/pops
          }
         if(number_tags_switch==4) // tags by area is a random number based on uniform dist with mean==((max+min)/2)
          {
           //number of tags is random number, but all areas receive tag releases (unless max and min set to 0 then no tags)
           ntags_region(i,n,x)=ntags_RN(i,n,x)*(max_tags(i,n)-min_tags(i,n))+min_tags(i,n); //scale uniform number ([0,1)) to interval that want
          }
         if(number_tags_switch==5) //number of tags is input and distributed proportional to survey abundance AND probability of tagging in a given AREA is random number (i.e., meant to simulate completely opportunistic tagging)
          {
           ntags_region(i,n,x)=input_total_tags(x)*(sum(true_survey_population_abundance(xx,i))/sum(true_survey_total_abundance(xx)))*(sum(true_survey_region_abundance(i,xx,n))/sum(true_survey_population_abundance(xx,i)))*prob_tag(i,n,x); //if prob_tag==0 then no tagging in that release year/area
          }
         if(number_tags_switch==6) //number of tags is input and distributed proportional to survey abundance AND probability of tagging in a given YEAR is random number
          {
           ntags_region(i,n,x)=input_total_tags(x)*(sum(true_survey_population_abundance(xx,i))/sum(true_survey_total_abundance(xx)))*(sum(true_survey_region_abundance(i,xx,n))/sum(true_survey_population_abundance(xx,i)))*prob_tag_year(x); //if prob_tag==0 then no tagging in that release year/area
          }
         if(number_tags_switch==7) //number of tags is random number AND probability of tagging in a given YEAR is random number AND probability of tagging in a given AREA is a random number
          {
           ntags_region(i,n,x)=(ntags_RN(i,n,x)*(max_tags(i,n)-min_tags(i,n))+min_tags(i,n))*prob_tag_year(x)*prob_tag(i,n,x); //if prob_tag==0 then no tagging in that release year/area
          }
         ntags_total_temp(i,n)=ntags_region(i,n,x);
          }
        }
         ntags_total(x)=sum(ntags_total_temp);
       }

 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
       for (int a=1;a<=nages;a++) //release age 
        {
         if(age_dist_tags==0)  //distribute tags evenly across ages
          {
           ntags(i,n,x,a)=ntags_region(i,n,x)/(nages);
          }
         if(age_dist_tags==1)  //distribute tags across ages based on regional survey proportions at age (i.e., based on the proportions you would expect to catch fish in a survey...use primary survey in a region)
          {
           ntags(i,n,x,a)=ntags_region(i,n,x)*survey_at_age_fleet_prop(i,n,xx,1,a);
          }
         if(age_dist_tags==2)  //distribute tags across ages based on regional catch proportions at age (i.e., based on the proportions you would expect to catch fish in the fishery...use average across fleets in a region)
          {
           ntags(i,n,x,a)=ntags_region(i,n,x)*catch_at_age_region_prop(i,n,xx,a);
          }
        }
       }
      }
     }


////////////////////////////////////////////////////////////////////////////////////////////
 // for the EM .dat when have mismatches
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age 
        {
        ntags_temp1(i,x,a,n)=ntags(i,n,x,a);
        ntags_population(i,x,a)=sum(ntags_temp1(i,x,a));
        ntags_pan_temp(x,a,i)=ntags_population(i,x,a);
        ntags_pan(x,a)=sum(ntags_pan_temp(x,a));
        }
       }
      }
     }
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////     
//tag dynamics
//////////////////////////////////////////////////

 //assume tags released in natal population
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age //because accounting for release age, don't need to account for recap age, just adjust mortality and T, to use plus group value if recapture age exceeds max age
        {
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
              if(y==1) //year of release for a cohort
               {
                if(sim_tag_mixing_switch==0) //assume complete mixing of tagged and untagged fish
                 {
                  tags_avail(i,n,x,a,y,j,r)=ntags(i,n,x,a)*T(i,n,xx,a,j,r); 
                  recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,xx,a)*(1.-mfexp(-(F(j,r,xx,a)+M(j,r,xx,a))))/(F(j,r,xx,a)+(M(j,r,xx,a)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)

               if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                {
                 tags_avail(i,n,x,a,y,j,r)=ntags(i,n,x,a)*T(i,n,xx,tag_age_sel,j,r); 
                 recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,xx,tag_age_sel)*(1.-mfexp(-(F(j,r,xx,tag_age_sel)+M(j,r,xx,tag_age_sel))))/(F(j,r,xx,tag_age_sel)+(M(j,r,xx,tag_age_sel)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                                
                }
                
                 }
                if(sim_tag_mixing_switch==1) //assume incomplete mixing of tagged and untagged fish
                 {
                 if(sim_tag_mixing_T_switch==0) //assume incomplete mixing of tagged and untagged fish only for F
                  {
                   tags_avail(i,n,x,a,y,j,r)=ntags(i,n,x,a)*T(i,n,xx,a,j,r); 
                  }
                 if(sim_tag_mixing_T_switch==1) //assume incomplete mixing of tagged and untagged fish 
                  {                 
                   if(i==j && n==r)
                    {
                      T_tag(i,n,x,a,j,r)=T_tag_res(x);
                    }
                   if(i!=j || n!=r)
                    {
                     if(move_switch==8)
                     {
                      T_tag(i,n,x,a,j,r)=(1-T_tag_res(x))*rel_bio(i,n,xx,a,j,r);
                     }
                     if(move_switch!=8)
                     {
                      T_tag(i,n,x,a,j,r)=(1-T_tag_res(x))/(sum(nregions)-1); //split movement evenly across remaining regions
                     }
                    }
                      if(T_tag(i,n,x,a,j,r)>1)
                       {
                        T_tag(i,n,x,a,j,r)=1;
                       }
                      if(T_tag(i,n,x,a,j,r)<0)
                       {
                        T_tag(i,n,x,a,j,r)=0;
                       }
                    tags_avail(i,n,x,a,y,j,r)=ntags(i,n,x,a)*T_tag(i,n,x,a,j,r);                        
                   }
                  F_tag(j,r,x,a)=F_tag_scalar(x)*F(j,r,xx,a);
                  recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*(F_tag_scalar(x)*F(j,r,xx,a))*(1.-mfexp(-((F_tag_scalar(x)*F(j,r,xx,a))+M(j,r,xx,a))))/((F_tag_scalar(x)*F(j,r,xx,a))+(M(j,r,xx,a)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)
                 }
               }
            if(y==2) // must account for the maximum age so use min function to ensure that not exceeding the max age
             {
             if(sim_tag_mixing_switch==0) //assume complete mixing of tagged and untagged fish
              {
               tags_avail_temp=0;
               for(int p=1;p<=npops;p++)
               {
                for (int s=1;s<=nregions(p);s++)
                {
                 if(natal_homing_switch==0) //if no natal homing
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),min((a+y),nages),j,r)*mfexp(-(F(p,s,(xx+y-2),min(((a+y)-1),nages))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)

                   if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                    {
                     tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),tag_age_sel,j,r)*mfexp(-(F(p,s,(xx+y-2),tag_age_sel)+(M(p,s,(xx+y-2),tag_age_sel)))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)
                    }
                 }                        
                //#####################################################################################################
                //  TRUE NATAL HOMING  T(n,x,a,y,j) becomes T(i,x,a,y,j) because need to maintain your natal origin
                //  movement values so T doesn't depend on current population only origin population and destination population
                //########################################################################################################              
                 if(natal_homing_switch==1) //if natal homing 
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(i,n,(xx+y-1),min((a+y),nages),j,r)*mfexp(-(F(p,s,(xx+y-2),min(((a+y)-1),nages))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)

                   if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                    {
                     tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),tag_age_sel,j,r)*mfexp(-(F(p,s,(xx+y-2),tag_age_sel)+(M(p,s,(xx+y-2),tag_age_sel)))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)
                    }
                 }               
                }
               }
                 tags_avail(i,n,x,a,y,j,r)=sum(tags_avail_temp); //sum across all pops/regs of tags that moved into pop j reg r
                 recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,(xx+y-1),min((a+y),nages))*(1.-mfexp(-(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))))))/(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)   

                   if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                    {
                     recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,(xx+y-1),tag_age_sel)*(1.-mfexp(-(F(j,r,(xx+y-1),tag_age_sel)+(M(j,r,(xx+y-1),tag_age_sel)))))/(F(j,r,(xx+y-1),tag_age_sel)+(M(j,r,(xx+y-1),tag_age_sel)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                 
                    }    
               }
             if(sim_tag_mixing_switch==1) //assume incomplete mixing of tagged and untagged fish
              {
               tags_avail_temp=0;
               for(int p=1;p<=npops;p++)
               {
                for (int s=1;s<=nregions(p);s++)
                {
                 if(natal_homing_switch==0) //if no natal homing
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),min((a+y),nages),j,r)*mfexp(-((F_tag_scalar(x)*F(p,s,(xx+y-2),min(((a+y)-1),nages)))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)
                 }                        
                //#####################################################################################################
                //  TRUE NATAL HOMING  T(n,x,a,y,j) becomes T(i,x,a,y,j) because need to maintain your natal origin
                //  movement values so T doesn't depend on current population only origin population and destination population
                //########################################################################################################              
                 if(natal_homing_switch==1) //if natal homing 
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(i,n,(xx+y-1),min((a+y),nages),j,r)*mfexp(-((F_tag_scalar(x)*F(p,s,(xx+y-2),min(((a+y)-1),nages)))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)             
                 }               
                }
               }
                 tags_avail(i,n,x,a,y,j,r)=sum(tags_avail_temp); //sum across all pops/regs of tags that moved into pop j reg r
                 recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,(xx+y-1),min((a+y),nages))*(1.-mfexp(-(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))))))/(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                 
               }
              }
              if(y>2) // must account for the maximum age so use min function to ensure that not exceeding the max age
              {
               tags_avail_temp=0;
               for(int p=1;p<=npops;p++)
               {
                for (int s=1;s<=nregions(p);s++)
                {
                 if(natal_homing_switch==0) //if no natal homing
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),min((a+y),nages),j,r)*mfexp(-(F(p,s,(xx+y-2),min(((a+y)-1),nages))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)

                   if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                    {
                     tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),tag_age_sel,j,r)*mfexp(-(F(p,s,(xx+y-2),tag_age_sel)+(M(p,s,(xx+y-2),tag_age_sel)))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)
                    }
                 }                        
                //#####################################################################################################
                //  TRUE NATAL HOMING  T(n,x,a,y,j) becomes T(i,x,a,y,j) because need to maintain your natal origin
                //  movement values so T doesn't depend on current population only origin population and destination population
                //########################################################################################################              
                 if(natal_homing_switch==1) //if natal homing 
                 {
                  tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(i,n,(xx+y-1),min((a+y),nages),j,r)*mfexp(-(F(p,s,(xx+y-2),min(((a+y)-1),nages))+(M(p,s,(xx+y-2),min(((a+y)-1),nages))))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)

                   if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                    {
                     tags_avail_temp(p,s)=tags_avail(i,n,x,a,y-1,p,s)*T(p,s,(xx+y-1),tag_age_sel,j,r)*mfexp(-(F(p,s,(xx+y-2),tag_age_sel)+(M(p,s,(xx+y-2),tag_age_sel)))); //tags_temp holds all tags moving into population j,region r; min function takes min of true age and max age allowed (plus group)
                    }
                 }               
                }
               }
                 tags_avail(i,n,x,a,y,j,r)=sum(tags_avail_temp); //sum across all pops/regs of tags that moved into pop j reg r
                 recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,(xx+y-1),min((a+y),nages))*(1.-mfexp(-(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))))))/(F(j,r,(xx+y-1),min((a+y),nages))+(M(j,r,(xx+y-1),min((a+y),nages))));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)

                   if(tag_fit_ages_switch_OM==1 && tag_fit_ages_switch==1) //need to determine age of full selectivity for tagging data if assuming don't know age of tags
                    {
                     recaps(i,n,x,a,y,j,r)=report_rate_TRUE(j,x,r)*tags_avail(i,n,x,a,y,j,r)*F(j,r,(xx+y-1),tag_age_sel)*(1.-mfexp(-(F(j,r,(xx+y-1),tag_age_sel)+(M(j,r,(xx+y-1),tag_age_sel)))))/(F(j,r,(xx+y-1),tag_age_sel)+(M(j,r,(xx+y-1),tag_age_sel)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                 
                    }  
               }
             }
            }
           }
          }
         }
        }
       }
  
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age //because accounting for release age, don't need to account for recap age, just adjust mortality and T, to use plus group value if recapture age exceeds max age
        {
        total_recap_temp.initialize();
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
              total_recap_temp(j,r,y)=recaps(i,n,x,a,y,j,r);
             }
            }
           }
             total_rec(i,n,x,a)=sum(total_recap_temp);
             not_rec(i,n,x,a)=ntags(i,n,x,a)-total_rec(i,n,x,a);  //for ntags  at a given age all entries represent all tags released so can just use any of the entries (hence the i,x,a,1 subscripts)
           }
          }
         }
        }
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age //because accounting for release age, don't need to account for recap age, just adjust mortality and T, to use plus group value if recapture age exceeds max age
        {
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
              if(ntags(i,n,x,a)>0)
               {
                tag_prop(i,n,x,a,y,j,r)=recaps(i,n,x,a,y,j,r)/ntags(i,n,x,a);
                tag_prop_not_rec(i,n,x,a)=not_rec(i,n,x,a)/ntags(i,n,x,a);
               }
              if(ntags(i,n,x,a)==0)
               {                   
                tag_prop(i,n,x,a,y,j,r)=0;
                tag_prop_not_rec(i,n,x,a)=0;
               }
             } 
            }
           }
          }
         }
        }
       }

  
 //in order to use the multinomial RNG, need to have a vector of probabilities
 //this essentially requires stacking the tag_prop array into vectors for each release cohort covering all recap states (recap year, population, region)
 //need to extract each recap prob vector and store, then combine into array where can extract the last index (essentially the columns of the array)
 //can then use the fill.multinomial with that vector of probabilities
 if(tag_fit_ages_switch==0) // fit cohorts by region and age
 {
 nreg_temp=rowsum(nregions_temp);
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
       for (int a=1;a<=nages;a++) //release age 
        {
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
              tag_prop_temp(j,y,r)=0; //for whatever reason ADMB won't let a 3darray=double...this is workaround for total_recap_temp=0;, ie setting temp 3darray to 0 after loops through all years 
             }
            }
           }
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++) //recap region
              {
               tag_prop_temp(j,y,r)=tag_prop(i,n,x,a,y,j,r); //fill temp array with a single release cohort's recapture probabilities (excluding not recaptured)
              }
            }
           }
         //tag_prop_temp2=0;
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
              {
               tag_prop_temp2(i,n,x,a,((y-1)*sum(nregions)+nreg_temp(j)+r))=tag_prop_temp(j,y,r); //stack array into single vector by year, population, region
              }
             }
            }
             for(int s=1;s<=(max_life_tags*sum(nregions)+1);s++) //create temp array that has columns of recap prob for each release cohort and add not recap probability to final entry of temp array
              {
              if(s<(min(max_life_tags,nyrs-xx+1)*sum(nregions)+1))
              {
               tag_prop_final(i,n,x,a,s)=tag_prop_temp2(i,n,x,a,s);
              }
              if(s==(min(max_life_tags,nyrs-xx+1)*sum(nregions)+1) && max_life_tags<=(nyrs-xx+1)) //add not recap probability to final entry of temp array
              {
               tag_prop_final(i,n,x,a,s)=tag_prop_not_rec(i,n,x,a);  //for estimation model will use this version of tag_prop in likelihood
              }
              if(s==(min(max_life_tags,nyrs-xx+1)*sum(nregions)+1) && max_life_tags>(nyrs-xx+1)) //add not recap probability to final entry of temp array with adjustment for release events where model ends before max_life_tags (so NR remains in last state)
              {
               tag_prop_final(i,n,x,a,(max_life_tags*sum(nregions)+1))=tag_prop_not_rec(i,n,x,a);  //for estimation model will use this version of tag_prop in likelihood
              }
            }
           }
          }
         }
        }
       }
 if(tag_fit_ages_switch==1) //only fit cohorts by region not age
 {
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age //because accounting for release age, don't need to account for recap age, just adjust mortality and T, to use plus group value if recapture age exceeds max age
        {
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
               tag_recap_no_age_temp(i,n,x,y,j,r,a)=recaps(i,n,x,a,y,j,r);
               total_rec_no_age(i,n,x)=sum(total_rec(i,n,x));
               not_rec_no_age(i,n,x)=sum(not_rec(i,n,x));
             }
            }
           }
          }
        ntags_no_age(i,n,x)=sum(ntags(i,n,x));
       }
      }
     }

 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
      for (int a=1;a<=nages;a++) //release age //because accounting for release age, don't need to account for recap age, just adjust mortality and T, to use plus group value if recapture age exceeds max age
        {
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
              if(ntags_no_age(i,n,x)>0)
               {
                tag_prop_no_age(i,n,x,y,j,r)=sum(tag_recap_no_age_temp(i,n,x,y,j,r))/ntags_no_age(i,n,x);
                tag_prop_not_rec_no_age(i,n,x)=not_rec_no_age(i,n,x)/ntags_no_age(i,n,x);
               }
              if(ntags_no_age(i,n,x)==0)
               {                   
                tag_prop_no_age(i,n,x,y,j,r)=0;
                tag_prop_not_rec_no_age(i,n,x)=0;
               }
             } 
            }
           }
          }
         }
        }
       }
 nreg_temp=rowsum(nregions_temp);
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
             {
              tag_prop_temp(j,y,r)=0; //for whatever reason ADMB won't let a 3darray=double...this is workaround for total_recap_temp=0;, ie setting temp 3darray to 0 after loops through all years 
             }
            }
           }
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++) //recap region
              {
               tag_prop_temp(j,y,r)=tag_prop_no_age(i,n,x,y,j,r); //fill temp array with a single release cohort's recapture probabilities (excluding not recaptured)
              }
            }
           }
         for(int y=1;y<=min(max_life_tags,nyrs-xx+1);y++)  //recap year
          {
           for(int j=1;j<=npops;j++) //recap stock
            {
             for (int r=1;r<=nregions(j);r++)
              {
               tag_prop_temp2_no_age(i,n,x,((y-1)*sum(nregions)+nreg_temp(j)+r))=tag_prop_temp(j,y,r); //stack array into single vector by year, population, region
              }
             }
            }
             for(int s=1;s<=(max_life_tags*sum(nregions)+1);s++) //create temp array that has columns of recap prob for each release cohort and add not recap probability to final entry of temp array
              {
              if(s<(min(max_life_tags,nyrs-xx+1)*sum(nregions)+1))
              {
               tag_prop_final_no_age(i,n,x,s)=tag_prop_temp2_no_age(i,n,x,s);
              }
              if(s==(min(max_life_tags,nyrs-xx+1)*sum(nregions)+1) && max_life_tags<=(nyrs-xx+1)) //add not recap probability to final entry of temp array
              {
               tag_prop_final_no_age(i,n,x,s)=tag_prop_not_rec_no_age(i,n,x);  //for estimation model will use this version of tag_prop in likelihood
              }
              if(s==(min(max_life_tags,nyrs-xx+1)*sum(nregions)+1) && max_life_tags>(nyrs-xx+1)) //add not recap probability to final entry of temp array with adjustment for release events where model ends before max_life_tags (so NR remains in last state)
              {
               tag_prop_final_no_age(i,n,x,(max_life_tags*sum(nregions)+1))=tag_prop_not_rec_no_age(i,n,x);  //for estimation model will use this version of tag_prop in likelihood
              }
            }
           }
          }
         }
        } //end tag_age_switch_loop
       } //end_do_tag loop
FUNCTION get_observed_tag_recaptures
 if(do_tag==1)
  {
  
 random_number_generator myrand_tag(myseed_tag);
 
 if(tag_fit_ages_switch==0) // fit cohorts by region and age
 {
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
       for (int a=1;a<=nages;a++) //release age 
        {
         rand_tag_prop_temp(i,n,x,a).fill_multinomial(myrand_tag,value(tag_prop_final(i,n,x,a))); //fill multinomial using recap prop for each cohort
        }
       }
      }
     }
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
       for (int a=1;a<=nages;a++) //release age 
        {
          rand_tag_prop_temp2=0;
            for(int s=1;s<=SIM_ntag;s++) /// look into changing this so can have ntag change by year (ie different sample sizes for beginning and end of timeseries)
             {
               rand_tag_prop_temp2(value(rand_tag_prop_temp(i,n,x,a,s)))+= 1.0; //count up number in each recap state
              }
           SIM_tag_prop(i,n,x,a)=rand_tag_prop_temp2; //put number in each recap state into appropriate release cohort
        }
       }
      }
     }

 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
      xx=yrs_releases(x);
       for (int a=1;a<=nages;a++) //release age 
        {
         for(int s=1;s<=(max_life_tags*sum(nregions)+1);s++) //create temp array that has columns of recap prob for each release cohort and add not recap probability to final entry of temp array
          {           
           OBS_tag_prop_final(i,n,x,a,s)=SIM_tag_prop(i,n,x,a,s)/SIM_ntag;
              if(ntags(i,n,x,a)==0)
               {                   
                OBS_tag_prop_final(i,n,x,a,s)=0; 
               }
           //summarize by population for panmictic EM
             OBS_tag_prop_population_temp(i,x,a,s,n)=SIM_tag_prop(i,n,x,a,s);
             OBS_tag_prop_population_temp2(i,x,a,s)=sum(OBS_tag_prop_population_temp(i,x,a,s));
             OBS_tag_prop_population_final(i,x,a,s)=OBS_tag_prop_population_temp2(i,x,a,s)/(nregions(i)*SIM_ntag);
           //  OBS_tag_prop_pan_temp(x,a,s,i)=OBS_tag_prop_population_temp2(i,x,a,s);
         //    OBS_tag_prop_pan_temp2(x,a,s)=sum(OBS_tag_prop_pan_temp(x,a,s));

     //        if(OM_structure>0 && EM_structure==0){
      //            for (int k=1;k<=max_life_tags+1;k++)
       //               {
       //                OBS_tag_prop_pan_final_temp(x,a,k)=OBS_tag_prop_pan_temp2(x,a,k)+OBS_tag_prop_pan_temp2(x,a,k+max_life_tags);
       //                OBS_tag_prop_pan_final_temp(x,a,max_life_tags+1)=OBS_tag_prop_pan_temp2(x,a,max_life_tags*sum(nregions)+1);// this is a nightmare
        //               OBS_tag_prop_pan_final(x,a,k)=OBS_tag_prop_pan_final_temp(x,a,k)/(SIM_ntag*sum(nregions));//this probably wont work with 3 areas
                       //OBS_tag_prop_pan_final(x,a,k)=1; //for place holding only
       //                }
        //              }
                    }
                    }
                   }
                  }
                 }
                }

 if(tag_fit_ages_switch==1) // fit cohorts by region
 {
  for (int i=1;i<=npops;i++)
   {
    for (int n=1;n<=nregions(i);n++)
     {
      for(int x=1; x<=nyrs_release; x++)
       {
        rand_tag_prop_temp_no_age(i,n,x).fill_multinomial(myrand_tag,value(tag_prop_final_no_age(i,n,x))); //fill multinomial using recap prop for each cohort
       }
     }
     }
     
 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
    for(int x=1; x<=nyrs_release; x++)
     {
          rand_tag_prop_temp2=0;
            for(int s=1;s<=SIM_ntag;s++) /// look into changing this so can have ntag change by year (ie different sample sizes for beginning and end of timeseries)
             {
               rand_tag_prop_temp2(value(rand_tag_prop_temp_no_age(i,n,x,s)))+= 1.0; //count up number in each recap state
              }
           SIM_tag_prop_no_age(i,n,x)=rand_tag_prop_temp2; //put number in each recap state into appropriate release cohort
        }
       }
      }
     

 for (int i=1;i<=npops;i++)
  {
   for (int n=1;n<=nregions(i);n++)
    {
   for(int x=1; x<=nyrs_release; x++)
    {
     xx=yrs_releases(x);
         for(int s=1;s<=(max_life_tags*sum(nregions)+1);s++) //create temp array that has columns of recap prob for each release cohort and add not recap probability to final entry of temp array
          {           
           OBS_tag_prop_final_no_age(i,n,x,s)=SIM_tag_prop_no_age(i,n,x,s)/SIM_ntag;
              if(ntags_no_age(i,n,x)==0)
               {                   
                OBS_tag_prop_final_no_age(i,n,x,s)=0; 
              }
           //summarize by population for panmictic EM
//             OBS_tag_prop_population_temp_no_age(i,x,s,n)=SIM_tag_prop_no_age(i,n,x,s);
//             OBS_tag_prop_population_temp2_no_age(i,x,s)=sum(OBS_tag_prop_population_temp_no_age(i,x,s));
 //            OBS_tag_prop_population_final_no_age(i,x,s)=OBS_tag_prop_population_temp2_no_age(i,x,s)/(nregions(i)*SIM_ntag);
           //  OBS_tag_prop_pan_temp(x,a,s,i)=OBS_tag_prop_population_temp2(i,x,a,s);
         //    OBS_tag_prop_pan_temp2(x,a,s)=sum(OBS_tag_prop_pan_temp(x,a,s));

     //        if(OM_structure>0 && EM_structure==0){
      //            for (int k=1;k<=max_life_tags+1;k++)
       //               {
       //                OBS_tag_prop_pan_final_temp(x,a,k)=OBS_tag_prop_pan_temp2(x,a,k)+OBS_tag_prop_pan_temp2(x,a,k+max_life_tags);
       //                OBS_tag_prop_pan_final_temp(x,a,max_life_tags+1)=OBS_tag_prop_pan_temp2(x,a,max_life_tags*sum(nregions)+1);// this is a nightmare
        //               OBS_tag_prop_pan_final(x,a,k)=OBS_tag_prop_pan_final_temp(x,a,k)/(SIM_ntag*sum(nregions));//this probably wont work with 3 areas
                       //OBS_tag_prop_pan_final(x,a,k)=1; //for place holding only
       //                }
                     // }
                    }
                    }
                   }
                  }
                 
                } //end do_age_tag_switch loop
  } //end do_tag loop

 //EOUT(recaps(1))
FUNCTION evaluate_the_objective_function
   f=0.0;

 if(active(log_F_est)) //penalties only apply if F is being estimated
 {
 if(MSY_model_type_switch==1 || MSY_model_type_switch==2) //ensure input harvest/TAC matches realized TAC but only does anything if F_est active
  {
   for (int y=1;y<=ny;y++)
    {
     for (int j=1;j<=npops;j++)
      {
       for (int r=1;r<=nregions(j);r++)
        {
         for (int z=1;z<=nfleets(j);z++)
          {
           res_TAC(j,r,z,y)=input_TAC(j,r,z)-yield_fleet(j,r,y,z);
          }
         res_u(j,r,y)=sum(input_u(j,r))-harvest_rate_region_bio(j,r,y);  //sum input_u because it is fleet specific
        }
       }
      }
     if(MSY_model_type_switch==1)
      {
       f=sum(res_TAC);
      }
     if(MSY_model_type_switch==2)
      {
       f=sum(res_u);
      }
    }
          
 if(MSY_model_type_switch==3) // search for F_MSY (max yield across entire model domain, on average, across last 10 years of sim)
  {
   for (int y=1;y<=nyrs_quasi_equil;y++)
    {
     yield_MSY_temp(y)=yield_total(nyrs-nyrs_quasi_equil+y);
    }
   f+=norm2(F_est); //penalize large Fs
   f+=-1.0*(sum(yield_MSY_temp)/nyrs_quasi_equil);   //maximize yield on average over last X years of sim (quasi-equilibrium...instead of just terminal year yield maximization)
  }
 }
REPORT_SECTION

 //for the mismatch calculate the abundance fraction by region/population
  for (int y=1;y<=nyrs;y++)
   {
   for (int p=1;p<=npops;p++)
    {
    for (int r=1;r<=nregions(p);r++)
      {
      for (int z=1;z<=nfleets(p);z++)
        {
         for(int a=1;a<=nages;a++)
          {
        abund_frac_age_region(p,r,y,a)=abundance_at_age_AM(p,r,y,a)/abundance_total(y,a);
        abund_frac_region_year(p,r,y)=sum(abundance_at_age_AM(p,r,y))/sum(abundance_total(y));
        abund_frac_region(p,r)=sum(abund_frac_region_year(p,r))/nyrs;//average of all years
        }}}}}


 //Aggregating OBS values and vitals for panmictic EM
 //if(EM_structure==0 && OM_structure>0){ 
 // for (int y=1;y<=nyrs;y++)
 //  {
 //  for (int p=1;p<=npops;p++)
 //   {
 //   for (int r=1;r<=nregions(p);r++)
 //     {
 //     for (int z=1;z<=nfleets(p);z++)
 //       {
 //        for(int a=1;a<=nages;a++)
 //         {
       //aggregating weight at age
 //       input_weight_region_temp(p,a,r)=input_weight(p,r,a)*abund_frac_region(p,r);//rearrange to summarize and weight for output
 //       input_weight_region(p,a)=sum(input_weight_region_temp(p,a));
 //       input_weight_population_temp(a,p)=input_weight_region(p,a);
 //       input_weight_population(a)=sum(input_weight_population_temp(a));

        //aggregating catch weight at age
 //       input_catch_weight_region_temp(p,a,r)=input_catch_weight(p,r,a)*abund_frac_region(p,r);//sum by region
 //       input_catch_weight_region(p,a)=sum(input_catch_weight_region_temp(p,a));
 //       input_catch_weight_population_temp(a,p)=input_catch_weight_region(p,a);
 //       input_catch_weight_population(a)=sum(input_catch_weight_population_temp(a));

        //aggregating fecundity
 //       fecundity_region_temp(p,a,r)=fecundity(p,r,a)*abund_frac_region(p,r);//sum by region
 ///       fecundity_region(p,a)=sum(fecundity_region_temp(p,a));
 //       fecundity_population_temp(a,p)=fecundity_region(p,a);
 //       fecundity_population(a)=sum(fecundity_population_temp(a));

        //aggregating maturity
 //       maturity_region_temp(p,a,r)=maturity(p,r,a)*abund_frac_region(p,r);//sum by region
 //       maturity_region(p,a)=sum(maturity_region_temp(p,a));
 //       maturity_population_temp(a,p)=maturity_region(p,a);
  //      maturity_population(a)=sum(maturity_population_temp(a));

        //aggregating selectivity
 //       selectivity_age_temp(p,a,z,r)=selectivity_age(p,r,a,z)*abund_frac_region(p,r);
 //       selectivity_age_pop(p,a,z)=sum(selectivity_age_temp(p,a,z));
 //       survey_selectivity_age_temp(p,a,z,r)=survey_selectivity_age(p,r,a,z)*abund_frac_region(p,r);
 //       survey_selectivity_age_pop(p,a,z)=sum(survey_selectivity_age_temp(p,a,z));

        //aggregating the age comps

        //survey
 //        OBS_survey_prop_temp(p,r,y,a,z)=OBS_survey_prop(p,r,z,y,a);
 //        OBS_survey_prop_temp2(p,r,y,a)=sum(OBS_survey_prop_temp(p,r,y,a));
 //       OBS_survey_prop_temp3(p,r,y,a)=OBS_survey_prop_temp2(p,r,y,a)*abund_frac_age_region(p,r,y,a);
 //       OBS_survey_prop_temp4(p,y,a,r)=OBS_survey_prop_temp3(p,r,y,a);
 //       OBS_survey_prop_population(p,y,a)=sum(OBS_survey_prop_temp4(p,y,a));
 //       OBS_survey_prop_pan_temp(y,a,p)= OBS_survey_prop_population(p,y,a);
 //       OBS_survey_prop_pan(y,a)=sum(OBS_survey_prop_pan_temp(y,a));

        //catch
 //       OBS_catch_prop_temp(p,r,y,a,z)= OBS_catch_prop(p,r,z,y,a);
 //       OBS_catch_prop_temp2(p,r,y,a)=sum(OBS_catch_prop_temp(p,r,y,a));
 //       OBS_catch_prop_temp3(p,r,y,a)= OBS_catch_prop_temp2(p,r,y,a)*abund_frac_age_region(p,r,y,a);
 //       OBS_catch_prop_temp4(p,y,a,r)= OBS_catch_prop_temp3(p,r,y,a);
 //       OBS_catch_prop_population(p,y,a)=sum(OBS_catch_prop_temp4(p,y,a));
 //       OBS_catch_prop_pan_temp(y,a,p)= OBS_catch_prop_population(p,y,a);
 //       OBS_catch_prop_pan(y,a)=sum(OBS_catch_prop_pan_temp(y,a));


 //       } //end age loop
 //       } //end fleets loop

        //proportion female
 //       prop_fem_temp(p,r)= prop_fem(p,r)*abund_frac_region(p,r); 
 //       prop_fem_pan=sum(prop_fem_temp);
 //       rec_index_temp(p,y,r)=rec_index_BM(p,r,y)*abund_frac_region_year(p,r,y); //rearrange and weight for summing

 //      } //end reg loop

       //rec index
 //       rec_index_BM_population(p,y)=sum(rec_index_temp(p,y));// combined by region
 //       rec_index_temp2(y,p)=rec_index_BM_population(p,y);
 //       rec_index_pan(y)=sum(rec_index_temp2(y));//npops; combined by populations
        
 //     } //end pop loop          
 //    } //end year loop
 // }


//Additional model structure parameters
  report<<"#nages"<<endl;
  report<<nages<<endl;
  report<<"#nyrs"<<endl;
  report<<nyrs<<endl;

//EM structure 
  report<<"#npops_EM"<<endl;
  report<<npops_EM<<endl;
  report<<"#nregions_EM"<<endl;
  report<<nregions_EM<<endl;
  report<<"#nfleets_EM"<<endl;
  report<<nfleets_EM<<endl;
  report<<"#nfleets_survey_EM"<<endl;
  report<<nfleets_survey_EM<<endl;
  
//OM structure 
  report<<"#npops_OM"<<endl;
  report<<npops<<endl;
  report<<"#nregions_OM"<<endl;
  report<<nregions<<endl;
  report<<"#nfleets_OM"<<endl;
  report<<nfleets<<endl;
  report<<"#nfleets_survey_OM"<<endl;
  report<<nfleets_survey<<endl;

//EM parameters input from OM .dat
  report<<"#tsurvey_EM"<<endl;
  report<<tsurvey_EM<<endl;
  report<<"#diagnostics_switch"<<endl;
  report<<diagnostics_switch<<endl;
  report<<"#move_switch"<<endl;
  report<<move_switch_EM<<endl;
  report<<"#report_rate_switch"<<endl;
  report<<report_rate_switch_EM<<endl;
  report<<"#natal_homing_switch"<<endl;
  report<<natal_homing_switch_EM<<endl;
  report<<"#spawn_return_switch"<<endl;
  report<<spawn_return_switch_EM<<endl;
  report<<"#select_switch"<<endl;
  report<<select_switch_EM<<endl;
  report<<"#select_switch_survey"<<endl;
  report<<select_switch_survey_EM<<endl;
  report<<"#maturity_switch_equil"<<endl;
  report<<maturity_switch_equil_EM<<endl;
  report<<"#SSB_type"<<endl;
  report<<SSB_type_EM<<endl;
  report<<"#Rec_type"<<endl;
  report<<Rec_type_EM<<endl;
  report<<"#apportionment_type"<<endl;
  report<<apportionment_type_EM<<endl;
  report<<"#use_stock_comp_info_survey"<<endl;
  report<<use_stock_comp_info_survey_EM<<endl;
  report<<"#use_stock_comp_info_catch"<<endl;
  report<<use_stock_comp_info_catch_EM<<endl;
  report<<"#F_switch"<<endl;
  report<<F_switch_EM<<endl;
  report<<"#M_switch"<<endl;
  report<<M_switch_EM<<endl;
  report<<"#recruit_devs_switch"<<endl;
  report<<recruit_devs_switch_EM<<endl;
  report<<"#recruit_randwalk_switch"<<endl;
  report<<recruit_randwalk_switch_EM<<endl;
  report<<"#init_abund_switch"<<endl;
  report<<init_abund_switch_EM<<endl;
  report<<"#est_dist_init_abund"<<endl;
  report<<est_dist_init_abund_EM<<endl;
  report<<"#tspawn_EM"<<endl;
  report<<tspawn_EM<<endl;
  report<<"#return_age"<<endl;
  report<<return_age<<endl;
  report<<"#return_probability"<<endl;
  report<<return_probability_EM<<endl;
  report<<"#spawn_return_prob"<<endl;
  report<<spawn_return_prob_EM<<endl;
  report<<"#do_tag"<<endl;
  report<<do_tag_EM<<endl;
  report<<"#fit_tag_age_switch"<<endl;
  report<<tag_fit_ages_switch<<endl;
  report<<"#do_tag_mult"<<endl;
  report<<do_tag_mult<<endl;
  report<<"#est_tag_mixing_switch"<<endl;
  report<<est_tag_mixing_switch<<endl;

  report<<"#sigma_recruit_EM"<<endl;
  report<<sigma_recruit_EM<<endl;

  report<<"#ph_lmr"<<endl;
  report<<ph_lmr<<endl;
  report<<"#Rave_start"<<endl;
  report<<Rave_start<<endl;
  report<<"#lb_R_ave"<<endl;
  report<<lb_R_ave<<endl;
  report<<"#ub_R_ave"<<endl;
  report<<ub_R_ave<<endl;
  report<<"#ph_rec"<<endl;
  report<<ph_rec<<endl;
  report<<"#lb_rec_devs"<<endl;
  report<<lb_rec_devs<<endl;
  report<<"#ub_rec_devs"<<endl;
  report<<ub_rec_devs<<endl;
  report<<"#Rdevs_start"<<endl;
  report<<Rdevs_start<<endl;

  report<<"#ph_rec_app_CNST"<<endl;
  report<<ph_rec_app_CNST<<endl;
  report<<"#ph_rec_app_YR"<<endl;
  report<<ph_rec_app_YR<<endl;
  report<<"#lb_rec_app"<<endl;
  report<<lb_rec_app<<endl;
  report<<"#ub_rec_app"<<endl;
  report<<ub_rec_app<<endl;
  report<<"#Rapp_start"<<endl;
  report<<Rapp_start<<endl;

  report<<"#ph_init_abund"<<endl;
  report<<ph_init_abund<<endl;
  report<<"#ph_init_abund_no_ag1"<<endl;
  report<<ph_init_abund_no_ag1<<endl;
  report<<"#N_start"<<endl;
  report<<N_start<<endl;
  report<<"#init_dist_start"<<endl;
  report<<init_dist_start<<endl;
  report<<"#ph_reg_init"<<endl;
  report<<ph_reg_init<<endl;
  report<<"#ph_reg_age_init"<<endl;
  report<<ph_reg_age_init<<endl;
  report<<"#ph_non_natal_init"<<endl;
  report<<ph_non_natal_init<<endl;
  report<<"#ph_non_natal_age_init"<<endl;
  report<<ph_non_natal_age_init<<endl;
  report<<"#lb_init_dist"<<endl;
  report<<lb_init_dist<<endl;
  report<<"#ub_init_dist"<<endl;
  report<<ub_init_dist<<endl;
  report<<"#lb_init_abund"<<endl;
  report<<lb_init_abund<<endl;
  report<<"#ub_init_abund"<<endl;
  report<<ub_init_abund<<endl;
  report<<"#ph_F"<<endl;
  report<<ph_F<<endl;
  report<<"#lb_F"<<endl;
  report<<lb_F<<endl;
  report<<"#ub_F"<<endl;
  report<<ub_F<<endl;
  report<<"#F_start"<<endl;
  report<<F_start<<endl;
  report<<"#ph_steep"<<endl;
  report<<ph_steep<<endl;
  report<<"#lb_steep"<<endl;
  report<<lb_steep<<endl;
  report<<"#ub_steep"<<endl;
  report<<ub_steep<<endl;
  report<<"#steep_start"<<endl;
  report<<steep_start<<endl;
  report<<"#ph_M_CNST"<<endl;
  report<<ph_M_CNST<<endl;
  report<<"#ph_M_pop_CNST"<<endl;
  report<<ph_M_pop_CNST<<endl;
  report<<"#ph_M_age_CNST"<<endl;
  report<<ph_M_age_CNST<<endl;
  report<<"#ph_M_pop_age"<<endl;
  report<<ph_M_pop_age<<endl;
  report<<"#lb_M"<<endl;
  report<<lb_M<<endl;
  report<<"#ub_M"<<endl;
  report<<ub_M<<endl;
  report<<"#M_start"<<endl;
  report<<M_start<<endl;
  report<<"#ph_sel_log"<<endl;
  report<<ph_sel_log<<endl;
  report<<"#lb_sel_beta1"<<endl;
  report<<lb_sel_beta1<<endl;
  report<<"#ub_sel_beta1"<<endl;
  report<<ub_sel_beta1<<endl;
  report<<"#sel_beta1_start"<<endl;
  report<<sel_beta1_start<<endl;
  report<<"#lb_sel_beta2"<<endl;
  report<<lb_sel_beta2<<endl;
  report<<"#ub_sel_beta2"<<endl;
  report<<ub_sel_beta2<<endl;
  report<<"#sel_beta2_start"<<endl;
  report<<sel_beta2_start<<endl;
  report<<"#lb_sel_beta3"<<endl;
  report<<lb_sel_beta3<<endl;
  report<<"#ub_sel_beta3"<<endl;
  report<<ub_sel_beta3<<endl;
  report<<"#sel_beta3_start"<<endl;
  report<<sel_beta3_start<<endl;
  report<<"#lb_sel_beta4"<<endl;
  report<<lb_sel_beta4<<endl;
  report<<"#ub_sel_beta4"<<endl;
  report<<ub_sel_beta4<<endl;
  report<<"#sel_beta4_start"<<endl;
  report<<sel_beta4_start<<endl;
  report<<"#lb_sel_beta1_surv"<<endl;
  report<<lb_sel_beta1_surv<<endl;
  report<<"#ub_sel_beta1_surv"<<endl;
  report<<ub_sel_beta1_surv<<endl;
  report<<"#sel_beta1_surv_start"<<endl;
  report<<sel_beta1_surv_start<<endl;
  report<<"#lb_sel_beta2_surv"<<endl;
  report<<lb_sel_beta2_surv<<endl;
  report<<"#ub_sel_beta2_surv"<<endl;
  report<<ub_sel_beta2_surv<<endl;
  report<<"#sel_beta2_surv_start"<<endl;
  report<<sel_beta2_surv_start<<endl;
  report<<"#lb_sel_beta3_surv"<<endl;
  report<<lb_sel_beta3_surv<<endl;
  report<<"#ub_sel_beta3_surv"<<endl;
  report<<ub_sel_beta3_surv<<endl;
  report<<"#sel_beta3_surv_start"<<endl;
  report<<sel_beta3_surv_start<<endl;
  report<<"#lb_sel_beta4_surv"<<endl;
  report<<lb_sel_beta4_surv<<endl;
  report<<"#ub_sel_beta4_surv"<<endl;
  report<<ub_sel_beta4_surv<<endl;
  report<<"#sel_beta4_surv_start"<<endl;
  report<<sel_beta4_surv_start<<endl;
  report<<"#ph_sel_log_surv"<<endl;
  report<<ph_sel_log_surv<<endl;
  report<<"#ph_sel_dubl"<<endl;
  report<<ph_sel_dubl<<endl;
  report<<"#ph_sel_dubl_surv"<<endl;
  report<<ph_sel_dubl_surv<<endl;
  report<<"#ph_q"<<endl;
  report<<ph_q<<endl;
  report<<"#lb_q"<<endl;
  report<<lb_q<<endl;
  report<<"#ub_q"<<endl;
  report<<ub_q<<endl;
  report<<"#q_start"<<endl;
  report<<q_start<<endl;
  report<<"#ph_F_rho"<<endl;
  report<<ph_F_rho<<endl;
  report<<"#lb_F_rho"<<endl;
  report<<lb_F_rho<<endl;
  report<<"#ub_F_rho"<<endl;
  report<<ub_F_rho<<endl;
  report<<"#Frho_start"<<endl;
  report<<Frho_start<<endl;
  report<<"#phase_T_YR"<<endl;
  report<<phase_T_YR<<endl;
  report<<"#phase_T_YR_ALT_FREQ"<<endl;
  report<<phase_T_YR_ALT_FREQ<<endl;
  report<<"#T_est_freq"<<endl;
  report<<T_est_freq<<endl;
  report<<"#phase_T_YR_AGE_ALT_FREQ"<<endl;
  report<<phase_T_YR_AGE_ALT_FREQ<<endl;
  report<<"#T_est_age_freq"<<endl;
  report<<T_est_age_freq<<endl;
  report<<"#juv_age"<<endl;
  report<<juv_age<<endl;  
  report<<"#phase_T_CNST"<<endl;
  report<<phase_T_CNST<<endl;
  report<<"#phase_T_CNST_AGE"<<endl;
  report<<phase_T_CNST_AGE<<endl;
  report<<"#phase_T_YR_AGE"<<endl;
  report<<phase_T_YR_AGE<<endl;
  report<<"#phase_T_CNST_AGE_no_AG1"<<endl;
  report<<phase_T_CNST_AGE_no_AG1<<endl;
  report<<"#phase_T_YR_AGE_no_AG1"<<endl;
  report<<phase_T_YR_AGE_no_AG1<<endl;
  report<<"#phase_T_YR_AGE_ALT_FREQ_no_AG1"<<endl;
  report<<phase_T_YR_AGE_ALT_FREQ_no_AG1<<endl;
  report<<"#lb_T"<<endl;
  report<<lb_T<<endl;
  report<<"#ub_T"<<endl;
  report<<ub_T<<endl;
  report<<"#T_start"<<endl;
  report<<T_start<<endl;
  report<<"#phase_rep_rate_YR"<<endl;
  report<<phase_rep_rate_YR<<endl;
  report<<"#phase_rep_rate_CNST"<<endl;
  report<<phase_rep_rate_CNST<<endl;
  report<<"#lb_B"<<endl;
  report<<lb_B<<endl;
  report<<"#ub_B"<<endl;
  report<<ub_B<<endl;
  report<<"#B_start"<<endl;
  report<<B_start<<endl;
  
  report<<"#ph_T_tag"<<endl;
  report<<ph_T_tag<<endl;
  report<<"#ph_F_tag"<<endl;
  report<<ph_F_tag<<endl;
  report<<"#lb_scalar_T"<<endl;
  report<<lb_scalar_T<<endl;
  report<<"#ub_scalar_T"<<endl;
  report<<ub_scalar_T<<endl;
  report<<"#lb_scalar_F"<<endl;
  report<<lb_scalar_F<<endl;
  report<<"#ub_scalar_F"<<endl;
  report<<ub_scalar_F<<endl;
  report<<"#scalar_T_start"<<endl;
  report<<scalar_T_start<<endl;
  report<<"#scalar_F_start"<<endl;
  report<<scalar_F_start<<endl;
  
  report<<"#ph_dummy"<<endl;
  report<<ph_dummy_EM<<endl;
  report<<"#wt_surv"<<endl;
  report<<wt_surv<<endl;
  report<<"#wt_catch"<<endl;
  report<<wt_catch<<endl;
  report<<"#wt_fish_age"<<endl;
  report<<wt_fish_age<<endl;
  report<<"#wt_srv_age"<<endl;
  report<<wt_srv_age<<endl;
  report<<"#wt_rec"<<endl;
  report<<wt_rec<<endl;
  report<<"#wt_tag"<<endl;
  report<<wt_tag<<endl;
  report<<"#wt_F_pen"<<endl;
  report<<wt_F_pen<<endl;
  report<<"#wt_M_pen"<<endl;
  report<<wt_M_pen<<endl;
  report<<"#wt_B_pen"<<endl;
  report<<wt_B_pen<<endl;
  report<<"#report_rate_sigma"<<endl;
  report<<report_rate_sigma<<endl;
  report<<"#report_rate_ave"<<endl;
  report<<report_rate_ave<<endl;

  report<<"#abund_pen_switch"<<endl;
  report<<abund_pen_switch<<endl;
  report<<"#wt_abund_pen"<<endl;
  report<<wt_abund_pen<<endl;
  report<<"#Mean_N"<<endl;
  report<<Mean_N<<endl;
  report<<"#abund_dist_pen_switch"<<endl;
  report<<abund_dist_pen_switch<<endl;
  report<<"#wt_abund_dist_pen"<<endl;
  report<<wt_abund_dist_pen<<endl;
  report<<"#Mean_N_dist"<<endl;
  report<<Mean_N_dist<<endl;
  report<<"#move_pen_switch"<<endl;
  report<<move_pen_switch<<endl;
  report<<"#wt_move_pen"<<endl;
  report<<wt_move_pen<<endl;
  report<<"#Tpen"<<endl;
  report<<Tpen<<endl;
  report<<"#sigma_Tpen_EM"<<endl;
  report<<sigma_Tpen_EM<<endl;
  report<<"#Rave_pen_switch"<<endl;
  report<<Rave_pen_switch<<endl;
  report<<"#wt_Rave_pen"<<endl;
  report<<wt_Rave_pen<<endl;
  report<<"#Rave_mean"<<endl;
  report<<Rave_mean<<endl;
  report<<"#R_app_pen_switch"<<endl;
  report<<R_app_pen_switch<<endl;
  report<<"#wt_R_app_pen"<<endl;
  report<<wt_R_app_pen<<endl;
  report<<"#R_app_pen"<<endl;
  report<<R_app_pen<<endl;
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
///REPORTING THE CORRECT EM PARAMETERS///////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

//Spatial to panmictic EM inputs
//   if(EM_structure==0 && OM_structure>=0){ 
//      report<<"#input_weight"<<endl;
//      report<<input_weight_population<<endl;
//      report<<"#input_catch_weight"<<endl;
//      report<<input_catch_weight_population<<endl;
//      report<<"#fecundity"<<endl;
//      report<<fecundity_population<<endl;
//      report<<"#maturity"<<endl;
//      report<<maturity_population<<endl;
//      report<<"#prop_fem"<<endl; 
//      report<<prop_fem_pan<<endl;
//      report<<"#OBS_rec_index_BM"<<endl;
//      report<<rec_index_pan<<endl;

 //for straight panmictic EM 
 //   if(sum(nfleets_EM)==1){
 //     report<<"#OBS_survey_fleet"<<endl;
 //     report<<OBS_survey_total_bio<<endl;
 //     report<<"#OBS_survey_fleet_bio_se_EM"<<endl;
 //     report<<OBS_survey_fleet_bio_se_EM<<endl;
 //     report<<"#OBS_survey_prop"<<endl;
 //     report<<OBS_survey_prop_pan<<endl;
 //     report<<"#OBS_survey_prop_N_EM"<<endl;
 //     report<<OBS_survey_prop_N_EM<<endl;
 //     report<<"#OBS_yield_fleet"<<endl;
 //     report<<OBS_yield_total<<endl;
 //     report<<"#OBS_yield_fleet_se_EM"<<endl;
 //     report<<OBS_yield_fleet_se_EM<<endl;    
 //     report<<"#OBS_catch_prop"<<endl;
 //     report<<OBS_catch_prop_pan<<endl; 
 //     report<<"#OBS_catch_prop_N_EM"<<endl;
 //     report<<OBS_catch_prop_N_EM<<endl;
      
//tagging information
 //     report<<"#nyrs_release"<<endl;
 //     report<<nyrs_release<<endl;
 //     report<<"#years_of_tag_releases "<<endl;
 //     report<<yrs_releases<<endl;
 //     report<<"#max_life_tags"<<endl;
 //     report<<max_life_tags<<endl;
 //     report<<"#age_full_selection"<<endl;
 //     report<<age_full_selection<<endl;
 //     report<<"#input_report_rate_EM"<<endl;
 //     report<<input_report_rate_EM<<endl;
 //     report<<"#ntags"<<endl;
 //     report<<ntags_pan<<endl;
 //     report<<"#ntags_total"<<endl;
 //     report<<ntags_total<<endl;
 //     report<<"#tag_N_EM"<<endl;
 //     report<<tag_N_EM<<endl;
 //     report<<"#input_T_EM"<<endl;
 //     report<<input_T_EM<<endl;
 //    report<<"#OBS_tag_prop_pan_final"<<endl;
 //     report<<OBS_tag_prop_pan_final<<endl;
     // report<<"#OBS_tag_prop_pan_final_no_age"<<endl;
     // report<<OBS_tag_prop_pan_final_no_age<<endl;
   //   }


 //for fleets-as-areas approach
 //     if(sum(nfleets_EM)>1){
      //fleet specific outputs for fishery, panmictic for survey
 //     report<<"#OBS_survey_fleet"<<endl;
 //     report<<OBS_survey_total_bio<<endl;
 //     report<<"#OBS_survey_fleet_bio_se_EM"<<endl;
 //     report<<OBS_survey_fleet_bio_se_EM<<endl;
 //     report<<"#OBS_survey_prop"<<endl;
 //     report<<OBS_survey_prop_pan<<endl;
 //     report<<"#OBS_survey_prop_N_EM"<<endl;
 //    report<<OBS_survey_prop_N_EM<<endl;
 //     report<<"#OBS_yield_fleet"<<endl;
 //     report<<OBS_yield_fleet<<endl;
 //     report<<"#OBS_yield_fleet_se_EM"<<endl;
 //     report<<OBS_yield_fleet_se_EM<<endl;
 //     report<<"#OBS_catch_prop"<<endl;
 //     report<<OBS_catch_prop<<endl;
 //     report<<"#OBS_catch_prop_N_EM"<<endl;
 //     report<<OBS_catch_prop_N_EM<<endl;
      
//tagging information
 //     report<<"#nyrs_release"<<endl;
 //     report<<nyrs_release<<endl;
 //     report<<"#years_of_tag_releases "<<endl;
 //     report<<yrs_releases<<endl;
 //     report<<"#max_life_tags"<<endl;
 //     report<<max_life_tags<<endl;
 //     report<<"#age_full_selection"<<endl;
 //     report<<age_full_selection<<endl;
 //     report<<"#input_report_rate_EM"<<endl;
 //     report<<input_report_rate_EM<<endl;
 //     report<<"#ntags"<<endl;
 //     report<<ntags_pan<<endl;
 //     report<<"#ntags_total"<<endl;
 //     report<<ntags_total<<endl;
 //     report<<"#tag_N_EM"<<endl;
 //     report<<tag_N_EM<<endl;
 //     report<<"#input_T_EM"<<endl;
 //     report<<input_T_EM<<endl;
 //     report<<"#OBS_tag_prop_pan_final"<<endl;
 //     report<<OBS_tag_prop_pan_final<<endl;
     // report<<"#OBS_tag_prop_pan_final_no_age"<<endl;
     // report<<OBS_tag_prop_pan_final_no_age<<endl;
 //     }
 //     }

//spatial to spatial EM inputs or panmictic matching - no aggregation needed
     if((EM_structure>0 && OM_structure>0) || (EM_structure==0 && OM_structure==0)){
      report<<"#input_weight"<<endl;
      report<<input_weight<<endl;
      report<<"#input_catch_weight"<<endl;
      report<<input_catch_weight<<endl;
      report<<"#fecundity"<<endl;
      report<<fecundity<<endl;
      report<<"#maturity"<<endl;
      report<<maturity<<endl;
      report<<"#prop_fem"<<endl; 
      report<<prop_fem<<endl;
      report<<"#OBS_rec_index_BM"<<endl;
      report<<rec_index_BM<<endl;
      report<<"#OBS_survey_fleet_bio"<<endl;
      report<<OBS_survey_fleet_bio<<endl;
      report<<"#OBS_survey_fleet_bio_se_EM"<<endl;
      report<<OBS_survey_fleet_bio_se_EM<<endl;
      report<<"#OBS_survey_prop"<<endl;
      report<<OBS_survey_prop<<endl;
      report<<"#OBS_survey_prop_N_EM"<<endl;
      report<<OBS_survey_prop_N_EM<<endl;
      report<<"#OBS_yield_fleet"<<endl;
      report<<OBS_yield_fleet<<endl;
      report<<"#OBS_yield_fleet_se_EM"<<endl;
      report<<OBS_yield_fleet_se_EM<<endl;
      report<<"#OBS_catch_prop"<<endl;
      report<<OBS_catch_prop<<endl;
      report<<"#OBS_catch_prop_N_EM"<<endl;
      report<<OBS_catch_prop_N_EM<<endl;

//tagging information
      report<<"#nyrs_release"<<endl;
      report<<nyrs_release<<endl;
      report<<"#years_of_tag_releases "<<endl;
      report<<yrs_releases<<endl;
      report<<"#max_life_tags"<<endl;
      report<<max_life_tags<<endl;
      report<<"#age_full_selection"<<endl;
      report<<age_full_selection<<endl;
      report<<"#input_report_rate_EM"<<endl;
      report<<input_report_rate_EM<<endl;
      report<<"#ntags"<<endl;
      report<<ntags<<endl;
      report<<"#ntags_total"<<endl;
      report<<ntags_total<<endl;
      report<<"#tag_N_EM"<<endl;
      report<<tag_N_EM<<endl;
      report<<"#input_T_EM"<<endl;
      report<<input_T_EM<<endl;
      report<<"#OBS_tag_prop_final"<<endl;
      report<<OBS_tag_prop_final<<endl;
      report<<"#OBS_tag_prop_final_no_age"<<endl;
      report<<OBS_tag_prop_final_no_age<<endl;
      }


///////////
// panmictic to spatial will go here eventually
/////////

//Additional inputs for EM specified in OM .dat
  report<<"#input_M_EM"<<endl;
  report<<input_M_EM<<endl;
  report<<"#input_Rec_Prop_EM"<<endl;
  report<<input_rec_prop_EM<<endl;
  report<<"#input_selectivity_EM"<<endl;
  report<<input_selectivity_EM<<endl;
  report<<"#input_survey_selectivity_EM"<<endl;
  report<<input_survey_selectivity_EM<<endl;
  report<<"#input_dist_init_abund"<<endl;
  report<<input_dist_init_abund_EM<<endl;
  report<<"#init_abund_EM"<<endl;
  report<<init_abund_EM<<endl;
/// TRUE VALUES FROM OM
  report<<"#frac_natal_true"<<endl;
  report<<frac_natal_true<<endl;
  report<<"#frac_natal_true_age"<<endl;
  report<<frac_natal_true_age<<endl;  
  report<<"#input_M_TRUE"<<endl;
  report<<input_M_TRUE<<endl;
  report<<"#init_abund_TRUE"<<endl;
  report<<init_abund<<endl;
  report<<"#q_survey"<<endl;
  report<<q_survey<<endl;
  report<<"#sel_beta1"<<endl;
  report<<sel_beta1<<endl;
  report<<"#sel_beta2"<<endl;
  report<<sel_beta2<<endl;
  report<<"#sel_beta3"<<endl;
  report<<sel_beta3<<endl;
  report<<"#sel_beta4"<<endl;
  report<<sel_beta4<<endl;
  report<<"#sel_beta1_survey"<<endl;
  report<<sel_beta1_survey<<endl;
  report<<"#sel_beta2_survey"<<endl;
  report<<sel_beta2_survey<<endl;
  report<<"#sel_beta3_survey"<<endl;
  report<<sel_beta3_survey<<endl;
  report<<"#sel_beta4_survey"<<endl;
  report<<sel_beta4_survey<<endl;
  report<<"#steep"<<endl;
  report<<steep<<endl;
  report<<"#R_ave"<<endl;
  report<<R_ave<<endl;
  report<<"#SSB_zero"<<endl;
  report<<SSB_zero<<endl;
  report<<"#rec_devs"<<endl;
  report<<rec_devs<<endl;
  report<<"#Rec_Prop"<<endl;
  report<<Rec_Prop<<endl;
  report<<"#recruits_BM"<<endl;
  report<<recruits_BM<<endl;
  report<<"#F"<<endl;
  report<<F<<endl;
  report<<"#Fyear"<<endl;
  report<<F_year<<endl;
  report<<"#biomass_AM"<<endl;
  report<<biomass_AM<<endl;
  report<<"#biomass_population"<<endl;
  report<<biomass_population<<endl;

//add to EM input
  report<<"#catch_at_age_fleet_prop"<<endl;
  report<<catch_at_age_fleet_prop<<endl;
  report<<"#yield_fleet"<<endl;
  report<<yield_fleet<<endl;
  report<<"#survey_at_age_fleet_prop"<<endl;
  report<<survey_at_age_fleet_prop<<endl;
  report<<"#true_survey_fleet_bio"<<endl;
  report<<true_survey_fleet_bio<<endl;
//

  report<<"#harvest_rate_region_bio"<<endl;
  report<<harvest_rate_region_bio<<endl;
  report<<"#depletion_region"<<endl;
  report<<depletion_region<<endl;
  report<<"#SSB_region"<<endl;
  report<<SSB_region<<endl;
  report<<"#Bratio_population"<<endl;
  report<<Bratio_population<<endl;
  report<<"#T_year"<<endl;
  report<<T_year<<endl;

//reporting true aggregated selectivity if fleets ==1
  if(EM_structure==0 && sum(nfleets_EM)==1)
  {
  report<<"#selectivity_age"<<endl;
  report<<selectivity_age_pop<<endl;
  report<<"#survey_selectivity_age"<<endl;
  report<<survey_selectivity_age_pop<<endl;
  }
  
  else{
  report<<"#selectivity_age"<<endl;
  report<<selectivity_age<<endl;
  report<<"#survey_selectivity_age"<<endl;
  report<<survey_selectivity_age<<endl;
  }
  
 //true tag information
  report<<"#TRUE_tag_prop"<<endl;
  report<<tag_prop_final<<endl;
  report<<"#TRUE_tag_prop_no_age"<<endl;
  report<<tag_prop_final_no_age<<endl;
  
  report<<"#T"<<endl;
  report<<T<<endl;
  report<<"#report_rate_TRUE"<<endl;
  report<<report_rate_TRUE<<endl;
  
  report<<"#move_switch_OM"<<endl;
  report<<move_switch<<endl;
  report<<"#DD_move_age_switch_OM"<<endl;
  report<<DD_move_age_switch<<endl;
  
//report the abundance_frac for later calcs if needed
  report<<"#abund_frac_age_region"<<endl;
  report<<abund_frac_age_region<<endl;
  report<<"#abund_frac_year"<<endl;
  report<<abund_frac_region_year<<endl;
  report<<"#abund_frac_region"<<endl;
  report<<abund_frac_region<<endl;

  report<<"#F_tag_scalar"<<endl;
  report<<F_tag_scalar<<endl;
  report<<"#T_tag_res"<<endl;
  report<<T_tag_res<<endl;
  report<<"#sim_tag_mixing_switch"<<endl;
  report<<sim_tag_mixing_switch<<endl;
  
  report<<"#debug"<<endl;
  report<<debug<<endl;


//to get abundance after movement for init abund
// report<<"#Abund_BM"<<endl;
// report<<abundance_at_age_BM<<endl;
// report<<"#Abund_AM"<<endl;
// report<<abundance_at_age_AM<<endl;

 // Some stuff for looking at tag calculations 
  /*
  report<<"#OBS_tag_prop_final"<<endl;
  report<<OBS_tag_prop_final(1,1,1)<<endl;
  report<<"nt"<<endl;
  report<<"true_survey_population_bio"<<true_survey_population_bio<<endl;
  report<<"true_survey_region_bio"<<true_survey_region_bio<<endl;
  report<<"abundance_total"<<abundance_total<<endl<<"abundance_population"<<abundance_population<<endl;
  report<<"true fleet_age"<<true_survey_fleet_age<<endl;
  report<<"survey_abundance"<<endl<<true_survey_region_abundance<<endl;
  report<<"temp_fleet_age"<<endl<<true_survey_fleet_age_temp<<endl;


  report<<"ntags(1,1,1,1)"<<endl<<ntags(1,1,1,9)<<endl<<"ntags_total(x)"<<endl<<ntags_total(1)<<endl;
  report<<"true_survey_population_abundance(xx,i,a)"<<endl<<true_survey_population_abundance(1,1,9)<<endl;
  report<<"true_survey_total_abundance(xx,a)"<<endl<<sum(true_survey_total_abundance(1))<<endl;
  report<<"true_survey_region_abundance(i,xx,n,a)"<<endl<<true_survey_region_abundance(1,1,1,9)<<endl;
  report<<"true_survey_population_abundance(xx,i,a)"<<endl<<true_survey_population_abundance(1,1,9)<<endl;
  report<<"survey_selectivity(i,n,xx,a,1)"<<endl<<survey_selectivity(1,1,1,9,1)<<endl;
  report<<"sum(survey_selectivity_temp(i,n,xx,1)))"<<endl<<sum(survey_selectivity_temp(1,1,1,1))<<endl;
 */

RUNTIME_SECTION
  convergence_criteria .001,.0001, 1.0e-4, 1.0e-7
  maximum_function_evaluations 100000
  


