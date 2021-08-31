#V3.30.17.00;_2021_06_11;_safe;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.3
#Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States.
#Foreign copyrights may apply. See copyright.txt for more information.
#_user_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_user_info_available_at:https://vlab.noaa.gov/group/stock-synthesis
#C growth parameters are estimated
#C spawner-recruitment bias adjustment Not tuned For optimality
#_data_and_control_files: data.ss // control.ss
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns (Growth Patterns, Morphs, Bio Patterns, GP are terms used interchangeably in SS)
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Platoon_within/between_stdev_ratio (no read if N_platoons=1)
#_Cond  1 #vector_platoon_dist_(-1_in_first_val_gives_normal_approx)
#
4 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
1 #_Nblock_Patterns
 1 #_blocks_per_pattern 
# begin and end years of blocks
 1970 1970
#
# controls for all timevary parameters 
1 #_time-vary parm bound check (1=warn relative to base parm bounds; 3=no bound check); Also see env (3) and dev (5) options to constrain with base bounds
#
# AUTOGEN
 0 0 0 0 0 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen time-varying parms of this category; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=P_base*exp(TVP*env(y));  2: P(y)=P_base+TVP*env(y);  3: P(y)=f(TVP,env_Zscore) w/ logit to stay in min-max;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=dev(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho;  5: like 4 with logit transform to stay in base min-max
#_DevLinks(more):  21-25 keep last dev for rest of years
#
#_Prior_codes:  0=none; 6=normal; 1=symmetric beta; 2=CASAL's beta; 3=lognormal; 4=lognormal with biascorr; 5=gamma
#
# setup for M, growth, wt-len, maturity, fecundity, (hermaphro), recr_distr, cohort_grow, (movement), (age error), (catch_mult), sex ratio 
#_NATMORT
3 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
 #_Age_natmort_by sex x growthpattern
 0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075	0.075
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0 #_Age(post-settlement)_for_L1;linear growth below this 
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
3 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
#_Age_Maturity by growth pattern
 0	0	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1 #_First_Mature_Age
3 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach for M, G, CV_G:  1- direct, no offset**; 2- male=fem_parm*exp(male_parm); 3: male=female*exp(parm) then old=young*exp(parm)
#_** in option 1, any male parameter with value = 0.0 and phase <0 is set equal to female parameter
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
# 0.05 0.15 0.075 0.1 0.8 0 -3 0 0 0 0 0 0 0 # NatM_uniform_Fem_GP_1
# Sex: 1  BioPattern: 1  Growth
 1 45 7.832525 36 10 6 2 0 0 0 0 0 0 0 # L_at_Amin_Fem_GP_1
 40 90 57.4842 70 10 6 4 0 0 0 0 0 0 0 # L_at_Amax_Fem_GP_1
 0.05 0.25 0.146478 0.15 0.8 6 4 0 0 0 0 0 0 0 # VonBert_K_Fem_GP_1
 0.05 0.25 0.08 0.08 0.8 0 -3 0 0 0 0 0 0 0 # CV_young_Fem_GP_1
 0.05 0.25 0.08 0.08 0.8 0 -3 0 0 0 0 0 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 -3 3 4.467e-05 4.467e-05 0.8 0 -3 0 0 0 0 0 0 0 # Wtlen_1_Fem_GP_1
 -3 4 2.793 2.793 0.8 0 -3 0 0 0 0 0 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 50 60 55 55 0.8 0 -3 0 0 0 0 0 0 0 # Mat50%_Fem_GP_1
 -3 3 -0.25 -0.25 0.8 0 -3 0 0 0 0 0 0 0 # Mat_slope_Fem_GP_1
 -3 3 1 1 0.8 0 -3 0 0 0 0 0 0 0 # Eggs/kg_inter_Fem_GP_1
 -3 3 1 1 0.8 0 -3 0 0 0 0 0 0 0 # Eggs/kg_slope_wt_Fem_GP_1
# Hermaphroditism
#  Recruitment Distribution  
# 0 0 0 0 0 0 -4 0 0 0 0 0 0 0 # RecrDist_GP_1
# 0 0 0 0 0 0 -4 0 0 0 0 0 0 0 # RecrDist_Area_1
# 0 0 0 0 0 0 -4 0 0 0 0 0 0 0 # RecrDist_month_1
#  Cohort growth dev base
 0.001 10 1 1 1 0 -1 0 0 0 0 0 0 0 # CohortGrowDev
#  Movement
#  Age Error from parameters
#  catch multiplier
#  fraction female, by GP
 1e-06 0.999999 0.5 0.5 0.5 0 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; Options: 1=NA; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
             3            31       8.59784          10.3            10             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
           0.2             1      		 1           0.7          0.05             1          4          0          0          0          0          0          0          0 # SR_BH_steep
             0             2           0.6           0.8           0.8             0         -4          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0             1             0         -4          0          0          0          0          0          0          0 # SR_regime
             0             0             0             0             0             0        -99          0          0          0          0          0          0          0 # SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1931 # first year of main recr_devs; early devs can preceed this era
2018 # last year of main recr_devs; forecast devs start in following year
2 #_recdev phase 
1 # (0/1) to read 13 advanced options
 0 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 -4 #_recdev_early_phase
 0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1900 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1900 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2020 #_last_yr_fullbias_adj_in_MPD
 2020 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
 -1 #_max_bias_adj_in_MPD (typical ~0.8; -3 sets all years to 0.0; -2 sets all non-forecast yrs w/ estimated recdevs to 1.0; -1 sets biasadj=1.0 for all yrs w/ recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -5 #min rec_dev
 5 #max rec_dev
 90 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#  1971R 1972R 1973R 1974R 1975R 1976R 1977R 1978R 1979R 1980R 1981R 1982R 1983R 1984R 1985R 1986R 1987R 1988R 1989R 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002F 2003F 2004F 2005F 2006F 2007F 2008F 2009F 2010F 2011F
#  0.1268 -0.0629442 0.0998014 -0.174095 0.0306484 0.714818 -0.0228752 0.00379775 0.261267 0.173626 0.0900049 -0.226622 -0.439888 -0.312088 0.393112 0.551725 0.21891 0.154932 -0.384786 0.596744 -0.682432 -0.273424 -0.829665 0.365024 -0.605267 0.455103 1.11072 -0.546499 -0.656469 0.171606 -0.301581 0 0 0 0 0 0 0 0 0 0
#1931R 1932R 1933R 1934R 1935R 1936R 1937R 1938R 1939R 1940R 1941R 1942R 1943R 1944R 1945R 1946R 1947R 1948R 1949R 1950R 1951R 1952R 1953R 1954R 1955R 1956R 1957R 1958R 1959R 1960R 1961R 1962R 1963R 1964R 1965R 1966R 1967R 1968R 1969R 1970R 1971R 1972R 1973R 1974R 1975R 1976R 1977R 1978R 1979R 1980R 1981R 1982R 1983R 1984R 1985R 1986R 1987R 1988R 1989R 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013R 2014R 2015R 2016R 2017R 2018R 2019R 2020R
1931R 0.18 1932R 0.18 1933R 0.18 1934R 0.18 1935R 0.18 1936R 0.18 1937R 0.18 1938R 0.18 1939R 0.18 1940R 0.18 1941R 0.18 1942R 0.18 1943R 0.18 1944R 0.18 1945R 0.18 1946R 0.18 1947R 0.18 1948R 0.18 1949R 0.18 1950R 0.18 1951R 0.18 1952R 0.18 1953R 0.18 1954R 0.18 1955R 0.18 1956R 0.18 1957R 0.18 1958R -0.0039228381609285 1959R -0.172398387171472 1960R 0.00326282149994594 1961R -0.198336440719912 1962R -0.253864582629862 1963R -0.243120043346885 1964R -0.223467105445491 1965R -0.344248644098131 1966R -0.0931219211204512 1967R -0.34932909533055 1968R -0.0135847490726654 1969R 0.156731373060646 1970R -0.178104536748327 1971R -0.162490308946776 1972R 0.0441802768574651 1973R 0.240153922819747 1974R -0.371647618286246 1975R -0.0633462586317292 1976R -0.0369130015635736 1977R -0.179536176219765 1978R -0.359568092631645 1979R -0.258504962186365 1980R -0.0787707289573609 1981R -0.678021823750179 1982R 0.201761491781513 1983R 0.607878702945065 1984R 0.312781111233818 1985R -0.199797361359587 1986R -1.61576749062559 1987R -0.199797361359587 1988R -1.55727128394399 1989R -0.0406466711156225 1990R -0.822393430927567 1991R -0.0356715364755087 1992R -0.470087691099498 1993R 0.048751713390046 1994R 0.437738196078709 1995R -1.15941077522104 1996R 0.583463105437491 1997R -0.136081546973479 1998R -0.166724613085564 1999R -0.210084006069862 2000R 0.144372822356849 2001R 0.152628803203868 2002R 0.0812840270608423 2003R 0.534171813720614 2004R 1.11883469630905 2005R 0.47416103854949 2006R 0.600025259439494 2007R 0.744745465755421 2008R 0.636158222423683 2009R 0.21825871211709 2010R 0.327557564357615 2011R 0.15775439105268 2012R 0.65374661552457 2013R 0.598052223613636 2014R 1.55624402526639 2015R 0.852434138962404 2016R 0.937060506303598 2017R 0.666123011125619 2018R 0.18 2019R 0.18 2020R 0.18
#
#Fishing Mortality info 
0.3 # F ballpark value in units of annual_F
-2001 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
2.95 # max F or harvest rate, depends on F_Method
# no additional F input needed for Fmethod 1
# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
# if Fmethod=3; read N iterations for tuning for Fmethod 3
4  # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms; for each fleet x season that has init_catch; nest season in fleet; count = 0
#_for unconstrained init_F, use an arbitrary initial catch and set lambda=0 for its logL
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
#
# F rates by fleet x season
# Yr:  1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# FISHERY 0 0.00211081 0.010609 0.0107037 0.0217063 0.0333334 0.0459509 0.0599453 0.0757167 0.107737 0.146876 0.162531 0.180868 0.202893 0.230365 0.266192 0.314644 0.338215 0.354481 0.356016 0.338877 0.238035 0.242891 0.250688 0.26355 0.283377 0.227156 0.238194 0.247552 0.252337 0.253174 0.0129829 0.0279253 0.038022 0.0447387 0.0493313 0.0527091 0.0554663 0.0579281 0.0602317 0.0624094
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm; 4=mirror with offset, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         6         1         0         0         0         1  #  CPUE_BT
         7         1         0         0         0         1  #  CPUE_BPT
         8         1         0         0         0         1  #  TrawSurvey
         9        1         0         0         0         1  #  TrawSurvey2yr
         10        1         0         0         0         1  #  TrawSurvey3yr
         11        1         0         0         0         1  #  TrawSurvey4yr
         12        1         0         0         0         1  #  TrawSurvey5yr
         13        1         0         0         0         1  #  TrawSurvey6_14yr
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -25            25      -9.2664              0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CPUE_BT(5)
           -25            25      -10.1797             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CPUE_BPT(6)
           -25            25      -8.37891             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_TrawSurveyAge2(7)
           -25            25      -8.37891             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_TrawSurveyAge2(7)
           -25            25      -8.37891             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_TrawSurveyAge2(7)
           -25            25      -8.37891             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_TrawSurveyAge2(7)
           -25            25      -8.37891             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_TrawSurveyAge2(7)
           -25            25      -8.37891             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_TrawSurveyAge2(7)
#_no timevary Q parameters
#
#_size_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for all sizes
#Pattern:_1; parm=2; logistic; with 95% width specification
#Pattern:_5; parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6; parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8; parm=8; New doublelogistic with smooth transitions and constant above Linf option
#Pattern:_9; parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in size
#Pattern:_27; parm=3+special; cubic spline 
#Pattern:_42; parm=2+special+3; // like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 0 0 0 0 # 1 FISHERY_BT
 0 0 0 0 # 2 FISHERY_BPT
 0 0 0 0 # 3 FISHERY_JP  approx LL select 5-30 = 1.0
 0 0 0 0 # 4 FISHERY_Rec1
 0 0 0 0 # 5 FISHERY_Rec2
 15 0 0 1 # 6 CPUE_BT
 15 0 0 2 # 7 CPUE_BPT
 0 0 0 0 # 8 TrawSurveyTotal
 0 0 0 0 # 9 TrawSurvey2yr
 0 0 0 0 # 10 TrawSurvey3yr
 0 0 0 0 # 11 TrawSurvey4yr
 0 0 0 0 # 12 TrawSurvey5yr
 0 0 0 0 # 13 TrawSurvey6_14yr

#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age
#Pattern:_42; parm=2+nages+1; // cubic spline; with 2 additional param for scaling (average over bin range)
#_Pattern Discard Male Special
 20 0 0 0 # 1 FISHERY_BT
 15 0 0 1 # 2 FISHERY_BPT
 11 0 0 0 # 3 FISHERY_JP approx LL select 5-30 = 1.0
 12 0 0 0 # 4 FISHERY_Rec1
 12 0 0 0 # 5 FISHERY_Rec2
 15 0 0 1 # 6 CPUE_BT
 15 0 0 1 # 7 CPUE_BPT
 20 0 0 0 # 8 TrawSurveyTotal
 11 0 0 0 # 9 TrawSurvey2yr
 11 0 0 0 # 10 TrawSurvey3yr
 11 0 0 0 # 11 TrawSurvey4yr
 11 0 0 0 # 12 TrawSurvey5yr
 11 0 0 0 # 13 TrawSurvey6_14yr
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   FISHERY_BT LenSelex
# full select at age 4 yr approx 32 cm
#            25            45            32          31.3             5             0          1          0          0          0          0          0          0          0  #  SizeSel_P1_FISHERY_BT(1)
#           -12             0           -10           -10             0             0         -5          0          0          0          0          0          0          0  #  SizeSel_P2_FISHERY_BT(1)
#           -10             9       3.39183             4             5             6          4          0          0          0          0          0          0          0  #  SizeSel_P3_FISHERY_BT(1)
#            -2            12       6.38008           6.1            10             6          4          0          0          0          0          0          0          0  #  SizeSel_P4_FISHERY_BT(1)
#           -10             9            -8            -8          1000             0         -5          0          0          0          0          0          0          0  #  SizeSel_P5_FISHERY_BT(1)
#           -16            10      -7.98845            -8            10             6          5          0          0          0          0          0          0          0  #  SizeSel_P6_FISHERY_BT(1)

# 2   FISHERY_BPT LenSelex
# full select at age 4 yr approx 32 cm
#            25            45            32          31.3             5             0          3          0          0          0          0          0          0          0  #  SizeSel_P1_FISHERY_BT(1)
#           -12             0           -10           -10             0             0         -5          0          0          0          0          0          0          0  #  SizeSel_P2_FISHERY_BT(1)
#           -10             9       3.39183             4             5             6          4          0          0          0          0          0          0          0  #  SizeSel_P3_FISHERY_BT(1)
#            -2            12       6.38008           6.1            10             6          4          0          0          0          0          0          0          0  #  SizeSel_P4_FISHERY_BT(1)
#           -10             9            -8            -8          1000             0         -5          0          0          0          0          0          0          0  #  SizeSel_P5_FISHERY_BT(1)
#           -16            10      -7.98845            -8            10             6          5          0          0          0          0          0          0          0  #  SizeSel_P6_FISHERY_BT(1)
## full select
#           -16            10             10            10            10             0          -5          0          0          0          0          0          0          0  #  SizeSel_P6_FISHERY_BT(1)

# 3   FISHERY_Rec1 Outside Harbours
#            25            45            32            32             5             0          3          0          0          0          0          0          0          0  #  SizeSel_P1_FISHERY_BT(1)
#           -12             0           -8             -8             0             0          -5          0          0          0          0          0          0          0  #  SizeSel_P2_FISHERY_BT(1)
#           -10             9            3              3             2             6          4          0          0          0          0          0          0          0  #  SizeSel_P3_FISHERY_BT(1)
#            -2            12            8              8             3             6          4          0          0          0          0          0          0          0  #  SizeSel_P4_FISHERY_BT(1)
#           -10             9            -8            -8          1000             0         -5          0          0          0          0          0          0          0  #  SizeSel_P5_FISHERY_BT(1)
#           -16            10      -7.98845            -8             3             6          5          0          0          0          0          0          0          0  #  SizeSel_P6_FISHERY_BT(1)
# 4   FISHERY_Rec2 Inside Harbours
#            25            45            28            28             3             6          3          0          0          0          0          0          0          0  #  SizeSel_P1_FISHERY_BT(1)
#           -12             0           -10           -10             0             0         -5          0          0          0          0          0          0          0  #  SizeSel_P2_FISHERY_BT(1)
#           -10             9             1             1             1             6          4          0          0          0          0          0          0          0  #  SizeSel_P3_FISHERY_BT(1)
#            -2            12           5.3           5.3             1             6          4          0          0          0          0          0          0          0  #  SizeSel_P4_FISHERY_BT(1)
#           -10             9            -8            -8          1000             0         -5          0          0          0          0          0          0          0  #  SizeSel_P5_FISHERY_BT(1)
#           -16            10           -10            -8            10             0         -5          0          0          0          0          0          0          0  #  SizeSel_P6_FISHERY_BT(1)

## FIXED PARAMETERS
# 3   FISHERY_Rec1 Outside Harbours
#            25            45       33.2189            32             5             0          -3          0          0          0          0          0          0          0  #  SizeSel_P1_FISHERY_BT(1)
#           -12             0      -8                  -8             0             0          -5          0          0          0          0          0          0          0  #  SizeSel_P2_FISHERY_BT(1)
#           -10             9       3.1213              3             2             0          -4          0          0          0          0          0          0          0  #  SizeSel_P3_FISHERY_BT(1)
#            -2            12       8.67436             8             3             0          -4          0          0          0          0          0          0          0  #  SizeSel_P4_FISHERY_BT(1)
#           -10             9      -8                  -8          1000             0          -5          0          0          0          0          0          0          0  #  SizeSel_P5_FISHERY_BT(1)
#           -16            10      -7.99597            -8             3             0          -5          0          0          0          0          0          0          0  #  SizeSel_P6_FISHERY_BT(1)
# 4   FISHERY_Rec2 Inside Harbours
#            25            45       27.7695            28             3             0          -3          0          0          0          0          0          0          0  #  SizeSel_P1_FISHERY_BT(1)
#           -12             0      -10                -10             0             0          -5          0          0          0          0          0          0          0  #  SizeSel_P2_FISHERY_BT(1)
#           -10             9        1.1773             1             1             0          -4          0          0          0          0          0          0          0  #  SizeSel_P3_FISHERY_BT(1)
#            -2            12        5.2451             5.3           1             0          -4          0          0          0          0          0          0          0  #  SizeSel_P4_FISHERY_BT(1)
#           -10             9       -8                 -8          1000             0          -5          0          0          0          0          0          0          0  #  SizeSel_P5_FISHERY_BT(1)
#           -16            10      -10                 -8            10             0          -5          0          0          0          0          0          0          0  #  SizeSel_P6_FISHERY_BT(1)

# Trawl survey
#            13            45            20            20             3             0          3          0          0          0          0          0          0          0  #  SizeSel_P1_FISHERY_BT(1)
#           -12             0           -10           -10             0             0          5          0          0          0          0          0          0          0  #  SizeSel_P2_FISHERY_BT(1)
#           -10             9       3.39183             4             5             6          4          0          0          0          0          0          0          0  #  SizeSel_P3_FISHERY_BT(1)
#            -2            12       6.38008           6.1            10             6          4          0          0          0          0          0          0          0  #  SizeSel_P4_FISHERY_BT(1)
#           -10             9            -8            -8            10             6         -6          0          0          0          0          0          0          0  #  SizeSel_P5_FISHERY_BT(1)
#           -16            10      -7.98845            -8            10             6          7          0          0          0          0          0          0          0  #  SizeSel_P6_FISHERY_BT(1)

# 5   CPUE_BT LenSelex
# 6   CPUE_BPT LenSelex
# 7   TrawSurveyAge2 LenSelex
# 8   TrawSurveyAge3 LenSelex
# 1   FISHERY_BT AgeSelex
# 2   FISHERY_BPT AgeSelex
# 3   FISHERY_Rec AgeSelex

# 1   BT and BPT AgeSelex
             3             8             4.49355             5          1              6          1          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
           -10             5              -5           -5          99              0          -5          0          0          0          0          0          0          0  #  AgeSel_P2_Survey_Tag(4)
             0             5             0.49028          0.5          1              6          4          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
             1             10            2.69634              3          2            6          6          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
           -15             10            -10           -10         99              0         -1          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
## constrain
#           -10             10             0             0          99              0          3          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
           -10             10             -0.0283879            -0.313788          2              6          3          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)


## 3 JP catch - approximate LL select
# 4   JP AgeSelex fix at 1.0 for ages 5-30
             0             5             5             5           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
             0             5            100            100           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P2_Survey_Tag(4)
# 3   FISHERY_Rec1 Outside Harbours
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
             0.1           20             3.4            5           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
             0.1           20            1.2            100           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P2_Survey_Tag(4)
# 4   FISHERY_Rec2 Inside Harbours
             0.1           20             2.7            5           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
             0.1           20            2.2            100           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P2_Survey_Tag(4)
# 5   CPUE_BT AgeSelex
# 6   CPUE_BPT AgeSelex
# Trawl Survey - all 
             1             15             3             3          99              0          1          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
           -10             5              -5           -5          99              0          -5          0          0          0          0          0          0          0  #  AgeSel_P2_Survey_Tag(4)
             0             10             2             2          99              0          4          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
             1             10           6.5          6.5           99              0          6          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
           -15             10            -10           -10         99              0         -1          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
           -10             10             0             0          99              0          3          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)

# Trawl Survey 2yr 
             0             6             2            2           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
             0             6             2            2           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P2_Survey_Tag(4)

# Trawl Survey 3yr 
             0             6             3            3           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
             0             6             3            3           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P2_Survey_Tag(4)
# Trawl Survey 4yr 
             0             6             4            4           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
             0             6             4            4           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P2_Survey_Tag(4)
# Trawl Survey 5yr 
             0             6             5            5           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
             0             6             5            5           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P2_Survey_Tag(4)
# Trawl Survey 5yr 
             0             20             6            6           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P1_Survey_Tag(4)
             0             20             14           14           0.5             0         -1          0          0          0          0          0          0          0  #  AgeSel_P2_Survey_Tag(4)
#_No_Dirichlet parameters
#_no timevary selex parameters
#
0   #  use 2D_AR1 selectivity(0/1)
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read and autogen if tag data exist; 1=read
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# no timevary parameters
#
#
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
 -9999   1    0  # terminator
#
4 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 3 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
 1 2 2 1 1
 4 2 2 1 1
 4 2 3 1 1
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  0 0 0 0 #_CPUE/survey:_1
#  1 1 1 1 #_CPUE/survey:_2
#  1 1 1 1 #_CPUE/survey:_3
#  1 1 1 1 #_lencomp:_1
#  1 1 1 1 #_lencomp:_2
#  0 0 0 0 #_lencomp:_3
#  1 1 1 1 #_agecomp:_1
#  1 1 1 1 #_agecomp:_2
#  0 0 0 0 #_agecomp:_3
#  1 1 1 1 #_size-age:_1
#  1 1 1 1 #_size-age:_2
#  0 0 0 0 #_size-age:_3
#  1 1 1 1 #_init_equ_catch1
#  1 1 1 1 #_init_equ_catch2
#  1 1 1 1 #_init_equ_catch3
#  1 1 1 1 #_recruitments
#  1 1 1 1 #_parameter-priors
#  1 1 1 1 #_parameter-dev-vectors
#  1 1 1 1 #_crashPenLambda
#  0 0 0 0 # F_ballpark_lambda
1 # (0/1/2) read specs for more stddev reporting: 0 = skip, 1 = read specs for reporting stdev for selectivity, size, and numbers, 2 = add options for M,Dyn. Bzero, SmryBio
 1 1 -1 5 # Selectivity: (1) 0 to skip or fleet, (2) 1=len/2=age/3=combined, (3) year, (4) N selex bins; NOTE: combined reports in age bins
 1 5 # Growth: (1) 0 to skip or growth pattern, (2) growth ages; NOTE: does each sex
 1 -1 5 # Numbers-at-age: (1) 0 or area(-1 for all), (2) year, (3) N ages;  NOTE: sums across morphs
 5 15 25 35 43 # vector with selex std bins (-1 in first bin to self-generate)
 1 2 14 26 40 # vector with growth std ages picks (-1 in first bin to self-generate)
 1 2 14 26 40 # vector with NatAge std ages (-1 in first bin to self-generate)
999

