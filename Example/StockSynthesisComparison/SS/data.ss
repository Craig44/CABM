#V3.30.17.00;_2021_06_11;_safe;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.3
#Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States.
#Foreign copyrights may apply. See copyright.txt for more information.
#_user_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_user_info_available_at:https://vlab.noaa.gov/group/stock-synthesis
#_Start_time: Fri Jun 11 08:38:55 2021
#_Number_of_datafiles: 1
#C data file for simple example
#_observed data: 
#V3.30.17.00;_2021_06_11;_safe;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.3
#Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States.
#Foreign copyrights may apply. See copyright.txt for more information.
1931 #_StartYr
2020 #_EndYr
1 #_Nseas
12 #_months/season
2 #_Nsubseasons (even number, minimum is 2)
1 #_spawn_month
1 #_Ngenders: 1, 2, -1  (use -1 for 1 sex setup with SSB multiplied by female_frac parameter)
80 #_Nages=accumulator age, first age is always age 0
1 #_Nareas
13 #_Nfleets (including surveys)
#_fleet_type: 1=catch fleet; 2=bycatch only fleet; 3=survey; 4=ignore 
#_sample_timing: -1 for fishing fleet to use season-long catch-at-age for observations, or 1 to use observation month;  (always 1 for surveys)
#_fleet_area:  area the fleet/survey operates in 
#_units of catch:  1=bio; 2=num (ignored for surveys; their units read later)
#_catch_mult: 0=no; 1=yes
#_rows are fleets
#_fleet_type fishery_timing area catch_units need_catch_mult fleetname
 1 -1 1 1 0 FISHERY_BT  # 1
 1 -1 1 1 0 FISHERY_BPT  # 2
 1 -1 1 1 0 FISHERY_JP  # 3
 1 -1 1 1 0 FISHERY_Rec1  # 4  Oustide Harbour catch
 1 -1 1 1 0 FISHERY_Rec2  # 5  Harbour catch
 3 1 1 2 0 CPUE_BT  # 6
 3 1 1 2 0 CPUE_BPT   # 7
 3 1 1 2 0 TrawSurveyTotal  # 8
 3 1 1 2 0 TrawSurvey2yr  # 9
 3 1 1 2 0 TrawSurvey3yr  # 10
 3 1 1 2 0 TrawSurvey4yr  # 11
 3 1 1 2 0 TrawSurvey5yr  # 12
 3 1 1 2 0 TrawSurvey6_14yr  # 13
#Bycatch_fleet_input_goes_next
#a:  fleet index
#b:  1=include dead bycatch in total dead catch for F0.1 and MSY optimizations and forecast ABC; 2=omit from total catch for these purposes (but still include the mortality)
#c:  1=Fmult scales with other fleets; 2=bycatch F constant at input value; 3=bycatch F from range of years
#d:  F or first year of range
#e:  last year of range
#f:  not used
# a   b   c   d   e   f 
#_Catch data: yr, seas, fleet, catch, catch_se
#_catch_se:  standard error of log(catch)
#_NOTE:  catch data is ignored for survey fleets
## BT
-999	1	1	0	0.01
1931	1	1	168.0	0.01
1932	1	1	190.8	0.01
1933	1	1	255.6	0.01
1934	1	1	228.0	0.01
1935	1	1	129.6	0.01
1936	1	1	123.6	0.01
1937	1	1	102.0	0.01
1938	1	1	106.8	0.01
1939	1	1	85.2	0.01
1940	1	1	91.2	0.01
1941	1	1	74.4	0.01
1942	1	1	68.4	0.01
1943	1	1	90.0	0.01
1944	1	1	82.8	0.01
1945	1	1	148.8	0.01
1946	1	1	292.8	0.01
1947	1	1	301.2	0.01
1948	1	1	258.0	0.01
1949	1	1	332.4	0.01
1950	1	1	381.6	0.01
1951	1	1	436.8	0.01
1952	1	1	433.2	0.01
1953	1	1	1348.8	0.01
1954	1	1	1311.6	0.01
1955	1	1	1442.4	0.01
1956	1	1	1395.6	0.01
1957	1	1	1766.4	0.01
1958	1	1	1353.6	0.01
1959	1	1	1336.8	0.01
1960	1	1	1442.4	0.01
1961	1	1	1413.6	0.01
1962	1	1	1622.4	0.01
1963	1	1	1747.2	0.01
1964	1	1	1531.2	0.01
1965	1	1	1418.4	0.01
1966	1	1	2197.2	0.01
1967	1	1	1772.4	0.01
1968	1	1	1789.2	0.01
1969	1	1	1612.8	0.01
1970	1	1	1905.6	0.01
1971	1	1	2222.4	0.01
1972	1	1	1176.6	0.01
1973	1	1	1822.8	0.01
1974	1	1	2604.0	0.01
1975	1	1	2530.2	0.01
1976	1	1	2556.5	0.01
1977	1	1	1418.8	0.01
1978	1	1	1041.6	0.01
1979	1	1	831.4	0.01
1980	1	1	794.2	0.01
1981	1	1	756.7	0.01
1982	1	1	632.6	0.01
1983	1	1	435.4	0.01
1984	1	1	414.0	0.01
1985	1	1	371.0	0.01
1986	1	1	438.7	0.01
1987	1	1	294.7	0.01
1988	1	1	616.4	0.01
1989	1	1	839.9	0.01
1990	1	1	1010.2	0.01
1991	1	1	1091.4	0.01
1992	1	1	1100.6	0.01
1993	1	1	1382.1	0.01
1994	1	1	1162.5	0.01
1995	1	1	970.8	0.01
1996	1	1	1128.3	0.01
1997	1	1	1484.8	0.01
1998	1	1	1505.1	0.01
1999	1	1	1314.8	0.01
2000	1	1	1241.3	0.01
2001	1	1	1235.1	0.01
2002	1	1	1389.2	0.01
2003	1	1	1413.6	0.01
2004	1	1	1592.1	0.01
2005	1	1	1255.2	0.01
2006	1	1	1154.8	0.01
2007	1	1	1294.6	0.01
2008	1	1	1230.8	0.01
2009	1	1	1331.6	0.01
2010	1	1	1223.0	0.01
2011	1	1	1289.1	0.01
2012	1	1	1487.6	0.01
2013	1	1	1464.1	0.01
2014	1	1	1402.5	0.01
2015	1	1	1399.2	0.01
2016	1	1	1460.8	0.01
2017	1	1	1467.4	0.01
2018	1	1	1415.7	0.01
2019	1	1	1421.2	0.01
2020	1	1	1480.6	0.01 ## E:\trophia\SNA8_2019\SSmodel2021\Catch history v3 2021 update.xlsx
##BPT
-999	1	2	0	0.01
1931	1	2	0.0	0.01
1932	1	2	0.0	0.01
1933	1	2	0.0	0.01
1934	1	2	0.0	0.01
1935	1	2	0.0	0.01
1936	1	2	0.0	0.01
1937	1	2	0.0	0.01
1938	1	2	0.0	0.01
1939	1	2	0.0	0.01
1940	1	2	0.0	0.01
1941	1	2	0.0	0.01
1942	1	2	0.0	0.01
1943	1	2	0.0	0.01
1944	1	2	0.0	0.01
1945	1	2	0.0	0.01
1946	1	2	0.0	0.01
1947	1	2	0.0	0.01
1948	1	2	0.0	0.01
1949	1	2	0.0	0.01
1950	1	2	0.0	0.01
1951	1	2	0.0	0.01
1952	1	2	0.0	0.01
1953	1	2	0.0	0.01
1954	1	2	0.0	0.01
1955	1	2	0.0	0.01
1956	1	2	0.0	0.01
1957	1	2	0.0	0.01
1958	1	2	0.0	0.01
1959	1	2	0.0	0.01
1960	1	2	0.0	0.01
1961	1	2	0.0	0.01
1962	1	2	0.0	0.01
1963	1	2	0.0	0.01
1964	1	2	0.0	0.01
1965	1	2	0.0	0.01
1966	1	2	0.0	0.01
1967	1	2	0.0	0.01
1968	1	2	0.0	0.01
1969	1	2	0.0	0.01
1970	1	2	0.0	0.01
1971	1	2	0.0	0.01
1972	1	2	1176.6	0.01
1973	1	2	1822.8	0.01
1974	1	2	2604.0	0.01
1975	1	2	2530.2	0.01
1976	1	2	3834.7	0.01
1977	1	2	3310.4	0.01
1978	1	2	4166.4	0.01
1979	1	2	3325.4	0.01
1980	1	2	3176.6	0.01
1981	1	2	3026.9	0.01
1982	1	2	2530.6	0.01
1983	1	2	1741.4	0.01
1984	1	2	1656.0	0.01
1985	1	2	1484.2	0.01
1986	1	2	1754.9	0.01
1987	1	2	687.6	0.01
1988	1	2	924.7	0.01
1989	1	2	839.9	0.01
1990	1	2	695.9	0.01
1991	1	2	733.5	0.01
1992	1	2	504.3	0.01
1993	1	2	315.2	0.01
1994	1	2	533.7	0.01
1995	1	2	608.8	0.01
1996	1	2	585.5	0.01
1997	1	2	289.5	0.01
1998	1	2	242.8	0.01
1999	1	2	484.8	0.01
2000	1	2	523.1	0.01
2001	1	2	559.0	0.01
2002	1	2	345.5	0.01
2003	1	2	300.2	0.01
2004	1	2	241.6	0.01
2005	1	2	574.1	0.01
2006	1	2	422.6	0.01
2007	1	2	165.1	0.01
2008	1	2	203.6	0.01
2009	1	2	147.9	0.01
2010	1	2	185.0	0.01
2011	1	2	155.2	0.01
2012	1	2	8.4	0.01
2013	1	2	0.0	0.01
2014	1	2	0.0	0.01
2015	1	2	0.0	0.01
2016	1	2	0.0	0.01
2017	1	2	0.0	0.01
2018	1	2	0.0	0.01
2019	1	2	0.0	0.01
2020	1	2	0.0	0.01
## JP catch (LL not reported)
-999	1	3	0	0.01
1931	1	3	0.0	0.01
1932	1	3	0.0	0.01
1933	1	3	0.0	0.01
1934	1	3	0.0	0.01
1935	1	3	0.0	0.01
1936	1	3	0.0	0.01
1937	1	3	0.0	0.01
1938	1	3	0.0	0.01
1939	1	3	0.0	0.01
1940	1	3	0.0	0.01
1941	1	3	0.0	0.01
1942	1	3	0.0	0.01
1943	1	3	0.0	0.01
1944	1	3	0.0	0.01
1945	1	3	0.0	0.01
1946	1	3	0.0	0.01
1947	1	3	0.0	0.01
1948	1	3	0.0	0.01
1949	1	3	0.0	0.01
1950	1	3	0.0	0.01
1951	1	3	0.0	0.01
1952	1	3	0.0	0.01
1953	1	3	0.0	0.01
1954	1	3	0.0	0.01
1955	1	3	0.0	0.01
1956	1	3	0.0	0.01
1957	1	3	0.0	0.01
1958	1	3	0.0	0.01
1959	1	3	0.0	0.01
1960	1	3	0.0	0.01
1961	1	3	0.0	0.01
1962	1	3	0.0	0.01
1963	1	3	0.0	0.01
1964	1	3	0.0	0.01
1965	1	3	2000.0	0.01
1966	1	3	2000.0	0.01
1967	1	3	2000.0	0.01
1968	1	3	2000.0	0.01
1969	1	3	2000.0	0.01
1970	1	3	2000.0	0.01
1971	1	3	2000.0	0.01
1972	1	3	2000.0	0.01
1973	1	3	2000.0	0.01
1974	1	3	2000.0	0.01
1975	1	3	749.0	0.01
1976	1	3	1127.0	0.01
1977	1	3	1104.0	0.01
1978	1	3	0.0	0.01
1979	1	3	0.0	0.01
1980	1	3	0.0	0.01
1981	1	3	0.0	0.01
1982	1	3	0.0	0.01
1983	1	3	0.0	0.01
1984	1	3	0.0	0.01
1985	1	3	0.0	0.01
1986	1	3	0.0	0.01
1987	1	3	0.0	0.01
1988	1	3	0.0	0.01
1989	1	3	0.0	0.01
1990	1	3	0.0	0.01
1991	1	3	0.0	0.01
1992	1	3	0.0	0.01
1993	1	3	0.0	0.01
1994	1	3	0.0	0.01
1995	1	3	0.0	0.01
1996	1	3	0.0	0.01
1997	1	3	0.0	0.01
1998	1	3	0.0	0.01
1999	1	3	0.0	0.01
2000	1	3	0.0	0.01
2001	1	3	0.0	0.01
2002	1	3	0.0	0.01
2003	1	3	0.0	0.01
2004	1	3	0.0	0.01
2005	1	3	0.0	0.01
2006	1	3	0.0	0.01
2007	1	3	0.0	0.01
2008	1	3	0.0	0.01
2009	1	3	0.0	0.01
2010	1	3	0.0	0.01
2011	1	3	0.0	0.01
2012	1	3	0.0	0.01
2013	1	3	0.0	0.01
2014	1	3	0.0	0.01
2015	1	3	0.0	0.01
2016	1	3	0.0	0.01
2017	1	3	0.0	0.01
2018	1	3	0.0	0.01
2019	1	3	0.0	0.01
2020	1	3	0.0	0.01
## 
-999 1 4 0 0.01
1931	1	4	45	0.01
1932	1	4	47	0.01
1933	1	4	50	0.01
1934	1	4	52	0.01
1935	1	4	54	0.01
1936	1	4	56	0.01
1937	1	4	59	0.01
1938	1	4	61	0.01
1939	1	4	63	0.01
1940	1	4	65	0.01
1941	1	4	68	0.01
1942	1	4	70	0.01
1943	1	4	72	0.01
1944	1	4	75	0.01
1945	1	4	77	0.01
1946	1	4	79	0.01
1947	1	4	81	0.01
1948	1	4	84	0.01
1949	1	4	86	0.01
1950	1	4	88	0.01
1951	1	4	91	0.01
1952	1	4	93	0.01
1953	1	4	95	0.01
1954	1	4	97	0.01
1955	1	4	100	0.01
1956	1	4	102	0.01
1957	1	4	104	0.01
1958	1	4	106	0.01
1959	1	4	109	0.01
1960	1	4	111	0.01
1961	1	4	113	0.01
1962	1	4	116	0.01
1963	1	4	118	0.01
1964	1	4	120	0.01
1965	1	4	122	0.01
1966	1	4	125	0.01
1967	1	4	127	0.01
1968	1	4	129	0.01
1969	1	4	131	0.01
1970	1	4	134	0.01
1971	1	4	136	0.01
1972	1	4	138	0.01
1973	1	4	141	0.01
1974	1	4	143	0.01
1975	1	4	145	0.01
1976	1	4	147	0.01
1977	1	4	150	0.01
1978	1	4	152	0.01
1979	1	4	154	0.01
1980	1	4	156	0.01
1981	1	4	159	0.01
1982	1	4	161	0.01
1983	1	4	163	0.01
1984	1	4	166	0.01
1985	1	4	168	0.01
1986	1	4	170	0.01
1987	1	4	172	0.01
1988	1	4	175	0.01
1989	1	4	177	0.01
1990	1	4	179	0.01
1991	1	4	179	0.01
1992	1	4	179	0.01
1993	1	4	179	0.01
1994	1	4	179	0.01
1995	1	4	179	0.01
1996	1	4	179	0.01
1997	1	4	173	0.01
1998	1	4	173	0.01
1999	1	4	181	0.01
2000	1	4	178	0.01
2001	1	4	186	0.01
2002	1	4	188	0.01
2003	1	4	190	0.01
2004	1	4	191	0.01
2005	1	4	191	0.01
2006	1	4	191	0.01
2007	1	4	195	0.01
2008	1	4	227	0.01
2009	1	4	281	0.01
2010	1	4	336	0.01
2011	1	4	400	0.01
2012	1	4	473	0.01
2013	1	4	518	0.01
2014	1	4	556	0.01
2015	1	4	591	0.01
2016	1	4	620	0.01
2017	1	4	647	0.01
2018	1	4	669	0.01
2019	1	4	701	0.01
2020	1	4	701	0.01 ## equivalent to 2019 
-999 1 5 0 0.01
1931	1	5	15	0.01
1932	1	5	16	0.01
1933	1	5	17	0.01
1934	1	5	17	0.01
1935	1	5	18	0.01
1936	1	5	19	0.01
1937	1	5	20	0.01
1938	1	5	20	0.01
1939	1	5	21	0.01
1940	1	5	22	0.01
1941	1	5	23	0.01
1942	1	5	23	0.01
1943	1	5	24	0.01
1944	1	5	25	0.01
1945	1	5	26	0.01
1946	1	5	26	0.01
1947	1	5	27	0.01
1948	1	5	28	0.01
1949	1	5	29	0.01
1950	1	5	29	0.01
1951	1	5	30	0.01
1952	1	5	31	0.01
1953	1	5	32	0.01
1954	1	5	32	0.01
1955	1	5	33	0.01
1956	1	5	34	0.01
1957	1	5	35	0.01
1958	1	5	35	0.01
1959	1	5	36	0.01
1960	1	5	37	0.01
1961	1	5	38	0.01
1962	1	5	39	0.01
1963	1	5	39	0.01
1964	1	5	40	0.01
1965	1	5	41	0.01
1966	1	5	42	0.01
1967	1	5	42	0.01
1968	1	5	43	0.01
1969	1	5	44	0.01
1970	1	5	45	0.01
1971	1	5	45	0.01
1972	1	5	46	0.01
1973	1	5	47	0.01
1974	1	5	48	0.01
1975	1	5	48	0.01
1976	1	5	49	0.01
1977	1	5	50	0.01
1978	1	5	51	0.01
1979	1	5	51	0.01
1980	1	5	52	0.01
1981	1	5	53	0.01
1982	1	5	54	0.01
1983	1	5	54	0.01
1984	1	5	55	0.01
1985	1	5	56	0.01
1986	1	5	57	0.01
1987	1	5	57	0.01
1988	1	5	58	0.01
1989	1	5	59	0.01
1990	1	5	60	0.01
1991	1	5	60	0.01
1992	1	5	60	0.01
1993	1	5	60	0.01
1994	1	5	60	0.01
1995	1	5	60	0.01
1996	1	5	60	0.01
1997	1	5	58	0.01
1998	1	5	58	0.01
1999	1	5	60	0.01
2000	1	5	59	0.01
2001	1	5	62	0.01
2002	1	5	63	0.01
2003	1	5	63	0.01
2004	1	5	64	0.01
2005	1	5	64	0.01
2006	1	5	64	0.01
2007	1	5	65	0.01
2008	1	5	76	0.01
2009	1	5	94	0.01
2010	1	5	112	0.01
2011	1	5	133	0.01
2012	1	5	158	0.01
2013	1	5	173	0.01
2014	1	5	185	0.01
2015	1	5	197	0.01
2016	1	5	207	0.01
2017	1	5	216	0.01
2018	1	5	223	0.01
2019	1	5	234	0.01
2020	1	5	234	0.01 ## equivalent to 2019 
-9999 0 0 0 0
#
 #_CPUE_and_surveyabundance_observations
#_Units:  0=numbers; 1=biomass; 2=F; 30=spawnbio; 31=recdev; 32=spawnbio*recdev; 33=recruitment; 34=depletion(&see Qsetup); 35=parm_dev(&see Qsetup)
#_Errtype:  -1=normal; 0=lognormal; >0=T
#_SD_Report: 0=no sdreport; 1=enable sdreport
#_Fleet Units Errtype SD_Report
1 1 0 0 # FISHERY_BT
2 1 0 0 # FISHERY_BPT
3 1 0 0 # FISHERY_BPT
4 1 0 0 # FISHERY_Rec1
5 1 0 0 # FISHERY_Rec2
6 1 0 1 # CPUE_BT
7 0 0 0 # CPUE_BPT
8 1 0 0 # TrawSurvey Total biomass
9 0 0 0 # TrawSurvey 2yr number
10 0 0 0 # TrawSurvey 3yr number
11 0 0 0 # TrawSurvey 4yr number
12 0 0 0 # TrawSurvey 5yr number
13 0 0 0 # TrawSurvey 6_14yr number
#_yr month fleet obs stderr
1990 3.4 6 9505 0.18 #_ Survey_Tag
2002 3.4 6 10441 0.12 #_ Survey_Tag
# updated recent CPUE indices and CV
# updated to include 2019/20 
# E:\trophia\SNA8_2019\CPUE2021\CPUETop5\tables, E:\trophia\SNA8_2019\SSmodel2021\cpue indices
1997	1	6	0.412	0.163
1998	1	6	0.376	0.156
1999	1	6	0.422	0.137
2000	1	6	0.513	0.124
2001	1	6	0.370	0.149
2002	1	6	0.545	0.127
2003	1	6	0.515	0.133
2004	1	6	0.390	0.15
2005	1	6	0.535	0.137
2006	1	6	0.668	0.176
2007	1	6	0.452	0.15
2008	1	6	0.773	0.138
2009	1	6	0.693	0.128
2010	1	6	0.931	0.138
2011	1	6	1.310	0.142
2012	1	6	1.374	0.133
2013	1	6	1.348	0.132
2014	1	6	1.647	0.126
2015	1	6	1.897	0.136
2016	1	6	1.848	0.139
2017	1	6	1.649	0.159
2018	1	6	1.839	0.193
2019	1	6	2.084	0.191
2020	1	6	1.411	0.171  ## drop in 2020 linked to single vessel
## includes process error = 0.2
1974 3 7 1 0.33 #_ CPUE_BPT
1975 3 7 0.85 0.329 #_ CPUE_BPT
1976 3 7 0.86 0.328 #_ CPUE_BPT
1977 3 7 0.46 0.33 #_ CPUE_BPT
1978 3 7 0.59 0.336 #_ CPUE_BPT
1979 3 7 0.37 0.335 #_ CPUE_BPT
1980 3 7 0.25 0.32 #_ CPUE_BPT
1981 3 7 0.21 0.343 #_ CPUE_BPT
1982 3 7 0.2 0.35 #_ CPUE_BPT
1983 3 7 0.15 0.333 #_ CPUE_BPT
1984 3 7 0.11 0.382 #_ CPUE_BPT
1985 3 7 0.13 0.354 #_ CPUE_BPT
1986 3 7 0.07 0.343 #_ CPUE_BPT
1987 3 7 0.13 0.354 #_ CPUE_BPT
1988 3 7 0.16 0.45 #_ CPUE_BPT
1989 3 7 0.25 0.4 #_ CPUE_BPT
1990 3 7 0.67 0.499 #_ CPUE_BPT
1991 3 7 0.36 0.422 #_ CPUE_BPT
## trawl survey biomass
## updated biomass indices from E:\trophia\SNA8_2019\data\TrawlSurvey\WCNI SNA core area time series for Adam ALL YEARS 2020.txt
1990	3.5	8	1234	0.278	#_	TrawSurveyTotal
1992	3.5	8	576.5	0.304	#_	TrawSurveyTotal
1995	3.5	8	409.8	0.243	#_	TrawSurveyTotal
1997	3.5	8	920.4	0.282	#_	TrawSurveyTotal
2000	3.5	8	869.9	0.167	#_	TrawSurveyTotal
2019	3.5	8	3821	0.223	#_	TrawSurveyTotal
2020	3.5	8	7999	0.197	#_	TrawSurveyTotal
2021	3.5	8	11796.7	0.203	#_	TrawSurveyTotal
## trawl survey number (1000s) age 2yr = 1+
1990	3.5	9	50.1	0.353	#_	TrawlSurvey2yr
1992	3.5	9	36.3	0.384	#_	TrawlSurvey2yr
1995	3.5	9	344.9	0.314	#_	TrawlSurvey2yr
1997	3.5	9	996	0.388	#_	TrawlSurvey2yr
2000	3.5	9	876.3	0.341	#_	TrawlSurvey2yr
2019	3.5	9	658.2	0.321	#_	TrawlSurvey2yr
2020	3.5	9	902	0.357	#_	TrawlSurvey2yr
2021	3.5	9	820.4	0.233	#_	TrawlSurvey2yr
## trawl survey number (1000s) age 3yr = 2+						
1990	3.5	10	407.7	0.306	#_	TrawlSurvey3yr
1992	3.5	10	344.9	0.192	#_	TrawlSurvey3yr
1995	3.5	10	192.7	0.319	#_	TrawlSurvey3yr
1997	3.5	10	244	0.204	#_	TrawlSurvey3yr
2000	3.5	10	109.2	0.199	#_	TrawlSurvey3yr
2019	3.5	10	1141.9	0.211	#_	TrawlSurvey3yr
2020	3.5	10	727.6	0.217	#_	TrawlSurvey3yr
2021	3.5	10	1188.6	0.213	#_	TrawlSurvey3yr
## trawl survey number (1000s) age 4yr = 3+						
1990	3.5	11	562	0.346	#_	TrawlSurvey4yr	
1992	3.5	11	38.7	0.393	#_	TrawlSurvey4yr	
1995	3.5	11	247.3	0.242	#_	TrawlSurvey4yr	
1997	3.5	11	262.6	0.251	#_	TrawlSurvey4yr	
2000	3.5	11	601.8	0.15	#_	TrawlSurvey4yr	
2019	3.5	11	731.6	0.18	#_	TrawlSurvey4yr	
2020	3.5	11	2353.1	0.129	#_	TrawlSurvey4yr	
2021	3.5	11	1410.7	0.156	#_	TrawlSurvey4yr	
## trawl survey number (1000s) age 5yr = 4+						
1990	3.5	12	343.9	0.285	#_	TrawlSurvey5yr
1992	3.5	12	136.6	0.329	#_	TrawlSurvey5yr
1995	3.5	12	36.9	0.305	#_	TrawlSurvey5yr
1997	3.5	12	60.7	0.454	#_	TrawlSurvey5yr
2000	3.5	12	157.2	0.259	#_	TrawlSurvey5yr
2019	3.5	12	683.2	0.182	#_	TrawlSurvey5yr
2020	3.5	12	809.4	0.204	#_	TrawlSurvey5yr
2021	3.5	12	2404.4	0.174	#_	TrawlSurvey5yr
## trawl survey number (1000s) age 6-14yr
1990	3.5	13	255.6	0.381	#_	TrawlSurvey6_14yr
1992	3.5	13	215.4	0.353	#_	TrawlSurvey6_14yr
1995	3.5	13	57.7	0.461	#_	TrawlSurvey6_14yr
1997	3.5	13	241.1	0.525	#_	TrawlSurvey6_14yr
2000	3.5	13	113.9	0.433	#_	TrawlSurvey6_14yr
2019	3.5	13	925	0.372	#_	TrawlSurvey6_14yr
2020	3.5	13	3246.8	0.279	#_	TrawlSurvey6_14yr
2021	3.5	13	3802.8	0.279	#_	TrawlSurvey6_14yr
-9999 1 1 1 1 # terminator for survey observations 
#
0 #_N_fleets_with_discard
#_discard_units (1=same_as_catchunits(bio/num); 2=fraction; 3=numbers)
#_discard_errtype:  >0 for DF of T-dist(read CV below); 0 for normal with CV; -1 for normal with se; -2 for lognormal; -3 for trunc normal with CV
# note: only enter units and errtype for fleets with discard 
# note: discard data is the total for an entire season, so input of month here must be to a month in that season
#_Fleet units errtype
# -9999 0 0 0.0 0.0 # terminator for discard data 
#
0 #_use meanbodysize_data (0/1)
#_COND_0 #_DF_for_meanbodysize_T-distribution_like
# note:  type=1 for mean length; type=2 for mean body weight 
#_yr month fleet part type obs stderr
#  -9999 0 0 0 0 0 0 # terminator for mean body size data 
#
# set up population length bin structure (note - irrelevant if not using size data and using empirical wtatage
2 # length bin method: 1=use databins; 2=generate from binwidth,min,max below; 3=read vector
1 # binwidth for population size comp 
1 # minimum size in the population (lower edge of first bin and size at age 0.00) 
80 # maximum size in the population (lower edge of last bin) 
1 # use length composition data (0/1)
#_mintailcomp: upper and lower distribution for females and males separately are accumulated until exceeding this level.
#_addtocomp:  after accumulation of tails; this value added to all bins
#_combM+F: males and females treated as combined gender below this bin number 
#_compressbins: accumulate upper tail by this number of bins; acts simultaneous with mintailcomp; set=0 for no forced accumulation
#_Comp_Error:  0=multinomial, 1=dirichlet
#_ParmSelect:  parm number for dirichlet
#_minsamplesize: minimum sample size; set to 1 to match 3.24, minimum value is 0.001
#
#_mintailcomp addtocomp combM+F CompressBins CompError ParmSelect minsamplesize
0.0005 1e-07 0 0 0 0 1 #_fleet:1_FISHERY_BT
0.0005 1e-07 0 0 0 0 1 #_fleet:2_FISHERY_BPT
0.0005 1e-07 0 0 0 0 1 #_fleet:3_FISHERY_JP
0.0005 1e-07 0 0 0 0 1 #_fleet:4_FISHERY_Rec1
0.0005 1e-07 0 0 0 0 1 #_fleet:5_FISHERY_Rec2
0.0005 1e-07 0 0 0 0 1 #_fleet:6_CPUE_BT
0.0005 1e-07 0 0 0 0 1 #_fleet:7_CPUE_BPT
0.0005 1e-07 0 0 0 0 1 #_fleet:8_TrawSurvey
0.0005 1e-07 0 0 0 0 1 #_fleet:9_TrawSurvey2yr
0.0005 1e-07 0 0 0 0 1 #_fleet:10_TrawSurvey3yr
0.0005 1e-07 0 0 0 0 1 #_fleet:11_TrawSurvey4yr
0.0005 1e-07 0 0 0 0 1 #_fleet:12_TrawSurvey5yr
0.0005 1e-07 0 0 0 0 1 #_fleet:13_TrawSurvey6_14yr# sex codes:  0=combined; 1=use female only; 2=use male only; 3=use both as joint sexxlength distribution
# partition codes:  (0=combined; 1=discard; 2=retained
25 #_N_LengthBins; then enter lower edge of each length bin
 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 68 72 76 78 80 
#_yr month fleet sex part Nsamp datavector(female-male)
 1971 7 1 3 0 125 0	0	0	0	0	0	0	0	0	4	1	1	2	4	1	5	6	2	3	11	8	4	5	0	0	
 1971 3 1 3 0 125 0	0	0	0	0	0	0	0	0	4	1	1	2	4	1	5	6	2	3	11	8	4	5	0	0	
 1971 12 1 3 0 125 0	0	0	0	0	0	0	0	0	4	1	1	2	4	1	5	6	2	3	11	8	4	5	0	0	
-9999	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	
#
30 #_N_age_bins
 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
2 #_N_ageerror_definitions
0.5	1.5	2.5	3.5	4.5	5.5	6.5	7.5	8.5	9.5	10.5	11.5	12.5	13.5	14.5	15.5	16.5	17.5	18.5	19.5	20.5	21.5	22.5	23.5	24.5	25.5	26.5	27.5	28.5	29.5	30.5	31.5	32.5	33.5	34.5	35.5	36.5	37.5	38.5	39.5	40.5	41.5	42.5	43.5	44.5	45.5	46.5	47.5	48.5	49.5	50.5	51.5	52.5	53.5	54.5	55.5	56.5	57.5	58.5	59.5	60.5	61.5	62.5	63.5	64.5	65.5	66.5	67.5	68.5	69.5	70.5	71.5	72.5	73.5	74.5	75.5	76.5	77.5	78.5	79.5	80.5
0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001
0.5	1.5	2.5	3.5	4.5	5.5	6.5	7.5	8.5	9.5	10.5	11.5	12.5	13.5	14.5	15.5	16.5	17.5	18.5	19.5	20.5	21.5	22.5	23.5	24.5	25.5	26.5	27.5	28.5	29.5	30.5	31.5	32.5	33.5	34.5	35.5	36.5	37.5	38.5	39.5	40.5	41.5	42.5	43.5	44.5	45.5	46.5	47.5	48.5	49.5	50.5	51.5	52.5	53.5	54.5	55.5	56.5	57.5	58.5	59.5	60.5	61.5	62.5	63.5	64.5	65.5	66.5	67.5	68.5	69.5	70.5	71.5	72.5	73.5	74.5	75.5	76.5	77.5	78.5	79.5	80.5
0.5	0.65	0.67	0.7	0.73	0.76	0.8	0.84	0.88	0.92	0.97	1.03	1.09	1.16	1.23	1.32	1.41	1.51	1.62	1.75	1.89	2.05	2.23	2.45	2.71	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3
#_mintailcomp: upper and lower distribution for females and males separately are accumulated until exceeding this level.
#_addtocomp:  after accumulation of tails; this value added to all bins
#_combM+F: males and females treated as combined gender below this bin number 
#_compressbins: accumulate upper tail by this number of bins; acts simultaneous with mintailcomp; set=0 for no forced accumulation
#_Comp_Error:  0=multinomial, 1=dirichlet
#_ParmSelect:  parm number for dirichlet
#_minsamplesize: minimum sample size; set to 1 to match 3.24, minimum value is 0.001
#
#_mintailcomp addtocomp combM+F CompressBins CompError ParmSelect minsamplesize
0.0005 1e-07 0 0 0 0 1 #_fleet:1_FISHERY_BT
0.0005 1e-07 0 0 0 0 1 #_fleet:2_FISHERY_BPT
0.0005 1e-07 0 0 0 0 1 #_fleet:3_FISHERY_JP
0.0005 1e-07 0 0 0 0 1 #_fleet:4_FISHERY_Rec1
0.0005 1e-07 0 0 0 0 1 #_fleet:5_FISHERY_Rec2
0.0005 1e-07 0 0 0 0 1 #_fleet:6_CPUE_BT
0.0005 1e-07 0 0 0 0 1 #_fleet:7_CPUE_BPT
0.0005 1e-07 0 0 0 0 1 #_fleet:8_TrawSurvey
0.0005 1e-07 0 0 0 0 1 #_fleet:9_TrawSurvey2yr
0.0005 1e-07 0 0 0 0 1 #_fleet:10_TrawSurvey3yr
0.0005 1e-07 0 0 0 0 1 #_fleet:11_TrawSurvey4yr
0.0005 1e-07 0 0 0 0 1 #_fleet:12_TrawSurvey5yr
0.0005 1e-07 0 0 0 0 1 #_fleet:13_TrawSurvey6_14yr
1 #_Lbin_method_for_Age_Data: 1=poplenbins; 2=datalenbins; 3=lengths
# sex codes:  0=combined; 1=use female only; 2=use male only; 3=use both as joint sexxlength distribution
# partition codes:  (0=combined; 1=discard; 2=retained
#_yr month fleet sex part ageerr Lbin_lo Lbin_hi Nsamp datavector(female-male)
 1971 7 1 3 0 2 1 -1 75 0 0 0 0 3 1 1 4 2 1 0 1 2 2 13 2 3 0 0 4 2 1 1 2 1 2 2 1 2 1
 1971 1 1 3 0 2 1 -1 75 0 0 0 0 3 1 1 4 2 1 0 1 2 2 13 2 3 0 0 4 2 1 1 2 1 2 2 1 2 1
 1971 9 1 3 0 2 1 -1 75 0 0 0 0 3 1 1 4 2 1 0 1 2 2 13 2 3 0 0 4 2 1 1 2 1 2 2 1 2 1
-9999  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
#
0 #_Use_MeanSize-at-Age_obs (0/1)
# sex codes:  0=combined; 1=use female only; 2=use male only; 3=use both as joint sexxlength distribution
# partition codes:  (0=combined; 1=discard; 2=retained
# ageerr codes:  positive means mean length-at-age; negative means mean bodywt_at_age
#_yr month fleet sex part ageerr ignore datavector(female-male)
#                                          samplesize(female-male)
#
0 #_N_environ_variables
# -2 in yr will subtract mean for that env_var; -1 will subtract mean and divide by stddev (e.g. Z-score)
#Yr Variable Value
#
0 # N sizefreq methods to read 
#
0 # do tags (0/1)
#
0 #    morphcomp data(0/1) 
#  Nobs, Nmorphs, mincomp
#  yr, seas, type, partition, Nsamp, datavector_by_Nmorphs
#
0  #  Do dataread for selectivity priors(0/1)
# Yr, Seas, Fleet,  Age/Size,  Bin,  selex_prior,  prior_sd
# feature not yet implemented
#
999

ENDDATA
